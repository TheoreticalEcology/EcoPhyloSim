# Runs the simulations of design.R with one model version and saves the
# indicators of each simulation (indicators.R) in a cache folder, one .rds file
# per simulation, and its raw output (species map, phylogeny, trait maps of the first
# replicate) in raw/; finished simulations are skipped (failed ones are rerun), so an
# interrupted run resumes. The cache is turned into an archive with archive.R.
# See README.md.
#
# Run from the repository root:
#   Rscript model-behavior/run.R <mode> <cores> <version> <cache folder> lib=<libraries>
#     [<blocks>] [part=i/n]
# <mode>: main or test (design.R); <version>: a name in `versions` of adapters.R.
# lib=: the library into which installVersion.R installed the simulated commit,
# followed by the libraries of its dependencies (separated by ";" on Windows and
# ":" otherwise); each simulation runs in its own R process with that package.
# The commit and package name are taken from INSTALLINFO.txt of that library.
# The current EcoPhyloSim must be installed in the default library (indicators).
# Environment variable MB_MIN_FREE_GB (default 30): free memory to wait for
# before each simulation (checked on Windows only).

if (!file.exists("model-behavior/design.R")) stop("run from the repository root")
library(EcoPhyloSim)
sources <- file.path("model-behavior", c("design.R", "indicators.R", "adapters.R"))
for (f in sources) source(f)
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
nCores <- as.integer(args[2])
version <- args[3]
cacheDir <- args[4]
rest <- args[-(1:4)]
lib <- sub("^lib=", "", grep("^lib=", rest, value = TRUE))
part <- grep("^part=", rest, value = TRUE)
blocks <- rest[!grepl("^(lib|part)=", rest)]
if (is.null(versions[[version]])) stop("unknown version ", version)
if (!length(lib)) stop("lib=<libraries> is needed (install the commit with installVersion.R)")
# commit and package of the simulated version, recorded by installVersion.R
info <- file.path(strsplit(lib, .Platform$path.sep)[[1]][1], "INSTALLINFO.txt")
if (!file.exists(info)) stop(info, " not found; install the commit with installVersion.R")
info <- readLines(info)
commit <- sub("^commit: ", "", grep("^commit: ", info, value = TRUE))
package <- strsplit(sub("^package: ", "", grep("^package: ", info, value = TRUE)), " ")[[1]][1]
if (length(commit) != 1 || is.na(package)) stop("commit or package missing in INSTALLINFO.txt")

# commit of the scripts, or "unknown" outside a git checkout
gitDescribe <- function() {
  d <- suppressWarnings(tryCatch(system("git describe --always --dirty", intern = TRUE), error = function(e) ""))
  if (length(d) == 1 && nzchar(d)) d else "unknown"
}

runTask <- function(task) {
  minFree <- as.numeric(Sys.getenv("MB_MIN_FREE_GB", "30"))
  freeGB <- function() if (.Platform$OS.type != "windows") NA else
    tryCatch(as.numeric(system2("powershell", c("-NoProfile", "-Command",
      "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"), stdout = TRUE)) / 1024^2,
      error = function(e) NA)
  if (isTRUE(freeGB() < minFree)) message("waiting for ", minFree, " GB of free memory (MB_MIN_FREE_GB)")
  while (isTRUE(freeGB() < minFree)) Sys.sleep(60)
  running <- sub("\\.rds$", ".running", task$file)
  file.create(running)
  on.exit(unlink(running))
  out <- tryCatch({
    par <- adapt(task, version)
    if (is.null(par)) list(unsupported = TRUE) else {
      # simulate in a separate R process with the package of the simulated version;
      # system2(env = ) does not work on Windows, and each worker is its own process
      tmp <- tempfile(c("par", "sim"), fileext = ".rds")
      on.exit(unlink(tmp), add = TRUE)
      saveRDS(par, tmp[1])
      Sys.setenv(R_LIBS = lib)
      log <- sub("\\.rds$", ".log", task$file)
      status <- system2(file.path(R.home("bin"), "Rscript"),
                        shQuote(c("model-behavior/simulate.R", package, tmp[1], tmp[2])),
                        stdout = log, stderr = log)
      if (status != 0 || !file.exists(tmp[2])) stop("simulation failed (exit status ", status, "), see ", log)
      unlink(log)
      sim <- readRDS(tmp[2])
      set.seed(task$seed) # R-side randomness of the indicators (sac(), nullModel())
      r <- behavior(sim, censusTimes(task$generations))
      # raw output of the last census, kept only locally (not in the archive tables): species map
      # and phylogeny of each replicate, trait maps of the first replicate
      final <- sim$Output[[length(sim$Output)]]
      keep <- c("specMat", "phylogeny", if (task$rep == 1) c("traitMat", "compMat", "neutMat"))
      dir.create(file.path(dirname(task$file), "raw"), showWarnings = FALSE)
      saveRDS(final[intersect(keep, names(final))], file.path(dirname(task$file), "raw", basename(task$file)),
              compress = "xz")
      r$specMat <- NULL
      # parameters actually used, including the defaults of the version
      model <- sim$Model[vapply(sim$Model, function(v) is.atomic(v) && length(v) == 1, logical(1))]
      c(r, list(runtime = sim$Model$runtime, package = sim$packageVersion,
                parameters = paste(names(model), model, sep = "=", collapse = "; ")))
    }
  }, error = function(e) list(error = conditionMessage(e)))
  # write to a temporary file first, so that an interrupted write leaves no broken file
  tmpFile <- paste0(task$file, ".tmp")
  for (attempt in 1:5) { # retry if the file system is briefly unavailable
    if (!inherits(try(saveRDS(c(list(task = task), out), tmpFile), silent = TRUE), "try-error") &&
        file.rename(tmpFile, task$file)) break
    if (attempt == 5) warning("could not save ", task$file)
    Sys.sleep(30)
  }
  invisible(NULL)
}

tasks <- taskTable(mode, cacheDir)
if (length(blocks)) tasks <- tasks[tasks$block %in% blocks, ]
if (length(part)) {
  p <- as.integer(strsplit(sub("part=", "", part), "/")[[1]])
  tasks <- tasks[seq_len(nrow(tasks)) %% p[2] == p[1] - 1, ]
}
# finished simulations are skipped, but only if they were run with the same
# parameters and seed; failed simulations and unreadable files are rerun
done <- file.exists(tasks$file)
for (i in which(done)) {
  old <- tryCatch(readRDS(tasks$file[i]), error = function(e) NULL)
  if (is.null(old)) { done[i] <- FALSE; next }
  keys <- c(names(base), "generations", "seed")
  if (!all(keys %in% names(old$task)) ||
      !isTRUE(all.equal(old$task[keys], tasks[i, keys], check.attributes = FALSE)))
    stop(tasks$file[i], " was run with other parameters, generations or seed; use a new cache folder")
  if (!is.null(old$error)) done[i] <- FALSE
}
tasks <- tasks[!done, ]
# runs of record only from committed scripts in a git checkout, so that RUNINFO.txt names them
if (mode == "main" && (gitDescribe() == "unknown" ||
    length(suppressWarnings(system("git status --porcelain -- model-behavior/*.R", intern = TRUE)))))
  stop("mode main runs only from a git checkout with committed scripts (model-behavior/*.R)")
dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)
message(nrow(tasks), " simulations to run on ", nCores, " cores")
# record the run and keep copies of the scripts in the cache (archived by archive.R)
stamp <- format(Sys.time(), "%Y%m%d-%H%M")
for (f in file.path("model-behavior", c("design.R", "indicators.R", "adapters.R", "run.R", "simulate.R")))
  file.copy(f, file.path(cacheDir, sub("\\.R$", paste0("-", stamp, ".R"), basename(f))))
cat(sprintf("%s: run %s, version %s, %d simulations on %d cores, simulation libraries %s, indicators with EcoPhyloSim %s\n",
            format(Sys.time()), mode, version, nrow(tasks), nCores, lib,
            packageVersion("EcoPhyloSim")),
    sprintf("simulated package: %s, commit %s (INSTALLINFO.txt)\n", package, commit),
    sprintf("scripts commit: %s (git describe at the start)\n",
            gitDescribe()),
    file = file.path(cacheDir, "RUNINFO.txt"), append = TRUE, sep = "")
# provenance: command, R, platform, and the packages of all libraries used
libs <- c(strsplit(lib, .Platform$path.sep)[[1]], .libPaths())
installInfo <- file.path(libs, "INSTALLINFO.txt")
writeLines(c(paste("command: Rscript model-behavior/run.R", paste(commandArgs(trailingOnly = TRUE), collapse = " ")),
             paste("platform:", R.version$platform, "-", paste(Sys.info()[c("sysname", "release")], collapse = " ")),
             "", "Session of the indicators:", capture.output(sessionInfo()), "",
             "Installations of the libraries used (INSTALLINFO.txt of installVersion.R):",
             unlist(lapply(installInfo[file.exists(installInfo)], readLines)), "",
             "Packages in the libraries (searched in this order by the simulations):",
             capture.output(print(installed.packages(libs)[, c("Package", "Version", "LibPath")],
                                  quote = FALSE, row.names = FALSE))),
           file.path(cacheDir, paste0("SESSION-", stamp, ".txt")))
cl <- parallel::makeCluster(nCores)
parallel::clusterExport(cl, c("sources", "version", "lib", "package"))
parallel::clusterEvalQ(cl, { library(EcoPhyloSim); for (f in sources) source(f) })
# longest tasks first, for better load balancing (rough cost: cells x dispersal
# kernel size x competition neighborhood size)
kernel <- ifelse(tasks$dispersal == "global", 0.1, 1 + (2 * suppressWarnings(as.numeric(tasks$dispersal)))^2)
neighborhood <- ifelse(tasks$density > 0, 1 + tasks$densityCut^2, 1)
tasks <- tasks[order(-tasks$generations * tasks$x * tasks$y * kernel * neighborhood), ]
# chunk.size = 1: each worker takes the next single task
invisible(parallel::parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))), runTask, chunk.size = 1))
parallel::stopCluster(cl)
