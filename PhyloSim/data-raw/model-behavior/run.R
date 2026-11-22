# Runs the simulations of design.R with one model version and saves the
# indicators of each simulation (indicators.R) in a cache folder, one .rds file
# per simulation, and its raw output (maps of the last six censuses, phylogeny) in
# raw/; finished simulations are skipped (failed ones are rerun), so an
# interrupted run resumes. A complete run of mode main without failed simulations
# is archived with archive.R. See README.md.
#
# Run from the repository root:
#   Rscript PhyloSim/data-raw/model-behavior/run.R <mode> <cores> <version> <commit or tag> [<blocks>] [part=i/n]
# <mode>: main or test (design.R); <version>: a name in `versions` of adapters.R.
# The commit and the EcoPhyloSim that computes the indicators (`indicatorCommit` of
# design.R) are installed into lib/ if missing; the cache is
# cache/<mode>_<version>_<commit>/ (folders in the main checkout, see README.md).
# Each simulation runs in its own R process with the package of the commit.
# Environment variable MB_MIN_FREE_GB (default 16): free memory to wait for
# before each simulation (checked on Windows only).

source("PhyloSim/data-raw/model-behavior/common.R")
sources <- file.path(mb, c("common.R", "design.R", "adapters.R"))
for (f in sources) source(f)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: run.R <mode> <cores> <version> <commit or tag> [<blocks>] [part=i/n]")
mode <- args[1]
nCores <- as.integer(args[2])
version <- args[3]
rest <- args[-(1:4)]
part <- grep("^part=", rest, value = TRUE)
blocks <- rest[!grepl("^part=", rest)]
if (is.null(modes[[mode]])) stop("unknown mode ", mode)
if (is.null(versions[[version]])) stop("unknown version ", version)

# the EcoPhyloSim of the indicators, the same for all versions; R_LIBS is set
# before the workers start, which do not see .libPaths() of this session (#190)
indLib <- libraryOf(indicatorCommit)
Sys.setenv(R_LIBS = indLib)
.libPaths(c(indLib, .libPaths()))
library(EcoPhyloSim, lib.loc = indLib)
sources <- c(sources, file.path(mb, "indicators.R"))
source(sources[4])
indicators <- installInfo(indLib, "commit")
# the simulated commit, its package and the name of its archive (release number,
# or the package version for a development version)
lib <- libraryOf(args[4], isTRUE(versions[[version]]$oldDeps))
commit <- installInfo(lib, "commit")
pkg <- strsplit(installInfo(lib, "package"), " ")[[1]]
package <- pkg[1]
archiveName <- paste0(if (mode != "main") paste0(mode, "-"), if (version == "current") pkg[2] else version,
                      "_", substr(commit, 1, 7))
cacheDir <- cacheOf(mode, version, commit)

runTask <- function(task) {
  minFree <- as.numeric(Sys.getenv("MB_MIN_FREE_GB", "16"))
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
      status <- system2(Rscript, shQuote(c(simScript, package, tmp[1], tmp[2])), stdout = log, stderr = log)
      if (status != 0 || !file.exists(tmp[2])) stop("simulation failed (exit status ", status, "), see ", log)
      unlink(log)
      sim <- readRDS(tmp[2])
      set.seed(task$seed) # R-side randomness of the indicators (sac(), nullModel())
      r <- behavior(sim, censusTimes(task$generations))
      # raw output, kept only locally (not in the archive tables), enough to recompute the
      # indicators at the end: all maps of the last six censuses (CNDD) and the last phylogeny
      last <- length(sim$Output)
      lastSix <- max(1, last - 5):last
      maps <- c("specMat", "traitMat", "envMat", "compMat", "neutMat")
      raw <- list(runs = sim$Model$runs[lastSix], Output = lapply(lastSix, function(i)
        sim$Output[[i]][intersect(c(maps, if (i == last) "phylogeny"), names(sim$Output[[i]]))]))
      dir.create(file.path(dirname(task$file), "raw"), showWarnings = FALSE)
      saveRDS(raw, file.path(dirname(task$file), "raw", basename(task$file)), compress = "xz")
      r$specMat <- NULL
      # parameters actually used, including the defaults of the version
      model <- sim$Model[vapply(sim$Model, function(v) is.atomic(v) && length(v) == 1, logical(1))]
      c(r, list(runtime = sim$Model$runtime, package = sim$packageVersion, indicators = indicators,
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

allTasks <- tasks <- taskTable(mode, cacheDir)
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
    stop(tasks$file[i], " was run with other parameters, generations or seed; delete the cache folder")
  if (!is.null(old$error)) done[i] <- FALSE
}
tasks <- tasks[!done, ]
scripts <- file.path(mb, c("common.R", "design.R", "indicators.R", "adapters.R", "run.R", "simulate.R",
                           "archive.R", "installVersion.R", "installOldDeps.R"))
if (mode == "main") {
  # runs of record only from committed scripts in a git checkout, so that RUNINFO.txt names them
  if (gitDescribe() == "unknown" || length(suppressWarnings(system(
      paste("git status --porcelain --", paste(scripts, collapse = " ")), intern = TRUE))))
    stop("mode main runs only from a git checkout with committed scripts (", mb, ")")
  # a resumed run uses the scripts of its first start (copies in the cache)
  first <- sort(list.files(cacheDir, "^run-.*\\.R$"))[1]
  if (!is.na(first)) {
    copies <- file.path(cacheDir, sub("\\.R$", sub("^run", "", first), basename(scripts)))
    if (!all(file.exists(copies)) || any(tools::md5sum(copies) != tools::md5sum(scripts)))
      stop("the scripts differ from those of the first start of the run in ", cacheDir, "; delete the cache folder")
  }
}
dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)
message(nrow(tasks), " simulations to run on ", nCores, " cores in ", cacheDir)
# record the run and keep copies of the scripts in the cache (archived by archive.R)
stamp <- format(Sys.time(), "%Y%m%d-%H%M")
for (f in scripts) file.copy(f, file.path(cacheDir, sub("\\.R$", paste0("-", stamp, ".R"), basename(f))))
simScript <- file.path(cacheDir, paste0("simulate-", stamp, ".R")) # the copy, not the checkout
cat(sprintf("%s: run %s, version %s, %d simulations on %d cores, simulation libraries %s\n",
            format(Sys.time()), mode, version, nrow(tasks), nCores, lib),
    sprintf("simulated package: %s, commit %s (INSTALLINFO.txt)\n", paste(pkg, collapse = " "), commit),
    sprintf("indicators: EcoPhyloSim %s, commit %s\n", packageVersion("EcoPhyloSim"), indicators),
    sprintf("scripts commit: %s (git describe at the start)\n", gitDescribe()),
    sprintf("archive: %s\n", archiveName),
    file = file.path(cacheDir, "RUNINFO.txt"), append = TRUE, sep = "")
# provenance: command, R, platform, and the packages of all libraries used
libs <- unique(c(strsplit(lib, .Platform$path.sep)[[1]], .libPaths()))
info <- file.path(libs, "INSTALLINFO.txt")
writeLines(c(paste("command: Rscript", file.path(mb, "run.R"), paste(args, collapse = " ")),
             paste("platform:", R.version$platform, "-", paste(Sys.info()[c("sysname", "release")], collapse = " ")),
             "", "Session of the indicators:", capture.output(sessionInfo()), "",
             "Installations of the libraries used (INSTALLINFO.txt of installVersion.R):",
             unlist(lapply(info[file.exists(info)], readLines)), "",
             "Packages in the libraries (searched in this order by the simulations):",
             capture.output(print(installed.packages(libs)[, c("Package", "Version", "LibPath")],
                                  quote = FALSE, row.names = FALSE))),
           file.path(cacheDir, paste0("SESSION-", stamp, ".txt")))
cl <- parallel::makeCluster(nCores)
parallel::clusterExport(cl, c("sources", "version", "lib", "package", "simScript", "indicators"))
parallel::clusterEvalQ(cl, { library(EcoPhyloSim); for (f in sources) source(f) })
# longest tasks first, for better load balancing (rough cost: cells x dispersal
# kernel size x competition neighborhood size)
kernel <- ifelse(tasks$dispersal == "global", 0.1, 1 + (2 * suppressWarnings(as.numeric(tasks$dispersal)))^2)
neighborhood <- ifelse(tasks$density > 0, 1 + tasks$densityCut^2, 1)
tasks <- tasks[order(-tasks$generations * tasks$x * tasks$y * kernel * neighborhood), ]
# chunk.size = 1: each worker takes the next single task
invisible(parallel::parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))), runTask, chunk.size = 1))
parallel::stopCluster(cl)

# archive a complete run of mode main without failed simulations (archive.R as at the start)
if (mode == "main" && !length(blocks) && !length(part)) {
  failed <- sum(vapply(allTasks$file, function(f) !file.exists(f) || !is.null(tryCatch(readRDS(f)$error, error = function(e) "unreadable")), logical(1)))
  if (failed) message(failed, " simulations failed; rerun run.R to retry them, or archive the run as it is with\n",
                      "  Rscript ", file.path(mb, "archive.R"), " ", cacheDir) else
    system2(Rscript, shQuote(c(file.path(mb, "archive.R"), cacheDir)))
}
