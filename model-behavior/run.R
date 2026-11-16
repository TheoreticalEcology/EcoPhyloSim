# Runs the simulations of design.R with one model version and saves the
# indicators of each simulation (indicators.R) in a cache folder, one .rds file
# per simulation; finished simulations are skipped, so an interrupted run resumes.
# The cache is turned into an archive with archive.R. See README.md.
#
# Run from the repository root:
#   Rscript model-behavior/run.R <mode> <cores> <version> <cache folder> [lib=<libraries>]
#     [commit=<commit of the simulated package>] [<blocks>] [part=i/n]
# <mode>: pilot, main or test (design.R); <version>: a name in `versions` of
# adapters.R. The current EcoPhyloSim must be installed (indicators). Without
# lib=, the simulations also use it; with lib= (the library of the version,
# followed by the libraries of its dependencies, separated by ";" on Windows and
# ":" otherwise), each simulation runs in its own R process with that package.
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
commit <- sub("^commit=", "", grep("^commit=", rest, value = TRUE))
part <- grep("^part=", rest, value = TRUE)
blocks <- rest[!grepl("^(lib|commit|part)=", rest)]
if (is.null(versions[[version]])) stop("unknown version ", version)
if (version != "current" && !length(lib)) stop("versions other than current need lib=<libraries>")

runTask <- function(task) {
  minFree <- as.numeric(Sys.getenv("MB_MIN_FREE_GB", "30"))
  freeGB <- function() if (.Platform$OS.type != "windows") NA else
    tryCatch(as.numeric(system2("powershell", c("-NoProfile", "-Command",
      "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"), stdout = TRUE)) / 1024^2,
      error = function(e) NA)
  while (isTRUE(freeGB() < minFree)) Sys.sleep(60)
  running <- sub("\\.rds$", ".running", task$file)
  file.create(running)
  on.exit(unlink(running))
  out <- tryCatch({
    par <- adapt(task, version)
    if (is.null(par)) list(unsupported = TRUE) else {
      if (length(lib)) {
        # simulate in a separate R process with the package of the other version;
        # system2(env = ) does not work on Windows, and each worker is its own process
        tmp <- tempfile(c("par", "sim"), fileext = ".rds")
        on.exit(unlink(tmp), add = TRUE)
        saveRDS(par, tmp[1])
        Sys.setenv(R_LIBS = lib)
        log <- sub("\\.rds$", ".log", task$file)
        status <- system2(file.path(R.home("bin"), "Rscript"),
                          c("model-behavior/simulate.R", versions[[version]]$package, tmp[1], tmp[2]),
                          stdout = log, stderr = log)
        if (status != 0 || !file.exists(tmp[2])) stop("simulation failed (exit status ", status, "), see ", log)
        unlink(log)
        sim <- readRDS(tmp[2])
      } else {
        capture.output(sim <- suppressWarnings(runSimulation(do.call(createCompletePar, par))))
        sim$packageVersion <- paste("EcoPhyloSim", packageVersion("EcoPhyloSim"))
      }
      set.seed(task$seed) # R-side randomness of the indicators (sac(), nullModel())
      r <- behavior(sim, censusTimes(task$generations))
      if (task$rep != 1) r$specMat <- NULL # species map of the first replicate only
      # parameters actually used, including the defaults of the version
      model <- sim$Model[vapply(sim$Model, function(v) is.atomic(v) && length(v) == 1, logical(1))]
      c(r, list(runtime = sim$Model$runtime, package = sim$packageVersion,
                parameters = paste(names(model), model, sep = "=", collapse = "; ")))
    }
  }, error = function(e) list(error = conditionMessage(e)))
  for (attempt in 1:5) { # retry if the file system is briefly unavailable
    if (!inherits(try(saveRDS(c(list(task = task), out), task$file), silent = TRUE), "try-error")) break
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
# finished simulations are skipped, but only if they were run with the same parameters
done <- file.exists(tasks$file)
for (i in which(done)) {
  old <- readRDS(tasks$file[i])$task
  if (!isTRUE(all.equal(old[names(base)], tasks[i, names(base)], check.attributes = FALSE)))
    stop(tasks$file[i], " was run with other parameters; use a new cache folder")
}
tasks <- tasks[!done, ]
dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)
message(nrow(tasks), " simulations to run on ", nCores, " cores")
# record the run and keep copies of the scripts in the cache (archived by archive.R)
stamp <- format(Sys.time(), "%Y%m%d-%H%M")
for (f in file.path("model-behavior", c("design.R", "indicators.R", "adapters.R", "run.R", "simulate.R")))
  file.copy(f, file.path(cacheDir, sub("\\.R$", paste0("-", stamp, ".R"), basename(f))))
cat(sprintf("%s: run %s, version %s, %d simulations on %d cores, simulation libraries %s, indicators with EcoPhyloSim %s\n",
            format(Sys.time()), mode, version, nrow(tasks), nCores, if (length(lib)) lib else "default",
            packageVersion("EcoPhyloSim")),
    sprintf("simulated package commit: %s\n", if (length(commit)) commit else "not given (default library)"),
    sprintf("scripts commit: %s (git describe at the start)\n",
            tryCatch(system("git describe --always --dirty", intern = TRUE), error = function(e) "unknown")),
    file = file.path(cacheDir, "RUNINFO.txt"), append = TRUE, sep = "")
cl <- parallel::makeCluster(nCores)
parallel::clusterExport(cl, c("sources", "version", "lib"))
parallel::clusterEvalQ(cl, { library(EcoPhyloSim); for (f in sources) source(f) })
# longest tasks first, for better load balancing (rough cost: cells x dispersal
# kernel size x competition neighborhood size)
kernel <- ifelse(tasks$dispersal == "global", 0.1, 1 + (2 * suppressWarnings(as.numeric(tasks$dispersal)))^2)
neighborhood <- ifelse(tasks$density > 0, 1 + tasks$densityCut^2, 1)
tasks <- tasks[order(-tasks$generations * tasks$x * tasks$y * kernel * neighborhood), ]
# chunk.size = 1: each worker takes the next single task
invisible(parallel::parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))), runTask, chunk.size = 1))
parallel::stopCluster(cl)
