# Runs the simulations of design.R with one model version and saves the
# indicators of each simulation (indicators.R) in a cache folder, one .rds file
# per simulation; finished simulations are skipped, so an interrupted run resumes.
# The cache is turned into an archive with archive.R.
#
# Run from the repository root:
#   Rscript model-behavior/run.R <mode> <cores> <version> <cache folder> [lib=<library>]
#     [commit=<commit of the simulated package>] [<blocks>] [part=i/n]
# <mode>: pilot, main or test (design.R); <version>: 0.2, 0.3, 0.3.1 or current
# (see adapters.R). The current EcoPhyloSim must be installed (indicators).
# Without lib=, the simulations also use it; with lib=, each simulation runs in
# its own R process with the package of that library (PhyloSim for 0.2 and 0.3,
# EcoPhyloSim otherwise). Use a new cache folder for each version: finished
# simulations are recognized by their file names only.

library(EcoPhyloSim)
for (f in c("design.R", "indicators.R", "adapters.R")) source(file.path("model-behavior", f))
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
if (version != "current" && !length(lib)) stop("versions other than current need lib=<library>")

runTask <- function(task) {
  # wait until at least 30 GB of memory are free (Windows only, otherwise no check)
  freeGB <- function() tryCatch(as.numeric(system2("powershell", c("-NoProfile", "-Command",
    "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"), stdout = TRUE)) / 1024^2,
    error = function(e) NA)
  while (isTRUE(freeGB() < 30)) Sys.sleep(60)
  running <- sub("\\.rds$", ".running", task$file)
  file.create(running)
  on.exit(unlink(running))
  par <- adapt(task, version)
  out <- if (is.null(par)) list(unsupported = TRUE) else tryCatch({
    if (length(lib)) {
      # simulate in a separate R process with the package of the other version
      tmp <- tempfile(c("par", "sim"), fileext = ".rds")
      on.exit(unlink(tmp), add = TRUE)
      saveRDS(par, tmp[1])
      pkg <- if (version %in% c("0.2", "0.3")) "PhyloSim" else "EcoPhyloSim"
      # system2(env = ) does not work on Windows; each worker is its own process
      Sys.setenv(R_LIBS = lib)
      log <- sub("\\.rds$", ".log", task$file)
      status <- system2(file.path(R.home("bin"), "Rscript"),
                        c("model-behavior/simulate.R", pkg, tmp[1], tmp[2]), stdout = log, stderr = log)
      if (status != 0 || !file.exists(tmp[2])) stop("simulation failed (exit status ", status, "), see ", log)
      unlink(log)
      sim <- readRDS(tmp[2])
    } else {
      capture.output(sim <- suppressWarnings(runSimulation(do.call(createCompletePar, par))))
    }
    set.seed(task$seed) # R-side randomness of the indicators (sac(), nullModel())
    r <- behavior(sim, censusTimes(task$generations))
    if (task$rep != 1) r$specMat <- NULL # species map of the first replicate only
    c(r, list(runtime = sim$Model$runtime))
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
tasks <- tasks[!file.exists(tasks$file), ]
dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)
message(nrow(tasks), " simulations to run on ", nCores, " cores")
# record the run and keep copies of the scripts in the cache (archived by archive.R)
stamp <- format(Sys.time(), "%Y%m%d-%H%M")
for (f in c("design.R", "indicators.R", "adapters.R", "run.R", "simulate.R"))
  file.copy(file.path("model-behavior", f), file.path(cacheDir, sub("\\.R$", paste0("-", stamp, ".R"), f)))
cat(sprintf("%s: run %s, version %s, %d simulations on %d cores, simulation library %s, indicators with EcoPhyloSim %s\n",
            format(Sys.time()), mode, version, nrow(tasks), nCores, if (length(lib)) lib else "default",
            packageVersion("EcoPhyloSim")),
    sprintf("simulated package commit: %s\n", if (length(commit)) commit else "not given (default library)"),
    sprintf("scripts commit: %s (git HEAD at the start)\n",
            tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) "unknown")),
    file = file.path(cacheDir, "RUNINFO.txt"), append = TRUE, sep = "")
cl <- parallel::makeCluster(nCores)
parallel::clusterExport(cl, c("version", "lib", "base", "censusTimes", "adapt", "behavior", "censusStats",
                              "distanceDecay", "cndd"))
parallel::clusterEvalQ(cl, library(EcoPhyloSim))
# longest tasks first, for better load balancing (rough cost: cells x dispersal
# kernel size x competition neighborhood size)
kernel <- ifelse(tasks$dispersal == "global", 0.1, 1 + (2 * suppressWarnings(as.numeric(tasks$dispersal)))^2)
neighborhood <- ifelse(tasks$density > 0, 1 + tasks$densityCut^2, 1)
tasks <- tasks[order(-tasks$generations * tasks$x * tasks$y * kernel * neighborhood), ]
# chunk.size = 1: each worker takes the next single task
invisible(parallel::parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))), runTask, chunk.size = 1))
parallel::stopCluster(cl)
