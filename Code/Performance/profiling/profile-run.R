# Simulation run to be profiled by profile.ps1: writes its process id, waits 4 s for the profiler, then runs.
# usage: Rscript profile-run.R <setting from settings.R> <pid file> <library folder> [generations]
args <- commandArgs(TRUE)
library(EcoPhyloSim, lib.loc = args[3])
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "settings.R"))
gens <- if (length(args) >= 4) as.integer(args[4]) else 2000
par <- settingPar(args[1], runs = gens, seed = 1)
writeLines(as.character(Sys.getpid()), args[2])
Sys.sleep(4)  # time for the profiler to attach
t <- system.time(sim <- runSimulation(par))
cat(sprintf("RESULT %s gens=%d elapsed=%.2f\n", args[1], gens, t[3]))
