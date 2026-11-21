# Runs one setting with a library built with prototypes/counters.patch, which prints the counters (events, deaths,
# parent draws, recruitment fallbacks, density updates, localDissimilarity and recruitFitness calls) when the
# model object is deleted. Not instrumented: the neutral branch with global dispersal (all counters 0).
# usage: Rscript count-run.R <setting from settings.R> <generations> <library folder>
a <- commandArgs(TRUE)
library(EcoPhyloSim, lib.loc = a[3])
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "settings.R"))
par <- settingPar(a[1], runs = as.integer(a[2]), seed = 1)
t <- system.time(invisible(capture.output(sim <- runSimulation(par), type = "message")))
rm(sim); invisible(gc())  # deletes the model object, which prints the counters
cat(sprintf("RESULT %s cpu=%.2f\n", a[1], t[["user.self"]]))
