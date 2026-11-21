# Runs short simulations covering the main code paths and saves their complete output, so that two builds can be
# compared with identity-compare.R (bit-identical results for the same seed).
# usage: Rscript identity.R <library folder> <output.rds> [long]
# Default: 16 settings x 2 seeds, censuses at 150 and 300 generations; "long": 7 density-dependence settings,
# seed 3, censuses at 1000 and 2000 generations.
a <- commandArgs(TRUE)
suppressMessages(library(EcoPhyloSim, lib.loc = a[1]))
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "settings.R"))
long <- length(a) >= 3 && a[3] == "long"
names <- if (long) c("dd2", "dd2species", "posDD", "dd2recruit", "envDDGlobal", "fissionDD", "reproDD") else
  c("neutralLocal1", "neutralLocal3", "neutralGlobal", "env", "dd2", "dd2gauss", "dd2species", "envDD2", "posDD",
    "dd2recruit", "envRecruit2", "repro", "ddGlobal", "fission1", "protracted", "redQueen")
out <- list()
for (n in names) for (seed in if (long) 3 else 1:2) {
  sim <- suppressWarnings(runSimulation(settingPar(n, runs = if (long) c(1000, 2000) else c(150, 300), seed = seed)))
  out[[paste(n, seed)]] <- lapply(sim$Output, function(o) o[setdiff(names(o), "runtime")])
}
saveRDS(out, a[2])
