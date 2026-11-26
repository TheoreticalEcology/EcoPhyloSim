# #141: results with global dispersal and with dispersal at multiples of 0.5
# must not change. Runs 14 settings (seed 3, censuses at 50, 150 and 300, longRun 500 and 2000
# generations) and saves species, traits and phylogeny of every census. Run
# once with the version before the fix and once with the fix, each installed in
# its own library, then compare:
#   Rscript identity.R <libraryBefore> before.rds
#   Rscript identity.R <libraryAfter> after.rds
#   Rscript -e 'mapply(identical, readRDS("before.rds"), readRDS("after.rds"))'
# Checked with R 4.6.2, master 1c778ea against branch i141-dispersal: all TRUE.
# Single thread, about a minute.
args <- commandArgs(TRUE)
library(EcoPhyloSim, lib.loc = args[1])
base <- list(x = 24, y = 18, runs = c(50, 150, 300), seed = 3, specRate = 3)
settings <- list(
  global = list(),
  globalEnv = list(environment = 1, negativeDD = 0.5),
  d0.5 = list(dispersal = 0.5),
  d1 = list(dispersal = 1),
  d1.5 = list(dispersal = 1.5),
  d2 = list(dispersal = 2),
  d4 = list(dispersal = 4),
  dd = list(dispersal = 1.5, negativeDD = 1, negativeDDCut = 2),
  env = list(dispersal = 2.5, environment = 1),
  repro = list(dispersal = 1, environment = 0.5, negativeDD = 0.5, fitnessActsOn = "both"),
  recr = list(dispersal = 3, environment = 1, negativeDD = 0.5, fitnessActsOn = c("mortality", "recruitment")),
  fission = list(dispersal = 0.5, fission = 2),
  posDD = list(dispersal = 1, positiveDD = 0.1, positiveDDCut = 2),
  longRun = list(x = 30, y = 30, runs = c(500, 2000), specRate = 2, dispersal = 1, seed = 7)
)
res <- lapply(settings, function(s) {
  p <- modifyList(base, s)
  invisible(capture.output(simu <- suppressMessages(suppressWarnings(
    runSimulation(do.call(createCompletePar, p), progress = FALSE)))))
  lapply(simu$Output, function(o) o[c("specMat", "traitMat", "envMat", "compMat", "neutMat", "phyloTXT")])
})
saveRDS(res, args[2])
