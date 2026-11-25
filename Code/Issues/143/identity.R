# #143: results with protracted = 0 must not change. Runs 22 settings with
# protracted = 0 (seed 3, censuses at 50, 150 and 300 generations; longRun
# seed 7) and saves species, traits and phylogeny of every census. Run once
# with the version before the fix and once with the fix, each installed in
# its own library, then compare:
#   Rscript identity.R <libraryBefore> before.rds
#   Rscript identity.R <libraryAfter> after.rds
#   Rscript -e 'mapply(identical, readRDS("before.rds"), readRDS("after.rds"))'
# Checked with R 4.6.2, master 116f180 against its merge with 7d72cdf: all TRUE.
# Single thread, a few minutes.
args <- commandArgs(TRUE)
library(EcoPhyloSim, lib.loc = args[1])
base <- list(x = 24, y = 18, runs = c(50, 150, 300), seed = 3, specRate = 3)
settings <- list(
  global = list(),
  local1 = list(dispersal = 1),
  local2 = list(dispersal = 2),
  fission1 = list(dispersal = 1, fission = 1),
  fission2 = list(dispersal = 1, fission = 2),
  fission1global = list(fission = 1),
  fission2global = list(fission = 2, specRate = 20),
  dd = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2),
  ddFission = list(dispersal = 1, negativeDD = 1, fission = 1, negativeDDKernel = "species"),
  env = list(dispersal = 1, environment = 1),
  envRecr = list(dispersal = 1, environment = 1, negativeDD = 0.5, fitnessActsOn = c("mortality", "recruitment")),
  jump = list(dispersal = 1, traitSpeciationJump = 0.2, environment = 0.5),
  jumpFission = list(dispersal = 1, traitSpeciationJump = c(0.1, 0.2, 0.3), fission = 2),
  redQueen = list(dispersal = 1, redQueen = 0.1, redQueenStrength = 1),
  redQueenFission = list(redQueen = 0.1, redQueenStrength = 5, fission = 1),
  highRate = list(specRate = 60),
  highRateLocal = list(dispersal = 1, specRate = 60, traitSpeciationJump = 0.1),
  highRateSmall = list(x = 4, y = 4, specRate = 12, fission = 0),
  fullTree = list(dispersal = 1, prunePhylogeny = FALSE, specRate = 5),
  fullTreeFission = list(prunePhylogeny = FALSE, fission = 2, specRate = 5),
  posDD = list(dispersal = 1, positiveDD = 0.1, positiveDDCut = 2),
  longRun = list(x = 30, y = 30, runs = c(500, 2000), specRate = 2, dispersal = 1, seed = 7)
)
res <- lapply(settings, function(s) {
  p <- modifyList(base, s)
  p$protracted <- 0
  invisible(capture.output(simu <- suppressMessages(suppressWarnings(
    runSimulation(do.call(createCompletePar, p))))))
  lapply(simu$Output, function(o) o[c("specMat", "traitMat", "envMat", "compMat", "neutMat", "phyloTXT")])
})
saveRDS(res, args[2])
