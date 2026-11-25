# #143: richness with protracted speciation, before and after offspring inherit
# the incipient state. Neutral model, global dispersal, 30 x 30 cells,
# specRate = 4, 3000 generations, seeds 1:10, one run per seed; fission 0, 1
# and 2. Compared with the Ewens expectation for theta = specRate * exp(-p),
# i.e. if only the founding individual can complete speciation.
# Results: see README.md. Run with R 4.6.2 and EcoPhyloSim at commit 619b469
# (before) and at commit 7d72cdf (after); the library
# with the installed version is the first command line argument, e.g.
# Rscript richness.R <library>. Single thread, about 10 minutes.
lib <- commandArgs(TRUE)[1]
library(EcoPhyloSim, lib.loc = lib)
ewens <- function(theta, n) sum(theta / (theta + 0:(n - 1)))
res <- expand.grid(protracted = c(0, 1, 2, 5), fission = 0:2)
res[c("mean", "sd")] <- t(apply(res, 1, function(s) {
  rich <- sapply(1:10, function(seed) {
    par <- createCompletePar(x = 30, y = 30, specRate = 4, runs = 3000, seed = seed,
                             protracted = s[["protracted"]], fission = s[["fission"]])
    invisible(capture.output(simu <- runSimulation(par)))
    length(unique(c(simu$Output[[1]]$specMat)))
  })
  c(mean(rich), sd(rich))
}))
res$ewensFounder <- sapply(res$protracted, function(p) ewens(4 * exp(-p), 900))
print(res, digits = 3)
