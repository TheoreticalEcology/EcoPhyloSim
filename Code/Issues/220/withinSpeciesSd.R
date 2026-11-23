# #220: checks the approximate standard deviation of a trait within species,
# traitEvolution / sqrt(3 * traitPull * (2 - traitPull)), which the help of
# createCompletePar() states and checkPar() uses for the message on a small
# Gaussian density width. Neutral model, global dispersal, 50 x 50 cells,
# 2000 generations, seeds 1:3; SD of the neutral trait within species with at
# least 50 individuals, pooled over species. Results: see README.md.
# Run with EcoPhyloSim at commit f16f2f5 (branch i220-traitEvolution) or later,
# R 4.6.2; single thread, about 2 minutes.
library(EcoPhyloSim)
set <- list(c(0.01, 0.2), c(0.05, 0.5), c(0.05, 0.1), c(0.1, 1))
for (s in set) {
  sds <- sapply(1:3, function(seed) {
    par <- createCompletePar(x = 50, y = 50, runs = 2000, specRate = 2, seed = seed,
                             traitEvolution = s[1], traitPull = s[2])
    invisible(capture.output(simu <- runSimulation(par)))
    o <- simu$Output[[1]]
    ab <- table(o$specMat)
    big <- names(ab)[ab >= 50]
    sqrt(mean(tapply(c(o$neutMat), c(o$specMat), var)[big]))
  })
  cat(sprintf("w = %.2f, pull = %.1f: predicted %.4f, simulated %s\n", s[1], s[2],
              s[1] / sqrt(3 * s[2] * (2 - s[2])), paste(round(sds, 4), collapse = " ")))
}
