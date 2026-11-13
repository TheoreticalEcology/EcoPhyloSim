# Hand-built PhyloSim object for tests of the analysis functions without a simulation.
# Each element of 'mats' is the species matrix of one recorded generation.
fakeSimulation <- function(mats) {
  if (is.matrix(mats)) mats <- list(mats)
  structure(list(Model = list(x = ncol(mats[[1]]), y = nrow(mats[[1]]), runs = seq_along(mats)),
                 Output = lapply(mats, function(m) list(specMat = m, traitMat = m, envMat = m, compMat = m))),
            class = "PhyloSim")
}

# Simulation without the console output of runSimulation()
quietSimulation <- function(...) {
  out <- NULL
  invisible(capture.output(out <- runSimulation(createCompletePar(...))))
  out
}
