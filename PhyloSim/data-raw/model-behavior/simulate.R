# Runs one simulation with the package of another model version (called by
# run.R with R_LIBS set to the library of that version, followed by the
# libraries of its dependencies):
#   Rscript PhyloSim/data-raw/model-behavior/simulate.R <package> <parameter file> <output file>
# The package is loaded from these libraries only, never from the default
# library, and its version is saved with the result.
args <- commandArgs(trailingOnly = TRUE)
libs <- strsplit(Sys.getenv("R_LIBS"), .Platform$path.sep)[[1]]
suppressPackageStartupMessages(library(args[1], character.only = TRUE, lib.loc = libs))
par <- do.call(createCompletePar, readRDS(args[2]))
set.seed(par$seed) # R-side randomness of versions before #174 (random resolution of polytomies)
capture.output(sim <- suppressWarnings(runSimulation(par)))
sim$packageVersion <- paste(args[1], utils::packageDescription(args[1], lib.loc = libs)$Version)
saveRDS(sim, args[3])
