# Runs one simulation with the package of another model version (called by
# run.R with R_LIBS set to the library of that version):
#   Rscript model-behavior/simulate.R <package> <parameter file> <output file>
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(args[1], character.only = TRUE))
par <- do.call(createCompletePar, readRDS(args[2]))
capture.output(sim <- suppressWarnings(runSimulation(par)))
saveRDS(sim, args[3])
