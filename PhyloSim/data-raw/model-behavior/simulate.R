# Runs one simulation with the package of another model version (called by
# run.R with R_LIBS set to the library of that version, followed by the
# libraries of its dependencies):
#   Rscript PhyloSim/data-raw/model-behavior/simulate.R <package> <parameter file> <output file>
# The package is loaded from these libraries only, never from the default
# library, and its version is saved with the result.
args <- commandArgs(trailingOnly = TRUE)
libs <- strsplit(Sys.getenv("R_LIBS"), .Platform$path.sep)[[1]]
suppressPackageStartupMessages(library(args[1], character.only = TRUE, lib.loc = libs))
# v0.2 stops when it cannot read a phylogeny with a single tip (one species left); as
# in v0.3, its phylogeny is then 0 (R wrapper only, the model code is unchanged). Applied
# to package PhyloSim (v0.2 and v0.3; v0.2 also has version number 0.3), as v0.3 never
# reads such a tree.
if (args[1] == "PhyloSim") {
  readTree <- ape::read.tree
  utils::assignInNamespace("read.tree", function(file = "", text = NULL, ...)
    if (!is.null(text) && !grepl(",", text)) 0 else readTree(file = file, text = text, ...), "ape")
}
par <- do.call(createCompletePar, readRDS(args[2]))
set.seed(par$seed) # R-side randomness of versions before #174 (random resolution of polytomies)
capture.output(sim <- suppressWarnings(runSimulation(par)))
if (is.null(sim$Model$runs)) sim$Model$runs <- sim$Model$generations # renamed in #199
sim$packageVersion <- paste(args[1], utils::packageDescription(args[1], lib.loc = libs)$Version)
saveRDS(sim, args[3])
