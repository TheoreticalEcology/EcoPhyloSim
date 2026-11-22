# Installs the archived CRAN packages that the old versions (package PhyloSim,
# v0.2 and v0.3) import into one library. Their other dependencies come from the
# default library.
#
# Run from the repository root:
#   Rscript PhyloSim/data-raw/model-behavior/installOldDeps.R <library>

lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) stop("usage: Rscript PhyloSim/data-raw/model-behavior/installOldDeps.R <library>")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))
cran <- "https://cloud.r-project.org"
install.packages("cubature", lib = lib, repos = cran) # needed by apTreeshape
archived <- c(kulife = "0.1-14", apTreeshape = "1.5-0", PhyloMeasures = "2.1")
for (p in names(archived))
  install.packages(sprintf("%s/src/contrib/Archive/%s/%s_%s.tar.gz", cran, p, p, archived[[p]]),
                   lib = lib, repos = NULL, type = "source")
missing <- setdiff(c("cubature", names(archived)), rownames(installed.packages(lib)))
if (length(missing)) stop("not installed: ", paste(missing, collapse = ", "))
