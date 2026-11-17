# Installs the package of one commit or tag of this repository into its own
# library, with the build fixes that old versions need (the model code is not
# changed):
# - src/RcppExports.cpp is missing (v0.2, v0.3; #102): generated with
#   Rcpp::compileAttributes();
# - v0.2: the lines with callLeipzig, a function without C++ code, are then
#   removed from the generated src/RcppExports.cpp.
# Writes <library>/INSTALLINFO.txt with the commit, the fixes and the R session.
#
# Run from the repository root:
#   Rscript PhyloSim/data-raw/model-behavior/installVersion.R <commit> <library> [<library of old dependencies>]
# The old versions (package PhyloSim) need the library of installOldDeps.R.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript PhyloSim/data-raw/model-behavior/installVersion.R <commit> <library> [<library>]")
commit <- system2("git", c("rev-parse", "--verify", paste0(args[1], "^{commit}")), stdout = TRUE)
lib <- args[2]
if (length(args) > 2) .libPaths(c(args[3], .libPaths()))
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
src <- tempfile("src")
dir.create(src)
tarFile <- file.path(src, "src.tar")
if (system2("git", c("archive", "-o", tarFile, commit, "PhyloSim")) != 0) stop("git archive failed")
untar(tarFile, exdir = src)
pkg <- file.path(src, "PhyloSim")
name <- read.dcf(file.path(pkg, "DESCRIPTION"), "Package")[1, 1]
fixes <- character(0)
rcpp <- file.path(pkg, "src", "RcppExports.cpp")
if (!file.exists(rcpp)) {
  Rcpp::compileAttributes(pkg)
  fixes <- "src/RcppExports.cpp generated with Rcpp::compileAttributes() (#102)"
  l <- readLines(rcpp)
  if (any(grepl("callLeipzig", l))) {
    writeLines(grep("callLeipzig", l, invert = TRUE, value = TRUE), rcpp)
    fixes <- c(fixes, "lines with callLeipzig removed from src/RcppExports.cpp")
  }
}
install.packages(pkg, lib = lib, repos = NULL, type = "source")
if (!file.exists(file.path(lib, name, "DESCRIPTION"))) stop(name, " was not installed")
writeLines(c(paste("commit:", commit), paste("package:", name, read.dcf(file.path(pkg, "DESCRIPTION"), "Version")[1, 1]),
             paste("installed:", format(Sys.time())),
             paste("build fixes:", if (length(fixes)) paste(fixes, collapse = "; ") else "none"),
             "", capture.output(sessionInfo())),
           file.path(lib, "INSTALLINFO.txt"))
unlink(src, recursive = TRUE)
