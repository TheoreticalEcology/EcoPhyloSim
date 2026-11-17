# Helpers of the scripts in this folder (sourced from the repository root).

mb <- file.path("PhyloSim", "data-raw", "model-behavior")
if (!file.exists(file.path(mb, "design.R"))) stop("run the scripts from the repository root")
Rscript <- file.path(R.home("bin"), "Rscript")

# git describe of the checkout, or "unknown" outside a git checkout
gitDescribe <- function() {
  d <- suppressWarnings(tryCatch(system("git describe --always --dirty", intern = TRUE), error = function(e) ""))
  if (length(d) == 1 && nzchar(d)) d else "unknown"
}

# the folders not in git (archives/, cache/, lib/) are in the main checkout, also
# when a script runs in a git worktree, which can be removed together with its files
dataDir <- function(...) {
  g <- suppressWarnings(tryCatch(system("git rev-parse --path-format=absolute --git-common-dir", intern = TRUE),
                                 error = function(e) character(0)))
  file.path(if (length(g) == 1) dirname(g) else ".", mb, ...)
}

# cache of a run of run.R
cacheOf <- function(mode, version, commit) dataDir("cache", paste(mode, version, substr(commit, 1, 7), sep = "_"))

# library with the package of a commit or tag (installed with installVersion.R into
# lib/<commit> if missing), followed by the library of the archived dependencies
# of the old versions (installOldDeps.R) if `oldDeps`; separated by ";" or ":"
libraryOf <- function(ref, oldDeps = FALSE) {
  h <- suppressWarnings(system2("git", c("rev-parse", "--short=7", paste0(ref, "^{commit}")), stdout = TRUE))
  if (length(h) != 1) stop("unknown commit or tag ", ref)
  deps <- if (oldDeps) dataDir("lib", "olddeps")
  if (oldDeps && !all(c("cubature", "kulife", "apTreeshape", "PhyloMeasures") %in% rownames(installed.packages(deps))))
    system2(Rscript, shQuote(c(file.path(mb, "installOldDeps.R"), deps)))
  lib <- dataDir("lib", h)
  if (!file.exists(file.path(lib, "INSTALLINFO.txt")) &&
      (system2(Rscript, shQuote(c(file.path(mb, "installVersion.R"), ref, lib, deps))) != 0 ||
       !file.exists(file.path(lib, "INSTALLINFO.txt"))))
    stop("installation of ", ref, " failed")
  paste(c(lib, deps), collapse = .Platform$path.sep)
}

# a field of INSTALLINFO.txt of a library (installVersion.R)
installInfo <- function(lib, field) {
  l <- readLines(file.path(strsplit(lib, .Platform$path.sep)[[1]][1], "INSTALLINFO.txt"))
  sub(paste0("^", field, ": "), "", grep(paste0("^", field, ": "), l, value = TRUE))
}

# rbind of data frames with different columns (e.g. parameters that a version lacks)
bind <- function(l) {
  l <- Filter(Negate(is.null), l)
  if (length(l) == 0) return(data.frame(simulation = character(0)))
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}
