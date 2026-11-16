# Attributes differences between releases to single changes of the model code
# (attribution.csv: each change with its parent commit, or the start of a range
# of commits, the adapters of both in adapters.R, the package name and the
# blocks for stage 2).
# Stage 1: mode test (2 replicates, 300 generations) of the whole design for the
# change and its parent with the same seeds; census indicators identical or not.
# Stage 2: mode main (30 replicates) of the blocks that differ in stage 1; the
# last census compared with Welch tests (Holm correction per change).
#
# Run from the repository root:
#   Rscript model-behavior/attribution.R run <stage> <cores> <work folder> [<commit> ...]
#   Rscript model-behavior/attribution.R compare <stage> <work folder>
# run installs each commit (installVersion.R, archived dependencies with
# installOldDeps.R) into <work folder>/lib-<commit> and runs run.R into
# <work folder>/stage<stage>-<commit>; finished simulations are skipped. compare
# writes <work folder>/stage<stage>.csv.

args <- commandArgs(trailingOnly = TRUE)
if (!file.exists("model-behavior/attribution.csv")) stop("run from the repository root")
a <- read.csv("model-behavior/attribution.csv", colClasses = "character")
Rscript <- file.path(R.home("bin"), "Rscript")
cacheOf <- function(stage, h) file.path(work, paste0("stage", stage, "-", h))

if (args[1] == "run") {
  stage <- args[2]
  cores <- args[3]
  work <- args[4]
  if (length(args) > 4) a <- a[a$commit %in% args[-(1:4)], ]
  if (stage == "2") a <- a[a$stage2 != "", ]
  deps <- file.path(work, "lib-olddeps")
  if (!all(c("cubature", "kulife", "apTreeshape", "PhyloMeasures") %in% rownames(installed.packages(deps))))
    system2(Rscript, shQuote(c("model-behavior/installOldDeps.R", deps)))
  for (i in seq_len(nrow(a))) for (h in c(a$commit[i], a$parent[i])) {
    lib <- file.path(work, paste0("lib-", h))
    # a failed install would fill the cache with failed simulations, which are not rerun
    if (!file.exists(file.path(lib, a$package[i], "DESCRIPTION")) &&
        (system2(Rscript, shQuote(c("model-behavior/installVersion.R", h, lib, deps))) != 0 ||
         !file.exists(file.path(lib, a$package[i], "DESCRIPTION"))))
      stop("installation of ", h, " failed")
    adapter <- if (h == a$commit[i]) a$adapter[i] else a$parentAdapter[i]
    blocks <- if (stage == "2") strsplit(a$stage2[i], " ")[[1]]
    status <- system2(Rscript, shQuote(c("model-behavior/run.R", if (stage == "1") "test" else "main", cores,
                                         adapter, cacheOf(stage, h), paste0("lib=", lib, .Platform$path.sep, deps),
                                         paste0("commit=", h), paste0("package=", a$package[i]), blocks)))
    if (status != 0) message("run.R failed for ", h, " (exit status ", status, ")")
  }
} else if (args[1] == "compare") {
  stage <- args[2]
  work <- args[3]
  comm <- c("richness", "shannon", "evenness", "dominance", "conspecificNeighbors",
            "corTraitEnv", "sdCompetitionTrait", "sdNeutralTrait")
  groups <- list(communities = comm, speciesAge = c("meanSpeciesAge", "corAgeAbundance"),
                 gammaTreeHeight = c("gamma", "treeHeight"), colless = "colless")
  # indicators.R used by a run (copies in the cache); differences would be reported as model differences
  indicatorsOf <- function(dir) unique(tools::md5sum(list.files(dir, "^indicators-.*\\.R$", full.names = TRUE)))
  out <- do.call(rbind, lapply(seq_len(nrow(a)), function(i) {
    new <- cacheOf(stage, a$commit[i])
    old <- cacheOf(stage, a$parent[i])
    f <- intersect(list.files(old, "rds$"), list.files(new, "rds$"))
    if (!length(f)) return(NULL)
    if (!setequal(indicatorsOf(old), indicatorsOf(new)))
      warning("the runs of ", a$commit[i], " and ", a$parent[i], " used different versions of indicators.R")
    r <- lapply(f, function(x) list(readRDS(file.path(old, x)), readRDS(file.path(new, x))))
    task <- do.call(rbind, lapply(r, function(z) z[[1]]$task[c("block", "setting", "rep")]))
    has <- sapply(1:2, function(k) vapply(r, function(z) !is.null(z[[k]]$census), logical(1)))
    ok <- has[, 1] & has[, 2]
    common <- function(z) setdiff(intersect(names(z[[1]]$census), names(z[[2]]$census)), "generation")
    if (stage == "1") {
      # per simulation: are the indicators of all censuses identical?
      d <- data.frame(task, simulations = 1, failedOld = !has[, 1], failedNew = !has[, 2], compared = ok)
      for (g in names(groups)) d[[paste0(g, "Differ")]] <- vapply(r, function(z) {
        if (is.null(z[[1]]$census) || is.null(z[[2]]$census)) return(NA)
        v <- intersect(groups[[g]], common(z))
        !isTRUE(all.equal(z[[1]]$census[v], z[[2]]$census[v], check.attributes = FALSE))
      }, logical(1))
      d <- aggregate(d[setdiff(names(d), c("block", "setting", "rep"))], d["block"], sum, na.rm = TRUE)
    } else {
      if (!any(ok)) return(NULL)
      # per setting and indicator: Welch test of the last census over the replicates
      v <- lapply(1:2, function(k) do.call(rbind, lapply(r[ok], function(z) {
        x <- z[[k]]$census
        x[nrow(x), common(z), drop = FALSE]
      })))
      key <- paste(task$block[ok], task$setting[ok])
      d <- do.call(rbind, lapply(split(seq_along(key), key), function(j) do.call(rbind, lapply(names(v[[1]]), function(n)
        data.frame(block = task$block[ok][j[1]], setting = task$setting[ok][j[1]], indicator = n, n = length(j),
                   mean.old = mean(v[[1]][j, n], na.rm = TRUE), mean.new = mean(v[[2]][j, n], na.rm = TRUE),
                   p = tryCatch(t.test(v[[1]][j, n], v[[2]][j, n])$p.value, error = function(e) NA_real_))))))
      d$p.adj <- p.adjust(d$p, "holm")
    }
    data.frame(commit = a$commit[i], parent = a$parent[i], d)
  }))
  write.csv(out, file.path(work, paste0("stage", stage, ".csv")), row.names = FALSE)
  print(out)
} else stop("first argument: run or compare")
