# Attributes differences between releases to single changes of the model code.
# changes.csv gives for each tested change its parent commit, the adapters of
# both in adapters.R, the blocks for stage 2 and, if stage 2 compares another
# commit with the parent (e.g. when the change itself crashes), that commit.
# Stage 1: mode test (2 replicates, 300 generations) of the whole design for the
# change and its parent with the same seeds; census indicators identical or not.
# Stage 2: mode main (30 replicates) of the blocks that differ in stage 1; the
# last census compared with Welch tests (false discovery rate per change).
#
# Run from the repository root:
#   Rscript model-behavior/attribution.R run <stage> <cores> <work folder> [<commit> ...]
#   Rscript model-behavior/attribution.R compare <stage> <work folder>
# run installs each commit (installVersion.R, archived dependencies with
# installOldDeps.R) into <work folder>/lib-<commit> and runs run.R into
# <work folder>/stage<stage>-<commit>; finished simulations are skipped. compare
# writes <work folder>/stage<stage>.csv; enter the results in the column
# evidence of changes.csv.

args <- commandArgs(trailingOnly = TRUE)
if (!file.exists("model-behavior/changes.csv")) stop("run from the repository root")
a <- read.csv("model-behavior/changes.csv", colClasses = "character")
a <- a[a$parent != "", ]
stage <- args[2]
if (args[1] == "run" && length(args) > 4) a <- a[a$commit %in% args[-(1:4)] | a$stage2Commit %in% args[-(1:4)], ]
if (stage == "2") {
  a <- a[a$stage2 != "", ]
  a$commit <- ifelse(a$stage2Commit != "", a$stage2Commit, a$commit)
}
Rscript <- file.path(R.home("bin"), "Rscript")
cacheOf <- function(h) file.path(work, paste0("stage", stage, "-", h))
installed <- function(lib) file.exists(file.path(lib, "INSTALLINFO.txt"))

if (args[1] == "run") {
  cores <- args[3]
  work <- args[4]
  deps <- file.path(work, "lib-olddeps")
  if (!all(c("cubature", "kulife", "apTreeshape", "PhyloMeasures") %in% rownames(installed.packages(deps))))
    system2(Rscript, shQuote(c("model-behavior/installOldDeps.R", deps)))
  for (i in seq_len(nrow(a))) for (h in c(a$commit[i], a$parent[i])) {
    lib <- file.path(work, paste0("lib-", h))
    # a failed install would let every simulation fail
    if (!installed(lib) && (system2(Rscript, shQuote(c("model-behavior/installVersion.R", h, lib, deps))) != 0 ||
                            !installed(lib)))
      stop("installation of ", h, " failed")
    adapter <- if (h == a$commit[i]) a$adapter[i] else a$parentAdapter[i]
    blocks <- if (stage == "2") strsplit(a$stage2[i], " ")[[1]]
    status <- system2(Rscript, shQuote(c("model-behavior/run.R", if (stage == "1") "test" else "main", cores,
                                         adapter, cacheOf(h), paste0("lib=", lib, .Platform$path.sep, deps),
                                         blocks)))
    if (status != 0) message("run.R failed for ", h, " (exit status ", status, ")")
  }
} else if (args[1] == "compare") {
  work <- args[3]
  comm <- c("richness", "shannon", "evenness", "dominance", "conspecificNeighbors",
            "corTraitEnv", "sdCompetitionTrait", "sdNeutralTrait")
  groups <- list(communities = comm, speciesAge = c("meanSpeciesAge", "corAgeAbundance"),
                 gammaTreeHeight = c("gamma", "treeHeight"), colless = "colless")
  # indicators.R used by a run (copies in the cache); differences would be reported as model differences
  indicatorsOf <- function(dir) unique(tools::md5sum(list.files(dir, "^indicators-.*\\.R$", full.names = TRUE)))
  out <- do.call(rbind, lapply(seq_len(nrow(a)), function(i) {
    new <- cacheOf(a$commit[i])
    old <- cacheOf(a$parent[i])
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
      # per simulation: are the indicators of all censuses identical? (missing: failed or unsupported)
      d <- data.frame(task, simulations = 1, missingOld = !has[, 1], missingNew = !has[, 2], compared = ok)
      for (g in names(groups)) d[[paste0(g, "Differ")]] <- vapply(r, function(z) {
        if (is.null(z[[1]]$census) || is.null(z[[2]]$census)) return(NA)
        v <- intersect(groups[[g]], common(z))
        !isTRUE(all.equal(z[[1]]$census[v], z[[2]]$census[v], check.attributes = FALSE))
      }, logical(1))
      d <- aggregate(d[setdiff(names(d), c("block", "setting", "rep"))], d["block"], sum, na.rm = TRUE)
    } else {
      # per setting and indicator: Welch test of the last census over the replicates
      key <- paste(task$block, task$setting)
      v <- lapply(1:2, function(k) do.call(rbind, lapply(seq_along(r), function(j) {
        x <- r[[j]][[k]]$census
        if (!ok[j]) return(NULL)
        data.frame(key = key[j], x[nrow(x), common(r[[j]]), drop = FALSE])
      })))
      if (is.null(v[[1]])) return(NULL)
      d <- do.call(rbind, lapply(split(seq_along(key), key), function(j) {
        o <- v[[1]][v[[1]]$key == key[j[1]], -1, drop = FALSE]
        w <- v[[2]][v[[2]]$key == key[j[1]], -1, drop = FALSE]
        do.call(rbind, lapply(names(o), function(n) data.frame(
          block = task$block[j[1]], setting = task$setting[j[1]], indicator = n,
          missingOld = sum(!has[j, 1]), missingNew = sum(!has[j, 2]),
          n = sum(!is.na(o[[n]]) & !is.na(w[[n]])), # replicates with results in both
          mean.old = mean(o[[n]], na.rm = TRUE), mean.new = mean(w[[n]], na.rm = TRUE),
          p = tryCatch(t.test(o[[n]], w[[n]])$p.value, error = function(e) NA_real_))))
      }))
      d$p.adj <- p.adjust(d$p, "BH")
    }
    data.frame(commit = a$commit[i], parent = a$parent[i], d)
  }))
  write.csv(out, file.path(work, paste0("stage", stage, ".csv")), row.names = FALSE)
  print(out)
} else stop("first argument: run or compare")
