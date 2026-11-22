# Attributes differences between releases to single changes of the model code.
# changes.csv gives for each tested change its parent commit, the adapters of
# both in adapters.R, the blocks for stage 2 and, if stage 2 compares another
# commit with the parent (e.g. when the change itself crashes), that commit.
# Stage 1: mode test (2 replicates, 300 generations) of the whole design for the
# change and its parent with the same seeds; census indicators identical or not.
# Stage 2: mode main (30 replicates) of the blocks that differ in stage 1; the
# last census compared with paired t-tests over the replicates, which have the same
# seeds in both runs (false discovery rate per change).
#
# Run from the repository root:
#   Rscript PhyloSim/data-raw/model-behavior/attribution.R run <stage> <cores> [<commit> ...]
#   Rscript PhyloSim/data-raw/model-behavior/attribution.R compare <stage>
# run runs run.R for each change and its parent (installed into lib/ if missing;
# caches in cache/, see README.md); finished simulations are skipped. compare
# writes cache/attribution-stage<stage>.csv; enter the results in the column
# evidence of changes.csv.

source("PhyloSim/data-raw/model-behavior/common.R")
args <- commandArgs(trailingOnly = TRUE)
a <- read.csv(file.path(mb, "changes.csv"), colClasses = "character")
a <- a[a$parent != "", ]
stage <- args[2]
mode <- if (stage == "1") "test" else "main"
if (args[1] == "run" && length(args) > 3) a <- a[a$commit %in% args[-(1:3)] | a$stage2Commit %in% args[-(1:3)], ]
if (stage == "2") {
  a <- a[a$stage2 != "", ]
  a$commit <- ifelse(a$stage2Commit != "", a$stage2Commit, a$commit)
}

if (args[1] == "run") {
  for (i in seq_len(nrow(a))) for (h in c(a$commit[i], a$parent[i])) {
    adapter <- if (h == a$commit[i]) a$adapter[i] else a$parentAdapter[i]
    blocks <- if (stage == "2") strsplit(a$stage2[i], " ")[[1]]
    status <- system2(Rscript, shQuote(c(file.path(mb, "run.R"), mode, args[3], adapter, h, blocks)))
    if (status != 0) message("run.R failed for ", h, " (exit status ", status, ")")
  }
} else if (args[1] == "compare") {
  comm <- c("richness", "shannon", "evenness", "dominance", "conspecificNeighbors",
            "corTraitEnv", "sdCompetitionTrait", "sdNeutralTrait")
  groups <- list(communities = comm, speciesAge = c("meanSpeciesAge", "corAgeAbundance"),
                 gammaTreeHeight = c("gamma", "treeHeight"), colless = "colless")
  # indicators.R used by a run (copies in the cache); differences would be reported as model differences
  indicatorsOf <- function(dir) unique(tools::md5sum(list.files(dir, "^indicators-.*\\.R$", full.names = TRUE)))
  out <- do.call(rbind, lapply(seq_len(nrow(a)), function(i) {
    new <- cacheOf(mode, a$adapter[i], a$commit[i])
    old <- cacheOf(mode, a$parentAdapter[i], a$parent[i])
    f <- intersect(list.files(old, "rds$"), list.files(new, "rds$"))
    if (!length(f)) { message("no common simulations in ", old, " and ", new); return(NULL) }
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
      # per setting and indicator: paired t-test of the last census over the replicates
      key <- paste(task$block, task$setting)
      v <- lapply(1:2, function(k) do.call(rbind, lapply(seq_along(r), function(j) {
        x <- r[[j]][[k]]$census
        if (!ok[j]) return(NULL)
        data.frame(key = key[j], x[nrow(x), common(r[[j]]), drop = FALSE])
      })))
      if (is.null(v[[1]])) return(NULL)
      d <- do.call(rbind, lapply(split(seq_along(key), key), function(j) {
        # rows of both runs in the same order of replicates (same seeds)
        o <- v[[1]][v[[1]]$key == key[j[1]], -1, drop = FALSE]
        w <- v[[2]][v[[2]]$key == key[j[1]], -1, drop = FALSE]
        do.call(rbind, lapply(names(o), function(n) data.frame(
          block = task$block[j[1]], setting = task$setting[j[1]], indicator = n,
          missingOld = sum(!has[j, 1]), missingNew = sum(!has[j, 2]),
          n = sum(!is.na(o[[n]]) & !is.na(w[[n]])), # replicates with results in both
          mean.old = mean(o[[n]], na.rm = TRUE), mean.new = mean(w[[n]], na.rm = TRUE),
          p = tryCatch(t.test(w[[n]], o[[n]], paired = TRUE)$p.value, error = function(e) NA_real_))))
      }))
      d$p.adj <- p.adjust(d$p, "BH")
    }
    data.frame(commit = a$commit[i], parent = a$parent[i], d)
  }))
  write.csv(out, dataDir("cache", paste0("attribution-stage", stage, ".csv")), row.names = FALSE)
  print(out)
} else stop("first argument: run or compare")
