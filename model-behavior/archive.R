# Converts the per-simulation results of a packageBehavior run (the cache of
# PhyloSim/data-raw/packageBehavior.R: all .rds files in the folder, plus
# RUNINFO.txt and the script copies) into the archive of model behavior for one
# version: per-replicate tables that can be read without R. The archive gets a
# README.md that describes the files. For archived runs, extract cache.tar.gz first.
#
# Run from the repository root:
#   Rscript model-behavior/archive.R <cache folder> <output folder> <label> <commit>
# e.g. Rscript model-behavior/archive.R PhyloSim/data-raw/cache ../EcoPhyloSim-runs 0.3.1 852c61e
# writes <output folder>/<label>_<commit>/.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
cacheDir <- args[1]
out <- file.path(args[2], paste0(args[3], "_", args[4]))
files <- list.files(cacheDir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("no .rds files in ", cacheDir)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
res <- lapply(files, readRDS)
id <- sub("\\.rds$", "", basename(files))

# rbind of data frames with different columns (e.g. parameters that a version lacks)
bind <- function(l) {
  l <- Filter(Negate(is.null), l)
  if (length(l) == 0) return(data.frame(simulation = character(0)))
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}
perSimulation <- function(f) bind(lapply(seq_along(res), function(i) {
  d <- f(res[[i]])
  if (is.null(d) || nrow(d) == 0) NULL else data.frame(simulation = id[i], d, row.names = NULL)
}))

simulations <- perSimulation(function(r) data.frame(
  mode = NA, r$task[names(r$task) != "file"],
  status = if (is.null(r$error)) "ok" else "error",
  error = if (is.null(r$error)) NA else r$error,
  runtime = if (is.null(r$runtime)) NA else r$runtime))
simulations$mode <- sub("_.*", "", simulations$simulation)

census <- perSimulation(function(r) if (!is.null(r$census)) {
  ind <- setdiff(names(r$census), "generation")
  data.frame(generation = rep(r$census$generation, length(ind)),
             indicator = rep(ind, each = nrow(r$census)),
             value = as.numeric(unlist(r$census[ind])))
})

curves <- perSimulation(function(r) rbind(
  if (!is.null(r$sac)) data.frame(curve = "sac", x = r$sac$size, value = r$sac$sr.Mean),
  if (!is.null(r$rac)) data.frame(curve = "rac", x = seq_along(r$rac), value = r$rac),
  if (!is.null(r$decay)) data.frame(curve = "decay", x = r$decay$distance, value = r$decay$sameSpecies)))

cndd <- perSimulation(function(r) rbind(
  if (!is.null(r$cndd)) data.frame(scope = "pooled", abundance = NA, slope = r$cndd),
  if (NROW(r$cnddSpecies) > 0) data.frame(scope = "species", abundance = r$cnddSpecies$abundance,
                                          slope = r$cnddSpecies$slope)))

write.csv(simulations, file.path(out, "simulations.csv"), row.names = FALSE)
for (n in c("census", "curves", "cndd")) {
  con <- gzfile(file.path(out, paste0(n, ".csv.gz")), "w")
  write.csv(get(n), con, row.names = FALSE)
  close(con)
}
file.copy(list.files(cacheDir, "(RUNINFO\\.txt|\\.R)$", full.names = TRUE), out, overwrite = TRUE)
writeLines(c(
  paste0("# Model behavior of EcoPhyloSim ", args[3], " (commit ", args[4], ")"), "",
  "Per-replicate results of the simulations of `PhyloSim/data-raw/packageBehavior.R` (vignette",
  "chapter \"Emerging patterns\"), created with `model-behavior/archive.R`. `RUNINFO.txt` and the",
  "copies of the run script (`*.R`) record the package version, the script versions used and",
  "notes on the run; the indicators are defined in `censusStats()` of the script.", "",
  "- `simulations.csv`: one row per simulation (`simulation` = `<mode>_<block>_<setting>_<rep>`;",
  "  `mode` pilot: few replicates and many censuses to check the time to equilibrium, main: more",
  "  replicates), with its parameters (arguments of `createCompletePar()`, `y` equals `x`), number",
  "  of generations, seed, status, error message and runtime of the core simulation in seconds.",
  "- `census.csv.gz`: value of each indicator at each census of each simulation.",
  "- `curves.csv.gz`: species-area curve (`sac`: mean richness in subplots of `x` cells),",
  "  rank-abundance curve (`rac`: abundance at rank `x`) and distance decay (`decay`: share of",
  "  cells at distance `x` with the same species) at the end of each simulation.",
  "- `cndd.csv.gz`: coefficient (`slope`) of the number of conspecifics within radius 2 in a",
  "  logistic regression of death over the last five generations: pooled with species intercepts",
  "  (`abundance` NA) and per species with at least 50 observations (`abundance`: mean number of",
  "  individuals per census).", "",
  "Species maps, trait maps and phylogenies are not stored."),
  file.path(out, "README.md"))
message("archived ", nrow(simulations), " simulations in ", out)
