# Converts the per-simulation results of a run (the cache folder of run.R, or of
# the earlier PhyloSim/data-raw/packageBehavior.R: all .rds files in the folder,
# plus RUNINFO.txt and the script copies) into the archive of model behavior for
# one version: per-replicate tables that can be read without R, with a README.md
# that describes the files.
#
# Run from the repository root:
#   Rscript model-behavior/archive.R <cache folder> <output folder> <version> <commit>
# e.g. Rscript model-behavior/archive.R ../cache-current ../EcoPhyloSim-runs 0.3.1.9000 7385ca7
# writes <output folder>/<version>_<commit>/. <version> must be a version number
# (e.g. 0.3.1.9000 for a development version), as makeData.R sorts by it.

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
  r$task[names(r$task) != "file"],
  status = if (!is.null(r$error)) "error" else if (isTRUE(r$unsupported)) "unsupported" else "ok",
  error = if (is.null(r$error)) NA else r$error,
  runtime = if (is.null(r$runtime)) NA else r$runtime,
  package = if (is.null(r$package)) NA else r$package,
  parameters = if (is.null(r$parameters)) NA else r$parameters))
simulations <- data.frame(simulation = simulations$simulation,
                          mode = sub("_.*", "", simulations$simulation), simulations[-1])

# indicators at each census; indicators computed only at the end (mpdLocal) at the last census
census <- perSimulation(function(r) if (!is.null(r$census)) {
  ind <- setdiff(names(r$census), "generation")
  rbind(data.frame(generation = rep(r$census$generation, length(ind)),
                   indicator = rep(ind, each = nrow(r$census)),
                   value = as.numeric(unlist(r$census[ind]))),
        if (!is.null(r$mpdLocal)) data.frame(generation = max(r$census$generation),
                                              indicator = "mpdLocal", value = r$mpdLocal))
})

curves <- perSimulation(function(r) rbind(
  if (!is.null(r$sac)) data.frame(curve = "sac", x = r$sac$size, value = r$sac$sr.Mean),
  if (!is.null(r$rac)) data.frame(curve = "rac", x = seq_along(r$rac), value = r$rac),
  if (!is.null(r$decay)) data.frame(curve = "decay", x = r$decay$distance, value = r$decay$sameSpecies)))

cndd <- perSimulation(function(r) rbind(
  if (!is.null(r$cndd)) data.frame(scope = "pooled", abundance = NA, slope = r$cndd),
  if (NROW(r$cnddSpecies) > 0) data.frame(scope = "species", abundance = r$cnddSpecies$abundance,
                                          slope = r$cnddSpecies$slope)))

maps <- perSimulation(function(r) if (!is.null(r$specMat))
  data.frame(row = c(row(r$specMat)), col = c(col(r$specMat)), species = c(r$specMat)))

write.csv(simulations, file.path(out, "simulations.csv"), row.names = FALSE)
for (n in c("census", "curves", "cndd", "maps")) {
  con <- gzfile(file.path(out, paste0(n, ".csv.gz")), "w")
  write.csv(get(n), con, row.names = FALSE)
  close(con)
}
file.copy(list.files(cacheDir, "((RUNINFO|SESSION-.*)\\.txt|\\.R)$", full.names = TRUE), out, overwrite = TRUE)
writeLines(c(
  paste0("# Model behavior of EcoPhyloSim ", args[3], " (commit ", args[4], ")"), "",
  "Per-replicate results of the simulations that document the behavior of the model (vignette",
  "chapter \"Emerging patterns\"; `model-behavior/run.R`, earlier `PhyloSim/data-raw/packageBehavior.R`),",
  "created with `model-behavior/archive.R`. `RUNINFO.txt` and the copies of the scripts (`*.R`)",
  "record the versions used and notes on the run; the indicators are defined in `indicators.R`",
  "(earlier in `censusStats()` of `packageBehavior.R`). `SESSION-*.txt` record the command, the",
  "platform, the installation of the simulated version and the versions of all packages. How to",
  "regenerate the archive is described in `model-behavior/README.md` of the repository.", "",
  "- `simulations.csv`: one row per simulation (`simulation` = `<mode>_<block>_<setting>_<rep>`;",
  "  `mode` main: the run of record, test: short check of the scripts; older archives also have",
  "  pilot: few replicates to check the time to equilibrium), main: more",
  "  replicates, test: short check of the scripts), with its parameters (as in `design.R`; `y` equals `x` if absent), number of",
  "  generations, seed, status (`ok`, `error`, or `unsupported` if the model version cannot run",
  "  the setting), error message, runtime of the core simulation in seconds, the simulated",
  "  package and version, and all parameters actually used (`Model` of the result, incl. defaults).",
  "- `census.csv.gz`: value of each indicator at each census of each simulation (`mpdLocal`, the",
  "  phylogenetic dispersion of local communities, only at the last census).",
  "- `maps.csv.gz`: species of each cell (`row`, `col`) at the end of the first replicate of each setting.",
  "- `curves.csv.gz`: species-area curve (`sac`: mean richness in subplots of `x` cells),",
  "  rank-abundance curve (`rac`: abundance at rank `x`) and distance decay (`decay`: share of",
  "  cells at distance `x` with the same species) at the end of each simulation.",
  "- `cndd.csv.gz`: coefficient (`slope`) of the number of conspecifics within radius 2 in a",
  "  logistic regression of death over the last five generations: pooled with species intercepts",
  "  (`abundance` NA) and per species with at least 50 observations (`abundance`: mean number of",
  "  individuals per census).", "",
  "Trait maps and phylogenies are not stored."),
  file.path(out, "README.md"))
message("archived ", nrow(simulations), " simulations in ", out)
