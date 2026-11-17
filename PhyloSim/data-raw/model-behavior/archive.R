# Converts the cache of a run of run.R (one .rds file per simulation, RUNINFO.txt,
# SESSION-*.txt, script copies, raw/) into the archive of model behavior of one
# version, archives/<version>_<commit>/ in the main checkout: per-replicate tables
# that can be read without R, with a README.md that describes the files, and the
# raw outputs (raw/). It also writes the summary of mode main for makeData.R to
# runs/<version>_<commit>.csv: replicate mean, sd and n of each indicator at the
# last census per setting, the runtime and the share of failed simulations.
# Existing archives are not overwritten. run.R calls it at the end of a complete
# run of mode main without failed simulations.
#
# Run from the repository root:
#   Rscript PhyloSim/data-raw/model-behavior/archive.R <cache folder> [partial]
# The archive of a run of mode test is named test-<version>_<commit>.
# It stops if simulations of the design are missing in the cache, unless `partial`
# is given (e.g. for runs of some blocks).

source("PhyloSim/data-raw/model-behavior/common.R")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 1:2 || (length(args) == 2 && args[2] != "partial"))
  stop("usage: archive.R <cache folder> [partial]")
cacheDir <- args[1]
runInfo <- file.path(cacheDir, "RUNINFO.txt")
name <- if (file.exists(runInfo)) unique(sub("^archive: ", "", grep("^archive: ", readLines(runInfo), value = TRUE)))
if (length(name) != 1) stop("no unique line 'archive:' in ", runInfo, " (cache of run.R?)")
out <- dataDir("archives", name)
if (dir.exists(out)) stop(out, " exists; archives are not overwritten")
files <- list.files(cacheDir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("no .rds files in ", cacheDir)
# simulations of the design of the run (its copy of design.R) missing in the cache (not counted as failed)
source(tail(sort(list.files(cacheDir, "^design-.*\\.R$", full.names = TRUE)), 1))
modesRun <- intersect(unique(sub("_.*", "", basename(files))), names(modes))
absent <- setdiff(unlist(lapply(modesRun, function(m) basename(taskTable(m, cacheDir)$file))), basename(files))
if (length(absent)) (if (length(args) == 2) message else stop)(
  length(absent), " simulations of the design are missing in the cache (rerun run.R)")
res <- lapply(files, readRDS)
ind <- unique(unlist(lapply(res, function(r) if (is.null(r$error) && !isTRUE(r$unsupported)) r$indicators)))
if (length(ind) > 1) stop("the indicators were computed with several installations: ", paste(ind, collapse = ", "))
id <- sub("\\.rds$", "", basename(files))

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
  indicators = if (is.null(r$indicators)) NA else r$indicators,
  parameters = if (is.null(r$parameters)) NA else r$parameters))
simulations <- data.frame(simulation = simulations$simulation,
                          mode = sub("_.*", "", simulations$simulation), simulations[-1])

# indicators at each census, one column per indicator; indicators computed only at the end
# (mpdLocal, cndd) at the last census
census <- perSimulation(function(r) if (!is.null(r$census)) {
  d <- r$census
  for (n in c("mpdLocal", "cndd")) if (!is.null(r[[n]])) d[[n]] <- ifelse(d$generation == max(d$generation), r[[n]], NA)
  d
})

curves <- perSimulation(function(r) rbind(
  if (!is.null(r$sac)) data.frame(curve = "sac", x = r$sac$size, value = r$sac$sr.Mean),
  if (!is.null(r$rac)) data.frame(curve = "rac", x = seq_along(r$rac), value = r$rac),
  if (!is.null(r$decay)) data.frame(curve = "decay", x = r$decay$distance, value = r$decay$sameSpecies)))

dir.create(out, recursive = TRUE)
write.csv(simulations, file.path(out, "simulations.csv"), row.names = FALSE)
# 6 significant digits keep the archive of a release small enough to attach it to the release
for (n in c("census", "curves")) {
  d <- get(n)
  num <- vapply(d, is.double, logical(1))
  d[num] <- lapply(d[num], signif, 6)
  con <- gzfile(file.path(out, paste0(n, ".csv.gz")), "w")
  write.csv(d, con, row.names = FALSE)
  close(con)
}
# raw outputs: kept with the archive, but not attached to the release
if (dir.exists(file.path(cacheDir, "raw"))) {
  if (!file.copy(file.path(cacheDir, "raw"), out, recursive = TRUE)) warning("raw/ could not be copied")
  nRaw <- length(list.files(file.path(out, "raw"), "\\.rds$"))
  if (nRaw < sum(simulations$status == "ok")) message("raw/ has outputs of ", nRaw, " of ",
                                                       sum(simulations$status == "ok"), " simulations")
}
invisible(file.copy(list.files(cacheDir, "((RUNINFO|SESSION-.*)\\.txt|\\.R)$", full.names = TRUE), out, overwrite = TRUE))
writeLines(c(
  paste0("# Model behavior of EcoPhyloSim ", sub("_", " (commit ", name), ")"), "",
  "Per-replicate results of the simulations that document the behavior of the model (vignettes",
  "\"The EcoPhyloSim R Package\", chapter \"Emerging patterns\", and \"Model history\"), created with",
  "`PhyloSim/data-raw/model-behavior/run.R` and `archive.R` of the package repository, whose README.md",
  "describes how to regenerate it. `RUNINFO.txt` and the copies of the scripts (`*.R`) record the",
  "versions used; the indicators are defined in `indicators.R`. `SESSION-*.txt` record the command,",
  "the platform, the installations of the simulated version and of the indicators, and the versions",
  "of all packages.", "",
  "- `simulations.csv`: one row per simulation (`simulation` = `<mode>_<block>_<setting>_<rep>`;",
  "  `mode` main: the run of record, test: short check of the scripts), with its parameters (as in",
  "  `design.R`), number of generations, seed, status (`ok`, `error`, or `unsupported` if the model",
  "  version cannot run the setting), error message, runtime of the core simulation in seconds, the",
  "  simulated package and version, the commit of the EcoPhyloSim that computed the indicators, and",
  "  all parameters actually used (`Model` of the result, incl. defaults).",
  "- `census.csv.gz`: one row per census of each simulation, one column per indicator; `mpdLocal`",
  "  (phylogenetic dispersion of local communities) and `cndd` (coefficient of the number of",
  "  conspecifics within radius 2 in a logistic regression of death over the last five generations,",
  "  with species intercepts) only at the last census.",
  "- `curves.csv.gz`: species-area curve (`sac`: mean richness in subplots of `x` cells),",
  "  rank-abundance curve (`rac`: abundance at rank `x`) and distance decay (`decay`: share of",
  "  cells at distance `x` with the same species) at the end of each simulation.",
  "- `raw/` (not attached to the GitHub release): for each simulation that ran, an .rds file with",
  "  the generations (`runs`) and all maps (`Output`) of the last six censuses, and the phylogeny",
  "  at the last census; enough to recompute the indicators at the end.", "",
  "Numbers in the `.csv.gz` tables are rounded to 6 significant digits."),
  file.path(out, "README.md"))
message("archived ", nrow(simulations), " simulations in ", out)

# summary of mode main for makeData.R (data set modelHistory)
sims <- simulations[simulations$mode == "main" & simulations$status != "unsupported", ]
if (nrow(sims)) {
  last <- merge(census, sims[c("simulation", "generations")])
  last <- last[last$generation == last$generations, setdiff(names(last), c("generation", "generations"))]
  s <- stack(last[-1])
  ok <- sims$status == "ok"
  values <- rbind(data.frame(simulation = last$simulation, indicator = as.character(s$ind), value = s$values),
                  data.frame(simulation = sims$simulation[ok], indicator = "runtime", value = sims$runtime[ok]),
                  data.frame(simulation = sims$simulation, indicator = "failed", value = as.numeric(!ok)))
  params <- c("block", "setting", names(base), "generations")
  d <- merge(sims[c("simulation", params)], values)
  d <- d[!is.na(d$value), ]
  a <- aggregate(d["value"], d[c(params, "indicator")], function(v) c(mean = mean(v), sd = sd(v), n = length(v)))
  summary <- data.frame(a[c(params, "indicator")], signif(a$value, 6), indicators = ind)
  summary <- summary[do.call(order, summary[c("block", "setting", "indicator")]), ]
  write.csv(summary, file.path(mb, "runs", paste0(name, ".csv")), row.names = FALSE)
  message("summary written to ", file.path(mb, "runs", paste0(name, ".csv")))
}
