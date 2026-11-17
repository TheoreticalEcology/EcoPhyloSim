# Builds the data sets of the package from the runs of PhyloSim/data-raw/model-behavior/:
# - modelHistory (vignette "Model history"): the summaries in runs/ (archive.R) of
#   all versions and the table of changes, changes.csv; all runs must have computed
#   the indicators with `indicatorCommit` of design.R, one run per version;
# - packageBehavior (chapter "Emerging patterns" of the main vignette), only if
#   asked for with the archive of the run of the current version (a folder name in
#   archives/ or a path): replicate mean, sd and n of the indicators at all
#   censuses, of the curves and of the CNDD coefficient, and the runtime.
#
# Run from the repository root:
#   Rscript PhyloSim/data-raw/model-behavior/makeData.R [packageBehavior <archive>]
# The new format of packageBehavior (all censuses in `main`, no `pilot` and `cnddSpecies`, generations
# per block instead of a "default") needs the adapted chapter "Emerging patterns"
# (#180); until then, build only modelHistory.

source("PhyloSim/data-raw/model-behavior/common.R")
source(file.path(mb, "design.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && (args[1] != "packageBehavior" || length(args) != 2))
  stop("usage: makeData.R [packageBehavior <archive>]")
scripts <- gitDescribe()

files <- list.files(file.path(mb, "runs"), pattern = "\\.csv$", full.names = TRUE)
files <- files[order(numeric_version(sub("_.*", "", basename(files))))]
v <- sub("_.*", "", basename(files))
if (anyDuplicated(v)) stop("several runs of version ", paste(unique(v[duplicated(v)]), collapse = ", "),
                           " in runs/; keep only the run of record")
summary <- bind(lapply(files, function(f) {
  run <- sub("\\.csv$", "", basename(f))
  data.frame(version = sub("_.*", "", run), commit = sub(".*_", "", run), read.csv(f))
}))
# all versions are compared with the same indicators (run.R)
if (is.null(summary$indicators)) summary$indicators <- NA
other <- unique(summary$version[is.na(summary$indicators) | !startsWith(summary$indicators, indicatorCommit)])
if (length(other)) stop("the runs of version ", paste(other, collapse = ", "), " did not compute the indicators with ",
                        indicatorCommit, " (indicatorCommit of design.R); run them again")
changes <- read.csv(file.path(mb, "changes.csv"), colClasses = "character")
changes$date <- as.Date(ifelse(changes$date == "", NA, changes$date))
modelHistory <- list(summary = summary[names(summary) != "indicators"],
                     changes = changes[c("date", "kind", "version", "commit", "reference",
                                         "type", "change", "outputs", "evidence")],
                     info = list(runs = basename(files), indicators = unique(summary$indicators),
                                 scripts = scripts, R = R.version.string))
save(modelHistory, file = "PhyloSim/data/modelHistory.rda", compress = "xz")
message("modelHistory: ", nrow(summary), " rows")

archive <- args[2]
if (!is.na(archive)) {
  if (!dir.exists(archive)) archive <- dataDir("archives", archive)
  sims <- read.csv(file.path(archive, "simulations.csv"))
  # replicates and generations per block of each mode, as run (not as in the current design.R)
  modes <- lapply(split(sims, sims$mode), function(d)
    list(reps = max(d$rep), generations = c(tapply(d$generations, d$block, max))))
  sims <- sims[sims$mode == "main", ]
  keys <- c("block", "setting")
  settings <- unique(sims[setdiff(names(sims), c("simulation", "mode", "rep", "seed", "file", "status",
                                                  "error", "runtime", "package", "indicators", "parameters"))])
  # replicate mean, sd and n of `value` per group
  stats <- function(d, by) {
    d <- merge(sims[c("simulation", keys)], d)
    d <- d[!is.na(d$value), ]
    a <- aggregate(d["value"], d[c(keys, by)], function(v) c(mean = mean(v), sd = sd(v), n = length(v)))
    data.frame(a[c(keys, by)], a$value)
  }
  census <- read.csv(file.path(archive, "census.csv.gz")) # one column per indicator (archive.R)
  if (!"cndd" %in% names(census)) stop("archive in an earlier format (archive.R before #180)")
  ind <- setdiff(names(census), c("simulation", "generation", "cndd"))
  s <- stack(census[ind])
  curves <- read.csv(file.path(archive, "curves.csv.gz"))
  curve <- function(type, x) { d <- stats(curves[curves$curve == type, ], "x"); names(d)[names(d) == "x"] <- x; d }
  packageBehavior <- list(
    main = list(settings = settings[order(settings$block, settings$setting), ],
                census = stats(data.frame(census[c("simulation", "generation")], indicator = as.character(s$ind),
                                          value = s$values), c("generation", "indicator")),
                sac = curve("sac", "size"), rac = curve("rac", "rank"), decay = curve("decay", "distance"),
                cndd = stats(data.frame(simulation = census$simulation, value = census$cndd), NULL),
                runtime = stats(data.frame(simulation = sims$simulation, value = sims$runtime), NULL)),
    info = list(archive = basename(normalizePath(archive)), commit = sub(".*_", "", basename(normalizePath(archive))),
                modes = modes, scripts = scripts,
                run = readLines(file.path(archive, "RUNINFO.txt")),
                runs = nrow(sims), errors = sum(sims$status == "error"),
                simulationHours = sum(sims$runtime, na.rm = TRUE) / 3600))
  save(packageBehavior, file = "PhyloSim/data/packageBehavior.rda", compress = "xz")
  message("packageBehavior: ", nrow(sims), " simulations of ", basename(normalizePath(archive)))
}
