# Builds the data sets of the package from the runs of PhyloSim/data-raw/model-behavior/:
# - modelHistory (vignette "Model history"): the summaries in PhyloSim/data-raw/model-behavior/runs/
#   (addRun.R) of all versions and the table of changes, changes.csv;
# - packageBehavior (chapter "Emerging patterns" of the main vignette), only if
#   asked for with the archive of the run of the current version (archive.R): replicate
#   mean, sd and n of the indicators at all censuses, of the curves and of the
#   CNDD coefficient, and the runtime.
#
# Run from the repository root:
#   Rscript PhyloSim/data-raw/model-behavior/makeData.R [packageBehavior <archive of the current version>]
# The new format of packageBehavior (all censuses in `main`, no `pilot` and no `cnddSpecies`, generations
# per block instead of a "default") needs the adapted chapter "Emerging patterns"
# (#180); until then, build only modelHistory.

if (!file.exists("PhyloSim/data-raw/model-behavior/changes.csv")) stop("run from the repository root")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && (args[1] != "packageBehavior" || length(args) != 2))
  stop("usage: makeData.R [packageBehavior <archive>]")
scripts <- suppressWarnings(tryCatch(system("git describe --always --dirty", intern = TRUE), error = function(e) ""))
if (length(scripts) != 1 || !nzchar(scripts)) scripts <- "unknown"

files <- list.files("PhyloSim/data-raw/model-behavior/runs", pattern = "\\.csv$", full.names = TRUE)
files <- files[order(numeric_version(sub("_.*", "", basename(files))))]
runs <- lapply(files, function(f) {
  run <- sub("\\.csv$", "", basename(f))
  data.frame(version = sub("_.*", "", run), commit = sub(".*_", "", run), read.csv(f))
})
cols <- unique(unlist(lapply(runs, names)))
summary <- do.call(rbind, lapply(runs, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
changes <- read.csv("PhyloSim/data-raw/model-behavior/changes.csv", colClasses = "character")
changes$date <- as.Date(ifelse(changes$date == "", NA, changes$date))
modelHistory <- list(summary = summary, changes = changes[c("date", "kind", "version", "commit", "reference",
                                                             "type", "change", "outputs", "evidence")],
                     info = list(runs = basename(files), scripts = scripts, R = R.version.string))
save(modelHistory, file = "PhyloSim/data/modelHistory.rda", compress = "xz")
message("modelHistory: ", nrow(summary), " rows")

archive <- args[2]
if (!is.na(archive)) {
  sims <- read.csv(file.path(archive, "simulations.csv"))
  # replicates and generations per block of each mode, as run (not as in the current design.R)
  modes <- lapply(split(sims, sims$mode), function(d)
    list(reps = max(d$rep), generations = c(tapply(d$generations, d$block, max))))
  sims <- sims[sims$mode == "main", ]
  keys <- c("block", "setting")
  settings <- unique(sims[setdiff(names(sims), c("simulation", "mode", "rep", "seed", "file", "status",
                                                  "error", "runtime", "package", "parameters"))])
  # replicate mean, sd and n of `value` per group
  stats <- function(d, by) {
    d <- merge(sims[c("simulation", keys)], d)
    d <- d[!is.na(d$value), ]
    a <- aggregate(d["value"], d[c(keys, by)], function(v) c(mean = mean(v), sd = sd(v), n = length(v)))
    data.frame(a[c(keys, by)], a$value)
  }
  census <- read.csv(file.path(archive, "census.csv.gz")) # one column per indicator (archive.R)
  if ("indicator" %in% names(census)) stop("archive in the earlier format (one row per indicator)")
  ind <- setdiff(names(census), c("simulation", "generation"))
  s <- stack(census[ind])
  census <- data.frame(census[c("simulation", "generation")], indicator = as.character(s$ind), value = s$values)
  curves <- read.csv(file.path(archive, "curves.csv.gz"))
  cndd <- read.csv(file.path(archive, "cndd.csv.gz"))
  curve <- function(type, x) { d <- stats(curves[curves$curve == type, ], "x"); names(d)[names(d) == "x"] <- x; d }
  packageBehavior <- list(
    main = list(settings = settings[order(settings$block, settings$setting), ],
                census = stats(census, c("generation", "indicator")),
                sac = curve("sac", "size"), rac = curve("rac", "rank"), decay = curve("decay", "distance"),
                cndd = stats(data.frame(simulation = cndd$simulation, value = cndd$slope)[cndd$scope == "pooled", ], NULL),
                runtime = stats(data.frame(simulation = sims$simulation, value = sims$runtime), NULL)),
    info = list(archive = basename(normalizePath(archive)), commit = sub(".*_", "", basename(normalizePath(archive))),
                modes = modes, scripts = scripts,
                run = readLines(file.path(archive, "RUNINFO.txt")),
                runs = nrow(sims), errors = sum(sims$status == "error"),
                simulationHours = sum(sims$runtime, na.rm = TRUE) / 3600))
  save(packageBehavior, file = "PhyloSim/data/packageBehavior.rda", compress = "xz")
  message("packageBehavior: ", nrow(sims), " simulations of ", basename(normalizePath(archive)))
}
