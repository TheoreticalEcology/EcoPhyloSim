# Builds the data sets of the package from the runs of model-behavior/:
# - modelHistory (vignette "Model history"): the summaries in model-behavior/runs/
#   (addRun.R) of all versions and the table of changes, changes.csv;
# - packageBehavior (chapter "Emerging patterns" of the main vignette), if the
#   archive of the run of the current version (archive.R) is given: replicate
#   mean, sd and n of the indicators at all censuses, of the curves and of the
#   CNDD coefficient, the CNDD coefficients of single species, and the runtime.
#
# Run from the repository root:
#   Rscript model-behavior/makeData.R [<archive of the current version>]

if (!file.exists("model-behavior/changes.csv")) stop("run from the repository root")
scripts <- tryCatch(system("git describe --always --dirty", intern = TRUE), error = function(e) "unknown")

files <- list.files("model-behavior/runs", pattern = "\\.csv$", full.names = TRUE)
files <- files[order(numeric_version(sub("_.*", "", basename(files))))]
runs <- lapply(files, function(f) {
  run <- sub("\\.csv$", "", basename(f))
  data.frame(version = sub("_.*", "", run), commit = sub(".*_", "", run), read.csv(f))
})
cols <- unique(unlist(lapply(runs, names)))
summary <- do.call(rbind, lapply(runs, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
changes <- read.csv("model-behavior/changes.csv", colClasses = "character")
changes$date <- as.Date(ifelse(changes$date == "", NA, changes$date))
modelHistory <- list(summary = summary, changes = changes[c("date", "kind", "version", "commit", "reference",
                                                             "type", "change", "outputs", "evidence")],
                     info = list(runs = basename(files), scripts = scripts, R = R.version.string))
save(modelHistory, file = "PhyloSim/data/modelHistory.rda", compress = "xz")
message("modelHistory: ", nrow(summary), " rows")

archive <- commandArgs(trailingOnly = TRUE)[1]
if (!is.na(archive)) {
  sims <- read.csv(file.path(archive, "simulations.csv"))
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
  census <- read.csv(file.path(archive, "census.csv.gz"))
  curves <- read.csv(file.path(archive, "curves.csv.gz"))
  cndd <- read.csv(file.path(archive, "cndd.csv.gz"))
  curve <- function(type, x) { d <- stats(curves[curves$curve == type, ], "x"); names(d)[names(d) == "x"] <- x; d }
  species <- merge(sims[c("simulation", keys, "rep")], cndd[cndd$scope == "species", c("simulation", "abundance", "slope")])
  packageBehavior <- list(
    main = list(settings = settings[order(settings$block, settings$setting), ],
                census = stats(census, c("generation", "indicator")),
                sac = curve("sac", "size"), rac = curve("rac", "rank"), decay = curve("decay", "distance"),
                cndd = stats(data.frame(simulation = cndd$simulation, value = cndd$slope)[cndd$scope == "pooled", ], NULL),
                cnddSpecies = species[c(keys, "rep", "abundance", "slope")],
                runtime = stats(data.frame(simulation = sims$simulation, value = sims$runtime), NULL)),
    info = list(archive = basename(normalizePath(archive)), commit = sub(".*_", "", basename(normalizePath(archive))),
                modes = local({ source("model-behavior/design.R", local = TRUE); modes }), scripts = scripts,
                run = readLines(file.path(archive, "RUNINFO.txt")),
                runs = nrow(sims), errors = sum(sims$status == "error"),
                simulationHours = sum(sims$runtime, na.rm = TRUE) / 3600))
  save(packageBehavior, file = "PhyloSim/data/packageBehavior.rda", compress = "xz")
  message("packageBehavior: ", nrow(sims), " simulations of ", basename(normalizePath(archive)))
}
