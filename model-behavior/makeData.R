# Builds the data set modelHistory of the package (vignette "Model history")
# from the runs in model-behavior/runs/ and changes.csv:
# replicate mean, sd and n of each indicator at the last census per version and
# setting, and the table of changes.
#
# Run from the repository root:
#   Rscript model-behavior/makeData.R

if (!file.exists("model-behavior/changes.csv")) stop("run from the repository root")
files <- list.files("model-behavior/runs", pattern = "csv.gz$", full.names = TRUE)
files <- files[order(numeric_version(sub("_.*", "", basename(files))))]
summaries <- lapply(files, function(f) {
  d <- read.csv(f)
  d <- d[!is.na(d$value), ]
  keys <- setdiff(names(d), c("simulation", "rep", "value"))
  # missing parameter values as "NA", so that aggregate() keeps their settings
  d[keys] <- lapply(d[keys], function(x) if (anyNA(x)) ifelse(is.na(x), "NA", x) else x)
  a <- aggregate(d["value"], d[keys], function(v) c(mean = mean(v), sd = sd(v), n = length(v)))
  run <- sub(".csv.gz", "", basename(f), fixed = TRUE)
  data.frame(version = sub("_.*", "", run), commit = sub(".*_", "", run), a[keys],
             mean = a$value[, "mean"], sd = a$value[, "sd"], n = a$value[, "n"])
})
cols <- unique(unlist(lapply(summaries, names)))
summary <- do.call(rbind, lapply(summaries, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
changes <- read.csv("model-behavior/changes.csv", colClasses = "character")
changes$date <- as.Date(ifelse(changes$date == "", NA, changes$date))
modelHistory <- list(summary = summary, changes = changes,
                      info = list(created = Sys.time(), runs = basename(files),
                                  scripts = system("git describe --always --dirty", intern = TRUE),
                                  R = R.version.string))
save(modelHistory, file = "PhyloSim/data/modelHistory.rda", compress = "xz")
message("saved ", nrow(summary), " rows, ", round(file.size("PhyloSim/data/modelHistory.rda") / 1024), " KB")
