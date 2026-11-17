# Summarizes the run of one version for makeData.R (data set modelHistory,
# vignette "Model history"): replicate mean, sd and number of replicates of each
# indicator at the last census per setting, plus the pooled CNDD coefficient, the
# runtime and the share of failed simulations ("failed"), from an archive created
# with archive.R.
# Run from the repository root:
#   Rscript model-behavior/addRun.R <archive folder> [<mode>]
# e.g. Rscript model-behavior/addRun.R ../EcoPhyloSim-runs/0.3.1_852c61e
# writes model-behavior/runs/<archive folder name>.csv (mode main by default).
# Parameters that are not in the archive (older designs) get their default from
# design.R, so that the settings of all versions can be matched.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 1:2)
a <- args[1]
mode <- if (length(args) == 2) args[2] else "main"
source("model-behavior/design.R")
sims <- read.csv(file.path(a, "simulations.csv"))
sims <- sims[sims$mode == mode & sims$status != "unsupported", ]
if (nrow(sims) == 0) stop("no simulations of mode ", mode, " in ", a)
for (n in setdiff(names(base), names(sims))) sims[[n]] <- if (n == "y") sims$x else base[[n]]
census <- read.csv(file.path(a, "census.csv.gz"))
if (!"indicator" %in% names(census)) { # one column per indicator (archive.R since #181)
  s <- stack(census[setdiff(names(census), c("simulation", "generation"))])
  census <- data.frame(census[c("simulation", "generation")], indicator = as.character(s$ind), value = s$values)
}
census <- merge(census, sims[c("simulation", "generations")])
census <- census[census$generation == census$generations, c("simulation", "indicator", "value")]
cndd <- read.csv(file.path(a, "cndd.csv.gz"))
cndd <- data.frame(simulation = cndd$simulation[cndd$scope == "pooled"], indicator = "cndd",
                   value = cndd$slope[cndd$scope == "pooled"])
ok <- sims$status == "ok"
values <- rbind(census, cndd,
                data.frame(simulation = sims$simulation[ok], indicator = "runtime", value = sims$runtime[ok]),
                data.frame(simulation = sims$simulation, indicator = "failed", value = as.numeric(!ok)))
params <- c("block", "setting", names(base), "generations")
d <- merge(sims[c("simulation", params)], values)
d <- d[!is.na(d$value), ]
s <- aggregate(d["value"], d[c(params, "indicator")], function(v) c(mean = mean(v), sd = sd(v), n = length(v)))
out <- data.frame(s[c(params, "indicator")], signif(s$value, 6))
out <- out[do.call(order, out[c("block", "setting", "indicator")]), ]
write.csv(out, file.path("model-behavior", "runs", paste0(basename(a), ".csv")), row.names = FALSE)
message("written: ", basename(a), " (", nrow(out), " rows)")
