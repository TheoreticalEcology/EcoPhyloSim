# Exports the per-replicate results of one version for makeData.R (data set
# modelHistory, vignette "Model history"): the value of each indicator at the last census of
# each simulation of the main run, plus the pooled CNDD coefficient and the
# runtime, from an archive created with archive.R.
# Run from the repository root:
#   Rscript model-behavior/addRun.R <archive folder> [<mode>]
# e.g. Rscript model-behavior/addRun.R ../EcoPhyloSim-runs/0.3.1_852c61e
# writes model-behavior/runs/<archive folder name>.csv.gz (mode main by default).
# Parameters that are not in the archive (older designs) get their default from
# design.R, so that the settings of all versions can be matched.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 1:2)
a <- args[1]
mode <- if (length(args) == 2) args[2] else "main"
source("model-behavior/design.R")
sims <- read.csv(file.path(a, "simulations.csv"))
sims <- sims[sims$mode == mode & sims$status == "ok", ]
if (nrow(sims) == 0) stop("no simulations of mode ", mode, " with status ok in ", a)
for (n in setdiff(names(base), names(sims))) sims[[n]] <- if (n == "y") sims$x else base[[n]]
census <- read.csv(file.path(a, "census.csv.gz"))
census <- merge(census, sims[c("simulation", "generations")])
census <- census[census$generation == census$generations, c("simulation", "indicator", "value")]
cndd <- read.csv(file.path(a, "cndd.csv.gz"))
cndd <- data.frame(simulation = cndd$simulation[cndd$scope == "pooled"], indicator = "cndd",
                   value = cndd$slope[cndd$scope == "pooled"])
runtime <- data.frame(simulation = sims$simulation, indicator = "runtime", value = sims$runtime)
params <- setdiff(names(sims), c("mode", "status", "error", "runtime", "seed", "package", "parameters"))
out <- merge(sims[params], rbind(census, cndd, runtime))
con <- gzfile(file.path("model-behavior", "runs", paste0(basename(a), ".csv.gz")), "w")
write.csv(out, con, row.names = FALSE)
close(con)
message("written: ", basename(a))
