# Exports the per-replicate results of one version for model-behavior.qmd, which
# compares the model versions: the value of each indicator at the last census of
# each simulation of the main run, plus the pooled CNDD coefficient and the
# runtime, from an archive created with archive.R.
# Run from the repository root:
#   Rscript model-behavior/addRun.R <archive folder>
# e.g. Rscript model-behavior/addRun.R ../EcoPhyloSim-runs/0.3.1_852c61e
# writes model-behavior/runs/<archive folder name>.csv.gz.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
a <- args[1]
sims <- read.csv(file.path(a, "simulations.csv"))
sims <- sims[sims$mode == "main" & sims$status == "ok", ]
census <- read.csv(file.path(a, "census.csv.gz"))
census <- merge(census, sims[c("simulation", "generations")])
census <- census[census$generation == census$generations, c("simulation", "indicator", "value")]
cndd <- read.csv(file.path(a, "cndd.csv.gz"))
cndd <- data.frame(simulation = cndd$simulation[cndd$scope == "pooled"], indicator = "cndd",
                   value = cndd$slope[cndd$scope == "pooled"])
runtime <- data.frame(simulation = sims$simulation, indicator = "runtime", value = sims$runtime)
params <- setdiff(names(sims), c("mode", "status", "error", "runtime", "seed"))
out <- merge(sims[params], rbind(census, cndd, runtime))
con <- gzfile(file.path("model-behavior", "runs", paste0(basename(a), ".csv.gz")), "w")
write.csv(out, con, row.names = FALSE)
close(con)
message("written: ", basename(a))
