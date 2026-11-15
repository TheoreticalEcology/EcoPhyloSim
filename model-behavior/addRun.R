# Exports the summaries of one run of PhyloSim/data-raw/packageBehavior.R as
# CSV files for model-behavior.qmd, which compares the model versions.
# Run from the repository root after combining the run:
#   Rscript model-behavior/addRun.R [path to packageBehavior.rda]
args <- commandArgs(trailingOnly = TRUE)
load(if (length(args)) args[1] else "PhyloSim/data/packageBehavior.rda")
m <- packageBehavior$main
id <- paste(format(as.Date(packageBehavior$info$combined)), packageBehavior$info$commit, sep = "_")
census <- merge(m$settings, m$census)
cndd <- merge(m$settings, transform(m$cndd, indicator = "cndd", generation = NA))
out <- rbind(census, cndd[names(census)])
out$mean <- signif(out$mean, 4)
out$sd <- signif(out$sd, 4)
write.csv(out, file.path("model-behavior", "runs", paste0(id, ".csv")), row.names = FALSE)
message("written: ", id)
