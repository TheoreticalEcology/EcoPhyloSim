# Summary statistics of runs with global dispersal and fitness on recruitment, to compare two builds whose
# results differ for a seed but should have the same distribution (binned parent proposals, #207).
# usage: Rscript recruitment-check.R <library folder> <tag> <results file> <first seed> <last seed> [settings]
# settings: comma-separated names from the list below (default: all); appends one row per run (semicolon-separated)
a <- commandArgs(TRUE)
suppressMessages(library(EcoPhyloSim, lib.loc = a[1]))
MR <- c("mortality", "recruitment")
settings <- list(env = list(environment = 1, fitnessActsOn = "recruitment"),
                 envDD = list(environment = 1, negativeDD = 1, fitnessActsOn = MR),
                 narrow = list(environment = 1, nicheWidth = 0.01, fitnessActsOn = "recruitment"))
names <- if (length(a) >= 6) strsplit(a[6], ",")[[1]] else names(settings)
for (n in names) for (seed in as.integer(a[4]):as.integer(a[5])) {
  p <- do.call(createCompletePar, c(list(x = 30, y = 30, dispersal = "global", runs = 300, specRate = 2, seed = seed),
                                    settings[[n]]))
  t <- system.time(invisible(capture.output(sim <- suppressWarnings(runSimulation(p)))))[["user.self"]]
  o <- sim$Output[[1]]
  ab <- table(o$specMat) / length(o$specMat)
  row <- data.frame(tag = a[2], setting = n, seed = seed, richness = length(ab), shannon = -sum(ab * log(ab)),
                    mismatch = mean(abs(o$traitMat - o$envMat)), meanTrait = mean(o$traitMat),
                    sdTrait = sd(c(o$traitMat)), meanComp = mean(o$compMat), cpu_s = t)
  write.table(row, a[3], sep = ";", row.names = FALSE, col.names = !file.exists(a[3]), append = file.exists(a[3]))
}
