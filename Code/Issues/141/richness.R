# #141: species richness against dispersal, before and after the fix. Neutral
# model with local dispersal, 30 x 30 cells, specRate = 2, 1000 generations,
# dispersal from 0.5 to 3 in steps of 0.1, seeds 1:10. Before the fix, values
# between multiples of 0.5 give the same result as the next lower multiple.
#   Rscript richness.R <libraryBefore> before.csv
#   Rscript richness.R <libraryAfter> after.csv
# Results: see README.md. Run with R 4.6.2, master 1c778ea (before) and branch
# i141-dispersal (after); 3 parallel workers, a few minutes each.
args <- commandArgs(TRUE)
lib <- args[1]
grid <- expand.grid(dispersal = seq(0.5, 3, by = 0.1), seed = 1:10)
cl <- parallel::makeCluster(3)
invisible(parallel::clusterCall(cl, function(lib) library(EcoPhyloSim, lib.loc = lib), lib))
grid$richness <- parallel::parApply(cl, grid, 1, function(g) {
  par <- createCompletePar(x = 30, y = 30, dispersal = g[["dispersal"]], runs = 1000, specRate = 2,
                           seed = g[["seed"]])
  length(unique(as.vector(runSimulation(par, progress = FALSE)$Output[[1]]$specMat)))
})
parallel::stopCluster(cl)
write.csv(grid, args[2], row.names = FALSE)
