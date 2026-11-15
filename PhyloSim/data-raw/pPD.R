# Creates the example data sets pPD.pvalues and pPD.positions for
# plotPhylogeneticDispersion(): 28 scenarios (4 dispersal modes x 7 combinations
# of density dependence and environmental filtering) on a 50 x 50 grid, and the
# p-values of the "PhylMeta" null model for 200 subplots of each size 2 to 10.
#
# Run from the repository root, with the package installed:
#   Rscript PhyloSim/data-raw/pPD.R [cores] [generations]
# The results are illustrative only (one run per scenario).

library(EcoPhyloSim)
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 1
generations <- if (length(args) >= 2) as.integer(args[2]) else 50000

# scenario order and labels as in the original data set: within each dispersal
# mode, "neutral", density 0 to 1 together with environment 1 to 0, and "both"
fit <- data.frame(density = c(0, seq(0, 1, len = 5), 1),
                  environment = c(0, seq(1, 0, len = 5), 1))
dispersal <- c("global", "0.5", "2", "4")
scen <- merge(fit, data.frame(dispersal = dispersal), sort = FALSE)
scen <- scen[order(match(scen$dispersal, dispersal)), ]

pars <- lapply(seq_len(nrow(scen)), function(i)
  createCompletePar(x = 50, y = 50, runs = generations, specRate = 2,
                    dispersal = if (scen$dispersal[i] == "global") "global"
                                else as.numeric(scen$dispersal[i]),
                    density = scen$density[i], environment = scen$environment[i],
                    scenario = paste0("disp=", scen$dispersal[i], " dens=", scen$density[i],
                                      " env=", scen$environment[i]),
                    seed = 1000 + i))

pvalues <- function(p) {
  set.seed(p$seed)
  simu <- runSimulation(p)
  calculatePhylogeneticDispersion(simu, plotlengths = 2:10, plots = 200,
                                  replicates = 500, type = "PhylMeta")[[1]][[1]]
}

if (cores > 1) {
  cl <- parallel::makeCluster(cores)
  parallel::clusterEvalQ(cl, library(EcoPhyloSim))
  res <- parallel::parLapplyLB(cl, pars, pvalues, chunk.size = 1)
  parallel::stopCluster(cl)
} else res <- lapply(pars, pvalues)

pPD.pvalues <- list(res)
pPD.positions <- list(x = 1:7, y = c(4, 1, 2, 3),
                      zname = c("neutral", "1", "0.75", "0.5", "0.25", "0", "both"),
                      yname = c("global", "0.5", "2", "4"),
                      xname = c("neutral", "0", "0.25", "0.5", "0.75", "1", "both"))
save(pPD.pvalues, file = "PhyloSim/data/pPD.pvalues.rda", compress = "xz")
save(pPD.positions, file = "PhyloSim/data/pPD.positions.rda", compress = "xz")
