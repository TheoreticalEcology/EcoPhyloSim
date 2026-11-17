# Creates the example data sets pPD.pvalues and pPD.positions for
# plotPhylogeneticDispersion(): 28 scenarios (4 dispersal modes x 7 combinations
# of density dependence and environmental filtering) on a 50 x 50 grid, and the
# p-values of the "PhylMeta" null model for 200 subplots of each size 2 to 10.
#
# Run at each release from the repository root, with the release commit checked
# out (its package is installed into lib/ if missing):
#   Rscript PhyloSim/data-raw/model-behavior/pPD.R [cores] [generations]
# The results are illustrative only (one run per scenario). The attribute
# "provenance" of pPD.pvalues records the run (#180).

source("PhyloSim/data-raw/model-behavior/common.R")
# the package of the checked-out commit (installed into lib/ if missing), also for the
# parallel workers (R_LIBS before they start, #190)
lib <- libraryOf("HEAD")
Sys.setenv(R_LIBS = lib)
library(EcoPhyloSim, lib.loc = lib)
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
                    negativeDD = scen$density[i], environment = scen$environment[i],
                    scenario = paste0("disp=", scen$dispersal[i], " dens=", scen$density[i],
                                      " env=", scen$environment[i]),
                    seed = 1000 + i))

# a failing scenario returns its error message instead of stopping the whole run (#180)
pvalues <- function(p) tryCatch({
  set.seed(p$seed)
  simu <- runSimulation(p)
  calculatePhylogeneticDispersion(simu, plotlengths = 2:10, plots = 200,
                                  replicates = 500, type = "PhylMeta")[[1]][[1]]
}, error = function(e) conditionMessage(e))
runAll <- function(pars) {
  if (cores <= 1) return(lapply(pars, pvalues))
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, library(EcoPhyloSim))
  parallel::parLapplyLB(cl, pars, pvalues, chunk.size = 1)
}
res <- runAll(pars)
# failed scenarios are run once more with another seed
failed <- which(vapply(res, is.character, logical(1)))
if (length(failed)) {
  message("rerunning failed scenarios ", paste(failed, collapse = ", "), " with other seeds")
  res[failed] <- runAll(lapply(pars[failed], function(p) { p$seed <- p$seed + 1e6; p }))
  attr(res, "reseeded") <- failed # scenarios run with seed + 1e6
}
failed <- which(vapply(res, is.character, logical(1)))
if (length(failed)) {
  partial <- dataDir("cache", "pPD-partial.rds")
  dir.create(dirname(partial), showWarnings = FALSE)
  saveRDS(res, partial)
  stop("scenarios ", paste(failed, collapse = ", "), " failed (", paste(unique(unlist(res[failed])), collapse = "; "),
       "); results of the others are in ", partial)
}

pPD.pvalues <- list(res)
attr(pPD.pvalues, "provenance") <- c(
  date = format(Sys.time(), "%Y-%m-%d %H:%M"), commit = gitDescribe(),
  package = paste("EcoPhyloSim", packageVersion("EcoPhyloSim", lib.loc = lib)),
  R = R.version.string, platform = R.version$platform, generations = generations,
  reseeded = paste(attr(res, "reseeded"), collapse = " "))
pPD.positions <- list(x = 1:7, y = c(4, 1, 2, 3),
                      zname = c("neutral", "1", "0.75", "0.5", "0.25", "0", "both"),
                      yname = c("global", "0.5", "2", "4"),
                      xname = c("neutral", "0", "0.25", "0.5", "0.75", "1", "both"))
save(pPD.pvalues, file = "PhyloSim/data/pPD.pvalues.rda", compress = "xz")
save(pPD.positions, file = "PhyloSim/data/pPD.positions.rda", compress = "xz")
