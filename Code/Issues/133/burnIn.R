# #133: how many generations of burn-in negative density dependence needs before
# it acts on mortality, with the trait kernels ("absolute", "gaussian") and the
# "species" kernel, with and without a jump of the competition trait at
# speciation (traitSpeciationJump) and with faster trait evolution
# (traitEvolution). Setting of the issue: 50 x 50 cells, dispersal = 1,
# negativeDDCut = 2, specRate = 1 (default), 6 consecutive censuses
# (g + 0:5) at g = 100, 300, 1000, 3000, 5000, 10000, 30000; 10 seeds per setting.
# Per census block: mortality of individuals with 0-3 and 9-12 conspecific
# neighbors (radius 2, getMortalityData()), logit slope of mortality per
# conspecific neighbor, sd of the competition trait, species richness.
# Results: see README.md; raw per-run tables go to results/ (not committed),
# the summary table summary.csv is committed.
# Run with EcoPhyloSim at commit 1c778ea (0.3.2.9000), R 4.6.2, Windows 11;
# 6 parallel workers, about 12 minutes. Run from this folder:
#   Rscript burnIn.R
library(parallel)

gens <- c(100, 300, 1000, 3000, 5000, 10000, 30000)
seeds <- 1:10
settings <- data.frame(
  name = c("none", "abs_w0", "abs_w0.1", "abs_w0.3", "abs0.5_w0", "abs0.5_w0.3",
           "abs_w0_evol0.05", "gauss_w0", "gauss_w0.3", "species"),
  negativeDD = c(0, 1, 1, 1, 0.5, 0.5, 1, 1, 1, 1),
  kernel = c("absolute", "absolute", "absolute", "absolute", "absolute", "absolute",
             "absolute", "gaussian", "gaussian", "species"),
  jump = c(0, 0, 0.1, 0.3, 0, 0.3, 0, 0, 0.3, 0),
  evolution = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05, 0.01, 0.01, 0.01))
jobs <- expand.grid(setting = seq_len(nrow(settings)), seed = seeds)

runOne <- function(i) {
  s <- settings[jobs$setting[i], ]; seed <- jobs$seed[i]
  file <- file.path("results", sprintf("%s_seed%02d.rds", s$name, seed))
  if (file.exists(file)) return(readRDS(file))
  par <- EcoPhyloSim::createCompletePar(
    x = 50, y = 50, dispersal = 1, negativeDD = s$negativeDD, negativeDDCut = 2,
    negativeDDKernel = s$kernel, traitSpeciationJump = c(0, s$jump, 0),
    traitEvolution = c(0.01, s$evolution, 0.01), runs = c(outer(0:5, gens, "+")),
    seed = seed)
  invisible(capture.output(simu <- EcoPhyloSim::runSimulation(par)))
  out <- do.call(rbind, lapply(gens, function(g) {
    sub <- simu
    idx <- which(par$runs %in% (g + 0:5))
    sub$Output <- simu$Output[idx]; sub$Model$runs <- par$runs[idx]
    d <- EcoPhyloSim::getMortalityData(sub)
    o <- simu$Output[[idx[1]]]
    data.frame(setting = s$name, seed = seed, generation = g,
               mortFew = mean(d$died[d$con <= 3]), mortMany = mean(d$died[d$con >= 9]),
               nMany = sum(d$con >= 9),
               slope = unname(coef(glm(died ~ con, family = binomial, data = d))[2]),
               sdComp = sd(c(o$compMat)), richness = length(unique(c(o$specMat))))
  }))
  saveRDS(out, file)
  out
}

dir.create("results", showWarnings = FALSE)
cl <- makeCluster(6)
clusterExport(cl, c("settings", "jobs", "gens"))
res <- do.call(rbind, parLapply(cl, seq_len(nrow(jobs)), runOne))
stopCluster(cl)

# mean and sd over seeds per setting and generation
# (mortMany is NaN for runs without individuals with >= 9 conspecific neighbors)
agg <- aggregate(cbind(mortFew, mortMany, nMany, slope, sdComp, richness) ~ setting + generation,
                 data = res, FUN = mean, na.rm = TRUE, na.action = na.pass)
sdSlope <- aggregate(slope ~ setting + generation, data = res, FUN = sd)
agg$slopeSd <- sdSlope$slope[match(paste(agg$setting, agg$generation),
                                   paste(sdSlope$setting, sdSlope$generation))]
agg <- agg[order(match(agg$setting, settings$name), agg$generation), ]
agg[, -(1:2)] <- signif(agg[, -(1:2)], 3)
write.csv(agg, "summary.csv", row.names = FALSE)
print(agg, row.names = FALSE)
