# #124: expected species richness es() of the neutral model with global dispersal.
# 1. Checks that theta = specRate * (J - 1) / (J - specRate), used by es() since #124, equals
#    theta = (J - 1) * nu / (1 - nu) of Vallade & Houchmandzadeh (2003, Phys. Rev. E 68, 061902,
#    eqn 6), with nu = specRate / J, and compares es() with the old theta = 2 * specRate.
# 2. Compares es() with the mean richness of runSimulation() (global dispersal, recorded once at
#    the end of runs of 5 * J generations, 100 seeds per setting).
# Run from the repository root: Rscript Code/Issues/124/esTheta.R (about 5 minutes, 1 thread).
# Results are written to Code/Issues/124/esTheta.csv; see README.md for the numbers and decisions.
# EcoPhyloSim 0.3.2.9000 (branch i124-esTheta), R 4.6.2; seeds 1 to 100.

devtools::load_all("PhyloSim", quiet = TRUE)

# 1. theta
for (J in c(6, 100, 2500)) for (specRate in c(0.5, 2, 4)) {
  nu <- specRate / J
  stopifnot(all.equal(specRate * (J - 1) / (J - specRate), (J - 1) * nu / (1 - nu)))
}
oldEs <- function(specRate, n) sum(2 * specRate / (2 * specRate + 0:(n^2 - 1)))
cat("es() with theta = 2 * specRate relative to es() since #124 (30 x 30 cells):\n")
print(sapply(c(1, 2, 4), function(r) round(oldEs(r, 30) / es(r, 30), 2)))

# 2. simulation
settings <- expand.grid(n = c(10, 20, 30), specRate = c(1, 2, 4))
res <- do.call(rbind, lapply(seq_len(nrow(settings)), function(i) {
  n <- settings$n[i]; specRate <- settings$specRate[i]
  rich <- sapply(1:100, function(seed) {
    invisible(capture.output(s <- runSimulation(createCompletePar(
      x = n, y = n, dispersal = "global", specRate = specRate, runs = 5 * n^2, seed = seed))))
    specRich(s)
  })
  data.frame(cells = n^2, specRate = specRate, mean = mean(rich), se = sd(rich) / sqrt(length(rich)),
             es = es(specRate, n), esOld = oldEs(specRate, n))
}))
res$ratio <- res$mean / res$es
res$excess <- res$mean - res$es
print(res, digits = 3)
write.csv(res, "Code/Issues/124/esTheta.csv", row.names = FALSE)
