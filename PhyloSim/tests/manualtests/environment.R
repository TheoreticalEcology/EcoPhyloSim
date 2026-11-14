# params: a,b,c like: a*x + b*y + c
# w, h: width and height
# scale: whether the gradient schould be scaled from 0 to 1
createLinearGradientEnv <- function (params, w, h, scale=TRUE) {
  out <- matrix(NA, ncol = w, nrow = h)
  for (x in 1:w) {
    for (y in 1:h) {
      out[y, x] <- x * params[1] + y * params[2] + params[3]
    }
  }
  if (scale == TRUE) out <- ( (out - min(out)) / (max(out) - min(out)) )
  return(out)
}


library(EcoPhyloSim)
w <- h <- 50
p <- c(1.5, 1, 0)
envMat <- createLinearGradientEnv(p, w, h, T)

nRuns <- 1000
dispersal <- 0
envStrength <- 1 # strong influence
competition <- 0

par1 <- createCompletePar(w, h, runs = nRuns,
                         dispersal = dispersal,
                         environment = envStrength, airmat = envMat,
                         density = competition)
sim1 <- runSimulation(par1)
plotTraitDistribution(sim1, type = "all")

# --------------------------------
envStrength <- 0.1 # weak influence

par2 <- createCompletePar(w, h, runs = nRuns,
                          dispersal = dispersal,
                          environment = envStrength, airmat = envMat,
                          density = competition)
sim2 <- runSimulation(par2)
plotTraitDistribution(sim2, type = "all")


# --------------------------------
envStrength <- 0 # no influence

par3 <- createCompletePar(w, h, runs = nRuns,
                          dispersal = dispersal,
                          environment = envStrength, airmat = envMat,
                          density = competition)
sim3 <- runSimulation(par3)
plotTraitDistribution(sim3, type = "all")
