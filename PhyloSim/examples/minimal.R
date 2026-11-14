
library(EcoPhyloSim)

par <- createCompletePar(x = 50, y = 50, dispersal = 1, runs = c(500,1000), density = 1, environment = 0.5, specRate = 1)

system.time(simu <- runSimulation(par))

plotSpatialPhylo(simu, plot = "both")

plotTraitDistribution(simu)
