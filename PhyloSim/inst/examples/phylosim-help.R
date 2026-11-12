data(simu.envcom.local)
plotSpatialPhylo(simu.envcom.local)

par <- createCompletePar(x = 50, y = 50, dispersal = "global",
                         runs = c(500,1000), density = 1,
                         environment = 0.5, specRate = 1)
new.simulation <- runSimulation(par)
plotSpatialPhylo(new.simulation)