

library(EcoPhyloSim)
 

par <- createCompletePar(x = 50, y = 50, dispersal = 1, runs = c(500,1000), density = 1, environment = 0.5, specRate = 1)

system.time(simu <- runSimulation(par))

plotSpatialPhylo(simu, plot = "both")

plotTraitDistribution(simu)

dis <- calculatePhylogeneticDispersion(simu)


plot(simu[[1]]$phylogeny)
  
  unique(as.vector(simu[[1]]$specMat))

plot(out, time = 1)

par(mfrow = c(1,2))
rac(simu$specMat)
sac(simu$specMat)


par <- createCompletePar(runs = c(100,200,300))
simu <- runSimulation(par)
str(simu)



simu <- NeutralMod()
str(simu)
