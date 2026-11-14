

library(EcoPhyloSim)

# Run the model

# dispersal = 1 
# dispersal = 3
# dispersal = 2 

simu <- fullMod(x = 100, y = 100, dispersal = 1, runs = round(seq(100,100000, len = 20)), specRate = 1.0, density = F, environment = T, fitnessBaseMortalityRatio = 5)

for (i in 1:20){
  plotSpatialPhylo(simu[[i]], plot = "both", main = i)  
}

specRich(simu$specMat)


par(mfrow = c(1,2))
rac(simu$specMat)
sac(simu$specMat)

plotTraitDistribution(simu)


