#Load data
data(simu.neutral.global)
data(simu.neutral.local)
data(simu.envcom.local)



#neutral model with global dispersal
plotTraitDistribution(simu=simu.neutral.global)
plotTraitDistribution(simu=simu.neutral.global, type ="all")


#neutral model with local dispersal
plotTraitDistribution(simu=simu.neutral.local)
plotTraitDistribution(simu=simu.neutral.local, type = "all")



#environment and competition model with local dispersal
plotTraitDistribution(simu=simu.envcom.local)
plotTraitDistribution(simu=simu.envcom.local, type = "all")

