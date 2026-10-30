
#Load data

data(simu.neutral.global)
data(simu.neutral.local)
data(simu.envcom.local)

#neutral model with global dispersal
simu <- simu.neutral.global
   #phylogeny and traits
   plotSpatialPhylo(simu, plot = "phylogeny")

   #phylogeny and spatial distribution
   plotSpatialPhylo(simu, plot = "both", plotTraits = FALSE)

   #all three
   plotSpatialPhylo(simu, plot = "both", plotTraits = TRUE)


 #neutral model with local dispersal
 simu <- simu.neutral.local
   #phylogeny and traits
   plotSpatialPhylo(simu, plot = "phylogeny")

   #phylogeny and spatial distribution
   plotSpatialPhylo(simu, plot = "both", plotTraits = FALSE)

   #all three
   plotSpatialPhylo(simu, plot = "both", plotTraits = TRUE)


 #evironment and competition model with local dispersal
 simu <- simu.envcom.local
   #phylogeny and traits
   plotSpatialPhylo(simu, plot = "phylogeny")

   #phylogeny and spatial distribution
   plotSpatialPhylo(simu, plot = "both", plotTraits = FALSE)

   #all three
   plotSpatialPhylo(simu, plot = "both", plotTraits = TRUE)

#Set simulation run with which.result
#Here, two simulation runs are stored in simu (one with 500 timesteps, one with 1000)
  par <- createCompletePar(
     x = 50, y = 50, dispersal = 1,
     runs = c(500,1000), density = 1,
     environment = 0.5, specRate = 1)
simu <- runSimulation(par)

#Choose to take the first simulation run (default is the last)
plotSpatialPhylo(simu, which.result = 1)


