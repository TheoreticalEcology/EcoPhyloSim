

library(EcoPhyloSim)
par <- createCompletePar(x = 50, y = 50, dispersal = "global" , runs = 10)
for (i in 1:1000) simu <- runSimulation(par)
plot(simu)

# before approx 500 MB additional memory, 1 MB per run

library(EcoPhyloSim)
par <- createCompletePar(x = 10, y = 10, dispersal = "global" , runs = 10)
for (i in 1:500) simu <- runSimulation(par)

# before approx 160 MB additional memory, 0.32 MB per run 

library(EcoPhyloSim)
par <- createCompletePar(x = 10, y = 10, dispersal = "global" , runs = 1000)
for (i in 1:500) simu <- runSimulation(par)

# before approx 410 MB additional memory, 0.8 MB per run 


