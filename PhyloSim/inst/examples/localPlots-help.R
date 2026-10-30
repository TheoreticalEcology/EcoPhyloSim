library(PhyloSim)
xSize <- 30
ySize <- 30
nGenerations <- 300

#dispersal <- c(0, seq(1, 10.5, 0.5))
#specRate <- seq(0, 20, 1)
#environment <- seq(0, 10, 0.5)
#densitiy <- seq(0, 10, 0.5)

par <- createCompletePar(x = xSize,
                         y = ySize,
                         runs = nGenerations)

sim <- runSimulation(par = par)

localPlots(simu = sim, size = 2, n = 3, plot = TRUE)

