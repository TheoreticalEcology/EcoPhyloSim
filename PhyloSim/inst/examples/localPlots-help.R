library(EcoPhyloSim)
xSize <- 30
ySize <- 30
nGenerations <- 300

par <- createCompletePar(x = xSize,
                         y = ySize,
                         runs = nGenerations)

sim <- runSimulation(par = par)

localPlots(simu = sim, size = 2, n = 3, plot = TRUE)

