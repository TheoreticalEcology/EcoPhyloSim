require(EcoPhyloSim)

# SET WD TO SOURCE FILE LOCATION !!!!

setwd("C:/Users/fhartig/Desktop/phylosim/PhyloSim/examples/PhylogeneticDispersionConfigurations")

runs = seq(25000, 50000, len = 50)

dispOptions = 4
fitOptions = 7

dispersal =  rep(c(0,0.5, 2,4), each = fitOptions)
density = rep(c(0,seq(0,1,len = fitOptions-2),1),dispOptions)
environment = rep(c(0,seq(1,0,len = fitOptions-2),1),dispOptions)

modes <- ifelse(dispersal == 0, "global", "local")
scenarios <- paste(modes, " dens=", density, " env=", environment, sep = "")


pars <- list()

for (i in 1:length(scenarios)){pars[[i]] = createCompletePar(x = 100,y = 100, runs = runs, modes = modes[i],scenario = scenarios[i], dispersal = dispersal[i],specRate = 1,density = density[i],environment = environment[i], fitnessBaseMortalityRatio = 5, densityCut = 1,seed = 1000)}


system.time(simulationOut <- fullModBatch(pars, parallel = T)) 

save(simulationOut, pars, file = "results100x100.Rdata" )
