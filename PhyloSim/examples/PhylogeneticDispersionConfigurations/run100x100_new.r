require(EcoPhyloSim)

# SET WD TO SOURCE FILE LOCATION !!!!

#old settings:
#runs = seq(25000, 50000, len = 50)
#setwd("C:/Users/fhartig/Desktop/phylosim/PhyloSim/examples/PhylogeneticDispersionConfigurations")


setwd("C:/Users/Stefan/Desktop/HIWI/nullModel tests")


runs = 2000

dispOptions = 4
fitOptions = 7

dispersal =  rep(c(0,0.5, 2,4), each = fitOptions)
density = rep(c(0,seq(0,1,len = fitOptions-2),1),dispOptions)
environment = rep(c(0,seq(1,0,len = fitOptions-2),1),dispOptions)


scenarios <- paste(modes, " dens=", density, " env=", environment, sep = "")


pars <- list()

#for (i in 1:length(scenarios)){pars[[i]] = createCompletePar(x = 100,y = 100, runs = runs,scenario = scenarios[i], dispersal = dispersal[i],specRate = 1,density = density[i],environment = environment[i], fitnessBaseMortalityRatio = 5, densityCut = 1,seed = 1000)}
for (i in 1:length(scenarios)){pars[[i]] = createCompletePar(x = 20,y = 20, runs = runs,scenario = scenarios[i], dispersal = dispersal[i],specRate = 2,density = density[i],environment = environment[i], fitnessBaseMortalityRatio = 5, densityCut = 1,seed = 1000)}



system.time(simulationOut <- runSimulationBatch(pars, parallel = "auto"))

save(simulationOut, pars, file = "results100x100.Rdata" )
