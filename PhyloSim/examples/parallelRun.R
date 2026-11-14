require(EcoPhyloSim)

# SET WD TO SOURCE FILE LOCATION !!!!

modes <- c(rep("global",4), rep("local",16))
modesDists <- paste(modes, c(rep("0",4),as.character(sort(rep(c(0.5,1,2,4),4)))),sep = "")
modesDistsMechs <- paste(modesDists, rep(c("Envi", "Neut", "Dens", "Hybr"),4), sep = "")

dispersal = dispersal = c(rep(0,4), rep(0.5,4), rep(1,4), rep(2,4), rep(4,4))[i]

density = rep(c(F,F,T,T),5)
environment = rep(c(T,F,F,T),5)

pars <- list()
n <- length(modes)

for (i in 1:20){pars[[i]] = createCompletePar(x = 25,y = 25,runs = round(seq(250, 300, len = 10)),modes = modes[i],scenario = modesDistsMechs[i], dispersal = dispersal[i],specRate = 1,density = density[i],environment = environment[i], fitnessBaseMortalityRatio = 10, densityCut = 1,seed = 1000)}





system.time(simulationOut <- fullModBatch(pars, parallel = T)) 

save(simulationOut, pars, file = "results3.Rdata" )
