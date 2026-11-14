library(EcoPhyloSim)

runs = 3

timings <- matrix(nrow = 8, ncol = runs)

pars = list()
pars[[1]] = createCompletePar(x = 100, y = 100, runs = 2000, scenario = "Global Neutral")
pars[[2]] = createCompletePar(x = 100, y = 100, runs = 2000, density = 1, scenario = "Global Density")
pars[[3]] = createCompletePar(x = 100, y = 100, runs = 2000, environment = 1, scenario = "Global Environment")
pars[[4]] = createCompletePar(x = 100, y = 100, runs = 2000, density = 1, environment = 1, scenario = "Global Both")
pars[[5]] = createCompletePar(x = 50, y = 50, runs = 500, dispersal = 4, scenario = "Local Neutral")
pars[[6]] = createCompletePar(x = 50, y = 50, runs = 500, dispersal = 4, density = 1, scenario = "Local Density")
pars[[7]] = createCompletePar(x = 50, y = 50, runs = 500, dispersal = 4, environment = 1, scenario = "Local Environment")
pars[[8]] = createCompletePar(x = 50, y = 50, runs = 500, dispersal = 4, density = 1, environment = 1, scenario = "Local Both")

for (j in 1:length(pars)){
  par = pars[[j]]
  for (i in 1:ncol(timings)){
    timings[j,i] = system.time(out <- runSimulation(par))[3] 
  }
}

oldpar <- par(mfrow = c(2,4))
apply(timings,1,hist)
par(oldpar)

save(timings, file = "referenceTimes.RData")



