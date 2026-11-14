library(EcoPhyloSim)
library(picante)

getSummaries <- function(param){
  par <- createCompletePar(x = 30, y = 30, dispersal = 1, runs = round(param[1]), density = param[2], environment = param[3])
  simu <- runSimulation(par)
  #s1 = collessImbalance(simu$phylogeny)
  dat <- rac(simu$specMat, plot = F)
  s2 = coef(lm(Abundance ~ Rank + I(Rank^2), data = dat))
  return(s2)
}

# Reference

pars = matrix(rep(c(300, 0.5,0.5),40), ncol = 3, byrow = T)
test <- apply(pars, 1, getSummaries)

testdata <- apply(test, 1, median)
variance <- apply(test, 1, median)

pairs(t(test))

# Samples 

samples = 100

pars = data.frame(runs= sample(50:1000, samples), density = runif(samples, 0,1), environment = runif(samples, 0,1))

res <- apply(pars, 1, getSummaries)

acceptance <- sqrt(apply((res - testdata)^2/(variance^2), 2, sum)) < 0.5

oldpar <- par(mfrow = c(1,3))
hist(pars[acceptance, 1], xlim = c(50,1000))
hist(pars[acceptance, 2], xlim = c(0,1))
hist(pars[acceptance, 3], xlim = c(0,1))
par(oldpar)



