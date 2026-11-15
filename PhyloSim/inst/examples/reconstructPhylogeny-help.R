# Run a simulation
par <- createCompletePar(x = 50, y = 50, dispersal = "global" , runs = 500,
        negativeDD = 1, environment = 0.5, specRate = 1)

simu <- runSimulation(par)

## Reconstruct phylogeny
rePhyl <- phyloReconstruct(simu)

## Compare to the real phylogeny (requires package 'ape')
par(mfrow=c(1,2))
plot(ape::drop.fossil(simu$Output[[1]]$phylogeny))
plot(rePhyl)