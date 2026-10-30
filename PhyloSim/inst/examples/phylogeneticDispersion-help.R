
## Create a set of parameters
par1 <- createCompletePar(x = 50, y = 50, dispersal = 0 , runs = c(500,1000),
        density = 0, environment = 0.5, specRate = 1)
par2 <- createCompletePar(x = 50, y = 50, dispersal = 0 , runs = c(500,1000),
        density = 1, environment = 0.5, specRate = 1)

## Merge the parameter sets.
par <- list(par1,par2)

## Run the model
simu <- runSimulationBatch(par, parallel = 2)

## Calculate null model and compare the observed results against
pValues <- calculatePhylogeneticDispersion(simu, plotlength = 20, plots = 20,
                                           replicates = 20, type = "PhylMeta")