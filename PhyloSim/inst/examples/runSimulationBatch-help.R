## Define two parameter sets
par1 <- createCompletePar(x = 50, y = 50, dispersal = 0 , runs = c(500,1000),
        negativeDD = 0, environment = 0.5, specRate = 1)
par2 <- createCompletePar(x = 50, y = 50, dispersal = 0 , runs = c(500,1000),
        negativeDD = 1, environment = 0.5, specRate = 1)

## Merge the parameter sets. It is important to note that the function
## needs a list of parameter sets.
par <- list(par1,par2)

## Run the model
simu <- runSimulationBatch(par, parallel = 2)


## Compare the results, here the SAC is shown as an example.
par(mfrow=c(1,2))
sac(simu[[1]])
sac(simu[[2]])
