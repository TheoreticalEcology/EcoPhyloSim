
### Usage with runSimulation
## Define a parameter set
par <- createCompletePar(x = 50, y = 50, dispersal = 0,
                         runs = c(500,1000), negativeDD = 1,
                         environment = 0.5, specRate = 1)

## Run Model
simu <- runSimulation(par)

## Compare to null model
nullModel(simu, abundance= FALSE, localPlotSize = 10, numberOfPlots = 10, repetitions = 10)

## Usage with runSimulationBatch
## Define two (or multiple) parameter sets
par1 <- createCompletePar(x = 50, y = 50, dispersal = 0 , runs = c(500,1000),
        negativeDD = 1, environment = 0.5, specRate = 1)
par2 <- createCompletePar(x = 50, y = 50, dispersal = 0 , runs = c(500,1000),
        negativeDD = 0, environment = 0.5, specRate = 1)

## Merge the parameter sets. It is important to note that the function
## needs a list of parameter sets.
par <- list(par1,par2)

## Run the model
simubatch <- runSimulationBatch(par, parallel = FALSE)

## Compare results to null model
nullModel(simubatch[[1]],localPlotSize = 10, numberOfPlots = 10, repetitions = 10)
nullModel(simubatch[[2]],localPlotSize = 10, numberOfPlots = 10, repetitions = 10)
