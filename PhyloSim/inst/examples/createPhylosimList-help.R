## Create two objects of type Phylosim
par1 <- createCompletePar(x = 50, y = 50, dispersal = "global" , runs = 500,
        negativeDD = 1, environment = 0.5, specRate = 1)
par2 <- createCompletePar(x = 50, y = 50, dispersal = "global" , runs = 500,
        negativeDD = 1, environment = 0.5, specRate = 2)

simu1 <- runSimulation(par1)
simu2 <- runSimulation(par2)

## Use the function
simuList <- createPhylosimList(list(simu1,simu2))

## The object can now be handed to functions designed for objects of type "PhylosimList".
