
#Load data
data(simu.neutral.global)
data(simu.neutral.local)
data(simu.envcom.local)



#Plot clades for neutral model with global dispersal with three subplots
par(mfrow=c(1,3))
highlightLocalPhylo(simu.neutral.global, size = 50, n = 3)



#Plot clades for neutral model with local dispersal with three subplots
par(mfrow=c(1,3))
highlightLocalPhylo(simu.neutral.local, size = 50, n = 3)



#Plot clades for environment and competition model with global dispersal with three subplots
par(mfrow=c(1,3))
highlightLocalPhylo(simu.envcom.local, size = 50, n = 3)



#increasing plot size
par(mfrow=c(2,2))
highlightLocalPhylo(simu.envcom.local, size = 50, n =2)
highlightLocalPhylo(simu.envcom.local, size = 50, n =2)