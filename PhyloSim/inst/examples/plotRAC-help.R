
#Load data
data(simu.neutral.global)
data(simu.neutral.local)
data(simu.envcom.local)



#Plot RAC curve for neutral model and global dispersal

rac(simu = simu.neutral.global)
rac(simu = simu.neutral.global, plot ="bar")



#Plot RAC curve for neutral model and local dispersal

rac(simu = simu.neutral.local)
rac(simu = simu.neutral.local, plot ="bar")



#Plot RAC curve for environment and competition model and local dispersal

rac(simu = simu.envcom.local)
rac(simu = simu.envcom.local, plot ="bar")


