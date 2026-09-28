
#Load data
data(simu.neutral.global)
data(simu.neutral.local)
data(simu.envcom.local)



#Plot SAC curve for neutral model and global dispersion

sac(simu = simu.neutral.global)



#Plot SAC curve for neutral model and local dispersion

sac(simu = simu.neutral.local)



#Plot SAC curve for environment and competition model and local dispersion

sac(simu = simu.envcom.local)



#Plot SAC curve with random plotsize

sac(simu=simu.envcom.local, size = sort(sample(c(10:100), size = 10)))



#Plot SAC curve with different repititions
simu <- simu.envcom.local
par(mfrow=c(3,1))
sac(simu=simu.envcom.local, rep = 3)
sac(simu=simu.envcom.local, rep = 30)
sac(simu=simu.envcom.local, rep = 30)
