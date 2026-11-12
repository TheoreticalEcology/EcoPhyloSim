context("run Simulation")
test_that("Model outputs are  as expected",{
  skip("reference data PhyloSimTestData is missing")
  # Here are the model runs that created the .RData files
#   
#   sacVal1<-sacVal2<-sacVal3<-sacVal4<-sacVal5<-sacVal6<-matrix(0,10,100)
#   
#   for(i in 1:100){
#     #Neutral Models
#   par1<-createCompletePar(x=100,y=100, dispersal="global", runs=1000, seed=i)
#   par2<-createCompletePar(x=100,y=100, dispersal=1, runs=1000, seed=i)
#   
#     #Other Model runs
#   par3<-createCompletePar(x=100,y=100, dispersal="global", runs=1000,
#                          density=0.5, environment=0.5, seed=i)
#   par4<-createCompletePar(x=100,y=100, dispersal=1, runs=1000,
#                          density=0.5, environment=0.5, seed=i) 
#   
#   par5<-createCompletePar(x=100,y=100, dispersal="global", runs=1000,
#                          density=1, environment=1, seed=i)
#   par6<-createCompletePar(x=100,y=100, dispersal=1, runs=1000,
#                          density=1, environment=1, seed=i) 
#   
#   parList<-list(par1,par2,par3,par4,par5,par6)
#   
#   simu<-runSimulationBatch(parList, parallel = 15)
#   
#   sacVal1[,i]<- sac(simu[[1]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
#   sacVal2[,i]<- sac(simu[[2]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
#   sacVal3[,i]<- sac(simu[[3]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
#   sacVal4[,i]<- sac(simu[[4]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
#   sacVal5[,i]<- sac(simu[[5]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
#   sacVal6[,i]<- sac(simu[[6]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
#   }
#   
#   
#   ## Save the data as a list in inst/extdata
#   simulations<-list(sacVal1,sacVal2,sacVal3,sacVal4,sacVal5,sacVal6)
#   save(simulations, file="PhyloSimTestData.RData")

  ## Now the model is tested
  sac1<-sac2<-sac3<-sac4<-sac5<-sac6<-matrix(0,10,10)
  
  for(i in 1:10){
    #Neutral Models
    par1<-createCompletePar(x=100,y=100, dispersal="global", runs=1000, seed=i)
    par2<-createCompletePar(x=100,y=100, dispersal=1, runs=1000, seed=i)
    
    #Other Model runs
    par3<-createCompletePar(x=100,y=100, dispersal="global", runs=1000,
                            density=0.5, environment=0.5, seed=i)
    par4<-createCompletePar(x=100,y=100, dispersal=1, runs=1000,
                            density=0.5, environment=0.5, seed=i) 
    
    par5<-createCompletePar(x=100,y=100, dispersal="global", runs=1000,
                            density=1, environment=1, seed=i)
    par6<-createCompletePar(x=100,y=100, dispersal=1, runs=1000,
                            density=1, environment=1, seed=i) 
    
    parList<-list(par1,par2,par3,par4,par5,par6)
    
    simu<-runSimulationBatch(parList, parallel = 2)
    
    
    
    ## Need to change to new format --> replace are with respective size.
    sac1[,i]<- sac(simu[[1]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
    sac2[,i]<- sac(simu[[2]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
    sac3[,i]<- sac(simu[[3]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
    sac4[,i]<- sac(simu[[4]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
    sac5[,i]<- sac(simu[[5]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
    sac6[,i]<- sac(simu[[6]],size=c(2,3,4,6,8,11,16,23,32,45), plot=FALSE)$sr.Mean
  }
  
  sac_run1<-sac_run2<-sac_run3<-sac_run4<-sac_run5<-sac_run6<-numeric()
  
  for(i in 1:nrow(sac1)){
    sac_run1[i] <- mean(sac1[i,])
    sac_run2[i] <- mean(sac2[i,])
    sac_run3[i] <- mean(sac3[i,])
    sac_run4[i] <- mean(sac4[i,])
    sac_run5[i] <- mean(sac5[i,])
    sac_run6[i] <- mean(sac6[i,])
  }
  

  ##Get the data
  data(PhyloSimTestData)
  
  sacVal1<-simulations[[1]]
  sacVal2<-simulations[[2]]
  sacVal3<-simulations[[3]]
  sacVal4<-simulations[[4]]
  sacVal5<-simulations[[5]]
  sacVal6<-simulations[[6]]
  
  for(i in 1:10){
    expect_true(sac_run1[i]<=quantile(sacVal1[i,],probs=0.975)& sac_run1[i]>=quantile(sacVal1[i,],probs=0.025))
    expect_true(sac_run2[i]<=quantile(sacVal2[i,],probs=0.975)& sac_run2[i]>=quantile(sacVal2[i,],probs=0.025))
    expect_true(sac_run3[i]<=quantile(sacVal3[i,],probs=0.975)& sac_run3[i]>=quantile(sacVal3[i,],probs=0.025))
    expect_true(sac_run4[i]<=quantile(sacVal4[i,],probs=0.975)& sac_run4[i]>=quantile(sacVal4[i,],probs=0.025))
    expect_true(sac_run5[i]<=quantile(sacVal5[i,],probs=0.975)& sac_run5[i]>=quantile(sacVal5[i,],probs=0.025))
    expect_true(sac_run6[i]<=quantile(sacVal6[i,],probs=0.975)& sac_run6[i]>=quantile(sacVal6[i,],probs=0.025))
  }
  
  
})
