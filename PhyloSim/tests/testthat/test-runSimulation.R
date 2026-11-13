context("run Simulation")

test_that("Simulation is running correctly",{
  par<-createCompletePar(x=50, y=50, dispersal="global",density=1)
  parRneutral<-createCompletePar(x=50, y=50, dispersal="global",density=1, type="Rneutral")
  simu<-runSimulation(par)
  simuRneutral<-runSimulation(parRneutral)
  
  expect_output(str(par), "List of 22")
  expect_is(simu, "PhyloSim")
  expect_is(simuRneutral, "PhyloSim")
  
  })

test_that("SimulationBatch is running correctly",{
  par1<-createCompletePar(x=50, y=50, dispersal="global",density=1)
  par2<-createCompletePar(x=50, y=50, dispersal="global",density=0.5)
  parlist<-list(par1,par2)
  simuBatch<-runSimulationBatch(parlist)

  expect_is(simuBatch, "PhylosimList")
  expect_is(simuBatch[[1]], "PhyloSim")
  expect_is(simuBatch[[2]], "PhyloSim")
  
})




# 
# test_that("Neutral Model is working as expected",{
#   testthat::skip_on_cran("Too slow")
#   
#   sacVal<-matrix(0,10,100)
# 
#   richness<-numeric()
# 
#   for(i in 1:100){
#     par<-createCompletePar(x=50,y=50,dispersal = "global",
#                            density=0, environment=0, runs=5e+05, seed = i)
#     simu<-runSimulation(par)
#     sacVal[,i]<- sac(simu,area=c(4,8,16,32,64,128,256,512,1024,2048), plot=FALSE)$sr.Mean
#     richness[i]<-specRich(simu)
#   }
#   
#   parRNeutral<-createCompletePar(x=50,y=50,dispersal = "global",
#                                  density=0, environment=0, type="Rneutral",runs=1000)
#   simuRNeutral<-runSimulation(parRNeutral)
#   sacR<-sac(simuRNeutral,area=c(4,8,16,32,64,128,256,512,1024,2048), plot=FALSE)$sr.Mean
#   richnessR<-specRich(simuRNeutral)
# 
#  expect_true(richnessR<=quantile(richness, probs=0.975) & richnessR>=quantile(richness,probs=0.025))
#   for(i in 1:10){
#     expect_true(sacR[i]<=quantile(sacVal[i,],probs=0.975)& sacR[i]>=quantile(sacVal[i,],probs=0.025))
#   }
#   })

