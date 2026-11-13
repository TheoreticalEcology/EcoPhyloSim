context("run Simulation")

test_that("Simulation is running correctly",{
  par<-createCompletePar(x=50, y=50, dispersal="global",density=1)
  simu<-runSimulation(par)
  
  expect_output(str(par), "List of 22")
  expect_is(simu, "PhyloSim")
  
  })

test_that("SimulationBatch is running correctly",{
  par1<-createCompletePar(x=50, y=50, dispersal="global",density=1)
  par2<-createCompletePar(x=50, y=50, dispersal="global",density=2)
  parlist<-list(par1,par2)
  simuBatch<-runSimulationBatch(parlist)

  expect_is(simuBatch, "PhylosimList")
  expect_is(simuBatch[[1]], "PhyloSim")
  expect_is(simuBatch[[2]], "PhyloSim")
  
})
