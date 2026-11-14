context("run Simulation")

test_that("Simulation is running correctly",{
  par<-createCompletePar(x=50, y=50, dispersal="global",density=1)
  simu<-runSimulation(par)
  
  expect_output(str(par), "List of 26")
  expect_is(simu, "PhyloSim")
  
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

test_that("runSimulation warns if only one species is left", {
  expect_warning(simu <- quietSimulation(x = 10, y = 10, runs = 1, specRate = 0),
                 "Only one species in the community at generation 1")
  expect_identical(simu$Output[[1]]$phylogeny, 0)
})
