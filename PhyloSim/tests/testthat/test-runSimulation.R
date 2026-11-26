context("run Simulation")

test_that("Simulation is running correctly",{
  par<-createCompletePar(x=50, y=50, dispersal="global",negativeDD=1)
  simu<-runSimulation(par)
  
  expect_output(str(par), "List of 33")
  expect_is(simu, "PhyloSim")
  
  })

test_that("SimulationBatch is running correctly",{
  par1<-createCompletePar(x=50, y=50, dispersal="global",negativeDD=1)
  par2<-createCompletePar(x=50, y=50, dispersal="global",negativeDD=0.5)
  parlist<-list(par1,par2)
  simuBatch<-runSimulationBatch(parlist)

  expect_is(simuBatch, "PhylosimList")
  expect_is(simuBatch[[1]], "PhyloSim")
  expect_is(simuBatch[[2]], "PhyloSim")

})

test_that("dispersal between multiples of 0.5 has an effect (#141)", {
  f <- function(d) runSimulation(createCompletePar(x = 20, y = 20, dispersal = d, runs = 100, specRate = 2,
                                                   seed = 1))$Output[[1]]$specMat
  expect_false(identical(f(0.5), f(0.75))) # more cells within 2 * dispersal
  expect_false(identical(f(1), f(1.05)))   # same cells, other kernel scale
  expect_identical(f(1.5), f(0.3 / 0.2)) # 1.4999999999999998 gives the same cells
})

test_that("runSimulation warns if only one species is left", {
  expect_warning(simu <- quietSimulation(x = 10, y = 10, runs = 1, specRate = 0),
                 "Only one species in the community at generation 1")
  expect_identical(simu$Output[[1]]$phylogeny, 0)
})

test_that("runSimulation shows progress only on request (#221)", {
  par <- createCompletePar(x = 20, y = 20, runs = 100, specRate = 2, seed = 1)
  expect_silent(runSimulation(par)) # Rscript and testthat are not interactive
  expect_error(runSimulation(par, progress = NA), "progress")
  expect_output(callModel(x = 10, y = 10, dispersal = 1, runs = c(5, 10), specRate = 1, dens = FALSE, env = FALSE,
                          neutral = TRUE, mort = FALSE, mortStrength = 1, repro = FALSE, dispersalCutoff = 1,
                          densityCutoff = 1, envStrength = 1, compStrength = 1, fission = 0, redQueen = 0,
                          redQueenStrength = 0, protracted = 0, airmatR = 1, soilmatR = 1, prunePhylogeny = TRUE,
                          nicheWidth = 0.03, progress = TRUE),
                "100%  gen 10/10, finished after")
})

test_that("runSimulation stops early at a time limit (#221)", {
  skip_on_cran()
  par <- createCompletePar(x = 20, y = 20, runs = 100, specRate = 2, seed = 1)
  simu <- runSimulation(par)
  # the C++ run checks for interrupts and time limits, so it stops early with the time limit error
  set.seed(2)
  seed <- .Random.seed
  time <- system.time(res <- tryCatch({
    setTimeLimit(elapsed = 0.5, transient = TRUE)
    runSimulation(createCompletePar(x = 100, y = 100, runs = 20000, seed = 1))
  }, error = function(e) conditionMessage(e), finally = setTimeLimit()))[["elapsed"]]
  expect_identical(res, gettext("reached elapsed time limit", domain = "R"))
  expect_lt(time, 3)
  expect_identical(.Random.seed, seed)
  expect_identical(runSimulation(par)$Output, simu$Output)
})
