smallPar <- function(...) {
  do.call(createCompletePar, utils::modifyList(list(x = 10, y = 10, runs = 2, seed = 1), list(...)))
}

test_that("valid settings pass", {
  expect_silent(smallPar())
  expect_silent(smallPar(airmat = matrix(0.5, 10, 10), dispersal = 1, density = 1,
                         environment = 0.5, runs = c(1, 1, 2)))
  expect_silent(createCompletePar(x = 20, y = 30, airmat = matrix(0.5, 30, 20)))
})

test_that("createCompletePar() rejects invalid single parameters", {
  expect_error(smallPar(x = 0), "'x'")
  expect_error(smallPar(runs = c(5, 2, 10)), "runs")
  expect_error(smallPar(runs = -1), "runs")
  expect_error(smallPar(dispersal = 0.3), "dispersal")
  expect_error(smallPar(dispersal = "local"), "dispersal")
  expect_error(smallPar(specRate = -1), "specRate")
  expect_error(smallPar(density = 1.2), "density") # #96
  expect_error(smallPar(environment = -0.1), "environment")
  expect_error(smallPar(fitnessActsOn = "birth"), "fitnessActsOn")
  expect_error(smallPar(fitnessBaseMortalityRatio = 0), "fitnessBaseMortalityRatio") # #94
  expect_error(smallPar(densityCut = 0), "densityCut")
  expect_error(smallPar(fission = 3), "fission")
  expect_error(smallPar(protracted = -1), "protracted")
  expect_error(smallPar(redQueen = -1), "redQueen")
  expect_error(smallPar(prunePhylogeny = NA), "prunePhylogeny")
})

test_that("createCompletePar() checks the size and values of airmat", {
  expect_error(createCompletePar(x = 20, y = 30, airmat = rep(0.5, 10)), "airmat")
  expect_error(createCompletePar(x = 20, y = 30, airmat = matrix(0.5, 20, 30)), "airmat")
  expect_error(smallPar(airmat = matrix(1.5, 10, 10)), "airmat")
  expect_error(smallPar(airmat = matrix(NA_real_, 10, 10)), "airmat")
})

test_that("createCompletePar() checks combinations of parameters", {
  expect_error(smallPar(dispersal = 6), "dispersal")
  expect_error(smallPar(densityCut = 11), "densityCut")
})

test_that("fitnessActsOn = 'reproduction' is rejected only where it is not implemented (#98)", {
  for (f in c("reproduction", "both")) {
    expect_error(smallPar(fitnessActsOn = f, environment = 0.5), "not implemented")
    expect_error(smallPar(fitnessActsOn = f, density = 0.5), "not implemented")
    expect_error(smallPar(fitnessActsOn = f, redQueen = 0.5), "not implemented")
    expect_silent(smallPar(fitnessActsOn = f))
    p <- smallPar(fitnessActsOn = f, dispersal = 1, environment = 0.5, density = 0.5)
    expect_s3_class(suppressWarnings(runSimulation(p)), "PhyloSim")
  }
})

test_that("runSimulation() checks parameter lists changed after createCompletePar()", {
  p <- smallPar()
  # wrong airmat lengths read memory outside the vector, length 0 crashed R (#118)
  for (len in c(0, 10, 101)) {
    p2 <- p
    p2$airmatR <- rep(0.5, len)
    expect_error(runSimulation(p2), "airmat")
  }
  p2 <- p; p2$x <- 0;               expect_error(runSimulation(p2), "'x'")
  p2 <- p; p2$runs <- c(5, 2);      expect_error(runSimulation(p2), "runs")
  p2 <- p; p2$runs <- numeric(0);   expect_error(runSimulation(p2), "runs")
  p2 <- p; p2$fitnessBaseMortalityRatio <- 0; expect_error(runSimulation(p2), "fitnessBaseMortalityRatio")
  p2 <- p; p2$density <- 2;         expect_error(runSimulation(p2), "density")
})

test_that("callModel() stops on input that would crash R (#13)", {
  args <- list(x = 10, y = 10, dispersal = 3, runs = 2, specRate = 1, dens = TRUE, env = FALSE,
               neutral = FALSE, mort = TRUE, mortStrength = 10, repro = FALSE, dispersalCutoff = 2,
               densityCutoff = 1, seed = 1, envStrength = 1, compStrength = 0.5, fission = 0,
               redQueen = 0, redQueenStrength = 0, protracted = 0, airmatR = 1, soilmatR = 1,
               prunePhylogeny = TRUE)
  call <- function(...) utils::capture.output(do.call(callModel, utils::modifyList(args, list(...))))
  expect_error(call(), NA)
  expect_error(call(x = 0), "x and y")
  expect_error(call(dispersal = 0), "dispersal")
  expect_error(call(dispersalCutoff = 0), "dispersalCutoff")
  expect_error(call(dispersalCutoff = 11), "dispersalCutoff")
  expect_error(call(densityCutoff = 11), "densityCutoff")
  expect_error(call(mortStrength = 0), "mortStrength")
  expect_error(call(specRate = NaN), "specRate")
  expect_error(call(runs = integer(0)), "runs")
  expect_error(call(runs = c(5, 2)), "runs")
  expect_error(call(airmatR = numeric(0)), "airmatR")
  expect_error(call(airmatR = rep(0.5, 7)), "airmatR")
  expect_error(call(soilmatR = rep(0.5, 100)), "soilmatR")
})
