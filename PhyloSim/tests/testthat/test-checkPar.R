smallPar <- function(...) {
  do.call(createCompletePar, utils::modifyList(list(x = 10, y = 10, runs = 2, seed = 1), list(...)))
}

test_that("valid settings pass", {
  expect_silent(smallPar())
  expect_silent(smallPar(airmat = matrix(0.5, 10, 10), dispersal = 1, negativeDD = 1,
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
  expect_error(smallPar(negativeDD = 1.2), "negativeDD") # #96
  expect_error(smallPar(environment = -0.1), "environment")
  expect_error(smallPar(nicheWidth = 0), "'nicheWidth'")
  expect_error(smallPar(fitnessActsOn = "birth"), "fitnessActsOn")
  expect_error(smallPar(fitnessBaseMortalityRatio = 0), "fitnessBaseMortalityRatio") # #94
  expect_error(smallPar(negativeDDCut = 0), "negativeDDCut")
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
  expect_silent(smallPar(specRate = 100, dispersal = 2.4, negativeDD = 0.5, negativeDDCut = 4))
  expect_silent(smallPar(negativeDDCut = 10))
  expect_error(smallPar(negativeDDCut = 11), "negativeDDCut")
  expect_error(smallPar(x = 1, y = 1), "at least 2 cells") # #81
  expect_message(smallPar(nicheWidth = 5, environment = 0.5), "outside")
  expect_message(smallPar(nicheWidth = 0.001, environment = 0.5), "outside")
  expect_error(smallPar(specRate = 101), "specRate") # #150
  # kernels must fit into the grid without wrapping onto the same cells (#152)
  expect_error(smallPar(dispersal = 2.5), "dispersal")
  expect_error(smallPar(negativeDD = 0.5, negativeDDCut = 5), "negativeDDCut")
})

test_that("fitnessActsOn = 'reproduction' is rejected only where it is not implemented (#98)", {
  for (f in c("reproduction", "both")) {
    expect_error(smallPar(fitnessActsOn = f, environment = 0.5), "not implemented")
    expect_error(smallPar(fitnessActsOn = f, negativeDD = 0.5), "not implemented")
    expect_error(smallPar(fitnessActsOn = f, redQueen = 0.5), "not implemented")
    expect_silent(smallPar(fitnessActsOn = f))
    p <- smallPar(fitnessActsOn = f, dispersal = 1, environment = 0.5, negativeDD = 0.5)
    expect_s3_class(suppressWarnings(runSimulation(p)), "PhyloSim")
  }
})

test_that("density dependence options and fitnessActsOn stages are checked (#144, #165)", {
  expect_error(smallPar(negativeDDKernel = "trait"), "negativeDDKernel")
  expect_error(smallPar(negativeDDWidth = 0), "negativeDDWidth")
  expect_error(smallPar(negativeDDResponse = "exponential"), "negativeDDResponse")
  expect_error(smallPar(negativeDD = 0.5, negativeDDResponse = function(S) 2 - S), "between 0 and 1")
  expect_message(smallPar(negativeDD = 0.5, negativeDDResponse = function(S) S), "increases")
  expect_message(smallPar(negativeDDResponse = function(S) 1 - S), "no effect")
  expect_message(smallPar(negativeDD = 0.5, negativeDDKernel = "gaussian", negativeDDWidth = 2), "negativeDDWidth")
  expect_error(smallPar(fitnessActsOn = "recruitment", redQueen = 0.5), "no effect")
  expect_error(smallPar(fitnessActsOn = character(0)), "fitnessActsOn")
  expect_silent(smallPar(fitnessActsOn = c("mortality", "recruitment"), negativeDD = 0.5, environment = 0.5))
  expect_silent(smallPar(fitnessActsOn = "recruitment", dispersal = "global", negativeDD = 0.5))
  # parameter lists from earlier versions get the default options
  p <- smallPar()
  p$negativeDDKernel <- p$negativeDDWidth <- p$negativeDDResponse <- p$positiveDD <- NULL
  p$density <- p$negativeDD; p$negativeDD <- NULL
  expect_equal(EcoPhyloSim:::checkPar(p)[names(smallPar())], smallPar())
  expect_null(EcoPhyloSim:::checkPar(p)$density)
})

test_that("positive density dependence is checked, and old argument names still work", {
  expect_error(smallPar(positiveDD = 1.5), "positiveDD")
  expect_error(smallPar(positiveDDKernel = "trait"), "positiveDDKernel")
  expect_error(smallPar(positiveDD = 0.5, positiveDDCut = 5), "positiveDDCut")
  expect_error(smallPar(positiveDD = 0.5, fitnessActsOn = "reproduction"), "not implemented")
  expect_message(smallPar(negativeDD = 0.5, positiveDD = 0.5, positiveDDKernel = "absolute"), "offset")
  expect_silent(smallPar(negativeDD = 0.5, negativeDDCut = 2, positiveDD = 0.5))
  expect_warning(p <- smallPar(density = 0.5, densityCut = 2), "deprecated")
  expect_equal(p, smallPar(negativeDD = 0.5, negativeDDCut = 2))
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
  p2 <- p; p2$negativeDD <- 2;      expect_error(runSimulation(p2), "negativeDD")
})

test_that("callModel() stops on invalid input (#13)", {
  args <- list(x = 10, y = 10, dispersal = 3, runs = 2, specRate = 1, dens = TRUE, env = FALSE,
               neutral = FALSE, mort = TRUE, mortStrength = 10, repro = FALSE, dispersalCutoff = 2,
               densityCutoff = 1, seed = 1, envStrength = 1, compStrength = 0.5, fission = 0,
               redQueen = 0, redQueenStrength = 0, protracted = 0, airmatR = 1, soilmatR = 1,
               prunePhylogeny = TRUE, nicheWidth = 0.03659906)
  call <- function(...) utils::capture.output(do.call(callModel, utils::modifyList(args, list(...))))
  expect_error(call(), NA)
  expect_error(call(x = 0), "x and y")
  expect_error(call(x = 1, y = 1), "x and y") # #81
  expect_error(call(dispersal = 0), "dispersal")
  expect_error(call(dispersalCutoff = 0), "dispersalCutoff")
  expect_error(call(dispersalCutoff = 4, densityCutoff = 4, specRate = 100), NA)
  expect_error(call(dispersalCutoff = 5), "dispersalCutoff")
  expect_error(call(densityCutoff = 5), "densityCutoff")
  expect_error(call(densityCutoff = 0), "densityCutoff")
  expect_error(call(dens = FALSE, densityCutoff = 10), NA)
  expect_error(call(dens = FALSE, densityCutoff = 11), "densityCutoff")
  expect_error(call(mortStrength = 0), "mortStrength")
  expect_error(call(specRate = NaN), "specRate")
  expect_error(call(nicheWidth = 0), "nicheWidth must be")
  expect_error(call(specRate = 101), "specRate")
  expect_error(call(runs = integer(0)), "runs")
  expect_error(call(runs = c(5, 2)), "runs")
  expect_error(call(airmatR = numeric(0)), "airmatR")
  expect_error(call(airmatR = rep(0.5, 7)), "airmatR")
  expect_error(call(airmatR = rep(1.5, 100)), "airmatR")
  expect_error(call(soilmatR = rep(0.5, 100)), "soilmatR")
  expect_error(call(densityKernel = 4), "densityKernel")
  expect_error(call(densityWidth = 0), "densityWidth")
  expect_error(call(positiveStrength = 2), "positiveStrength")
  expect_error(call(positiveStrength = 0.5, positiveCutoff = 5), "positiveCutoff")
  expect_error(call(positiveKernel = 0), "positiveKernel")
  expect_error(call(densityTable = 0.5), "densityTable")
  expect_error(call(densityTable = c(0, 2)), "densityTable")
  expect_error(call(recruit = TRUE, mort = FALSE, densityKernel = 3, densityTable = c(1, 0.5, 0.3, 0.2, 0)), NA)
})
