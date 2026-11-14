# Behaviour of the fitness combination, density kernels, density response and fitness on
# recruitment (#144, #165). Small runs; thresholds are far from the values observed with
# several seeds (e.g. CNDD slope about 0.45 with the species kernel, about 0 without).

conspecificShare <- function(m) {
  right <- m[, c(2:ncol(m), 1)]
  down <- m[c(2:nrow(m), 1), ]
  mean(c(m == right, m == down))
}

test_that("environment and density dependence together are no longer neutral (#144)", {
  simu <- quietSimulation(x = 30, y = 30, dispersal = 1, runs = 500, specRate = 2,
                          environment = 0.5, density = 0.5, seed = 1)
  out <- simu$Output[[1]]
  expect_gt(cor(c(out$traitMat), c(out$envMat)), 0.5)
})

test_that("a density response function equal to the linear one gives identical results", {
  run <- function(...) quietSimulation(x = 20, y = 20, dispersal = 1, runs = 100, specRate = 2,
                                       density = 0.5, densityKernel = "species", seed = 1, ...)
  expect_identical(run()$Output[[1]]$specMat,
                   run(densityResponse = function(S) 1 - 0.5 * S)$Output[[1]]$specMat)
})

test_that("the species kernel gives conspecific density dependence in mortality", {
  slope <- function(density) {
    simu <- quietSimulation(x = 30, y = 30, dispersal = 1, runs = 300:305, specRate = 2,
                            density = density, densityKernel = "species", seed = 1)
    coef(glm(died ~ con, binomial(link = "cloglog"), data = getMortalityData(simu)))[2]
  }
  expect_gt(slope(1), 0.2)
  expect_lt(abs(slope(0)), 0.1)
})

test_that("density dependence on recruitment reduces conspecific neighbors, also with global dispersal", {
  share <- function(...) conspecificShare(quietSimulation(x = 30, y = 30, runs = 300, specRate = 2,
                                                          seed = 1, ...)$Output[[1]]$specMat)
  expect_lt(share(dispersal = 1, density = 1, densityKernel = "species", fitnessActsOn = "recruitment"),
            share(dispersal = 1) - 0.2)
  simu <- quietSimulation(x = 20, y = 20, dispersal = "global", runs = 100, specRate = 2, density = 1,
                          densityKernel = "species", fitnessActsOn = "recruitment", seed = 1)
  expect_s3_class(simu, "PhyloSim")
})

test_that("environment on recruitment filters traits, and all stages can be combined", {
  out <- quietSimulation(x = 20, y = 20, dispersal = 1, runs = 150, specRate = 2, environment = 1,
                         fitnessActsOn = "recruitment", seed = 1)$Output[[1]]
  expect_gt(cor(c(out$traitMat), c(out$envMat)), 0.5)
  simu <- quietSimulation(x = 20, y = 20, dispersal = 1, runs = 100, specRate = 2, environment = 1,
                          density = 1, densityKernel = "gaussian", densityWidth = 0.05,
                          fitnessActsOn = c("mortality", "reproduction", "recruitment"), seed = 1)
  expect_s3_class(simu, "PhyloSim")
})
