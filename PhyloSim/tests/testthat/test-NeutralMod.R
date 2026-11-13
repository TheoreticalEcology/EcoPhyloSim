test_that("NeutralMod returns a PhyloSim object with a species matrix", {
  simu <- NeutralMod(x = 12, y = 8, runs = c(5, 10), specRate = 2, seed = 1, verbose = FALSE)
  expect_s3_class(simu, "PhyloSim")
  expect_length(simu$Output, 2)
  expect_equal(dim(simu$Output[[2]]$specMat), c(8, 12))
  expect_equal(simu$Model$type, "NeutralMod")
  expect_identical(simu$Output, NeutralMod(x = 12, y = 8, runs = c(5, 10), specRate = 2,
                                           seed = 1, verbose = FALSE)$Output)
  expect_equal(specRich(simu), length(unique(c(simu$Output[[2]]$specMat))))
  square <- NeutralMod(x = 10, y = 10, runs = 5, specRate = 2, seed = 1, verbose = FALSE)
  expect_silent(sac(square, plot = FALSE))
})

test_that("NeutralMod without speciation keeps a single species", {
  simu <- NeutralMod(x = 5, y = 5, runs = 3, specRate = 0, seed = 1, verbose = FALSE)
  expect_equal(specRich(simu), 1)
})

test_that("NeutralMod runs with local dispersal and density dependence", {
  simu <- NeutralMod(x = 10, y = 10, runs = 5, specRate = 2, seed = 1, verbose = FALSE,
                     dispersal = function(n) stats::rexp(n, 1/2), density = 0.5, densityCut = 2)
  expect_s3_class(simu, "PhyloSim")
  expect_error(NeutralMod(dispersal = 2, verbose = FALSE), "dispersal")
  expect_error(NeutralMod(x = 5, y = 5, runs = 1, dispersal = function(n) rep(0.1, n),
                          verbose = FALSE), "dead cell")
})

test_that("runSimulation no longer accepts type = 'Rneutral'", {
  expect_error(runSimulation(createCompletePar(type = "Rneutral")), "use NeutralMod")
  expect_error(runSimulation(createCompletePar(type = "Leipzig")), "unknown model type")
})

test_that("NeutralMod richness matches the Ewens sampling formula with theta = specRate", {
  skip_on_cran()
  # mean and variance of the number of species among J = 100 individuals
  theta <- 2
  i <- 0:99
  expected <- sum(theta / (theta + i))
  se <- sqrt(sum(theta * i / (theta + i)^2) / 100)
  richness <- sapply(1:100, function(s) specRich(NeutralMod(x = 10, y = 10, runs = 500,
                                    specRate = theta, seed = s, verbose = FALSE)))
  expect_lt(abs(mean(richness) - expected), 3 * se)
})
