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
  expect_error(NeutralMod(x = 1, y = 1, verbose = FALSE), "at least 2 cells")
})

test_that("NeutralMod without speciation keeps a single species", {
  simu <- NeutralMod(x = 5, y = 5, runs = 3, specRate = 0, seed = 1, verbose = FALSE)
  expect_equal(specRich(simu), 1)
})

test_that("NeutralMod runs with local dispersal and density dependence", {
  simu <- NeutralMod(x = 10, y = 10, runs = 5, specRate = 2, seed = 1, verbose = FALSE,
                     dispersal = function(n) stats::rexp(n, 1/2), negativeDD = 0.5, negativeDDCut = 2)
  expect_s3_class(simu, "PhyloSim")
  expect_error(NeutralMod(dispersal = 2, verbose = FALSE), "dispersal")
  expect_error(NeutralMod(x = 5, y = 5, runs = 1, dispersal = function(n) rep(0.1, n),
                          verbose = FALSE), "dead cell")
  expect_error(NeutralMod(x = 5, y = 5, runs = 1, dispersal = function(n) NA,
                          verbose = FALSE), "dispersal\\(1\\)")
  expect_error(NeutralMod(x = 5, y = 5, negativeDD = 0.5, negativeDDCut = 3, verbose = FALSE),
               "negativeDDCut")
  expect_error(NeutralMod(x = 5, y = 5, positiveDD = 0.5, positiveDDCut = 3, verbose = FALSE),
               "positiveDDCut")
})

test_that("NeutralMod draws the parent directly after 1000 rejected recruits", {
  # with one species and negativeDD = 1, every recruit is rejected
  simu <- NeutralMod(x = 5, y = 5, runs = 2, specRate = 0, negativeDD = 1, seed = 1, verbose = FALSE)
  expect_equal(specRich(simu), 1)
})

test_that("NeutralMod: density dependence increases species richness", {
  skip_on_cran()
  richness <- function(density) mean(sapply(1:10, function(s) specRich(NeutralMod(
    x = 10, y = 10, runs = 100, specRate = 2, negativeDD = density, seed = s, verbose = FALSE))))
  expect_gt(richness(1), richness(0))
})

test_that("NeutralMod restores the random number stream of the session", {
  set.seed(123)
  NeutralMod(x = 5, y = 5, runs = 1, seed = 7, verbose = FALSE)
  after <- stats::runif(1)
  set.seed(123)
  expect_equal(after, stats::runif(1))
})

test_that("runSimulation no longer accepts type = 'Rneutral'", {
  expect_error(runSimulation(createCompletePar(type = "Rneutral")), "use NeutralMod")
  expect_error(runSimulation(createCompletePar(type = "Leipzig")), "unknown model type")
})

test_that("NeutralMod richness matches the Ewens sampling formula", {
  skip_on_cran()
  # mean and variance of the number of species among J = 6 individuals, with
  # theta = specRate * (J - 1) / (J - specRate) as documented in NeutralMod(); on this
  # small grid, 500 runs separate it from theta = specRate / (1 - specRate / J), which
  # holds if the dead cell can be its own parent (#81)
  specRate <- 1.5
  J <- 6
  theta <- specRate * (J - 1) / (J - specRate)
  i <- 0:(J - 1)
  expected <- sum(theta / (theta + i))
  se <- sqrt(sum(theta * i / (theta + i)^2) / 500)
  richness <- sapply(1:500, function(s) specRich(NeutralMod(x = 3, y = 2, runs = 60,
                                    specRate = specRate, seed = s, verbose = FALSE)))
  expect_lt(abs(mean(richness) - expected), 3 * se)
})

test_that("NeutralMod: positive density dependence reduces species richness", {
  skip_on_cran()
  richness <- function(positiveDD) mean(sapply(1:5, function(s) specRich(NeutralMod(
    x = 10, y = 10, runs = 100, specRate = 2, positiveDD = positiveDD, seed = s, verbose = FALSE))))
  expect_lt(richness(1), richness(0))
})
