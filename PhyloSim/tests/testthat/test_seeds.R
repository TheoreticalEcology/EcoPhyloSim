context("seeds")

testthat::test_that("simulations with the same seed return same summary statistics", {
  seeds <- c(1, 2)
  for (seed in seeds) {
    set.seed(seed)
    simu1 <- runSimulation(createCompletePar(specRate = 2, calculateSummaries = T))
    s1 <- simu1$Output[[1]]$summaries
    set.seed(seed)
    simu2 <- runSimulation(createCompletePar(specRate = 2, calculateSummaries = T))
    s2 <- simu2$Output[[1]]$summaries
    testthat::expect_equal(s1, s2)
  }
})
