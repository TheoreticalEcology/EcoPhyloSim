test_that("simulations with the same seed return the same summary statistics", {
  skip_on_cran() # two calls of calculateSummaryStatistics(), slow
  run <- function() {
    set.seed(1) # R-side randomness, see #135
    quietSimulation(x = 20, y = 20, specRate = 2, seed = 1, calculateSummaries = TRUE)$Output[[1]]$summaries
  }
  expect_equal(run(), run())
})
