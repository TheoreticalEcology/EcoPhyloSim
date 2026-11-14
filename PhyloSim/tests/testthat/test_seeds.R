test_that("simulations with the same seed return the same summary statistics (#135)", {
  skip_on_cran() # two calls of calculateSummaryStatistics(), slow
  run <- function() quietSimulation(x = 20, y = 20, specRate = 2, seed = 1,
                                    calculateSummaries = TRUE)$Output[[1]]$summaries
  set.seed(2)
  s1 <- run()
  r1 <- stats::runif(1)
  set.seed(2)
  expect_identical(stats::runif(1), r1) # runSimulation() restores R's random number stream
  expect_identical(run(), s1) # with another state of R's random number generator
})
