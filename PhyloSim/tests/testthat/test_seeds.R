test_that("simulations with the same seed return the same summary statistics (#135)", {
  skip_on_cran() # two calls of calculateSummaryStatistics(), slow
  run <- function() quietSimulation(x = 20, y = 20, specRate = 2, seed = 1,
                                    calculateSummaries = TRUE)$Output[[1]]$summaries
  s1 <- run()
  set.seed(2)
  expect_identical(run(), s1) # with another state of R's random number generator
})
