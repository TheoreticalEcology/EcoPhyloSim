test_that("calculateSummaryStatistics returns all statistics", {
  set.seed(1)
  simu <- quietSimulation(x = 20, y = 20, specRate = 2, runs = 200, seed = 1, density = 0.5, environment = 0.5)
  summStat <- calculateSummaryStatistics(simu, strict = FALSE)
  expect_length(summStat, 16)
  expect_false(anyNA(unlist(summStat)))
  expect_equal(summStat$gammaDiversity, specRich(simu))
})

test_that("calculateSummaryStatistics uses the last result if runs has several values", {
  simu <- quietSimulation(x = 20, y = 20, specRate = 2, runs = c(100, 200), seed = 1)
  expect_equal(calculateSummaryStatistics(simu, strict = FALSE)$gammaDiversity, specRich(simu))
  expect_error(calculateSummaryStatistics(simu, which.result = 3), "which.result")
})

test_that("the gamma statistic resolves polytomies (#169)", {
  skip_on_cran() # calculateSummaryStatistics(), slow
  simu <- quietSimulation(x = 30, y = 30, specRate = 20, runs = 200, seed = 1,
                          convertToBinaryTree = FALSE)
  phy <- simu$Output[[1]]$phylogeny
  expect_false(ape::is.binary(phy))
  gamma <- ape::gammaStat(ape::multi2di(phy, random = FALSE))
  expect_equal(calculateSummaryStatistics(simu, strict = FALSE)$gammaStatistics, gamma)
  expect_output(summary(simu), paste("Gamma: ", format(gamma)), fixed = TRUE)
})
