test_that("calculateSummaryStatistics returns all statistics", {
  set.seed(1)
  simu <- quietSimulation(x = 20, y = 20, specRate = 2, runs = 200, seed = 1, negativeDD = 0.5, environment = 0.5)
  summStat <- calculateSummaryStatistics(simu, strict = FALSE)
  expect_length(summStat, 14)
  expect_false(anyNA(unlist(summStat)))
  expect_equal(summStat$gammaDiversity, specRich(simu))
})

test_that("calculateSummaryStatistics uses the last result if runs has several values", {
  simu <- quietSimulation(x = 20, y = 20, specRate = 2, runs = c(100, 200), seed = 1)
  expect_equal(calculateSummaryStatistics(simu, strict = FALSE)$gammaDiversity, specRich(simu))
  expect_error(calculateSummaryStatistics(simu, which.result = 3), "which.result")
})

test_that("clarkEvans() places cells correctly on wide, tall and square grids (#163)", {
  set.seed(1)
  for (d in list(c(20, 50), c(50, 20), c(30, 30))) {
    mat <- matrix(sample(5, prod(d), replace = TRUE), nrow = d[1], ncol = d[2])
    expected <- mean(sapply(1:5, function(s) {
      cells <- which(mat == s, arr.ind = TRUE)
      pp <- spatstat.geom::ppp(x = cells[, "col"], y = cells[, "row"],
                               window = spatstat.geom::owin(c(0, d[2]), c(0, d[1])))
      spatstat.explore::clarkevans(pp, correction = "Donnelly")
    }))
    expect_no_warning(ce <- EcoPhyloSim:::clarkEvans(mat, correction = "Donnelly"))
    expect_equal(ce, expected)
  }
})

test_that("envSorting() is defined under strong environmental filtering (#126)", {
  data(simu.envcom.local, simu.neutral.global, envir = environment())
  expect_gt(EcoPhyloSim:::envSorting(simu.envcom.local$Output[[1]]), 0.9)
  expect_lt(EcoPhyloSim:::envSorting(simu.neutral.global$Output[[1]]), 0.05)
  # all cells count, not only the first row of the environment
  expect_equal(EcoPhyloSim:::envSorting(list(envMat = matrix(c(0, 1), 2, 2), specMat = matrix(c(1, 2), 2, 2))), 1)
  expect_identical(EcoPhyloSim:::envSorting(list(envMat = matrix(0.5, 2, 2), specMat = matrix(1:4, 2, 2))), NA_real_)
})

test_that("clarkEvans() returns NA if no species has a finite index (#126)", {
  # each species occupies a single cell
  expect_identical(suppressWarnings(EcoPhyloSim:::clarkEvans(matrix(1:4, 2, 2), correction = "Donnelly")),
                   NA_real_)
})
