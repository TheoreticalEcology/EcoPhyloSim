test_that("multinomialDraw() draws index i with probability w_i / sum(w)", {
  set.seed(1)
  w <- c(0, 1, 2, 0, 3, 0.5, 0)
  n <- 1e5
  counts <- multinomialDrawCounts(w, n)
  expect_equal(sum(counts), n)
  expect_equal(counts[w == 0], c(0L, 0L, 0L))
  expect_gt(stats::chisq.test(counts[w > 0], p = w[w > 0] / sum(w))$p.value, 0.001)
  expect_equal(multinomialDrawCounts(5, 10), 10L)
})

test_that("runSimulation() leaves the random number state and kind of the R session unchanged", {
  oldKind <- RNGkind()
  on.exit(RNGkind(oldKind[1], oldKind[2], oldKind[3]))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(2)
  state <- .Random.seed
  a <- quietSimulation(x = 10, y = 10, runs = 20, seed = 1)
  expect_identical(.Random.seed, state)
  RNGkind("Mersenne-Twister")
  b <- quietSimulation(x = 10, y = 10, runs = 20, seed = 1)
  expect_identical(a$Output[[1]]$specMat, b$Output[[1]]$specMat) # independent of RNGkind()
})
