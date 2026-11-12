context("collessImbalance")

test_that("internal Colless index gives known values", {
  caterpillar <- ape::stree(8, type = "left")
  balanced <- ape::stree(8, type = "balanced")

  # caterpillar tree with n tips: sum of |left - right| = (n - 1)(n - 2) / 2
  expect_equal(PhyloSim:::.collessIndex(caterpillar), 21)
  expect_equal(PhyloSim:::.collessIndex(balanced), 0)

  # normalizations as in apTreeshape::colless(): yule = (21 - E[Yule]) / 8, pda = 21 / 8^1.5
  expect_equal(PhyloSim:::.collessIndex(caterpillar, norm = "yule"), 1.661490, tolerance = 1e-6)
  expect_equal(PhyloSim:::.collessIndex(caterpillar, norm = "pda"), 0.928078, tolerance = 1e-6)
  expect_error(PhyloSim:::.collessIndex(caterpillar, norm = "other"))
})
