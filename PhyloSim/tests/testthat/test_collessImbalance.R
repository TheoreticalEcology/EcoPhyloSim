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

test_that("collessImbalance() uses only extant species by default", {
  simu <- runSimulation(createCompletePar(x = 20, y = 20, runs = 1000, seed = 1,
                                          prunePhylogeny = FALSE))
  phylo <- simu$Output[[1]]$phylogeny
  expect_gt(ape::Ntip(phylo), ape::Ntip(ape::drop.fossil(phylo)))
  expect_equal(collessImbalance(simu),
               PhyloSim:::.collessIndex(ape::drop.fossil(phylo)))
})