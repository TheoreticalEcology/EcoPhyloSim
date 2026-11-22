# Numeric regression test: results of small simulations are compared with the snapshots in
# _snaps/reproducibility.md (for intended changes, see "Tests and checks" in CONTRIBUTING.md).

local_edition(3)

snapshotSimulation <- function(...) {
  simu <- quietSimulation(x = 30, y = 30, runs = 200, seed = 1, ...)
  out <- simu$Output[[1]]
  set.seed(1) # random plots in sac()
  list(abundances = sort(as.vector(table(out$specMat)), decreasing = TRUE),
       nTips = ape::Ntip(out$phylogeny),
       sac = sac(simu, size = c(2, 5, 10, 20), plot = FALSE)$sr.Mean)
}

test_that("simulation results are numerically reproducible", {
  skip_on_cran()
  expect_snapshot_value(snapshotSimulation(specRate = 2), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, environment = 0.5, negativeDD = 0.5),
                        style = "deparse", tolerance = 1e-6)
})

test_that("results with only one fitness mechanism are numerically reproducible", {
  skip_on_cran()
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, environment = 1),
                        style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, negativeDD = 1, negativeDDCut = 2),
                        style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, negativeDD = 1, fitnessActsOn = "both"),
                        style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, redQueen = 0.1, redQueenStrength = 1),
                        style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, positiveDD = 0.05),
                        style = "deparse", tolerance = 1e-6)
})

# traits and the speciation paths fission and protracted, which the snapshots above do not cover;
# with the default speciationTraitJump = 0 they must not change (#114)
test_that("traits with fission and protracted speciation are numerically reproducible", {
  skip_on_cran()
  snapshotTraits <- function(...) {
    out <- quietSimulation(x = 30, y = 30, runs = 200, seed = 1, specRate = 2, dispersal = 1, ...)$Output[[1]]
    list(abundances = sort(as.vector(table(out$specMat)), decreasing = TRUE),
         traitSums = c(sum(out$traitMat), sum(out$compMat), sum(out$neutMat)))
  }
  expect_snapshot_value(snapshotTraits(environment = 0.5), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotTraits(fission = 1), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotTraits(fission = 2), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotTraits(protracted = 1), style = "deparse", tolerance = 1e-6)
})
