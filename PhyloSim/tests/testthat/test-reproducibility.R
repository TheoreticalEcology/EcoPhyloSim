# Numeric regression test: results of small simulations are compared with the snapshots in
# _snaps/reproducibility.md (for intended changes, see "Tests and checks" in CONTRIBUTING.md).
# The C++ random number generator (std::default_random_engine) is implementation-defined, so
# results on macOS (libc++) may differ; snapshots are only compared on Windows and Linux.

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
  skip_on_os("mac")
  expect_snapshot_value(snapshotSimulation(specRate = 2), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, environment = 0.5, density = 0.5),
                        style = "deparse", tolerance = 1e-6)
})
