# Properties of small simulations that hold on every platform. The C++ random number
# generator is implementation-defined, so exact values are tested only in test-reproducibility.R.

test_that("output has one result per recorded generation and grids with y rows and x columns", {
  simu <- quietSimulation(x = 30, y = 20, runs = c(50, 100), seed = 1)
  expect_s3_class(simu, "PhyloSim")
  expect_length(simu$Output, 2)
  for (out in simu$Output) {
    for (mat in out[c("specMat", "traitMat", "envMat", "compMat")]) expect_equal(dim(mat), c(20, 30))
  }
  expect_equal(sum(rac(simu, plot = FALSE)$Abundance), 30 * 20)
})

test_that("the seed controls the species matrix", {
  run <- function(seed) quietSimulation(x = 20, y = 20, runs = 200, specRate = 2, seed = seed)$Output[[1]]$specMat
  expect_identical(run(1), run(1))
  expect_false(identical(run(1), run(2)))
})

test_that("the seed alone controls the phylogeny", {
  skip("the seed does not control the R-side randomness yet, see #135")
  run <- function() quietSimulation(x = 20, y = 20, runs = 300, specRate = 10, seed = 1)$Output[[1]]$phylogeny # with polytomies
  expect_equal(run(), run())
})

test_that("environmental filtering correlates traits with the environment", {
  simu <- quietSimulation(x = 30, y = 30, runs = 1000, environment = 1, dispersal = 1, seed = 1)
  expect_gt(cor(c(simu$Output[[1]]$traitMat), c(simu$Output[[1]]$envMat)), 0.5)
})
