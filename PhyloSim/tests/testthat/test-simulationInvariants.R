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

test_that("extant species end at the present and no branch length is negative (#81)", {
  # on small grids, global dispersal used to draw the dying individual as its own parent
  for (args in list(list(), list(environment = 0.5), list(dispersal = 0.5))) {
    for (seed in 1:20) {
      par <- c(list(x = 3, y = 3, specRate = 1, runs = 2000, seed = seed, convertToBinaryTree = FALSE), args)
      out <- suppressWarnings(do.call(quietSimulation, par))$Output[[1]] # warning for results with only one species (#153)
      if (!inherits(out$phylogeny, "phylo")) next # a single species gives no tree
      expect_true(all(out$phylogeny$edge.length >= 0))
      expect_length(setdiff(paste0("s", unique(c(out$specMat))), ape::drop.fossil(out$phylogeny)$tip.label), 0)
    }
  }
})

test_that("environmental filtering correlates traits with the environment, less so with a very narrow niche", {
  corTraitEnv <- function(...) {
    out <- quietSimulation(x = 30, y = 30, runs = 1000, environment = 1, dispersal = 1, seed = 1, ...)$Output[[1]]
    cor(c(out$traitMat), c(out$envMat))
  }
  r <- corTraitEnv()
  expect_gt(r, 0.5)
  expect_lt(suppressMessages(corTraitEnv(nicheWidth = 0.001)), r - 0.3)
})
