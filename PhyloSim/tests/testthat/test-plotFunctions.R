simu <- quietSimulation(x = 20, y = 20, dispersal = 1, runs = 200, specRate = 2, seed = 1)

test_that("plotTraitDistribution() does not change the random number stream", {
  pdf(NULL)
  on.exit(dev.off())
  set.seed(1)
  expected <- runif(1)
  set.seed(1)
  plotTraitDistribution(simu)
  expect_equal(runif(1), expected)
})

test_that("speciesColors() gives one distinct color per extant species", {
  cols <- speciesColors(ape::drop.fossil(simu$Output[[1]]$phylogeny))
  expect_setequal(names(cols), paste0("s", unique(c(simu$Output[[1]]$specMat))))
  expect_false(anyDuplicated(cols) > 0)
})

test_that("plotPhylogeneticDispersion() fills all squares, also for NA p-values and p = 1 (#195)", {
  pdf(NULL)
  on.exit(dev.off())
  cols <- list()
  local_mocked_bindings(rect = function(..., col = NULL) cols[[length(cols) + 1]] <<- col, .package = "EcoPhyloSim")
  pvalues <- list(list(list("2" = c(NA, NA), "3" = c(0.2, 0.4)), # one plot size with only NA
                       list("2" = c(1, 1), "3" = c(1, 1)),       # p = 1
                       list("2" = c(NA, NA), "3" = c(NA, NA))))  # only NA
  plotPhylogeneticDispersion(pvalues, positions = list(x = 1:3, y = 1, xname = 1:3, yname = 1, zname = 1:3))
  squares <- unlist(cols[lengths(cols) == 1])
  expect_length(squares, 3)
  expect_false(anyNA(squares))
  expect_equal(squares[3], "grey")
})
