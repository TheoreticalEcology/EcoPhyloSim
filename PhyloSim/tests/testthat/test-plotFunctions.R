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
