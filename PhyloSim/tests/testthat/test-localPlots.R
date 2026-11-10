test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("localPlots returns the requested number of plots", {
  simu <- runSimulation(createCompletePar(x = 30, y = 30, runs = 100, seed = 1))
  
  for (n in c(1, 2)) {
    out <- localPlots(simu, size = 5, n = n, community = TRUE)
    
    # every output should contain exactly n plots
    expect_length(out$subPlots, n)
    expect_length(out$envPlots, n)
    expect_length(out$compPlots, n)
    expect_length(out$positions, n)
    expect_equal(nrow(out$communityTable), n)
    
    # each plot has 5 x 5 = 25 cells, so species counts must add up to 25
    expect_true(all(rowSums(out$communityTable) == 25))
  }
})

test_that("localPlots works without a community table", {
  simu <- runSimulation(createCompletePar(x = 30, y = 30, runs = 100, seed = 1))
  
  out <- localPlots(simu, size = 5, n = 1, community = FALSE)
  expect_length(out$subPlots, 1)
  expect_null(out$communityTable)
})

test_that("localPlots rejects invalid n", {
  simu <- runSimulation(createCompletePar(x = 30, y = 30, runs = 100, seed = 1))
  
  expect_error(localPlots(simu, size = 5, n = 0))
  expect_error(localPlots(simu, size = 5, n = 1.5))
})