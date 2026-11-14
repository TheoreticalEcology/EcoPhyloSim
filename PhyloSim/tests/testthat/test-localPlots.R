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
  
  expect_error(localPlots(simu, size = 5, n = 0), "Assertion on 'n' failed")
  expect_error(localPlots(simu, size = 5, n = 1.5), "Assertion on 'n' failed")
})

test_that("localPlots works on non-square grids", {
  # 4 rows (y) x 6 columns (x); each cell stores its row + 100 * column
  m <- outer(1:4, 1:6, function(r, c) r + 100 * c)
  fake <- fakeSimulation(m)
  set.seed(1)
  out <- localPlots(fake, size = c(1, 3, 4), n = 20)
  expect_identical(out$envPlots, out$subPlots)
  for (i in seq_along(out$subPlots)) {
    p <- out$subPlots[[i]]
    k <- c(1, 3, 4)[(i - 1) %% 3 + 1]
    expect_equal(dim(as.matrix(p)), c(k, k))
    # rows and columns of the plot are neighboring rows and columns of m (wrapped)
    rows <- as.matrix(p)[, 1] %% 100
    cols <- as.matrix(p)[1, ] %/% 100
    expect_equal(rows, (rows[1] - 1 + 0:(k - 1)) %% 4 + 1)
    expect_equal(cols, (cols[1] - 1 + 0:(k - 1)) %% 6 + 1)
    # positions are (x, y), i.e. (column, row)
    expect_equal(out$positions[[i]][1], PhyloSim:::wrap_coordinates(cols[1] + k / 2 - 0.5, 6))
    expect_equal(out$positions[[i]][2], PhyloSim:::wrap_coordinates(rows[1] + k / 2 - 0.5, 4))
  }
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(localPlots(fake, size = 3, n = 2, plot = TRUE))
  expect_error(localPlots(fake, size = 5, n = 1), "Assertion on 'size' failed")
})
