test_that("nicheWidth changes the results", {
  a <- runSimulation(createCompletePar(x = 20, y = 20, runs = 100, environment = 0.5, seed = 1))
  b <- runSimulation(createCompletePar(x = 20, y = 20, runs = 100, environment = 0.5, seed = 1,
                                       nicheWidth = 0.3))
  expect_false(identical(a$Output[[1]]$specMat, b$Output[[1]]$specMat))
})

test_that("nicheWidth must be a positive number", {
  expect_error(createCompletePar(nicheWidth = 0), "nicheWidth")
  expect_error(createCompletePar(nicheWidth = -1), "nicheWidth")
})