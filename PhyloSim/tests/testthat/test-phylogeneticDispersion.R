test_that("calculatePhylogeneticDispersion accepts several results", {
  simu <- quietSimulation(x = 20, y = 20, dispersal = "global", runs = c(50, 100), specRate = 10, seed = 1)
  out <- calculatePhylogeneticDispersion(simu, plotlengths = 3, plots = 5, replicates = 5,
                                         which.result = c(1, 2)) # #131
  expect_length(out, 2)
  for (t in 1:2) expect_length(out[[t]][[1]][["3"]], 5)
})
