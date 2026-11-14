test_that("nullModel draws null communities of the size of the observed plots", {
  sizes <- list()
  local_mocked_bindings(mpd = function(samp, dis, abundance.weighted = FALSE) {
    sizes[[length(sizes) + 1]] <<- unname(rowSums(samp))
    rep(0.5, nrow(samp))
  }, .package = "picante")
  nullModel(simu.neutral.local, localPlotSize = 3, numberOfPlots = 2, repetitions = 5)
  expect_equal(sizes[[1]], rep(9, 2)) # observed plots
  expect_equal(sizes[[2]], rep(9, 5)) # null communities, #132
  expect_error(nullModel(simu.neutral.local, localPlotSize = 0, numberOfPlots = 2, repetitions = 5), "localPlotSize")
})
