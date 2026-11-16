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

test_that("nullModel(fun = 'pd') works for a phylogeny with a polytomy at the root (#176)", {
  simu <- simu.neutral.local
  phy <- simu$Output[[1]]$phylogeny
  phy$root.edge <- NULL # as in simulations since #160
  # collapse an internal edge below the root into a root polytomy
  phy$edge.length[which(phy$edge[, 1] == ape::Ntip(phy) + 1 & phy$edge[, 2] > ape::Ntip(phy))[1]] <- 0
  simu$Output[[1]]$phylogeny <- ape::di2multi(phy)
  expect_false(ape::is.rooted(simu$Output[[1]]$phylogeny))
  expect_length(nullModel(simu, localPlotSize = 3, numberOfPlots = 4, repetitions = 20, fun = "pd"), 4)
})
