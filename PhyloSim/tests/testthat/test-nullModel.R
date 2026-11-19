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

test_that("nullModel gives roughly uniform p-values in a neutral model with global dispersal (#186)", {
  skip_on_cran()
  simus <- lapply(1:3, function(s) quietSimulation(x = 50, y = 50, dispersal = "global", runs = 5000,
                                                   specRate = 5, seed = s))
  for (abundance in c(FALSE, TRUE)) {
    set.seed(1)
    p <- unlist(lapply(simus, nullModel, abundance = abundance, localPlotSize = 5,
                       numberOfPlots = 100, repetitions = 500))
    expect_false(anyNA(p)) # no plot with a single species
    expect_gt(mean(p), 0.4)
    expect_lt(mean(p), 0.6)
    expect_lt(mean(p < 0.05), 0.12) # plots overlap and share one null distribution, so no strict test
  }
})

test_that("nullModel gives NA, not p = 0, for plots with a single species (#186)", {
  # same seed, so that nullModel() draws the same plots (localPlots() is its first use of random numbers)
  set.seed(1)
  single <- rowSums(localPlots(simu.neutral.global, size = 2, n = 100, community = TRUE)$communityTable > 0) == 1
  expect_true(any(single))
  set.seed(1)
  p <- nullModel(simu.neutral.global, localPlotSize = 2, numberOfPlots = 100, repetitions = 100)
  expect_equal(is.na(p), unname(single))
})

test_that("nullModel gives NA p-values if all null communities have a single species (#196)", {
  for (abundance in c(FALSE, TRUE))
    expect_equal(nullModel(simu.neutral.global, abundance = abundance, localPlotSize = 1, numberOfPlots = 5, repetitions = 20), rep(NA_real_, 5))
  p <- calculatePhylogeneticDispersion(simu.neutral.global, plotlengths = 1, plots = 5, replicates = 20)
  expect_equal(p[[1]][[1]][["1"]], rep(NA_real_, 5))
})
