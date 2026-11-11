test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("runSimulation works with prunePhylogeny = FALSE", {
  simu <- runSimulation(createCompletePar(x = 20, y = 20, runs = 1000, seed = 1,
                                          prunePhylogeny = FALSE))
  expect_s3_class(simu$Output[[1]]$phylogeny, "phylo")
})

test_that("the full phylogeny has at least as many tips as the pruned one", {
  full   <- runSimulation(createCompletePar(x = 20, y = 20, runs = 1000, seed = 1,
                                            prunePhylogeny = FALSE))
  pruned <- runSimulation(createCompletePar(x = 20, y = 20, runs = 1000, seed = 1,
                                            prunePhylogeny = TRUE))
  expect_gte(ape::Ntip(full$Output[[1]]$phylogeny),
             ape::Ntip(pruned$Output[[1]]$phylogeny))
})

