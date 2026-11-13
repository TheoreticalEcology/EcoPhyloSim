test_that("the full phylogeny has at least as many tips as the pruned one", {
  full   <- runSimulation(createCompletePar(x = 20, y = 20, runs = 1000, seed = 1,
                                            prunePhylogeny = FALSE))
  pruned <- runSimulation(createCompletePar(x = 20, y = 20, runs = 1000, seed = 1,
                                            prunePhylogeny = TRUE))
  expect_gte(ape::Ntip(full$Output[[1]]$phylogeny),
             ape::Ntip(pruned$Output[[1]]$phylogeny))
  extant <- paste0("s", unique(c(full$Output[[1]]$specMat)))
  expect_true(all(extant %in% full$Output[[1]]$phylogeny$tip.label))
})

test_that("phylogenies are complete for several runs and early time points", {
  for (prune in c(TRUE, FALSE)) {
    simu <- runSimulation(createCompletePar(x = 20, y = 20, specRate = 2, seed = 1,
                                            runs = c(3, 1000), prunePhylogeny = prune))
    for (out in simu$Output) expect_s3_class(out$phylogeny, "phylo")
  }
})

test_that("phylogenies with and without pruning agree on the extant species", {
  # #76 (pruned and full phylogeny agree); polytomies are not resolved randomly (#135)
  run <- function(prune) quietSimulation(x = 20, y = 20, runs = 1000, seed = 3, prunePhylogeny = prune,
                                         convertToBinaryTree = FALSE)$Output[[1]]
  full <- run(FALSE)
  pruned <- run(TRUE)
  expect_identical(full$specMat, pruned$specMat)
  extant <- paste0("s", unique(c(full$specMat)))
  expect_true(ape::all.equal.phylo(ape::keep.tip(full$phylogeny, extant),
                                   ape::keep.tip(pruned$phylogeny, extant)))
})
