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

test_that("pruning removes extinct lineages and leaves the community unchanged (#140)", {
  runs <- c(1000, 3000)
  sims <- lapply(c(TRUE, FALSE), function(prune)
    runSimulation(createCompletePar(x = 20, y = 20, specRate = 8, runs = runs, seed = 1,
                                    convertToBinaryTree = FALSE, prunePhylogeny = prune)))
  for (i in seq_along(runs)) {
    pruned <- sims[[1]]$Output[[i]]
    full <- sims[[2]]$Output[[i]]
    expect_identical(pruned$specMat, full$specMat)
    expect_identical(pruned$traitMat, full$traitMat)
    # the extant species are related in the same way in both trees
    extant <- paste0("s", unique(c(pruned$specMat)))
    expect_true(all(extant %in% pruned$phylogeny$tip.label))
    expect_equal(ape::cophenetic.phylo(ape::keep.tip(pruned$phylogeny, extant))[extant, extant],
                 ape::cophenetic.phylo(ape::keep.tip(full$phylogeny, extant))[extant, extant])
    # extinct species remain only as ancestors, i.e. also as internal nodes
    extinct <- setdiff(pruned$phylogeny$tip.label, extant)
    expect_true(all(extinct %in% pruned$phylogeny$node.label))
  }
})

test_that("the pruned phylogeny does not depend on the number of censuses (#140)", {
  last <- function(runs) {
    simu <- runSimulation(createCompletePar(x = 20, y = 20, specRate = 8, runs = runs, seed = 1))
    simu$Output[[length(runs)]]
  }
  one <- last(5000)
  many <- last(seq(500, 5000, 500))
  expect_identical(one$specMat, many$specMat)
  expect_identical(one$phyloTXT, many$phyloTXT)
})
