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

test_that("pruning removes extinct lineages and leaves the community unchanged (#140, #147)", {
  runs <- c(1000, 3000)
  sims <- lapply(c(TRUE, FALSE), function(prune)
    runSimulation(createCompletePar(x = 20, y = 20, specRate = 8, runs = runs, seed = 1,
                                    convertToBinaryTree = FALSE, prunePhylogeny = prune)))
  for (i in seq_along(runs)) {
    pruned <- sims[[1]]$Output[[i]]
    full <- sims[[2]]$Output[[i]]
    expect_identical(pruned$specMat, full$specMat)
    expect_identical(pruned$traitMat, full$traitMat)
    # the pruned tree contains only the extant species (#147), related as in the full tree
    extant <- paste0("s", unique(c(pruned$specMat)))
    expect_setequal(pruned$phylogeny$tip.label, extant)
    expect_true(ape::is.ultrametric(pruned$phylogeny))
    expect_equal(ape::cophenetic.phylo(pruned$phylogeny)[extant, extant],
                 ape::cophenetic.phylo(ape::keep.tip(full$phylogeny, extant))[extant, extant])
    # simultaneous speciation events form polytomies, not zero-length edges or single-child nodes
    for (phy in list(pruned$phylogeny, full$phylogeny)) {
      expect_false(ape::has.singles(phy))
      expect_true(all(phy$edge.length[phy$edge[, 2] > ape::Ntip(phy)] > 0))
      expect_null(phy$root.edge)
    }
  }
})

test_that("species that died out in the last generation are removed by drop.fossil() (#147)", {
  for (seed in 1:10) for (prune in c(TRUE, FALSE)) {
    # suppressWarnings: some runs end with a single species (phylogeny set to 0)
    o <- suppressWarnings(quietSimulation(x = 3, y = 3, specRate = 1, runs = 2000, seed = seed,
                         convertToBinaryTree = FALSE, prunePhylogeny = prune))$Output[[1]]
    extant <- paste0("s", unique(c(o$specMat)))
    if (length(extant) < 2) next
    expect_setequal(ape::drop.fossil(o$phylogeny)$tip.label, extant)
    if (prune) expect_setequal(o$phylogeny$tip.label, extant)
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

test_that(".dropFossil() gives the same trees as ape::drop.fossil() (#168)", {
  set.seed(1)
  for (n in c(3, 10, 50, 200)) for (poly in c(FALSE, TRUE)) {
    phy <- ape::rcoal(n)
    if (poly) phy <- ape::di2multi(phy, tol = stats::median(phy$edge.length))
    tips <- which(phy$edge[, 2] <= n)
    short <- sample(tips, n %/% 2)
    phy$edge.length[short] <- phy$edge.length[short] * stats::runif(length(short), 0, 0.9)
    expect_equal(EcoPhyloSim:::.dropFossil(phy), ape::drop.fossil(phy))
  }
  # simulated tree with extinct species
  phy <- quietSimulation(x = 20, y = 20, specRate = 8, runs = 1000, seed = 1,
                         convertToBinaryTree = FALSE, prunePhylogeny = FALSE)$Output[[1]]$phylogeny
  expect_lt(ape::Ntip(EcoPhyloSim:::.dropFossil(phy)), ape::Ntip(phy))
  expect_equal(EcoPhyloSim:::.dropFossil(phy), ape::drop.fossil(phy))
})
