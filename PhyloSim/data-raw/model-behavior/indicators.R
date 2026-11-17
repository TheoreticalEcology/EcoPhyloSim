# Indicators of model behavior, computed with the current EcoPhyloSim on the
# output of runSimulation() of any version (list with Model$runs and Output[[i]]
# with specMat, traitMat, envMat, compMat, neutMat, phylogeny).

# share of cells at distance r (along rows and columns, torus) with the same
# species; r = 1 gives the share of the 4 neighbors of the same species
distanceDecay <- function(m, r = unique(pmin(c(1, 2, 4, 8, 16, 32), floor(min(dim(m)) / 2)))) {
  shift <- function(m, dr, dc) m[(seq_len(nrow(m)) + dr - 1) %% nrow(m) + 1,
                                 (seq_len(ncol(m)) + dc - 1) %% ncol(m) + 1]
  data.frame(distance = r, sameSpecies = sapply(r, function(d)
    mean(c(m == shift(m, d, 0), m == shift(m, 0, d)))))
}

# indicators of census i
censusStats <- function(sim, i) {
  o <- sim$Output[[i]]
  ab <- as.vector(table(o$specMat))
  p <- ab / sum(ab)
  # tree of the extant species (older versions also keep extinct species as
  # tips); multi2di() for polytomies, on which gammaStat() gives wrong values
  phylo <- tryCatch(ape::multi2di(ape::keep.tip(o$phylogeny, paste0("s", unique(c(o$specMat)))),
                                  random = FALSE), error = function(e) NULL)
  enoughTips <- !is.null(phylo) && ape::Ntip(phylo) >= 3
  # species age from the species ID: IDs are counted up at each speciation, at
  # a roughly constant rate (largest ID / generation); NA if only one species is left
  G <- sim$Model$runs[i]
  id <- as.integer(names(table(o$specMat)))
  age <- if (inherits(o$phylogeny, "phylo"))
    G * (1 - id / max(as.integer(sub("s", "", o$phylogeny$tip.label)), id)) else NA
  shannon <- -sum(p * log(p))
  data.frame(
    generation = G,
    richness = length(ab),
    shannon = shannon,
    evenness = if (length(ab) > 1) shannon / log(length(ab)) else NA,
    dominance = max(p),
    conspecificNeighbors = distanceDecay(o$specMat, 1)$sameSpecies,
    corTraitEnv = suppressWarnings(cor(c(o$traitMat), c(o$envMat))),
    sdCompetitionTrait = sd(c(o$compMat)),
    sdNeutralTrait = sd(c(o$neutMat)),
    colless = if (enoughTips) tryCatch(EcoPhyloSim::collessImbalance(
      list(Output = list(list(phylogeny = phylo))), norm = "yule", dropFossils = FALSE),
      error = function(e) NA) else NA,
    gamma = if (enoughTips) tryCatch(ape::gammaStat(phylo), error = function(e) NA) else NA,
    treeHeight = if (enoughTips) max(ape::node.depth.edgelength(phylo)) else NA,
    meanSpeciesAge = mean(age),
    corAgeAbundance = if (length(ab) >= 3) suppressWarnings(cor(age, ab, method = "spearman")) else NA)
}

# conspecific density dependence of mortality over the last consecutive
# generations: coefficient of the number of conspecifics within radius 2 in a
# logistic regression of death, with a species effect
cndd <- function(sim, last) {
  sub <- sim
  sub$Output <- sim$Output[last]
  sub$Model$runs <- sim$Model$runs[last]
  d <- EcoPhyloSim::getMortalityData(sub, radius = 2)
  d <- d[d$interval == 1, ]
  if (length(unique(d$species)) < 2) NA else
    tryCatch(coef(glm(died ~ con + factor(species), binomial, data = d))[["con"]], error = function(e) NA)
}

# all indicators of one simulation: time series at the censuses and, at the
# end, curves, CNDD and the phylogenetic dispersion of local communities
# (mean p-value of the MPD of 20 subplots of 10 x 10 cells against random
# communities drawn from the whole grid; small: clustered, large: overdispersed)
behavior <- function(sim, censuses) {
  last <- length(sim$Output)
  final <- sim$Output[[last]]
  list(census = do.call(rbind, lapply(which(sim$Model$runs %in% censuses), censusStats, sim = sim)),
       decay = distanceDecay(final$specMat),
       cndd = cndd(sim, last - 5:0),
       mpdLocal = tryCatch(mean(EcoPhyloSim::nullModel(sim, which.result = last, localPlotSize = 10,
                                                       numberOfPlots = 20, repetitions = 99), na.rm = TRUE),
                           error = function(e) NA),
       sac = EcoPhyloSim::sac(sim, which.result = last,
                              size = unique(round(2^seq(0, log2(min(dim(final$specMat))), length.out = 8))),
                              rep = 20, plot = FALSE)[c("size", "sr.Mean")],
       rac = sort(as.vector(table(final$specMat)), decreasing = TRUE),
       specMat = final$specMat)
}
