# PhyloSim 0.3.1

## New features

* New pipeline for large simulation experiments: `runPhylosimAnalysis()` runs many simulations with fixed and prior-drawn parameters and saves them in chunks; `chunkwiseCalculateSummaries()`, `chunkwiseCleanup()`, `combineChunks()`, `prepareForAnalysis()` and `extractParams()` process the results (#85).
* New functions `rejectionSample()` and `rejectionSamplePhylosimAnalysis()` for ABC rejection sampling on simulation outputs (#85).
* New functions `getSpecies()`, `getSpeciesTrait()` and `getSpeciesTraits()` extract the extant species and their traits from a simulation (#85).
* New function `getMortalityData()` returns, for each individual and census interval, its species, number of conspecific neighbors and whether it died before the next recorded generation, for analyses of conspecific negative density dependence (@AndyEnvFr, #100).
* `calculateSummaryStatistics()` returns additional statistics (Clark-Evans index, 8-neighbor index, environmental response coefficients) (#85).
* `createCompletePar()`: new argument `nicheWidth` sets the width (standard deviation) of the environmental niche, which was fixed before (@AndyEnvFr, #103, #138).

## Changes

* `calculateSummaryStatistics()` returns `NA` if any statistic cannot be calculated; set `strict = FALSE` to get the previous behavior (#85).
* `createCompletePar()` and `runSimulation()` stop with an error if `fitnessActsOn` is `"reproduction"` or `"both"` with global dispersal and `environment`, `density`, `redQueen` or `redQueenStrength` > 0, because this is not implemented (#98, #101).
* `createCompletePar()` and `runSimulation()` check their parameters with checkmate, also in parameter lists changed by hand, as described in the help of `createCompletePar()`; newly rejected are, among others, `density` > 1 (#96), decreasing `runs`, negative `redQueen` and `redQueenStrength`, and non-integer values of integer parameters (#85, #118).
* The `Model` element of the output of `runSimulation()` keeps the numeric values of `density` and `environment` instead of `TRUE`/`FALSE` (#91, #118).
* `runSimulationBatch()`: errors in single simulations are no longer silently turned into `NULL` during parallel runs, and `backup = TRUE` saves each result as `<scenario>.rds` (#106).
* Default plot sizes in `sac()` are log-spaced, as documented (#115).
* `localPlots()` stops with an error if `n` is not a single positive whole number (checkmate was added to Imports for such checks), and the species columns of its community table are sorted by species number and contain integer counts (#118, #119).
* The archived packages kulife, apTreeshape and PhyloMeasures are no longer needed: the Colless index is computed internally, and the phylogenetic dispersion uses `picante::ses.mpd()` with 999 randomizations, so the `dispersion` statistic of `calculateSummaryStatistics()` is stochastic and slower than before. The `useApTreeshape` argument of `collessImbalance()` is deprecated; `useApTreeshape = FALSE` no longer returns the normalized index (@melina-leite, #89).
* spatstat.geom, spatstat.explore and vegan were added to Imports; rmarkdown, testthat, abctools, ranger, data.table and gtools to Suggests (the last four are only needed for `runPhylosimAnalysis()`, its helper functions and the rejection sampling functions) (#113).
* License changed from GPL (>= 2) to GPL (>= 3).
* PhyloSim now requires R (>= 4.1.0), as required by its dependencies vegan and MESS (#113).
* The `plot()`, `print()` and `summary()` methods for `PhyloSim` and `PhylosimList` objects use the argument names of the generics (`x`, `object`), so calls such as `plot(simu = s)` no longer work (#113).
* `collessImbalance()` and `summary()` use only extant species by default (`dropFossils = TRUE`), as `calculateSummaryStatistics()` already did (#117).
* The non-functional Leipzig model (`type = "Leipzig"`) and its C++ code were removed, and `runSimulation()` no longer accepts `type = "Rneutral"`. Unknown values of `type` now give an error instead of returning `NULL` (#112).
* `NeutralMod()` was rewritten: argument names as in `createCompletePar()` (`xdim`/`ydim` renamed to `x`/`y`), new defaults (`x = 50`, `y = 50`, `runs = 100`, `specRate = 1`), local dispersal, density dependence, speciation at single births, and a `PhyloSim` object as output; the dead cell is never its own parent, as in `runSimulation()` (#112, #81).
* `plotSpatialPhylo()` and `plotTraitDistribution()` color the species with the colorblind-safe viridis palette, ordered along the phylogeny, so that a species has the same color in both plots (#156).

## Bug fixes

* Speciation by fission (`fission = 1` or `2`) assigned the wrong individuals to the new species (#85).
* `runSimulationBatch()`: `strip = "summaries"` was ignored without parallelization, and the results were not named after their scenarios (#106).
* Installation from GitHub failed because `src/RcppExports.cpp` was missing (#102).
* `createCompletePar()` and `runSimulation()`: a matrix passed as `airmat` was assigned to the wrong cells (and memory outside of it was read) if the grid is not square (`x != y`) (#122).
* `calculateSummaryStatistics()`: `varPart_2` returned the same value as `varPart_1` instead of the fraction explained by competition (#113).
* `rac()`: a vector in `which.result` gave an error, and the bar plot had no rank labels (#113).
* `sac()` with a single result (including the default, the last result) computed the curve from the first result; this also affected the `sacAuc` statistic of `calculateSummaryStatistics()` (#113).
* `summary()` of a `PhylosimList` printed the column headers for the number of species and the number of runs in the wrong order (#113).
* `getAverageTraits()` did not accept the output of `runSimulation()` (#113); the bundled data sets are still affected (#128).
* `localPlots()` with `n = 1` returned two subplots (#116).
* `runSimulation()` read memory outside of `airmat`, or crashed R, if its length did not match the grid; `callModel()` stops with an error instead of crashing R for invalid input such as `mortStrength < 1` (#13, #94, #118).
* `runSimulation()`: the phylogeny was cut off, so that reading it failed, with `prunePhylogeny = FALSE`, for early values in `runs`, or with `fission`; with `prunePhylogeny = FALSE` and several values in `runs`, R crashed (#82, #97, #105).
* `runSimulation()`: with `prunePhylogeny = TRUE`, the returned phylogeny contained extinct lineages without living descendants, and with either setting, each value in `runs` leaked a copy of all species ever created, so that memory use grew with the run length and the number of values in `runs`. Such lineages are now removed during the run; the simulated communities are unchanged. With `prunePhylogeny = FALSE`, the branches of extant species in intermediate results no longer extend to the last value in `runs` (#140).
* `runSimulation()` with global dispersal: the dying individual could be its own parent, so that living species could get an extinction date (removed by `ape::drop.fossil()`, negative branch lengths, wrong species trait means); results with global dispersal change for a given seed. Grids must have at least 2 cells (#81).
* `plotTraitDistribution()` changed the random number stream of the user, and its histograms were cut off above log(count) = 7 (#156).
* `runSimulation()` crashed R with a stack overflow on large grids with global dispersal and `environment`, `density`, `redQueen` or `redQueenStrength` > 0, or with large local dispersal kernels (#149).

## Documentation

* The vignette was converted from Sweave (PDF) to R Markdown (HTML).
* Package help converted to roxygen; examples moved to `inst/examples/` (#87).
* The help of `createCompletePar()` explains vectors in `runs`, the valid ranges of `environment`, `airmat` and `fitnessBaseMortalityRatio` (#85), `density` (#96), `dispersal`, `densityCut`, `redQueen` and `redQueenStrength` (#118), and the dimensions of `airmat` (#122).
* Help pages and vignette revised: corrected trait-evolution formula in the vignette and return value of `sac()`, documented parameters of `callModel()`, faster examples, corrected spelling and several incorrect argument and return value descriptions (#113).
