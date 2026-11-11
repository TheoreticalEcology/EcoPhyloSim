# PhyloSim 0.3.1

## New features

* New pipeline for large simulation experiments: `runPhylosimAnalysis()` runs many simulations with fixed and prior-drawn parameters and saves them in chunks; `chunkwiseCalculateSummaries()`, `chunkwiseCleanup()`, `combineChunks()`, `prepareForAnalysis()` and `extractParams()` process the results (#85).
* New functions `rejectionSample()` and `rejectionSamplePhylosimAnalysis()` for ABC rejection sampling on simulation outputs (#85).
* New functions `getSpecies()`, `getSpeciesTrait()` and `getSpeciesTraits()` extract the extant species and their traits from a simulation (#85).
* `calculateSummaryStatistics()` returns additional statistics (Clark-Evans index, 8-neighbor index, environmental response coefficients) (#85).

## Changes

* `calculateSummaryStatistics()` returns `NA` if any statistic cannot be calculated; set `strict = FALSE` to get the previous behavior (#85).
* `createCompletePar()` and `runSimulation()` stop with an error if `environment` or the values of `airmat` are outside [0, 1], or if `fitnessBaseMortalityRatio < 1` (#85).
* `createCompletePar()` stops with an error if `fitnessActsOn` is `"reproduction"` or `"both"`, because these settings crashed R (#98, #101).
* `runSimulationBatch()`: errors in single simulations are no longer silently turned into `NULL` during parallel runs, and `backup = TRUE` saves each result as `<scenario>.rds` (#106).
* Default plot sizes in `sac()` are log-spaced, as documented (#115).
* The archived packages kulife, apTreeshape and PhyloMeasures are no longer needed: the Colless index is computed internally, and the phylogenetic dispersion uses `picante::ses.mpd()` with 999 randomizations, so the `dispersion` statistic of `calculateSummaryStatistics()` is stochastic and slower than before. The `useApTreeshape` argument of `collessImbalance()` is deprecated; `useApTreeshape = FALSE` no longer returns the normalized index (#89).
* spatstat.geom, spatstat.explore and vegan were added to Imports; rmarkdown, testthat, abctools, ranger, data.table and gtools to Suggests (the last four are only needed for `runPhylosimAnalysis()`, its helper functions and the rejection sampling functions) (#113).
* License changed from GPL (>= 2) to GPL (>= 3).

## Bug fixes

* Speciation by fission (`fission = 1` or `2`) assigned the wrong individuals to the new species (#85).
* `runSimulationBatch()`: `strip = "summaries"` was ignored without parallelization, and the results were not named after their scenarios (#106).
* Installation from GitHub failed because `src/RcppExports.cpp` was missing (#102).
* `createCompletePar()` and `runSimulation()`: a matrix passed as `airmat` was assigned to the wrong cells (and memory outside of it was read) if the grid is not square (`x != y`) (#122).

## Documentation

* The vignette was converted from Sweave (PDF) to R Markdown (HTML).
* Package help converted to roxygen; examples moved to `inst/examples/` (#87).
* The help of `createCompletePar()` explains vectors in `runs`, the valid ranges of `environment`, `airmat` and `fitnessBaseMortalityRatio` (#85), and the dimensions of `airmat` (#122).
* Help pages and vignette revised: corrected trait-evolution formula in the vignette and return value of `sac()`, documented parameters of `callModel()`, and faster examples (#113).
