# PhyloSim 0.3.1

* `combineCunks()` was renamed to `combineChunks()`.
* The archived packages apTreeshape and PhyloMeasures are no longer needed: the Colless index is computed internally, and the phylogenetic dispersion uses `picante::ses.mpd()`. The `useApTreeshape` argument of `collessImbalance()` is deprecated; `useApTreeshape = FALSE` no longer returns the normalised index.
* The `dispersion` statistic of `calculateSummaryStatistics()` now uses 999 randomisations, so it is stochastic and slower than before.
* abctools, ranger, data.table and gtools moved from Imports to Suggests; they are only needed for `runPhylosimAnalysis()`, its helper functions and the rejection sampling functions.
