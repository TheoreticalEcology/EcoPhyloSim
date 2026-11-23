# EcoPhyloSim (development version)

## New features

* `createCompletePar()`: new argument `traitSpeciationJump` (one value or one per trait, default 0) changes the traits of a new species by a random jump at speciation, shared by all its individuals (#114).
* `createCompletePar()`: new arguments `traitEvolution` and `traitPull` (one value or one per trait) set the half-width of the random trait change at each birth and the weight of the pull to the species mean; the defaults 0.01 and 0.2 are the values fixed so far (#220).
* `runSimulation()`: new argument `progress` (default `TRUE`) shows a progress bar with the estimated remaining time in interactive sessions; runs can be interrupted with Esc or Ctrl+C, and the end messages are no longer printed (#221).

## Changes

* `runSimulation()` is faster: about 5 to 60 times with local dispersal (more with larger `dispersal`; less with fitness acting on reproduction), 5 to 25 times with density dependence, about 2 times in global neutral runs and up to 7 times with global dispersal and fitness on recruitment. Results for a given seed are unchanged, except with global dispersal and the environment acting on recruitment, where parents are now proposed from bins of their environmental trait: each parent is drawn with the same probability as before, but with other random numbers (#207).

## Bug fixes

* `getAverageTraits(orderPhylo = FALSE)` returns one row per species instead of repeating species, usually once per grid cell (#200).
* `plotSpatialPhylo()`, `plotTraitDistribution(type = "all")` and `localPlots(plot = TRUE)` draw the maps with the x-axis of the grid horizontally and in the aspect ratio of the grid; before, x and y were swapped and the maps were stretched to a square or to the plot panel (#211).
* `plotSpatialPhylo(plotTraits = FALSE)` shows the panel titles, which were drawn above the top edge of the plot (#217).
* `rejectionSample()` and `rejectionSamplePhylosimAnalysis()` work with a single parameter or a single summary statistic instead of stopping with an error (#201).
* `runSimulation()` with `fission = 1` or `2`: without protracted speciation, each speciation event creates a new species. Before, two events on the same species in one generation created only one, and with `fission = 2`, events on species within one or two adjacent columns created none. `fission = 2` now splits a species by rank along the axis on which it extends farther, measured on the torus, instead of at its mean x coordinate. With `protracted > 0`, a new event no longer restarts the incipient stage of individuals marked by an earlier event. Results with fission change (#208).

# EcoPhyloSim 0.3.2

Results differ from 0.3.1 for the same parameters and seed, because the C++ model uses R's random number generator (#14); with `negativeDD > 0`, the distributions of results also change (#144, #165).

## New features

* `createCompletePar()`: new arguments `negativeDDKernel` and `negativeDDWidth` measure the similarity to the neighbors by the competition trait (absolute difference or Gaussian kernel) or by species identity (conspecific density dependence without burn-in) (@AndyEnvFr, #104, #133, #144).
* `createCompletePar()`: new argument `negativeDDResponse` (experimental) sets the term of negative density dependence as a function of the similarity to the neighbors (stored as its values) or as a vector of such values (#144).
* `createCompletePar()`: new arguments `positiveDD`, `positiveDDCut`, `positiveDDKernel` and `positiveDDWidth` (experimental) for positive density dependence with its own strength, radius and similarity measure; `NeutralMod()` gets `positiveDD` and `positiveDDCut` (@AndyEnvFr, #104, #108).
* `createCompletePar()`: `fitnessActsOn` accepts `"recruitment"` (fitness of the offspring in the free cell) and combinations such as `c("mortality", "recruitment")` (#165).
* `runSimulationBatch()`: new argument `onError = "NA"` returns NA for failed simulations instead of stopping (#158).
* `highlightLocalPhylo()`: new argument `title` replaces the default title of the plots (e.g. `title = ""` for none) (#180).

## Changes

* `density` and `densityCut` in `createCompletePar()` and `NeutralMod()` were renamed to `negativeDD` and `negativeDDCut`; the old names still work with a deprecation warning, and parameter lists with the old names are converted (#144). Results store the new names in `$Model`, so code that reads `$Model$density` or extracts it with `extractParams()` must use `negativeDD`; `runPhylosimAnalysis()` renames the old names with a deprecation warning.
* Environment and density dependence now combine multiplicatively (their death rates for mortality, their fitness terms for reproduction and recruitment) instead of adding their fitness terms; this changes all results with `environment > 0` and `negativeDD > 0`, which were neutral with `fitnessActsOn = "mortality"` if `environment + negativeDD <= 1` (#144). Models without `negativeDD` are not affected by this change.
* The default `negativeDDKernel = "absolute"` averages the trait differences over the `k` neighbors instead of dividing their sum by `k + 1`, which had reduced their effect by `k / (k + 1)`; this changes results with `negativeDD > 0` (#144).
* `NeutralMod()` with `negativeDD > 0` draws parents until a recruit is accepted, as in `runSimulation()` with `fitnessActsOn = "recruitment"`: with global dispersal, after `min(1000, x * y - 1)` rejections the parent is drawn directly with fitness-weighted probabilities; with a dispersal function, the 1000th recruit is accepted. Before, the 10th recruit was accepted without a test. This changes results with `negativeDD > 0` (#165).
* The C++ simulation uses R's random number generator, which `runSimulation()` sets from `seed` and restores afterwards, so that results for a given seed no longer depend on the compiler's random number library (but differ from earlier versions); this also seeds the summary statistics and replaces their separate seed handling from 0.3.1 (#135). `callModel()` has no `seed` argument any more, and `createCompletePar()` draws a missing `seed` from a larger range (#14).
* `callModel()`, the interface to the C++ simulation that `runSimulation()` calls, is no longer exported; use `runSimulation()` instead (#74).
* `citation("EcoPhyloSim")` is generated from DESCRIPTION instead of the removed `inst/CITATION` file, so that it gives the CRAN URL and DOI after a CRAN release (#30).

## Bug fixes

* `plotSpatialPhylo()` and `highlightLocalPhylo()` scale the tip labels down for phylogenies with many species, where they overlapped (#180).
* Functions that need the phylogeny (e.g. `plot()`, `summary()`, `collessImbalance()`, `nullModel()`, `calculatePhylogeneticDispersion()`) stop with a clear message for results with a single species, whose phylogeny is 0, instead of the error "$ operator is invalid for atomic vectors" (#180).
* `nullModel()` with `fun = "mpd"` and `calculatePhylogeneticDispersion()` with `type = "PhylMeta"` or `"PhylSample"` give p = NA instead of p = 0 for plots with a single species (#186), and instead of an error in `ecdf()` if all null communities have a single species, e.g. with a plot size of 1 (#196).
* `plotPhylogeneticDispersion()` fills the square of a scenario also if a plot size has only NA p-values (grey if all are NA) or the mean p-value is 1, where the square was left empty, and ignores NA values in the means over simulations (`means > 1`) (#195).
* `nullModel()` and `calculatePhylogeneticDispersion()` stop with an error, and `calculateSummaryStatistics()` gives `dispersion = NA` with a message, for trees with more than 23170 extant species (in a binary tree) on Windows, where `cophenetic()` crashed R (#168).
* `nullModel(fun = "pd")` works for phylogenies with a polytomy at the root, which failed because `picante::pd()` treated them as unrooted (#176).
* `extractParams()` returns one row per parameter set also if a parameter has several values (e.g. `runs` or `fitnessActsOn`); such values are combined into one string.
* `calculateSummaryStatistics()` and `runSimulation(calculateSummaries = TRUE)` no longer fail on grids with fewer than 15 cells in x or y; the statistics from plots are NA on such grids (#157).
* `runPhylosimAnalysis()`: a failed simulation gives NA at its position instead of replacing the results of its whole chunk by a single NA, so that `chunkwiseCleanup()` removes only that simulation (#158).
* `runSimulationBatch()`, `runPhylosimAnalysis()`, `chunkwiseCalculateSummaries()`, `chunkwiseCleanup()`, `rejectionSample()` and `rejectionSamplePhylosimAnalysis()`: the parallel workers load EcoPhyloSim from the library paths of the session; before, they could load another installed version from their default library (#190).
* `chunkwiseCleanup()` gives a message with the number of removed simulations; before, it removed them silently, e.g. all simulations on grids with fewer than 15 cells in x or y (#191).

## Documentation

* The vignette is split into "EcoPhyloSim user guide" (`vignette("PhyloSimManual")`: running the model, plots and summary statistics, fitting it to data) and "EcoPhyloSim model" (`vignette("ModelDescription")`: model overview, processes and community patterns). The user guide has new chapters "Quick start" and "Fitting the model to data" (approximate Bayesian computation and neural simulation-based inference) and explanations of the plots and summary statistics; the model vignette has three chapters: an overview, the processes with their equations, parameters and specific patterns, and the patterns of the community that arise from the processes together (#180).
* New vignette "Model history" and data set `modelHistory`: which changes of the model code alter simulation results, and how the patterns of the vignette "EcoPhyloSim model" compare between releases; the scripts that produce the data are in the folder `PhyloSim/data-raw/model-behavior` of the repository (#180, #181).
* The data sets `packageBehavior` and `pPD.pvalues` were simulated again with 0.3.2; `packageBehavior` has new indicators, the runtime and all censuses, and no longer `pilot` and `cnddSpecies` (#180).
* The help of `createCompletePar()` describes `fitnessBaseMortalityRatio` as the ratio of death rates per mechanism (equivalent to the previous definition for a single mechanism) and the neighborhoods of density dependence, and says that with `protracted > 0` the offspring of an incipient individual do not inherit its incipient state; the help of `getMortalityData()` explains why effects estimated from census data are smaller than in the model (#92, #143, #144).
* New help page `?PhylosimList` for the class returned by `runSimulationBatch()` and `createPhylosimList()`; the help of `rac()` and `sac()` says that `which.result` accepts several results (#92).

# EcoPhyloSim 0.3.1

## New features

* New pipeline for large simulation experiments: `runPhylosimAnalysis()` runs many simulations with fixed and prior-drawn parameters and saves them in chunks; `chunkwiseCalculateSummaries()`, `chunkwiseCleanup()`, `combineChunks()`, `prepareForAnalysis()` and `extractParams()` process the results (#85).
* New functions `rejectionSample()` and `rejectionSamplePhylosimAnalysis()` for ABC rejection sampling on simulation outputs (#85).
* New functions `getSpecies()`, `getSpeciesTrait()` and `getSpeciesTraits()` extract the extant species and their traits from a simulation (#85).
* New function `getMortalityData()` returns, for each individual and census interval, its species, number of conspecific neighbors and whether it died before the next recorded generation, for analyses of conspecific negative density dependence (@AndyEnvFr, #100).
* `calculateSummaryStatistics()` returns additional statistics (Clark-Evans index, 8-neighbor index, environmental sorting of species) (#85, #126).
* `createCompletePar()`: new argument `nicheWidth` sets the width (standard deviation) of the environmental niche, which was fixed before (@AndyEnvFr, #103, #138).
* The package has a CITATION file; `citation("EcoPhyloSim")` now also gives the URL of the GitHub repository, and `DESCRIPTION` lists the repository and its issue tracker (#30).
* New data set `packageBehavior` with summaries of simulation experiments that show how the outputs of the model respond to its parameters (#136).

## Changes

* The package was renamed from PhyloSim to EcoPhyloSim because of a name clash with the archived CRAN package phylosim; class and function names are unchanged (#64).
* `calculateSummaryStatistics()`: new argument `strict` (default `FALSE`); with `strict = TRUE`, the result is `NA` if any statistic cannot be calculated (#85, #126).
* `createCompletePar()` and `runSimulation()` stop with an error if `fitnessActsOn` is `"reproduction"` or `"both"` with global dispersal and `environment`, `density`, `redQueen` or `redQueenStrength` > 0, because this is not implemented (#98, #101).
* `createCompletePar()` and `runSimulation()` check their parameters with checkmate, also in parameter lists changed by hand, as described in the help of `createCompletePar()`; newly rejected are, among others, `density` > 1 (#96), decreasing `runs`, negative `redQueen` and `redQueenStrength`, and non-integer values of integer parameters (#85, #118).
* `createCompletePar()` and `runSimulation()` reject `specRate` > `x * y` (values above about 2.1e9 made the simulation hang, #150), and dispersal or `densityCut` kernels wider than the grid, which counted cells more than once: `2 * floor(2 * dispersal) + 1` and, if `density` > 0, `2 * densityCut + 1` must not exceed `min(x, y)` (#152).
* The `Model` element of the output of `runSimulation()` keeps the numeric values of `density` and `environment` instead of `TRUE`/`FALSE` (#91, #118).
* `runSimulationBatch()`: errors in single simulations are no longer silently turned into `NULL` during parallel runs, and `backup = TRUE` saves each result as `<scenario>.rds` (#106).
* Default plot sizes in `sac()` are log-spaced, as documented (#115).
* `localPlots()` stops with an error if `n` is not a single positive whole number (checkmate was added to Imports for such checks), and the species columns of its community table are sorted by species number and contain integer counts (#118, #119).
* The archived packages kulife, apTreeshape and PhyloMeasures are no longer needed: the Colless index is computed internally, and the phylogenetic dispersion uses `picante::ses.mpd()` with 999 randomizations, so the `dispersion` statistic of `calculateSummaryStatistics()` is stochastic and slower than before. The `useApTreeshape` argument of `collessImbalance()` is deprecated; `useApTreeshape = FALSE` no longer returns the normalized index (@melina-leite, #89).
* spatstat.geom, spatstat.explore and vegan were added to Imports; rmarkdown, testthat, abctools, ranger, data.table and gtools to Suggests (the last four are only needed for `runPhylosimAnalysis()`, its helper functions and the rejection sampling functions) (#113).
* License changed from GPL (>= 2) to GPL (>= 3).
* EcoPhyloSim now requires R (>= 4.1.0), as required by its dependency vegan (#113).
* adephylo, fields, MESS and shape were removed from Imports: adephylo was not used, and the few functions used from the others were replaced by base R code; `calculateSummaryStatistics()` gives the same values, and `plotPhylogeneticDispersion()` looks the same (#167).
* The unused C++ code for reading parameters from XML files and the bundled pugixml library were removed, and with them the pugixml authors in `Authors@R` (#170).
* The `plot()`, `print()` and `summary()` methods for `PhyloSim` and `PhylosimList` objects use the argument names of the generics (`x`, `object`), so calls such as `plot(simu = s)` no longer work (#113).
* `collessImbalance()` and `summary()` use only extant species by default (`dropFossils = TRUE`), as `calculateSummaryStatistics()` already did (#117).
* The non-functional Leipzig model (`type = "Leipzig"`) and its C++ code were removed, and `runSimulation()` no longer accepts `type = "Rneutral"`. Unknown values of `type` now give an error instead of returning `NULL` (#112).
* `NeutralMod()` was rewritten: argument names as in `createCompletePar()` (`xdim`/`ydim` renamed to `x`/`y`), new defaults (`x = 50`, `y = 50`, `runs = 100`, `specRate = 1`), local dispersal, density dependence, speciation at single births, and a `PhyloSim` object as output; the dead cell is never its own parent, as in `runSimulation()` (#112, #81).
* `plotSpatialPhylo()` and `plotTraitDistribution()` color the species with the colorblind-safe viridis palette, ordered along the phylogeny, so that a species has the same color in both plots (#156).
* `localPlots()` stops with an error if `size` is not a whole number between 1 and the shorter side of the grid, and the default plot sizes of `sac()` go up to the shorter side of the grid instead of `sqrt(x * y)` (#151).
* `runSimulation()`: the warning for results with only one species says that the phylogeny is set to 0 and how to avoid this, instead of "Cannot build Phylogeny" (#153).
* `runSimulation()`: with `prunePhylogeny = TRUE` (the default), the phylogeny contains only the extant species, so that its size no longer grows with the run length; extinct ancestors become internal nodes. Each census is placed half a generation after the events of its generation, and the tree has no root edge (see `?runSimulation`), so tree statistics change slightly. Writing large phylogenies is faster (#147).
* Bingxue Wu was added as author and further people as contributors in `DESCRIPTION`; the file `CONTRIBUTORS.md` in the repository lists all contributions (#136).

## Bug fixes

* Speciation by fission (`fission = 1` or `2`) assigned the wrong individuals to the new species (#85).
* `runSimulationBatch()`: `strip = "summaries"` was ignored without parallelization, and the results were not named after their scenarios (#106).
* Installation from GitHub failed because `src/RcppExports.cpp` was missing (#102).
* `createCompletePar()` and `runSimulation()`: a matrix passed as `airmat` was assigned to the wrong cells (and memory outside of it was read) if the grid is not square (`x != y`) (#122).
* `calculateSummaryStatistics()`: `varPart_2` returned the same value as `varPart_1` instead of the fraction explained by competition (#113).
* `rac()`: a vector in `which.result` gave an error, and the bar plot had no rank labels (#113).
* `sac()` with a single result (including the default, the last result) computed the curve from the first result; this also affected the `sacAuc` statistic of `calculateSummaryStatistics()` (#113).
* `summary()` of a `PhylosimList` printed the column headers for the number of species and the number of runs in the wrong order (#113).
* `getAverageTraits()` did not accept the output of `runSimulation()` (#113).
* `localPlots()` with `n = 1` returned two subplots (#116).
* `runSimulation()` read memory outside of `airmat`, or crashed R, if its length did not match the grid; `callModel()` stops with an error instead of crashing R for invalid input such as `mortStrength < 1` (#13, #94, #118).
* `runSimulation()`: the phylogeny was cut off, so that reading it failed, with `prunePhylogeny = FALSE`, for early values in `runs`, or with `fission`; with `prunePhylogeny = FALSE` and several values in `runs`, R crashed (#82, #97, #105).
* `runSimulation()`: with `prunePhylogeny = TRUE`, the returned phylogeny contained extinct lineages without living descendants, and with either setting, each value in `runs` leaked a copy of all species ever created, so that memory use grew with the run length and the number of values in `runs`. Such lineages are now removed during the run; the simulated communities are unchanged. With `prunePhylogeny = FALSE`, the branches of extant species in intermediate results no longer extend to the last value in `runs` (#140).
* `runSimulation()` with global dispersal: the dying individual could be its own parent, so that living species could get an extinction date (removed by `ape::drop.fossil()`, negative branch lengths, wrong species trait means); results with global dispersal change for a given seed. Grids must have at least 2 cells (#81).
* `plotTraitDistribution()` changed the random number stream of the user, and its histograms were cut off above log(count) = 7 (#156).
* The bundled data sets `simu.envcom.local`, `simu.neutral.local` and `simu.neutral.global` had the class `c("list", "Phylosim")` instead of `PhyloSim`, so that `print()`, `summary()`, `plot()` and `calculatePhylogeneticDispersion()` did not recognize them (#128).
* `calculatePhylogeneticDispersion()`: a vector in `which.result` gave an error (#131).
* `nullModel()`, and `calculatePhylogeneticDispersion()` with `type = "PhylMeta"` or `"PhylSample"`, drew null communities with `localPlotSize` instead of `localPlotSize^2` individuals, so that the p-values compared communities of different size; `localPlotSize` must now be a positive whole number. The example data set `pPD.pvalues` was recreated with the corrected null model (#132, #136).
* `calculateSummaryStatistics()` mixed the first and last result and could fail if `runs` had several values, also via `runSimulation(calculateSummaries = TRUE)`; all statistics now use the last result, or the one given by the new argument `which.result` (#134).
* `summary()` of a `PhyloSim` object printed the class name as `Phylosim` (#155).
* `runSimulation()` crashed R with a stack overflow on large grids with global dispersal and `environment`, `density`, `redQueen` or `redQueenStrength` > 0, or with large local dispersal kernels (#149).
* `localPlots()` and `sac()` failed with "subscript out of bounds" if the grid is not square (`x != y`); on square grids, the `positions` of the subplots had x and y swapped. The fix changes `sac()` results for a fixed seed (#151).
* `sac()`: a vector in `which.result` gave an error (#151).
* `runSimulation()`: species that went extinct in the last generation before a census were kept by `ape::drop.fossil()`; with `prunePhylogeny = TRUE`, the species pool of the `dispersion` statistic of `calculateSummaryStatistics()` contained extinct ancestors (#147).
* `runSimulation()`: without `airmat`, the default environmental gradient exceeded 1 for an odd number of columns `x` (e.g., 1.2 for `x = 5`); results for odd `x` change (#142).
* `calculateSummaryStatistics()`: the Clark-Evans index `ce` used a wrong row for each cell, so that cells were moved to other rows or, on grids with more columns than rows, dropped; `ce` changes on all grids (#163).
* `calculateSummaryStatistics()`: `envIntercept`, `env` and `envSqr` (binomial GLM coefficients) were `NaN` under strong environmental filtering because the fits failed, and used only the first row of the environment; they are replaced by `envSorting`; `ce` is `NA` instead of `NaN` if no species has a finite Clark-Evans index (#126).
* R crashed on Windows in `calculateSummaryStatistics()`, `summary()`, `collessImbalance()`, `getAverageTraits()`, `nullModel()`, `calculatePhylogeneticDispersion()`, `plotSpatialPhylo()`, `plotHighlightLocalPhylo()` and `plotTraitDistribution()` if the phylogeny had more than 23170 tips (possible with `prunePhylogeny = FALSE`), because `ape::drop.fossil()` and `cophenetic()` build a matrix of all node distances. Extinct species are now removed in linear time, and the `dispersion` statistic of `calculateSummaryStatistics()` uses only the extant species as species pool also with `prunePhylogeny = FALSE`, which changes its results for this setting. More than 23170 extant species still crash `cophenetic()` (#168).
* `calculateSummaryStatistics()` and `summary()`: the gamma statistic was wrong for phylogenies with polytomies (with `convertToBinaryTree = FALSE`), and `summary()` suppressed the warning; polytomies are now resolved with `ape::multi2di(random = FALSE)` before the gamma statistic is computed (#169).
* `runSimulation()`: the same `seed` gave different phylogenies and summary statistics, because polytomies were resolved at random and the summary statistics use R's random number generator. Polytomies are now resolved deterministically, and with `calculateSummaries = TRUE`, R's seed is set from `seed` and restored afterwards. For a given seed, the phylogeny (with `convertToBinaryTree = TRUE`), `collessImbalance()` and the summary statistics change; with `convertToBinaryTree = TRUE`, the gamma statistic does not (#135).

## Documentation

* The vignette was converted from Sweave (PDF) to R Markdown (HTML).
* Package help converted to roxygen; examples moved to `inst/examples/` (#87).
* The help of `createCompletePar()` explains vectors in `runs`, the valid ranges of `environment`, `airmat` and `fitnessBaseMortalityRatio` (#85), `density` (#96), `dispersal`, `densityCut`, `redQueen` and `redQueenStrength` (#118), and the dimensions of `airmat` (#122); it also defines `fitnessBaseMortalityRatio` and explains that only multiples of 0.5 of `dispersal` have an effect (#92, #141).
* Help pages and vignette revised: corrected trait-evolution formula in the vignette and return value of `sac()`, documented parameters of `callModel()`, faster examples, corrected spelling and several incorrect argument and return value descriptions (#113).
* The vignette was restructured: an abstract and a quick start; a model description with a schematic of a generation, parameter tables with defaults and ranges, and known limitations; a description of the output; and a new chapter "Emerging patterns" that compares how the outputs respond to the parameters with expectations from theory and the literature, using the data set `packageBehavior` (#77, #99, #136).
