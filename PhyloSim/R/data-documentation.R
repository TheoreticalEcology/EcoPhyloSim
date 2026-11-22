#' Simulation output of community assembly with local dispersal under the influence of environment and competition
#'
#' @name simu.envcom.local
#' @format A list containing two lists: The model output and the model input. The model output is described here
#' \describe{
#'   \item{specMat}{species matrix, a matrix containing the species ID of the individual in each cell}
#'   \item{traitMat}{trait matrix, a matrix containing the environmental trait of the individual in each cell (corresponding to the species matrix)}
#'   \item{envMat}{environment matrix, contains the environmental conditions for each cell}
#'   \item{compMat}{competition matrix, contains the competition attribute for each cell}
#'   \item{neutMat}{neutral matrix, contains the neutral traits}
#'   \item{phylogeny}{phylogeny of the community, an object of class 'phylo'}
#'   \item{phyloTXT}{the phylogeny in Newick format}
#'    }
#' @source created with \code{\link{runSimulation}}
NULL





#' Simulation output of community assembly with local dispersal under neutral conditions
#'
#' @name simu.neutral.local
#' @format A list containing two lists: The model output and the model input. The model output is described here
#' \describe{
#'   \item{specMat}{species matrix, a matrix containing the species ID of the individual in each cell}
#'   \item{traitMat}{trait matrix, a matrix containing the environmental trait of the individual in each cell (corresponding to the species matrix)}
#'   \item{envMat}{environment matrix, contains the environmental conditions for each cell}
#'   \item{compMat}{competition matrix, contains the competition attribute for each cell}
#'   \item{neutMat}{neutral matrix, contains the neutral traits}
#'   \item{phylogeny}{phylogeny of the community, an object of class 'phylo'}
#'   \item{phyloTXT}{the phylogeny in Newick format}
#'    }
#' @source created with \code{\link{runSimulation}}
NULL







#' Simulation output of community assembly with global dispersal under neutral conditions
#'
#' @name simu.neutral.global
#' @format A list containing two lists: The model output and the model input. The model output is described here
#' \describe{
#'   \item{specMat}{species matrix, a matrix containing the species ID of the individual in each cell}
#'   \item{traitMat}{trait matrix, a matrix containing the environmental trait of the individual in each cell (corresponding to the species matrix)}
#'   \item{envMat}{environment matrix, contains the environmental conditions for each cell}
#'   \item{compMat}{competition matrix, contains the competition attribute for each cell}
#'   \item{neutMat}{neutral matrix, contains the neutral traits}
#'   \item{phylogeny}{phylogeny of the community, an object of class 'phylo'}
#'   \item{phyloTXT}{the phylogeny in Newick format}
#'    }
#' @source created with \code{\link{runSimulation}}
NULL

#' Example p-values for the phylogenetic dispersion null model, used to demonstrate \code{\link{plotPhylogeneticDispersion}}
#'
#' @name pPD.pvalues
#' @format A list of length 1, wrapping a list of 28 scenarios (one for each combination of the 7 x-axis and 4 y-axis positions in \code{\link{pPD.positions}}). Each scenario is itself a list of 9 elements, named "2" to "10", corresponding to the sample plot sizes (argument \code{plotlengths}) used when the data was generated with \code{\link{calculatePhylogeneticDispersion}}. Each of these holds a numeric vector of p-values, one for each of the 200 sample plots drawn for that scenario/plot-size combination.
#' @source created with the script \code{PhyloSim/data-raw/model-behavior/pPD.R} in the package
#' repository at each release: one run of 50,000 generations per scenario on a 50 x 50
#' grid, and the "PhylMeta" null model of \code{\link{calculatePhylogeneticDispersion}}
#' with 500 replicates; \code{attr(pPD.pvalues, "provenance")} records the run. The
#' data only illustrate the plot.
NULL


#' Example axis positions and labels for \code{\link{plotPhylogeneticDispersion}}, matching \code{\link{pPD.pvalues}}
#'
#' @name pPD.positions
#' @format A list with 5 entries, determining where the scenarios in \code{\link{pPD.pvalues}} should be plotted and how the axes should be labeled
#' \describe{
#'   \item{x}{integer vector of length 7, the positions along the x-axis}
#'   \item{y}{numeric vector of length 4, the positions along the y-axis}
#'   \item{xname}{character vector of length 7, the labels for the lower x-axis (density dependence strength)}
#'   \item{zname}{character vector of length 7, the labels for the upper x-axis (environmental trait strength)}
#'   \item{yname}{character vector of length 4, the labels for the y-axis (dispersal mode)}
#'    }
#' @source created together with \code{\link{pPD.pvalues}}
NULL


#' Summaries of simulations on how the model outputs depend on the parameters
#'
#' Used in the vignette "EcoPhyloSim model" (chapters "Overview", "Processes" and "Community patterns") to show
#' how the outputs of the model respond to its parameters, how fast they reach equilibrium and how long the simulations take.
#'
#' @name packageBehavior
#' @format A list with two entries, \code{main} and \code{info}. \code{main}
#' contains:
#' \describe{
#'   \item{settings}{data frame with one row per parameter setting: \code{block}
#'   (group of settings in which one or two parameters are varied), \code{setting}
#'   (number within the block), the parameters (names of release 0.3.1, e.g.
#'   \code{density} and \code{densityCut} for \code{negativeDD} and
#'   \code{negativeDDCut}, and the current names for parameters added later;
#'   \code{dispersal} as character) and the number of \code{generations}}
#'   \item{census}{data frame with the replicate \code{mean}, \code{sd} and number of
#'   replicates \code{n} (without \code{NA} values, e.g. tree statistics of runs
#'   with fewer than 3 species) of each \code{indicator} per setting and census
#'   (\code{generation}): species richness, Shannon diversity, evenness (Shannon
#'   diversity / log(richness)), dominance (share of the most abundant species),
#'   share of the 4 neighbors that belong to the same species, correlation between
#'   environmental trait and environment, sd of the competition and of the neutral
#'   trait, Colless index (Yule-normalized, \code{\link{collessImbalance}}), gamma
#'   statistic (\code{\link[ape]{gammaStat}}) and height of the tree of the extant
#'   species, mean species age and rank correlation between species age and
#'   abundance (\code{meanSpeciesAge}, \code{corAgeAbundance}; the age is estimated
#'   from the species ID, which is counted up at each speciation), and, at the last
#'   census only, the mean p-value of the phylogenetic dispersion of 20 local plots
#'   of 10 x 10 cells (\code{mpdLocal}, \code{\link{nullModel}})}
#'   \item{sac}{data frame with the replicate mean, sd and n of the species richness
#'   per plot size (\code{size}, in cells) at the last census, from \code{\link{sac}}}
#'   \item{rac}{data frame with the replicate mean, sd and n of the abundance per rank
#'   at the last census}
#'   \item{decay}{data frame with the replicate mean, sd and n of the share of cells
#'   at a given \code{distance} (along rows and columns) that contain the same
#'   species, at the last census}
#'   \item{cndd}{data frame with the replicate mean, sd and n of the coefficient of
#'   the number of conspecifics within radius 2 in a logistic regression of death
#'   (with one intercept per species) over the last five generations, from
#'   \code{\link{getMortalityData}}}
#'   \item{runtime}{data frame with the replicate mean, sd and n of the runtime of
#'   the core simulation in seconds}
#'    }
#' \code{info} contains metadata for reproducibility: the name of the archive of
#' the run and the simulated commit (\code{archive}, \code{commit}), the replicates
#' and generations per block (\code{modes}), the version of the scripts
#' (\code{scripts}), the provenance of the run (\code{run}, from
#' \code{RUNINFO.txt} of the archive), and the number of simulations, failed
#' simulations and hours of simulation time (\code{runs}, \code{errors},
#' \code{simulationHours}; sum of the runtimes of the core simulations).
#' @source simulated and built with the scripts in the folder
#' \code{PhyloSim/data-raw/model-behavior} of the package repository
#' (\url{https://github.com/TheoreticalEcology/EcoPhyloSim}), from the run of the
#' current release (archive named in \code{info$archive}); its \code{README.md}
#' describes how to regenerate the data.
NULL


#' Model history: simulation results across releases
#'
#' Summaries of the same simulation design run with each release of the package,
#' and the changes of the model code that alter simulation results. Used in the
#' vignette "Model history".
#'
#' @name modelHistory
#' @format A list with three entries:
#' \describe{
#'   \item{summary}{data frame with one row per run of a release (\code{version},
#'   \code{commit}), setting and indicator: \code{block} (group of settings in which
#'   one or two parameters are varied), \code{setting} (number within the block),
#'   the parameters of the setting (names of release 0.3.1, e.g. \code{density} for
#'   \code{negativeDD}, and the current names for parameters added later;
#'   \code{dispersal} as character), the number of
#'   \code{generations}, the \code{indicator} and the replicate \code{mean},
#'   \code{sd} and number of replicates \code{n} (without \code{NA} values) at the
#'   last census. Indicators as in \code{\link{packageBehavior}} (richness,
#'   shannon, conspecificNeighbors, corTraitEnv, sdCompetitionTrait, colless,
#'   gamma, meanSpeciesAge, corAgeAbundance, cndd) and the
#'   \code{runtime} of the core simulation in seconds, the share of failed
#'   simulations (\code{failed}), treeHeight, evenness, dominance,
#'   sdNeutralTrait and mpdLocal (see the vignette "Model history").}
#'   \item{changes}{data frame of the releases and the changes of the model code
#'   that alter results: \code{date}, \code{kind} (release or change), release
#'   \code{version}, \code{commit}, issue \code{reference}, \code{type} (N: results
#'   per seed change, but no change of their distributions was detected; B: patterns change; S: parameters renamed or redefined;
#'   S+B: both), the \code{change}, its effect on the \code{outputs} and the
#'   \code{evidence}.}
#'   \item{info}{the runs used, the commit of the EcoPhyloSim that computed the
#'   indicators of all runs (\code{indicators}), the version of the scripts and of R.}
#' }
#' @source Simulated and built with the scripts in the folder
#' \code{PhyloSim/data-raw/model-behavior} of the package repository
#' (\url{https://github.com/TheoreticalEcology/EcoPhyloSim}); its
#' \code{README.md} describes how to regenerate the data.
NULL
