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
#' @source created with \code{\link{calculatePhylogeneticDispersion}}
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
#' Used in the vignette (chapter "Emerging patterns") to show how the outputs of
#' the model respond to its parameters and how fast they reach equilibrium.
#'
#' @name packageBehavior
#' @format A list with three entries: \code{pilot} (few replicates, long runs, used
#' to check the time to equilibrium), \code{main} (more replicates, used for the
#' responses to the parameters) and \code{info}. \code{pilot} contains
#' \code{settings} and \code{census}, \code{main} all of the following:
#' \describe{
#'   \item{settings}{data frame with one row per parameter setting: \code{block}
#'   (group of settings in which one or two parameters are varied), \code{setting}
#'   (number within the block), the arguments passed to
#'   \code{\link{createCompletePar}} and the number of generations}
#'   \item{census}{data frame with the replicate \code{mean}, \code{sd} and number of
#'   replicates \code{n} of each indicator per setting and census (\code{generation};
#'   in \code{main} only the last census):
#'   species richness, Shannon diversity, share of the 4 neighbors that belong to the
#'   same species, correlation between environmental trait and environment, sd of the
#'   competition trait, Colless index (Yule-normalized, \code{\link{collessImbalance}}),
#'   gamma statistic (\code{\link[ape]{gammaStat}}), mean species age and rank
#'   correlation between species age and abundance (\code{meanSpeciesAge},
#'   \code{corAgeAbundance}; the age is estimated from the species ID, which is
#'   counted up at each speciation)}
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
#'   \item{cnddSpecies}{data frame with the \code{slope} of the same regression for
#'   single species (at least 50 observations) and their \code{abundance}, for each
#'   replicate (\code{rep})}
#'    }
#' \code{info} contains metadata for reproducibility: the package commit and a
#' checksum of the script used (\code{run}), the date, package and R versions,
#' platform, operating system, processor, versions of the loaded packages, the
#' number of replicates and generations (\code{modes}), and the number of runs,
#' failed runs and the hours of simulation time (sum of the runtimes of the core
#' simulations) per mode.
#' @source created with the script \code{data-raw/packageBehavior.R} in the
#' package repository (\url{https://github.com/TheoreticalEcology/EcoPhyloSim})
NULL
