#' Simulation output of community assembly with local dispersion under the influence of environment and competition
#'
#' @name simu.envcom.local
#' @format A list containing two lists: The model output and the model input. The model output is described here
#' \describe{
#'   \item{specMat}{species matrix, a matrix representing the species abundance}
#'   \item{traitMat}{trait matrix, a matrix containing the environmental trait for each species (corresponding to the species matrix)}
#'   \item{envMat}{environment matrix, contains the environmental conditions for each cell}
#'   \item{compMat}{competition matrix, contains the competition attribute for each cell}
#'   \item{neutMat}{neutral matrix, contains the neutral traits}
#'   \item{phylogeny}{phylogeny of the community stored in a list}
#'    }
#' @source created with \code{\link{runSimulation}}
NULL





#' Simulation output of community assembly with local dispersion under neutral conditions
#'
#' @name simu.neutral.local
#' @format A list containing two lists: The model output and the model input. The model output is described here
#' \describe{
#'   \item{specMat}{species matrix, a matrix representing the species abundance}
#'   \item{traitMat}{trait matrix, a matrix containing the environmental trait for each species (corresponding to the species matrix)}
#'   \item{envMat}{environment matrix, contains the environmental conditions for each cell}
#'   \item{compMat}{competition matrix, contains the competition attribute for each cell}
#'   \item{neutMat}{neutral matrix, contains the neutral traits}
#'   \item{phylogeny}{phylogeny of the community stored in a list}
#'    }
#' @source created with \code{\link{runSimulation}}
NULL







#' Simulation output of community assembly with global dispersion under neutral conditions
#'
#' @name simu.neutral.global
#' @format A list containing two lists: The model output and the model input. The model output is described here
#' \describe{
#'   \item{specMat}{species matrix, a matrix representing the species abundance}
#'   \item{traitMat}{trait matrix, a matrix containing the environmental trait for each species (corresponding to the species matrix)}
#'   \item{envMat}{environment matrix, contains the environmental conditions for each cell}
#'   \item{compMat}{competition matrix, contains the competition attribute for each cell}
#'   \item{neutMat}{neutral matrix, contains the neutral traits}
#'   \item{phylogeny}{phylogeny of the community stored in a list}
#'    }
#' @source created with \code{\link{runSimulation}}
NULL

#' Example p-values for the phylogenetic dispersion null model, used to demonstrate \code{\link{plotPhylogeneticDispersion}}
#'
#' @name pPD.pvalues
#' @format A list of length 1, wrapping a list of 28 scenarios (one for each combination of the 7 x-axis and 4 y-axis positions in \code{\link{pPD.positions}}). Each scenario is itself a list of 9 elements, named "2" to "10", corresponding to the sample plot sizes (plotlengths) used when the data was generated with \code{\link{calculatePhylogeneticDispersion}}. Each of these holds a numeric vector of p-values, one for each of the 200 sample plots drawn for that scenario/plot-size combination.
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