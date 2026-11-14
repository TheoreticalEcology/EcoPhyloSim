#' @title  Parameter Generator 
#' @description Function to create a list of parameters for biogeographical simulations with \code{\link{runSimulation}} or \code{\link{runSimulationBatch}}
#' @param x  Integer, Dimension of the model landscape in x-direction; the grid must have at least 2 cells (x * y >= 2)
#' @param y  Integer, Dimension of the model landscape in y-direction
#' @param dispersal Numeric or character. Type 0 or "global" for global dispersal. For local dispersal, a number >= 0.5 sets the dispersal distance; 2 * floor(2 * dispersal) + 1 must not exceed min(x, y).
#' @param runs  Integer or vector of Integers, Number of generations or sequence of generations the model runs over (see Details); values must not decrease.
#' @param specRate Numeric, speciation rate: mean number of speciation events per generation (the realized number is drawn from a Poisson distribution); must not exceed x * y
#' @param negativeDD Numeric between 0 and 1, strength of negative density dependence (formerly \code{density}). By default (0) there is none. Its fitness term is \code{1 - negativeDD * S}, where \code{S} is the similarity to the neighbors within \code{negativeDDCut} (see \code{negativeDDKernel}); with a function as \code{negativeDDResponse}, \code{negativeDD > 0} only switches it on.
#' @param environment Numeric, determining whether or how strong the environment influences the model.  By default (environment=0) there is no influence of the environment. The higher the value of the parameter, the stronger is the influence of the environment. Environment value must be between 0 and 1.
#' @param nicheWidth Numeric > 0, width (standard deviation) of the Gaussian environmental niche (no effect on the simulation if environment = 0). The default 0.03659906 is the value that was fixed in earlier versions. Values much smaller or larger than the spread of the environment (standard deviation about 0.29 for the default gradient) weaken environmental filtering; values below 0.003 or above 3 give a message.
#' @param fitnessActsOn Character vector, the life stages on which fitness acts: "mortality" (default; individuals with low fitness in their cell die more often), "reproduction" (parents with high fitness in their own cell produce more offspring), "recruitment" (the offspring of a parent is established in a free cell with a probability given by its fitness in that cell, evaluated with the traits of the parent), or a combination, e.g. \code{c("mortality", "recruitment")}; "both" is short for \code{c("mortality", "reproduction")}. With several stages, environment and density dependence act at each of them. With global dispersal, "reproduction" is only implemented if environment, negativeDD, positiveDD, redQueen and redQueenStrength are 0. redQueen and redQueenStrength do not act on recruitment. With "recruitment", parents are drawn until a recruit is accepted, which can be slow if most recruits have low fitness, e.g. with global dispersal at the start of a run, when all individuals belong to one species.
#' @param fitnessBaseMortalityRatio Integer >= 1, the largest ratio of death rates between individuals for each mechanism (environment, negative and positive density dependence) acting on mortality: an individual with fitness term \code{w} between 0 and 1 dies at a rate proportional to \code{1 + (fitnessBaseMortalityRatio - 1) * (1 - w)}. 1 gives no fitness differences in mortality. No effect if fitness does not act on mortality. See Details.
#' @param negativeDDCut Integer >= 1, radius of the neighborhood for negative density dependence (formerly \code{densityCut}; no effect if negativeDD = 0); must not exceed min(x, y), and 2 * negativeDDCut + 1 must not exceed min(x, y) if negativeDD > 0
#' @param fission Integer, determining which fission type should be used. Options are 0 (none = default), 1 (every second individual becomes part of new species) and 2 (population is geographically split in two parts).
#' @param protracted Integer, determining the time span in generations a new species stays 'incipient' before turning into a 'good' species. Default is 0.
#' @param redQueenStrength Numeric, determining the strength of the Red Queen effect. A value > 0 means a new species gets a fitness boost due to its novelty. Must not be negative.
#' @param redQueen Numeric, determining the strength of the fitness decline of an aging species. Must not be negative.
#' @param airmat matrix determining the environment of the simulation, with `y` rows and `x` columns (the grid dimensions); values must be between 0 and 1. A single value (the default) gives an environmental gradient along `x` that rises from 0 to 1 and falls back.
#' @param seed Integer, sets the random seed. If NULL (default), a random seed is drawn.
#' @param type Character, determining which model should be used. Currently only "base" (the default) is available. For a simple neutral model in R, see \code{\link{NeutralMod}}.
#' @param scenario String, further information you want to add to the parameter set in order to refer to a model run more conveniently. 
#' @param calculateSummaries Logical, determining whether summary statistics should be calculated
#' @param convertToBinaryTree Logical, determining if the phylogeny should be converted into a binary tree
#' @param negativeDDKernel Character, how the similarity \code{S} of an individual to its neighbors is measured for negative density dependence: "absolute" (default) uses the absolute difference in the competition trait, \code{S = 1 - D}, where \code{D} is the sum of the differences to the neighbors divided by the number of cells including the focal cell (so \code{S >= 1 / (k + 1)} for \code{k} neighbors); "gaussian" uses the mean of \code{exp(-0.5 * (difference / negativeDDWidth)^2)} over the neighbors; "species" uses the proportion of conspecifics among the neighbors, i.e., conspecific density dependence that acts from the first speciation on. The trait-based kernels act only once the competition traits have diverged, which takes thousands of generations, because all individuals start with the same traits.
#' @param negativeDDWidth Numeric > 0, width of the Gaussian kernel for \code{negativeDDKernel = "gaussian"}; values above 1 give a message, as the competition trait lies between 0 and 1.
#' @param negativeDDResponse "linear" (default) for the density term \code{1 - negativeDD * S}, or a function of the similarity \code{S} between 0 and 1 that returns the density term (between 0 and 1), e.g. \code{function(S) exp(-a * k * S)} for a survival (or establishment) factor \code{exp(-a)} per conspecific neighbor with \code{negativeDDKernel = "species"} and \code{k} neighbors. The function is evaluated before the run, exactly for each number of conspecifics with \code{negativeDDKernel = "species"}, otherwise at 1001 values with linear interpolation.
#' @param positiveDD Numeric between 0 and 1, strength of positive density dependence (e.g. facilitation by mutualists): individuals with dissimilar neighbors within \code{positiveDDCut} are disadvantaged, with the fitness term \code{1 - positiveDD * (1 - S)}. As the grid is always full, this is the same as an advantage of individuals with similar neighbors. By default (0) there is none.
#' @param positiveDDCut,positiveDDKernel,positiveDDWidth Radius, similarity measure (default "species") and Gaussian width for positive density dependence, as \code{negativeDDCut}, \code{negativeDDKernel} and \code{negativeDDWidth}.
#' @param density,densityCut Deprecated names of \code{negativeDD} and \code{negativeDDCut}.
#' @param prunePhylogeny Logical, whether extinct lineages without living descendants are removed from the phylogeny during the run, which keeps the tree and memory use small; the phylogeny then contains only the extant species (see \code{\link{runSimulation}}). If FALSE, all species are kept.
#' @details Fitness: environment and density dependence give fitness terms between 0 and 1, \code{w_env = environment * exp(-0.5 * ((E - t) / nicheWidth)^2) + 1 - environment} and the terms of negative and positive density dependence (see \code{negativeDD} and \code{positiveDD}). If several act on mortality, their relative death rates multiply (proportional hazards): an individual dies at a rate proportional to the product of \code{1 + (m - 1) * (1 - w)} over these terms, with \code{m = fitnessBaseMortalityRatio}. As the number of deaths per generation is fixed (\code{x * y}), only relative rates matter, and a term that is equal for all individuals has no effect. For reproduction and recruitment, the terms multiply directly.
#'
#' If runs is a sequence of generations, intermediate results are saved. E.g. when runs is c(500, 600, 700), the simulation runs 700 generations in total, and the intermediate results at generations 500 and 600 are saved additionally. The intermediate and end results are saved in the output of \code{runSimulation}.
#' @return A List with parameters
#' @example /inst/examples/parCreator-help.R
#' @export

createCompletePar <- function(x = 50, y = 50, dispersal = "global", runs = 100, specRate = 1.0, negativeDD = 0, environment = 0, nicheWidth = 0.03659906, fitnessActsOn = "mortality" , fitnessBaseMortalityRatio = 10, negativeDDCut = 1, seed=NULL,  type = "base", fission = 0, redQueen = 0, redQueenStrength = 0, protracted = 0, airmat = 1, scenario = NULL, calculateSummaries = FALSE, convertToBinaryTree = TRUE, prunePhylogeny = TRUE, negativeDDKernel = "absolute", negativeDDWidth = 0.1, negativeDDResponse = "linear", positiveDD = 0, positiveDDCut = 1, positiveDDKernel = "species", positiveDDWidth = 0.1, density, densityCut){

  # old names, before negative and positive density dependence were separated
  if (!missing(density)) {
    if (!missing(negativeDD)) stop("Give only 'negativeDD', not also the deprecated 'density'.", call. = FALSE)
    warning("Argument 'density' is deprecated; use 'negativeDD' instead.", call. = FALSE)
    negativeDD <- density
  }
  if (!missing(densityCut)) {
    if (!missing(negativeDDCut)) stop("Give only 'negativeDDCut', not also the deprecated 'densityCut'.", call. = FALSE)
    warning("Argument 'densityCut' is deprecated; use 'negativeDDCut' instead.", call. = FALSE)
    negativeDDCut <- densityCut
  }
    
  soilmat <- 1 # Needs to be defined here as long as it is not
               # implemented in the model
  if (length(airmat) != 1) checkmate::assert_matrix(airmat, nrows = y, ncols = x)

  if(length(soilmat) != 1){
    if((nrow(soilmat) != y) | (ncol(soilmat)!= x)) stop("Environment and matrix size do not match")
  }

  airmat <- as.numeric(airmat)
  
  # The soilmat is not implemented in the model so far.
  soilmat = 1
  soilmat <- as.numeric(soilmat)
  
  if (is.null(seed)) seed = sample(1:10000,1)
  
  par = list(x=x,y=y,dispersal = dispersal, runs = runs, specRate = specRate, 
             negativeDD = negativeDD, environment = environment, nicheWidth = nicheWidth, fitnessActsOn=fitnessActsOn,
             fitnessBaseMortalityRatio=fitnessBaseMortalityRatio, negativeDDCut = negativeDDCut, 
             seed = seed, type = type, scenario = scenario, fission = fission, 
             redQueen = redQueen, redQueenStrength = redQueenStrength, protracted = protracted,
             airmatR = airmat, soilmatR=soilmat, calculateSummaries=calculateSummaries,
             convertToBinaryTree = convertToBinaryTree, prunePhylogeny = prunePhylogeny,
             negativeDDKernel = negativeDDKernel, negativeDDWidth = negativeDDWidth,
             negativeDDResponse = negativeDDResponse, positiveDD = positiveDD, positiveDDCut = positiveDDCut,
             positiveDDKernel = positiveDDKernel, positiveDDWidth = positiveDDWidth)
  
  checkPar(par)

  return(par)
  
}

