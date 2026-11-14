#' @title  Parameter Generator 
#' @description Function to create a list of parameters for biogeographical simulations with \code{\link{runSimulation}} or \code{\link{runSimulationBatch}}
#' @param x  Integer, Dimension of the model landscape in x-direction; the grid must have at least 2 cells (x * y >= 2)
#' @param y  Integer, Dimension of the model landscape in y-direction
#' @param dispersal Numeric or character. Type 0 or "global" for global dispersal. For local dispersal, a number >= 0.5 sets the dispersal distance; 2 * floor(2 * dispersal) + 1 must not exceed min(x, y).
#' @param runs  Integer or vector of Integers, Number of generations or sequence of generations the model runs over (see Details); values must not decrease.
#' @param specRate Numeric, speciation rate: mean number of speciation events per generation (the realized number is drawn from a Poisson distribution); must not exceed x * y
#' @param density Numeric, determining whether or how strong the density dependence influences the model. By default (density=0) there is no density dependence. The higher the value of the parameter, the stronger is the density dependence. Density value must be between 0 and 1.
#' @param environment Numeric, determining whether or how strong the environment influences the model.  By default (environment=0) there is no influence of the environment. The higher the value of the parameter, the stronger is the influence of the environment. Environment value must be between 0 and 1.
#' @param nicheWidth Numeric > 0, width (standard deviation) of the Gaussian environmental niche (no effect on the simulation if environment = 0). The default 0.03659906 is the value that was fixed in earlier versions. Values much smaller or larger than the spread of the environment (standard deviation about 0.29 for the default gradient) weaken environmental filtering; values below 0.003 or above 3 give a message.
#' @param fitnessActsOn Character, determining how the fitness influences the individuals. Possible inputs are "mortality" (default), "reproduction" or "both". With global dispersal, "reproduction" and "both" are only implemented if environment, density, redQueen and redQueenStrength are 0.
#' @param fitnessBaseMortalityRatio Integer, determines the fitness based mortality ratio. Must be greater than or equal to 1.
#' @param densityCut Integer >= 1, defines the effective range of the competition (no effect on the simulation if density = 0); must not exceed min(x, y), and 2 * densityCut + 1 must not exceed min(x, y) if density > 0
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
#' @param prunePhylogeny Logical, whether extinct lineages without living descendants are removed from the phylogeny during the run, which keeps the tree and memory use small (see \code{\link{callModel}}). If FALSE, all species are kept.
#' @details If runs is a sequence of generations, intermediate results are saved. E.g. when runs is c(500, 600, 700), the simulation runs 700 generations in total, and the intermediate results at generations 500 and 600 are saved additionally. The intermediate and end results are saved in the output of \code{runSimulation}.
#' @return A List with parameters
#' @example /inst/examples/parCreator-help.R
#' @export

createCompletePar <- function(x = 50, y = 50, dispersal = "global", runs = 100, specRate = 1.0, density = 0, environment = 0, nicheWidth = 0.03659906, fitnessActsOn = "mortality" , fitnessBaseMortalityRatio = 10, densityCut = 1, seed=NULL,  type = "base", fission = 0, redQueen = 0, redQueenStrength = 0, protracted = 0, airmat = 1, scenario = NULL, calculateSummaries = FALSE, convertToBinaryTree = TRUE, prunePhylogeny = TRUE){
    
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
             density = density, environment = environment, nicheWidth = nicheWidth, fitnessActsOn=fitnessActsOn,
             fitnessBaseMortalityRatio=fitnessBaseMortalityRatio, densityCut = densityCut, 
             seed = seed, type = type, scenario = scenario, fission = fission, 
             redQueen = redQueen, redQueenStrength = redQueenStrength, protracted = protracted,
             airmatR = airmat, soilmatR=soilmat, calculateSummaries=calculateSummaries,
             convertToBinaryTree = convertToBinaryTree, prunePhylogeny = prunePhylogeny)
  
  checkPar(par)

  return(par)
  
}

