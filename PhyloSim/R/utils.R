#' @title Get Recent Species
#' @description Function to extract the recent species of a Phylosim object
#' @param sim Phylosim object
#' @details The function extracts the IDs of the recent species of the grid in the Phylosim object.
#' @return List of vectors containing the species IDs.
#' @author Tankred Ott
#' @examples 
#' ## Define a parameter set
#' par <- createCompletePar(x = 50, y = 50, dispersal = FALSE , runs = c(500,1000),
#'         density = 1, environment = 0.5, specRate = 1)
#'
#' ## Run the model
#' simu <- runSimulation(par)
#' 
#' ## Extract species IDs
#' speciesIDs <- getSpecies(simu)
#' 
#' @export
getSpecies <- function (sim) {
  nRuns <- length(sim$Output)
  
  out <- list(rep(NA, nRuns))
  for (i in 1:nRuns) {
    out[[i]] <- unique(as.vector(sim$Output[[i]]$specMat))
  }
  return(out)
}

#' @title Get Species Trait
#' @description Function to extract the traits of a species from a Phylosim object.
#' @param sim Phylosim object
#' @param speciesID Integer, species ID
#' @param trait String, determining which trait should be returned ("evironment", "competition", "neutral").
#' @param whichRun Integer, specifying from which run the data should be collected.
#' @details The function extracts the environmental, competition or neutral trait of the species with speciesID from the Phylosim object. If there are multiple runs in whichRun determines which run should be used.
#' @return Vector, containing the trait values
#' @author Tankred Ott
#' @examples 
#' ## Define a parameter set
#' par <- createCompletePar(x = 50, y = 50, runs = c(500,1000),
#'         density = 1, environment = 0.5, specRate = 1)
#'
#' ## Run the model
#' simu <- runSimulation(par)
#' 
#' ## Extract species IDs
#' speciesIDs <- getSpecies(simu)
#' 
#' ## Get Species Traits of a species
#' envTraits <- getSpeciesTrait(simu, speciesIDs[[1]][1], "environment", 1)
#' 
#' hist(envTraits)
#' 
#' @export
getSpeciesTrait <- function (sim, speciesID, trait="neutral", whichRun = 1) {
  out <- NA
  specIndices <- getSpeciesIndices(sim, speciesID, whichRun)
  
  if (trait == "neutral") {
    out <- sim$Output[[whichRun]]$neutMat[specIndices]
  } else if (trait == "competition") {
    out <- sim$Output[[whichRun]]$compMat[specIndices]
  } else if (trait == "environment") {
    out <- sim$Output[[whichRun]]$traitMat[specIndices]
  } else {
    stop("Trait must be 'neutral', 'competition' or 'environment'")
  }
  
  return(out)
}

#' @title Get Species Traits
#' @description Function to extract the traits of a species from a Phylosim object.
#' @param sim Phylosim object
#' @param speciesID Integer, species ID
#' @param whichRun Integer, specifying from which run the data should be collected.
#' @details The function extracts the environmental, competition and neutral trait of the species with speciesID from the Phylosim object. If there are multiple runs in whichRun determines which run should be used.
#' @return List of trait types containing vectors with trait values.
#' @author Tankred Ott
#' @examples 
#' ## Define a parameter set
#' par <- createCompletePar(x = 50, y = 50, runs = c(500,1000),
#'         density = 1, environment = 0.5, specRate = 1)
#'
#' ## Run the model
#' simu <- runSimulation(par)
#' 
#' ## Extract species IDs
#' speciesIDs <- getSpecies(simu)
#' 
#' ## Get Species Traits of a species
#' traits <- getSpeciesTraits(simu, speciesIDs[[1]][1])
#' 
#' #hist(traits$environment)
#' 
#' @export
getSpeciesTraits <- function (sim, speciesID, whichRun = 1) {
  
  out <- list(environment=NA, competition=NA, neutral=NA)
  
  # specIndices <- which(sim$Output[[whichRun]]$specMat == speciesID)
  specIndices <- getSpeciesIndices(sim, speciesID, whichRun)
  
  # out$environment <- sim$Output[[whichRun]]$traitMat[specIndices]
  # out$competition <- sim$Output[[whichRun]]$compMat[specIndices]
  # out$neutral <- sim$Output[[whichRun]]$neutMat[specIndices]
  out$environment <- getSpeciesTrait(sim, speciesID, "environment", whichRun)
  out$competition <- getSpeciesTrait(sim, speciesID, "competition", whichRun)
  out$neutral <- getSpeciesTrait(sim, speciesID, "neutral", whichRun)

  return(out)
}


#' @title Get Species Indices
#' @description Function to get the indices (positions) of a species on the Phylosim grid
#' @param sim Phylosim object
#' @param speciesID Integer, species ID
#' @param whichRun Integer, specifying from which run the data should be collected.
#' @details The function extracts the environmental, competition and neutral trait of the species with speciesID from the Phylosim object. If there are multiple runs in 
#' @return Vector of indices
#' @author Tankred Ott
getSpeciesIndices <- function (sim, speciesID, whichRun=1) {
  return(which(sim$Output[[whichRun]]$specMat == speciesID))  
}


