#' @title Null Model
#' @description Generates a null model and compares the observed results against it.
#' @param simu an object of class "PhyloSim" as created by \code{\link{runSimulation}} or \code{\link{runSimulationBatch}}
#' @param which.result Integer, determines which result should be used. This argument is only useful if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used.
#' @param abundance Logical, determining whether the null model should be based on the abundance distribution of the sample plots. Default is FALSE.
#' @param localPlotSize edge length (in grid cells) of the square plots drawn within the metacommunity
#' @param numberOfPlots number of plots drawn within the metacommunity
#' @param repetitions number of generated null model plots
#' @param fun String, phylogenetic measure that is used. Possible inputs are "mpd" (\code{\link[picante]{mpd}}) and "pd" (\code{\link[picante]{pd}}). 
#' @return A numeric vector with pValues for each plot in the observed metacommunity
#' @example /inst/examples/nullModel-help.R
 
#' @export
nullModel <- function(simu, which.result=NULL, abundance = FALSE, localPlotSize, numberOfPlots, repetitions, fun="mpd"){
  
  comMat <- localPlots(simu=simu, which.result=which.result, size = localPlotSize, n = numberOfPlots, community = T)$communityTable # create community matrix from local communities (requires PhyloSim)
  
  
  if (is.null(which.result)) which.result = length(simu$Output) 
  simu <- simu$Output[[which.result]]
  phylogeny <- simu$phylogeny
  speciesMatrix <- simu$specMat
  
  extantPhyloCophen <- stats::cophenetic(phylogeny) # create distance matrix from phylogeny (requires picante)
  
  # calculate fun for "real matrix"
  ifelse (fun == "mpd", observedMPD <- picante::mpd(samp = comMat, dis = extantPhyloCophen, abundance.weighted = T),
          ifelse(fun == "pd", observedMPD <- picante::pd(samp = comMat, phylogeny)$PD, stop("wrong dist argument")))
  
  
  # calculate abundances of species in metacommunity
  
  ifelse(abundance==FALSE, metaabundance <- table(speciesMatrix),
         ifelse(abundance==TRUE, metaabundance <- colSums(comMat), stop("wrong abundance argument") ))
  
  
  if(abundance==FALSE) names(metaabundance) <- paste("s", names(metaabundance), sep="") # set names equal to those in the local communities and the distance matrix
  
  
  distribution <- rmultinom(n = repetitions, size=localPlotSize,  prob=metaabundance) # create multinomial distributed matrix
  
  if(fun == "mpd"){
  NullDistribution <- picante::mpd(samp = t(distribution), dis = extantPhyloCophen, abundance.weighted = T) # calculate mpd for null matrix
  }
  else{
   NullDistribution<- picante::pd(samp=t(distribution), tree= phylogeny)$PD
  }
   
  observedMPD[which(is.na(observedMPD))] <- 0
  
  CummulativeDensityFunction <- ecdf(NullDistribution) # create cumulative density function from null matrix
  
  
  pValues <- CummulativeDensityFunction(observedMPD) # use  cumulative density function on "real" mpd
  
  return(pValues)
}

