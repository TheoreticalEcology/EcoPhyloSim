#' @title Species richness
#' @description Calculate the total species richness of a community
#' @param simu Object of class "PhyloSim", usually consisting of several lists. Needs at least the spatial distribution of the species stored in a matrix ($specMat)
#' @param which.result Integer, determines which result should be used. This argument is only useful if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used.
#' @return An integer value for the species richness
#' @export

specRich <-function(simu,which.result=NULL){
  
  if (is.null(which.result)) which.result = length(simu$Output) 
  simu <- simu$Output[[which.result]]
  
  matrix <- simu$specMat
  
  sr <- length(unique(c(matrix)))
  return(speciesRichness = sr)
}

