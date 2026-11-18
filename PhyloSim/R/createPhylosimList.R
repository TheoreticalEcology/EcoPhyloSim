#' @title PhyloSimList creator
#' @description Creates an object of class "PhyloSimList" from a list of "PhyloSim" objects.
#' @param simu List of objects of type PhyloSim as created by \code{\link{runSimulation}}
#' @return An object of class \code{\link{PhyloSimList}}
#' @details The function is used to make functions for objects of type "PhyloSimList" usable for objects of class "PhyloSim". The function needs at least two objects to work.  
#' @example /inst/examples/createPhylosimList-help.R
#' @export


createPhyloSimList<-function(simu){
  simuBatch<-simu
  if(length(simuBatch) < 2) stop("Invalid argument. List of objects too short.")
  class(simuBatch)<-"PhyloSimList"
  return(simuBatch)
}


#' @rdname createPhyloSimList
#' @details \code{createPhylosimList()} is the deprecated name of \code{createPhyloSimList()}.
#' @export
createPhylosimList <- function(simu){
  .Deprecated("createPhyloSimList")
  createPhyloSimList(simu)
}
