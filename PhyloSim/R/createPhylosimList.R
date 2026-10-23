#' @title PhylosimList creator
#' @description Creates an object of class "PhylosimList" from a list of "Phylosim" objects.
#' @param simu List of objects of type Phylosim as created by \code{\link{runSimulation}}
#' @return An object of type PhylosimList
#' @details The function is used to make functions for objets of type "PhylosimList" usable for objects of class "Phylosim". The function needs at least two objects to work.  
#' @example /inst/examples/createPhylosimList-help.R
#' @export


createPhylosimList<-function(simu){
  simuBatch<-simu
  if(length(simuBatch) < 2) stop("Invalid argument. List of objects too short.")
  class(simuBatch)<-"PhylosimList"
  return(simuBatch)
}

