#' @title Phylogeny Reconstruction
#' @description reconstructs a phylogeny from a neutral trait matrix as returned by \code{\link{runSimulation}}
#' @param simu An object of class PhyloSim
#' @param which.result Integer, determines which result should be used. This argument is only useful if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used.
#' @param ...  Additional arguments to be passed to function \code{\link{hclust}}
#' @return An object of class 'phylo'
#' @example /inst/examples/reconstructPhylogeny-help.R
#' @export
#' 
phyloReconstruct <- function(simu, which.result = NULL, ...)
{
  if (is.null(which.result)) which.result = length(simu$Output) 
  simu <- simu$Output[[which.result]]
  
  speciesTable <- table(simu$specMat)
  traitMatrix <- simu$neutMat
  
  meanNeutralTraits <-numeric()  
  for(i in 1:length(speciesTable))
  {
    meanNeutralTraits[i] <- mean(traitMatrix[which(simu$specMat == names(speciesTable)[i])])
  }
  names(meanNeutralTraits) = names(speciesTable)
 
  distNeutral <- data.frame(x = names(speciesTable))
  columnX <- numeric()
  for(i in 1:length(meanNeutralTraits))
  {
    for(j in 1: length(meanNeutralTraits))
    {
      columnX[j] <- abs(meanNeutralTraits[i] - meanNeutralTraits[j] )
    }
    distNeutral <- cbind(distNeutral, columnX)
  }
  distNeutral <- distNeutral[,-1]
  colnames(distNeutral) <- paste("s",names(meanNeutralTraits),  sep = "")
  rownames(distNeutral) <- paste("s",names(meanNeutralTraits),  sep = "")
  
  distMat <- dist(distNeutral, diag = TRUE)

  phyloclust <- ape::as.phylo(stats::hclust(distMat, method = ...))
  return(phyloclust)
}


