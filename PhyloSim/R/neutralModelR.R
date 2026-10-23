#' @title Model of neutral community assembly
#' @description Runs a simple model of species community assembly under neutral conditions as defined by Hubbell (2001). Usually called by \code{\link{runSimulation}}.
#' @param xdim Integer, Dimension of the model landscape in x-direction
#' @param ydim Integer, Dimension of the model landscape in y-direction
#' @param specRate The speciation rate in total speciation events per generation
#' @param seed The random seed (do not change this option if you previously set the seed!)
#' @param runs The number of generations the model runs through
#' @references Hubbell, Stephen P. The unified neutral theory of biodiversity and biogeography (MPB-32). Vol. 32. Princeton University Press, 2001.
#' @details This model is to be seen only for test and teaching purpose. However, if you want to use the results for further analysis please call the model with \code{\link{runSimulation}}. \cr\cr  Be careful with the dimensions you choose. Large grids and many generations may take very long to compute. You can assume to have reached the equilibrium state of the model after dim^2/2.  
#' @example /inst/examples/neutralModelR-help.R

#' @export
NeutralMod <- function(  xdim = 100, ydim=100, specRate = 2,  seed = NULL, runs = 500 ){
  
  ptm <- proc.time()
  if (! is.null(seed)) set.seed(seed)
  
  counter <- 1
  mat  <- matrix(counter,xdim, ydim)
  
  for(i in 1:runs){
    
    for (j in 1:(xdim*ydim)){
      coordx = sample.int(n=xdim,2, T )
      coordy = sample.int(n=ydim,2, T )
      mat[coordx[1], coordy[1]] = mat[coordx[2], coordy[2]]
    }
    for (j in 1:specRate){
      SpecCoordx = sample.int(n=xdim,1, T )
      SpecCoordy = sample.int(n=ydim,1, T )
      mat[SpecCoordx, SpecCoordy] <- counter +1
      counter <- counter +1  
    }
    if(i %% 100 == 0) print(i)
  }
  print (paste("Finished after",floor(((proc.time() - ptm)[3])/60), "minute(s) and", ((proc.time() - ptm)[3])%%60, "second(s)."))
  return (mat)
}
