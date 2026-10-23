## Class functions for class Phylosim


#plot.phylosim

#' @title Plot objects of class PhyloSim
#' @description Plots the phylogeny tree ,the spatial abundance of the evolved species and their traits.
#' @param x Simulation output of \code{\link{runSimulation}}
#' @param ... further arguments, currently not used

#' @details The tiplabels of the phylogeny tree are colored in correspondency to the spatial abundance plot

#' @export

plot.PhyloSim <- function(x, ...){
  plotSpatialPhylo(x)
}


#summary.phylosim

#' @title Summary for objects of class PhyloSim
#' @description Gives summary about parameter settings, runtime and species Richness.
#' @param object Simulation output of \code{\link{runSimulation}}
#' @param ... further arguments, currently not used
#' @export

summary.PhyloSim <- function(object, ...){
  
  if(length(object$Model$runs)==1){
    runsarg<-object$Model$runs
  } else{
    runsarg<- "c("
    for(i in 1:(length(object$Model$runs)-1)){
      runsarg<- paste(runsarg,object$Model$runs[[i]],",")
    }
    runsarg<- paste(runsarg, object$Model$runs[[length(object$Model$runs)]], ")" )
  }
  
  cat("Summary for class Phylosim",  "\n")
  cat("# Call: runSimulation(x = ", object$Model$x, ", y = ", object$Model$y, ", dispersal = '", 
      object$Model$dispersal,"'", ", runs =", runsarg,"\n", ", density = ", object$Model$density, 
      ", environment = ", object$Model$environment, ", specRate = ",object$Model$specRate,")", "\n", sep="")
  cat("------------------------------------", "\n", sep="")
  cat("Number of Species in Metacommunity: ", specRich(object), "\n", sep="")
  cat("Colless' Imbalance: ", collessImbalance(object), "\n", sep="")
  cat("Gamma: ", suppressWarnings(ape::gammaStat(ape::drop.fossil(object$Output[[length(object$Output)]]$phylogeny))), "\n" )
  cat("Time elapsed: ", object$Model$runtime,"seconds", "\n", sep="")
}



#print.phylosim 
#' @title Prints summary for objects of class PhyloSim
#' @description Prints parameter settings, runtime and species Richness.
#' @param x Simulation output of \code{\link{runSimulation}}
#' @param ... further arguments, currently not used
#' @export

print.PhyloSim <- function(x, ...){
  
  if(length(x$Model$runs)==1){
    runsarg<-x$Model$runs
  } else{
    runsarg<- "c("
    for(i in 1:(length(x$Model$runs)-1)){
      runsarg<- paste(runsarg,x$Model$runs[[i]],",")
    }
    runsarg<- paste(runsarg, x$Model$runs[[length(x$Model$runs)]], ")" )
  }
  
  cat("Summary for class PhyloSim",  "\n")
  cat("# Call: runSimulation(x = ", x$Model$x, ", y = ", x$Model$y, ", dispersal = '", 
      x$Model$dispersal,"'", ", runs =", runsarg, ", density = ", x$Model$density, 
      ", environment = ", x$Model$environment, ", specRate = ",x$Model$specRate,")", "\n", sep="")
  cat("------------------------------------", "\n", sep="")
  cat("Time elapsed: ", x$Model$runtime,"seconds", "\n", sep="")
  
}