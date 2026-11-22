## Class functions for class Phylosim


#plot.phylosim

#' @title Plot objects of class PhyloSim
#' @description Plots the phylogeny tree, the spatial abundance of the evolved species and their traits.
#' @param x Simulation output of \code{\link{runSimulation}}
#' @param ... further arguments, currently not used

#' @details The tip labels of the phylogeny tree are colored according to the spatial abundance plot.

#' @export

plot.PhyloSim <- function(x, ...){
  plotSpatialPhylo(x)
}


#summary.phylosim

#' @title Summary for objects of class PhyloSim
#' @description Gives a summary of parameter settings, runtime and species richness.
#' @param object Simulation output of \code{\link{runSimulation}}
#' @param ... further arguments, currently not used
#' @export

summary.PhyloSim <- function(object, ...){
  
  gens <- modelGenerations(object$Model)
  if(length(gens)==1){
    runsarg<-gens
  } else{
    runsarg<- "c("
    for(i in 1:(length(gens)-1)){
      runsarg<- paste(runsarg,gens[[i]],",")
    }
    runsarg<- paste(runsarg, gens[[length(gens)]], ")" )
  }
  
  cat("Summary for class PhyloSim",  "\n")
  cat("# Call: runSimulation(x = ", object$Model$x, ", y = ", object$Model$y, ", dispersal = '", 
      object$Model$dispersal,"'", ", generations = ", runsarg, ", negativeDD = ", c(object$Model$negativeDD, object$Model$density)[1], 
      ", environment = ", object$Model$environment, ", specRate = ",object$Model$specRate,")", "\n", sep="")
  cat("------------------------------------", "\n", sep="")
  cat("Number of Species in Metacommunity: ", specRich(object), "\n", sep="")
  cat("Colless' Imbalance: ", collessImbalance(object), "\n", sep="")
  cat("Gamma: ", .gammaStat(object$Output[[length(object$Output)]]$phylogeny), "\n" )
  cat("Time elapsed: ", object$Model$runtime,"seconds", "\n", sep="")
}



#print.phylosim 
#' @title Prints summary for objects of class PhyloSim
#' @description Prints parameter settings and runtime.
#' @param x Simulation output of \code{\link{runSimulation}}
#' @param ... further arguments, currently not used
#' @export

print.PhyloSim <- function(x, ...){
  
  gens <- modelGenerations(x$Model)
  if(length(gens)==1){
    runsarg<-gens
  } else{
    runsarg<- "c("
    for(i in 1:(length(gens)-1)){
      runsarg<- paste(runsarg,gens[[i]],",")
    }
    runsarg<- paste(runsarg, gens[[length(gens)]], ")" )
  }
  
  cat("Summary for class PhyloSim",  "\n")
  cat("# Call: runSimulation(x = ", x$Model$x, ", y = ", x$Model$y, ", dispersal = '", 
      x$Model$dispersal,"'", ", generations =", runsarg, ", negativeDD = ", c(x$Model$negativeDD, x$Model$density)[1], 
      ", environment = ", x$Model$environment, ", specRate = ",x$Model$specRate,")", "\n", sep="")
  cat("------------------------------------", "\n", sep="")
  cat("Time elapsed: ", x$Model$runtime,"seconds", "\n", sep="")
  
}
