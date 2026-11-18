#' @title Class PhyloSimList
#' @name PhyloSimList
#' @aliases PhylosimList
#' @description An object of class "PhyloSimList" is a list of simulation outputs of class "PhyloSim" (see \code{\link{runSimulation}}), one for each parameter set. It is returned by \code{\link{runSimulationBatch}} and created from a list of "PhyloSim" objects by \code{\link{createPhyloSimList}}. Before version 0.3.2, the class was called "PhylosimList"; such objects keep working. Its elements are accessed by indexing, e.g. \code{simuList[[2]]}, and can be passed to all functions for "PhyloSim" objects; with \code{onError = "NA"} in \code{\link{runSimulationBatch}}, failed simulations are NA.
#' @details Methods for "PhyloSimList" objects: \code{\link[=plot.PhyloSimList]{plot}}, \code{\link[=print.PhyloSimList]{print}} and \code{\link[=summary.PhyloSimList]{summary}}.
NULL

#' @title Plot all objects in an object of class PhyloSimList
#' @description Plots the phylogeny tree, the spatial abundance of the evolved species and their traits.
#' @param x Simulation output of \code{\link{runSimulationBatch}}
#' @param ... further arguments, currently not used

#' @details The tip labels of the phylogeny tree are colored according to the spatial abundance plot.

#' @export

plot.PhyloSimList <- function(x, ...){
  for(i in 1:(length(x)-1)){
    plotSpatialPhylo(x[[i]])
    cat ("Press [enter] to continue")
    line <- readline()
  } 
  plotSpatialPhylo(x[[length(x)]])
}


#print.phylosimlist 
#' @title Prints summary for objects of class PhyloSimList
#' @description Prints parameter settings, species richness and number of generations for each simulation.
#' @param x Simulation output of \code{\link{runSimulationBatch}}
#' @param ... further arguments, currently not used
#' @export

print.PhyloSimList <- function(x, ...){
  cat("Summary for class PhyloSimList",  "\n")
  cat("------------------------------------", "\n","------------------------------------", "\n", sep="")
  cat("# Simu", "X ", "Y", "Dispersal", "Density", "Environment", "#Species","#Generations","\n")
  cat("------------------------------------", "\n")
  cat()
  for(i in 1:length(x)){
    cat("#",i," ",x[[i]]$Model$x,"", x[[i]]$Model$y," ",x[[i]]$Model$dispersal,"     ",
        as.numeric(c(x[[i]]$Model$negativeDD, x[[i]]$Model$density)[1]),"      ", 
        as.numeric(x[[i]]$Model$environment),"      ", specRich(x[[i]]),"   ", tail(modelGenerations(x[[i]]$Model), 1), "\n")
  }
}


#' @title Summary for objects of class PhyloSimList
#' @description Gives a summary of parameter settings, species richness, Colless' imbalance and Gamma for each simulation.
#' @param object Simulation output of \code{\link{runSimulationBatch}}
#' @param ... further arguments, currently not used
#' @export

summary.PhyloSimList <- function(object, ...){
  cat("Summary for class PhyloSimList",  "\n")
  cat("------------------------------------", "\n","------------------------------------", "\n", sep="")
  cat("# Simu", "X ", "Y", "Dispersal", "Density", "Environment","#Species","#Generations","Colless'","Gamma","\n")
  cat("------------------------------------", "\n")
  cat()
  for(i in 1:length(object)){
    cat("#",i," ",object[[i]]$Model$x,"", object[[i]]$Model$y," ",object[[i]]$Model$dispersal,"     ",
        as.numeric(c(object[[i]]$Model$negativeDD, object[[i]]$Model$density)[1]),"      ", 
        as.numeric(object[[i]]$Model$environment),"      ", specRich(object[[i]]),"   ", 
        tail(modelGenerations(object[[i]]$Model), 1),
        "  ", round(collessImbalance(object[[i]]),3), "  ",
        round(.gammaStat(object[[i]]$Output[[length(object[[i]]$Output)]]$phylogeny),3),"\n")
  }
}

# objects created before version 0.3.2 have class PhylosimList (#92)
#' @rawNamespace S3method(plot, PhylosimList, plot.PhyloSimList)
#' @rawNamespace S3method(print, PhylosimList, print.PhyloSimList)
#' @rawNamespace S3method(summary, PhylosimList, summary.PhyloSimList)
NULL
