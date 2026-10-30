
#' @title Plot all objects in an object of class PhylosimList
#' @description Plots the phylogeny tree ,the spatial abundance of the evolved species and their traits.
#' @param x Simulation output of \code{\link{runSimulationBatch}}
#' @param ... further arguments, currently not used

#' @details The tiplabels of the phylogeny tree are colored in correspondency to the spatial abundance plot

#' @export

plot.PhylosimList <- function(x, ...){
  for(i in 1:(length(x)-1)){
    plotSpatialPhylo(x[[i]])
    cat ("Press [enter] to continue")
    line <- readline()
  } 
  plotSpatialPhylo(x[[length(x)]])
}


#print.phylosimlist 
#' @title Prints summary for objects of class PhylosimList
#' @description Prints parameter settings, runtime and species Richness.
#' @param x Simulation output of \code{\link{runSimulation}}
#' @param ... further arguments, currently not used
#' @export

print.PhylosimList <- function(x, ...){
  cat("Summary for class PhylosimList",  "\n")
  cat("------------------------------------", "\n","------------------------------------", "\n", sep="")
  cat("# Simu", "X ", "Y", "Dispersal", "Density", "Environment", "#Species","#Runs","\n")
  cat("------------------------------------", "\n")
  cat()
  for(i in 1:length(x)){
    cat("#",i," ",x[[i]]$Model$x,"", x[[i]]$Model$y," ",x[[i]]$Model$dispersal,"     ",
        as.numeric(x[[i]]$Model$density),"      ", 
        as.numeric(x[[i]]$Model$environment),"      ", specRich(x[[i]]),"   ", x[[i]]$Model$runs[length(x[[i]]$Model$runs)], "\n")
  }
}


#' @title Summary for objects of class PhylosimList
#' @description Gives summary about parameter settings, runtime and species Richness.
#' @param object Simulation output of \code{\link{runSimulationBatch}}
#' @param ... further arguments, currently not used
#' @export

summary.PhylosimList <- function(object, ...){
  cat("Summary for class PhylosimList",  "\n")
  cat("------------------------------------", "\n","------------------------------------", "\n", sep="")
  cat("# Simu", "X ", "Y", "Dispersal", "Density", "Environment","#Species","#Runs","Colless'","Gamma","\n")
  cat("------------------------------------", "\n")
  cat()
  for(i in 1:length(object)){
    cat("#",i," ",object[[i]]$Model$x,"", object[[i]]$Model$y," ",object[[i]]$Model$dispersal,"     ",
        as.numeric(object[[i]]$Model$density),"      ", 
        as.numeric(object[[i]]$Model$environment),"      ", specRich(object[[i]]),"   ", 
        object[[i]]$Model$runs[length(object[[i]]$Model$runs)],
        "  ", round(collessImbalance(object[[i]]),3), "  ",
        suppressWarnings(round(ape::gammaStat(ape::drop.fossil(object[[i]]$Output[[length(object[[i]]$Output)]]$phylogeny)),3)),"\n")
  }
}