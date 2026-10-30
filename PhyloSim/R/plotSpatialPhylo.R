#' @title Phylogeny And Spatial Abundance Plot
#' @description Plots the phylogeny tree, the spatial abundance of the evolved species and their traits.
#' @param simu Simulation output of the class "Phylosim", usually consisting out of several lists. Needs to contain at least the spatial species matrix ($specMat) and the phylogeny ($phylogeny)   
#' @param plot defines what to plot."spatial" plots the spatial abundance, "phylogeny" the phylogeny tree. "both" plots both of them side by side 
#' @param plotTraits logical, defines if the traits should be plotted 
#' @param which.result Integer, determines which result should be used. This argument is only useful if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used.
#' @details The tiplabels of the phylogeny tree are colored in correspondency to the spatial abundance plot
#' @example /inst/examples/plotSpatialPhylo-help.R
#' @import ape
#' 
#' @export
plotSpatialPhylo <- function (simu, plot = "both", plotTraits = T, which.result = NULL){
  
  oldpar<-par()
  
  if (is.null(which.result)) which.result = length(simu$Output) 
  simu <- simu$Output[[which.result]]
  
  phylogeny = simu$phylogeny
  
  extantPhylogeny <- ape::drop.fossil(phylogeny)
  
  nSpecies = length(extantPhylogeny$tip.label)
  
  landscape = simu$specMat
  
  if (plotTraits == T){
    #library(shape)
    traits = getAverageTraits(simu)  
    #rootDist <- adephylo::distRoot(extantPhylogeny, 1)
  }
  
  cols <- rainbow(nSpecies,start=0)
  
  tipNumbers <- as.numeric(sapply(list(extantPhylogeny$tip.label), substring, 2))
  
  landscape <- matrix(as.numeric(factor(landscape, levels = tipNumbers, labels = 1:nSpecies)), nrow = nrow(landscape))
  
  if(plot!="both"){
    at <- 0.5
    if(plotTraits==T){
      layout(mat = matrix(1:4, nrow = 2, byrow = T), widths = c(0.5,1), heights = c(0.15,0.8))
      par(mar = c(0,0,0,0), oma = c(2,2,2,2))
      plot(NULL,frame = F, yaxt='n', xaxt='n', ann=FALSE, xlim = c(0,7), ylim = c(0,3))
      text(1:4 - 0.5, 0, colnames(traits), srt = 45, pos = 4)
      plot.new()
      size = dev.size()[2]/3
      plot(rep(1:4, each = nSpecies), rep(1:nSpecies, 4), cex = size*traits, pch = 16, frame = F, yaxt='n', xaxt='n', ann=FALSE, xlim = c(0,7))
      mtext("Traits", at=0.1, outer=T, cex=1.5)
      at <- 0.6
      
    }
    
    if(plotTraits==F){
      par(mfrow=c(1,1),mar = c(1,1,1,1))
    }
    if(plot=="spatial"){
      image(landscape, col = cols, yaxt='n', xaxt='n', asp=1, bty="n")
      mtext(outer =T ,"Spatial distribution", cex=1.2, at =at)
      
    }
    else{
      if(plot=="phylogeny"){
        plot(extantPhylogeny, tip.color = cols, cex = 1.3)
        mtext(outer =T ,"Phylogeny", cex=1.2, at = at)
        
        
      }
      else stop("wrong plot argument")
    }
  }
  else{
    if(plotTraits==T){
      layout(mat = matrix(1:6, nrow = 2, byrow = T), widths = c(1,0.5,1), heights = c(0.15,0.8))
      par(mar = c(0,0,0,0), oma = c(2,2,2,2))
      plot.new()
      plot(NULL,frame = F, yaxt='n', xaxt='n', ann=FALSE, xlim = c(0,7), ylim = c(0,3))
      text(1:4 - 0.5, 0, colnames(traits), srt = 45, pos = 4)
      plot.new()
      plot(extantPhylogeny, tip.color = cols, cex = 1.3)
      size = dev.size()[2]/3
      plot(rep(1:4, each = nSpecies), rep(1:nSpecies, 4), cex = size*traits, pch = 16, frame = F, yaxt='n', xaxt='n', ann=FALSE, xlim = c(0,7))
      par(mar = c(1,1,1,1))
      image(landscape, col = cols, yaxt='n', xaxt='n', bty="n", asp=1)
      mtext("Traits", at=0.48, outer=T, cex=1.2)
      mtext(outer =T ,"Spatial distribution",at =0.8, cex=1.2)
      mtext(outer =T ,"Phylogeny",at =0.2, cex=1.2)
    }
    else{
      par(mfrow=c(1,2))
      plot(extantPhylogeny, tip.color = cols, cex = 1.3)
      image(landscape, col = cols, yaxt='n', xaxt='n')
      mtext(outer =T ,"Spatial distribution",at =0.75, cex=1.2)
      mtext(outer =T ,"Phylogeny",at =0.25, cex=1.2)
    }
  }
  par(mfrow=oldpar$mfrow, mar=oldpar$mar, oma=oldpar$oma)
}


