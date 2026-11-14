#' @title Trait plots
#' @description Plots trait-histograms, trait ~ Environment relationship and the spatial distribution of the different traits
#' @param simu Simulation output of the class "PhyloSim", usually consisting of several lists. Needs to contain at least the three trait matrices ($traitMat;$compMat;$neutMat) and the environment matrix ($envMat) 
#' @param which.result Integer, determines which result should be used. This argument is only useful if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used.
#' @param type defines whether the histograms should be plotted standalone or with the trait ~ Environment relationship and the spatial distribution. The latter is done by type = "all". Default is "hist"
#' @example /inst/examples/plotTraitDistribution-help.R
#' 

#' @export

plotTraitDistribution <- function (simu, type = "hist", which.result = NULL){
  
  oldpar<-par()
  
  if (is.null(which.result)) which.result = length(simu$Output) 
  dat <- simu$Output[[which.result]]
  
  names <- unique(as.numeric(dat$specMat))
  # same species colors as in plotSpatialPhylo(), without random numbers (#156)
  cols <- unname(speciesColors(ape::drop.fossil(dat$phylogeny), alpha = 0.5)[paste0("s", names)])
  cols[is.na(cols)] <- grey(0.5, alpha = 0.5)
  ymax <- max(1, log(max(table(dat$specMat)))) # no histogram bin can exceed the largest abundance
  colmat <- dat$specMat
  for(i in 1:length(names)){
    colmat[colmat==names[i]] <- cols[i]}
  
  if(type == "all"){
    par(mfrow=c(3,3))
    par(mar = c(5, 4, 2, 2))
    
    plot(dat$traitMat~dat$envMat, col = colmat, ylab="Trait", xlab="Environment"  )
    title("Environmental Trait")
    plot(dat$compMat~dat$envMat, col = colmat, ylab= "", xlab="Environment")
    title("Competition Trait")
    plot(dat$neutMat~dat$envMat, col = colmat, ylab= "", xlab="Environment")
    title("Neutral Trait")
    
    image(dat$traitMat, xaxt="n", yaxt="n", ylab="Spatial Distribution",useRaster =T,col = grey(seq(0, 1, length = 256)), asp=1, bty="n")
    image(dat$compMat, xaxt="n", yaxt="n",useRaster =T,col = grey(seq(0, 1, length = 256)), asp=1, bty="n")
    image(dat$neutMat, xaxt="n", yaxt="n",useRaster =T,col = grey(seq(0, 1, length = 256)), asp=1, bty="n")
    
    }else{par(mfrow=c(1,3))}
   
  xmax <- max(dat$traitMat) 
  xmin <- min(dat$traitMat)
  
  f <- hist(dat$traitMat[dat$specMat==names[1]], plot=F)
  f$counts <- log(f$counts)
  plot(f,ylim = c(0, ymax), xlim=c(xmin,xmax), col =cols[1], ylab="log(Frequency)" , xlab="Value", main = "Environment Histogram", border =0)
  
  for(i in 2:length(names)){
    f <- hist(dat$traitMat[dat$specMat==names[i]], plot=F)
    f$counts <- log(f$counts)
    plot(f, add=T, col =cols[i], border =0)
  }
  
  
  xmax <- max(dat$compMat) 
  xmin <- min(dat$compMat)
  
  f <- hist(dat$compMat[dat$specMat==names[1]], plot=F)
  f$counts <- log(f$counts)
  plot(f,ylim = c(0, ymax), xlim=c(xmin,xmax), col =cols[1], ylab="log(Frequency)" , xlab="Value", main = "Competition Histogram", border =0)
  
  for(i in 2:length(names)){
    f <- hist(dat$compMat[dat$specMat==names[i]], plot=F)
    f$counts <- log(f$counts)
    plot(f, add=T, col =cols[i], border =0)
  }
  
  
  xmax <- max(dat$neutMat) 
  xmin <- min(dat$neutMat)
  
  f <- hist(dat$neutMat[dat$specMat==names[1]], plot=F)
  f$counts <- log(f$counts)
  plot(f,ylim = c(0, ymax), xlim=c(xmin,xmax), col =cols[1], ylab="log(Frequency)" , xlab="Value", main = "Neutral Histogram", border =0)
  
  for(i in 2:length(names)){
    f <- hist(dat$neutMat[dat$specMat==names[i]], plot=F)
    f$counts <- log(f$counts)
    plot(f, add=T, col =cols[i], border =0)
    
  }
  par(mfrow=oldpar$mfrow, mar=oldpar$mar)
}


