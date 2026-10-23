#' @title Phylogenetic clades
#' @description Visualises the phylogenetic clades of a local community within a metacommunity. 
#' @param simu Simulation output of the class "Phylosim", usually consisting out of several lists. Needs to contain at least the phylogeny ($phylogeny) and species matrix ($specMat). The species matrix is used as metacommunity
#' @param which.result Integer, determines which result should be used. This argument is only usefull if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used.
#' @param size size of the local plots.
#' @param n The number of local communities. For each community a new plot is created. Needs to be greater than 1
#' @param subplot if no simu object is given, a communitytable of a local community can be given here
#' @param phylogeny if no simu object is given, an object of the class 'phylo' can be given here to describe the community's phylogeny
#' @details This function works with a single simu object (class "phylosim"). If no simu object is given, it is also possible to calculate the clades from a given phylogeny and local plots. Here the phylogeny represents the metacommunity within which the local plots create the clades
#' @example /inst/examples/plotHighlightLocalPhylo-help.R
#' @export
highlightLocalPhylo <- function(simu,which.result=NULL, size, n, subplot = NULL, phylogeny=NULL){
  
  
  if(is.null(subplot)){
  lp <- localPlots(simu = simu,size = size, n=n,community=T , which.result =which.result)
  community <- lp$communityTable} else community <- subplot
  
  if(is.null(phylogeny)){
    if (is.null(which.result)) which.result = length(simu$Output) 
    simu <- simu$Output[[which.result]]
  
  phylo <- simu$phylogeny} else phylo <- phylogeny
  phylo <- ape::drop.fossil(phylo)
  
  for( i in 1:n){
    x <- phylo$tip.label[is.element(phylo$tip.label,names(community[which(community[i,] > 0)]))]
    y <- phylo$tip.label[!is.element(phylo$tip.label,names(community[which(community[i,] > 0)]))]
    xy <- list(x,y)
    
    group <- ape::which.edge(phylo, x)
    
    cladecolor<-rep("darkgrey", dim(phylo$edge)[1])
    cladecolor[group]<-"black"
    
    cladewidth<-rep(1, dim(phylo$edge)[1])
    cladewidth[group]<-2
    
    title <- paste("Local Plots (size =",size,")","Within Metacommunity \n Local Community", i,"/",n)
    plot(phylo,  direction="downwards", show.tip.label=F, show.node.label=F, edge.color=cladecolor, edge.width=cladewidth, main = title)
    
    place<- numeric(length(x))
    for(k in 1:length(x)){
      place[k] <- which(phylo$tip.label==x[k])
    }
    
    ape::tiplabels(text = x,tip = place, frame ="none", srt = 270, adj=0)
    ape::add.scale.bar()
  }
}



