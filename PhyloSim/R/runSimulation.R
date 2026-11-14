#' @title  Species Community Simulation
#' @description A model of species community assembly under different assembly mechanisms.  
#' @param par a list of model parameters created with \link{createCompletePar}
#' @return An object of class \code{PhyloSim}. This object contains the species matrix, the trait matrix, the environmental matrix, the competition matrix and the neutral matrix, as well as the phylogeny and the parameter settings of the simulation. 
#' @details If the runs argument of your parameter set contains more than one element, each interim step is saved in the PhyloSim object. \cr\cr For larger simulations consider \link{runSimulationBatch} to make use of parallel computing. \cr\cr It is possible that more than one new species arises per generation. This leads to a multifurcated phylogeny, yet many packages such as "ape" can only work with bifurcating trees. If `convertToBinaryTree = TRUE` (the default), the phylogeny is converted to a bifurcating one using `multi2di()` from the "ape" package. \cr\cr Branch lengths of the phylogeny are in generations. Speciation and extinction events of generation g are placed at time g, and the census of generation g (an entry in \code{runs}) at time g + 0.5, after all events of that generation. Tips of extant species therefore end half a generation after those of species that went extinct in the last generation, so that \code{ape::drop.fossil()} removes the latter. Internal nodes are labeled with the species from which the descendants emerged (nodes added by \code{convertToBinaryTree} have no label); the tree has no root edge. With \code{prunePhylogeny = TRUE} (the default), the phylogeny contains only the extant species; with \code{prunePhylogeny = FALSE}, it also contains all extinct species as tips. 
#' @importFrom adephylo distTips
#' @import Rcpp
#' @useDynLib PhyloSim, .registration = TRUE
#' @export
#' @example /inst/examples/runSimulation-help.R


runSimulation <- function(par)
{
  checkPar(par)
  
  # GENERAL SETTINGS
  
  if (par$type == "base"){
    
    if (par$dispersal == "global" | par$dispersal == 0){
      dispersal = 1
      dispersalCut = 1
    }else if (is.numeric(par$dispersal)){
      dispersalCut = 2*par$dispersal
      dispersal = 3
    }else stop("wrong dispersal parameter") 
    
    if(par$fitnessActsOn == "mortality"){
      reproductiveFitness = F
      mortalityFitness = T    
    }else if(par$fitnessActsOn == "reproduction"){
      reproductiveFitness = T
      mortalityFitness = F    
    }else if (par$fitnessActsOn == "both"){
      reproductiveFitness = T
      mortalityFitness = T    
    }else stop("wrong fitness parameters")
    
    
    # density and environment stay numeric, so that the stored $Model passes checkPar() when rerun (#118)
    par$compStrength <- if (par$density == 0) 1 else par$density
    par$envStrength <- if (par$environment == 0) 1 else par$environment
    neutral <- par$density == 0 && par$environment == 0
    
    ptm <- proc.time()
    
    out <- callModel( x = par$x, 
                      y = par$y, 
                      dispersal = dispersal, 
                      runs = round(par$runs), 
                      specRate = par$specRate, 
                      dens = par$density > 0, 
                      env = par$environment > 0, 
                      neutral = neutral, 
                      mort = mortalityFitness, 
                      mortStrength = par$fitnessBaseMortalityRatio, 
                      repro = reproductiveFitness, 
                      dispersalCutoff = dispersalCut, 
                      densityCutoff = par$densityCut, 
                      seed = par$seed, 
                      envStrength = par$envStrength, 
                      compStrength = par$compStrength, 
                      fission = par$fission, 
                      redQueen = par$redQueen, 
                      redQueenStrength = par$redQueenStrength, 
                      protracted = par$protracted, 
                      airmatR = par[["airmatR"]], 
                      soilmatR = par[["soilmatR"]],
                      prunePhylogeny = par$prunePhylogeny,
                      nicheWidth = par$nicheWidth)
    
    runtime <- as.numeric((proc.time() - ptm)[3])
    
    print (paste("Core simulation finished after",floor(((proc.time() - ptm)[3])/60), "minute(s) and", ((proc.time() - ptm)[3])%%60, "second(s). Converting data"))
    
    #     Rcpp::List listResults = Rcpp::List::create(Rcpp::Named("Species") = specOut,
    #                                                 Rcpp::Named("EnvTrait") = traitOut,
    #                                                 Rcpp::Named("NeutralTrait") = neutralOut,
    #                                                 Rcpp::Named("CompetitionTrait") = compOut,
    #                                                 Rcpp::Named("Environment") = envOut,
    #                                                 Rcpp::Named("Phylogeny") = phyloOut);
    
    # The following code reshuffles the output to the structure 
    # of a phylosim class
    output<-list()
    for (i in 1:length(par$runs)){
      
      specMati = matrix(out[[i]]$Species,ncol=par$x, nrow=par$y)
      if(length(unique(c(specMati)))==1){
        phyloi <- 0
        warning("Cannot build Phylogeny")
      }else{
        phyloi <- ape::read.tree(text= out[[i]]$Phylogeny)
      }
       
      output$Output[[i]] = list(
        specMat = matrix(out[[i]]$Species,ncol=par$x, nrow=par$y), 
        traitMat= matrix(out[[i]]$EnvTrait,ncol=par$x, nrow=par$y), 
        envMat = matrix(out[[i]]$Environment,ncol=par$x, nrow=par$y), 
        compMat = matrix(out[[i]]$CompetitionTrait,ncol=par$x, nrow=par$y), 
        neutMat = matrix(out[[i]]$NeutralTrait,ncol=par$x, nrow=par$y),
        phylogeny = if (par$convertToBinaryTree && !is.double(phyloi)) ape::multi2di(phyloi) else phyloi, 
        phyloTXT = out[[i]]$Phylogeny,
        summaries = NA)
    }
    cat("done! \n")
    
    par$runtime <- runtime
    
    output$Model <- par
    
    class(output) <- "PhyloSim"
    
    if(par$calculateSummaries) {
      output$Output[[1]]$summaries = calculateSummaryStatistics(output)
    }
    
    return(output)
    
  }else stop("unknown model type: ", par$type, "; for a neutral model in R, use NeutralMod()")
}


# neutral = T
# x = 50
# y = 50
# dispersal = 0.5
# runs = c(100,200)
# specRate = 1.0
# density = F
# environment = F
# seed = 1
# fitnessBaseMortalityRatio = 5
# fitnessActsOn = "mortality" 
# densityCut = 1
# saveTimes = "last"
# 
# 
# if (dispersal == "global" | dispersal == 0){
#   dispersal = 1
#   dispersalCut = 1
# }else if (is.numeric(dispersal)){
#   if (dispersal < 0.5) stop("dispersal parameter too small") 
#   dispersalCut = 2*dispersal
#   dispersal = 3
# }else stop("wrong dispersal parameter") 
# 
# if(fitnessActsOn == "mortality"){
#   reproductiveFitness = F
#   mortalityFitness = T    
# }else if(fitnessActsOn == "reproduction"){
#   reproductiveFitness = T
#   mortalityFitness = F    
# }else if (fitnessActsOn == "both"){
#   reproductiveFitness = T
#   mortalityFitness = T    
# }else stop("wrong fitness parameters")
# 
# if(density == FALSE & environment == FALSE)
# {
#   neutral = TRUE
# }else{
#   neutral = F
# }
# if (length(runs)>1){
#   if (any(runs[-length(runs)] > runs[-1])) stop( "wrong argument to runs")
# }
# 
# 
# out <- callModel( x,  y,  dispersal,  runs,  specRate, density, environment, neutral, mortalityFitness, fitnessBaseMortalityRatio, reproductiveFitness, dispersalCut, densityCut, seed)  


