#' @title  Species Community Simulation
#' @description A model of species community assembly under different assembly mechanisms.  
#' @param par a list of model parameters created with \link{createCompletePar}
#' @return An object of class \code{PhyloSim}. This object contains the species matrix, the trait matrix, the environmental matrix, the competition matrix and the neutral matrix, as well as the phylogeny and the parameter settings of the simulation. 
#' @details If the runs argument of your parameter set contains more than one element, each interim step is saved in the PhyloSim object. \cr\cr For larger simulations consider \link{runSimulationBatch} to make use of parallel computing. \cr\cr It is possible that more than one new species arises per generation. This leads to a multifurcated phylogeny, yet many packages such as "ape" can only work with bifurcating trees. If `convertToBinaryTree = TRUE` (the default), the phylogeny is converted to a bifurcating one using `multi2di(random = FALSE)` from the "ape" package, so that the same `seed` gives the same phylogeny. \cr\cr Branch lengths of the phylogeny are in generations. Speciation and extinction events of generation g are placed at time g, and the census of generation g (an entry in \code{runs}) at time g + 0.5, after all events of that generation. Tips of extant species therefore end half a generation after those of species that went extinct in the last generation, so that \code{ape::drop.fossil()} removes the latter. Internal nodes are labeled with the species from which the descendants emerged (nodes added by \code{convertToBinaryTree} have no label); the tree has no root edge. With \code{prunePhylogeny = TRUE} (the default), the phylogeny contains only the extant species; with \code{prunePhylogeny = FALSE}, it also contains all extinct species as tips. \cr\cr If a result contains only one species, no phylogeny can be built, and its `phylogeny` element is set to 0 with a warning. This is common in short runs, because all individuals start as one species and new species arise at a rate of `specRate` per generation. \cr\cr The simulation uses R's random number generator, set with \code{set.seed(par$seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")}; the random number state of the R session is restored afterwards. Results for a given seed therefore do not depend on the compiler's random number library.
#' @import Rcpp
#' @useDynLib EcoPhyloSim, .registration = TRUE
#' @export
#' @example /inst/examples/runSimulation-help.R


runSimulation <- function(par)
{
  par <- checkPar(par)
  localSeed(par$seed) # also seeds the summary statistics (#135)
  
  # GENERAL SETTINGS
  
  if (par$type == "base"){
    
    if (par$dispersal == "global" | par$dispersal == 0){
      dispersal = 1
      dispersalCut = 1
    }else if (is.numeric(par$dispersal)){
      dispersalCut = 2*par$dispersal
      dispersal = 3
    }else stop("wrong dispersal parameter") 
    
    actsOn <- fitnessStages(par$fitnessActsOn)
    mortalityFitness <- "mortality" %in% actsOn
    reproductiveFitness <- "reproduction" %in% actsOn

    kernels <- c("absolute", "gaussian", "species")
    # user-defined density response, stored by checkPar() as values over the similarity S in [0, 1]
    densityTable <- if (is.numeric(par$negativeDDResponse)) par$negativeDDResponse else numeric(0)
    
    # negativeDD and environment stay numeric, so that the stored $Model passes checkPar() when rerun (#118)
    par$compStrength <- if (par$negativeDD == 0) 1 else par$negativeDD
    par$envStrength <- if (par$environment == 0) 1 else par$environment
    neutral <- par$negativeDD == 0 && par$positiveDD == 0 && par$environment == 0
    
    ptm <- proc.time()
    
    out <- callModel( x = par$x, 
                      y = par$y, 
                      dispersal = dispersal, 
                      runs = round(par$runs), 
                      specRate = par$specRate, 
                      dens = par$negativeDD > 0, 
                      env = par$environment > 0, 
                      neutral = neutral, 
                      mort = mortalityFitness, 
                      mortStrength = par$fitnessBaseMortalityRatio, 
                      repro = reproductiveFitness, 
                      dispersalCutoff = dispersalCut, 
                      densityCutoff = par$negativeDDCut, 
                      envStrength = par$envStrength, 
                      compStrength = par$compStrength, 
                      fission = par$fission, 
                      redQueen = par$redQueen, 
                      redQueenStrength = par$redQueenStrength, 
                      protracted = par$protracted, 
                      airmatR = par[["airmatR"]], 
                      soilmatR = par[["soilmatR"]],
                      prunePhylogeny = par$prunePhylogeny,
                      nicheWidth = par$nicheWidth,
                      recruit = "recruitment" %in% actsOn,
                      densityKernel = match(par$negativeDDKernel, kernels),
                      densityWidth = par$negativeDDWidth,
                      densityTable = densityTable,
                      positiveStrength = par$positiveDD,
                      positiveCutoff = par$positiveDDCut,
                      positiveKernel = match(par$positiveDDKernel, kernels),
                      positiveWidth = par$positiveDDWidth,
                      speciationTraitJump = rep_len(par$speciationTraitJump, 3))
    
    runtime <- as.numeric((proc.time() - ptm)[3])
    
    print (paste("Core simulation finished after",floor(((proc.time() - ptm)[3])/60), "minute(s) and", ((proc.time() - ptm)[3])%%60, "second(s). Converting data"))
    
    # The following code reshuffles the output to the structure 
    # of a phylosim class
    output<-list()
    for (i in 1:length(par$runs)){
      
      specMati = matrix(out[[i]]$Species,ncol=par$x, nrow=par$y)
      if(length(unique(c(specMati)))==1){
        phyloi <- 0
        warning("Only one species in the community at generation ", par$runs[i], "; phylogeny set to 0. Increase runs or specRate.")
      }else{
        phyloi <- ape::read.tree(text= out[[i]]$Phylogeny)
      }
       
      output$Output[[i]] = list(
        specMat = matrix(out[[i]]$Species,ncol=par$x, nrow=par$y), 
        traitMat= matrix(out[[i]]$EnvTrait,ncol=par$x, nrow=par$y), 
        envMat = matrix(out[[i]]$Environment,ncol=par$x, nrow=par$y), 
        compMat = matrix(out[[i]]$CompetitionTrait,ncol=par$x, nrow=par$y), 
        neutMat = matrix(out[[i]]$NeutralTrait,ncol=par$x, nrow=par$y),
        phylogeny = if (par$convertToBinaryTree && !is.double(phyloi)) ape::multi2di(phyloi, random = FALSE) else phyloi, 
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
