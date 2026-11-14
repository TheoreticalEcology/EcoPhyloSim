#' @title Null Model Batch Extension
#' @description Generates null models for a list of PhyloSim objects as created by \code{\link{runSimulationBatch}} and compares the observed results against them.
#' @param simu A PhyloSim object or a list of PhyloSim objects
#' @param plotlengths Integer or vector of integers, edge length(s) of the sample plots used in the null models.
#' @param plots Integer, number of sample plots
#' @param replicates Integer, number of replicates
#' @param type String, determines which null model should be created. Possible inputs are "PhylMeta" (default), "PhylSample", "PhylPool", "SamplePool" (see details).
#' @param which.result Integer, determines which result should be used. This argument is only useful if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used. You can choose one or several results (e.g. which.result = c(1, 2)), or use all results with which.result = "all".
#' @param reduce Logical, currently not used by this function (reserved for a future implementation; has no effect on the output).
#' @param fun String, currently not used by this function; MPD is always used.
#' @details In the type argument "PhylMeta" is equivalent to \code{\link{nullModel}} with abundance = FALSE, "PhylSample" is equivalent to the same function with abundance = TRUE, whereas the two other options use \code{\link[picante]{ses.mpd}}. "PhylPool" uses the argument null.model = "phylogeny.pool", "SamplePool" is setting null.model to "sample.pool". \cr\cr For more information see \code{\link{nullModel}}. \cr\cr If you choose which.result = NULL (default), the length of the output corresponds to the length of the 'runs' argument in your parameter definition. Otherwise the list has length max(which.result); elements of results not requested are NULL. \cr\cr The structure of the output is organized as follows: with `output[[t]]` you will get the results for result (time step) 't'. With `output[[t]][[n]]` you will get the results for scenario 'n' at 't'. `output[[t]][[n]]` is a list that contains one vector with p values for each element in the plotlengths argument.
#' @return A list with pValues for each plot in the observed metacommunity. If you have calculated multiple scenarios, the results can be visualized by \code{\link{plotPhylogeneticDispersion}}.
#' @example /inst/examples/phylogeneticDispersion-help.R
#'
#' @export

# Version that accepts multiple types
calculatePhylogeneticDispersion <- function(simu, plotlengths = 10,  plots = 200, replicates = 500, type = "PhylMeta",  which.result = NULL, reduce= FALSE, fun="mpd"){
  
  if ("PhyloSim" %in% class(simu) ==T) simu = list(simu)
  
  ## TODO implement function for only 1 scenario
  
  nScenarios = length(simu) 
  nRuns = length(simu[[1]]$Output)
  parIndex = simu[[1]] 
  times <- which.result
  if(is.null(times)){
    times = nRuns
  }
  
  if(is.null(times) == FALSE){
    if(identical(times, "all")){
      times = c(1:nRuns)     
    }
  }
  # type / scenario / plotlength / repetition
  
  pValues = list()
  
  
  
  for(i in times){
    pValues[[i]] = list()
    timeindex = i
    
    for(k in 1:nScenarios){
      pValues[[i]][[k]] = list()
      
      
      phyl <- simu[[k]]$Output[[timeindex]]$phylogeny
      extantPhylo <- ape::drop.fossil(phy = phyl)
      extantPhyloDis <- cophenetic(extantPhylo)
      
      pval = list()
      for (j in 1:length(plotlengths)){
        
        plotsize = plotlengths[j]
        
        
        
        if (type == "PhylMeta"){
          pval[[j]] <- nullModel(simu[[k]], localPlotSize = plotsize, numberOfPlots = plots,
                                 which.result = timeindex, repetitions = replicates,
                                 fun = "mpd")
        
          }else if (type == "PhylSample"){
          pval[[j]] <- nullModel(simu[[k]], localPlotSize = plotsize, abundance = TRUE,
                                 numberOfPlots = plots, which.result = timeindex, 
                                 repetitions = replicates)    
        
          }else if (type == "PhylPool"){
          comMat <- localPlots(simu[[k]],size = plotsize, n = plots, 
                               which.result = timeindex, community = T)$communityTable
          
          pval[[j]] <- picante::ses.mpd(samp = comMat, dis = extantPhyloDis, 
                                        null.model = "phylogeny.pool", abundance.weighted = TRUE,
                                        runs = replicates)$mpd.obs.p      
        
          }else if (type == "SamplePool"){
          comMat <- localPlots(simu[[k]], size = plotsize, n = plots, 
                               which.result = timeindex, community = T)$communityTable
         
          pval[[j]] <- picante::ses.mpd(samp = comMat, dis = extantPhyloDis, 
                                        null.model = "sample.pool", abundance.weighted = TRUE, 
                                        runs = replicates)$mpd.obs.p
        }
        
      }
      names(pval)<-plotlengths
      
      pValues[[i]][[k]] <- pval
   
    }
    
  }
  
  return(pValues)
}




