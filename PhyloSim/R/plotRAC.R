#' @title Rank Abundance Curve
#' @description Plots the Rank Abundance Curve for a given community. 
#' @param simu Simulation output of the class "Phylosim", usually consisting out of several lists. Needs at least the spatial distribution of the species stored in a matrix ($specMat)
#' @param which.result Integer, determines which result should be used. This argument is only usefull if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used. If you choose "all" all results are shown in one plot (see 'Details).
#' @param plot determining whether to plot the RAC as "line"(default) or "bar".
#' @param title String, determining the title of the plot.
#' @details Each species is given a rank according to their abundance (highest = rank 1). Then the species' abundance is plotted in dependency of their rank. It can be used as an indicator for the amount of equally abundant species a community can support. \cr\cr If which.result = "all" all intermediate results are shown in one plot. The colors of the lines are plotted as a gradient from blue (first results) to red (end result).
#' @return A dataframe containing the ranked abundances, sorted by ascending rank. If which.result = "all" only the plot will be returned.
#' @example inst/examples/plotRAC-help.R
#' @export

rac <- function(simu,which.result=NULL, plot="line", title="RAC"){
  
  if(is.null(which.result)) which.result = length(simu$Output) 
  
  if(is.null(which.result) == FALSE){
    if(identical(which.result, "all")){
      if(plot=="bar")stop("Argument 'bar' not possible for which.result='all'")
      simulations <- c(1:length(simu$Output))
    }else{
      simulations<-which.result
    } 
  }
  
  colfunc <- colorRampPalette(c("blue", "red"))
  cols<-colfunc(length(simulations))
  
  RAC<-list()
  max_abundance<-0
  max_rank<-0
  
  for(i in simulations){
    simu_t <- simu$Output[[i]]
    matrix <- simu_t$specMat
    
    Abundances <- as.data.frame(table(matrix))
    sel <- order(Abundances$Freq, decreasing=T)
    RAC[[i]] <- data.frame(Rank = seq(1,length(Abundances$Freq),1), 
                  Abundance = Abundances$Freq[sel], Species = Abundances$matrix[sel])

    if(max(RAC[[i]]$Abundance)>max_abundance) max_abundance<-max(RAC[[i]]$Abundance)
    if(max(RAC[[i]]$Rank)>max_rank) max_rank<-max(RAC[[i]]$Rank)
    
  }
  
  if(plot=="bar"){
    barplot(RAC[[i]]$Abundance, log="y",ylab="Log Abundance", xlab="Rank", 
             main=title, names.arg = RAC[[i]]$Rank)
  }
  if(plot=="line"){
    if(length(simulations) == 1){
      plot(RAC[[i]]$Rank, RAC[[i]]$Abundance, type="l",log="y",ylab="Log Abundance",
            xlab="Rank", main=title, lwd=2)
    }else{
      for(i in simulations){  
        if(i ==1){
          plot(RAC[[i]]$Rank, RAC[[i]]$Abundance, type="l", log="y", ylab="Log Abundance", 
                xlab="Rank", main=title, col=cols[i],xlim=c(0, max_rank), ylim=c(1, max_abundance)) 
        }else{
          lines(RAC[[i]]$Rank, RAC[[i]]$Abundance, type="l", main=title, col=cols[i])
        }
      }
    }
  }
  if(length(simulations)==1) return(RAC[[simulations]]) 
}




