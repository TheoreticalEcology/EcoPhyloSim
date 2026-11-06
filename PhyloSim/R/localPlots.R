#' @title Create Subplots
#' @description Generates subplots from a metacommunity, simulating field conditions
#' @param simu Simulation output of the class "Phylosim", usually consisting out of several lists. Needs to contain at least a species matrix ($specMat)
#' @param which.result Integer, determines which result should be used. This argument is only usefull if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default, the last result is used.
#' @param size A single value (or a vector) determining the edge length(s) of the generated subplots
#' @param n The number of subplots to be generated
#' @param community Logical, determining whether to generate a community table or not. default is FALSE
#' @param plot Logical, determining whether the results should be shown graphically. Default is FALSE.
#' @param nested Logical, determining whether the subplots should be nested (see Details).
#' @return A list of subplots and (if coomunity = T) a community table with plots in rows and species in collumns
#' @details If nested == TRUE the subplots are poduced in a sequential order. That means if your 'size' argument has the length 5 and your 'n' argument has the size 2 you will get ten subplots. The first five will be one group of nested plots and the last five the second group of nested plots.
#' @example inst/examples/localPlots-help.R
#' @export
localPlots <- function(simu, which.result=NULL,size, n, community = FALSE, plot = FALSE, nested = FALSE){
  
  checkmate::assert_count(n, positive = TRUE)
  if(community == T && length(size) != 1) stop("Community table only possible for same plot sizes")
  
  if (is.null(which.result)) which.result = length(simu$Output) 
  simu <- simu$Output[[which.result]]
  
  
  matrix <- simu$specMat
  env <- simu$envMat
  comp <- simu$compMat
  
  if(plot == T) plotmat <-matrix(0, ncol(matrix), nrow(matrix))
  
  
  subPlots <- list()
  envPlots <- list()
  compPlots <- list()
  positions <- list()
  count <-0

  for(i in 1:n)
  {           
    if(nested == TRUE){
      x = sample(x=1:length(matrix[1,]),size=1)
      y = sample(x=1:length(matrix[,1]),size=1)
    }
    
    for(k in 1:length(size)){
      count <- count + 1
      edge <- size[k] 
      
      if(nested == FALSE){
        x = sample(x=1:length(matrix[1,]),size=1)
        y = sample(x=1:length(matrix[,1]),size=1)      
      }
      
      # Check for boundary issues and implemant warped bounaries if necessary
      if(x+(edge-1) > length(matrix[1,]) && y+(edge-1) <= length(matrix[,1])){
        rsel = c(seq(x,length(matrix[1,])), seq(1,(edge-length(seq(x,length(matrix[1,]))))))
        csel = seq(y,y+(edge-1))
      }
      else if(x+(edge-1) <= length(matrix[1,]) && y+(edge-1) > length(matrix[,1])){
        rsel = seq(x,x+(edge-1)) 
        csel = c(seq(y,length(matrix[,1])),seq(1,(edge-length(seq(y,length(matrix[,1]))))))
      }
      else if(x+(edge-1) > length(matrix[1,]) && y+(edge-1) > length(matrix[,1])){
        rsel = c(seq(x,length(matrix[1,])), seq(1,(edge-length(seq(x,length(matrix[1,]))))))
        csel = c(seq(y,length(matrix[,1])), seq(1,(edge-length(seq(y,length(matrix[,1]))))))
      }
      else{
        rsel = seq(x,x+(edge-1))
        csel = seq(y,y+(edge-1))
      }
      
      
      subPlots[[count]]  <- matrix[rsel, csel]
      envPlots[[count]] <- env[rsel, csel]
      compPlots[[count]] <- comp[rsel, csel]
      
      # TODO: Is this the correct way to calculate the center of a rectangle in a matrix?
      position <- c(rsel[1] + size[k] / 2 - 0.5, csel[1] + size[k] / 2 - 0.5)
      positions[[count]] <- c(wrap_coordinates(position[1], limit=length(matrix[1,])),
                              wrap_coordinates(position[2], limit=length(matrix[,1])))
      
      
      if(plot == T){
        plotmat[rsel, csel] <- 1
        
      }
    }
  }
  
  
  if(community == T)
  {
    plotID <- rep(seq_along(subPlots), times = lengths(subPlots))
    communityTable <- as.data.frame.matrix(table(
      factor(plotID, levels = seq_along(subPlots), labels = paste0("plot", seq_along(subPlots))),
      factor(unlist(subPlots))))
    names(communityTable) <- paste0("s", names(communityTable))
  }
  
  if(plot == T){
    oldpar <- par()$mar
    par(mar=c(1,1,1,1))
    image(matrix,  yaxt='n', xaxt='n', asp=1, bty="n")
    image(plotmat, col=c(rgb(0,0,0,0),rgb(0,0,1,0.7)), add=T)
    par(mar=oldpar)
  } 
  
  if(community == T){
    return(list(subPlots=subPlots, communityTable=communityTable, envPlots=envPlots, compPlots=compPlots, positions=positions))
  } else return(list(subPlots=subPlots, envPlots=envPlots, compPlots=compPlots, positions=positions))
  
}

wrap_coordinates <- function(x, limit) if (x > limit) x - limit else if (x <= 0) limit + x else x
