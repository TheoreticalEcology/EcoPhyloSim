#' @title Get mortality of individuals
#' @param out object of type Phylosim or PhylosimList 
#' @description sdfg
#' @details Isdfg
#' @export 
getMortality <- function(out) {
  UseMethod("getMortality")
}


getMortality.Phylosim <- function(out){ # core function
  
  # computes mortality only if generations are consecutive
  if (!(1 %in% diff(out$Model$runs))) {
    stop("No consecutive generations. Can't calculate mortality")
  } 
  
  # create empty mortality matrix for all generations
  for (i in 1:length(out$Output)) {
    out$Output[[i]]$mortMat <- matrix(data = NA, nrow = out$Model$x, ncol = out$Model$y)
  }
  
  # calculates mortalities based on the unique traits per individuum
  for (i in 1:(length(out$Output) -1)) {
    out$Output[[i+1]]$mortMat <-
      ifelse((out$Output[[i]]$traitMat - out$Output[[i+1]]$traitMat) == 0, FALSE, TRUE)
  }
  
  return(out)
}


getMortality.PhylosimBatch <- function(out){ # wrapper function for batch
  for (runs in seq_along(out)) {
    out[[runs]] <- getMortality(out[[runs]])
  }
  return(out)
}


