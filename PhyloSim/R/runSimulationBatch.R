#' @title Batch runner
#' @description A model of species community assembly under different assembly mechanisms, using parallel computing to make use of multi core CPUs and clusters in order to reduce computation time. The function is an extension to \link{runSimulation} in order to accelerate the simulation of multiple scenarios.
#' @param pars A list of parameter sets as created by \code{\link{createCompletePar}} 
#' @param parallel FALSE (default) for sequential runs, or an integer giving the number of cores used for parallel computing. If parallel = TRUE or "auto", n-1 cores are used.
#' @param backup Logical, determining whether the results of the individual scenario runs should be saved with saveRDS() as \emph{scenario}.rds in the working directory (advised if the simulation takes a long time, or if the individual scenarios vary greatly in runtime. Default is FALSE)
#' @param strip NULL or "summaries", determining whether the whole simulation or only the summary statistics should be returned. Summary statistics are only available if \code{calculateSummaries = TRUE} in \code{\link{createCompletePar}}.
#' @param onError "stop" (default) to stop at the first simulation that fails, or "NA" to return NA for each failed simulation, with the error message as attribute "error" and a warning that lists the failed simulations.
#' @return An object of class \code{\link{PhyloSimList}}, or, if strip = "summaries", a plain list of the summary statistics of each scenario. With onError = "NA", failed simulations are NA, so that methods for "PhyloSimList" may fail.
#' @details This function uses the \pkg{foreach} and \pkg{doParallel} packages to compute the model scenarios in parallel on several cores. \cr\cr The phylogeny is passed to R in the newick format and parsed to an object of class "phylo" with the function \code{\link[ape]{read.tree}} from the \pkg{ape} package.
#' @example /inst/examples/runSimulationBatch-help.R
#' @importFrom foreach foreach %dopar% %do%
#' @export
#' 
runSimulationBatch <- function(pars, parallel = FALSE, backup = FALSE, strip = NULL, onError = c("stop", "NA")){
  #start timing
  ptm <- proc.time() 
  onError <- match.arg(onError)
  errorHandling <- if (onError == "NA") "pass" else "stop" # fixes #158
  if (!is.null(strip) && !identical(strip, "summaries")) stop("EcoPhyloSim::runSimulationBatch unrecognized argument to strip")
  
  if (parallel != FALSE){
    cat("running", length(pars), "batch simulations with parallelization\n")
    
    if (parallel == TRUE | parallel == "auto") cores <- parallel::detectCores() - 1
    if (is.numeric(parallel)) cores <- parallel
    cl <- startCluster(cores)
    
    out <- foreach(i=1:length(pars), .packages = c("EcoPhyloSim"), .errorhandling = errorHandling) %dopar%{
      
      OUT <- runSimulation(pars[[i]])
      
      if(backup == TRUE){
        name <- paste0(pars[[i]]$scenario, ".rds")
        saveRDS(OUT, file = name)
      }
      
      if(!is.null(strip) && !is.null(OUT)) {
        if(strip == "summaries") {
          OUT <- OUT$Output[[1]]$summaries
          class(OUT) <- "list"
        }
      }
      OUT
    }
    parallel::stopCluster(cl)
  } else {
    cat("running", length(pars), "batch simulations without parallelization\n")
    out <- foreach(i=1:length(pars), .packages = c("EcoPhyloSim"), .errorhandling = errorHandling) %do%{
      
      cat("running parameter", i, "\n")
      
      OUT <- runSimulation(pars[[i]])
      
      if(backup == TRUE){

        name <- paste0(pars[[i]]$scenario, ".rds")
        saveRDS(OUT, file = name)
      }
      if(!is.null(strip)) {
        if(strip == "summaries") {
          OUT <- OUT$Output[[1]]$summaries
          class(OUT) <- "list"
        }
      }
      OUT
    }
  }
  
  # Fix the naming assignment
  for (i in 1:length(pars)) names(out)[i] <- pars[[i]]$scenario
  
  failed <- which(vapply(out, inherits, logical(1), what = "error"))
  if (length(failed) > 0) {
    warning(length(failed), " of ", length(out), " simulations failed and are NA: ",
            paste0(failed, " (", vapply(out[failed], conditionMessage, ""), ")", collapse = "; "), call. = FALSE)
    out[failed] <- lapply(out[failed], function(e) structure(NA, error = conditionMessage(e)))
  }
  
  if(is.null(strip)) {
    class(out) <- "PhyloSimList"
  } else {
    class(out) <- "list"
  }
  
  # Fix timing output
  time <- proc.time() - ptm
  cat("Finished after", floor(time[3]/60), "minute(s) and", round(time[3]%%60, 2), "second(s).\n")
  
  return(out)
}
