#' @title Run PhyloSim Analysis
#' @description This function is used to perform many (>1000) PhyloSim simulations in parallel with only one function call. The user can specify parameters that are fixed across all simulations and parameters that are variable. The results are saved to a specified folder. 
#' @param variableParams A named list, where the names are parameter names and the values are priors (functions). A call to the function must return a single value. See details.
#' @param fixedParams A named list, where the names are parameter names and the values are parameter values for the PhyloSim model.
#' @param nRuns Integer, specifying the number of simulations to run
#' @param outDir String, specifying the directory to which the output will be written
#' @param summariesOnly Logical. If TRUE, only the summary statistics are saved, to the files summary_FROM-TO instead of out_FROM-TO, so that \code{\link{chunkwiseCalculateSummaries}} is not needed. Requires calculateSummaries = TRUE in fixedParams.
#' @param nCores Integer or "auto", number of parallel cores to use for simulations. "auto" uses all detected cores minus one.
#' @param chunkSize Integer or "auto" (default), specifying the number of simulations per save file. Values above 10,000 give a warning, as they might deplete the RAM.
#' @details variableParams must be a list of functions. The idea behind this is to draw from a prior, so the functions should return a single value drawn from the corresponding parameter prior. However, the implementation allows the user to specify stateful functions instead of priors to perform analyses with variable parameters that are not drawn from a prior, but follow a sequence e.g. run simulations with varying specRate.
#' @examples 
#' \dontrun{
#' variableParams <- list(environment=function()return(runif(1)),
#'                        specRate=function()return(runif(1, 1, 10)))
#' fixedParams <- list(x=30, y=30,
#'                     runs=1000,
#'                     densityCut=1,
#'                     calculateSummaries=TRUE)
#' 
#' runPhylosimAnalysis(variableParams,
#'                     fixedParams,
#'                     1000,
#'                     "outDirName")
#' }
#' @return None. Writes paramOverview.txt and the files param_FROM-TO and out_FROM-TO (summary_FROM-TO if summariesOnly = TRUE) to outDir.
#' @author Tankred Ott
#' @export
# TODO: Rename this
runPhylosimAnalysis <- function(variableParams, fixedParams, nRuns, outDir, summariesOnly=FALSE, nCores="auto", chunkSize="auto") {
  # Get available parameters
  arguments <- formals(EcoPhyloSim::createCompletePar)
  
  # Validate input
  validateParams(names(variableParams), names(arguments))
  validateParams(names(fixedParams), names(arguments))
  validatePriors(variableParams)
  cat("ATTENTION: createCompletePar() stops with an error if a fixed parameter or a value drawn from a prior is invalid, which stops the analysis. Please make sure to read the help of EcoPhyloSim::createCompletePar to choose valid parameters.\n")
  
  # Find parameters that stay at default values
  defaultParams <- getDefaultParams(arguments, c(fixedParams, variableParams))
  
  # if output directory doesn't exist create it
  if(!dir.exists(outDir)){
    dir.create(outDir)
  }
  
  # write an first overview over the parameters to file
  writeParamOverview(defaultParams, fixedParams, variableParams, joinPath(c(outDir, "paramOverview.txt")))
  # # read it later in like this:
  # inp <- readLines("argsOut.txt")
  # strsplit(inp, split = " ")
  
  # get the number of parallel cores
  nCores <- getValidatedNCores(nCores)
  
  # determine chunk size
  chunkSize <- getValidatedChunkSize(chunkSize, nRuns)
  
  # Create a list of parameter that are fixed over all simulations
  fixedParamsList <- c(defaultParams, fixedParams)
  
  # Runs the simulation in chunks and saves the results to outDir
  runSimulationChunkwise(fixedParamsList, variableParams, nRuns, chunkSize, outDir, summariesOnly, nCores)  
}

#' @title Chunkwise summary statistics calculation
#' @description Calculates summary statistics for all "out_" files created with runPhylosimAnalysis and saves the results to "summary_" files in the same folder.
#' @param dir Directory from which the "out_" files should be read and to which the "summary_" files will be saved.
#' @param nCores "auto" or number of parallel cores.
#' @return None. Writes summary_FROM-TO files to dir.
#' @author Tankred Ott
#' @export
chunkwiseCalculateSummaries <- function (dir=".", nCores="auto") {
  simFiles <- getFileNames(path=dir)
  if (length(simFiles) == 0) stop(paste("No PhyloSim output files found in", if (dir == ".") getwd() else dir))
  
  # this is inefficient for large files, as the whole file is loaded into memory just to get the chunkSize
  # better might be to run the first summary calculation without parallelization and start the parallel runs with the 
  # second file
  chunkSize <- length(loadLocal(joinPath(c(dir, simFiles[1]))))
  
  nCores <- getValidatedNCores(nCores)
  cl <- parallel::makeCluster(nCores)
  parallel::clusterCall(cl, function(lp) .libPaths(lp), .libPaths())
  doParallel::registerDoParallel(cl)
  
  foreach(i=1:length(simFiles), .packages = c("EcoPhyloSim")) %dopar% {
    loadLocal <- function(path) return(local(get(load(path))))
    joinPath <- function(x, sysName = "auto") {
      if (sysName == "auto") sysName <- Sys.info()[["sysname"]]
      if(sysName == "Windows") return(paste(x, collapse = "\\"))
      else if(sysName %in% c("Linux", "Darwin")) return(paste(x, collapse = "/"))
      else stop(paste(c("System", sysName, "is not supported."), collapse = " "))
    }
    
    # Load simulation chunk
    simChunk <- loadLocal(joinPath(c(dir, simFiles[i])))
    
    # Calculate summaries
    tryCalcSummaries <- function(x,...)return(tryCatch(EcoPhyloSim::calculateSummaryStatistics(x,...),
                                                   error=function(x)return(NA)))
    # summaryChunk <- lapply(simChunk, EcoPhyloSim::calculateSummaryStatistics, strict=TRUE)
    summaryChunk <- lapply(simChunk, tryCalcSummaries, strict=TRUE)
    
    # Determine file name
    n <- length(simChunk)
    from <- as.integer((i-1) * chunkSize + 1)
    to <- as.integer(from + n - 1)
    suffix <- paste(c("_", from,"-", to), collapse = "")
    summaryFileName <-  paste(c("summary", suffix), collapse = "")
    # save file
    save(summaryChunk, file = joinPath(c(dir, summaryFileName)))
  }
  parallel::stopCluster(cl)
}

#' @title Clean up parameter and summary chunks
#' @description Reads in summary files created with chunkwiseCalculateSummaries or with runPhylosimAnalysis(summariesOnly = TRUE), and parameter files created with runPhylosimAnalysis and removes crashed simulations (those where summaries == NA). The resulting cleaned-up chunks are saved to the subfolder "cleaned_up".
#' @author Tankred Ott
#' @param dir Directory to which the results of runPhylosimAnalysis (and chunkwiseCalculateSummaries) were saved.
#' @param nCores "auto" or number of parallel cores to use
#' @return None. Writes summary_clean_FROM-TO and param_clean_FROM-TO files to dir/cleaned_up.
#' @export
chunkwiseCleanup <- function (dir=".", nCores="auto") {
  nCores <- getValidatedNCores(nCores)
  cl <- parallel::makeCluster(nCores)
  parallel::clusterCall(cl, function(lp) .libPaths(lp), .libPaths())
  doParallel::registerDoParallel(cl)
  
  sumFiles <- getFileNames(path=dir, pattern="^summary_[0-9]+-[0-9]+$")
  if (length(sumFiles) == 0) stop(paste("No PhyloSim summary files found in", if (dir == ".") getwd() else dir))
  
  parFiles <- getFileNames(path=dir, pattern="^param_[0-9]+-[0-9]+$")
  if (length(parFiles) == 0) stop(paste("No PhyloSim parameter files found in", if (dir == ".") getwd() else dir))
  
  chunkSize <- length(loadLocal(joinPath(c(dir, sumFiles[1]))))
  
  # Create output directory
  outDir <- joinPath(c(dir, "cleaned_up"))
  if(!dir.exists(outDir)){
    dir.create(outDir)
  }
  
  foreach(i=1:length(sumFiles), .packages = c("EcoPhyloSim")) %dopar% {
    # loadLocal <- function(path) return(local(get(load(path))))
    # joinPath <- function(x, sysName = "auto") {
    #   if (sysName == "auto") sysName <- Sys.info()[["sysname"]]
    #   if(sysName == "Windows") return(paste(x, collapse = "\\"))
    #   else if(sysName %in% c("Linux", "Darwin")) return(paste(x, collapse = "/"))
    #   else stop(paste(c("System", sysName, "is not supported."), collapse = " "))
    # }
    
    # Load simulation chunk and parameter chunk
    sumChunk <- loadLocal(joinPath(c(dir, sumFiles[i])))
    parChunk <- loadLocal(joinPath(c(dir, parFiles[i])))
    # Find indices where summaries are NA and remove them
    naIndices <- which(is.na(sumChunk))
    naIndices <- c(naIndices, which(sapply(sumChunk, function (x) sum(is.na(x)) > 0)))
    
    sumChunk[naIndices] <- NULL
    parChunk[naIndices] <- NULL
    
    # replace NULL with "NULL". Necessary for data.table::rbindlist in later analysis
    for(j in 1:length(parChunk)) {
      for(k in 1:length(parChunk[[j]])) {
        if (is.null(parChunk[[j]][[k]])) parChunk[[j]][[k]] <- "NULL"
      }
    }
    
    # Determine file name
    n <- length(sumChunk)
    from <- as.integer((i-1) * chunkSize + 1)
    to <- as.integer(from + n - 1)
    suffix <- paste(c("_", from,"-", to), collapse = "")
    sumFileName <-  paste(c("summary_clean", suffix), collapse = "")
    parFileName <-  paste(c("param_clean", suffix), collapse = "")
    # save file
    save(sumChunk, file = joinPath(c(outDir, sumFileName)))
    save(parChunk, file = joinPath(c(outDir, parFileName)))
  }
  parallel::stopCluster(cl)
}

#' @title combine chunks
#' @description Combines files with a common prefix into a combined file and saves it as combined_PREFIX.
#' @param dir Directory. Location of the files
#' @param prefix Common file prefix. Possible to use RegEx like syntax here.
#' @return None. Writes the combined list to the file combined_PREFIX in dir.
#' @author Tankred Ott
#' @export
combineChunks <- function(dir=".", prefix) {
  files <- getFileNames(path=dir, pattern=paste("^", prefix,".*", sep = ""))
  if (length(files) == 0) stop(paste("No files with prefix", prefix, "found in", if (dir == ".") getwd() else dir))
  out <- list()
  for (i in 1:length(files)) {
    chunk <- loadLocal(joinPath(c(dir, files[i])))
    out <- c(out, chunk)
  }
  fileName <- paste("combined_", prefix, sep = "")
  save(out, file = joinPath(c(dir, fileName)))
}

#' @title Extract params
#' @description Extracts parameters from a list of PhyloSim parameter objects.
#' @param parList List of phylosim parameter objects
#' @param params Vector of parameter names to extract
#' @param outFormat Type of output that should be returned. Only "df" (data frame) is implemented; "list" gives an error.
#' @return A data frame with one column per parameter in params.
#' @author Tankred Ott
#' @export
extractParams <- function(parList, params, outFormat="df") {
  # check if params in parList
  parNames <- names(parList[[1]])
  invalidParams <- which(!is.element(params, parNames))
  if (length(invalidParams != 0)) stop(paste(paste(params[invalidParams], collapse = ", "), "is/are invalid parameter/s"))

  out <- NA
  if (outFormat == "df") {
    # out <- as.data.frame(do.call(what = rbind, args = parList))
    checkSuggested("data.table")
    out <- as.data.frame(data.table::rbindlist(parList))
    out <- out[, params]
  } else if (outFormat == "list") {
    stop("no implemented")
  } else stop("Wrong argument for outFormat")
  return(out)
}

#' @author Tankred Ott
# This function devides the runs in chunks of chunkSize, generates a phylosim parameter object for each run by
# combining the fixed parameters and draws from the priors of the variable parameter, runs the phylosim model in
# parallel mode and saves the results and the parameters in outDir
runSimulationChunkwise <- function (fixedParamsList, variableParams, nRuns, chunkSize, outDir, summariesOnly, nCores) {
  if (summariesOnly && !isTRUE(fixedParamsList$calculateSummaries)) stop("summariesOnly is set to TRUE, but calculateSummaries is not TRUE in fixedParams")
  nChunks <- ceiling(nRuns / chunkSize)
  remainder <- nRuns %% chunkSize
  for (i in 1:nChunks) {
    cat(paste(c("Running chunk",i ,"of", nChunks, "\n"), collapse = " "))
    t0 <- proc.time()
    
    n <- chunkSize # number of simulations per chunk
    # if nRuns not divisible by chunkSize, we need to add an additional
    # smaller chunk at the end
    if (i == nChunks && (remainder != 0)) {
      n <- remainder
    }
    
    # create phyloSim parameter object
    params <- rep(list(NA), n)
    for (j in 1:n) {
      # create list of variable parameters drawn from prior
      varParamsList <- lapply(variableParams, FUN = do.call, list())
      params[[j]] <- do.call(EcoPhyloSim::createCompletePar, args = c(fixedParamsList, varParamsList))
    }
    
    # generate "unique" suffix for data files
    from <- (i-1) * chunkSize + 1
    to <- from + n - 1
    from <- format(from, scientific = FALSE)
    to <- format(to, scientific = FALSE)
    suffix <- paste(c("_", from,"-", to), collapse = "")
    paramFileName <- paste(c("param", suffix), collapse = "")
    outFileName <- paste(c(if (summariesOnly) "summary" else "out", suffix), collapse = "")
    
    save(params, file = joinPath(c(outDir, paramFileName)))
    
    # run simulations in parallel
    out <- tryCatch(runSimulationBatch(params, nCores, strip = if (summariesOnly) "summaries"), error=function(e)return(NA))
    
    save(out, file = joinPath(c(outDir, outFileName)))
    cat(paste(c("Finished chunk", i ,"of", paste0(nChunks, "."), "time elapsed:", proc.time()[3]-t0[3], "\n"), collapse = " "))
  }
  cat(paste(c("Output and parameters are saved to directory", c("'", outDir, "'")), collapse = " "))
}

#' @title Prepare for Analysis
#' @description Reads the files combined_summary_clean and combined_param_clean in the subfolder cleaned_up (created by combineChunks with prefix "summary_clean" and "param_clean") and creates data frames for later analysis.
#' @param dir root directory of the analysis (as specified in runPhylosimAnalysis)
#' @param variableParams A named list (or vector) whose names match the variable parameter names used in the corresponding \code{\link{runPhylosimAnalysis}} call. Only the names are used here to select which parameter columns to extract -- the values themselves do not need to match.
#' @return Create a subfolder "analysis_ready" and write the files "parameters" (variable parameters), "summaries" and "summaries_scaled" containing data.frames with the corresponding data.
#' @export
#' @author Tankred Ott


prepareForAnalysis <- function(dir=".", variableParams) {
  # load Data
  cleanedDir <- joinPath(c(dir,"cleaned_up"))
  summaries <- loadLocal(joinPath(c(cleanedDir, "combined_summary_clean")))
  pars <- loadLocal(joinPath(c(cleanedDir, "combined_param_clean")))
  
  # extract variable parameters' values
  pars_df <- EcoPhyloSim::extractParams(pars, names(variableParams))
  
  # scale summaries
  summaries_df <- as.data.frame(data.table::rbindlist(summaries))
  summaries_sd <- apply(X = summaries_df, MARGIN = 2, FUN=sd)
  summaries_df_scaled <- as.data.frame(scale(summaries_df, scale = summaries_sd, center = F))
  
  analysisDir <- joinPath(c(dir,"analysis_ready"))
  if(!dir.exists(analysisDir)){
    dir.create(analysisDir)
  }
  save(pars_df, file = joinPath(c(analysisDir, "parameters")))
  save(summaries_df, file = joinPath(c(analysisDir, "summaries")))
  save(summaries_df_scaled, file = joinPath(c(analysisDir, "summaries_scaled"))) 
}

#' @author Tankred Ott
# Checks if the input parameters are valid parameters for a function
# with the parameters arguments
validateParams <- function(input, arguments) {
  isE <- is.element(input, arguments)
  if (sum(isE) != length(input)) {
    errInputs <- input[!isE]
    multiple <- length(errInputs) > 1
    errMsg <- paste(c(paste(errInputs, collapse = ", ")),
                    if (multiple) "are" else "is","not", if (!multiple) "a",
                    "valid", if (multiple) "parameters" else "parameter")
    stop(errMsg)
  }
}

#' @author Tankred Ott
# Check if priors are behaving as expected
validatePriors <- function(variableParams) {
  for (i in 1:length(variableParams)) {
    if(!is.function(variableParams[[i]])) stop(paste(c("The prior for parameter",
                                                       names(variableParams)[i],
                                                       "must be a function"),
                                                     collapse = " "))
    # Call to prior should perform a single draw/return a single number
    # This is a little bite ineffecient, as we are only drawing 1 value at a time.
    # However, this should make no difference as the bottleneck is the simulation,
    # and the code becomes more readable and extensible this way
    x <- variableParams[[i]]()
    if(!(is.vector(x) && length(x) == 1)) stop(paste(c("Expected prior",
                                                       names(variableParams)[i], 
                                                       "to return a single value"),
                                                     collapse = " "))
  }
}

#' @author Tankred Ott
# expects a lists of parameter names and the corresponding values
# return a list of parameters and corresponding values that stay at the default value
getDefaultParams <- function (allParams, nonDefaultParams) {
  replacedIndices <- is.element(names(allParams), names(nonDefaultParams))
  return(allParams[!replacedIndices])
}

#' @author Tankred Ott
# Writes the parameters and corresponding values (for fixed and default, but not
# for variable as those are priors) to a file
writeParamOverview <- function(defaultParams, fixedParams, variableParams, f="param_overview.txt") {
  argsOut <- list(defaultParams=names(defaultParams),
                  defaultArgsParams=unlist(defaultParams, use.names = FALSE),
                  fixedParams=names(fixedParams),
                  fixedParamsVals=unlist(fixedParams, use.names = FALSE),
                  variableParams=names(variableParams))
  
  for(i in 1:length(argsOut)) {
    append <- if (i == 1) FALSE else TRUE
    out <- paste(c(names(argsOut)[i], argsOut[[i]]))
    write(file = f, x = out, ncolumns = length(out), append = append)
  }  
}

#' @author Tankred Ott
# checks the chunkSize parameter corresponding to nRuns
# Implements the automatic chunkSize adaption when chunkSize="auto"
getValidatedChunkSize<- function (chunkSize, nRuns, maxChunkSize = 10000) {
  if (is.character(chunkSize)) {
    if (chunkSize == "auto") {
      if (nRuns <= 100) chunkSize <- nRuns
      else if (nRuns <= 10000) chunkSize <- 100
      else chunkSize <- min(floor(nRuns / 100), maxChunkSize)
      cat(paste(c("chunkSize automatically set to", chunkSize, "\n"), collapse = " "))
    } else {
      stop(paste(c(chunkSize, "is not a valid input for chunkSize."), collapse = " "))
    }
  } else if (chunkSize > maxChunkSize) {
    warning(paste(c("chunkSize is greater than",
                    maxChunkSize,
                    ". There might not be enough RAM available."),
                  collapse = " "))
  }
  return(chunkSize)
}

#' @author Tankred Ott
# joins directory paths (windows and Linux)
joinPath <- function(x, sysName = "auto") {
  if (sysName == "auto") sysName <- Sys.info()[["sysname"]]
  if(sysName == "Windows") return(paste(x, collapse = "\\"))
  else if(sysName %in% c("Linux", "Darwin")) return(paste(x, collapse = "/"))
  else stop(paste(c("System", sysName, "is not supported."), collapse = " "))
}

# returns the number of parallel cores
getValidatedNCores <- function (nCores) {
  nMaxCores <- parallel::detectCores()
  if (is.character(nCores)) {
    if (nCores == "auto") {
      if (is.na(nMaxCores)) stop("Could not detect the number of available cores. Please enter the number of parallel cores to use manually.")
      else {
        nCores <- nMaxCores - 1
        cat(paste(c("nCores set to", nCores, "\n"), collapse = " "))
      }
    }
    # nCores is not a string
  } else {
    if (!is.na(nMaxCores)) {
      if (nCores > nMaxCores) stop("nCores may not be larger than the number of cores available")
      else if (nCores == nMaxCores) warning("All available cores are used for simulation. This will affect the performance of your computer.")
    }
  }
  return(nCores)
}

#' @author Tankred Ott
loadLocal <- function(path) return(local(get(load(path))))

# stops with an informative error if a package listed in Suggests is missing
checkSuggested <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required for this function. Please install it.", call. = FALSE)
}

#' @author Tankred Ott
getFileNames <- function(path=".", pattern="^out_[0-9]+-[0-9]+$") {
  files <- list.files(path = path, pattern = pattern)
  checkSuggested("gtools")
  files <- gtools::mixedsort(files)
  return(files)
}

#' @title Rejection sampling for PhyloSim
#' @description Applies rejection sampling using parameter and summary files generated by the PhyloSim analysis pipeline (runPhylosimAnalysis, etc.).
#' @param p Path to a file containing the parameters in form of a data.frame/matrix, or a matrix/data.frame containing the parameters.
#' @param s Path to a file containing the (scaled) summary statistics, or a matrix/data.frame containing the summaries.
#' @param nDraws Number of ABC rejection sampling draws
#' @param eQuantile Epsilon expressed as quantile. E.g. 0.025 are only the nearest 2.5%.
#' @param nCores Number of parallel cores, or "auto".
#' @param which vector of strings, determining which summary statistics reduction should be applied. Each supplied method will be applied separately. Possible options: "none", "saABC", "rf"
#' @param sumFun Function to summarize the accepted simulated parameters, e.g. mean or median
#' @param ... additional parameters passed to summaryRedRF (and to ranger::ranger)
#' @return A list of matrices, one for each method in \code{which}, each containing the true parameters postfixed with "_true" and the estimated parameters postfixed with "_estimated".
#' @author Tankred Ott
#' @export

rejectionSamplePhylosimAnalysis <- function (p, s, nDraws, eQuantile=0.025, nCores="auto", which=c("none", "saABC", "rf"), sumFun=median, ...) {
  pars <- NULL
  sums <- NULL
  
  if (is.character(p)) {
    pars <- loadLocal(p)  
  } else {
    pars <- p
  }
  
  if (is.character(s)) {
    sums <- loadLocal(s)  
  } else {
    sums <- s
  }
  
  pars <- as.matrix(pars)
  sums <- as.matrix(sums)
  
  rsResult <- vector("list", 0)
  
  for (type in which) {
    if (type == "none") {
      rsResult$none <- rejectionSample(pars, sums, nDraws, eQuantile, nCores, sumFun)
    } else if (type == "saABC") {
      predPar <- summaryRedsaABC(p = pars, s = sums, ...)
      rsResult$saABC <- rejectionSample(pars, predPar, nDraws, eQuantile, nCores, sumFun)
    } else if (type == "rf") {
      predPar <- summaryRedRF(p = pars, s = sums, nCores = nCores, ...)
      rsResult$rf <- rejectionSample(predPar$pTrue, predPar$pPred, nDraws, eQuantile, nCores, sumFun)
    } else {
      warning(paste("Could no recognize argument", type,"for 'which'"))
    }
  }
  
  return(rsResult)
}

#' @title Rejection Sampling
#' @description Applies rejection sampling.
#' @param params A data.frame of the parameters (different parameters in columns)
#' @param summaries A data.frame of the summary statistics (different summaries in columns)
#' @param nDraws Number of ABC rejection sampling draws
#' @param eQuantile Epsilon expressed as quantile. E.g. 0.025 are only the nearest 2.5%.
#' @param nCores Number of parallel cores, or "auto".
#' @param sumFun Function to summarize the accepted simulated parameters, e.g. mean or median
#' @return A matrix containing the true parameters postfixed with "_true" and the estimated parameters postfixed with "_estimated".
#' @author Tankred Ott
#' @export
rejectionSample <- function(params, summaries, nDraws = 10, eQuantile = 0.025, nCores="auto", sumFun=median){
  # draw random samples
  drawIndices <- sample(1:nrow(params), nDraws)
  
  # initialize parallel cores
  nCores <- getValidatedNCores(nCores)
  cl <- parallel::makeCluster(nCores)
  parallel::clusterCall(cl, function(lp) .libPaths(lp), .libPaths())
  doParallel::registerDoParallel(cl)
  
  # for each draw do rejection sampling against all remaining samples
  OUT <- foreach(i=1:nDraws, .combine = rbind) %dopar% {
    drawIndex <- drawIndices[i]
    
    obsSummaries <- summaries[drawIndex,]
    obsParams <- params[drawIndex,]
    
    simSummaries <- summaries[-drawIndex,]
    simParams <- params[-drawIndex,]
    
    distances <- apply(simSummaries, MARGIN = 1, function(x) return(sqrt(sum((x-obsSummaries)^2))))
    selection <- distances <= quantile(distances, probs = eQuantile)
    
    acceptedParams <- simParams[selection,]
    estimatedParams <- apply(acceptedParams, 2, sumFun)
    
    out <- c(obsParams, estimatedParams)
    namesParams <- colnames(params)
    names(out) <- c(paste(namesParams, "_true", sep = ""), paste(namesParams, "_estimated", sep = ""))
    
    out
  }
  parallel::stopCluster(cl)
  return(OUT)
}

#' @title random forest parameter prediction
#' @description Function to predict model parameters on the basis of summary statistics. This function is used to create a reduced set of summary statistics by fitting a random forest with p as response and s as predictors. The predicted parameters can then be used in the ABC rejection sampling function instead of the summaries.
#' @param p matrix or dataframe containing the parameters as columns
#' @param s matrix or dataframe containing the summary statistics as columns
#' @param nCores number of parallel cores
#' @param nTrees number of trees to construct
#' @param scale logical, determining whether the resulting predictions should be scaled by their standard deviation
#' @param ... further arguments passed to ranger::ranger
#' @return A list with two items. Item one contains the true parameters, item two contains the predicted parameters.
#' @author Tankred Ott
summaryRedRF <- function (p, s, nCores="auto", nTrees=500, scale = TRUE, ...) {
  checkSuggested("ranger")
  nCores <- getValidatedNCores(nCores)
  
  s <- as.data.frame(s)
  p <- as.data.frame(p)
  
  nSum <- ncol(s)
  nPar <- ncol(p)
  
  d <- cbind(s, p)
  d <- na.omit(d)
  nRows <- nrow(d)
  
  pNames <- colnames(p)
  sNames <- colnames(s)
  formulas <- paste(paste(pNames, " ~ "), paste(sNames, collapse = " + "))
  
  pPredicted <- matrix(NA, nrow = nRows, ncol = nPar)
  for (i in 1:length(formulas)) {
    rf <- ranger::ranger(formula = formulas[i], data = d, num.trees = nTrees, num.threads = nCores, ...)
    pPredicted[,i] <- rf$predictions
  }
  
  if (scale == TRUE) {
    pPredSd <- apply(X = pPredicted, MARGIN = 2, FUN=sd)
    pPredicted <- scale(pPredicted, scale = pPredSd, center = F)
  }
  
  return(list(pTrue = as.matrix(d[,(nSum + 1):(nSum+nPar)]),
              pPred = as.matrix(pPredicted)))
}


#' @title saABC parameter prediction
#' @description Function to predict model parameters on the basis of summary statistics. This function is used to create a reduced set of summary statistics using saABC with p as response and s as predictors. The predicted parameters can then be used in the ABC rejection sampling function instead of the summaries.
#' @param p matrix or dataframe containing the parameters as columns
#' @param s matrix or dataframe containing the summary statistics as columns
#' @param scale logical, determining whether the resulting predictions should be scaled by their standard deviation
#' @param ... currently ignored
#' @return Matrix of parameter predictions
#' @author Tankred Ott
summaryRedsaABC <- function (p, s, scale = TRUE, ...) {
  checkSuggested("abctools")
  weights <- abctools::saABC(p, s, plot = FALSE)$B
  x <- rep(TRUE, ncol(s)) # fixes #129 (x undefined if no weight is NA)
  if (any(is.na(weights))) {
    x <- apply(X = weights, MARGIN = 2, FUN = function(x) return(!any(is.na(x))))
    warning(paste("Linear regression problems: ill-conditioning? Ignored summaries leading to NAs! (", paste(names(x)[!x], collapse = ", "), ")"))
    weights <- abctools::saABC(p, s[, x], plot = FALSE)$B
  }
  
  pPred <- as.matrix(s[, x]) %*% t(weights)
  
  if (scale == TRUE) {
    pPredSd <- apply(X = pPred, MARGIN = 2, FUN=sd)
    pPred <- scale(pPred, scale = pPredSd, center = F)
  }
  
  return(pPred)
}
