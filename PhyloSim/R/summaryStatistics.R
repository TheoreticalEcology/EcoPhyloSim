
#' @title Summary Statistics
#' @description Calculates the summary statistics for a simulation of class PhyloSim.
#' @param simulation Simulation output of the class "PhyloSim"
#' @param strict Logical. If TRUE, the output of the function is NA if at least one summary statistic could not be calculated. Default is FALSE.
#' @param which.result Integer, determines which result should be used. This argument is only useful if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used.
#' @details Computes the following summary statistics from one result of a simulation run. If calculateSummaries = TRUE in \code{\link{createCompletePar}}, \code{\link{runSimulation}} stores the summary statistics of the last result in `Output[[1]]$summaries`:
#' \describe{
#'   \item{racAuc}{Area under the Rank Abundance Curve (see \code{\link{rac}}).}
#'   \item{sacAuc}{Area under the Species Area Curve (see \code{\link{sac}}).}
#'   \item{alphaDiversity}{Mean local (alpha) species richness.}
#'   \item{betaDiversity}{Mean beta diversity (see \code{\link[vegan]{betadiver}}).}
#'   \item{gammaDiversity}{Total species richness of the community.}
#'   \item{imbalance}{Yule-normalized Colless index of phylogenetic imbalance.}
#'   \item{dispersion}{Standardized phylogenetic dispersion (SES-MPD z-score).}
#'   \item{gammaStatistics}{Gamma statistic (see \code{\link[ape]{gammaStat}}) testing for deviation from a constant diversification rate.}
#'   \item{meanNodeAge}{Mean branch length of the phylogeny.}
#'   \item{varPart_1}{Variance in community composition explained by environment (see \code{\link[vegan]{varpart}}).}
#'   \item{varPart_2}{Variance in community composition explained by competition (see \code{\link[vegan]{varpart}}).}
#'   \item{ce}{Clark-Evans index of spatial clustering.}
#'   \item{nb8}{Index of local species clustering among neighboring grid cells.}
#'   \item{envSorting}{Share of the variance of the environment across all cells that is explained by species identity, 1 - sum((e - e_s)^2) / sum((e - mean(e))^2), where e is the environment of a cell and e_s the mean environment of the cells occupied by its species. 0 if all species occupy the same mix of environments, 1 if each species occupies a single environmental value; NA if the environment is constant. Without environmental sorting, the expected value is about (S - 1) / (N - 1) for S species in N cells, so it increases with the number of species; with local dispersal, spatial clustering alone gives values above 0.}
#' }
#' @return A list containing the summary statistics, or NA if strict = TRUE and at least one statistic could not be calculated. If the simulation contains only one species, a list of NA is returned.
#' @example /inst/examples/summaryStatistics-help.R
#' @export
#' 

calculateSummaryStatistics <- function(simulation, strict=FALSE, which.result=NULL) {
  summaryStatistics <- list(racAuc=NA,
                            sacAuc=NA,
                            alphaDiversity = NA,
                            betaDiversity = NA,
                            gammaDiversity = NA,
                            imbalance = NA,
                            dispersion = NA,
                            gammaStatistics = NA,
                            meanNodeAge = NA,
                            varPart_1 = NA,
                            varPart_2 = NA,
                            ce = NA,
                            nb8 = NA,
                            envSorting = NA)
  
  # use one result for all statistics (fixes #134)
  checkmate::assert_int(which.result, lower = 1, upper = length(simulation$Output), null.ok = TRUE)
  if (is.null(which.result)) which.result <- length(simulation$Output)
  result <- simulation$Output[[which.result]]
  
  # TODO: implement special case: only 1 species
  if (is.double(result$phylogeny)) {
    return(summaryStatistics)
  }
  
  # area under the curve by the trapezoidal rule, as MESS::auc() (#167)
  auc <- function(x, y) {
    if (length(x) < 2) return(NA)
    o <- order(x)
    sum(diff(x[o]) * (y[o][-1] + y[o][-length(y)])) / 2
  }

  # RAC
  rac <- rac(simulation, which.result = which.result, plot = F)
  summaryStatistics$racAuc <- auc(x = rac$Rank,
                                  y = rac$Abundance)
  
  # SAC
  sac <- sac(simulation, which.result = which.result, size = seq(1,15), plot = F)
  summaryStatistics$sacAuc <- auc(x = sac$size,
                                  y = sac$sr.Mean)
  
  
  # alpha and beta diversity
  
  # Draw random plots
  # TODO: let the user chose the plot size and number of plots
  nSubPlots <- 20
  sizeSubPlots <- 10
  
  plots <- localPlots(simulation,
                      which.result = which.result,
                      n = nSubPlots,
                      size = sizeSubPlots,
                      plot = F,
                      community = T)
  
  # alpha diversity
  # in this case alpha diversity = mean number of species per subplot
  summaryStatistics$alphaDiversity <- mean(rowSums(plots$communityTable>0))
  
  # beta diversity
  # uses vegan's beta diversity function
  # see betadiver(help=TRUE)
  # TODO: let the user chose the method
  method = 1
  summaryStatistics$betaDiversity <- mean(vegan::betadiver(plots$communityTable,
                                                           method=method))
  
  # gamma diversity
  summaryStatistics$gammaDiversity <- gammaDiversity(result)
  
  # phylogenetic imbalance
  # Colless index (Yule-normalized), computed internally (see collessImbalance.R).
  # Previously relied on apTreeshape, which was archived from CRAN.
  summaryStatistics$imbalance <- .collessIndex(
    ape::drop.fossil(result$phylogeny),
    norm = "yule")

  # phylogenetic dispersion
  # standardized mean pairwise distance (z-score, analogous to -NRI).
  # Previously used PhyloMeasures::mpd.query (archived from CRAN); now uses
  # picante::ses.mpd, which is maintained and already a dependency.
  # The species pool is all species in the phylogeny: only the extant species with
  # prunePhylogeny = TRUE, all species ever created with prunePhylogeny = FALSE (#147).
  sesMpd <- picante::ses.mpd(plots$communityTable,
                             stats::cophenetic(result$phylogeny),
                             null.model = "taxa.labels", runs = 999)
  summaryStatistics$dispersion <- mean(sesMpd$mpd.obs.z, na.rm = TRUE)
  
  # gamma statistics
  # uses ape's implementation of the gammaStatistics
  # attention: only works on ultrametric trees
  summaryStatistics$gammaStatistics <- ape::gammaStat(ape::drop.fossil(result$phylogeny))
  
  # mean node age
  #summaryStatistics$meanNodeAge <- ape::chronoMPL(simulation$Output[[1]]$phylogeny)
  summaryStatistics$meanNodeAge <- mean(ape::drop.fossil(result$phylogeny)$edge.length)
  
  
  # variation partitioning
  #distance_matrix <- distances(plots$positions, limits=c(nrow(simulation$Output[[1]]$specMat), ncol(simulation$Output[[1]]$specMat)))
  env <- data.frame(unlist(lapply(plots$envPlots, mean)))
  comp <- data.frame(unlist(lapply(plots$compPlots, mean)))
  
  #mod <- vegan::varpart(plots$communityTable, as.dist(distance_matrix), env, comp)
  mod <- vegan::varpart(plots$communityTable, env, comp)
  
  summaryStatistics$varPart_1 <- mod$part$fract$Adj.R.squared[1]
  summaryStatistics$varPart_2 <- mod$part$fract$Adj.R.squared[2]
  
  
  # Clark Evans
  summaryStatistics$ce <- clarkEvans(result$specMat, correction="Donnelly")
  
  # nb8
  summaryStatistics$nb8 <- nb8Idx(result$specMat)
  
  # environmental sorting of species (replaces GLM coefficients that were NaN under strong filtering, #126)
  summaryStatistics$envSorting <- envSorting(result)
  
  if (strict == TRUE) {
    if (sum(is.na(summaryStatistics)) > 0) {
      summaryStatistics <- NA
    }
  }
  
  return(summaryStatistics)
}

gammaDiversity <- function(result, q=0){
  spec_vector <- c(result$specMat)
  specs <- unique(spec_vector)
  n_specs <- length(specs)
  
  if (q == 0) return(n_specs)
  
  n_individuals <- length(spec_vector)
  
  prop_abundances <- rep(NA, length(specs))
  
  for(i in 1:length(specs)){
    prop_abundances[i] <- sum(spec_vector == specs[i]) / n_individuals
  }
  
  gamma <- prop_abundances * prop_abundances^(q-1)
  gamma <- sum(gamma)^(1/(q-1))
  gamma <- 1/gamma
  return(gamma)
}

distances <- function(positions, limits){
  distance <- function(a, b, limits){
    # sqrt(min(|x1 - x2|, w - |x1 - x2|)^2 + min(|y1 - y2|, h - |y1-y2|)^2)
    return(sqrt(min(abs(a[1]-b[1]), limits[1] - abs(a[1]-b[1]))^2 + min(abs(a[2]-b[2]), limits[2] - abs(a[2]-b[2]))^2))
  }
  dist_mat <- matrix(0, ncol=length(positions), nrow=length(positions))
  # TODO: maybe there is a better way to implement this instead of a for loop
  for(i in 1:nrow(dist_mat)){
    for(j in 1:ncol(dist_mat)){
      dist_mat[i,j] <- distance(positions[[i]], positions[[j]], limits=limits)
    }
  }
  return(dist_mat)
}

clarkEvans <- function (mat, ...) {
  nCol <- ncol(mat)
  nRow <- nrow(mat)
  
  ids <- unique(as.vector(mat))
  ceSum <- 0
  n <- 0
  
  for (i in 1:length(ids)) {
    pos <- which(mat == ids[i])
    col <- ceiling(pos / nRow)
    row <- (pos - 1) %% nRow + 1 # fixes #163 (wrong row index)
    pp <- spatstat.geom::ppp(x=col, y=row, window=spatstat.geom::owin(xrange = c(0, nCol), yrange = c(0, nRow)))
    ce <- spatstat.explore::clarkevans(pp, ...)
    if (!is.infinite(ce) && !is.na(ce) && !is.nan(ce)) {
      ceSum <- ceSum + ce
      n <- n + 1
    }
  }
  if (n == 0) return(NA_real_)
  return(ceSum/n)
}

getNeighbours8 <- function(x, y, mat, nCol=NULL, nRow=NULL, isTorus=TRUE) {
  if (is.null(nCol)) nCol <- ncol(mat)
  if (is.null(nRow)) nRow <- nrow(mat)
  
  neighbours <- vector("numeric")
  
  for (yOffset in c(-1, 0, 1)) {
    for (xOffset in c(-1, 0, 1)) {
      if (!(yOffset == 0 && xOffset == 0)) {
        xPos <- x + xOffset
        yPos <- y + yOffset
        
        if (xPos > nCol || xPos < 1 || yPos > nRow || yPos < 1) {
          if (isTorus != TRUE) next
        }
        
        if (xPos > nCol) {
          xPos <- xPos - nCol
        } else if (xPos < 1) {
          xPos <- nCol - xPos
        }
        
        if (yPos > nRow) {
          yPos <- yPos - nRow
        } else if (yPos < 1) {
          yPos <- nRow - yPos
        }
        
        neighbours <- append(neighbours, mat[yPos, xPos])
      }
    }
  }
  return(neighbours)
}

nb8Idx <- function (mat) {
  nCol <- ncol(mat)
  nRow <- nrow(mat)
  n <- nRow * nCol
  
  ids <- unique(as.vector(mat))
  
  id_freqs <- vector(mode = "numeric", length = length(ids))
  for (i in 1:length(ids)) {
    id_freqs[i] <- sum(mat == ids[i])
  }
  
  # print(nRow)
  # print(nCol)
  
  nbRates <- vector("numeric", length(ids))
  names(nbRates) <- ids
  
  for (row in 1:nRow) {
    for (col in 1:nCol) {
      id <- mat[row, col]
      nbrs <- getNeighbours8(col, row, mat, nRow = nRow, nCol = nCol)
      
      r <- sum(id == nbrs) / length(nbrs)
      # r <- r * (id_freqs[which(ids == id)] / n)
      # r <- r * (n / id_freqs[which(ids == id)])
      r <- r / id_freqs[which(ids == id)]
      nbRates[as.character(id)] <- nbRates[as.character(id)] + r
    }
  }
  
  # return(sum(nbRates) / length(nbRates))
  return(sum(nbRates))
}

envSorting <- function (result) {
  env <- as.vector(result$envMat)
  totalSS <- sum((env - mean(env))^2)
  if (length(unique(env)) < 2) return(NA_real_) # constant environment
  1 - sum((env - stats::ave(env, as.vector(result$specMat)))^2) / totalSS
}


