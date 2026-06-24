#' @title Colless' imbalance
#' @description Calculates the Colless' imbalance for a Phylogeny.
#' @param simu An object of type "PhyloSim"
#' @param which.result Integer, determines which result should be used. This argument is only usefull if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used.
#' @param useApTreeshape Deprecated and ignored. The Colless index is now computed internally; the (archived) apTreeshape package is no longer required. Kept only for backward compatibility.
#' @param norm Normalization: Either NULL (raw Colless index), "yule" or "pda".
#' @param dropFossils Boolean, if TRUE applies ape's drop.fossil on the phylogeny
#' @return A numeric value for the Colless' Imbalance
#' @details Multifurcations are randomly resolved with \code{\link[ape]{multi2di}} before the index is computed, because the Colless index is only defined for binary trees. If dropFossils == TRUE only extant species are included in the phylogeny.
#' @references Colless, D. H. "Review of phylogenetics: the theory and practice of phylogenetic systematics." Syst. Zool 31 (1982): 100-104.
#' @export

collessImbalance <- function(simu, which.result = NULL, useApTreeshape = FALSE, norm = NULL, dropFossils = FALSE){
  if (is.null(which.result)) which.result = length(simu$Output)
  phylo <- simu$Output[[which.result]]$phylogeny

  if (dropFossils == TRUE) {
    phylo <- ape::drop.fossil(phylo)
  }

  if (isTRUE(useApTreeshape)) {
    warning("'useApTreeshape' is deprecated: apTreeshape was archived from CRAN. ",
            "The Colless index is now computed internally.", call. = FALSE)
  }

  return(.collessIndex(phylo, norm = norm))
}


# Internal helper: Colless' imbalance index computed directly from an ape
# 'phylo' object, reproducing the definitions used by the (archived)
# apTreeshape::colless function (norm = NULL, "yule" or "pda").
.collessIndex <- function(phy, norm = NULL) {
  if (!inherits(phy, "phylo"))
    stop("object \"phy\" is not of class \"phylo\"")

  # The Colless index is only defined for binary trees: resolve polytomies.
  phy <- ape::multi2di(phy)

  n <- length(phy$tip.label)
  if (n < 2) return(NA_real_)

  # Number of descending tips for every node (tips themselves count as 1).
  nTips <- ape::node.depth(phy)

  internal <- (n + 1):(n + phy$Nnode)
  ic <- 0
  for (nd in internal) {
    children <- phy$edge[phy$edge[, 1] == nd, 2]
    if (length(children) == 2L) {
      ic <- ic + abs(nTips[children[1]] - nTips[children[2]])
    }
  }

  if (is.null(norm)) return(as.numeric(ic))

  if (identical(norm, "yule")) {
    # Expected Colless index under the Yule model (Euler-Mascheroni constant).
    EICN <- n * log(n) + (0.5772156649015329 - 1 - log(2)) * n
    return(as.numeric((ic - EICN) / n))
  }

  if (identical(norm, "pda")) {
    return(as.numeric(ic / n^(1.5)))
  }

  stop("'norm' must be NULL, \"yule\" or \"pda\"")
}
