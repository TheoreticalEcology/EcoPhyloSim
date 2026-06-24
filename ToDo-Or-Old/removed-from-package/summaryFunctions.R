
require(ape)
require(caper)
require(geiger)
require(picante)
require(apTreeshape)
# require(dispRity)
require(RPANDA)
require(nLTT)
require(stringr)
require(Rfast)
require(castor)


#' Calculate summary statistics from a phylogenetic tree
#' 
#' @param object an object of class phylo 
#' @param selection selection of metrics. 1:10 is fast, > 11 is slow. Max range is 1:13
#' @param dropFossils whether to apply ape::drop.fossil. Safer, but can be slow for large trees. 
#' 
#' @author Florian Hartig, modified from a function by Allen Hulbert
#' 
#' @details 
#' 
#' See also 
#' 
#' https://www.lirmm.fr/MEP05/talk/14_Mooers.pdf
#' 
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5790134/
#' 
#' 
treeMetrics = function(phy, selection = 1:10, dropFossils = FALSE) {
  
  if(class(phy) != "phylo") {
    stop("The object passed to treeMetrics() is not of class 'phylo'")
  }
  
  out = list()
  
  ##### Tree preparation #####
  
  # Drop root edge
  phy$root.edge = 0
  
  # Prune out extinct species
  if(dropFossils == TRUE){
    phy = ape::drop.fossil(phy)    
  }

  # Randomly resolve any existing polytomies
  tree = ape::multi2di(phy)
  
  ##### Calculation of Statistics #####
  
  # Absolute tree length, which is the diagonal of the vcv matrix
  v.matrix = ape::vcv(tree, corr=F)
  tree.length = diag(v.matrix)[1]
  
  # Tree with branch lengths scaled by total tree length
  tree.scaled = tree
  tree.scaled$edge.length = tree$edge.length/tree.length
  
  # Create treeshape objects for certain calculations
  tree.shape = apTreeshape::as.treeshape(tree)
  tree.shape.scaled = apTreeshape::as.treeshape(tree.scaled)
  
  # using castor package
  # --get_all_distances_to_root returns tip and node distances; only want tip distances which are returned first
  root.dist <- castor::get_all_distances_to_root(tree, as_edge_count = TRUE)[1:(tree$Nnode+1)]
  
  # Species richness
  if(1 %in% selection){
    out$richness = length(tree$tip.label) 
  }
  
  # Length
  if(1 %in% selection){
    out$length = length(tree$tip.label) 
  }
  
  # TODO: why not use pd(samp, tree, include.root=TRUE)
  # PD (phylogenetic diversity, summed branch lengths) on scaled tree
  if(3 %in% selection){
    out$PD = sum(tree.scaled$edge.length)
  }

  # Colless index  
  if(4 %in% selection){
    out$Colless = apTreeshape::colless(tree.shape)
  }

  # Sackin index  
  if(5 %in% selection){
    out$Sackin = apTreeshape::sackin(tree.shape)
  }

  #Pybus & Harvey (2000)'s gamma statistic  
  if(6 %in% selection){
    out$gamma.stat = ape::gammaStat(tree.scaled)
  }
  
  # logarithm of the ratio of the likelihoods under the Yule model and the PDA model
  if(7 %in% selection){
    out$Yule.PDA.ratio = apTreeshape::shape.statistic(tree.shape)
  }
  
  # Mean Root Distance (MRD, also abbreviated N-bar, Shao & Sokal 1990) (Agapow & Purvis 2002)
  if(8 %in% selection){
    out$MRD = mean(root.dist)
  }
  
  # Variance in root distances across the tips, Shao & Sokal 1990 (Agapow & Purvis 2002)
  if(9 %in% selection){
    out$MRD = var(root.dist)
  }
  
  
  # Calculate Helmus et al. 2007's Phylogenetic Species Variability (PSV) metric
  # based on the phylogenetic variance/covariance matrix
  # While Helmus et al.'s equation is as follows
  # PSV = (n*trace(C) - sum(C))/(n*(n-1)) 
  # this only holds for a tree with contemporaneous tips where C is the correlation matrix.
  # When the tree is not ultra-metric, Algar et al. 2009 (Ecol Lett) show the calculation
  # using the variance-covariance matrix, V:
  # PSV = (n*trace(V) - sum(V))/(trace(V)*(n-1)) 
  
  if(10 %in% selection){
    n = nrow(v.matrix)
    out$PSV = (n*sum(diag(v.matrix)) - sum(v.matrix))/(sum(diag(v.matrix))*(n-1))
  }
  
  # mean I' from Purvis et al. 2002. Evaluating phylogenetic tree shape: two modifications to Fusco & Cronk's method
  
  if(11 %in% selection){
    if(tree$Nnode >= 3) {
      tree$node.label = NULL
      fusco = caper::fusco.test(tree)
      out$mean.Iprime = fusco$mean.Iprime
    } else {
      out$mean.Iprime = NA
    }
  }
  

  if(12 %in% selection){  
    # Mean Pairwise Distance (in scaled tree)
    tryCatch({
      pairwise.dists = ape::dist.nodes(tree.scaled)
      sum.pairwise.dists = Rfast::upper_tri(pairwise.dists, diag = FALSE, suma = TRUE)
      MPD = sum.pairwise.dists/(ncol(pairwise.dists)*(ncol(pairwise.dists)-1)/2)
    }, error = function(e) {
      MPD = NA
    })
    out$MPD = MPD
  }
  
  
  
  # Variance Pairwise Distance (in scaled tree)
  if(13 %in% selection){
    tryCatch({
      tri.pairwise.dists = Rfast::upper_tri(pairwise.dists, diag = FALSE, suma = FALSE)
      VPD = var(as.vector(tri.pairwise.dists))
    }, error = function(e) {
      VPD = NA
    })
    out$VPD = VPD
  }
  
  
  # beta
  # Calculate Blum & Francois (2006)'s Beta metric of tree imbalance using apTreeshape package
  
  # tryCatch({
  #   beta.out = maxlik.betasplit.AH(tree.scaled)
  #   beta.stat = beta.out$max_lik
  #   
  # }, error = function(e) {
  #   beta.stat = NA
  # })
  

  

  
  # Some metrics are difficult to calculate on very large trees.
  # Also calculate RPANDA spectral density metrics (Lewitus & Morlon 2016)
  
  # RPANDA analyses seem to bonk on very large phylogenies, so only try calculating for fewer than 6000 species
  # if(S < 6000) {
  #   
  #   # RPANDA
  #   tryCatch({
  #     MGL = spectR(tree)
  #     MGL_principal_eigenvalue = MGL$principal_eigenvalue 
  #     MGL_asymmetry = MGL$asymmetry  
  #     MGL_peakedness = MGL$peakedness
  #     MGL_eigengap = MGL$eigengap
  #   }, error = function(e) {
  #     MGL_principal_eigenvalue = NA
  #     MGL_asymmetry = NA
  #     MGL_peakedness = NA
  #     MGL_eigengap = NA
  #   })
  #   
  # } else {
  #   MGL_principal_eigenvalue = NA
  #   MGL_asymmetry = NA  
  #   MGL_peakedness = NA
  #   MGL_eigengap = NA
  # }

  
  return(out)
}


