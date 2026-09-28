#' @keywords internal
#' @title Analysing phylogenetic and biogeographic patterns and inferring species community assembly 
#' @name PhyloSim
#' @useDynLib PhyloSim, .registration = TRUE
#' @description
#' This package is designed to simulate and analyse species community assembly processes. It allows the user to implement different mechanisms of community assembly such as environmental filtering, competitive exclusion or dispersal limitation. These mechanisms can be combined or activated seprately, which allows the user to identify biogeographic or phylogenetic patterns exclusively correlated with a certain mechanism or combination of mechanisms.
#' @details
#' The core of this package is basically a community assembly model, which can be run with \code{\link{runSimulation}} or \code{\link{runSimulationBatch}}. All other functions help analize and visualize the model outcome. For visualization the following functions can be useful:
#'\tabular{ll}{
#'  \code{\link{plotSpatialPhylo}} \tab Plots the phylogeny tree ,the spatially abundance of the evolved species and their traits.\cr
#'  \code{\link{rac}} \tab Plots the Rank Abundance Curve for a given community\cr
#'  \code{\link{sac}} \tab Plots the species area curve for a given community\cr
#'  \code{\link{highlightLocalPhylo}} \tab Visualises the phylogenetic clades of a local community within a metacommunity.\cr
#'  \code{\link{plotTraitDistribution}} \tab Plots trait-histograms, trait ~ Environment relationship and the spatial distribution of the different traits\cr
#'}
#'
#'
#'The function \code{\link{phyloReconstruct}} creates a phylogeny from given neutral traits of a community in order to compare them with the existing phylogeny. The observed patterns in the simulation output can be tested against null models using the function \code{\link{nullModel}}.
#'
#'
#' @author Paul Bauche, Florian Hartig
#' @importFrom grDevices colorRampPalette dev.size grey rainbow rgb
#' @importFrom graphics barplot hist image layout lines mtext par plot.new polygon text title
#' @importFrom stats binomial cophenetic density dist ecdf glm median na.omit quantile rmultinom sd
#' 
#' @example /inst/examples/phylosim-help.R
#'  
"_PACKAGE"
## usethis namespace: start
## usethis namespace: end
NULL