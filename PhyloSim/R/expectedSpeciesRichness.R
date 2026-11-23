#' @title Expected species richness
#' @description Calculates the expected species richness of the neutral model with global dispersal in equilibrium, from the Ewens sampling formula.
#' @param specRate the speciation rate in speciations per generation
#' @param dimensions the length of the grid edge of a square grid, or the grid dimensions \code{c(x, y)}
#' @details The function is only valid for neutral models with global dispersal. With \eqn{J} cells and \code{specRate} new species per generation of \eqn{J} death-birth events, the speciation probability per birth is \eqn{\nu = specRate / J}{nu = specRate / J}. For this Moran model, the expected number of species is \deqn{E(S) = \sum_{i=0}^{J-1} \frac{\theta}{\theta + i}, \quad \theta = \frac{(J - 1) \nu}{1 - \nu} = \frac{specRate (J - 1)}{J - specRate}}{E(S) = sum_{i = 0}^{J - 1} theta / (theta + i), theta = (J - 1) nu / (1 - nu) = specRate (J - 1) / (J - specRate)} (Vallade and Houchmandzadeh 2003, eqns 6 and 7). \eqn{\theta}{theta} is about \code{specRate}, half of the fundamental biodiversity number \eqn{2 J \nu}{2 J nu} of Hubbell (2001).
#'
#' Simulations with \code{\link{runSimulation}} have 7 to 15% (0.5 to 2.3 species) more species than this expectation for \code{specRate} = 1 to 4 on 10 x 10 to 30 x 30 cells, probably because the model adds the new species at the end of each generation, right before the species are recorded.
#' @references Hubbell, S. P. (2001) The unified neutral theory of biodiversity and biogeography. Princeton University Press.
#'
#' Vallade, M. and Houchmandzadeh, B. (2003) Analytical solution of a neutral model of biodiversity. Physical Review E 68, 061902. \doi{10.1103/PhysRevE.68.061902}
#' @return A float value for the expected species richness
#' @export
#' @example inst/examples/expectedSpeciesRichness-help.R
es <- function(specRate, dimensions){
  checkmate::assert_number(specRate, lower = 0, finite = TRUE)
  checkmate::assert_integerish(dimensions, lower = 1, any.missing = FALSE, min.len = 1, max.len = 2)
  J <- if (length(dimensions) == 2) prod(dimensions) else dimensions^2
  if (specRate >= J) return(J) # every birth is a new species (specRate > J is not a valid model)
  theta <- specRate * (J - 1) / (J - specRate) # Moran model, not 2 * specRate (fixes #124)
  if (theta == 0) return(1)
  sum(theta / (theta + 0:(J - 1)))
}
