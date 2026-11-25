#' @title Expected species richness
#' @description Calculates the expected species richness of the neutral model with global dispersal in equilibrium, from the Ewens sampling formula.
#' @param specRate the speciation rate in speciations per generation, at most the number of cells
#' @param dimensions the length of the grid edge of a square grid, or the grid dimensions \code{c(x, y)}
#' @details The function is only valid for neutral models with global dispersal and without protracted speciation. In the Moran model with \eqn{J} cells, each death-birth event replaces one individual, and the offspring founds a new species with probability \eqn{\nu}{nu}; with \code{specRate} new species per generation of \eqn{J} events, \eqn{\nu = specRate / J}{nu = specRate / J}. The expected number of species is \deqn{E(S) = \sum_{i=0}^{J-1} \frac{\theta}{\theta + i}, \quad \theta = \frac{(J - 1) \nu}{1 - \nu} = \frac{specRate (J - 1)}{J - specRate}}{E(S) = sum_{i = 0}^{J - 1} theta / (theta + i), theta = (J - 1) nu / (1 - nu) = specRate (J - 1) / (J - specRate)} (Vallade and Houchmandzadeh 2003, eqns 6 and 7). \eqn{\theta}{theta} is about \code{specRate}, half of the fundamental biodiversity number \eqn{2 J \nu}{2 J nu} of Hubbell (2001).
#'
#' \code{\link{runSimulation}} instead turns on average \code{specRate} random individuals into new species at the end of each generation. Apart from this timing, its \eqn{\theta}{theta} is \eqn{specRate (J - 1) / J}{specRate (J - 1) / J}, which gives up to 3% fewer species than \code{es()} for \code{specRate} = 1 to 4 on 10 x 10 to 30 x 30 cells. Simulations in these settings have 7 to 15% (0.5 to 2.3 species) more species than \code{es()}, probably because the species are recorded right after the new species of a generation are added.
#' @references Hubbell, S. P. (2001) The unified neutral theory of biodiversity and biogeography. Princeton University Press.
#'
#' Vallade, M. and Houchmandzadeh, B. (2003) Analytical solution of a neutral model of biodiversity. Physical Review E 68, 061902. \doi{10.1103/PhysRevE.68.061902}
#' @return the expected number of species
#' @export
#' @example inst/examples/expectedSpeciesRichness-help.R
es <- function(specRate, dimensions){
  checkmate::assert_integerish(dimensions, lower = 1, any.missing = FALSE, min.len = 1, max.len = 2)
  J <- prod(rep_len(dimensions, 2))
  checkmate::assert_number(specRate, lower = 0, upper = J)
  if (specRate == J) return(J) # every birth is a new species
  theta <- specRate * (J - 1) / (J - specRate) # Moran model, not 2 * specRate (fixes #124)
  1 + sum(theta / (theta + seq_len(J - 1)))
}
