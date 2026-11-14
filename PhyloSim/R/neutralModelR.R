#' @title Neutral model in R
#' @description A simple neutral model of community assembly (Hubbell 2001), written in R, for teaching and testing. It is independent of the C++ model of \code{\link{runSimulation}} and returns an object of class \code{PhyloSim} with a species matrix only.
#' @param x Integer, dimension of the landscape in x-direction; the grid must have at least 2 cells (x * y >= 2)
#' @param y Integer, dimension of the landscape in y-direction
#' @param dispersal \code{"global"}, or a function with one argument \code{n} that returns \code{n} non-negative dispersal distances in cells (see Details)
#' @param runs Integer or increasing vector of integers, number of generations. If a vector, the species matrix is also saved at the intermediate generations, as in \code{\link{createCompletePar}}.
#' @param specRate Numeric, mean number of speciation events per generation (at most \code{x * y})
#' @param density Numeric between 0 and 1, strength of conspecific negative density dependence. 0 (default) gives a neutral model.
#' @param densityCut Integer, radius of the neighborhood for density dependence. If \code{density > 0}, at most \code{(min(x, y) - 1) / 2}, so that no cell is counted twice.
#' @param seed Integer, random seed. If NULL (default), a seed is drawn. Apart from this draw, the random number stream of the R session is not changed.
#' @param verbose Logical, if TRUE, reports progress every 100 generations and the runtime
#' @details The landscape is a torus of \code{x} times \code{y} cells, initially occupied by a single species. One generation consists of \code{x * y} events. In each event, a random cell dies and is replaced by a recruit. With probability \code{specRate / (x * y)}, the recruit is a new species; otherwise it is the offspring of a parent cell. As in \code{\link{runSimulation}}, the parent is never the dead cell itself: it is a random other cell (\code{dispersal = "global"}), or the cell at a distance drawn from \code{dispersal} in a random direction, rounded to the grid, where draws that land on the dead cell are repeated. With global dispersal, the equilibrium species richness follows the Ewens sampling formula with \eqn{\theta = specRate (x y - 1) / (x y - specRate)}{theta = specRate * (x * y - 1) / (x * y - specRate)} (from the coalescent of this model), which is approximately \code{specRate} if \code{specRate} is much smaller than \code{x * y}.
#'
#' With \code{density > 0}, a recruit is rejected with probability \code{density} times the proportion of conspecifics among the cells within distance \code{densityCut}, and a new parent is drawn; at most 10 parents are drawn, and the last recruit is accepted. Note that \code{dispersal} and \code{density} are defined differently than in \code{\link{createCompletePar}}, where density dependence acts through the competition trait. Speciation also differs: \code{\link{runSimulation}} adds new species in a batch at the end of each generation.
#'
#' The model is slow and only suitable for small landscapes. The output contains the species matrix (\code{specMat}) and the parameters (\code{Model}), so that functions that only need the species matrix, such as \code{\link{specRich}} and \code{\link{rac}}, can be used (\code{\link{sac}} and \code{\link{localPlots}} currently only on square grids). There is no phylogeny and there are no traits, so functions that need them (e.g. \code{plot()} and \code{summary()}) do not work.
#' @references Hubbell, Stephen P. The unified neutral theory of biodiversity and biogeography (MPB-32). Vol. 32. Princeton University Press, 2001.
#' @return An object of class \code{PhyloSim}, see \code{\link{runSimulation}}. Each element of \code{Output} contains only \code{specMat}, a matrix of species IDs with \code{y} rows and \code{x} columns.
#' @example /inst/examples/neutralModelR-help.R
#' @export
NeutralMod <- function(x = 50, y = 50, dispersal = "global", runs = 100, specRate = 1,
                       density = 0, densityCut = 1, seed = NULL, verbose = TRUE) {

  checkmate::assert_count(x, positive = TRUE)
  checkmate::assert_count(y, positive = TRUE)
  if (x * y < 2) stop("The grid must have at least 2 cells (x * y >= 2).")
  if (!is.function(dispersal) && !identical(dispersal, "global"))
    stop("dispersal must be \"global\" or a function(n) returning n distances")
  checkmate::assert_integerish(runs, lower = 1, any.missing = FALSE, min.len = 1, sorted = TRUE, unique = TRUE)
  checkmate::assert_number(specRate, lower = 0, upper = x * y)
  checkmate::assert_number(density, lower = 0, upper = 1)
  checkmate::assert_count(densityCut, positive = TRUE)
  if (density > 0) checkmate::assert_int(densityCut, upper = (min(x, y) - 1) %/% 2)
  checkmate::assert_int(seed, null.ok = TRUE)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  oldSeed <- globalenv()$.Random.seed
  on.exit(if (is.null(oldSeed)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", oldSeed, envir = globalenv()))
  set.seed(seed)

  ptm <- proc.time()
  J <- x * y
  mat <- matrix(1L, y, x)
  counter <- 1L
  if (density > 0) neighbors <- as.matrix(getCircularOffsets(densityCut))
  maxTries <- 10

  # linear index of the cell (row, col) + (dr, dc) on the torus
  shift <- function(cell, dr, dc) {
    r <- (cell - 1) %% y
    c <- (cell - 1) %/% y
    ((r + dr) %% y) + ((c + dc) %% x) * y + 1
  }

  # NeutralMod() follows the rules of the C++ model of runSimulation() except for the
  # differences listed in the help (Details), e.g. the dead cell is never its own parent (#81)
  drawParent <- function(cell) {
    if (!is.function(dispersal)) {
      parent <- sample.int(J - 1, 1)
      return(parent + (parent >= cell))
    }
    for (k in 1:100) {
      d <- dispersal(1)
      checkmate::assert_number(d, lower = 0, finite = TRUE, .var.name = "dispersal(1)")
      a <- stats::runif(1, 0, 2 * pi)
      dr <- round(d * sin(a))
      dc <- round(d * cos(a))
      parent <- shift(cell, dr, dc)
      if (parent != cell) return(parent)
    }
    stop("dispersal keeps landing on the dead cell; use larger distances")
  }

  output <- list()
  for (i in 1:max(runs)) {
    for (cell in sample.int(J, J, replace = TRUE)) {
      if (stats::runif(1) < specRate / J) {
        counter <- counter + 1L
        mat[cell] <- counter
        next
      }
      for (k in 1:maxTries) {
        recruit <- mat[drawParent(cell)]
        if (density == 0 || k == maxTries) break
        conspecific <- mean(mat[shift(cell, neighbors[, 1], neighbors[, 2])] == recruit)
        if (stats::runif(1) > density * conspecific) break
      }
      mat[cell] <- recruit
    }
    if (i %in% runs) output[[length(output) + 1]] <- list(specMat = mat)
    if (verbose && i %% 100 == 0) message("Generation ", i)
  }

  runtime <- as.numeric((proc.time() - ptm)[3])
  if (verbose) message("Finished after ", round(runtime, 1), " seconds.")

  model <- list(x = x, y = y,
                dispersal = if (is.function(dispersal)) paste(deparse(dispersal), collapse = " ") else dispersal,
                runs = runs, specRate = specRate, density = density, densityCut = densityCut,
                environment = 0, seed = seed, type = "NeutralMod", runtime = runtime)
  out <- list(Output = output, Model = model)
  class(out) <- "PhyloSim"
  return(out)
}
