#' @title Neutral model in R
#' @description A simple neutral model of community assembly (Hubbell 2001), written in R, for teaching and testing. It is a separate implementation that follows the rules of the C++ model of \code{\link{runSimulation}}, with the differences listed in Details, and returns an object of class \code{PhyloSim} with a species matrix only.
#' @param x Integer, dimension of the landscape in x-direction; the grid must have at least 2 cells (x * y >= 2)
#' @param y Integer, dimension of the landscape in y-direction
#' @param dispersal \code{"global"}, or a function with one argument \code{n} that returns \code{n} non-negative dispersal distances in cells (see Details)
#' @param runs Integer or increasing vector of integers, number of generations. If a vector, the species matrix is also saved at the intermediate generations, as in \code{\link{createCompletePar}}.
#' @param specRate Numeric, mean number of speciation events per generation (at most \code{x * y})
#' @param negativeDD Numeric between 0 and 1, strength of conspecific negative density dependence; 0 (default) switches it off.
#' @param negativeDDCut Integer, radius of the neighborhood for negative density dependence. If \code{negativeDD > 0}, at most \code{(min(x, y) - 1) / 2}, so that no cell is counted twice.
#' @param positiveDD,positiveDDCut Strength (between 0 and 1) and radius of conspecific positive density dependence, as \code{negativeDD} and \code{negativeDDCut}.
#' @param seed Integer, random seed. If NULL (default), a seed is drawn. Apart from this draw, the random number stream of the R session is not changed.
#' @param verbose Logical, if TRUE, reports progress every 100 generations and the runtime
#' @details The landscape is a torus of \code{x} times \code{y} cells, initially occupied by a single species. One generation consists of \code{x * y} events. In each event, a random cell dies and is replaced by a recruit. With probability \code{specRate / (x * y)}, the recruit is a new species; otherwise it is the offspring of a parent cell. As in \code{\link{runSimulation}}, the parent is never the dead cell itself: it is a random other cell (\code{dispersal = "global"}), or the cell at a distance drawn from \code{dispersal} in a random direction, rounded to the grid, where draws that land on the dead cell are repeated. With global dispersal, the equilibrium species richness follows the Ewens sampling formula with \eqn{\theta = specRate (x y - 1) / (x y - specRate)}{theta = specRate * (x * y - 1) / (x * y - specRate)} (from the coalescent of this model), which is approximately \code{specRate} if \code{specRate} is much smaller than \code{x * y}.
#'
#' With \code{negativeDD > 0} or \code{positiveDD > 0}, a recruit of species \code{s} is accepted with probability \code{(1 - negativeDD * S_n) * (1 - positiveDD * (1 - S_p))}, where \code{S_n} and \code{S_p} are the proportions of species \code{s} among the cells within distance \code{negativeDDCut} and \code{positiveDDCut} of the dead cell; otherwise a new parent is drawn. The parent is therefore drawn with a probability proportional to its dispersal probability times this fitness. This is the same rule as in \code{\link{runSimulation}} with \code{fitnessActsOn = "recruitment"} and the \code{"species"} kernels. After 1000 rejections, the parent is drawn directly from these probabilities with global dispersal; with a dispersal function, the last recruit is accepted. Note that \code{dispersal} is defined differently than in \code{\link{createCompletePar}}. Speciation also differs: \code{\link{runSimulation}} adds new species in a batch at the end of each generation.
#'
#' The model is slow and only suitable for small landscapes. The output contains the species matrix (\code{specMat}) and the parameters (\code{Model}), so that functions that only need the species matrix, such as \code{\link{specRich}}, \code{\link{rac}}, \code{\link{sac}} and \code{\link{localPlots}}, can be used. There is no phylogeny and there are no traits, so functions that need them (e.g. \code{plot()} and \code{summary()}) do not work.
#' @references Hubbell, Stephen P. The unified neutral theory of biodiversity and biogeography (MPB-32). Vol. 32. Princeton University Press, 2001.
#' @return An object of class \code{PhyloSim}, see \code{\link{runSimulation}}. Each element of \code{Output} contains only \code{specMat}, a matrix of species IDs with \code{y} rows and \code{x} columns.
#' @example /inst/examples/neutralModelR-help.R
#' @export
NeutralMod <- function(x = 50, y = 50, dispersal = "global", runs = 100, specRate = 1,
                       negativeDD = 0, negativeDDCut = 1, positiveDD = 0, positiveDDCut = 1, seed = NULL,
                       verbose = TRUE) {

  checkmate::assert_count(x, positive = TRUE)
  checkmate::assert_count(y, positive = TRUE)
  if (x * y < 2) stop("The grid must have at least 2 cells (x * y >= 2).")
  if (!is.function(dispersal) && !identical(dispersal, "global"))
    stop("dispersal must be \"global\" or a function(n) returning n distances")
  checkmate::assert_integerish(runs, lower = 1, any.missing = FALSE, min.len = 1, sorted = TRUE, unique = TRUE)
  checkmate::assert_number(specRate, lower = 0, upper = x * y)
  checkmate::assert_number(negativeDD, lower = 0, upper = 1)
  checkmate::assert_count(negativeDDCut, positive = TRUE)
  if (negativeDD > 0) checkmate::assert_int(negativeDDCut, upper = (min(x, y) - 1) %/% 2)
  checkmate::assert_number(positiveDD, lower = 0, upper = 1)
  checkmate::assert_count(positiveDDCut, positive = TRUE)
  if (positiveDD > 0) checkmate::assert_int(positiveDDCut, upper = (min(x, y) - 1) %/% 2)
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
  negNb <- as.matrix(getCircularOffsets(negativeDDCut))
  posNb <- as.matrix(getCircularOffsets(positiveDDCut))
  maxTries <- 1000 # as in the C++ model (Landscape::maxRecruitTries)

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

  # fitness of a recruit of species s in the dead cell, as in the C++ model with the species kernels
  recruitFitness <- function(s, cell) {
    w <- 1
    if (negativeDD > 0) w <- w * (1 - negativeDD * mean(mat[shift(cell, negNb[, 1], negNb[, 2])] == s))
    if (positiveDD > 0) w <- w * (1 - positiveDD * (1 - mean(mat[shift(cell, posNb[, 1], posNb[, 2])] == s)))
    w
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
        if (negativeDD == 0 && positiveDD == 0) break
        if (stats::runif(1) < recruitFitness(recruit, cell)) break
        if (k == maxTries && !is.function(dispersal)) { # exact draw, as in the C++ model (#165)
          prob <- vapply(seq_len(counter), recruitFitness, numeric(1), cell = cell)[mat] + .Machine$double.xmin * 100
          prob[cell] <- 0
          recruit <- mat[sample.int(J, 1, prob = prob)]
        }
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
                runs = runs, specRate = specRate, negativeDD = negativeDD, negativeDDCut = negativeDDCut,
                positiveDD = positiveDD, positiveDDCut = positiveDDCut,
                environment = 0, seed = seed, type = "NeutralMod", runtime = runtime)
  out <- list(Output = output, Model = model)
  class(out) <- "PhyloSim"
  return(out)
}
