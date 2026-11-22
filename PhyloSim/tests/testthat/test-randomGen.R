test_that("multinomialDraw() draws index i with probability w_i / sum(w)", {
  set.seed(1)
  w <- c(0, 1, 2, 0, 3, 0.5, 0)
  n <- 1e5
  counts <- multinomialDrawCounts(w, n)
  expect_equal(sum(counts), n)
  expect_equal(counts[w == 0], c(0L, 0L, 0L))
  expect_gt(stats::chisq.test(counts[w > 0], p = w[w > 0] / sum(w))$p.value, 0.001)
  expect_equal(multinomialDrawCounts(5, 10), 10L)
})

test_that("the C++ copy of R's Mersenne-Twister gives the numbers of runif() and continues its stream", {
  k <- 2^16 # equal weights: the index drawn is floor(k * u), so the counts check practically every number
  n <- 2000 # more than one block of 624 numbers, starting inside a block
  set.seed(3)
  runif(100)
  counts <- multinomialDrawCounts(rep(1, k), n)
  state <- .Random.seed
  set.seed(3)
  runif(100)
  expect_identical(counts, tabulate(floor(runif(n) * k) + 1L, k))
  expect_identical(state, .Random.seed)
})

test_that("runSimulation() leaves the random number state and kind of the R session unchanged", {
  oldKind <- RNGkind()
  on.exit(RNGkind(oldKind[1], oldKind[2], oldKind[3]))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(2)
  state <- .Random.seed
  a <- quietSimulation(x = 10, y = 10, runs = 20, seed = 1)
  expect_identical(.Random.seed, state)
  RNGkind("Mersenne-Twister")
  b <- quietSimulation(x = 10, y = 10, runs = 20, seed = 1)
  expect_identical(a$Output[[1]]$specMat, b$Output[[1]]$specMat) # independent of RNGkind()
})

# test of draw counts against probabilities ~ p: cells with expected counts below 5 are pooled (with the next
# smallest cells while the pool is still below 5); chi-square test, or exact binomial test if two groups remain
expectProportional <- function(counts, p) {
  expect_true(all(counts[p == 0] == 0))
  e <- sum(counts) * p / sum(p)
  o <- order(e)
  pool <- e[o] < 5
  while (sum(e[o][pool]) < 5 && sum(!pool) > 1) pool[which(!pool)[1]] <- TRUE
  obs <- c(sum(counts[o][pool]), counts[o][!pool])
  ex <- c(sum(e[o][pool]), e[o][!pool])
  pValue <- if (length(obs) == 2) stats::binom.test(obs[1], sum(obs), ex[1] / sum(ex))$p.value else
    stats::chisq.test(obs, p = ex / sum(ex))$p.value
  expect_gt(pValue, 1e-4)
}

test_that("binned proposals draw the parent of a recruit with probability ~ its fitness (#207)", {
  set.seed(4)
  x <- 12
  y <- 10
  n <- x * y
  env <- runif(n)
  draw <- function(traits, envStrength = 1, compStrength = 0, nicheWidth = 0.05, target = 17) {
    r <- binnedParentCounts(x, y, env, traits, runif(n), envStrength, compStrength, nicheWidth, target, 2e5)
    expect_equal(r$failed, 0)
    expectProportional(r$counts, r$fitness)
  }
  draw(runif(n))
  nBins <- ceiling(4 / 0.05)
  draw(pmin(1, pmax(0, sample(0:nBins, n, TRUE) / nBins + sample(c(-1e-12, 0, 1e-12), n, TRUE)))) # bin edges
  draw(runif(n), compStrength = 1)
  draw(runif(n), envStrength = 0.5, compStrength = 0.7)
  draw(rep(0.3, n), target = which.min(env) - 1) # all fit badly: fitness about 1e-8
  draw(runif(n), nicheWidth = 0.003) # 512 bins
  # fitness 0 for all cells: no draw, the fallback in GlobalEnvironment::reproduce() takes over
  expect_equal(binnedParentCounts(x, y, env, rep(0.5, n), runif(n), 1, 0, 0.05, which.min(env) - 1, 10)$failed, 10)
})
