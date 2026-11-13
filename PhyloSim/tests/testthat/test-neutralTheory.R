# Neutral model with global dispersal against neutral theory (#125): at equilibrium, the
# number of species follows the Ewens sampling formula with theta = specRate (one generation
# has x * y death-birth events, specRate speciation events per generation; derivation in #124).
# Replace the formula by es() once #124 is fixed.

test_that("richness of the neutral model matches the Ewens sampling formula", {
  skip_on_cran()
  specRate <- 2
  J <- 20 * 20
  richness <- sapply(1:20, function(seed) {
    specRich(quietSimulation(x = 20, y = 20, dispersal = "global", specRate = specRate, runs = 3000, seed = seed))
  })
  i <- 0:(J - 1)
  expected <- sum(specRate / (specRate + i))
  se <- sqrt(sum(specRate * i / (specRate + i)^2) / length(richness))
  expect_lt(abs(mean(richness) - expected), 4 * se)
})
