# Neutral model with global dispersal against neutral theory (#125): at equilibrium, the
# number of species follows the Ewens sampling formula, calculated by es() (#124).

test_that("richness of the neutral model matches es()", {
  skip_on_cran()
  richness <- sapply(1:20, function(seed) {
    specRich(quietSimulation(x = 20, y = 20, dispersal = "global", specRate = 2, runs = 1000, seed = seed))
  })
  expected <- es(2, 20)
  se <- sd(richness) / sqrt(length(richness))
  # 7 to 15% more species than es(), probably because the species are recorded right after the new
  # species of a generation are added (Code/Issues/124); es() with theta = 2 * specRate gave 50%
  # more species than the simulations
  expect_gt(mean(richness), expected - 4 * se)
  expect_lt(mean(richness), 1.2 * expected)
})
