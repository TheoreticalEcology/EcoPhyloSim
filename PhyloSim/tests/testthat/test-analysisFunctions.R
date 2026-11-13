# Analysis functions on hand-built species matrices (no simulation)

m <- matrix(c(1, 1, 1, 2,
              2, 3, 1, 1,
              4, 1, 2, 1), nrow = 3, byrow = TRUE) # abundances: species 1: 7, 2: 3, 3: 1, 4: 1
fake <- fakeSimulation(list(matrix(1, 3, 4), m))

test_that("specRich counts the species of the selected result", {
  expect_equal(specRich(fake), 4)
  expect_equal(specRich(fake, which.result = 1), 1)
})

test_that("rac ranks species by abundance", {
  pdf(NULL)
  on.exit(dev.off())
  out <- rac(fake)
  expect_equal(out$Rank, 1:4)
  expect_equal(out$Abundance, c(7, 3, 1, 1))
  expect_equal(as.character(out$Species[1:2]), c("1", "2"))
  expect_equal(rac(fake, which.result = 1)$Abundance, 12)
  expect_error(rac(fake, which.result = "all", plot = "bar"), "not possible")
})

test_that("sac counts every cell as a new species if all cells differ", {
  sq <- fakeSimulation(matrix(1:16, 4, 4)) # every cell a different species
  out <- sac(sq, size = 1:4, plot = FALSE)
  expect_equal(out$size, (1:4)^2)
  expect_equal(out$sr.Mean, (1:4)^2)
})

test_that("es gives one species for one individual and increases with specRate", {
  expect_equal(es(3, 1), 1) # a single individual is a single species
  expect_gt(es(2, 10), es(1, 10))
  skip("es() uses theta = 2 * specRate, see #124")
  expect_equal(es(2, 3), sum(2 / (2 + 0:8)))
})

test_that("createCompletePar rejects invalid parameters", {
  expect_error(createCompletePar(runs = c(100, 50)), "largest")
  expect_error(createCompletePar(x = 10, y = 20, airmat = matrix(0.5, 10, 10)), "do not match")
  expect_error(createCompletePar(airmat = 2), "between 0 and 1")
  expect_error(createCompletePar(environment = 2), "between 0 and 1")
  expect_error(createCompletePar(fitnessBaseMortalityRatio = 0.5), "greater than or equal to 1")
})

test_that("runSimulation rejects invalid dispersal", {
  expect_error(runSimulation(createCompletePar(x = 10, y = 10, dispersal = 0.2)), "too small")
  expect_error(runSimulation(createCompletePar(x = 10, y = 10, dispersal = "other")), "wrong dispersal")
})
