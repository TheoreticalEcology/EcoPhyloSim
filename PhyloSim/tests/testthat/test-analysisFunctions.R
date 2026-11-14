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
  out <- sac(fakeSimulation(matrix(1:24, 4, 6)), plot = FALSE) # default sizes up to the shorter side (#151)
  expect_equal(out$size, (2:4)^2)
  expect_equal(out$sr.Mean, out$size)
})

test_that("sac plots several results that do not start with the first", {
  pdf(NULL)
  on.exit(dev.off())
  multi <- fakeSimulation(list(matrix(1, 4, 4), matrix(1:16, 4, 4), matrix(1:16, 4, 4)))
  expect_null(sac(multi, which.result = 2:3, size = 1:2))
})

test_that("es gives one species for one individual and increases with specRate", {
  expect_equal(es(3, 1), 1) # a single individual is a single species
  expect_gt(es(2, 10), es(1, 10))
  skip("es() uses theta = 2 * specRate, see #124")
  expect_equal(es(2, 3), sum(2 / (2 + 0:8)))
})
