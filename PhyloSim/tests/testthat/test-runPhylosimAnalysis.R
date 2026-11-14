test_that("summaryRedsaABC works if no weight is NA", {
  skip_if_not_installed("abctools")
  set.seed(1)
  p <- data.frame(a = runif(200), b = runif(200))
  s <- data.frame(s1 = p$a + rnorm(200, sd = 0.1), s2 = p$b + rnorm(200, sd = 0.1), s3 = rnorm(200))
  expect_equal(dim(summaryRedsaABC(p, s)), c(200, 2)) # #129
})
