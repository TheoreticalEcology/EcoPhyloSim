test_that("summaryRedsaABC works if no weight is NA", {
  skip_if_not_installed("abctools")
  set.seed(1)
  p <- data.frame(a = runif(200), b = runif(200))
  s <- data.frame(s1 = p$a + rnorm(200, sd = 0.1), s2 = p$b + rnorm(200, sd = 0.1), s3 = rnorm(200))
  expect_equal(dim(summaryRedsaABC(p, s)), c(200, 2)) # #129
})

test_that("runPhylosimAnalysis saves only the summaries with summariesOnly = TRUE", {
  dir <- tempfile()
  on.exit(unlink(dir, recursive = TRUE))
  run <- function(...) capture.output(runPhylosimAnalysis(variableParams = list(specRate = function() 2),
                                                          nRuns = 2, outDir = dir, summariesOnly = TRUE,
                                                          nCores = 1, chunkSize = 1, ...))
  expect_error(run(fixedParams = list(x = 20, y = 20, runs = 100)), "calculateSummaries")
  skip_on_cran() # starts a cluster for each chunk
  run(fixedParams = list(x = 20, y = 20, runs = 100, calculateSummaries = TRUE)) # 2 chunks of 1 run, #130
  expect_setequal(list.files(dir), c("paramOverview.txt", "param_1-1", "param_2-2", "summary_1-1", "summary_2-2"))
  summaries <- get(load(file.path(dir, "summary_2-2")))
  expect_length(summaries, 1)
  expect_true("racAuc" %in% names(summaries[[1]]))
})
