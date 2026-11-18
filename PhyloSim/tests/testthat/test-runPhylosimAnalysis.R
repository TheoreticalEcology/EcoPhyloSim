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
  expect_error(run(fixedParams = list(x = 20, y = 20, generations = 100)), "calculateSummaries")
  skip_on_cran() # starts a cluster for each chunk
  run(fixedParams = list(x = 20, y = 20, generations = 100, calculateSummaries = TRUE)) # 2 chunks of 1 run, #130
  expect_setequal(list.files(dir), c("paramOverview.txt", "param_1-1", "param_2-2", "summary_1-1", "summary_2-2"))
  summaries <- get(load(file.path(dir, "summary_2-2")))
  expect_length(summaries, 1)
  expect_true("racAuc" %in% names(summaries[[1]]))
})

test_that("runSimulationBatch(onError = \"NA\") keeps the other simulations if one fails (#158)", {
  p <- createCompletePar(x = 10, y = 10, generations = 10, seed = 1)
  bad <- p
  bad$x <- 0
  run <- function(...) {
    out <- NULL
    capture.output(out <- runSimulationBatch(list(p, bad, p), ...))
    out
  }
  expect_error(run(), "Must be >= 1")
  expect_warning(out <- run(onError = "NA"), "1 of 3 simulations failed")
  expect_length(out, 3)
  expect_s3_class(out[[1]], "PhyloSim")
  expect_true(is.na(out[[2]]))
  expect_match(attr(out[[2]], "error"), "Must be >= 1")
  skip_on_cran() # starts a cluster
  expect_warning(out <- run(parallel = 1, onError = "NA"), "1 of 3 simulations failed")
  expect_true(is.na(out[[2]]))
  expect_s3_class(out[[3]], "PhyloSim")
})

test_that("runPhylosimAnalysis and chunkwiseCleanup drop only the failed simulations (#158)", {
  dir <- tempfile()
  on.exit(unlink(dir, recursive = TRUE))
  # simulations 2 and 4 fail; chunk 4-4 has no successful simulation, which chunkwiseCleanup() must handle
  local_mocked_bindings(runSimulationBatch = function(pars, parallel, strip, onError = "stop") {
    if (onError != "NA") stop("simulation failed")
    out <- lapply(pars, function(p) list(specRate = p$specRate))
    if (length(pars) == 3) out[[2]] <- NA else out[[1]] <- NA
    out
  })
  i <- 0
  capture.output(runPhylosimAnalysis(variableParams = list(specRate = function() i <<- i + 1),
                                     fixedParams = list(x = 20, y = 20, generations = 100, calculateSummaries = TRUE),
                                     nRuns = 4, outDir = dir, summariesOnly = TRUE, nCores = 1, chunkSize = 3))
  expect_length(get(load(file.path(dir, "summary_1-3"))), 3)
  skip_on_cran() # starts a cluster
  expect_message(chunkwiseCleanup(dir, nCores = 1), "Removed 2 of 4 simulations")
  cleaned <- file.path(dir, "cleaned_up")
  expect_equal(sapply(get(load(file.path(cleaned, "summary_clean_1-2"))), `[[`, "specRate"), c(2, 4))
  expect_equal(sapply(get(load(file.path(cleaned, "param_clean_1-2"))), `[[`, "specRate"), c(2, 4))
})

test_that("chunkwiseCleanup gives a message with the number of removed simulations (#191)", {
  skip_on_cran() # starts a cluster
  dir <- tempfile()
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  summaryChunk <- list(list(gammaDiversity = 5, alphaDiversity = NA), list(gammaDiversity = 6, alphaDiversity = 2))
  paramChunk <- list(list(x = 10), list(x = 20))
  save(summaryChunk, file = file.path(dir, "summary_1-2"))
  save(paramChunk, file = file.path(dir, "param_1-2"))
  expect_message(chunkwiseCleanup(dir, nCores = 1), "Removed 1 of 2 simulations")
  expect_equal(get(load(file.path(dir, "cleaned_up", "param_clean_1-1"))), list(list(x = 20)))
})

test_that("cluster workers load EcoPhyloSim from the first library of the session (#190)", {
  skip_on_cran() # copies the installed package and starts a cluster
  pkg <- find.package("EcoPhyloSim")
  skip_if_not(file.exists(file.path(pkg, "Meta", "package.rds")), "EcoPhyloSim is not installed")
  lib <- tempfile("lib")
  dir.create(lib)
  oldLib <- .libPaths()
  on.exit({.libPaths(oldLib); unlink(lib, recursive = TRUE)})
  file.copy(pkg, lib, recursive = TRUE)
  .libPaths(c(lib, oldLib))
  cl <- startCluster(1)
  on.exit(parallel::stopCluster(cl), add = TRUE, after = FALSE)
  workerPath <- parallel::clusterEvalQ(cl, getNamespaceInfo(loadNamespace("EcoPhyloSim"), "path"))[[1]]
  expect_equal(normalizePath(workerPath), normalizePath(file.path(lib, "EcoPhyloSim")))
})

test_that("old parameter names replace the defaults of the new ones (#92)", {
  d <- getDefaultParams(formals(createCompletePar), list(runs = 5, density = 0.5))
  expect_false(any(c("generations", "negativeDD") %in% names(d)))
  expect_true("negativeDDCut" %in% names(d))
})
