# Numeric regression test: results of small simulations are compared with the snapshots in
# _snaps/reproducibility.md (for intended changes, see "Tests and checks" in CONTRIBUTING.md).

local_edition(3)

snapshotSimulation <- function(...) {
  simu <- quietSimulation(x = 30, y = 30, runs = 200, seed = 1, ...)
  out <- simu$Output[[1]]
  set.seed(1) # random plots in sac()
  list(abundances = sort(as.vector(table(out$specMat)), decreasing = TRUE),
       nTips = ape::Ntip(out$phylogeny),
       sac = sac(simu, size = c(2, 5, 10, 20), plot = FALSE)$sr.Mean)
}

test_that("simulation results are numerically reproducible", {
  skip_on_cran()
  expect_snapshot_value(snapshotSimulation(specRate = 2), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1), style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, environment = 0.5, negativeDD = 0.5),
                        style = "deparse", tolerance = 1e-6)
})

test_that("results with only one fitness mechanism are numerically reproducible", {
  skip_on_cran()
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, environment = 1),
                        style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, negativeDD = 1, negativeDDCut = 2),
                        style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, negativeDD = 1, fitnessActsOn = "both"),
                        style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, redQueen = 0.1, redQueenStrength = 1),
                        style = "deparse", tolerance = 1e-6)
  expect_snapshot_value(snapshotSimulation(specRate = 2, dispersal = 1, positiveDD = 0.05),
                        style = "deparse", tolerance = 1e-6)
})

# MD5 hashes of the complete output of the C++ model, one per element (all censuses), written as bytes in a fixed
# byte order; the phylogeny and the summaries are computed from these elements in R
outputHash <- function(...) {
  simu <- suppressMessages(suppressWarnings(quietSimulation(x = 24, y = 16, runs = c(100, 200), seed = 1,
                                                            specRate = 2, ...)))
  sapply(c("specMat", "traitMat", "envMat", "compMat", "neutMat", "phyloTXT"), function(element) {
    f <- tempfile()
    on.exit(unlink(f))
    con <- file(f, "wb")
    for (o in simu$Output) {
      if (element == "phyloTXT") writeBin(charToRaw(o$phyloTXT), con)
      else writeBin(as.vector(o[[element]]), con, endian = "little")
    }
    close(con)
    unname(tools::md5sum(f))
  })
}

test_that("simulation results are identical to the stored ones", {
  skip_on_cran()
  # local and global dispersal, fitness on mortality, recruitment and reproduction, all density kernels (#207)
  MR <- c("mortality", "recruitment")
  hashes <- c(
    global = outputHash(),
    local2 = outputHash(dispersal = 2),
    env = outputHash(dispersal = 1, environment = 1),
    dd = outputHash(dispersal = 1, negativeDD = 1, negativeDDCut = 2),
    ddGaussian = outputHash(dispersal = 1, negativeDD = 1, negativeDDCut = 2, negativeDDKernel = "gaussian"),
    ddSpeciesResponse = outputHash(dispersal = 1, negativeDD = 1, negativeDDCut = 2, negativeDDKernel = "species",
                                   negativeDDResponse = function(S) 1 - S^2),
    posDD = outputHash(dispersal = 1, negativeDD = 1, negativeDDCut = 2, positiveDD = 0.1, positiveDDCut = 3),
    envDDRecruitment = outputHash(dispersal = 1, environment = 1, negativeDD = 1, fitnessActsOn = MR),
    envRecruitmentNarrow = outputHash(dispersal = 1, environment = 1, nicheWidth = 0.01, fitnessActsOn = "recruitment"),
    envDDReproduction = outputHash(dispersal = 1, environment = 1, negativeDD = 1,
                                   fitnessActsOn = c("reproduction", "mortality")),
    globalDD = outputHash(negativeDD = 1, negativeDDCut = 2),
    globalEnvRecruitment = outputHash(environment = 1, fitnessActsOn = "recruitment"),
    fissionDD = outputHash(dispersal = 1, fission = 1, negativeDD = 1, negativeDDKernel = "species"))
  expect_snapshot_value(hashes, style = "deparse")
})
