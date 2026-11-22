# Regression test of the values returned by the analysis functions for one small simulation (#212).
# Only numbers are stored, without names or structure, so that interface changes (#199, #202)
# may adapt the calls and the extraction below (keeping the order of values), but not the snapshot.

local_edition(3)

# all numbers of an output, in order, without names; factors by their labels
values <- function(x) {
  if (is.list(x)) return(unlist(lapply(x, values), use.names = FALSE))
  if (is.factor(x)) x <- as.character(x)
  as.numeric(x)
}

test_that("outputs of the analysis functions are numerically unchanged", {
  skip_on_cran()
  # two recorded generations, each selected explicitly; local dispersal, because PR #210 changes
  # results with global dispersal and the environment acting on recruitment
  simu <- quietSimulation(x = 30, y = 30, runs = c(100, 200), seed = 1, specRate = 2, dispersal = 1,
                          environment = 0.5, negativeDD = 0.5)
  traits <- getAverageTraits(simu, which.result = 2, orderPhylo = TRUE) # orderPhylo = FALSE changes with #200
  seeded <- function(expr) {set.seed(1); values(expr)}
  out <- list(
    rac = values(rac(simu, which.result = 2, plot = FALSE)),
    racFirst = values(rac(simu, which.result = 1, plot = FALSE)),
    summaryStatistics = seeded(calculateSummaryStatistics(simu, which.result = 2)),
    colless = collessImbalance(simu, which.result = 2), # the other numbers of summary() are in summaryStatistics
    averageTraits = c(as.numeric(rownames(traits)), values(traits)),
    speciesTraits = values(getSpeciesTraits(simu, speciesID = 360, whichRun = 2)), # 13 individuals
    sac = seeded(sac(simu, which.result = 2, size = c(2, 5, 10), rep = 10, plot = FALSE)),
    nullModel = seeded(nullModel(simu, which.result = 2, localPlotSize = 5, numberOfPlots = 10,
                                 repetitions = 50)),
    nullModelPd = seeded(nullModel(simu, which.result = 2, localPlotSize = 5, numberOfPlots = 10,
                                   repetitions = 50, abundance = TRUE, fun = "pd")),
    phylogeneticDispersion = values(lapply(c("PhylPool", "SamplePool"), function(type)
      seeded(calculatePhylogeneticDispersion(simu, which.result = 2, plotlengths = 5, plots = 10,
                                             replicates = 50, type = type)))),
    es = es(specRate = 2, dimensions = 30)
  )
  expect_snapshot_value(out, style = "deparse", tolerance = 1e-6)
})
