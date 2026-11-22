# Checks that a change of the C++ model keeps its results: runs short simulations covering the main code paths with
# two builds and compares their complete output with identical() (bit-identical results for the same seed). This
# covers more settings and longer runs than the hash test in PhyloSim/tests/testthat/test-reproducibility.R.
# usage: Rscript identity.R run <library folder> <output.rds> [long]
#        Rscript identity.R compare <reference.rds> <new.rds>
# Default: 16 settings x 2 seeds, censuses at 150 and 300 generations; "long": 7 density-dependence settings,
# seed 3, censuses at 1000 and 2000 generations. Grid 50 x 50 and specRate = 2 unless a setting says otherwise.
a <- commandArgs(TRUE)
if (!length(a) || !a[1] %in% c("run", "compare")) stop("first argument must be \"run\" or \"compare\"")

if (a[1] == "compare") {
  ref <- readRDS(a[2]); new <- readRDS(a[3])
  same <- mapply(identical, ref, new[names(ref)])
  cat("identical:", sum(same), "of", length(same), "\n")
  if (any(!same)) print(names(same)[!same])
  quit(save = "no")
}

suppressMessages(library(EcoPhyloSim, lib.loc = a[2]))
MR <- c("mortality", "recruitment")
settings <- list(
  neutralLocal1 = list(dispersal = 1),
  neutralLocal3 = list(dispersal = 3),
  neutralGlobal = list(dispersal = "global"),
  env           = list(dispersal = 1, environment = 1),
  dd2           = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2),
  dd2gauss      = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, negativeDDKernel = "gaussian"),
  dd2species    = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, negativeDDKernel = "species"),
  envDD2        = list(dispersal = 1, environment = 1, negativeDD = 1, negativeDDCut = 2),
  posDD         = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, positiveDD = 0.1, positiveDDCut = 3),
  dd2recruit    = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, fitnessActsOn = MR),
  envRecruit2   = list(dispersal = 2, environment = 1, fitnessActsOn = "recruitment"),
  repro         = list(dispersal = 1, environment = 1, fitnessActsOn = "reproduction"),
  ddGlobal      = list(dispersal = "global", negativeDD = 1, negativeDDCut = 2),
  fission1      = list(dispersal = 1, fission = 1, x = 30, y = 30),
  protracted    = list(dispersal = 1, protracted = 5),
  redQueen      = list(dispersal = 1, redQueen = 0.5, redQueenStrength = 2),
  # only in the long runs
  envDDGlobal   = list(dispersal = "global", negativeDD = 1, negativeDDCut = 2, environment = 1),
  fissionDD     = list(dispersal = 1, fission = 1, negativeDD = 1, x = 30, y = 30),
  reproDD       = list(dispersal = 1, negativeDD = 1, environment = 1, fitnessActsOn = c("reproduction", "mortality"))
)
long <- length(a) >= 4 && a[4] == "long"
sel <- if (long) c("dd2", "dd2species", "posDD", "dd2recruit", "envDDGlobal", "fissionDD", "reproDD") else
  names(settings)[1:16]
out <- list()
for (n in sel) for (seed in if (long) 3 else 1:2) {
  par <- do.call(createCompletePar, modifyList(list(x = 50, y = 50, specRate = 2, seed = seed,
                                                    runs = if (long) c(1000, 2000) else c(150, 300)), settings[[n]]))
  sim <- suppressWarnings(runSimulation(par))
  out[[paste(n, seed)]] <- lapply(sim$Output, function(o) o[setdiff(names(o), "runtime")])
}
saveRDS(out, a[3])
