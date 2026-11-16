# Design of the simulations that document the behavior of the model (vignette
# chapter "Emerging patterns", vignette "Model history"): blocks of settings in which
# one or two parameters are varied around a reference setting, with replicates.
# Parameter names are those of EcoPhyloSim 0.3.1 (density, densityCut are the
# arguments negativeDD, negativeDDCut since #172); adapters.R maps them to the
# arguments of each version.

# replicates and generations per mode; generations per block ("default" for all others)
modes <- list(
  pilot = list(reps = 5, generations = c(default = 50000, density = 100000, densityDispersal = 100000)),
  main  = list(reps = 30, generations = c(default = 50000, density = 100000, densityDispersal = 100000,
                                          neutral = 200000)),
  test  = list(reps = 2, generations = c(default = 300)) # quick check of the scripts
)

# log-spaced censuses, plus a dense second half for the equilibrium check
censusTimes <- function(G) sort(unique(c(round(10^seq(log10(50), log10(G), length.out = 20)),
                                         round(seq(G / 2, G, length.out = 9)))))

# reference setting; y = NA means y = x
base <- list(x = 50, y = NA, dispersal = "1", specRate = 2, density = 0, densityCut = 1,
             environment = 0, nicheWidth = 0.03659906, fitnessActsOn = "mortality",
             fitnessBaseMortalityRatio = 10, fission = 0, protracted = 0, redQueen = 0,
             redQueenStrength = 0)

grid <- function(...) expand.grid(..., stringsAsFactors = FALSE)
design <- list(
  neutral     = grid(dispersal = "global", specRate = c(0.5, 1, 2, 4, 8), x = c(20, 50, 100)),
  dispersal   = grid(dispersal = c("global", "0.5", "1", "2", "3")),
  specRateLocal = grid(specRate = c(0.5, 1, 4, 8)),
  environment = grid(environment = c(0, 0.25, 0.5, 0.75, 1), dispersal = c("global", "1")),
  selection   = rbind(grid(fitnessBaseMortalityRatio = c(1, 2, 5, 10, 20, 50), environment = 1, density = 0),
                      grid(fitnessBaseMortalityRatio = c(1, 2, 5, 10, 20, 50), environment = 0, density = 1)),
  fitnessActsOn = rbind(grid(fitnessActsOn = c("mortality", "reproduction", "both", "recruitment"),
                             environment = 1, density = 0),
                        grid(fitnessActsOn = c("mortality", "reproduction", "both", "recruitment"),
                             environment = 0, density = 1)),
  density     = grid(density = c(0, 0.25, 0.5, 0.75, 1), densityCut = c(1, 2)),
  densityDispersal = grid(density = 1, dispersal = c("global", "0.5", "2", "3")),
  envDensity  = grid(environment = c(0, 0.5, 1), density = c(0, 0.5, 1)),
  assembly    = grid(environment = c(0, 1), density = c(0, 1), dispersal = c("global", "1")),
  nicheWidth  = grid(nicheWidth = c(0.01, 0.03659906, 0.1, 0.3), environment = 1,
                     dispersal = c("global", "1")),
  speciation  = rbind(grid(protracted = c(0, 1, 2, 5), redQueen = 0, redQueenStrength = 0),
                      grid(protracted = 0, redQueen = c(0.01, 0.1), redQueenStrength = 100)),
  fission     = grid(fission = c(1, 2)),
  # speciation mechanisms over speciation rates, to compare them at the same richness
  # (point mutation: block specRateLocal; specRate = 2: blocks fission and speciation)
  speciationRate = rbind(
    grid(specRate = c(0.5, 1, 4, 8), fission = c(1, 2), protracted = 0, redQueen = 0, redQueenStrength = 0),
    grid(specRate = c(0.5, 1, 4, 8), fission = 0, protracted = c(1, 2), redQueen = 0, redQueenStrength = 0),
    grid(specRate = c(0.5, 1, 4, 8), fission = 0, protracted = 0, redQueen = 0.01, redQueenStrength = 100)),
  gridSize    = grid(x = c(20, 50, 100)),
  # odd number of columns (#142) and a non-square grid, with environmental filtering
  gridShape   = data.frame(x = c(49, 50), y = c(49, 40), environment = 1)
)

# seed from a hash of mode, block, parameter values and replicate, so that seeds
# do not depend on the order of blocks or settings
hashSeed <- function(s) {
  h <- 0
  for (b in utf8ToInt(s)) h <- (h * 31 + b) %% 2147483647
  as.integer(h)
}

# one row per simulation: block, setting, replicate, all parameters, generations, seed, file
taskTable <- function(mode, cacheDir) {
  m <- modes[[mode]]
  do.call(rbind, lapply(names(design), function(b) {
    d <- design[[b]]
    for (n in names(base)) if (is.null(d[[n]])) d[[n]] <- base[[n]]
    d$y[is.na(d$y)] <- d$x[is.na(d$y)]
    d$block <- b
    d$setting <- seq_len(nrow(d))
    d <- d[rep(seq_len(nrow(d)), each = m$reps), ]
    d$rep <- rep(seq_len(m$reps), length.out = nrow(d))
    d$generations <- if (b %in% names(m$generations)) m$generations[[b]] else m$generations[["default"]]
    d$seed <- vapply(paste(mode, b, do.call(paste, c(d[names(base)], sep = "/")), d$rep, sep = "/"),
                     hashSeed, integer(1))
    d$file <- file.path(cacheDir, sprintf("%s_%s_%02d_%02d.rds", mode, b, d$setting, d$rep))
    d[c("block", "setting", "rep", names(base), "generations", "seed", "file")]
  }))
}
