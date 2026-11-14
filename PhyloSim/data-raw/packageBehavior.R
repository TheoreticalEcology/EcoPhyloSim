# Simulations for the vignette section on emerging patterns (#136).
#
# Run from the repository root, with PhyloSim installed:
#   Rscript PhyloSim/data-raw/packageBehavior.R run pilot 2           # all blocks, 2 cores
#   Rscript PhyloSim/data-raw/packageBehavior.R run main 8 density    # selected blocks
#   Rscript PhyloSim/data-raw/packageBehavior.R run main 8 part=1/3   # a third of the tasks
#   Rscript PhyloSim/data-raw/packageBehavior.R combine
#
# Each simulation is an independent task with a fixed seed. Its summary is saved
# at once in PhyloSim/data-raw/cache/ (not in git); finished tasks are skipped,
# so an interrupted run resumes where it stopped, and parts can run on different
# machines. "combine" aggregates the cache into PhyloSim/data/packageBehavior.rda.

library(PhyloSim)

cacheDir <- "PhyloSim/data-raw/cache"
dataFile <- "PhyloSim/data/packageBehavior.rda"

## Settings ----------------------------------------------------------------

# generations: length of runs per block ("default" for all others); the main run
# lengths are at least twice the time to equilibrium found in the pilot run
modes <- list(
  pilot = list(reps = 5, generations = c(default = 50000, density = 100000)),
  main  = list(reps = 30, generations = c(default = 50000, density = 100000, neutral = 200000))
)

# Log-spaced censuses, plus a dense second half for the equilibrium check
censusTimes <- function(G) sort(unique(c(round(10^seq(log10(50), log10(G), length.out = 20)),
                                         round(seq(G / 2, G, length.out = 9)))))

base <- list(x = 50, dispersal = "1", specRate = 2, density = 0, densityCut = 1,
             environment = 0, fitnessBaseMortalityRatio = 10, protracted = 0,
             redQueen = 0, redQueenStrength = 0)

grid <- function(...) expand.grid(..., stringsAsFactors = FALSE)
design <- list(
  neutral     = grid(dispersal = "global", specRate = c(0.5, 1, 2, 4, 8), x = c(20, 50, 100)),
  dispersal   = grid(dispersal = c("global", "0.5", "1", "2", "3")),
  environment = grid(environment = c(0, 0.25, 0.5, 0.75, 1), dispersal = c("global", "1")),
  selection   = rbind(grid(fitnessBaseMortalityRatio = c(1, 2, 5, 10, 20, 50), environment = 1, density = 0),
                      grid(fitnessBaseMortalityRatio = c(1, 2, 5, 10, 20, 50), environment = 0, density = 1)),
  density     = grid(density = c(0, 0.25, 0.5, 0.75, 1, 1.2), densityCut = c(1, 2)),
  envDensity  = grid(environment = c(0, 0.5, 1), density = c(0, 0.5, 1)),
  speciation  = rbind(grid(protracted = c(0, 1, 2, 5), redQueen = 0, redQueenStrength = 0),
                      grid(protracted = 0, redQueen = c(0.01, 0.1), redQueenStrength = 100)),
  gridSize    = grid(x = c(20, 50, 100)),
  # neutral, density dependence, environment and both, for global and local dispersal
  assembly    = grid(environment = c(0, 1), density = c(0, 1), dispersal = c("global", "1"))
)

# One row per simulation: block, setting, replicate, all parameters and seed
taskTable <- function(mode) {
  m <- modes[[mode]]
  do.call(rbind, lapply(names(design), function(b) {
    d <- design[[b]]
    for (n in names(base)) if (is.null(d[[n]])) d[[n]] <- base[[n]]
    d$block <- b
    d$setting <- seq_len(nrow(d))
    # density > 1 is not allowed (#145); removed after numbering, so that the
    # seeds and file names of the other settings stay unchanged
    d <- d[d$density <= 1, ]
    d <- d[rep(seq_len(nrow(d)), each = m$reps), ]
    d$rep <- rep(seq_len(m$reps), length.out = nrow(d))
    d$generations <- if (b %in% names(m$generations)) m$generations[[b]] else m$generations[["default"]]
    d$seed <- (match(mode, names(modes)) * 100 + match(b, names(design))) * 1e5 +
      d$setting * 100 + d$rep
    d$file <- file.path(cacheDir, sprintf("%s_%s_%02d_%02d.rds", mode, b, d$setting, d$rep))
    d[c("block", "setting", "rep", names(base), "generations", "seed", "file")]
  }))
}

## One simulation ----------------------------------------------------------

# share of cells at distance r (along rows and columns, torus) with the same
# species; r = 1 gives the share of the 4 neighbors of the same species
distanceDecay <- function(m, r = unique(pmin(c(1, 2, 4, 8, 16, 32), floor(min(dim(m)) / 2)))) {
  shift <- function(m, dr, dc) m[(seq_len(nrow(m)) + dr - 1) %% nrow(m) + 1,
                                 (seq_len(ncol(m)) + dc - 1) %% ncol(m) + 1]
  data.frame(distance = r, sameSpecies = sapply(r, function(d)
    mean(c(m == shift(m, d, 0), m == shift(m, 0, d)))))
}

censusStats <- function(sim, i) {
  o <- sim$Output[[i]]
  ab <- as.vector(table(o$specMat))
  p <- ab / sum(ab)
  # extant species only (the phylogeny also contains extinct ancestors)
  phylo <- tryCatch(ape::keep.tip(o$phylogeny, paste0("s", unique(c(o$specMat)))),
                    error = function(e) NULL)
  enoughTips <- !is.null(phylo) && ape::Ntip(phylo) >= 3
  # species age from the species ID: IDs are counted up at each speciation, at
  # a roughly constant rate (largest ID / generation); NA if only one species
  # is left (then runSimulation() returns no phylogeny)
  G <- sim$Model$runs[i]
  id <- as.integer(names(table(o$specMat)))
  age <- if (inherits(o$phylogeny, "phylo"))
    G * (1 - id / max(as.integer(sub("s", "", o$phylogeny$tip.label)), id)) else NA
  data.frame(
    generation = G,
    richness = length(ab),
    shannon = -sum(p * log(p)),
    conspecificNeighbors = distanceDecay(o$specMat, 1)$sameSpecies,
    corTraitEnv = suppressWarnings(cor(c(o$traitMat), c(o$envMat))),
    sdCompetitionTrait = sd(c(o$compMat)),
    colless = if (enoughTips) tryCatch(collessImbalance(list(Output = list(list(phylogeny = phylo))),
                                                        norm = "yule"),
                                       error = function(e) NA) else NA,
    gamma = if (enoughTips) tryCatch(ape::gammaStat(phylo), error = function(e) NA) else NA,
    meanSpeciesAge = mean(age),
    corAgeAbundance = if (length(ab) >= 3) suppressWarnings(cor(age, ab, method = "spearman")) else NA)
}

# conspecific density dependence of mortality over the last consecutive
# generations: logistic regression of death on the number of conspecifics
# within radius 2, with a species effect (pooled), and per species with at
# least 50 observations
cndd <- function(sim, last) {
  sub <- sim
  sub$Output <- sim$Output[last]
  sub$Model$runs <- sim$Model$runs[last]
  d <- getMortalityData(sub, radius = 2)
  d <- d[d$interval == 1, ]
  pooled <- tryCatch(coef(glm(died ~ con + factor(species), binomial, data = d))[["con"]],
                     error = function(e) NA)
  bySpecies <- do.call(rbind, lapply(split(d, d$species), function(s) {
    if (nrow(s) < 50 || length(unique(s$con)) < 2) return(NULL)
    data.frame(abundance = nrow(s) / length(unique(d$census)), slope = tryCatch(
      coef(glm(died ~ con, binomial, data = s))[["con"]], error = function(e) NA))
  }))
  list(pooled = pooled, bySpecies = bySpecies)
}

runTask <- function(task) {
  # long runs need several GB (the phylogeny keeps extinct species); wait until
  # at least 30 GB of memory are free (Windows only, otherwise no check)
  freeGB <- function() tryCatch(as.numeric(system2("powershell", c("-NoProfile", "-Command",
    "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"), stdout = TRUE)) / 1024^2,
    error = function(e) NA)
  while (isTRUE(freeGB() < 30)) Sys.sleep(60)
  # marker for tasks that crash R (removed when the task finishes)
  running <- sub("\\.rds$", ".running", task$file)
  file.create(running)
  on.exit(unlink(running))
  set.seed(task$seed) # R-side randomness (multi2di(), sac() subplots)
  par <- createCompletePar(
    x = task$x, y = task$x,
    dispersal = if (task$dispersal == "global") "global" else as.numeric(task$dispersal),
    # plus five consecutive generations at the end for the mortality analysis
    runs = sort(unique(c(censusTimes(task$generations), task$generations - 5:1))),
    specRate = task$specRate,
    density = task$density, densityCut = task$densityCut, environment = task$environment,
    fitnessBaseMortalityRatio = task$fitnessBaseMortalityRatio, protracted = task$protracted,
    redQueen = task$redQueen, redQueenStrength = task$redQueenStrength, seed = task$seed)
  out <- tryCatch({
    capture.output(sim <- suppressWarnings(runSimulation(par)))
    last <- length(sim$Output)
    mort <- cndd(sim, last - 5:0)
    list(census = do.call(rbind, lapply(which(sim$Model$runs %in% censusTimes(task$generations)),
                                        censusStats, sim = sim)),
         decay = distanceDecay(sim$Output[[last]]$specMat),
         cndd = mort$pooled,
         cnddSpecies = mort$bySpecies,
         sac = sac(sim, size = unique(round(2^seq(0, log2(task$x), length.out = 8))),
                   rep = 20, plot = FALSE)[c("size", "sr.Mean")],
         rac = sort(as.vector(table(sim$Output[[last]]$specMat)), decreasing = TRUE),
         specMat = if (task$rep == 1) sim$Output[[last]]$specMat,
         runtime = sim$Model$runtime)
  }, error = function(e) list(error = conditionMessage(e)))
  for (attempt in 1:5) { # retry if the file system is briefly unavailable
    if (!inherits(try(saveRDS(c(list(task = task), out), task$file), silent = TRUE),
                  "try-error")) break
    Sys.sleep(30)
  }
  invisible(NULL)
}

## Combine the cache into the package data set ------------------------------

meanSdN <- function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), n = sum(!is.na(x)))

combineMode <- function(mode) {
  files <- taskTable(mode)$file
  files <- files[file.exists(files)]
  if (length(files) == 0) return(NULL)
  res <- lapply(files, readRDS)
  ok <- !vapply(res, function(r) "error" %in% names(r), logical(1))
  params <- do.call(rbind, lapply(res, `[[`, "task"))
  params$file <- NULL
  keys <- c("block", "setting")
  settings <- unique(params[!names(params) %in% c("rep", "seed")])

  # replicate mean, sd and n of each indicator per setting and census
  long <- do.call(rbind, lapply(res[ok], function(r) cbind(r$task[keys], r$census)))
  indicators <- setdiff(names(long), c(keys, "generation"))
  census <- do.call(rbind, lapply(indicators, function(v) {
    a <- aggregate(long[[v]], long[c(keys, "generation")], meanSdN)
    data.frame(a[c(keys, "generation")], indicator = v, a$x)
  }))

  sacLong <- do.call(rbind, lapply(res[ok], function(r) cbind(r$task[keys], r$sac)))
  a <- aggregate(sacLong$sr.Mean, sacLong[c(keys, "size")], meanSdN)
  sacAgg <- data.frame(a[c(keys, "size")], a$x)

  racLong <- do.call(rbind, lapply(res[ok], function(r)
    cbind(r$task[keys], rank = seq_along(r$rac), abundance = r$rac)))
  a <- aggregate(racLong$abundance, racLong[c(keys, "rank")], meanSdN)
  racAgg <- data.frame(a[c(keys, "rank")], a$x)

  decayLong <- do.call(rbind, lapply(res[ok], function(r)
    if (!is.null(r$decay)) cbind(r$task[keys], r$decay)))
  decayAgg <- if (!is.null(decayLong)) {
    a <- aggregate(decayLong$sameSpecies, decayLong[c(keys, "distance")], meanSdN)
    data.frame(a[c(keys, "distance")], a$x)
  }
  # pooled CNDD coefficient per setting, and per-species slopes of all runs
  cnddLong <- do.call(rbind, lapply(res[ok], function(r)
    if (!is.null(r$cndd)) cbind(r$task[keys], cndd = r$cndd)))
  cnddAgg <- if (!is.null(cnddLong)) {
    a <- aggregate(cnddLong$cndd, cnddLong[keys], meanSdN)
    data.frame(a[keys], a$x)
  }
  cnddSpecies <- do.call(rbind, lapply(res[ok], function(r)
    if (!is.null(r$cnddSpecies)) cbind(r$task[c(keys, "rep")], r$cnddSpecies)))

  list(settings = settings, census = census, sac = sacAgg, rac = racAgg,
       decay = decayAgg, cndd = cnddAgg, cnddSpecies = cnddSpecies,
       runtime = data.frame(params[ok, c(keys, "rep")],
                            runtime = vapply(res[ok], `[[`, numeric(1), "runtime")),
       errors = data.frame(params[!ok, c(keys, "rep")],
                           error = vapply(res[!ok], `[[`, character(1), "error")),
       specMat = setNames(lapply(res[ok & params$rep == 1], `[[`, "specMat"),
                          with(params[ok & params$rep == 1, ], paste(block, setting, sep = "_"))))
}

## Main ----------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) && args[1] == "run") {
  mode <- args[2]
  nCores <- as.integer(args[3])
  tasks <- taskTable(mode)
  rest <- args[-(1:3)]
  part <- grep("^part=", rest, value = TRUE)
  blocks <- setdiff(rest, part)
  if (length(blocks)) tasks <- tasks[tasks$block %in% blocks, ]
  if (length(part)) {
    p <- as.integer(strsplit(sub("part=", "", part), "/")[[1]])
    tasks <- tasks[seq_len(nrow(tasks)) %% p[2] == p[1] - 1, ]
  }
  tasks <- tasks[!file.exists(tasks$file), ]
  dir.create(cacheDir, showWarnings = FALSE, recursive = TRUE)
  message(nrow(tasks), " simulations to run on ", nCores, " cores")
  cl <- parallel::makeCluster(nCores)
  parallel::clusterExport(cl, c("censusTimes", "censusStats", "distanceDecay", "cndd"))
  parallel::clusterEvalQ(cl, library(PhyloSim))
  # longest tasks first, for better load balancing (rough cost: cells x dispersal
  # kernel size x competition neighborhood size)
  kernel <- ifelse(tasks$dispersal == "global", 0.1,
                   1 + (2 * suppressWarnings(as.numeric(tasks$dispersal)))^2)
  neighborhood <- ifelse(tasks$density > 0, 1 + tasks$densityCut^2, 1)
  tasks <- tasks[order(-tasks$generations * tasks$x^2 * kernel * neighborhood), ]
  parallel::parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))), runTask)
  parallel::stopCluster(cl)
} else if (length(args) && args[1] == "combine") {
  packageBehavior <- c(lapply(setNames(names(modes), names(modes)), combineMode),
                       list(info = list(date = Sys.time(),
                                        version = as.character(packageVersion("PhyloSim")),
                                        commit = system("git rev-parse --short HEAD", intern = TRUE),
                                        modes = modes, sessionInfo = sessionInfo())))
  save(packageBehavior, file = dataFile, compress = "xz")
}
