# Model versions that can be run, and the mapping of a task of design.R to the
# arguments of createCompletePar() of a version. adapt() returns NULL if the
# version cannot run the task (a parameter or value it does not support).
# To add a release, see "Release" in README.md.

versions <- list(
  # oldDeps: needs the archived dependencies of installOldDeps.R
  "0.2"   = list(oldDeps = TRUE,
                 params = c("x", "y", "dispersal", "specRate", "density", "densityCut", "environment",
                            "fitnessActsOn", "fitnessBaseMortalityRatio")),
  # without fission: it assigns the wrong individuals in 0.3 (#85), so that the Newick string
  # breaks or runs of 50,000 generations take hours instead of minutes
  "0.3"   = list(oldDeps = TRUE,
                 params = c("x", "y", "dispersal", "specRate", "density", "densityCut", "environment",
                            "fitnessActsOn", "fitnessBaseMortalityRatio", "protracted",
                            "redQueen", "redQueenStrength")),
  "0.3.1" = list(
                 params = c("x", "y", "dispersal", "specRate", "density", "densityCut", "environment",
                            "nicheWidth", "fitnessActsOn", "fitnessBaseMortalityRatio", "fission",
                            "protracted", "redQueen", "redQueenStrength")),
  # since #172: density, densityCut are the arguments negativeDD, negativeDDCut;
  # fitness can act on recruitment; kernels and positive density dependence
  "0.3.2" = list(
                 params = c("x", "y", "dispersal", "specRate", "density", "densityCut", "environment",
                            "nicheWidth", "fitnessActsOn", "fitnessBaseMortalityRatio", "fission",
                            "protracted", "redQueen", "redQueenStrength", "negativeDDKernel",
                            "positiveDD", "positiveDDCut"),
                 rename = c(density = "negativeDD", densityCut = "negativeDDCut"),
                 fitnessActsOn = c("mortality", "reproduction", "both", "recruitment")),
  # since #114: trait jump at speciation
  current = list(
                 params = c("x", "y", "dispersal", "specRate", "density", "densityCut", "environment",
                            "nicheWidth", "fitnessActsOn", "fitnessBaseMortalityRatio", "fission",
                            "protracted", "redQueen", "redQueenStrength", "negativeDDKernel",
                            "positiveDD", "positiveDDCut", "speciationTraitJump"),
                 rename = c(density = "negativeDD", densityCut = "negativeDDCut"),
                 fitnessActsOn = c("mortality", "reproduction", "both", "recruitment"))
)

adapt <- function(task, version) {
  v <- versions[[version]]
  stopifnot(!is.null(v))
  # a parameter the version lacks must have its default value in the task
  missing <- setdiff(names(base), v$params)
  if (any(vapply(missing, function(n) !identical(task[[n]], base[[n]]), logical(1)))) return(NULL)
  if (!task$fitnessActsOn %in% c(v$fitnessActsOn, "mortality", "reproduction", "both")) return(NULL)
  # fitness on reproduction is not implemented with global dispersal and fitness differences (#98, #101)
  if (task$fitnessActsOn %in% c("reproduction", "both") && task$dispersal == "global" &&
      (task$environment > 0 || task$density > 0 || task$redQueen > 0 || task$redQueenStrength > 0)) return(NULL)
  args <- as.list(task[intersect(v$params, names(task))])
  args$dispersal <- if (task$dispersal == "global") "global" else as.numeric(task$dispersal)
  for (n in names(v$rename)) names(args)[names(args) == n] <- v$rename[[n]]
  args$runs <- sort(unique(c(censusTimes(task$generations), task$generations - 5:1)))
  args$seed <- task$seed
  args
}
