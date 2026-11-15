# Maps a task of design.R to the arguments of createCompletePar() of a model
# version, or returns NULL if the version cannot run the task (a parameter or
# value it does not support). Versions: "0.2", "0.3" (package PhyloSim),
# "0.3.1" (EcoPhyloSim before #172), "current" (since #172).

adapt <- function(task, version) {
  supported <- list(
    "0.2"   = c("x", "y", "dispersal", "specRate", "density", "densityCut", "environment",
                "fitnessActsOn", "fitnessBaseMortalityRatio"),
    "0.3"   = c("x", "y", "dispersal", "specRate", "density", "densityCut", "environment",
                "fitnessActsOn", "fitnessBaseMortalityRatio", "fission", "protracted", "redQueen",
                "redQueenStrength"),
    "0.3.1" = c("x", "y", "dispersal", "specRate", "density", "densityCut", "environment",
                "nicheWidth", "fitnessActsOn", "fitnessBaseMortalityRatio", "fission", "protracted",
                "redQueen", "redQueenStrength"))
  supported$current <- supported[["0.3.1"]]
  stopifnot(version %in% names(supported))
  # a parameter the version lacks must have its default value in the task
  missing <- setdiff(names(base), supported[[version]])
  if (any(vapply(missing, function(n) !identical(task[[n]], base[[n]]), logical(1)))) return(NULL)
  if (task$fitnessActsOn == "recruitment" && version != "current") return(NULL)
  # fitness on reproduction is not implemented with global dispersal and fitness differences (#98, #101)
  if (task$fitnessActsOn %in% c("reproduction", "both") && task$dispersal == "global" &&
      (task$environment > 0 || task$density > 0 || task$redQueen > 0 || task$redQueenStrength > 0)) return(NULL)
  args <- as.list(task[intersect(supported[[version]], names(task))])
  args$dispersal <- if (task$dispersal == "global") "global" else as.numeric(task$dispersal)
  if (version == "current") {
    names(args)[names(args) == "density"] <- "negativeDD"
    names(args)[names(args) == "densityCut"] <- "negativeDDCut"
  }
  args$runs <- sort(unique(c(censusTimes(task$generations), task$generations - 5:1)))
  args$seed <- task$seed
  args
}
