# Checks a parameter list as created by createCompletePar(). Called by
# createCompletePar() and runSimulation(), so that lists changed by hand are
# also checked before they reach the C++ code (#118).
checkPar <- function(par) {
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_count(par[["x"]], positive = TRUE, .var.name = "x", add = coll)
  checkmate::assert_count(par[["y"]], positive = TRUE, .var.name = "y", add = coll)
  checkmate::assert_integerish(par[["runs"]], lower = 0, any.missing = FALSE, min.len = 1,
                               sorted = TRUE, .var.name = "runs", add = coll)
  d <- par[["dispersal"]]
  if (!(identical(d, "global") || (checkmate::test_number(d, finite = TRUE) && (d == 0 || d >= 0.5)))) {
    coll$push("Variable 'dispersal': Must be 'global', 0 or a number >= 0.5.")
  }
  checkmate::assert_number(par[["specRate"]], lower = 0, finite = TRUE, .var.name = "specRate", add = coll)
  checkmate::assert_number(par[["density"]], lower = 0, upper = 1, .var.name = "density", add = coll)
  checkmate::assert_number(par[["environment"]], lower = 0, upper = 1, .var.name = "environment", add = coll)
  checkmate::assert_number(par[["nicheWidth"]], lower = .Machine$double.eps, .var.name = "nicheWidth", add = coll) # > 0
  checkmate::assert_choice(par[["fitnessActsOn"]], c("mortality", "reproduction", "both"),
                           .var.name = "fitnessActsOn", add = coll)
  checkmate::assert_int(par[["fitnessBaseMortalityRatio"]], lower = 1,
                        .var.name = "fitnessBaseMortalityRatio", add = coll)
  checkmate::assert_int(par[["densityCut"]], lower = 1, .var.name = "densityCut", add = coll)
  checkmate::assert_int(par[["seed"]], .var.name = "seed", add = coll)
  checkmate::assert_choice(par[["fission"]], 0:2, .var.name = "fission", add = coll)
  checkmate::assert_count(par[["protracted"]], .var.name = "protracted", add = coll)
  checkmate::assert_number(par[["redQueen"]], lower = 0, finite = TRUE, .var.name = "redQueen", add = coll)
  checkmate::assert_number(par[["redQueenStrength"]], lower = 0, finite = TRUE,
                           .var.name = "redQueenStrength", add = coll)
  checkmate::assert_numeric(par[["airmatR"]], lower = 0, upper = 1, any.missing = FALSE, min.len = 1,
                            .var.name = "airmat", add = coll)
  checkmate::assert_string(par[["scenario"]], null.ok = TRUE, .var.name = "scenario", add = coll)
  checkmate::assert_flag(par[["calculateSummaries"]], .var.name = "calculateSummaries", add = coll)
  checkmate::assert_flag(par[["convertToBinaryTree"]], .var.name = "convertToBinaryTree", add = coll)
  checkmate::assert_flag(par[["prunePhylogeny"]], .var.name = "prunePhylogeny", add = coll)
  checkmate::reportAssertions(coll)

  # Conditions that connect several parameters
  x <- par[["x"]]; y <- par[["y"]]
  if (!length(par[["airmatR"]]) %in% c(1, x * y)) {
    stop("airmat must be a single value or have x * y = ", x * y, " values, not ",
         length(par[["airmatR"]]), ".", call. = FALSE)
  }
  local <- is.numeric(par[["dispersal"]]) && par[["dispersal"]] > 0
  if (local && 2 * par[["dispersal"]] > min(x, y)) {
    stop("dispersal = ", par[["dispersal"]], " is too large for a grid of ", x, " x ", y,
         ": 2 * dispersal must not exceed min(x, y).", call. = FALSE)
  }
  if (par[["densityCut"]] > min(x, y)) {
    stop("densityCut = ", par[["densityCut"]], " must not exceed min(x, y) = ", min(x, y), ".",
         call. = FALSE)
  }
  # fixes #98 (reproduction fitness is not implemented in the non-neutral global model)
  if (!local && par[["fitnessActsOn"]] != "mortality" &&
      (par[["environment"]] > 0 || par[["density"]] > 0 || par[["redQueen"]] != 0 ||
       par[["redQueenStrength"]] != 0)) {
    stop("fitnessActsOn = '", par[["fitnessActsOn"]], "' is not implemented for global dispersal ",
         "with environment, density, redQueen or redQueenStrength > 0. ",
         "Use local dispersal (dispersal >= 0.5) or fitnessActsOn = 'mortality'.", call. = FALSE)
  }
  invisible(par)
}
