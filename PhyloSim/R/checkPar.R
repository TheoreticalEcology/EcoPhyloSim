# Checks a parameter list as created by createCompletePar(). Called by
# createCompletePar() and runSimulation(), so that lists changed by hand are
# also checked before they reach the C++ code (#118).
checkPar <- function(par) {
  # parameter lists created before negative and positive density dependence were separated
  for (old in c("density", "densityCut")) {
    new <- sub("density", "negativeDD", old)
    if (is.null(par[[new]])) par[[new]] <- par[[old]]
    else if (!is.null(par[[old]])) warning("'", old, "' is ignored, as '", new, "' is also given.", call. = FALSE)
    par[[old]] <- NULL
  }
  defaults <- list(negativeDDKernel = "absolute", negativeDDWidth = 0.1, negativeDDResponse = "linear",
                   positiveDD = 0, positiveDDCut = 1, positiveDDKernel = "species", positiveDDWidth = 0.1)
  for (n in names(defaults)) if (is.null(par[[n]])) par[[n]] <- defaults[[n]]
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
  checkmate::assert_number(par[["negativeDD"]], lower = 0, upper = 1, .var.name = "negativeDD", add = coll)
  checkmate::assert_number(par[["positiveDD"]], lower = 0, upper = 1, .var.name = "positiveDD", add = coll)
  checkmate::assert_number(par[["environment"]], lower = 0, upper = 1, .var.name = "environment", add = coll)
  checkmate::assert_number(par[["nicheWidth"]], lower = .Machine$double.eps, finite = TRUE,
                           .var.name = "nicheWidth", add = coll) # > 0
  checkmate::assert_subset(par[["fitnessActsOn"]], c("mortality", "reproduction", "recruitment", "both"),
                           empty.ok = FALSE, .var.name = "fitnessActsOn", add = coll)
  for (dd in c("negativeDD", "positiveDD")) {
    checkmate::assert_choice(par[[paste0(dd, "Kernel")]], c("absolute", "gaussian", "species"),
                             .var.name = paste0(dd, "Kernel"), add = coll)
    checkmate::assert_number(par[[paste0(dd, "Width")]], lower = .Machine$double.eps, finite = TRUE,
                             .var.name = paste0(dd, "Width"), add = coll) # > 0
    checkmate::assert_int(par[[paste0(dd, "Cut")]], lower = 1, .var.name = paste0(dd, "Cut"), add = coll)
  }
  if (!is.function(par[["negativeDDResponse"]]) && !identical(par[["negativeDDResponse"]], "linear")) {
    coll$push("Variable 'negativeDDResponse': Must be 'linear' or a function.")
  }
  checkmate::assert_int(par[["fitnessBaseMortalityRatio"]], lower = 1,
                        .var.name = "fitnessBaseMortalityRatio", add = coll)
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
  if (x * y < 2) stop("The grid must have at least 2 cells (x * y >= 2).", call. = FALSE)
  if (!length(par[["airmatR"]]) %in% c(1, x * y)) {
    stop("airmat must be a single value or have x * y = ", x * y, " values, not ",
         length(par[["airmatR"]]), ".", call. = FALSE)
  }
  if (par[["specRate"]] > x * y) {
    stop("specRate = ", par[["specRate"]], " must not exceed x * y = ", x * y, ".", call. = FALSE)
  }
  # fixes #152 (wider kernels wrap around the torus and count cells repeatedly)
  local <- is.numeric(par[["dispersal"]]) && par[["dispersal"]] > 0
  if (local && 2 * floor(2 * par[["dispersal"]]) + 1 > min(x, y)) {
    stop("dispersal = ", par[["dispersal"]], " is too large for a grid of ", x, " x ", y,
         ": 2 * floor(2 * dispersal) + 1 must not exceed min(x, y).", call. = FALSE)
  }
  for (dd in c("negativeDD", "positiveDD")) {
    cut <- par[[paste0(dd, "Cut")]]
    maxCut <- if (par[[dd]] > 0) (min(x, y) - 1) %/% 2 else min(x, y)
    if (cut > maxCut) {
      stop(dd, "Cut = ", cut, " is too large for a grid of ", x, " x ", y, ": ", dd,
           "Cut must not exceed min(x, y), and 2 * ", dd, "Cut + 1 must not exceed min(x, y) if ", dd,
           " > 0.", call. = FALSE)
    }
  }
  # fixes #98 (reproduction fitness is not implemented in the non-neutral global model)
  actsOn <- fitnessStages(par[["fitnessActsOn"]])
  if (!local && "reproduction" %in% actsOn &&
      (par[["environment"]] > 0 || par[["negativeDD"]] > 0 || par[["positiveDD"]] > 0 ||
       par[["redQueen"]] != 0 || par[["redQueenStrength"]] != 0)) {
    stop("Fitness acting on reproduction is not implemented for global dispersal ",
         "with environment, negativeDD, positiveDD, redQueen or redQueenStrength > 0. ",
         "Use local dispersal (dispersal >= 0.5) or fitnessActsOn = 'mortality'.", call. = FALSE)
  }
  # Red Queen changes the fitness on mortality or reproduction only, not of recruits
  if (identical(actsOn, "recruitment") && (par[["redQueen"]] != 0 || par[["redQueenStrength"]] != 0)) {
    stop("redQueen and redQueenStrength have no effect with fitnessActsOn = 'recruitment'; ",
         "add 'mortality' or 'reproduction' to fitnessActsOn.", call. = FALSE)
  }
  if (is.function(par[["negativeDDResponse"]])) {
    tab <- densityResponseTable(par[["negativeDDResponse"]], 101)
    if (any(diff(tab) > 1e-12)) {
      message("negativeDDResponse increases with the similarity to the neighbors, i.e. similar neighbors ",
              "increase fitness (use positiveDD for positive density dependence).")
    }
    if (par[["negativeDD"]] == 0) message("negativeDDResponse has no effect, because negativeDD = 0.")
  }
  # Valid, but probably not intended
  w <- par[["nicheWidth"]]
  if (par[["environment"]] > 0 && (w < 0.003 || w > 3)) {
    message("nicheWidth = ", w, " is outside [0.003, 3], so environmental filtering is very weak.")
  }
  for (dd in c("negativeDD", "positiveDD")) {
    width <- par[[paste0(dd, "Width")]]
    if (par[[dd]] > 0 && par[[paste0(dd, "Kernel")]] == "gaussian" && width > 1) {
      message(dd, "Width = ", width, " is larger than the range of the competition trait [0, 1], ",
              "so all neighbors count as similar.")
    }
  }
  same <- function(s) par[[paste0("negativeDD", s)]] == par[[paste0("positiveDD", s)]]
  if (par[["negativeDD"]] > 0 && par[["positiveDD"]] > 0 && same("Kernel") && same("Cut") &&
      (par[["negativeDDKernel"]] != "gaussian" || same("Width"))) {
    message("negativeDD and positiveDD use the same neighbors and kernel, so they partly offset each other; ",
            "use different kernels, widths or cutoffs to separate them.")
  }
  invisible(par)
}

# fitnessActsOn as a set of stages; "both" means mortality and reproduction
fitnessStages <- function(fitnessActsOn) {
  unique(unlist(lapply(fitnessActsOn, function(f) if (f == "both") c("mortality", "reproduction") else f)))
}

# Values of a user-defined negativeDDResponse at n equally spaced similarities in [0, 1]
densityResponseTable <- function(f, n) {
  tab <- vapply(seq(0, 1, length.out = n), function(S) as.numeric(f(S))[1], numeric(1))
  if (anyNA(tab) || any(tab < 0 | tab > 1)) {
    stop("negativeDDResponse must return values between 0 and 1 for similarities between 0 and 1.", call. = FALSE)
  }
  tab
}
