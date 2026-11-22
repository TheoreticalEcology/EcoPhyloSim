# Compares the runtime of two builds of the package (e.g., before and after a change of the C++ model) over settings
# that cover the main costs of the model (each run takes about 1-2 s with the version of November 2026).
# Each simulation runs in a fresh R process; per setting and repetition, both builds run with the same seed, in the
# order A, B in odd and B, A in even repetitions. Reports per setting the median CPU time (user) of each build and the
# median, minimum and maximum of the ratio B / A over the repetitions.
# CPU times vary with other load on the machine: with the same build on both sides, 2 repetitions on a busy machine
# gave ratios between 0.85 and 1.21 (November 2026). This shows large effects reliably (e.g., the 1.2-60x speed-ups
# of #207); for differences of a few percent, use many repetitions, or pin both builds to the two logical CPUs of one
# physical core (pair.ps1 at commit 5efcd39 of PR #210, Code/Performance/runtime/scripts/).
# usage: Rscript benchmark.R <library A> <library B> [repetitions = 3] [settings, comma-separated; default: all]
a <- commandArgs(TRUE)

MR <- c("mortality", "recruitment")
settings <- list( # added to x = 50, y = 50, specRate = 2
  local1       = list(dispersal = 1, runs = 10000),
  local3       = list(dispersal = 3, runs = 9000),
  global       = list(dispersal = "global", x = 100, y = 100, runs = 4000),
  env          = list(dispersal = 1, environment = 1, runs = 4000),
  dd2          = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, runs = 5000),
  dd3          = list(dispersal = 1, negativeDD = 1, negativeDDCut = 3, runs = 3000),
  dd2gauss     = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, negativeDDKernel = "gaussian", runs = 1000),
  dd2species   = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, negativeDDKernel = "species", runs = 2000),
  posDD        = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, positiveDD = 0.1, positiveDDCut = 3, runs = 2000),
  dd2recruit   = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, fitnessActsOn = MR, runs = 1000),
  repro        = list(dispersal = 1, environment = 1, fitnessActsOn = "reproduction", runs = 600),
  globalRecruit = list(dispersal = "global", environment = 1, fitnessActsOn = "recruitment", runs = 300)
)

if (length(a) && a[1] == "--run") { # one simulation, called by the main part below: --run <library> <setting> <seed>
  suppressMessages(library(EcoPhyloSim, lib.loc = a[2]))
  par <- do.call(createCompletePar, modifyList(list(x = 50, y = 50, specRate = 2, seed = as.integer(a[4])),
                                               settings[[a[3]]]))
  t <- system.time(suppressWarnings(runSimulation(par)))
  cat(t[["user.self"]], "\n")
  quit(save = "no")
}

if (length(a) < 2) stop("usage: Rscript benchmark.R <library A> <library B> [repetitions] [settings]")
reps <- if (length(a) >= 3) as.integer(a[3]) else 3
sel <- if (length(a) >= 4) strsplit(a[4], ",")[[1]] else names(settings)
stopifnot(all(sel %in% names(settings)))
script <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
rscript <- file.path(R.home("bin"), "Rscript")
cpu <- function(lib, setting, seed) {
  out <- suppressWarnings(as.numeric(tail(system2(rscript, c(shQuote(script), "--run", shQuote(lib), setting, seed),
                                                  stdout = TRUE), 1)))
  if (length(out) != 1 || is.na(out)) stop("run failed: ", lib, " ", setting)
  out
}

res <- do.call(rbind, lapply(sel, function(s) {
  t <- t(sapply(seq_len(reps), function(i) {
    if (i %% 2) { A <- cpu(a[1], s, i); B <- cpu(a[2], s, i) } else { B <- cpu(a[2], s, i); A <- cpu(a[1], s, i) }
    c(A = A, B = B)
  }))
  r <- t[, "B"] / t[, "A"]
  row <- data.frame(setting = s, cpuA = median(t[, "A"]), cpuB = median(t[, "B"]), ratio = median(r),
                    min = min(r), max = max(r))
  print(row, row.names = FALSE, digits = 3)
  row
}))
cat("\nratio = CPU time B / A (< 1: B is faster)\n")
print(res, row.names = FALSE, digits = 3)
