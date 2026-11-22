# Recomputes the tables of README.md from the files in results/.
# usage: Rscript analyse.R [results folder]   (default: ../results relative to this script)
a <- commandArgs(TRUE)
here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
res <- if (length(a)) a[1] else file.path(here, "..", "results")
rd <- function(f) { d <- read.csv2(file.path(res, f), dec = ".", stringsAsFactors = FALSE); d$t <- as.numeric(as.POSIXct(d$time)); d }

# Each run of a pair is matched to the run of the other build with the same parameters and seed that ended
# closest in time (within window seconds, as both runs of a pair end within seconds); runs of B from other comparisons
# have no such partner and are dropped. Returns one row per pair with the ratio cpu(B) / cpu(A).
pairs <- function(d, A, B, window = 30) {
  b <- d[d$tag == B, ]; ref <- d[d$tag == A, ]
  m <- t(sapply(seq_len(nrow(b)), function(i) {
    r <- ref[ref$parameters == b$parameters[i] & ref$seed == b$seed[i], ]
    if (!nrow(r)) return(c(NA, Inf))
    j <- which.min(abs(r$t - b$t[i])); c(r$cpu_s[j], abs(r$t[j] - b$t[i])) }))
  b$cpuA <- m[, 1]; b <- b[m[, 2] <= window, ]
  b$ratio <- b$cpu_s / b$cpuA
  b
}
summ <- function(p, by) do.call(rbind, lapply(split(p, p[[by]]), function(q) data.frame(
  parameters = q$parameters[1], n = nrow(q), ratio = round(median(q$ratio), 3), min = round(min(q$ratio), 3),
  max = round(max(q$ratio), 3), cpuA = round(median(q$cpuA), 2), speedup = round(1 / median(q$ratio), 1))))

cat("== A/A test (pairs-aa.csv): ratios\n")
aa <- rd("pairs-aa.csv"); print(round(sort(pairs(aa, "m802a", "m802b")$ratio), 3))

steps <- function(file, vers) {
  d <- rd(file)
  do.call(rbind, lapply(seq_len(length(vers) - 1), function(k) {
    s <- summ(pairs(d, vers[k], vers[k + 1]), "parameters"); s$step <- paste(vers[k], "->", vers[k + 1]); s }))
}
cat("\n== Version history (pairs-history.csv)\n")
h <- steps("pairs-history.csv", c("v0.3", "v0.3.1", "b00c3d2", "dd2331720", "rng7385ca7", "m802e558"))
print(h[, c("step", "parameters", "n", "ratio", "min", "max")], row.names = FALSE)
cat("\nProduct of the step medians, 0.3 -> master:\n"); print(round(tapply(h$ratio, h$parameters, prod), 3))
cat("\n== Bisection 0.3 -> 0.3.1 (pairs-bisect.csv)\n")
print(steps("pairs-bisect.csv", c("v0.3", "46a21df", "0919b9b", "v0.3.1"))[, c("step", "parameters", "n", "ratio", "min", "max")], row.names = FALSE)

cat("\n== Prototypes against master (prototypes.csv): speed-up = 1 / median ratio; cpuA = median CPU s of the paired master runs\n")
speedups <- function(file, tags, window = 30) {
  d <- rd(file)
  for (tg in tags) {
    s <- summ(pairs(d, "m802", tg, window), "parameters")
    cat("\n", tg, "\n"); print(s[, c("parameters", "n", "cpuA", "speedup", "min", "max")], row.names = FALSE)
  }
}
speedups("prototypes.csv", c("proto1", "proto2", "proto3"))
cat("\n== Implementation in the package against master (implementation.csv), after each commit of #207\n")
# window 600 s: in the global recruitment settings, master ends minutes after the other side of a pair
speedups("implementation.csv", c("i207-c2", "i207-c3", "i207-c4", "i207-c5", "i207-c6"), window = 600)

cat("\n== Binned parent proposals (recruitment-check.csv): means before (i207-c5) and after (i207-c6), Wilcoxon p\n")
rc <- read.csv2(file.path(res, "recruitment-check.csv"), dec = ".")
for (s in unique(rc$setting)) {
  q <- rc[rc$setting == s, ]
  cat("\n", s, "(", sum(q$tag == "i207-c6"), "seeds )\n")
  print(signif(t(sapply(c("richness", "shannon", "mismatch", "meanTrait", "sdTrait", "meanComp"), function(k) {
    a <- q[q$tag == "i207-c5", k]; b <- q[q$tag == "i207-c6", k]
    c(before = mean(a), after = mean(b), p = suppressWarnings(wilcox.test(a, b)$p.value)) })), 3))
}

cat("\n== Scaling of master (scaling.csv): ns of CPU time per death, median\n")
sc <- rd("scaling.csv"); print(round(tapply(1e9 * sc$cpu_s / sc$deaths, sc$parameters, median)))
