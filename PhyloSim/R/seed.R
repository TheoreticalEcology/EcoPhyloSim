# Sets R's random number generator to 'seed', with fixed kinds so that results do not depend on
# RNGkind() (#14), and restores the caller's random number state when the calling function exits.
localSeed <- function(seed, envir = parent.frame()) {
  oldSeed <- globalenv()$.Random.seed
  oldKind <- RNGkind()
  restore <- function() {
    if (is.null(oldSeed)) {
      suppressWarnings(RNGkind(oldKind[1], oldKind[2], oldKind[3]))
      rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", oldSeed, envir = globalenv())
  }
  do.call(on.exit, list(as.call(list(restore)), add = TRUE), envir = envir)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
}
