# Archive: model version of the MSc thesis by AndyEnvFr (2025)

Model code used in the MSc thesis "How empirical patterns of stabilizing density
dependence emerge from complex feedbacks between phylogenetic effects and spatial
dynamics" (2025). Not part of the package history and not maintained.

## Origin

Files of the fork AndyEnvFr/EcoPhyloSim, branch master, commit 53948a3, without the
compiled files (`*.o`, `*.dll`, `*.so`) committed there and with `apTreeshape` and
`PhyloMeasures` removed from `Imports` in `PhyloSim/DESCRIPTION` (both archived on CRAN;
only used by summary statistics and `summary()`); the fork branched off at a25d281,
the parent of this commit.

## Model changes relative to a25d281

- Negative and positive density dependence (`negativeDensity`, `positiveDensity`) with
  their own radius (`nDensityCut`, `pDensityCut`) and specificity (`nDDNicheWidth`,
  `pDDNicheWidth`). With local dispersal, the specificity is used as lambda (0 < lambda
  <= 20) in the kernel `lambda / (1 - exp(-lambda)) / 20 * exp(-20 * lambda * d)` of the
  difference `d` in the competition trait; with global dispersal, a Gaussian kernel with
  this width is used.
- Fitness `1 + environment * exp(-0.5 * ((temp - mean) / envNicheWidth)^2)
  - negativeDensity * N + positiveDensity * P`, where `N` and `P` are the mean kernel
  values over the neighbours within the radius.
- Mortality: a drawn individual dies if `U(0, 3)` exceeds its fitness (every
  `fitnessBaseMortalityRatio`-th draw dies regardless); a generation has `x * y / 3`
  deaths.
- Phylogeny calculation and output switched off (macro `PHYL_OFF` in `src/debug.h`,
  proposed in #107).

Parts were later integrated into the package: the mortality data functions as
`getMortalityData()` (#100) and the environmental niche width as `nicheWidth` (#103,
#138). The density-dependence changes were proposed in #104 and #108; the redesign of
fitness and density dependence is #172.

## Building

`R CMD INSTALL PhyloSim` (the package is called `PhyloSim`). Built with R 4.6.2 on
Windows.
