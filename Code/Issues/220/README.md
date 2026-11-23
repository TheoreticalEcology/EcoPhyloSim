# #220: trait evolution within species

Decisions for the parameters `traitEvolution` and `traitPull` (PR #223) that rest on simulations.

## Spread within species and the message on small density widths

Question: from which Gaussian width (`negativeDDWidth`, `positiveDDWidth`) on does almost no neighbor count as similar, now that the trait change at birth can be set?

The old fixed threshold 0.001 was a tenth of the spread within species at the fixed values. With half-width `w` and pull `a`, a trait follows `t' = (1 - a) t + a m + r`, `r ~ U(-w, w)`, so its stationary standard deviation around a fixed species mean `m` is `w / sqrt(3 a (2 - a))`. `withinSpeciesSd.R` compares this with simulations (SD within species, 3 seeds):

| w | a | predicted | simulated |
|---|---|---|---|
| 0.01 | 0.2 | 0.0096 | 0.0095-0.0100 |
| 0.05 | 0.5 | 0.033 | 0.033-0.035 |
| 0.05 | 0.1 | 0.066 | 0.064-0.069 |
| 0.1 | 1 | 0.058 | 0.057-0.060 |

Decision: for a Gaussian kernel, `checkPar()` gives the message below a tenth of this value for the competition trait (second value of `traitEvolution` and `traitPull`) (0.00096 at the defaults instead of 0.001). The value is capped at 1 / sqrt(12) (about 0.29, the SD of a uniform trait in [0, 1]), since with large `w` and small `a` the formula exceeds what traits in [0, 1] allow. With `a = 0` there is no stable spread and no message.

## Vignette section "Trait evolution and niche conservatism"

Question: which small simulation shows how `traitEvolution` and `traitSpeciationJump` affect niche conservatism and the phylogenetic structure of local communities, with a build time of a few minutes?

- Setting: environmental filtering at recruitment with global dispersal (`environment = 1`, `fitnessActsOn = "recruitment"`, 50 x 50, `specRate = 2`, 3000 generations, 3 seeds), so that local communities are sorted by filtering rather than by local speciation.
- Clustering measure: `picante::ses.mpd()` with `null.model = "taxa.labels"` on the unweighted community matrix of 30 plots of 5 x 5 cells, share of plots with p <= 0.05 (99 repetitions). `nullModel()` was not used: its abundance-weighted MPD includes pairs of conspecific individuals, so it mostly measures spatial aggregation of species rather than phylogenetic relatedness. No numbers were recorded for this comparison.
- Niche conservatism: Blomberg's K of the species means of the environmental trait (`picante::Kcalc()`).

The code is the chunk `nicheConservatism` in `PhyloSim/vignettes/ModelDescription.Rmd`; it runs when the vignette is built (about 4 minutes on one thread). Means of 3 seeds (K / share of clustered plots):

| traitEvolution | jump 0 | jump 0.3 | jump 1 |
|---|---|---|---|
| 0.01 | 1.39 / 0.60 | 0.012 / 0.20 | 0.0023 / 0.02 |
| 0.03 | 0.66 / 0.36 | 0.028 / 0.21 | 0.0015 / 0.00 |
| 0.1 | 0.014 / 0.18 | 0.011 / 0.03 | 0.0034 / 0.04 |

K varies strongly between seeds (0.07 to 1.07 at `traitEvolution = 0.03` without jump).
