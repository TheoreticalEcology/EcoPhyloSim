# #124: theta of es()

`esTheta.R` checks the theta used by `es()` since #124 and compares `es()` with neutral simulations; `esTheta.csv` holds its results (regenerate with `Rscript Code/Issues/124/esTheta.R` from the repository root, about 5 minutes).

## Theta

In `runSimulation()` with global dispersal, one generation has J = x · y death–birth events, and the parent is a random cell other than the dead one (#81). `Landscape::speciation()` turns on average `specRate` random individuals into new species at the end of each generation. In the Moran model of Vallade & Houchmandzadeh (2003, Phys. Rev. E 68, 061902, doi:10.1103/PhysRevE.68.061902), instead, each offspring founds a new species with probability ν, here ν = `specRate` / J: their eqns 2 and 3 draw the parent from the J − 1 other individuals, eqn 6 gives θ = (J − 1) ν / (1 − ν) = `specRate` (J − 1) / (J − `specRate`), and eqn 7 gives E(S) = Σ θ / (θ + i), i = 0, …, J − 1. The paper notes that this θ equals Hubbell's fundamental biodiversity number (2 J ν) within a factor of 2 for large J and small ν. For the speciation of `runSimulation()`, apart from its timing, the same coalescent argument gives θ = `specRate` (J − 1) / J, i.e. up to 2.7% fewer species than `es()` in the simulated settings (10 × 10 cells, `specRate` = 4). `es()` used θ = 2 · `specRate`; for 30 × 30 cells, it gave 1.72 to 1.74 times the expected richness (`specRate` = 1, 2, 4).

## Simulations

Mean richness of `runSimulation()` (global dispersal, 100 seeds, recorded once after 5 · J generations) against `es()`:

| cells | `specRate` | mean (SE) | `es()` | old `es()` | mean / `es()` |
|---|---|---|---|---|---|
| 100 | 1 | 5.72 (0.22) | 5.19 | 8.39 | 1.10 |
| 400 | 1 | 7.27 (0.23) | 6.57 | 11.14 | 1.11 |
| 900 | 1 | 7.86 (0.23) | 7.38 | 12.76 | 1.07 |
| 100 | 2 | 9.24 (0.26) | 8.45 | 13.53 | 1.09 |
| 400 | 2 | 12.43 (0.32) | 11.17 | 18.98 | 1.11 |
| 900 | 2 | 14.58 (0.36) | 12.77 | 22.20 | 1.14 |
| 100 | 4 | 15.86 (0.35) | 13.82 | 21.29 | 1.15 |
| 400 | 4 | 20.77 (0.39) | 19.09 | 31.96 | 1.09 |
| 900 | 4 | 24.54 (0.41) | 22.26 | 38.36 | 1.10 |

The simulations have 7 to 15% (0.5 to 2.3 species) more species than `es()`. The likely cause (comment on #124) is that `Landscape::speciation()` adds the Poisson(`specRate`) new species at the end of each generation, right before the species are recorded, whereas the Moran model of the theory speciates at single births. `NeutralMod()`, which speciates at single births, matches `es()` (test in `test-NeutralMod.R`).

## Decisions

- `es()` uses the exact θ of the Moran model, the model of the theory and of `NeutralMod()`, rather than θ = `specRate` or the θ of `runSimulation()` above; its expected richness differs from that with θ = `specRate` by at most 2.1% (10 × 10 cells, `specRate` = 4) in the simulated settings, and the exact θ also holds for small grids such as the 3 × 2 grid of the `NeutralMod()` test.
- `dimensions` accepts `c(x, y)`, because the model allows non-square grids; one number still means a square grid.
- `test-neutralTheory.R` accepts simulated richness from `es()` − 4 SE to 1.2 · `es()` (seeds 1 to 20 give 12.45, i.e. 1.11 · `es()`), i.e., it allows the offset of the end-of-generation speciation, but fails for the old θ (lower bound 15.9).
- The end-of-generation speciation in the C++ model is not changed here: speciating at single births would change the results of all simulations.
