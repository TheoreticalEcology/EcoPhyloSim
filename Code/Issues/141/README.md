# #141: dispersal between multiples of 0.5

## Problem

`runSimulation()` passed `2 * dispersal` to the integer argument `dispersalCutoff` of `callModel()`, so it was truncated: parents were drawn from the cells within the radius `R = floor(2 * dispersal)`, with weights `exp(-d / delta)` and `delta = R / 2`. Values between multiples of 0.5 therefore gave the same simulation as the next lower multiple.

## Decision

**Rule.** Parents are drawn from the cells within the real-valued distance `2 * dispersal` (cells with `dx^2 + dy^2 <= (2 * dispersal)^2`, without the dead cell), with weights `exp(-d / dispersal)`. `callModel(dispersalCutoff)`, the Landscape constructors and `circleOffsets()` take a `double`; `runSimulation()` is unchanged (it already passed `2 * dispersal`). The size check stays `2 * floor(2 * dispersal) + 1 <= min(x, y)`.

**Options considered.** Three independent proposals were collected; all three recommended this rule.

- `R = ceiling(2 * dispersal)`, `delta = dispersal` (the solution suggested in the issue, implemented first). Rejected: just above each multiple of 0.5 a whole ring of cells enters, so the cutoff moves between 2 and almost 4 times `delta`. Up to 44% of the weight then lies beyond `2 * delta`, and the mean dispersal distance jumps by up to 0.43 cells (e.g. 1.316 at `dispersal = 1.00`, 1.666 at 1.01). It also made the size check stricter (e.g. `dispersal = 2.4` on 10 x 10 cells was rejected).
- A smooth taper at the edge or a wider cutoff (e.g. `3 * delta`). Rejected: it changes the weights, and thus the results, at multiples of 0.5, needs larger grids and is harder to document.

The chosen rule cuts the kernel always at `2 * delta`, as before. Results at multiples of 0.5 are the same as before, because the cell sets and `delta` are identical. The neighbourhood still grows in steps, because lattice distances take only the values 1, √2, 2, √5, ...; the largest step in mean dispersal distance between 0.5 and 3 is 0.22 cells (master 0.58; ceiling rule 0.43), and it never decreases with `dispersal`. Below `dispersal = √2 / 2 ≈ 0.71`, the neighbourhood holds only the four nearest cells with equal weights, so all these values give the same model; this is documented in the vignette. (Numbers from the proposals' R computations of the discrete kernel on a 0.001 grid.)

**Rounding.** `circleOffsets()` and the size checks in `checkPar()` and `callModel()` add a tolerance of `1e-9`, so that inputs with rounding errors keep the cells at exactly `2 * dispersal`: `0.3 / 0.2` is 1.4999999999999998 and gives the same simulation as 1.5 (without the tolerance, it would lose the cells at distance 3). For integer radii (density dependence), the tolerance changes nothing.

**Other arguments.** The other integer arguments of `callModel()` (`negativeDDCut`, `positiveDDCut`, `fitnessBaseMortalityRatio`, `fission`, `protracted`) are already checked as whole numbers by `checkPar()`, and real-valued ones (`nicheWidth`, the kernel widths) are passed as `double`, so none has the same problem.

## Checks

- `identity.R`: with global dispersal and with `dispersal` at multiples of 0.5, results for a given seed are identical before and after the fix (14 settings with environment, density dependence, fitness on reproduction and recruitment, fission, positive density dependence, two or three censuses each; R 4.6.2, master 1c778ea against branch `i141-dispersal`). The snapshots in `test-reproducibility.R` are unchanged.
- `richness.R`: species richness against `dispersal` (neutral, 30 x 30 cells, `specRate = 2`, 1000 generations, 10 seeds). Mean (sd) of selected values:

| dispersal | before | after |
|---|---|---|
| 0.5 | 19.7 (5.3) | 19.7 (5.3) |
| 0.7 | 19.7 (5.3) | 19.7 (5.3) |
| 0.8 | 19.7 (5.3) | 19.1 (3.3) |
| 1.0 | 18.3 (4.1) | 18.3 (4.1) |
| 1.3 | 18.3 (4.1) | 14.8 (3.4) |
| 1.5 | 15.2 (2.1) | 15.2 (2.1) |
| 2.0 | 14.8 (2.7) | 14.8 (2.7) |
| 2.1 | 14.8 (2.7) | 13.9 (2.7) |
| 3.0 | 14.4 (3.8) | 14.4 (3.8) |

Before the fix, the values are constant between multiples of 0.5; after it, they change with every step of 0.1 from 0.8 on (Spearman correlation of mean richness with `dispersal` over 0.5 to 3: -0.89; with 10 seeds, the means are noisy).
