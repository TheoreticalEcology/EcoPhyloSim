# Runtime of the simulation model: diagnosis and prototypes

Where the C++ model (`PhyloSim/src/`, master 802e558, November 2026) spends its time, how past changes affected it, and five changes that make it up to 60 times faster: three prototyped as patches in `prototypes/`, all five part of the package since #207 (see [Implementation](#implementation)); results are unchanged except with global dispersal and the environment acting on recruitment, where only the random numbers change. The benchmark tables can be recomputed from `results/` with `scripts/analyse.R`; profile shares and counters are taken from the summaries in `results/`. Issue: #207.

## Summary

- **Local dispersal**: for every death, the weights of all K_c cells of the dispersal kernel are recomputed with `exp()` and stored in newly allocated vectors, although they depend only on the offset and never change during a run.
- **Density dependence**: after every death, the densities of all K cells within the cutoff are recomputed, each over its K neighbours; only 1-4 of them are read before they change again.
- **Global dispersal with `fitnessActsOn = "recruitment"`**: a fallback that evaluates all N cells dominates when many parents are rejected; its cost per generation grows with N² (addressed by commit 5, see [Implementation](#implementation)).
- **History**: master is 12 % slower than 0.3 with local dispersal and 3 % faster with global dispersal; the local slowdown between 0.3 and 0.3.1 comes entirely from the #149 fix (0919b9b).
- **Prototypes**: a dispersal table computed once per run, an exact rewrite of the density loops, and densities recomputed only when read. Output identical to master in all 39 test simulations (52 for the third); 6.6x faster with `dispersal = 1`, 42x with `dispersal = 3`, 3.8-19x with density dependence.
- **Implementation** (#207): the three prototypes without their fallbacks to the old code, an inline copy of R's Mersenne-Twister, and binned parent proposals for global dispersal with recruitment; 1.2-60x faster in the benchmark settings, output identical to master except in the last setting (same distribution, other random numbers).

## Methods

**Builds.** All versions were installed with `R CMD INSTALL --debug` (the normal `-O2`, plus DWARF symbols, not stripped). Versions before #164 are package PhyloSim and use the argument names `density` and `densityCut`; versions before #113 need the archived dependencies of that time.

**Paired CPU-time benchmark** (`scripts/pair.ps1`, `scripts/run-one.R`). The machine was fully loaded by other simulation runs during all measurements; wall-clock times varied by +-20 %. Two builds therefore ran at the same time, pinned to the two logical CPUs of one physical core, at priority BelowNormal, with the sides swapped in every repetition; the CPU time of `runSimulation()` (`user.self`) was compared within each pair. In an A/A test (the same build on both sides, 6 pairs), the pair ratios were 0.999-1.023 (`results/pairs-aa.csv`). Absolute times are comparable only within a pair. All runs are 50x50 with `specRate = 2` unless stated otherwise.

**Profiles** (`scripts/profile.ps1`, `profile-run.R`, `profile-analyse.R`). Very Sleepy 0.91 sampled 20 s of the early phase of 14 settings (`scripts/settings.R`, seed 1). Addresses in the package DLL were resolved with `addr2line` from Rtools, including inlined functions; time in other DLLs is reported per DLL. Shares are relative to the time inside `callModel()`. Summaries: `results/profiles-802e558/`, `results/profiles-prototype3/`, and for 0.3 and 0.3.1 `results/profiles-history/` (`dispersal = 1`, 25 s from 10 s after the start).

**Counters** (`prototypes/counters.patch`, `scripts/count-run.R`). An instrumented build counted events, deaths, parent draws, recruitment fallbacks, density updates, and calls of `localDissimilarity()` and `recruitFitness()` over 500 generations (`results/counters.txt`).

**Identity of results** (`scripts/identity.R`, `identity-compare.R`). 16 settings x 2 seeds (censuses at 150 and 300 generations) covering local and global dispersal, environment, all density kernels, positive density dependence, `fitnessActsOn` = recruitment and reproduction, fission, protracted speciation and Red Queen; and 7 density-dependence settings over 2000 generations. All output elements (`specMat`, `traitMat`, `envMat`, `compMat`, `neutMat`, `phylogeny`, `phyloTXT`, `summaries`) were compared with `identical()`. In addition, 13 other combinations were checked for prototype 3 (30x30 and grids of 5x5 and 5x9, 300 generations), e.g. positive density dependence with a larger cutoff than the negative one, `dispersal = 0.5` and `0.75`, and global dispersal with recruitment.

## Where the time goes (master)

Shares of time in the profiles. "R.dll" is presumably R's `unif_rand()` (functions inside R.dll are not resolved); "heap" is ntdll + ucrtbase, presumably memory allocation; `exp` is `__expl_internal`.

| Setting | Main costs |
|---|---|
| global, neutral | R.dll 51 %, species bookkeeping (`reportDeath`, `updateMean`) 15 %, copying the individual 5 % |
| `dispersal = 1`, neutral | dispersal weights (`exp` 33 %, `kernel` 12 %), heap 15 %, R.dll 7 % |
| `dispersal = 3`, neutral | dispersal weights (`exp` 51 %, `kernel` 14 %) |
| `environment = 1` | dispersal weights (`exp` 33 %, `kernel` 11 %), heap 14 %, R.dll 9 % |
| `negativeDD = 1`, cutoff 1 / 2 / 3 | `localDissimilarity` 20 / 54 / 79 % |
| cutoff 2, kernel `"gaussian"` | `localDissimilarity` 86 % (of which `exp` 59 %) |
| cutoff 2, kernel `"species"` | `localDissimilarity` 56 % |
| `positiveDD = 0.1` (cutoff 3) and `negativeDD = 1` (cutoff 2) | `localDissimilarity` 80 % |
| `negativeDD = 1`, `fitnessActsOn = c("mortality", "recruitment")` | `localDissimilarity` 58 %; 20 % of the time is in `recruitFitness()`, mostly in its calls of `localDissimilarity` |

**Dispersal weights** (`Grid.cpp:458-492`, `Individual.cpp:103`): the kernel radius is `2 * dispersal` (`runSimulation.R:25`; truncated to an integer), so `dispersal` = 1, 2, 3, 5, 8 gives K_c = 12, 48, 112, 316, 796 kernel cells. The CPU time per death grows by about 72 ns per kernel cell (`results/scaling.csv`):

| `dispersal` | 1 | 2 | 3 | 5 | 8 |
|---|---|---|---|---|---|
| ns per death | 1342 | 4078 | 8736 | 23570 | 57525 |

**Density updates** (`Grid.cpp:531-584`): every death calls `densityUpdate()` once, which recomputes `localDissimilarity()` for the K cells within the cutoff (K = 5, 13, 29, 49, 113 for cutoff 1, 2, 3, 4, 6), each over K - 1 neighbours. The work per death therefore grows with K²; the measured extra CPU time over the neutral run grows about with K^1.65, as part of it is fixed overhead per call:

| `negativeDDCut` (`dispersal = 1`) | 1 | 2 | 3 | 4 | 6 |
|---|---|---|---|---|---|
| ns per death | 1866 | 3572 | 9344 | 20720 | 91267 |

At cutoff 6, density updates account for about 98 % of the time ((91267 - 1342) / 91267). A stored density is read only when its individual is tested for mortality (`getFitness()`, `deathProbability()`) or, with `fitnessActsOn = "reproduction"`, as a parent (`getSeedsTo()`): 1-4 times per death (the "events" below).

**Recruitment** (`Grid.cpp:332-356`, `497-514`): with `fitnessActsOn = "recruitment"`, parents are drawn until one is accepted; after min(K_c, 1000) rejections (local, `Grid.cpp:500`) or min(1000, N - 1) rejections (global, `Grid.cpp:335`), a fallback evaluates all K_c or N - 1 candidates.

Counters per death (`results/counters.txt`):

| Setting | Events | Parent draws | Fallbacks | `localDissimilarity` calls |
|---|---|---|---|---|
| `environment = 1` | 1.23 | 1 | - | - |
| `negativeDD = 1`, cutoff 1 / 2 | 1.02 / 1.01 | 1 | - | 5 / 13 |
| cutoff 2, `"species"` | 3.84 | 1 | - | 13 |
| `positiveDD` (cutoff 3) and `negativeDD` (cutoff 2) | 1.91 | 1 | - | 42 |
| `negativeDD`, cutoff 2, mortality and recruitment | 1.17 | 6.2 | 23 % | 22 |
| same, global dispersal | 1.26 | 8.2 | 0.01 % | 22 |
| `environment = 1`, recruitment only | 1 | 2.8 | 9 % | - |
| same, global dispersal | 1 | 69 | 5 % | - (192 `recruitFitness` calls) |

In the last setting, the fallback makes up 64 % of all fitness evaluations, at 50x50; with global dispersal its cost grows with N per death.

**Grid size** (`results/scaling.csv`): from 50x50 to 800x800, the CPU time per death grows from 1342 to 2509 ns with `dispersal = 1` and from 145 to 404 ns with global dispersal. This is probably due to CPU caches: the grid is allocated per column, and each `Individual` has 15 data members, a `RandomGen` object and a pointer to its species (not checked with hardware counters).

## Effect of past changes

Paired CPU time, ratio new / old, median of 6 pairs (`results/pairs-history.csv`, `results/pairs-bisect.csv`). Settings: local neutral (`dispersal = 1`, 4000 generations), local DD (`negativeDD = 1`, `negativeDDCut = 2`, 2000 generations), global neutral (100x100, 4000 generations).

| Step | Local neutral | Local DD | Global neutral |
|---|---|---|---|
| 0.3 (73ebaa3) -> 46a21df | 0.99 | 1.00 | |
| 46a21df -> 0919b9b (#149) | **1.12** | **1.07** | |
| 0919b9b -> 0.3.1 (269bb8d) | 1.00 | 1.00 | |
| 0.3 -> 0.3.1 | 1.12 | 1.07 | 1.01 |
| 0.3.1 -> b00c3d2 (before #172) | 1.00 | 1.00 | 1.01 |
| b00c3d2 -> 2331720 (#172) | 1.02 | 1.03 | 1.01 |
| 2331720 -> 7385ca7 (#175) | 0.98 | 1.01 | 0.95 |
| 7385ca7 -> 802e558 (master) | 1.00 | 1.00 | 1.00 |
| 0.3 -> master (product of the steps) | **1.12** | **1.12** | **0.97** |

All pair ratios were within +-3.3 % of the step median. 0919b9b (#149, merged with #161) replaced the variable-length arrays `weights` and `cumWeights` on the stack by `std::vector`; with local dispersal these are allocated for every death (in the non-neutral global branch once per generation). In profiles of `dispersal = 1`, ntdll grew from 3.9 % (0.3) to 9.8 % (0.3.1) of the time, and ntdll + ucrtbase from less than 4.9 % to 12.2 % (`results/profiles-history/`; these two profiles were taken 10 s after the start with an earlier version of `profile-run.R`).

In the model-behavior runs of the vignette (archives of 0.3 and 0.3.2 in `PhyloSim/data-raw/model-behavior/`, not in this folder; 14 parallel workers on 8 physical cores, with other load), 46 of 64 settings with local dispersal (without Red Queen or protracted speciation) took 8-20 % longer in 0.3.2 than in 0.3 (median of the per-setting medians 13 %, range -30 % to +23 %); this agrees with the measurements above within the noise of that setup. The only large difference, Red Queen with `specRate = 8` (median 5758 s in 0.3, 263 s in 0.3.2), comes from the phylogeny: at each of the 33 censuses, 0.3 copied all species that ever existed and built the Newick string with `insert(0, ...)`, which is quadratic in its length (20000 generations: 134 s with one census, 369 s with 33). This was fixed before 0.3.1 by a9aba4e (#140) and 852c61e (#160).

## Prototypes

The patches apply in this order to 802e558 (`patch -p1`). They are the prototypes as measured before the implementation, not finished code: the old loops are kept as fallback, and comments and tests are missing.

1. **`1-dispersal-table.patch`**: the offsets, weights and cumulative weights of the dispersal kernel are computed once per run, in the same order and with the same function as before; per death, a parent is drawn from the stored cumulative weights and only its coordinates are computed. The recruitment fallback uses a buffer allocated once. Used when fitness does not act on reproduction and the kernel radius is smaller than both grid dimensions (guaranteed by `checkPar.R:72`); otherwise the old code runs.
2. **`2-density-loops.patch`**: `densityUpdate()` and `localDissimilarity()` iterate over precomputed offset lists in the old order, instead of computing `floor(sqrt())` per row and two `%` per cell; the kernel type is a template parameter instead of a branch per neighbour; for the `"species"` kernel, species are compared by pointer instead of by ID (each species ID belongs to exactly one `Species` object).
3. **`3-lazy-density.patch`**: `densityUpdate()` only marks the cells whose density may have changed; a marked density is recomputed right before it is read. A stored density depends only on the current neighbourhood (the initial value 0 is exact, as all individuals start identical, and every change of a marker or species calls `densityUpdate()`), so the value read is the same as with immediate recomputation.

All three keep the same random draws and the same floating-point operations in the same order. Their output was identical to master in all identity checks (32 short and 7 long runs for each prototype, 13 additional runs for prototype 3; `results/identity.txt`).

Speed-up against master (1 / median ratio of 4-8 pairs, `results/prototypes.csv`; "master CPU s" is the median of the master runs paired with prototype 1):

| Setting | Master CPU s | 1 | 1+2 | 1+2+3 |
|---|---|---|---|---|
| `dispersal = 1`, 2000 generations | 6.99 | 6.4x | 6.4x | 6.6x |
| `dispersal = 3`, 300 generations | 6.62 | 38x | 40x | 42x |
| `environment = 1`, 1000 generations | 3.97 | 4.0x | 4.0x | 4.1x |
| `negativeDD = 1`, cutoff 2, 500 generations | 4.46 | 1.5x | 3.3x | 10x |
| `negativeDD = 1`, cutoff 3, 200 generations | 4.62 | 1.2x | 2.4x | 19x |
| cutoff 2, `"gaussian"`, 200 generations | 7.68 | 1.1x | 1.2x | 13x |
| cutoff 2, `"species"`, 500 generations | 5.50 | 1.3x | 2.8x | 5.5x |
| `positiveDD` (cutoff 3) and `negativeDD` (cutoff 2), 150 generations | 4.58 | 1.2x | 2.8x | 14x |
| `negativeDD`, cutoff 2, mortality and recruitment, 300 generations | 4.34 | 1.3x | 2.6x | 3.8x |

**After the prototypes** (`results/profiles-prototype3/`): with `dispersal = 1`, R.dll is the largest item (34 %), followed by copying the individual (10 %), the binary search for the parent (9 %) and the species bookkeeping in `reportDeath` (8 %). With density dependence (cutoff 2), recomputing the densities that are read (`ensureDensity()`, inclusive) takes 24 % (`"absolute"`), 34 % (`"species"`) and 68 % (`"gaussian"`, mostly `exp()`), marking them 3-9 %, and R.dll 8-24 %. With positive (cutoff 3) and negative density dependence, recomputing takes 50 % and marking 16 %. With selection at recruitment, `recruitFitness()` takes 34 %.

## Implementation

The three prototypes are implemented in the package, one commit each (#207), without the fallbacks to the old loops: `callModel()` (`PhylSimModelInterface.cpp`) and `checkPar()` guarantee `2 * cutoff + 1 <= min(x, y)` for the dispersal kernel and for active density dependence, so that no offset reaches the focal cell again and one addition or subtraction wraps a coordinate. The dispersal table is also used with `fitnessActsOn = "reproduction"`, with the weights of each death stored in buffers allocated once. After each commit, the output was identical to master in the identity checks (32 short and 7 long runs), in the reference test in `test-reproducibility.R` and in 54 further runs on small and non-square grids with the largest allowed cutoffs.

Two further commits (#207) address the remaining hot spots: (4) the random numbers come from an inline copy of R's Mersenne-Twister (`RandomGen.h`, `RandomGen.cpp`; same numbers, the state is taken from and written back to `.Random.seed`, `rpois()` still uses R's generator), and (5) with global dispersal and the environment acting on recruitment, parents are proposed from bins of their environmental trait after 10 rejected uniform proposals (`GlobalEnvironment::binnedParent()`). Without (5), early in a run, when few individuals fit the environment of a cell, many deaths needed 1000 rejected proposals and a fallback over all N cells, so that the cost per generation grew with N^2. (5) draws each parent with the same probability as before (~ fitness of its recruit), but with other random numbers, so that results for a given seed change in this setting only (one hash of the reference test). Checks of (5): the frequencies of drawn parents equal the exact probabilities in a test (`test-randomGen.R`, which detects a bound that is too small and a missing exclusion of the target cell); all other settings stay identical; and summary statistics agree (200 runs per build in one setting, 40 in two others; `scripts/recruitment-check.R`, `results/recruitment-check.csv`; e.g. species richness 20.6 and 20.7, all Wilcoxon p >= 0.11).

Speed-up against master after each commit (1 / median ratio of 4 pairs, `results/implementation.csv`, tags `i207-c2` to `i207-c6` for commits 1 to 5; "master CPU s" is the median of the master runs paired with the first commit in which the setting was measured; commit 2 was measured before a small cleanup outside the loops):

| Setting | Master CPU s | 1 | 1+2 | 1+2+3 | +4 | +5 |
|---|---|---|---|---|---|---|
| `dispersal = 1`, 2000 generations | 6.76 | 6.7x | 7.0x | 6.6x | 10x | 10x |
| `dispersal = 3`, 300 generations | 6.61 | 39x | 42x | 44x | 56x | 60x |
| `environment = 1`, 1000 generations | 3.87 | 4.1x | 4.1x | 4.0x | 5.4x | 5.2x |
| `negativeDD = 1`, cutoff 2, 500 generations | 4.40 | 1.5x | 3.6x | 11x | 12x | 12x |
| `negativeDD = 1`, cutoff 3, 200 generations | 4.74 | 1.2x | 2.6x | 21x | 27x | 25x |
| cutoff 2, `"gaussian"`, 200 generations | 7.64 | 1.1x | 1.2x | 13x | 13x | 13x |
| cutoff 2, `"species"`, 500 generations | 5.55 | 1.4x | 3.0x | 5.6x | 6.8x | 6.9x |
| `positiveDD` (cutoff 3) and `negativeDD` (cutoff 2), 150 generations | 4.74 | 1.2x | 2.8x | 16x | 17x | 18x |
| `negativeDD`, cutoff 2, mortality and recruitment, 300 generations | 4.20 | 1.3x | 2.8x | 4.1x | 4.8x | 4.8x |
| global, neutral, 100x100, 2000 generations | 3.56 | | | | 2.0x | 1.9x |
| global, `environment = 1`, recruitment, 300 generations | 15.39 | | | | | 4.5x |
| global, `environment = 1` and `negativeDD = 1`, mortality and recruitment, 200 generations | 78.10 | | | | | 6.9x |
| global, `environment = 1`, recruitment, 100x100, 30 generations | 107.92 | | | | | 5.8x |
| global, `environment = 0.5`, recruitment, 300 generations | 0.30 | | | | | 1.2x |

## Further options (not prototyped)

From three independent design reviews of the code and profiles; all would keep results unchanged unless noted.

- **Recruitment**: compute the neighbourhood of the target cell once per death instead of once per tried parent (the O(N) fallback with global dispersal is now mostly avoided by binned proposals, see [Implementation](#implementation)).
- **Incremental density sums**: update the stored sums of the K neighbours by the change of one cell, without recomputation. Exact for the `"species"` kernel (integer counts); for `"absolute"` and `"gaussian"` the summation order changes, so results change in the last bits and trajectories diverge. For `"gaussian"`, caching the K kernel values of each cell would save most `exp()` calls and keep results.
- **Memory layout for large grids**: one contiguous allocation for the grid; per-run constants (`m_nicheWidth`, `m_dispersalDistance`, `m_envStrength`, `m_compStrength`, `m_ddTable`, `m_posStrength`) moved out of `Individual`; species means computed when read instead of after every birth and death.
- **Global neutral runs**: one random index per cell instead of x and y changes results (an inline copy of R's Mersenne-Twister, which keeps them, is implemented, see [Implementation](#implementation)).
- Not worth it: alias tables or Gillespie-type sampling (changes results, little gain after prototype 1), threads within a run (the death-replacement chain is sequential), `-ffast-math` or `-march=native` (results would depend on the machine).

## Testing future changes

`tests/testthat/test-reproducibility.R` compares, besides summaries with a tolerance, an MD5 hash of the complete output of the C++ model for 13 short simulations (24x16, 200 generations, skipped on CRAN) with a snapshot; changes meant to keep results must pass it unchanged (the snapshot was made on Windows, x86-64). For more settings and longer runs, compare the full output of two builds with `scripts/identity.R` and `identity-compare.R`, and measure with `scripts/pair.ps1`. A change that is meant to change only the random numbers, not their distribution, needs a test of the distribution (as `test-randomGen.R` for the parents drawn) and a comparison of summary statistics over many seeds (as `scripts/recruitment-check.R`).

## Files

| Path | Content |
|---|---|
| `prototypes/` | the three prototypes and the counter instrumentation (patches against 802e558) |
| `scripts/` | `settings.R` (named settings), `run-one.R` and `pair.ps1` (paired benchmark; adjust `-cpus` and `-Rscript` to the machine), `identity.R` and `identity-compare.R`, `count-run.R`, `profile.ps1`, `profile-run.R` and `profile-analyse.R` (needs Very Sleepy and Rtools), `recruitment-check.R` (summary statistics for global recruitment), `analyse.R` (tables) |
| `results/` | `pairs-aa.csv`, `pairs-history.csv`, `pairs-bisect.csv`, `prototypes.csv`, `implementation.csv`, `scaling.csv` (one row per run: tag, parameters, seed, CPU, CPU s, elapsed s, deaths, time), `recruitment-check.csv` (one row per run: tag, setting, seed, statistics, CPU s), `counters.txt`, `identity.txt`, profile summaries |

The pair scripts of the history and bisection used named settings; `results/pairs-*.csv` give the equivalent `run-one.R` parameters. Raw profiles (`*.sleepy`), libraries and simulation output are not kept.
