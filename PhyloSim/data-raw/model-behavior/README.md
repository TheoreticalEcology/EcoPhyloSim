# Model behavior across versions

Scripts that document how the behavior of EcoPhyloSim changes between versions
(#180, #181).

## Overview

- **Design** (`design.R`): 18 blocks of settings around a reference setting, in
  which a few parameters are varied; 134 settings with 30 replicates each,
  so 4,020 simulations per version, of 50,000 generations (100,000 for the
  blocks `density`, `densityDispersal`, `densityKernel` and `positiveDD`,
  200,000 for `neutral`). Parameter names are those
  of release 0.3.1 (e.g. `density`), and the current names for parameters added
  later. Each simulation has its own seed, which depends only on the mode, the
  block, the parameters that differ from `base` and the replicate, so results do
  not depend on the number of cores, and seeds do not change when a parameter is
  added to `base`. `design.R` also sets the censuses and `indicatorCommit`.
- **Versions** (`versions` in `adapters.R`, whose names are called adapters):
  the releases 0.2, 0.3 and 0.3.1, and `current`, any commit since #172
  (2331720); earlier commits run with `0.3.1`. The indicators of all runs are computed with
  one fixed installation of EcoPhyloSim (`indicatorCommit` in `design.R`), so
  that they differ between versions only through the model.
- **Stages** (all commands from the repository root):
  1. `run.R` installs the version and `indicatorCommit` if needed
     (`installVersion.R`, `installOldDeps.R`), runs each simulation in its own R
     process (`simulate.R`, with `adapters.R`) and computes its indicators
     (`indicators.R`), into `cache/`.
  2. `archive.R` turns the cache into an archive and writes its summary to
     `runs/` (see Running a version).
  3. `makeData.R` builds `modelHistory` (vignette "Model history") from `runs/`
     and `changes.csv`, and, if asked for, `packageBehavior` (vignette "EcoPhyloSim
     model") from the archive of the
     current version.
  4. The vignettes only read these data sets; they run none of these simulations.

  Separate from these: `attribution.R` tests single changes against their parent
  commit (evidence in `changes.csv`), and `pPD.R` creates the example data sets
  `pPD.pvalues` and `pPD.positions` of `plotPhylogeneticDispersion()` (one run
  per scenario, provenance in an attribute). `common.R` holds helpers of all scripts.
- **Runtime**: a run of mode `main` of one version takes days; mode `test` a few
  minutes. Check every edited script with a run of mode `test`.

| Output | Written by | In git |
|---|---|---|
| `changes.csv`: releases and the changes of the model code that change results (type N: results per seed change, but no change of their distributions was detected; B: patterns change; S: parameters renamed or redefined; S+B: both), with their effects, the evidence and the settings for `attribution.R` | by hand | yes |
| `runs/<version>_<commit>.csv`: summary of a run (mean, sd and n per setting and indicator at the last census, runtime, share of failed simulations) | `archive.R` | releases only (`*.9000_*.csv` ignored) |
| `PhyloSim/data/`: `modelHistory.rda`, `packageBehavior.rda`; `pPD.pvalues.rda`, `pPD.positions.rda` | `makeData.R`; `pPD.R` | yes |
| `archives/<version>_<commit>/`: per-replicate tables, copies of the scripts, provenance and raw outputs (`raw/`: all maps of the last six censuses and the last phylogeny of each simulation, enough to recompute the indicators at the end) | `archive.R` | no; release archives are attached to the GitHub release without `raw/` (about 10 MB, see Release), other runs stay local; `raw/` (about 0.9 GB per run) is kept only in the main checkout of the maintainer, without a backup |
| `cache/<mode>_<adapter>_<commit>/`: cache of a run | `run.R` | no; can be deleted once archived |
| `lib/<commit>/`, `lib/olddeps/`: libraries of the installed versions and of the archived dependencies | `installVersion.R`, `installOldDeps.R` | no |

`<adapter>` is the name in `versions` (e.g. `current`), `<version>` the release,
or for `current` the version in `PhyloSim/DESCRIPTION` (e.g. `0.3.1.9000`; at a
release, the release version). `archives/`,
`cache/` and `lib/` (ignored in `.gitignore`) are in the main checkout also when
a script runs in a git worktree, so that they are kept when the worktree is
removed; `git clean -x` deletes them. `runs/` and `PhyloSim/data/` are written
in the checkout the script runs in.

## Requirements

Git 2.31 or newer in the PATH, a compiler for R packages (Rtools on Windows),
and in the default library the packages that the old and the current versions
import (e.g. `Rcpp`, `ape`, `picante`, `treestats`). The scripts install the
model from the committed state of a commit (`git archive`), so uncommitted
changes of the model code are not simulated.

## When a pull request changes the model

This complements the pull request workflow in
[CONTRIBUTING.md](../../../CONTRIBUTING.md#pull-request-workflow).

In the pull request: a change that alters simulation results or the meaning of
arguments gets a row in `changes.csv`, below the last row of the release in preparation
(if that release has no row yet, add one with `kind` release and empty `date`
and `commit`):
`kind` change, `version` the next release, `reference` the issue, `type` (N, B,
S or S+B), the `change` and its expected effect on the `outputs`. New or renamed
arguments are added to `design.R` and `adapters.R` (see Extending). Contributors
do not need to run the simulations.

After the merge, a maintainer fills in `date` and `commit` (the commit on
`master`), and for the attribution `parent` (its parent commit), `adapter` and
`parentAdapter` (the names in `versions` of `adapters.R` that run both, usually
`current`) and, after stage 1, `stage2` (the blocks that differ); before the
release, the maintainer runs the attribution below and fills in `evidence`.

## Running a version

    Rscript PhyloSim/data-raw/model-behavior/run.R main <cores> current HEAD
    Rscript PhyloSim/data-raw/model-behavior/makeData.R

`run.R <mode> <cores> <adapter> <commit or tag>` runs the design with a commit
on `<cores>` worker processes. Interrupted runs resume when the same command is
run again, and failed simulations are run again. A run of mode `main` of all
blocks (no `part=`) without failed simulations is archived automatically.
Otherwise (e.g. v0.2, where some simulations always crash, see Older versions)
archive it with `Rscript PhyloSim/data-raw/model-behavior/archive.R <cache folder>`,
which keeps failed simulations as such; `archive.R <cache folder> partial` also
archives runs of some blocks. Runs of mode `test` are not archived
automatically; `archive.R` names their archive `test-<version>_<commit>`.

- Modes (`design.R`): `main` (30 replicates, the runs of record) and `test`
  (2 replicates, 300 generations; a check of the scripts, e.g.
  `run.R test 4 current HEAD gridShape fission`). Blocks after the commit restrict
  the run to these blocks, `part=i/n` to every n-th simulation.
- Run only one `run.R` per cache folder at a time; to split a run between
  processes, give each its own `part=i/n` (there is no locking), and archive the
  run with `archive.R <cache folder>` once all parts have finished.
- Mode `main` runs only with committed scripts, so that `RUNINFO.txt` names
  their commit, and a resumed run only with the scripts of its first start.
  Commit before a run; do not edit the scripts during a run.
- `MB_MIN_FREE_GB` (default 16): free memory in GB that `run.R` waits for
  before each simulation (Windows only).
- Progress of a run, without starting simulations: at the top level of its
  cache folder (`.running` files of a killed run remain), each
  finished simulation has an `.rds` file (failed ones contain an element
  `error`), each running one a `.running` file; `RUNINFO.txt` gives the number of
  simulations still to run at each start.

## Regenerating all runs

    Rscript PhyloSim/data-raw/model-behavior/run.R main <cores> 0.2 v0.2
    Rscript PhyloSim/data-raw/model-behavior/archive.R PhyloSim/data-raw/model-behavior/cache/main_0.2_2742f3a
    Rscript PhyloSim/data-raw/model-behavior/run.R main <cores> 0.3 v0.3
    Rscript PhyloSim/data-raw/model-behavior/archive.R PhyloSim/data-raw/model-behavior/cache/main_0.3_73ebaa3
    Rscript PhyloSim/data-raw/model-behavior/run.R main <cores> 0.3.1 v0.3.1
    Rscript PhyloSim/data-raw/model-behavior/run.R main <cores> current v0.3.2
    Rscript PhyloSim/data-raw/model-behavior/pPD.R <cores>
    Rscript PhyloSim/data-raw/model-behavior/makeData.R packageBehavior 0.3.2_<commit of v0.3.2>

Runs with failed simulations (0.2, 0.3) are not archived automatically, hence
the `archive.R` lines; for 0.3.2, `current` uses the adapter of the current
version until 0.3.2 is in `versions` of `adapters.R`. The run of record of
0.3.2 used commit dbdf3fc of #209 (scripts at 8c23d06, data sets built at
3905334), whose model code equals that of the tag v0.3.2 (archive
`0.3.2_dbdf3fc`); these commits are not on `master` (`git fetch origin
pull/209/head`). `archive.R` does not overwrite archives, so move existing ones
away first and delete the old files in `runs/`; run `pPD.R` with the release
commit checked out.

`makeData.R` reads every file in `runs/` and accepts one run per version, all
with `indicatorCommit` (see Release, step 3). `makeData.R packageBehavior
<archive>` also builds `packageBehavior`. `pPD.R` does not resume; it reruns failed scenarios once
with another seed, and if some still fail, it writes no data sets and keeps the
results of the others in `cache/pPD-partial.rds`.

## What to rerun after a change

| Changed | Rerun |
|---|---|
| model code (C++ or R code that `runSimulation()` uses) | none in the pull request; attribution and the run of the next release (see When a pull request changes the model, Release) |
| `indicators.R`, or a package function it uses (`collessImbalance()`, `getMortalityData()`, `nullModel()`, `sac()`; then set a new `indicatorCommit`) | all versions (`makeData.R` detects only a new `indicatorCommit`; the raw outputs in the archives would allow recomputing the final indicators, but no script does this yet) |
| `design.R` (blocks, settings, replicates) | all versions (a run uses the scripts of its first start, and archives are not overwritten) |
| `calculatePhylogeneticDispersion()` | `pPD.R` (also run at each release) |

## Release

1. Freeze the model code of the release (all changes in `changes.csv` with
   their attribution, see "When a pull request changes the model") and set the
   release version in `PhyloSim/DESCRIPTION`, after which the archive is named
   (e.g. `0.3.2_<commit>`); commit. The run is made with this commit; the tag
   follows in step 5, on a commit with the same model code.
2. Run the design with this commit (`run.R main <cores> current HEAD`), and
   `pPD.R` (`Rscript PhyloSim/data-raw/model-behavior/pPD.R <cores>`).
3. Delete the summaries of development runs from `runs/` (`makeData.R` reads all
   files there). Fill in `date` (release date) and `commit` (the commit of step
   1) of the release row in `changes.csv`, so that the vignette does not show
   the release as planned, and the same date in `Date` of
   `PhyloSim/DESCRIPTION`. Commit, then run `makeData.R packageBehavior <archive>`.
4. Update `notes` and `inShort` per release in `PhyloSim/vignettes/ModelHistory.Rmd` and the
   description of the runs in the section Simulation design of
   `PhyloSim/vignettes/ModelDescription.Rmd`.
5. Commit `runs/` and the data sets, tag the release, add the release to
   `versions` in `adapters.R` (a copy of `current`) and a row for the next
   release to `changes.csv`, and attach the archive to the GitHub release (in
   `archives/`: `tar -czf modelBehavior-<archive>.tar.gz --exclude=raw <archive>`).

## Attribution of single changes

    Rscript PhyloSim/data-raw/model-behavior/attribution.R run 1 6
    Rscript PhyloSim/data-raw/model-behavior/attribution.R compare 1
    Rscript PhyloSim/data-raw/model-behavior/attribution.R run 2 6
    Rscript PhyloSim/data-raw/model-behavior/attribution.R compare 2

For each change in `changes.csv` with a `parent`, stage 1 runs mode `test` for
the change and its parent with the same seeds and reports per block how many
simulations have no result (failed or unsupported) and in how many the census
indicators of the communities, species ages, tree statistics and the Colless
index differ. Stage 2 runs mode `main` for the blocks in `stage2` (if the change
itself crashes, `stage2Commit` gives another commit to compare with the parent)
and compares the last census with paired t-tests over the replicates (same seeds;
false discovery rate per change). Stage 2 can show a change of the
distributions, but not prove that there is none. `compare` writes
`cache/attribution-stage<stage>.csv`; enter the results in the column `evidence`
of `changes.csv`. Commits after the run can be added to the command to run only
these changes.

## Older versions

- `installVersion.R` applies build fixes only: v0.2 and v0.3 lack
  `src/RcppExports.cpp` (#102), which is generated with `Rcpp::compileAttributes()`;
  for v0.2 the lines with the unused `callLeipzig` are then removed from it.
- v0.2 imports `kulife`, v0.3 also `apTreeshape` and `PhyloMeasures`, which are
  archived on CRAN; `installOldDeps.R` installs `kulife` 0.1-14, `apTreeshape`
  1.5-0 with `cubature`, and `PhyloMeasures` 2.1. Their other dependencies (e.g.
  `ape`, `picante`, `quantreg`) must be in the default library.
- The DESCRIPTION of v0.2 gives the version number 0.3; the archive name comes
  from the version in `adapters.R`.
- v0.2 stops when only one species is left, as it cannot read a phylogeny with a
  single tip; `simulate.R` then sets the phylogeny to 0, as v0.3 does. v0.2 also
  crashes without a message (exit status 5) in some simulations, most often in
  the neutral model on 20 x 20 cells at high speciation rates (29 of 30 with
  `specRate` 8, 8 of 30 with 4), but also 1 to 4 of 30 in 16 other settings,
  e.g. with environmental filtering or density dependence; these count as failed (59 of
  2,730 supported simulations in the run of record).
- v0.3 stops with "numbers of left and right parentheses in Newick string not
  equal" in many simulations with the Red Queen advantage (all 30 with
  `redQueen` 0.1 and 59 of 150 with 0.01; 89 in the run of record), and it is run without fission, which assigns
  the wrong individuals in 0.3 (#85) and does not finish within hours
  (`adapters.R`); failed simulations count as failed.

## Extending

- New parameter: add it to `base` (with the default of the current version) and,
  if varied, to a block of `design`; add it to `params` of the versions that
  support it in `adapters.R`. Seeds of existing settings do not change (see
  Overview).
- Renamed argument: keep the design name and add it to `rename` of the version.
- New indicator: add it to `censusStats()` (at each census) or `behavior()` (at
  the end; the name of a single value also goes into the list with `mpdLocal`
  and `cndd` in `archive.R`, which adds them to the last census) in
  `indicators.R`, describe it in the help of `packageBehavior`
  (`PhyloSim/R/data-documentation.R`) and, to show it in the key figure, add it
  to `indicators` in `PhyloSim/vignettes/ModelHistory.Rmd`. Rerun: see What to
  rerun after a change.
- New block: add it to `design` and, if it needs another run length, to
  `generations` in `modes`; rerun as above. The vignettes select
  blocks by name and settings by parameter values, so renaming a block or
  changing a value requires updating them.

## Provenance

Each run records in its cache, and `archive.R` copies to the archive:

- `RUNINFO.txt`: start, mode, version, number of simulations and cores, the
  simulated package and commit, the installation of the indicators, the commit
  of the scripts (`git describe --dirty`) and the archive name;
- `SESSION-<start>.txt`: the command, the platform, the R session of the
  indicators, `INSTALLINFO.txt` of the libraries used, and all packages with
  their versions in these libraries;
- copies of the scripts at the start of the run.

Since #175 (R's random number generator), results for a given seed no longer
depend on the random number library of the compiler; for older versions they may
depend on platform and compiler.
