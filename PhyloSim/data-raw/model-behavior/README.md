# Model behavior across versions

Scripts that document how the behavior of EcoPhyloSim changes between versions
(#180, #181). Every version runs the same simulation design (`design.R`), and the
indicators of all runs are computed with one fixed installation of EcoPhyloSim
(`indicatorCommit` in `design.R`). The runs give two data sets of the package:
`packageBehavior` (chapter "Emerging patterns" of the main vignette, from the run
of the current version) and `modelHistory` (vignette "Model history", from the
runs of all versions, with the changes in `changes.csv`). Run all commands from
the repository root; after editing a script, check it with a run of mode `test`
(a few minutes, see below).

## Files

- `design.R`: blocks of settings, replicates, run lengths, censuses, seeds and
  `indicatorCommit`. Parameter names are those of release 0.3.1 (e.g. `density`),
  and the current names for parameters added later.
- `indicators.R`: all indicators of one simulation (`behavior()`).
- `adapters.R`: the model versions (`versions`: package name, supported
  parameters, renamed arguments, archived dependencies) and `adapt()`, which maps
  a task to the arguments of `createCompletePar()` of a version.
- `run.R`, `simulate.R`: run the design with one commit, each simulation in its
  own R process; `run.R` installs the commit and the indicators if needed and
  archives a complete run.
- `archive.R`: turns the cache of a run into its archive and writes its summary
  to `runs/`.
- `makeData.R`: builds `modelHistory` from `runs/` and `changes.csv`, and, if
  asked for, `packageBehavior` from the archive of the current version.
- `attribution.R`: tests single changes against their parent commit.
- `installVersion.R`, `installOldDeps.R`: install the package of a commit or tag
  into its own library, and the archived dependencies of the old versions (called
  by `run.R`).
- `common.R`: helpers of the scripts.
- `pPD.R`: the example data sets `pPD.pvalues` and `pPD.positions` of
  `plotPhylogeneticDispersion()` (one run per scenario; the attribute
  `provenance` of `pPD.pvalues` records the run).
- `changes.csv`: releases and the changes of the model code that change results
  (type N: results per seed change, but no change of their distributions was
  detected; B: patterns change; S: parameters renamed or redefined; S+B: both),
  with their effects, the evidence and the settings for `attribution.R`.
- `runs/<version>_<commit>.csv`: summary of the run of a version (mean, sd and n
  per setting and indicator at the last census, runtime, share of failed
  simulations); only runs of releases are committed (summaries of development
  versions, `*.9000_*.csv`, are ignored by git).

Requirements: git in the PATH, a compiler for R packages (Rtools on Windows), and
in the default library the packages that the old and the current versions import
(e.g. `Rcpp`, `ape`, `picante`, `treestats`). The scripts install the model from
the committed state of a commit (`git archive`), so uncommitted changes of the
model code are not simulated.

## Folders not in git

`archives/`, `cache/` and `lib/` are ignored by git (`.gitignore`). The scripts
always use those of the main checkout of the repository, also when they run in a
git worktree, which can be removed together with its files. `git clean -x`
deletes them.

- `archives/<version>_<commit>/`: the archive of a run: per-replicate tables,
  copies of the scripts, provenance, and raw outputs (`raw/`: all maps of the
  last six censuses and the last phylogeny of each simulation, enough to
  recompute the indicators at the end). Release archives are attached to the
  GitHub release without `raw/`: `tar -czf <archive>.tar.gz --exclude=raw <archive>`
  (about 10 MB); development and attribution runs stay local.
- `cache/<mode>_<version>_<commit>/`: the cache of a run of `run.R`; can be
  deleted once archived.
- `lib/<commit>/`, `lib/olddeps/`: libraries of `installVersion.R` and
  `installOldDeps.R`.

Archives in the main checkout of the maintainer's machine:

| Archive | Content |
|---|---|
| `0.3.1_852c61e` | run of the pre-release commit of 0.3.1 (not the tag) with the earlier script; 3010 simulations; not published |
| `0.3.1.9000_2331720` | run of #172 with the earlier script; 3010 simulations; not published |
| `legacy/` | runs of the earlier script `packageBehavior.R` (with the cache of 852c61e and a partial pilot of a92899d) and their earlier index |

## Running a version

    Rscript PhyloSim/data-raw/model-behavior/run.R main 12 current HEAD
    Rscript PhyloSim/data-raw/model-behavior/run.R main 12 0.3.1 v0.3.1
    Rscript PhyloSim/data-raw/model-behavior/makeData.R

`run.R <mode> <cores> <version> <commit or tag>` runs the design with a commit;
`<version>` is a name in `versions` of `adapters.R` (the release, or `current`).
It installs the commit and `indicatorCommit` into `lib/` if they are missing and
runs into `cache/<mode>_<version>_<commit>/`. Interrupted runs resume, and failed
simulations are run again. A complete run of mode `main` without failed
simulations is archived as `archives/<version>_<commit>/` (development
versions: the package version, e.g. `0.3.1.9000`), and its summary is written
to `runs/`. Otherwise (e.g. v0.2, where some simulations always crash, see below)
archive it with `Rscript PhyloSim/data-raw/model-behavior/archive.R <cache folder>`,
which keeps failed simulations as such; `archive.R <cache folder> partial` also
archives runs of some blocks. Runs of mode `test` are archived as
`test-<version>_<commit>`.

- Modes (`design.R`): `main` (30 replicates, the runs of record) and `test`
  (2 replicates, 300 generations; a check of the scripts, e.g.
  `run.R test 4 current HEAD gridShape fission`). Blocks after the commit restrict
  the run to these blocks, `part=i/n` to every n-th simulation.
- Mode `main` runs only with committed scripts, so that `RUNINFO.txt` names
  their commit, and a resumed run only with the scripts of its first start.
  Commit before a run; do not edit the scripts during a run.
- A run of mode `main` takes days (about 20 hours on 10 cores for the earlier,
  smaller design). Results do not depend on the number of cores, as each
  simulation has its own seed.
- `MB_MIN_FREE_GB` (default 16): free memory in GB that `run.R` waits for
  before each simulation (Windows only).
- `makeData.R` accepts one run per version in `runs/`, all with `indicatorCommit`;
  `makeData.R packageBehavior <archive>` also builds `packageBehavior`.

## Release

1. Freeze the model code of the release (all changes in `changes.csv`, see below,
   with their attribution) and set the release version in `PhyloSim/DESCRIPTION`,
   after which the archive is named (e.g. `0.3.2_<commit>`); commit. The run is
   made with this commit; the tag follows in step 5, on a commit with the same
   model code.
2. Run the design with this commit (`run.R main <cores> current HEAD`), and
   `pPD.R` (`Rscript PhyloSim/data-raw/model-behavior/pPD.R 12`).
3. Delete the summaries of development runs from `runs/`, and run
   `makeData.R packageBehavior <archive>`.
4. Update the notes per release in `PhyloSim/vignettes/ModelHistory.Rmd` and the
   description of the runs in the chapter "Emerging patterns".
5. Commit `runs/` and the data sets, tag the release, add the release to
   `versions` in `adapters.R` (a copy of `current`) and to `changes.csv`, and
   attach the archive to the GitHub release (in `archives/`:
   `tar -czf <archive>.tar.gz --exclude=raw <archive>`).

A new `indicatorCommit` (e.g. a changed `nullModel()` or `sac()`) means new runs
of all versions; the raw outputs in the archives allow recomputing the
indicators at the end without new simulations, but no script does this yet.

## When a pull request changes the model

A change that alters simulation results or the meaning of arguments gets a row in
`changes.csv` in its pull request: `kind` change, `version` the next release,
`reference` the issue, `type` (N, B, S or S+B), the `change` and its expected
effect on the `outputs`; new or renamed arguments are added to `design.R` and
`adapters.R` (see Extending). After the merge, a maintainer fills in `date` and
`commit` (the commit on `master`), and for the attribution `parent` (its parent
commit), `adapter` and `parentAdapter` (the names in `versions` of `adapters.R`
that run both, usually `current`) and, after stage 1, `stage2` (the blocks that
differ); before the release, the maintainer runs the attribution below and fills
in `evidence`.

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
  crashes without a message (exit status 5) in some simulations with `specRate` 8
  (seen with grids of 20 and 50 cells); these count as failed.

## Extending

- New parameter: add it to `base` (with the default of the current version) and,
  if varied, to a block of `design`; add it to `params` of the versions that
  support it in `adapters.R`. Seeds of existing settings do not change, as they
  depend only on the mode, the block, the parameters that differ from `base` and
  the replicate.
- Renamed argument: keep the design name and add it to `rename` of the version.
- New indicator: add it to `censusStats()` (at each census) or `behavior()` (at
  the end; a single value also goes into the census table, as `mpdLocal` and
  `cndd` in `archive.R`) in `indicators.R`, and, to show it in the key figure, to
  `indicators` in `PhyloSim/vignettes/ModelHistory.Rmd`. Earlier runs need a new
  run to get it.
- New block: add it to `design`.

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

## Transition (until the runs of #180)

The runs in `runs/` and in the archives above were made with the earlier script
`PhyloSim/data-raw/packageBehavior.R` (see tag v0.3.1), whose indicators were
computed with the simulated version; `makeData.R` therefore does not accept them
next to runs of `run.R`, and they are replaced by new runs of all versions. The
data sets in `PhyloSim/data/` are committed only at releases; the chapter
"Emerging patterns" is adapted to the new format of `packageBehavior` (all
censuses in `main`, no `pilot` and `cnddSpecies`, generations per block, CNDD in
the census) with the next run of the current version. Texts to update then: the
data documentation of `packageBehavior` and `modelHistory`
(`PhyloSim/R/data-documentation.R`), the indicator table of
`PhyloSim/vignettes/ModelHistory.Rmd` ("runs made with `run.R`"), and this
section. The evidence of #172 in `changes.csv` was counted with an earlier
version of the comparison in the vignette.
