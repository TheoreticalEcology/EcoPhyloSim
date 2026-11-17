# Model behavior across versions

Scripts that document how the biological behavior of EcoPhyloSim changes between
versions (issues #180, #181). Each version runs the same simulation design; the
indicators are computed with the current package (in the runs of 0.3.1 and #172,
made with the earlier script, with the simulated version). The runs give the data
set `packageBehavior` (chapter "Emerging patterns" of the main vignette, from the
run of the current version) and the data set `modelHistory` (vignette "Model
history", from the runs of all versions). The `packageBehavior` of release 0.3.1
was still made with the earlier script `PhyloSim/data-raw/packageBehavior.R`
(see tag v0.3.1), which these scripts replace.

The per-replicate archives are not tracked in git (they are too large; `runs/`
keeps one summary per run), but everything needed to regenerate them is: the
scripts, the design, the commits of the simulated versions, and in each run and
archive the exact command, platform and package versions. The runs of 0.3.1
(852c61e) and of #172 (2331720) were still made with the earlier script; copies of
it are in their archives. The data sets in `PhyloSim/data/` are committed only at
releases. `makeData.R packageBehavior <archive>` builds `packageBehavior` in a
new format (all censuses in `main`, no `pilot` and no CNDD coefficients of single
species); the chapter "Emerging patterns"
is adapted to it with the next run of the current version (#180), so until then
build only `modelHistory`.

## Files

- `design.R`: blocks of settings, replicates, run lengths, censuses and seeds.
  Parameter names are those of release 0.3.1 (e.g. `density`), and the current
  names for parameters added later.
- `indicators.R`: all indicators of one simulation (`behavior()`).
- `adapters.R`: the model versions (`versions`: package name, supported
  parameters, renamed arguments) and `adapt()`, which maps a task to the
  arguments of `createCompletePar()` of a version.
- `installVersion.R`, `installOldDeps.R`: install the package of any commit or
  tag into its own library, and the archived dependencies of the old versions.
- `run.R`, `simulate.R`: run the design with one installed version, each
  simulation in its own R process.
- `archive.R`: turns the cache of a run into the archive (tables, raw outputs and
  README).
- `addRun.R`: writes the summary of an archive (mean, sd and n per setting and
  indicator at the last census, and the share of failed simulations) to `runs/`.
- `makeData.R`: builds `modelHistory` from `runs/` and `changes.csv`, and, if
  asked for, `packageBehavior` from the archive of the current version.
- `changes.csv`: releases and the changes of the model code that change results,
  with their effects, the evidence and the settings of the attribution runs.
- `attribution.R`: tests single changes against their parent commit (stages 1
  and 2, see below).

## Workflow (from the repository root)

    Rscript model-behavior/installVersion.R HEAD ../lib-current
    Rscript model-behavior/run.R main 12 current ../cache-current lib=../lib-current
    Rscript model-behavior/archive.R ../cache-current ../EcoPhyloSim-runs 0.3.1.9000 <commit>
    Rscript model-behavior/addRun.R ../EcoPhyloSim-runs/0.3.1.9000_<commit>
    Rscript model-behavior/makeData.R
    # only once the chapter "Emerging patterns" is adapted to the new format (#180):
    Rscript model-behavior/makeData.R packageBehavior ../EcoPhyloSim-runs/0.3.1.9000_<commit>

- Modes (`design.R`): `main` (30 replicates, the runs of record; its censuses
  also give the time to equilibrium) and `test` (2 replicates, 300 generations,
  a quick check of the scripts and stage 1).
- `run.R` takes the simulated commit and package name from `INSTALLINFO.txt` of
  the library (written by `installVersion.R`). The current EcoPhyloSim computes
  the indicators: install it in the default library, or, to record its commit,
  with `installVersion.R` into a library put first in `R_LIBS`; `SESSION-*.txt`
  then contains its `INSTALLINFO.txt`.
- A `main` run of the earlier, smaller design took about 20 hours on 10 cores of a
  Windows PC (2026, mostly the density blocks); the current design takes longer.
  Results do not depend on the number of cores, as each simulation has its own seed.
- Use a new cache folder per version and design (run.R stops if a cached
  simulation was run with other parameters, generations or seeds, e.g. in caches
  of earlier script versions). Interrupted runs resume; failed simulations are
  run again.
- Do not edit the scripts while a run is going: Rscript reads a script while it
  runs. For long runs, use a separate git checkout.
- `run.R` runs mode `main` only from a git checkout with committed scripts, so
  that `RUNINFO.txt` names their commit.
- The archive name must start with a version number; `<commit>` must match the
  simulated commit in `RUNINFO.txt`. Existing archives are not overwritten (delete
  the folder of a failed attempt).
- Only the archives of release runs (at the release tag) are attached to the
  GitHub release, without the raw outputs (`raw/`: species maps and phylogenies
  of all simulations, trait maps of the first replicates), which stay local:
  `tar -czf <archive>.tar.gz --exclude=raw <archive>`. Development and
  attribution runs stay local.

## Older versions and single commits

    Rscript model-behavior/installOldDeps.R ../lib-olddeps
    Rscript model-behavior/installVersion.R v0.3 ../lib-v0.3 ../lib-olddeps
    Rscript model-behavior/run.R main 12 0.3 ../cache-v0.3 "lib=../lib-v0.3;../lib-olddeps"

- `installVersion.R` applies build fixes only: v0.2 and v0.3 lack
  `src/RcppExports.cpp` (#102), which is generated with `Rcpp::compileAttributes()`;
  for v0.2 the lines with the unused `callLeipzig` are then removed from it.
- v0.2 imports `kulife`, v0.3 also `apTreeshape` and `PhyloMeasures`, which are
  archived on CRAN; `installOldDeps.R` installs `kulife` 0.1-14, `apTreeshape`
  1.5-0 with `cubature`, and `PhyloMeasures` 2.1. Their other dependencies (e.g.
  `ape`, `picante`, `quantreg`) must be in the default library.
- Separate the libraries with `;` on Windows and `:` otherwise.
- The version name (`0.3` above) selects the adapter in `adapters.R`.

## Attribution of changes

The columns `parent`, `adapter`, `parentAdapter`, `stage2` and `stage2Commit` of
`changes.csv` give, for each change that is tested, its parent commit, the
adapters of both, the blocks for stage 2 and, if stage 2 compares another commit
with the parent (e.g. when the change itself crashes), that commit.

    Rscript model-behavior/attribution.R run 1 6 ../attribution
    Rscript model-behavior/attribution.R compare 1 ../attribution
    Rscript model-behavior/attribution.R run 2 6 ../attribution
    Rscript model-behavior/attribution.R compare 2 ../attribution

Stage 1 runs mode `test` for each change and its parent with the same seeds and
reports per block how many simulations have no result (failed or unsupported) and
in how many the census indicators of the communities, species ages, tree
statistics and the Colless index differ. Stage 2 runs mode `main` for the blocks
in `stage2` and compares the last census with Welch tests (false discovery rate
per change). With 30 replicates, only differences of about 0.7 standard
deviations or more are detected reliably, so stage 2 can show a change of the
distributions but not prove that there is none. Enter the results in the column
`evidence` of `changes.csv`.

## Provenance

Each run records in its cache folder, and `archive.R` copies to the archive:

- `RUNINFO.txt`: start, mode, version, number of simulations and cores, libraries,
  simulated package and commit, and commit of the scripts (`git describe --dirty`);
- `SESSION-<start>.txt`: the command, the platform, the R session of the
  indicators, `INSTALLINFO.txt` of the libraries used, and all packages with
  their versions in these libraries;
- copies of the scripts at the start of the run.

Since #175 (R's random number generator), results for a given seed no longer
depend on the random number library of the compiler; for older versions they may
depend on platform and compiler. The archives are kept in `EcoPhyloSim-runs` next
to the repository.

## Extending

- New parameter: add it to `base` (with the default of the current version) and,
  if varied, to a block of `design`; add it to `params` of the versions that
  support it in `adapters.R`. Seeds of existing settings do not change, as they
  depend only on the mode, the block, the parameters that differ from `base` and
  the replicate. Older archives get the
  default in `addRun.R`.
- Renamed argument: keep the design name and add it to `rename` of the version.
- New indicator: add it to `censusStats()` or `behavior()` in `indicators.R`
  and, to show it in the key figure, to `indicators` in
  `PhyloSim/vignettes/ModelHistory.Rmd`. Older versions need a new run to get it.
- New block: add it to `design`.
- New release:
  1. Before the model changes: copy `current` in `versions` to the version number.
  2. In `changes.csv`, fill date and commit of the release, add its changes and,
     for the changes to test, `parent`, `adapter`, `parentAdapter` and `stage2`;
     run `attribution.R` and enter the results in `evidence`.
  3. Install the release tag with `installVersion.R`, run the design, archive it
     as `<release>_<commit>`, run `addRun.R`, delete the summaries of
     development runs of that release from `runs/`, and run `makeData.R`
     (with `packageBehavior <archive>` to rebuild `packageBehavior`).
  4. Update the notes per release in `PhyloSim/vignettes/ModelHistory.Rmd` and
     the description of the runs in the chapter "Emerging patterns".
  5. Commit `runs/` and the data sets, and attach the archive to the GitHub release
     (without `raw/`, see Workflow).
