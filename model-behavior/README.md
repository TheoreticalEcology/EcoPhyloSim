# Model behavior across versions

Scripts that document how the biological behavior of EcoPhyloSim changes between
versions (issues #180, #181). Each version runs the same simulation design; the
indicators are always computed with the current package; the per-replicate results
are archived outside the repository and compared in the vignette "Model history"
of the package. The chapter "Emerging patterns" of the main vignette is to use the
run of the current version (#180); so far it uses `packageBehavior`, made with
`PhyloSim/data-raw/packageBehavior.R`.

The full per-replicate archives are not tracked in git (they are too large; `runs/`
keeps only the last census per replicate), but everything needed to regenerate them is: the scripts, the design, the commits of
the simulated versions, and in each run and archive the exact command, platform and
package versions.

## Files

- `design.R`: blocks of settings, replicates, run lengths, censuses and seeds.
  Parameter names are those of the design (e.g. `density`), not necessarily the
  arguments of the current version.
- `indicators.R`: all indicators of one simulation (`behavior()`).
- `adapters.R`: the model versions (`versions`: package name, supported
  parameters, renamed arguments) and `adapt()`, which maps a task to the
  arguments of `createCompletePar()` of a version.
- `run.R`, `simulate.R`: run the design with one version (other versions in a
  separate R process with their own library).
- `installVersion.R`, `installOldDeps.R`: install the package of any commit or tag,
  and the archived dependencies of the old versions.
- `archive.R`: turns the cache of a run into the archive (tables and README).
- `addRun.R`: exports the last census per replicate of an archive to `runs/`.
- `makeData.R`: builds the package data set `modelHistory`
  (`PhyloSim/data/modelHistory.rda`) from `runs/` and `changes.csv`, which the
  vignette "Model history" (`PhyloSim/vignettes/ModelHistory.Rmd`)
  shows.
- `changes.csv`: releases and the changes of the model code that change results,
  with their effects and the evidence.
- `attribution.R`, `attribution.csv`: tests single changes against their parent
  commit (stages 1 and 2, see below).

## Workflow (from the repository root)

Run the design with the current package (installed in the default library):

    Rscript model-behavior/run.R main 12 current ../cache-current
    Rscript model-behavior/archive.R ../cache-current ../EcoPhyloSim-runs 0.3.1.9000 <commit>
    Rscript model-behavior/addRun.R ../EcoPhyloSim-runs/0.3.1.9000_<commit>
    Rscript model-behavior/makeData.R

- Modes (`design.R`): `main` (30 replicates, the runs of record), `pilot` (5
  replicates, many censuses, to check the time to equilibrium), `test` (2
  replicates, 300 generations, a quick check of the scripts and stage 1).
- A `main` run of the earlier, smaller design took about 20 hours on 12 cores of a
  Windows PC (2026, mostly the density blocks); the current design takes longer.
  Results do not depend on the number of cores, as each simulation has its own seed.
- Use a new cache folder per version and design (run.R stops if a cached
  simulation was run with other parameters). Interrupted runs resume.
- Do not edit the scripts while a run is going: Rscript reads a script while it
  runs. For long runs, run a copy of the scripts or a separate checkout.
- Commit the scripts before a run of record, so that `RUNINFO.txt` names a clean
  commit.
- The archive name must start with a version number.

## Older versions and single commits

    Rscript model-behavior/installOldDeps.R ../lib-olddeps
    Rscript model-behavior/installVersion.R v0.3 ../lib-v0.3 ../lib-olddeps
    Rscript model-behavior/run.R main 12 0.3 ../cache-v0.3 "lib=../lib-v0.3;../lib-olddeps" commit=v0.3

- `installVersion.R` installs the package of any commit or tag into its own
  library and writes `INSTALLINFO.txt` there. It applies build fixes only: v0.2
  and v0.3 lack `src/RcppExports.cpp` (#102), which is generated with
  `Rcpp::compileAttributes()`; for v0.2 the lines with the unused `callLeipzig`
  are then removed from it.
- v0.2 imports `kulife`, v0.3 also `apTreeshape` and `PhyloMeasures`, which are
  archived on CRAN; `installOldDeps.R` installs `kulife` 0.1-14, `apTreeshape`
  1.5-0 with `cubature`, and `PhyloMeasures` 2.1. Their other dependencies (e.g.
  `ape`, `picante`, `quantreg`) must be in the default library.
- Separate the libraries with `;` on Windows and `:` otherwise. For commits before
  the renaming to EcoPhyloSim (#164), add `package=PhyloSim`.
- The version name (`0.3` above) selects the adapter in `adapters.R`.

## Attribution of changes

`attribution.csv` lists the changes that may alter results, each with its parent
commit (or the start of a range of commits, see `note`), the adapters of both in
`adapters.R` and the blocks for stage 2.

    Rscript model-behavior/attribution.R run 1 6 ../attribution
    Rscript model-behavior/attribution.R compare 1 ../attribution
    Rscript model-behavior/attribution.R run 2 6 ../attribution
    Rscript model-behavior/attribution.R compare 2 ../attribution

Stage 1 runs mode `test` for each change and its parent with the same seeds and
reports per block how many simulations failed and in how many the census
indicators of the communities, species ages, tree statistics and the Colless index
differ. Stage 2 runs mode `main` for the blocks in `stage2`
and compares the last census with Welch tests (Holm correction). Enter the results
in the column `evidence` of `changes.csv`.

## Provenance

Each run records in its cache folder, and `archive.R` copies to the archive:

- `RUNINFO.txt`: start, mode, version, number of simulations and cores, libraries,
  commit of the simulated package and of the scripts (`git describe --dirty`);
- `SESSION-<start>.txt`: the command, the platform, the R session of the
  indicators, `INSTALLINFO.txt` of the simulated version, and all packages with
  their versions in the libraries used;
- copies of the scripts at the start of the run.

Since #175 (R's random number generator), results for a given seed no longer
depend on the random number library of the compiler; for older versions they may
depend on platform and compiler.
The archives are kept in `EcoPhyloSim-runs` next to the repository and are to be
attached to the GitHub releases.

## Extending

- New parameter: add it to `base` (with the default of the current version) and,
  if varied, to a block of `design`; add it to `params` of the versions that
  support it in `adapters.R`. Older archives get the default in `addRun.R`.
- Renamed argument: keep the design name and add it to `rename` of the version.
- New indicator: add it to `censusStats()` or `behavior()` in `indicators.R`
  and, to show it in the key figure, to `indicators` in `PhyloSim/vignettes/ModelHistory.Rmd`.
  Older versions need a new run to get it.
- New block: add it to `design`; seeds depend only on parameter values.
- New release: copy `current` in `versions` to the version number before the
  model changes, add the release and its changes to `changes.csv`, and run the
  design with the release.
