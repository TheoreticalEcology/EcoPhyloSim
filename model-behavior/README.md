# Model behavior across versions

Scripts that document how the biological behavior of EcoPhyloSim changes between
versions (issues #180, #181). Each version runs the same simulation design; the
indicators are always computed with the current package; the per-replicate results
are archived outside the repository and compared in `model-behavior.qmd`. The
vignette chapter "Emerging patterns" uses the run of the current version.

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
- `archive.R`: turns the cache of a run into the archive (tables and README).
- `addRun.R`: exports the last census per replicate of an archive to `runs/`.
- `model-behavior.qmd`: comparison of the versions and notes per version.

## Workflow (from the repository root)

    Rscript model-behavior/run.R main 12 current ../cache-current
    Rscript model-behavior/archive.R ../cache-current ../EcoPhyloSim-runs 0.3.1.9000 <commit>
    Rscript model-behavior/addRun.R ../EcoPhyloSim-runs/0.3.1.9000_<commit>
    quarto render model-behavior/model-behavior.qmd

Use `test` instead of `main` for a quick check of the scripts, and a new cache
folder per version and design (run.R stops if a cached simulation was run with
other parameters). The archive name must start with a version number.

## Older versions

Install the version from its tag into its own library, and its dependencies into
another one, e.g. for v0.3:

    git archive v0.3 PhyloSim | tar -x -C <folder>

- v0.2 and v0.3 (package PhyloSim) import `kulife` and `PhyloMeasures`, which
  are archived on CRAN (install `kulife_0.1-14`, `apTreeshape_1.5-0` with
  `cubature`, and `PhyloMeasures` from the CRAN archive).
- v0.3 lacks `src/RcppExports.cpp` (#102): run `Rcpp::compileAttributes()` first.
  v0.2 needs the same, and then the registration of the unused `callLeipzig` must
  be removed from the generated `src/RcppExports.cpp`. These are build fixes only.
- Run with `lib=<library of the version>;<library of the dependencies>` and
  `commit=<commit of the tag>`.

## Extending

- New parameter: add it to `base` (with the default of the current version) and,
  if varied, to a block of `design`; add it to `params` of the versions that
  support it in `adapters.R`. Older archives get the default in `addRun.R`.
- Renamed argument: keep the design name and add it to `rename` of the version.
- New indicator: add it to `censusStats()` or `behavior()` in `indicators.R`
  and, to show it in the key tables, to `indicators` in `model-behavior.qmd`.
  Older versions need a new run to get it.
- New block: add it to `design`; seeds depend only on parameter values.
- New release: copy `current` in `versions` to the version number before the
  model changes.
