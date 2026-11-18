# Scripts that create the data sets of the package

The data sets in `PhyloSim/data/` are created by the scripts in this folder
(except the three older example simulations in the table) and are not edited by
hand: to change one, change its script; a maintainer rebuilds the data set
(`modelHistory` and `packageBehavior` only with new runs of releases, see
[model-behavior/README.md](model-behavior/README.md)). This folder is not
included in the package that users install (it is excluded in
`PhyloSim/.Rbuildignore`).

| Data set | Content | Created by |
|---|---|---|
| `modelHistory` | behavior of all versions of the model and the changes between them (vignette "Model history") | `PhyloSim/data-raw/model-behavior/makeData.R` (not yet in `PhyloSim/data/`, see "Transition" in its README) |
| `packageBehavior` | patterns that emerge in the current version (chapter "Emerging patterns" of the vignette "EcoPhyloSim model") | `PhyloSim/data-raw/model-behavior/makeData.R` (the current file is still from `PhyloSim/data-raw/packageBehavior.R` of release 0.3.1, see "Transition" in its README) |
| `pPD.pvalues`, `pPD.positions` | example data of `plotPhylogeneticDispersion()` | `PhyloSim/data-raw/model-behavior/pPD.R` |
| `simu.envcom.local`, `simu.neutral.local`, `simu.neutral.global` | example simulations | no script in the repository (made with `runSimulation()`, see their help) |

## Folders

- [`model-behavior/`](model-behavior/README.md): simulations that document how
  the model behaves and how its behavior changes between versions, and the
  table of model changes (`changes.csv`) that pull requests changing results or
  the meaning of arguments must update; its `runs/` holds the committed
  summaries of the release runs.
- `cache/` (ignored by git) is left over from the earlier script
  `packageBehavior.R`.

## Adding a data set

Put its script into a folder here (a new one if it fits none), describe it in
that folder's `README.md` (what it creates, how to run it, how long it takes,
which outputs are not kept in git) and in the table above, and document the data
set in `PhyloSim/R/data-documentation.R` with `@source` naming the script. See
also the [pull request workflow](../../CONTRIBUTING.md#pull-request-workflow).
