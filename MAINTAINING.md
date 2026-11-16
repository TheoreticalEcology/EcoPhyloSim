# Maintaining EcoPhyloSim

This page is for the EcoPhyloSim team (maintainers and collaborators with write access to this repository). The [contributing guidelines](CONTRIBUTING.md) apply to the team as well; this page only adds the tasks that require write access.

## Issues and pull requests

- Set the target release as the milestone of the issue or pull request.
- Team members may merge their own pull requests. If in doubt, request a review from another team member first.
- As an exception to the pull request workflow, trivial changes and changes outside the package folder `PhyloSim/` (e.g., `.github/`; in `PhyloSim/` only `NEWS.md` and the version bump after a release) can be committed directly to `master`.

## Releases

- With every version change in `PhyloSim/DESCRIPTION`, also set its `Date` field to the current date: `citation("EcoPhyloSim")`, also in the vignette, takes the year from it unless the package is installed from CRAN, and `R CMD check --as-cran` notes a `Date` older than a month.

## Branches

- Existing branches created before the [naming convention](CONTRIBUTING.md#pull-request-workflow) keep their names.
- A branch deleted after a merge can be restored with "Restore branch" on the pull request page.
