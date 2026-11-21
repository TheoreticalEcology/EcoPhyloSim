# Code outside the package

Scripts and material that are not part of the package. They are not updated when the package changes, so older scripts may not run with the current version.

| Folder | Content |
|---|---|
| `Issues/` | reprexes and simulations for GitHub issues (created with its first file) |
| `Performance/` | runtime diagnosis and benchmarks of the simulation model |

## Adding files

Code for an issue goes to `Issues/<issue number>.R` (or a folder `Issues/<issue number>/` if it needs several files) and starts with a comment that says what it shows. Other code goes to the matching folder above. Raw simulation output and other files that the scripts regenerate are not committed (see the `.gitignore` and `README.md` of the subfolder that holds them, e.g. `Performance/runtime/`).
