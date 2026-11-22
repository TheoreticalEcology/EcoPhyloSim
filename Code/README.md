# Code outside the package

Scripts and material that are not part of the package. Except for the tools marked as maintained in `Development/` and `Performance/`, they are not updated when the package changes, so older scripts may not run with the current version.

| Folder | Content |
|---|---|
| `Development/` | tools for changes of the package, see the [README](Development/README.md) there |
| `Issues/` | reprexes and simulations for GitHub issues |
| `Performance/` | tools to measure the runtime of the simulation model and options for speed-ups, see the [README](Performance/README.md) there |

## Adding files

Code for an issue goes to `Issues/<issue number>.R` (or a folder `Issues/<issue number>/` if it needs several files) and starts with a comment that says what it shows. Other code goes to the matching folder above. Raw simulation output and other files that the scripts regenerate are not committed (see the `.gitignore` and `README.md` of the subfolder that holds them).
