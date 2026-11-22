# Development tools

Tools for changes of the package that are needed again in the future.

| File | Content | Maintained |
|---|---|---|
| `identity.R` | compares the complete output of two builds over many settings, for changes of the C++ model that are meant to keep results (e.g., to make it faster) | yes |
| `recruitment-check.R` | compares summary statistics of two builds over many seeds, for changes that alter the random numbers but should keep their distribution (written for the binned parent proposals of #207) | no |

Maintained tools run with the current version of the package, and a change of the package that breaks them also updates them. The others are templates that may need adjustments for later versions. Only tools that will be used again belong here, with a header that says what they do and how they are called. Prototypes, diagnoses and results of a specific change belong to its issue or pull request.
