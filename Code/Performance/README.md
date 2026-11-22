# Performance of the simulation model

Tools that are needed again when the C++ model changes. They are maintained: they must run with the current version of the package, and a change of the package that breaks them also updates them.

| File | Content |
|---|---|
| `identity.R` | compares the complete output of two builds over many settings, for changes that are meant to keep results (e.g., to make the model faster) |

## What belongs here

- Only tools that will be used again for future changes, one script per tool (settings, run and comparison in the same file), with a header that says what it checks and how it is called.
- Not: diagnoses, benchmarks, profiles and results of a specific change. They belong to its issue or pull request, e.g., the runtime diagnosis of #207 (profiles, paired benchmarks, counters) in the branch of PR #210 ([commit 5efcd39](https://github.com/TheoreticalEcology/EcoPhyloSim/tree/5efcd39/Code/Performance/runtime)).
