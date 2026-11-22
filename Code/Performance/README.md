# Performance of the simulation model

Tools to measure the runtime of the simulation model, and what is known about further speed-ups. Diagnoses and benchmark results of a specific change belong to its issue or pull request, e.g., the runtime diagnosis of #207 at [commit 5efcd39](https://github.com/TheoreticalEcology/EcoPhyloSim/tree/5efcd39/Code/Performance/runtime) of PR #210. Whether a change keeps results is checked with `Code/Development/identity.R`.

| File | Content | Maintained |
|---|---|---|
| `benchmark.R` | compares the CPU time of two builds over settings that cover the main costs of the model | yes |
| `profiling/` | sampling profiles of single settings: where the time goes, per function (Windows, Very Sleepy and Rtools; usage in the header of `profile.ps1`) | no |

Maintained tools run with the current version of the package, and a change of the package that breaks them also updates them. The others are templates that may need adjustments for later versions.

## Further options for speed-ups

Not implemented in #207; all would keep results unless noted.

- **Recruitment**: compute the neighbourhood of the target cell once per death instead of once per tried parent.
- **Fitness on reproduction**: the dispersal weights are still recomputed for every parent (only 1.2x faster in #207).
- **Incremental density sums**: update the stored sums of the neighbours by the change of one cell. Exact for the `"species"` kernel; for `"absolute"` and `"gaussian"` the summation order changes, so results change in the last bits. For `"gaussian"`, caching the kernel values of each cell would save most `exp()` calls and keep results.
- **Memory layout for large grids**: one contiguous allocation for the grid; per-run constants moved out of `Individual`; species means computed when read instead of after every birth and death.
- Not worth it: alias tables or Gillespie-type sampling (changes results, little gain), threads within a run (the death-replacement chain is sequential), `-ffast-math` or `-march=native` (results would depend on the machine).
