# #143: offspring inherit the incipient state

Decisions for protracted speciation (`protracted > 0`) and the numbers they rest on.

## Problem and reproduction

Before the fix, a speciation event marked one individual as incipient, and only this individual could complete speciation `protracted` generations later; its offspring were not incipient. `richness.R` runs the neutral model with global dispersal (30 x 30 cells, `specRate = 4`, 3000 generations, seeds 1:10) before (commit 619b469) and after the fix. Mean richness (sd):

| protracted | fission 0, before | fission 0, after | fission 1, before | fission 1, after | fission 2, before | fission 2, after |
|---|---|---|---|---|---|---|
| 0 | 24.2 (3.1) | 24.2 (3.1) | 67.2 (5.5) | 67.2 (5.5) | 68.0 (5.2) | 68.0 (5.2) |
| 1 | 10.8 (2.7) | 20.1 (2.2) | 58.3 (5.7) | 63.2 (5.3) | 58.9 (5.7) | 64.0 (4.8) |
| 2 | 5.1 (2.6) | 18.1 (1.9) | 44.7 (5.2) | 60.5 (5.6) | 44.4 (4.9) | 60.8 (4.8) |
| 5 | 1.1 (0.3) | 14.7 (2.0) | 12.4 (3.6) | 51.0 (5.0) | 11.5 (3.7) | 51.4 (5.1) |

Before the fix, richness without fission follows the Ewens expectation for a speciation rate reduced by the survival of a single founder, theta = 4 exp(-p): 22.2, 10.0, 4.6 and 1.2 species for p = 0, 1, 2 and 5. With fission, the marked individuals are many, so the loss was smaller.

## Decisions

The three design proposals made for this issue agreed on the main points below.

- **Inheritance.** Each event starts an incipient lineage with its own tag, stored in the individual (`Individual::m_incipTag`) and copied to its offspring; `Landscape::m_incipient` keeps the start generation of each lineage. The lineage becomes one new species `protracted` generations after the event if any of its members is alive, with one trait jump (#114) for all of them. The per-cell tag and the age counter of #208 were replaced by this, so that offspring born in a cell keep the tag of their parent.
- **Event on an incipient individual (without fission).** It starts a lineage nested in the one of the individual. When the outer lineage completes, the members of the nested lineage join the new species and keep their own tag; the nested lineage completes later as a daughter of it. If the founding branch of the outer lineage dies out but a nested lineage survives, the outer lineage still completes, since its descendants are alive. Ignoring such events would lose them (as in #208); letting the nested lineage split off from the old parent species would give a wrong topology. A second event on the same individual in the same generation has no effect, as without protracted speciation (#208); nesting it would give a new species without individuals. With fission, an event on an incipient individual still splits its lineage, and both parts keep the start of the lineage (#208).
- **Date of the split.** In the phylogeny, a new species splits off when its speciation starts (`Species::m_Date_of_Split`), not when it completes. In the protracted birth-death model, branching points lie at the times incipient species are produced (Etienne & Rosindell 2012, Syst. Biol. 61:204, Fig. 1, https://pmc.ncbi.nlm.nih.gov/articles/PMC3280041/). Dating at completion would shift every split by `protracted` generations towards the present. The age used by the Red Queen advantage still counts from the completion of speciation (`m_Date_of_Emergence`), since the advantage is meant for new species.
- **Unchanged results.** With `protracted = 0`, results for a given seed are identical to commit 619b469 (species, traits and phylogeny of 12 settings with global and local dispersal, fission 1 and 2, density dependence, environment, trait jump, Red Queen, many events per generation, several recorded generations, pruned and full trees); two events in one generation on the same cell still give one species there. Only the hash `protracted` in `test-reproducibility.R` changes.

Not verified here: whether Rosindell et al. (2010, Ecol. Lett. 13:716) let descendants of an incipient individual carry the incipient state (a comment on #143 says they do), and how they treat an event on an incipient individual; the full text was not accessible.
