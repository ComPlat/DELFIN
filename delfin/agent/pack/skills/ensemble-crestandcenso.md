---
domains: chemistry
---
# Conformer ensemble refinement (CENSO / weighted ensembles)
> Take a CREST/GOAT conformer ensemble, refine and rank it consistently at
> one level of theory, weight it at a stated temperature, and hand the
> best conformers to the DFT workflow.

This skill starts where the conformer-search skill
(`conformer-search.md`) ends: its charge/spin/topology rules and its
ensemble-reading section apply unchanged. What is added here: the
REFINEMENT side — CENSO levels and their cost, Boltzmann weighting at a
stated temperature, degenerate and duplicate conformers, and handing the
best structures onward.

## 1. What CENSO does, and what each level costs
CENSO ranks an ensemble at a composite-DFT level through sequential
phases; the number of structures per phase shrinks. In DELFIN's
`censo_sort` adapter (`delfin/tools/adapters/censo.py`):
- `prescreening` (default on): single points + cheap geometries on every
  conformer — the bulk of the cost is here, one pass over all structures.
- `screening` (default on): refined single points and constrained
  geometries on the surviving set (threshold, default 4.0 kcal/mol —
  DELFIN's `censo.py` adapter; the phase scheme and threshold follow
  the CENSO paper, Lit. [1]).
- `optimization` (default OFF): full optimizations per conformer — the
  most expensive phase; turn it on only when geometries themselves are
  the result (Raman, NMR, steric arguments).
- Functional default `r2scan-3c`; solvent via `solvent=` (ALPB).
  r2SCAN-3c: Grimme, Hansen, Ehlert, Mewes, J. Chem. Phys. 154, 064103
  (2021), DOI 10.1063/5.0040021, Lit. [3].

Say in the report which phases ran and how many conformers entered and
left each one. "The ensemble was refined" alone is not a result.

## 2. Boltzmann weighting at the stated temperature
- Weights come from ONE level of theory only. xtb energies weight an xtb
  ensemble; CENSO's r2scan-3c energies weight the refined ensemble.
  Never mix (see the trap task in this package's benchmark file).
- The temperature is a stated parameter in Kelvin — never a silent
  298.15 K default (DELFIN's `boltzmann_weights` takes T as a required
  argument, `crest_ensemble.py`): write "at 298.15 K" (or whatever was
  asked) next to the weights. p_i ~ exp(-dE_i/RT) is the Boltzmann
  distribution; the ensemble-weighting convention follows Lit. [1].
- Degeneracy enters explicitly: p_i ~ g_i exp(-dE_i/RT) with g_i the
  degeneracy of conformer i (mirror-image pairs g=2 — the
  symmetry-number convention of the ensemble-weighting formalism, Lit.
  [1]). `crest_ensemble.collapse_degeneracies`
  groups frames degenerate within a tolerance (default 0.01 kcal/mol —
  DELFIN's `collapse_degeneracies` in `crest_ensemble.py`)
  and sets each group's degeneracy to the GROUP SIZE (the number of
  frames in it — it does NOT preserve a per-frame g a search tool may
  have reported); `boltzmann_weights` multiplies by it. With bare
  Boltzmann factors over the un-collapsed frames the two entries of a
  degenerate pair give the same 2x weight, which is why g = group size
  is consistent.
- Sum of all weights is exactly 1; an ensemble whose frames carry
  energies tagged from two levels, or only PARTIALLY tagged, is refused
  by `boltzmann_weights` (ValueError) — that refusal is the correct
  answer, not something to route around by hand-weighting.

## 3. Duplicates and degenerate conformers
- Two frames with the SAME energy within tolerance are degenerate (often
  mirror images) — they are counted, not collapsed silently.
- Duplicates (identical geometry) are a separate thing: check
  connectivity/geometry before claiming an ensemble has N conformers;
  the RMSD filtering that CREST already did does not catch everything.

## 4. Handing the best conformers to the DFT workflow
- Take the lowest-energy structures that together carry the requested
  population mass (e.g. all with weight >= 0.05, mirroring DELFIN's
  OCCUPIER threshold convention in `delfin/occupier_auto.py:212`).
- Write each as its own .xyz, in the same charge/multiplicity as the
  search, and run the DFT workflow (CONTROL `method=` templates) per
  conformer; ranks and weights travel with the structures in the report.
- A rank from xtb is a PRESELECTION, never the final ordering: after DFT
  refinement the ordering can change — say which level each ordering
  came from.

## 5. When a CREST run has converged
- The search is done when the conformer count and the ensemble's
  composition stop changing across the search's iteration blocks
  (`crest.out` shows the final ensemble table) and the requested energy
  window (`ewin`, default 6 kcal/mol in DELFIN's
  `build_crest_cmd`, `delfin/tools/adapters/crest.py:24`; the CREST
  energy-window concept follows Lit. [2]) is covered.
- Constraints (CREST modes: default, quick, constrained, NCI): a
  constrained or NCI run explores a DIFFERENT space than the default;
  state which mode ran and what it excluded.

## Output format
**Request and system** (charge, unpaired electrons, solvent, temperature)
-> **Input** (skill reference: conformer-search.md; Literature: [1]
CENSO — Grimme et al., J. Phys. Chem. A 125, 4039 (2021), DOI
10.1021/acs.jpca.1c00971; [2] CREST — Pracht, Bohle, Grimme, Phys.
Chem. Chem. Phys. 22, 7169 (2020), DOI 10.1039/C9CP06869D; [3]
r2SCAN-3c — Grimme, Hansen, Ehlert, Mewes, J. Chem. Phys. 154, 064103
(2021), DOI 10.1063/5.0040021; [4] GFN2-xTB — Bannwarth, Ehlert,
Grimme, J. Chem. Theory Comput. 15 (3), 1652-1671 (2019), DOI
10.1021/acs.jctc.8b01176) -> **Refinement and weighting**
(phases run, counts in/out per phase, level per ensemble, Boltzmann
table with degeneracies at the stated T) -> **What it means and what it
does not** (xtb preselection vs DFT ordering) -> **Handed-over
structures** (paths, weights).
