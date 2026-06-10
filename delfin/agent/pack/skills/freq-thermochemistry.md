# Frequencies & thermochemistry
> Run/interpret frequency calculations: imaginary modes, thermochemistry, and physical-plausibility checks.

The user wants frequencies — for minima verification, transition states,
IR spectra, or Gibbs energies.

## 1. What the frequencies are FOR decides the protocol
- **Minimum check**: zero imaginary modes required.
- **Transition state**: EXACTLY ONE imaginary mode — and it must look like
  the reaction coordinate (describe the motion from the normal mode).
- **Thermochemistry**: frequencies must come from the SAME method/basis
  as the geometry — mixed levels invalidate ZPE/G corrections; flag it.

## 2. Ground the input in the manual
%freq keywords via `search_docs` (e.g. `scalfreq`, numerical vs analytic
options) — cite the section; never recite keyword names from memory.

## 3. Imaginary modes — the decision tree
- 1 small imaginary mode (<~50i cm⁻¹) at a "minimum": usually a floppy
  rotor / numerical noise → options: tighter opt, displaced restart along
  the mode, or document-and-accept for screening. State the trade-offs.
- Large imaginary mode at a "minimum": NOT a minimum — go back to
  geometry (use geometry-convergence-debug skill).
- TS with extra imaginary modes: displace along the spurious mode and
  re-optimise; do not just delete the mode from the analysis.

## 4. Thermochemistry checks (physical plausibility)
- T, p, standard state: state them explicitly (ORCA defaults vs 1M
  solution standard state — the conversion matters for ΔG in solution).
- ZPE/thermal corrections must be from the same level as stated above.
- Sanity: ΔG−ΔE difference typically a few kcal/mol; tens of kcal/mol →
  suspect a unit or geometry problem. Entropies of small molecules have
  known magnitudes — flag wild outliers.

## Safety rules
- Never report Gibbs energies from a structure with unexplained imaginary
  modes without an explicit warning in the answer.
- Keywords only with manual citation; no submits without confirmation.

## Output format
**Purpose** → **Input** (section cited) → **Mode analysis** →
**Thermochemistry with plausibility checks**.
