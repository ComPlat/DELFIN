# TDDFT excited states
> Set up and interpret TDDFT excitation calculations — roots, TDA, triplets, and what the spectrum means.

The user wants excited states (UV/Vis, phosphorescence, TADF screening).

## 1. Scope the physics
1. How many states and which kind? Absorption → singlets; ISC/TADF →
   also triplets; emission → optimized excited-state geometry needed.
2. Functional choice matters more than usual: warn about charge-transfer
   states with pure GGAs; recommend range-separated (e.g. CAM/ωB97 family)
   for CT — and say WHY in one sentence.

## 2. Ground the block in the manual
Look up the %tddft section via `search_docs` before writing it. Canonical
keywords: `nroots`, `tda`, `triplets`, `maxdim`, `iroot` — verify against
the manual; `nstates` does NOT exist in ORCA's %tddft (a common
hallucination — never use it).

## 3. Build the input
- `nroots`: at least 1.5× the states of interest (root flipping).
- `tda true/false`: TDA is more robust for triplets; full TDDFT for
  quantitative absorption — state which you chose and why.
- `triplets true` when ISC/phosphorescence is in scope.
- Excited-state optimization: `iroot` selects the state — verify the
  syntax in the manual before proposing it.

## 4. Interpretation checks
- Report energies in eV AND nm, oscillator strengths, dominant
  orbital transitions.
- Flag suspicious results: near-zero f for the "bright" state, S1 above
  experimental absorption by >0.5 eV, CT states with tiny orbital overlap.
- ΔE(S1–T1) < ~0.3 eV → mention TADF relevance explicitly.

## Safety rules
- Keywords only with a `search_docs` citation.
- No submits without explicit user confirmation.

## Output format
**Goal (physics)** → **Method/functional (justification)** → **Input block**
(manual section cited) → **Analysis checks**.
