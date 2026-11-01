---
name: orca-spectra-spin-occupier
description: >-
  Determine preferred spin states and broken-symmetry solutions through
  DELFIN's OCCUPIER workflow. Use when the user asks about spin states,
  high-spin vs low-spin, or spin contamination of a metal complex.
domains: code
---

# Spin states and broken symmetry via OCCUPIER

## Route

Spin-state comparisons run through DELFIN's OCCUPIER workflow (CONTROL
`OCCUPIER=auto` or the typed concept `occupier`; template
`delfin/tools/templates/_occupier.py`, entry point `occupier_stages`,
engine `delfin/occupier.py`).

**Before starting a run: propose the CONTROL, then confirm.** The agent
researches the system (charge states, candidate multiplicities, metal
basis need — see `metal_basis` / `NewGTO` handling in occupier.py:174-296),
writes the proposed CONTROL keys with one-line reasons, shows them to the
user, and starts OCCUPIER only after the user confirms. Never launch a
spin-state scan from an unconfirmed CONTROL.
What the workflow does:

- **Staged relaxation (FoB)**: the structure relaxes through fraction-
  of-bond stages (default fob_values = 0.0, 0.25, 0.5, 0.75, 1.0 in
  the template, each stage an `orca_opt` + `orca_freq` sub-pipeline
  built by `_make_fob_stage`); the full workflow builds an
  OCCUPIER_FOB job tree under the calc folder.
- **Preferred spin**: the workflow generates the candidate spin states
  (per charge/multiplicity you give) and compares their energies;
  `extract_preferred_spin` (delfin/copy_helpers.py, wired in cli.py:90)
  reads the preference.
- **Broken symmetry**: when spin contamination exceeds the threshold,
  the template's `check_contamination=True` branch inserts a broken-
  symmetry (BS) single point (`broken_sym=True`); the BS logic lives
  in delfin/bs_evolution.py. OCCUPIER scratch folders
  (`ox_step_*_OCCUPIER`, `red_step_*_OCCUPIER`, `initial_OCCUPIER`)
  are cleaned by delfin/cleanup.py:44-47.

### Broken-symmetry inputs: site order is not free

When OCCUPIER writes `%scf BrokenSym NA,NB` (delfin/occupier.py:1315-1319,
`broken_sym` comes from the `BS` key of the sequence entry), the ORCA
manual (6.1.1, §5.27) requires: **"site A" — the FIRST number — must be
the site with the LARGER number of unpaired electrons.** A `BrokenSym 1,2`
where site A holds 1 unpaired electron and site B holds 2 does not do what
the input says; it reconverges from a high-spin start built on the wrong
site. Check the atom order of the coordinates in the .inp: site A is
whatever the first listed metal center is. If the two centers carry equal
numbers, order does not matter.
For the `FlipSpin`/`FinalMs` mechanism (manual §5.27) the atom list counts
from **0**, and `FinalMs` MUST be given — ORCA cannot determine it itself.
Exchange couplings J from HS/BS pairs: use the Yamaguchi form
(Eq. 5.167, `J = -(E_HS - E_BS)/(<S²>HS - <S²>BS)`), which ORCA's manual
prefers as approximately valid over the whole coupling regime; report
<S²> of both solutions with it.

## What to report

- **The multiplicities compared** — the comparison is the result, not
  only the winner. Relative Gibbs energies in kcal/mol between the
  spin states.
- **S² and spin contamination** of every open-shell state (and of the
  BS solution). A heavily contaminated state (<S²> far from the ideal
  S(S+1)) is reported as such, not silently accepted.
- **The deciding stage**: which OCCUPIER stage / step produced the
  preference, so the number is traceable to a job.
- **Method, basis, solvation** with every value; a spin-state gap is
  only comparable within one level of theory.

## Honesty rules

- High-spin/low-spin gaps at hybrid-GGA DFT carry real method
  uncertainty (often several kcal/mol); a difference within ~2-3
  kcal/mol of the winner is "near-degenerate", not a definite order.
- Report imaginary frequencies of each state's Hessian; a
  "preferred" state that is not a minimum is not preferred.
- Deviations from OCCUPIER defaults (FoB schedule, threshold,
  metal basis NewGTO from occupier.py:174-296) need a reason in the
  same sentence.
- Never report a spin-state ordering from memory or a ligand-field
  argument alone — the comparison must come from the workflow's runs.

## Run `/check` before you report

`/check` (result_critic.py, docs/result_critic.md) reads each output
back and flags spin contamination (ground state via the UHF block;
per excited state where the TDDFT table header names an S^2 column),
imaginary modes on a supposedly-preferred state, and unfinished runs.
Report every ⚠️ verbatim; never report over a ❌ without saying why.
Report-only.
