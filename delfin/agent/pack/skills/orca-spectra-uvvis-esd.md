---
name: orca-spectra-uvvis-esd
description: >-
  Compute UV-Vis absorption through DELFIN's _esd workflow: TD-DFT
  states, oscillator strengths, solvation, n->pi* vs pi->pi* character.
  Use when the user asks for an absorption spectrum, excitation
  energies, or an excited state of a small molecule.
domains: code
---

# UV-Vis via TD-DFT and the _esd workflow

## Route

Absorption spectra run through DELFIN's ESD workflow (CONTROL
`method=esd`, template `delfin/tools/templates/_esd.py`). What the
template encodes:

- `esd_states(charge, mult, method, basis, states=[...], mode=...)`:
  `mode="tddft"` runs S0 as the trunk (opt, optional freq) and every
  excited state as a parallel branch with `orca_tddft` (nroots =
  highest requested state index, iroot per state, triplet flag from
  the T-prefix); `mode="deltascf"` builds a sequential chain.
- `include_freq=True` adds a frequency run per state;
  `include_isc=True` adds inter-system crossing rate steps between
  singlet/triplet pairs.
- For redox-coupled workflows the first S0 job always runs with FREQ
  (esd_input_generator.py:791) so the ground-state Gibbs energy
  exists for the cycle.

For a plain UV-Vis absorption of a small closed-shell molecule, the
S0 trunk + TD-DFT branch is what you want; you do not need ISC rates
unless the question is about phosphorescence or intersystem crossing.

## What to request and what to report

- **States**: request enough roots to see the visible band (typical
  small molecule: 5-10 singlet roots). State the number requested.
- **Report per state**: energy in eV AND nm (1 eV = 8065.54 cm-1,
  ORCA manual §10.3), oscillator strength f, and the character
  (which orbitals the excitation connects — n->pi* is dark, f ≈
  10^-3-10^-2; pi->pi* is bright, f ≈ 0.1-1).
- **Typed tools**: `extract_excited_states` returns the transition
  table plus first_bright / brightest_visible / strongest_visible with
  their definitions (bright = fosc >= 0.01, visible = 380-780 nm);
  `plot_uvvis_spectrum` draws the broadened curve.
- **Solvation**: solvated TD-DFT (CPCM etc.) belongs in the workflow
  defaults; state the solvent with the value. A gas-phase vertical
  excitation of a charge-transfer state can shift by >0.5 eV against
  the solvated one.

## Honesty rules

- A vertical TD-DFT excitation is not the experimental 0-0 band; label
  which one you report. TD-DFT typically overestimates n->pi*
  vertical energies of carbonyls.
- Unconverged TD-DFT (SCF or root convergence) is reported, never
  hidden.
- The E_00 energy and excited-state redox potentials that DELFIN's
  ESD banner mentions (cli_banner.py) come from the full ESD analysis
  — do not estimate them from vertical energies without saying so.
- Deviations from the workflow defaults need a reason in the same
  sentence (e.g. more roots than the default because the visible band
  starts at root 8).
- The companion skill tddft-excited-states covers the dashboard/TD-DFT
  mechanics; this skill covers the _esd workflow specifics. Use both
  when the task is a full spectrum.

## Run `/check` before you report the spectrum

`/check` (result_critic.py, docs/result_critic.md) reads the TDDFT
table back and flags: fewer roots than requested (truncated window),
per-state spin contamination (⟨S²⟩ deviation — only read where the
table header names the column), triplet states inside a singlet
spectrum, and a missing broadening width. Report every ⚠️ with the
user verbatim; never report over a ❌ without saying why. Report-only.
