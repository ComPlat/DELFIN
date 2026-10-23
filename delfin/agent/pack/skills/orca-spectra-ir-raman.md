---
name: orca-spectra-ir-raman
description: >-
  Compute and report IR (and Raman) spectra in DELFIN/ORCA: OPT+FREQ,
  scaling factors with their source, IR activity and point group.
  Use when the user asks for an IR spectrum, vibrational frequencies
  to compare with experiment, or Raman activities.
---

# IR / Raman spectra through DELFIN

## What computes what

IR intensities are computed automatically in every ORCA frequency run;
nothing needs to be switched on (ORCA manual §5.18.1 — "IR spectral
intensities are calculated automatically in frequency runs. Thus, there
is nothing to control by the user."). In DELFIN the frequencies come
from the classic opt+freq workflow (CONTROL `method=classic`, template
`delfin/tools/templates/_classic.py`) or a FREQ step inside a larger
workflow. Run the optimisation first; report frequencies of the
optimised geometry only.

## Workflow

1. **Optimise + frequency** through DELFIN's classic workflow. Do not
   run a bare Hessian on an unoptimised geometry.
2. **Check for imaginary frequencies.** Report `n_imag`. A minimum has
   none; a TS has exactly one. Never silently drop an imaginary mode.
3. **Read the modes and IR intensities** from the output (the typed
   tool `extract_vibrational_modes` gives mode index, frequency in
   cm-1, IR intensity in km/mol, and `is_imaginary`).
4. **Scale before comparing to experiment** (see below). Report both:
   the scaled value AND the factor with its source.
5. **Assign the modes**: highest-frequency X-H stretches, the
   characteristic double-bond stretch (C=O etc.), bending region.
   State the point group and which modes are IR-active (for C2v all
   three fundamentals of H2CO are IR-active; for centrosymmetric
   molecules the mutual-exclusion rule applies — g modes are Raman
   only, u modes IR).

## Scaling factors — with a source, always

Harmonic frequencies overestimate fundamentals; scale before comparing
to experiment. The factor depends on functional AND basis set.

- Source of the mechanism: ORCA manual §4.6.2 — edit
  `$frequency_scale_factor` in the `.hess` file and rerun `orca_vib`.
- Source of the numbers: NIST CCCBDB vibrational scaling factors
  (cccbdb.nist.gov/vibscalex.asp) and the Truhlar group database
  (comp.chem.umn.edu/freqscale); for double-hybrids and DFT:
  Kesharwani, Brauer & Martin, *J. Phys. Chem. A* 2015, 119, 1701–1714,
  DOI 10.1021/jp508422u.
- Typical order of magnitude (verify the exact pair against the
  databases above before quoting): B3LYP with a double-zeta basis
  ≈ 0.957–0.961; BP86 similar; TPSSh/hybrid-meta slightly higher;
  HF ≈ 0.90. A ZPE factor is a DIFFERENT number (~0.97–0.99) — do not
  use one for the other.

## Deviations need a reason

A deviation from DELFIN's workflow defaults (method, basis, solvation)
is only acceptable with a stated reason in the same sentence as the
deviation. Using an unscaled frequency in a comparison with experiment
is a deviation: name it ("raw harmonic, factor 1.0, because the
comparison is with another harmonic calculation").

## Honesty rules

- Every frequency carries cm-1; every intensity km/mol (IR) or Å⁴ (Raman).
- The scaled number is labelled scaled; the raw number raw.
- The scaling factor is quoted WITH its method/basis pair and its source.
- Imaginary modes are reported, never hidden (see the
  freq-thermochemistry skill for thermochemistry after a clean Hessian).
