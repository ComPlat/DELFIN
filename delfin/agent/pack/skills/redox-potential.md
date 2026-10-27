---
name: redox-potential
description: >-
  Compute redox potentials through DELFIN's redox workflow: the
  thermodynamic cycle, the reference electrode (E_ref), the M1/M2/M3
  methods, and what makes a reported potential scientifically valid.
  Use when the user asks for an oxidation, reduction, or redox potential.
domains: code
---

# Redox potentials through DELFIN's redox workflow

## Route

Redox potentials run through DELFIN's redox workflow (template
`delfin/tools/templates/_redox.py`, `redox_potential`; CLI
`delfin-app template redox_potential` writes the keys file,
cli_app.py:7). What the template encodes:

- **Trunk**: SMILES → xTB pre-opt → ORCA opt + freq for the neutral
  species (template steps `smiles_to_xyz`, `xtb_opt`, `orca_opt`,
  `orca_freq`, labels neutral_opt/neutral_freq).
- **Branches**: oxidation (+1 charge, multiplicity mult_ox) and
  reduction (-1 charge, multiplicity mult_red), each opt + freq.
- DELFIN collects the Gibbs energies of the redox states
  (delfin/energies.py:11-32, canonical 7-state map) and computes
  potentials by three methods (delfin/cli_calculations.py:8,
  calculate_redox_potentials → M1 average, M2 stepwise, M3 mixed;
  final selection priority 3 > 2 > 1, cli_calculations.py:98).

## The thermodynamic cycle (state it, always)

For a one-electron oxidation Ox/Red couple, the potential comes from
the reaction free energy ΔG°_ox = G°(Ox) − G°(Red) (including ZPE and
thermal corrections from the frequency runs; solvation model identical
on both sides) via

    E° = −ΔG°_ox / (n·F)     with n = 1, F = 23.0605 kcal/(mol·V)

referenced against the electrode potential E_ref (DELFIN's
`get_E_ref`, delfin/config.py:1027-1084). Deviations from the cycle
(an isodesmic correction, a different solvation model on one side)
need a stated reason — and a mixed solvation model across the couple
is never acceptable.

## The reference electrode (never omit it)

A potential without its reference electrode is not a result. DELFIN's
per-solvent E_ref values vs SHE (config.py:1074-1084): MeCN 4.745 V,
DMF 4.795 V, DCM 4.805 V, THF 4.905 V, DMSO 4.780 V, acetone 4.825 V,
DME 4.855 V, default 4.345 V (gas-phase absolute scale). ONIOM
(QM/QM2) runs use the adjusted −3.31 V branch (config.py:1052-1071) —
an ONIOM potential is on a different reference and must be labelled.

Report E with its electrode explicitly: "E_ox = +0.55 V vs SHE
(E_ref = 4.745 V, MeCN)". To convert to Fc/Fc+, subtract the Fc/Fc+
potential in that solvent — state the conversion value and its source.

## Honesty rules

- Same level of theory on BOTH sides of the couple: same functional,
  basis, dispersion, solvation model, geometry provenance. A potential
  computed across mixed levels is meaningless — refuse it and name
  each file's level instead (that is the wave-20 comparability rule).
- Multiplicities of the oxidised and reduced states must be stated
  (an electron removed changes the spin state).
- Report unconverged runs and imaginary frequencies; a Gibbs energy
  from a geometry that is not a minimum invalidates the cycle.
- The three DELFIN methods (M1/M2/M3) may disagree; report the
  selected one AND the spread, not only the winner.
- Deviations from the workflow defaults need a reason in the same
  sentence as the deviation.
