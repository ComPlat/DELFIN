# pKa module validation — SLURM run 2026-10-03

Isodesmic proton-transfer cycles (delfin/pka.py), reference acid acetic
(exp. pKa 4.756). Ten ORCA Opt+Freq jobs (SLURM 7435268-7435277,
operator-submitted), each species B3LYP/def2-SVP CPCM(water),
%pal nprocs 4, %maxcore 1500, <= 2 h walltime per job.
Structures from DELFIN's SMILES converter; Gibbs energies read with
delfin.energies.find_gibbs_energy via pka.read_cycle_gibbs_energies;
pKa anchored with pka.pka_from_cycle (1 atm -> 1 M correction included).

## Job status

All ten ORCA outputs end with "ORCA TERMINATED NORMALLY".
The parser's "failed (exit code 743526x)" outcome reads the SLURM job id
from the wrapper log as an exit code — the ORCA outputs themselves are
clean; no imaginary frequencies in any of the ten outputs.

## Gibbs energies (Hartree, CPCM(water))

| Species         | G (Eh)      | SLURM  |
|-----------------|-------------|--------|
| acetic_HA       | -228.77388592 | 7435269 |
| acetic_A        | -228.30084208 | 7435268 |
| formic_HA       | -189.52987223 | 7435275 |
| formic_A        | -189.06317920 | 7435274 |
| cyanoacetic_HA  | -320.90076837 | 7435273 |
| cyanoacetic_A   | -320.44399750 | 7435272 |
| benzoic_HA      | -420.21258742 | 7435271 |
| benzoic_A       | -419.74619293 | 7435270 |
| phenol_HA       | -306.99896471 | 7435277 |
| phenol_A        | -306.51777778 | 7435276 |

## Deviation table

| Acid         | computed pKa | experimental pKa | delta |
|--------------|--------------|------------------|-------|
| acetic (ref) | 4.76 (anchored) | 4.756 | 0.00 |
| formic       | 1.83         | 3.75             | -1.91 |
| cyanoacetic  | -2.73        | 2.45             | -5.18 |
| benzoic      | 1.70         | 4.20             | -2.50 |
| phenol       | 8.50         | 9.99             | -1.49 |

MAE = 2.77 pKa units (acetic as reference excluded from the mean).

## Interpretation (labeled as such)

The isodesmic cycle cancels the systematic errors of the functional and
the basis set only when proton-transfer reactions between HA/A pairs of
SIMILAR chemistry are involved. The large cyanoacetic deviation (-5.18)
is the expected consequence of a reference (acetic) that does not
resemble the target: the deprotonation Gibbs energy of the nitrile-
substituted acid contains strong through-inductive effects that an
acetic-anchored cycle cannot cancel. Phenol (-1.49) suffers the same
reference mismatch (aryl-OH vs alkyl-OH). Formic (-1.91) and benzoic
(-2.50) lie between. The deviations are systematically negative, i.e.
the computed deprotonation Gibbs energies are too positive relative to
experiment — consistent with CPCM missing specific solvation of the
deprotonated anions (hydrogen bonding to water is not in a continuum
model).

Known limitations of this validation, stated openly:
- def2-SVP without dispersion correction, single conformer per species
  (SMILES converter geometry, no conformer search).
- The thermochemistry uses ORCA's default 298.15 K / 1 atm standard
  state; the 1 atm -> 1 M correction (-1.89 kcal/mol at 298.15 K) is
  applied inside pka_from_cycle as documented in the module docstring.
- Experimental pKa values are the module's KNOWN_ACIDS constants
  (delfin/pka.py:71-99), quoted from standard tables at 298.15 K.
