# pKa module validation — SLURM run 2026-10-03

Isodesmic proton-transfer cycles (delfin/pka.py), reference acid acetic
(exp. pKa 4.756). Ten ORCA Opt+Freq jobs (SLURM 7435268-7435277,
operator-submitted), each species B3LYP/def2-SVP CPCM(water),
%pal nprocs 4, %maxcore 1500, <= 2 h walltime per job.
Structures from DELFIN's SMILES converter; Gibbs energies read with
delfin.energies.find_gibbs_energy via pka.read_cycle_gibbs_energies;
pKa anchored with pka.pka_from_cycle. No 1 atm -> 1 M standard-state
correction is applied, and none is needed: in the isodesmic reaction
HA + Ref- -> A- + RefH both sides carry two solutes, so the correction
cancels between target and reference (pka_from_cycle, pka.py:173-188,
computes only the deprotonation-Gibbs difference).

## Job status

All ten ORCA outputs end with "ORCA TERMINATED NORMALLY".
The parser's "failed (exit code 743526x)" outcome reads the SLURM job id
as an exit code: the run leaves a marker file named
`.exit_code_<jobid>` (content: the real code, here 0), and the reader
parses the job id out of the FILE NAME instead of the content
(delfin/doc_server/calc_indexer.py:176-179 names the glob and the
split). The ORCA outputs themselves are
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
(-2.50) lie between. The deviations are systematically negative,
i.e. the computed pKa values come out too LOW: the isodesmic
dG_deprot(target) - dG_deprot(acetic) is too NEGATIVE, so the target
anions come out too STABLE relative to acetate.
Interpretation, not a tested cause: candidates are CPCM missing
specific solvation of the anions (hydrogen bonding to water is not
in a continuum model), delocalised anion electronic structure that
the small basis handles unevenly, the single conformer per species,
and the missing dispersion correction. Which of them dominates is
not decided by this run.

Known limitations of this validation, stated openly:
- def2-SVP without dispersion correction, single conformer per species
  (SMILES converter geometry, no conformer search).
- ORCA's thermochemistry is at the 298.15 K / 1 atm standard state;
  no 1 atm -> 1 M correction is applied because the isodesmic reaction
  has two solutes on each side and the term cancels (see above).
- Experimental pKa values are the module's KNOWN_ACIDS constants
  (delfin/pka.py:71-99), quoted from standard tables at 298.15 K.
