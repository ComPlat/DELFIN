# Phase E capability tests — Haiku (2026-04-28)

5 new behavioral tests for the Phase E parsers (SCF convergence,
Mulliken charges, vibrational modes, DELFIN_data.json, multi-property
summary). Single sweep on Haiku.

**Result: 5/5 PASS in ~$0.74.**

## Per-cell verdicts

| Test | Verdict | Cost |
|---|---|---|
| SCF convergence question → `extract_scf_convergence` | ✓ tool | $0.109 |
| Mulliken charges → `extract_mulliken_charges` | ✓ tool | $0.184 |
| vibrational modes / IR → `extract_vibrational_modes` | ✓ tool | $0.164 |
| DELFIN_data.json question → `extract_delfin_json` | ✓ tool | $0.105 |
| multi-property summary → `extract_calc_summary_table` | ✓ named | $0.179 |

## What this proves

The dashboard agent's decision-tree extension (commit `e898538`) for
Phase E intents lands on the typed MCP tool 100% of the time on
Haiku. No iteration needed — the agent picks the right tool on
first try.

The 6 Phase E parsers + 3 plots add 9 new typed MCP tools to the
ops surface (total 65 → 74). All wired into the inline-emit hook
so plots auto-display in the dashboard chat.

## Tool-by-tool coverage now

After Phase E, the dashboard agent has end-to-end typed-tool
support for every common ORCA-output question:

- imag freq / minimum / TS → `extract_imaginary_frequencies`
- HOMO/LUMO / gap / orbitals → `extract_orbital_energies`
- UV/Vis / TDDFT → `extract_excited_states` + `plot_uvvis_spectrum`
- dipole moment → `extract_dipole`
- opt steps / convergence → `extract_optimization_trajectory` +
  `plot_optimization_convergence`
- **SCF iteration history → `extract_scf_convergence` +
  `plot_scf_convergence`** (new)
- Gibbs/SPE/ZPE one folder → `parse_orca_output`
- Gibbs/SPE many folders → `extract_energy_table`
- lowest/highest property → `find_calculation_extreme`
- Thermochem (T, P, H, S, G) → `extract_thermochem`
- ORCA errors / SCF / OOM → `find_orca_errors`
- compare two calcs → `compare_calculations`
- compare across functionals → `compare_across_functionals`
- **Mulliken/Loewdin charges → `extract_mulliken_charges` /
  `extract_loewdin_charges` + `plot_population_charges`** (new)
- **all vib modes + IR → `extract_vibrational_modes` +
  `plot_vibrational_spectrum`** (new)
- **DELFIN_data.json → `extract_delfin_json`** (new)
- **multi-property summary → `extract_calc_summary_table`** (new)

## How to run

```bash
DELFIN_AGENT_DRYRUN=1 python -m tests.agent_dryrun.report \\
  --models haiku \\
  --include "SCF convergence,Mulliken,vibrational modes,DELFIN_data,multi-property" \\
  --out /tmp/phaseE.md
```
