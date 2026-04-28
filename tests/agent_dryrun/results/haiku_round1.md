# Haiku 34-case dry-run sweep — Round 1 (2026-04-28)

First end-to-end behavioral sweep of every typed MCP tool the user
asked to verify. Run on `claude --print --permission-mode plan` with
the dashboard agent prompt that included the **first** version of the
parsing-first decision tree (commit `0938d37`). The decision tree
covered parsing tools but did NOT yet list workflow/lifecycle tools
(validate_orca_input, submit_calculation, list_ssh_transfer_jobs,
list_dashboard_widgets, …). Those got added in commit `39d4d5e`
afterwards — Round 2 reruns on the extended prompt.

Total: **20/34 PASS** ($3.025 spent, ~5 min/cell median).

## Per-cell verdicts

| # | Test | Verdict | Cost |
|---|---|---|---|
| 1 | imag-freq → extract_imaginary_frequencies | ✓ tool | $0.102 |
| 2 | functional comparison → compare_across_functionals | ✓ tool | $0.189 |
| 3 | rename calc folder → rename_calc_folder | ✓ named | $0.048 |
| 4 | archive intent → move_to_archive | ✓ named | $0.051 |
| 5 | cancel-all → kill_all_user_jobs | ✗ TIMED OUT | $0.000 |
| 6 | smart-recalc prep → prepare_recalc | ✗ TIMED OUT | $0.000 |
| 7 | HOMO/LUMO → extract_orbital_energies | ✓ tool | $0.137 |
| 8 | UV/Vis → plot_uvvis_spectrum / extract_excited_states | ✓ tool | $0.150 |
| 9 | opt convergence → extract_optimization_trajectory | ✓ tool | $0.076 |
| 10 | batch UX safety → no 'initial.xyz' folder filter | ✗ TIMED OUT | $0.000 |
| 11 | catalog discovery → list_tools / describe_tool | ✓ tool | $0.066 |
| 12 | thermochem → extract_thermochem | ✓ tool | $0.116 |
| 13 | dipole → extract_dipole | ✓ tool | $0.093 |
| 14 | multi-folder energy table → extract_energy_table | ✓ tool | $0.078 |
| 15 | ORCA error scan → find_orca_errors | ✓ tool | $0.124 |
| 16 | lowest-Gibbs → find_calculation_extreme | ✓ tool | $0.111 |
| 17 | two-folder diff → compare_calculations | ✓ tool | $0.209 |
| 18 | INP validation → validate_orca_input | ✗ missed | $0.075 |
| 19 | submit one folder → submit_calculation | ✗ missed | $0.073 |
| 20 | list active jobs → list_active_calculations | ✓ tool | $0.120 |
| 21 | ssh transfer status → list_ssh_transfer_jobs | ✗ missed | $0.056 |
| 22 | ask for batch recipe → get_dashboard_pattern | ✓ tool | $0.071 |
| 23 | widget catalog → list_dashboard_widgets | ✗ TIMED OUT | $0.000 |
| 24 | DELFIN concept → explain_delfin_feature | ✗ TIMED OUT | $0.000 |
| 25 | describe one tool → describe_tool | ✗ missed | $0.172 |
| 26 | ORCA-syntax → check_orca_manual_indexed | ✓ tool | $0.172 |
| 27 | ad-hoc PDF read → read_pdf | ✗ missed | $0.054 |
| 28 | options-menu lookup → list_calc_options | ✓ tool | $0.108 |
| 29 | create new folder → create_calc_folder | ✗ missed | $0.051 |
| 30 | move folder → move_calc_folder | ✗ TIMED OUT | $0.000 |
| 31 | delete folder → delete_calc_folder | ✓ named | $0.073 |
| 32 | energy distribution → plot_energy_distribution | ✗ TIMED OUT | $0.000 |
| 33 | energy correlation → plot_energy_correlation | ✗ missed | $0.290 |
| 34 | orbital diagram → plot_orbital_diagram | ✓ named | $0.162 |

## Failure breakdown

**Timeouts (7/14)** — 90 s per-cell budget too tight for multi-step
workflow reasoning:
- cancel-all (5), smart-recalc (6), batch UX (10), widget catalog (23),
  DELFIN concept (24), move folder (30), energy distribution (32).
- **Fix:** bumped per-cell timeout to 150 s (commit `5399725`).

**Misses (7/14)** — agent didn't pick the typed tool nor mention it:
- validate_orca_input (18), submit_calculation (19),
  list_ssh_transfer_jobs (21), describe_tool (25), read_pdf (27),
  create_calc_folder (29), plot_energy_correlation (33).
- **Fix:** extended decision tree in dashboard prompt to list
  workflow/lifecycle tools (commit `39d4d5e`). Includes:
  validate_orca_input, submit a folder, list_active_calculations,
  list_ssh_transfer_jobs, get_dashboard_pattern.

## Big wins (Round 0 → Round 1)

The first standalone test — before the parsing-first prompt change —
showed haiku falling back to Glob/Grep on `.out` files for every
parsing query. Round 1 (with the strengthened prompt) shows haiku
correctly picking the typed tool for **every** parsing question:

- imag-freq (1) ✓
- HOMO/LUMO (7) ✓
- UV/Vis (8) ✓
- opt convergence (9) ✓
- thermochem (12) ✓
- dipole (13) ✓
- energy table (14) ✓
- ORCA errors (15) ✓
- lowest-Gibbs (16) ✓
- two-folder diff (17) ✓

The hard "Typed MCP tool BEFORE Glob/Grep" rule + the per-question
decision tree closed the gap that the very first behavioral test
exposed.

Round 2 sweep (`--rerun-failed`) covers the 14 Round-1 failures with
the extended prompt + 150 s per-cell timeout — see haiku_round2.md
when complete.
