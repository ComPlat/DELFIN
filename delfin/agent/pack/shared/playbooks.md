# DELFIN Playbooks

Step-by-step strategies for common task types. Follow these when working on
the corresponding module — they encode validated approaches from prior sessions.

## build_up_complex.py (1531 lines)

1. Grep for the target function name
2. Read only the function body (offset+limit)
3. Check `_pso_fitness`, `_compute_clash_score`, `_prepare_force_field_terms` for side effects
4. Edit, then run: `pytest tests/ -x -q -k 'not browser'`
5. If touching VdW/clash logic: verify radii are VdW (not covalent), threshold 75% of sum

**Invariants:** PSO places ligands as rigid bodies (6 DOF) — never modify
intra-ligand topology. VdW radii for clash detection, NOT covalent. Procrustes
init for bidentate+ ligands (>=2 donors → SVD alignment).

## config.py (834 lines)

1. Grep for the relevant key/function
2. Check `_parse_control_file` (~line 353) for parsing flow
3. Check `validate_control_text` (~line 501) for validation rules
4. Run: `pytest tests/test_co2_control_overrides.py tests/test_functional_contracts.py -x -q`

**Invariants:** Config dicts are mutable — changes propagate. OCCUPIER_parser
is separate from main CONTROL parser.

## cli.py (2034 lines)

1. Grep for the subcommand name — each is a `_run_*_subcommand` function
2. Check `_load_full_cli_dependencies` (line 24) for lazy imports
3. Run: `pytest tests/test_cli_*.py -x -q`

**Invariants:** CLI uses lazy loading. Override logic:
`_parse_occupier_overrides`, `_apply_occupier_overrides`. Downstream cascade:
`_downstream_stages` invalidates dependent stages.

## orca_recovery.py (1728 lines)

1. `OrcaErrorType` enum (line 28) — all known error types
2. `OrcaErrorDetector` (line 62) — matches output patterns
3. `RecoveryStrategy` (line 275) — maps errors to fixes
4. `OrcaInputModifier` (line 684) — applies input changes
5. Run: `pytest tests/test_orca_workflow_contracts.py -x -q`

**Invariants:** `RetryStateTracker` (line 1633) prevents infinite loops.
Recovery must work for both local and SLURM. Never modify the original input.

## smiles_converter.py (18073 lines) — LARGEST MODULE

1. ALWAYS Grep first — never Read the whole file
2. Main classes: `_HybridHaptoFragment`, `_PrimaryOrganometalModule`
3. Entry point: `_try_multiple_strategies` (~line 627)
4. Metal handling: `_manual_metal_embed` (~line 835)
5. Hapto groups: `_find_hapto_groups` (~line 1944)

**Invariants:** Multiple fallback strategies (RDKit → OpenBabel → manual embed).
Metal bonds → dative bonds via `_convert_metal_bonds_to_dative`. Hapto
approximation is optional (`_hapto_approx_enabled` flag).

## dashboard/tab_agent.py (7435 lines)

1. ALWAYS Grep first
2. Entry point: `create_tab` (line 931)
3. Grep for specific widget or callback names
4. Run: `pytest tests/test_agent_*.py -x -q`
