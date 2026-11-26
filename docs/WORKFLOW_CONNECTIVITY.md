# WORKFLOW CONNECTIVITY — what a CONTROL file switches on, and how it is wired

<!-- GENERATED FILE - do not edit by hand.
     Source: delfin/workflow_connectivity.py (CONTROL_MODULES) plus
     delfin/define.py (TEMPLATE) and delfin/config.py (the alias table).
     Regenerate: python -m delfin.workflow_connectivity
     Pinned by: tests/test_the_control_connectivity_map_is_current.py -->

A CONTROL file is the only place a DELFIN run is configured. This map
answers, for every module one can switch on: which key switches it, the
code path from that key to the module, what the module is handed from the
step before it, and where its result goes. It is generated from the code
so it cannot drift: a renamed entry point or an unlisted template key
fails the test named above.

`delfin/define.TEMPLATE` writes **186 keys**. All of them are
central: the last table lists every one with the modules that read it.

## 1. Modules a CONTROL file switches on

| Module | Switch key | Values | Entry point the pipeline calls |
|---|---|---|---|
| SMILES / structure construction | `smiles_converter` | `QUICK`, `NORMAL`, `MANTA`, `ARCHITECTOR`, `MOLSIMPLIFY`, `MACE`, `RXEMBED`, `MOLASSEMBLER`, `METALLOGEN` | `delfin.workflows.pipeline:normalize_input_file` |
| global geometry optimisation (xTB/CREST) | `global_optimizer` | `GOAT`, `CREST` | `delfin.workflows.pipeline:run_classic_phase` |
| IMAG (imaginary-mode elimination) | `IMAG` | `yes` | `delfin.imag:run_IMAG` |
| classic redox phase | `method` | `classic` | `delfin.workflows.pipeline:run_classic_phase` |
| manually-written spin states | `method` | `manually` | `delfin.workflows.pipeline:run_manual_phase` |
| OCCUPIER (preferred spin state) | `method` | `OCCUPIER` | `delfin.workflows.pipeline:run_occuper_phase` |
| ESD (excited-state dynamics) | `ESD_modul` | `yes` | `delfin.esd_module:execute_esd_jobs` |
| hyperpolarisability (sTD-DFT-xTB) | `hyperpol_xTB` | `yes` | `delfin.workflows.pipeline:run_hyperpol_xtb_phase` |
| TADF xTB screening | `tadf_xTB` | `yes` | `delfin.workflows.pipeline:run_tadf_xtb_phase` |
| CO2 coordination | `co2_coordination` | `on`, `yes`, `true`, `1` | `delfin.cli:_run_co2_subcommand` |
| thermodynamics / stability constants | `stability_constant` | `yes` | `delfin.stability_constant:run_stability_constant_phase` |

## 2. What crosses each step boundary

| Module | Handed by the step before | Its result goes to |
|---|---|---|
| SMILES / structure construction | the SMILES line in input.txt (or CONTROL's SMILES key) | start.txt -- DELFIN coordinate format, the geometry every later step reads |
| global geometry optimisation (xTB/CREST) | start.txt | the optimised geometry written back to start.txt; XTB_OPT and XTB_preOPT select the pre-optimisation instead |
| IMAG (imaginary-mode elimination) | the step's .out/.inp and its .hess, in the step's own folder | a relaxed geometry copied back to the calling step (and, for OCCUPIER, to input_<step>_OCCUPIER.xyz) |
| classic redox phase | start.txt / initial.xyz, charge and multiplicity from CONTROL | initial/ox_step_N/red_step_N folders, each with its .inp/.out; the Gibbs energies the report is built from |
| manually-written spin states | multiplicity_N / BrokenSym_N per step from CONTROL | the same step folders as classic, using the stated spin states |
| OCCUPIER (preferred spin state) | start.txt, charge, and the OCCUPIER_sequence_profiles block | one *_OCCUPIER folder per step and OCCUPIER.txt, which names the preferred configuration every later step and the report use |
| ESD (excited-state dynamics) | start.txt / initial.xyz (the geometry the previous step left in the workspace), charge, solvent, metals and both basis sets | ESD/ with one input and output per state (S0, S1, T1, ...), plus ISC/IC/fluorescence/phosphorescence outputs the reports read |
| hyperpolarisability (sTD-DFT-xTB) | hyperpol_xTB_xyz (default start.txt), charge and multiplicity, in the CONTROL file's directory | hyperpol_xtb/ with the response run and hyperpol_xtb_summary.json |
| TADF xTB screening | tadf_xTB_xyz (default start.txt), charge and multiplicity | tadf_xtb/ with the singlet/triplet screening the TADF report reads |
| CO2 coordination | the finished run's preferred structures and co2_species_delta | CO2_coordination/ (its own workspace, entered with os.chdir), written back into the report |
| thermodynamics / stability constants | the initial step's Gibbs energy and FREQ output; embedded into the OCCUPIER scheduler when that method runs | stability_constant/ jobs and the logK the report reads |

## 3. Notes on individual modules

- **SMILES / structure construction** — MANTA is the current name for the GUPPY sampler; a legacy smiles_converter=GUPPY or a bare GUPPY=yes is read as MANTA (_apply_guppy_legacy, _resolve_smiles_converter).
- **global geometry optimisation (xTB/CREST)** — an empty global_optimizer means no global optimisation, which is what the two keys it replaced (XTB_GOAT=no, CREST=no) said.
- **IMAG (imaginary-mode elimination)** — a key read but with no step to act on is reported, not silently ignored (delfin/imag.py:_refuse_rates_on_a_saddle family).
- **ESD (excited-state dynamics)** — reached in parallel inside the classic scheduler when method is classic (add_esd_jobs_to_scheduler), and as its own phase otherwise (run_esd_phase).  It reads the workspace, not the process directory -- see tests/test_a_control_run_writes_into_its_workspace.py.
- **thermodynamics / stability constants** — the template spells the switch thermodynamics=, which the alias table maps onto stability_constant; thermodynamics_mode and thermodynamics_reaction are aliases too.

## 4. Which modules the agent's tools reach

The map above says what a CONTROL *run* reaches. An agent reaches a
module by a different route, and the two are not the same count:
`delfin/tool_reach.py` lists, per module, the registered building blocks
(`delfin.tools._registry`) and the agent tools (`delfin/api.py`, the
tool catalogue). Both are checked against the code by
`tests/test_the_tool_reach_names_real_routes.py`, so a renamed tool shows
up as a missing route rather than an inflated figure.

| Module | Building blocks | Agent tools | Reachable by a tool |
|---|---|---|---|
| ESD (excited-state dynamics) | `esd`, `esd_opi`, `orca_tddft` | `pipeline_run` | yes |
| IMAG (imaginary-mode elimination) | `imag_fix`, `imag_opi` | — | yes |
| OCCUPIER (preferred spin state) | `occupier`, `occupier_opi` | — | yes |
| SMILES / structure construction | `smiles_to_xyz`, `obabel_from_smiles`, `guppy_sampling` | `geometry_build_smiles` | yes |
| global geometry optimisation (xTB/CREST) | `crest_conformers`, `xtb_goat`, `orca_goat_opi` | — | yes |
| hyperpolarisability (sTD-DFT-xTB) | — | `hyperpol` | yes |
| TADF xTB screening | — | `tadf_xtb` | yes |
| CO2 coordination | `co2_coordinator_opi` | `co2` | yes |
| thermodynamics / stability constants | — | — | **no — CONTROL run only** |
| classic redox phase | `orca_opt_freq`, `orca_sp`, `orca_freq` | `pipeline_run`, `run_orca_input`, `submit_calculation` | yes |

Reachable by a CONTROL run but by **no tool at all**: `thermodynamics / stability constants`.

## 5. Every template key, and who reads it

`by alias` means the key reaches its reader under another name
(`delfin/config.py:_CONTROL_KEY_ALIASES`); a scan for the literal name
would miss it, which is how two of these were once reported as unused.

| Key | Read as | Modules that read it |
|---|---|---|
| `input_file` | own name | `delfin/agent/api_client.py`, `delfin/agent/calc_policy.py`, `delfin/analysis_tools/conformer_fragment_recognition.py`, `delfin/analysis_tools/conformer_shape_classification.py` (+34 more) |
| `NAME` | own name | `delfin/agent/hooks.py`, `delfin/cli.py`, `delfin/cli_report.py`, `delfin/common/control_allowed_values.py` (+4 more) |
| `SMILES` | own name | `delfin/agent/api_client.py`, `delfin/agent/calc_policy.py`, `delfin/cli_manta.py`, `delfin/common/control_allowed_values.py` (+8 more) |
| `charge` | own name | `delfin/agent/briefing.py`, `delfin/agent/calc_policy.py`, `delfin/agent/chat_plots.py`, `delfin/agent/cli_new_workflow.py` (+100 more) |
| `implicit_solvation_model` | own name | `delfin/cli.py`, `delfin/cli_fukui.py`, `delfin/co2/CO2_Coordinator6.py`, `delfin/co2/chain_setup.py` (+14 more) |
| `solvent` | own name | `delfin/agent/api_client.py`, `delfin/agent/briefing.py`, `delfin/agent/engine.py`, `delfin/agent/pack/benchmark/setup/a_small_calc_archive.py` (+59 more) |
| `XTB_SOLVATOR` | own name | `delfin/cli.py`, `delfin/cli_recalc.py`, `delfin/recalc_control.py`, `delfin/workflows/engine/occupier.py` (+2 more) |
| `number_explicit_solv_molecules` | own name | `delfin/cli.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/recalc_control.py` (+1 more) |
| `smiles_converter` | own name | `delfin/cli_fukui.py`, `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/dashboard/tab_submit.py` (+5 more) |
| `xTB_method` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/dashboard/goat_method.py`, `delfin/dashboard/tab_help.py` (+7 more) |
| `XTB_preOPT` | read as `XTB_OPT` (alias) | `delfin/cli.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py` (+4 more) |
| `global_optimizer` | own name | `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/dashboard/constants.py`, `delfin/dashboard/goat_method.py` (+5 more) |
| `multiplicity_global_opt` | own name | `delfin/cli.py`, `delfin/cli_fukui.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py` (+6 more) |
| `IMAG` | own name | `delfin/cli_imag.py`, `delfin/dashboard/tab_help.py`, `delfin/doc_server/calc_indexer.py`, `delfin/doc_server/calc_search.py` (+9 more) |
| `IMAG_scope` | own name | `delfin/cli_imag.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/imag.py` (+3 more) |
| `IMAG_option` | own name | `delfin/common/control_validator.py`, `delfin/tools/workflows.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/engine/occupier.py` |
| `allow_imaginary_freq` | own name | `delfin/cli_imag.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/imag.py` (+4 more) |
| `IMAG_sp_energy_window` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/imag.py`, `delfin/tools/workflows.py` (+1 more) |
| `IMAG_optimize_candidates` | own name | `delfin/imag.py`, `delfin/workflow_connectivity.py` |
| `IMAG_max_rounds` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/imag.py`, `delfin/tools/workflows.py` (+2 more) |
| `method` | own name | `delfin/agent/api_client.py`, `delfin/agent/briefing.py`, `delfin/agent/calc_budget.py`, `delfin/agent/calc_policy.py` (+90 more) |
| `calc_initial` | own name | `delfin/cli.py`, `delfin/cli_imag.py`, `delfin/cli_report.py`, `delfin/common/control_validator.py` (+7 more) |
| `oxidation_steps` | own name | `delfin/cli.py`, `delfin/cli_imag.py`, `delfin/cli_report.py`, `delfin/common/control_validator.py` (+7 more) |
| `reduction_steps` | own name | `delfin/cli.py`, `delfin/cli_imag.py`, `delfin/cli_report.py`, `delfin/common/control_validator.py` (+7 more) |
| `calc_potential_method` | own name | `delfin/cli_calculations.py`, `delfin/common/control_validator.py`, `delfin/recalc_control.py`, `delfin/tools/workflows.py` |
| `E_ref` | own name | `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/recalc_control.py`, `delfin/reporting/delfin_collector.py` (+2 more) |
| `multiplicity_0` | own name | `delfin/agent/calc_policy.py`, `delfin/cli.py`, `delfin/cli_fukui.py`, `delfin/cli_report.py` (+8 more) |
| `BrokenSym_0` | own name | `delfin/cli_imag.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py` (+1 more) |
| `multiplicity_ox1` | own name | `delfin/cli.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py` |
| `BrokenSym_ox1` | own name | `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py` |
| `multiplicity_ox2` | own name | `delfin/cli.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `BrokenSym_ox2` | own name | `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `multiplicity_ox3` | own name | `delfin/cli.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `BrokenSym_ox3` | own name | `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `multiplicity_red1` | own name | `delfin/cli.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `BrokenSym_red1` | own name | `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `multiplicity_red2` | own name | `delfin/cli.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `BrokenSym_red2` | own name | `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `multiplicity_red3` | own name | `delfin/cli.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `BrokenSym_red3` | own name | `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` |
| `OCCUPIER_method` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/copy_helpers.py`, `delfin/occupier.py` (+6 more) |
| `OWN_TREE_PURE_WINDOW` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/occupier_auto.py`, `delfin/workflow_connectivity.py` |
| `OWN_progressive_from` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/occupier_auto.py`, `delfin/workflow_connectivity.py` |
| `fob_equal_weights` | own name | `delfin/common/control_validator.py`, `delfin/occupier_flat_extraction.py`, `delfin/workflow_connectivity.py` |
| `OCCUPIER_compare` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config_manager.py`, `delfin/tools/workflows.py` (+1 more) |
| `occupier_selection` | own name | `delfin/common/control_validator.py`, `delfin/reporting/occupier_reports.py`, `delfin/reporting/occupier_selection.py`, `delfin/workflow_connectivity.py` |
| `occupier_precision` | own name | `delfin/common/control_validator.py`, `delfin/reporting/occupier_reports.py`, `delfin/reporting/occupier_selection.py`, `delfin/workflow_connectivity.py` |
| `occupier_epsilon` | own name | `delfin/common/control_validator.py`, `delfin/reporting/occupier_reports.py`, `delfin/reporting/occupier_selection.py`, `delfin/workflow_connectivity.py` |
| `clean_override_window_h` | own name | `delfin/common/control_validator.py`, `delfin/reporting/occupier_reports.py`, `delfin/reporting/occupier_selection.py` |
| `clean_quality_improvement` | own name | `delfin/common/control_validator.py`, `delfin/reporting/occupier_reports.py`, `delfin/reporting/occupier_selection.py` |
| `clean_quality_good` | own name | `delfin/common/control_validator.py`, `delfin/reporting/occupier_reports.py`, `delfin/reporting/occupier_selection.py` |
| `maxiter_occupier` | own name | `delfin/cli.py`, `delfin/common/control_validator.py`, `delfin/occupier.py`, `delfin/tools/occupier_opi.py` |
| `geom_opt_OCCUPIER` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/occupier.py`, `delfin/reporting/occupier_reports.py` (+3 more) |
| `pass_wavefunction` | own name | `delfin/common/control_validator.py`, `delfin/occupier.py`, `delfin/occupier_flat_extraction.py`, `delfin/tools/adapters/occupier_opi.py` (+1 more) |
| `approximate_spin_projection_APMethod` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/dashboard/tab_help.py`, `delfin/occupier.py` (+2 more) |
| `even_seq` | own name | `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/copy_helpers.py`, `delfin/occupier.py` (+4 more) |
| `odd_seq` | own name | `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/copy_helpers.py`, `delfin/occupier.py` (+4 more) |
| `keyword:basename` | own name | **none by name** — see the note below |
| `additions:basename` | own name | **none by name** — see the note below |
| `co2_coordination` | own name | `delfin/cli.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py` (+5 more) |
| `co2_species_delta` | own name | `delfin/cli.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py` (+4 more) |
| `calc_prop_of_interest` | own name | `delfin/doc_server/calc_indexer.py`, `delfin/esd_input_generator.py`, `delfin/recalc_control.py`, `delfin/workflow_connectivity.py` (+2 more) |
| `properties_of_interest` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/doc_server/calc_indexer.py`, `delfin/esd_input_generator.py` (+5 more) |
| `reorganisation_energy` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/dashboard/tab_help.py`, `delfin/doc_server/calc_indexer.py` (+5 more) |
| `ESD_modul` | own name | `delfin/cli.py`, `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/doc_server/calc_indexer.py` (+5 more) |
| `ESD_modus` | own name | `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/doc_server/calc_indexer.py`, `delfin/esd_input_generator.py` (+6 more) |
| `ESD_T1_opt` | own name | `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/esd_input_generator.py`, `delfin/tools/workflows.py` (+1 more) |
| `ESD_frequency` | own name | `delfin/esd_input_generator.py`, `delfin/esd_module.py`, `delfin/tools/esd_opi.py`, `delfin/workflow_connectivity.py` |
| `states` | own name | `delfin/agent/api_client.py`, `delfin/agent/chat_plots.py`, `delfin/agent/generate_fact_tasks.py`, `delfin/agent/office.py` (+13 more) |
| `ISCs` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/esd_module.py` (+3 more) |
| `ICs` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/esd_module.py` (+3 more) |
| `emission_rates` | own name | `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/dashboard/tab_help.py`, `delfin/esd_module.py` (+4 more) |
| `phosp_IROOT` | own name | `delfin/esd_input_generator.py`, `delfin/recalc_control.py`, `delfin/tools/esd_opi.py` |
| `phosp_keywords` | own name | `delfin/esd_input_generator.py`, `delfin/recalc_control.py`, `delfin/tools/esd_opi.py` |
| `fluor_keywords` | own name | `delfin/esd_input_generator.py`, `delfin/recalc_control.py`, `delfin/tools/esd_opi.py` |
| `TROOTSSL` | own name | `delfin/esd_input_generator.py`, `delfin/esd_module.py`, `delfin/esd_results.py`, `delfin/recalc_control.py` (+1 more) |
| `addition_S0` | own name | `delfin/esd_input_generator.py`, `delfin/recalc_control.py`, `delfin/tools/esd_opi.py` |
| `DOHT` | own name | `delfin/esd_input_generator.py`, `delfin/recalc_control.py`, `delfin/tools/esd_opi.py` |
| `ESD_LINES` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `ESD_LINEW` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `ESD_INLINEW` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `ESD_NPOINTS` | own name | `delfin/esd_input_generator.py` |
| `ESD_MAXTIME` | own name | `delfin/esd_input_generator.py` |
| `hybrid1_geom_MaxIter` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `elprop_Dipole` | own name | `delfin/esd_input_generator.py` |
| `elprop_Quadrupole` | own name | `delfin/esd_input_generator.py` |
| `elprop_Hyperpol` | own name | `delfin/esd_input_generator.py` |
| `elprop_Polar` | own name | `delfin/esd_input_generator.py` |
| `elprop_PolarVelocity` | own name | `delfin/esd_input_generator.py` |
| `elprop_PolarDipQuad` | own name | `delfin/esd_input_generator.py` |
| `elprop_PolarQuadQuad` | own name | `delfin/esd_input_generator.py` |
| `deltaSCF_DOMOM` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `deltaSCF_PMOM` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `deltaSCF_keepinitialref` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `deltaSCF_SOSCFHESSUP` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `deltaSCF_keywords` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `deltaSCF_maxiter` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `deltaSCF_SOSCFConvFactor` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `deltaSCF_SOSCFMaxStep` | own name | `delfin/esd_input_generator.py`, `delfin/tools/esd_opi.py` |
| `TDDFT_nroots` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/orca_input.py`, `delfin/common/tddft_settings.py` (+4 more) |
| `TDDFT_maxdim` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/orca_input.py`, `delfin/common/tddft_settings.py` (+1 more) |
| `TDDFT_maxiter` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/orca_input.py`, `delfin/common/tddft_settings.py` (+2 more) |
| `TDDFT_TDA` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/orca_input.py`, `delfin/common/tddft_settings.py` (+3 more) |
| `TDDFT_followiroot` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/orca_input.py`, `delfin/common/tddft_settings.py` (+2 more) |
| `TDDFT_SOC` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/orca_input.py`, `delfin/common/tddft_settings.py` (+3 more) |
| `TDDFT_additions` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/tddft_settings.py`, `delfin/config.py` (+1 more) |
| `functional` | own name | `delfin/agent/api_client.py`, `delfin/agent/briefing.py`, `delfin/agent/engine.py`, `delfin/agent/generate_fact_tasks.py` (+49 more) |
| `disp_corr` | own name | `delfin/cli_fukui.py`, `delfin/co2/CO2_Coordinator6.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_validator.py` (+11 more) |
| `ri_jkx` | own name | `delfin/cli_fukui.py`, `delfin/co2/CO2_Coordinator6.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_allowed_values.py` (+11 more) |
| `relativity` | own name | `delfin/agent/workflow_grader.py`, `delfin/cli.py`, `delfin/co2/CO2_Coordinator6.py`, `delfin/common/control_allowed_values.py` (+12 more) |
| `aux_jk` | own name | `delfin/co2/CO2_Coordinator6.py`, `delfin/co2/chain_setup.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py` (+3 more) |
| `aux_jk_rel` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/reporting/occupier_reports.py`, `delfin/utils.py` |
| `main_basisset` | own name | `delfin/agent/api_client.py`, `delfin/cli_fukui.py`, `delfin/cli_report.py`, `delfin/co2/CO2_Coordinator6.py` (+20 more) |
| `main_basisset_rel` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/doc_server/calc_indexer.py`, `delfin/utils.py` |
| `metal_basisset` | own name | `delfin/cli_fukui.py`, `delfin/cli_report.py`, `delfin/co2/CO2_Coordinator6.py`, `delfin/co2/chain_setup.py` (+21 more) |
| `metal_basisset_rel` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/doc_server/calc_indexer.py`, `delfin/reporting/occupier_reports.py` (+1 more) |
| `first_coordination_sphere_metal_basisset` | own name | `delfin/co2/CO2_Coordinator6.py`, `delfin/co2/chain_setup.py`, `delfin/esd_input_generator.py`, `delfin/occupier.py` (+3 more) |
| `first_coordination_sphere_scale` | own name | `delfin/co2/CO2_Coordinator6.py`, `delfin/co2/chain_setup.py`, `delfin/esd_input_generator.py`, `delfin/occupier.py` (+3 more) |
| `geom_opt` | own name | `delfin/cli_fukui.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/esd_input_generator.py` (+8 more) |
| `freq_type` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/esd_input_generator.py`, `delfin/recalc_control.py` (+3 more) |
| `initial_guess` | own name | `delfin/cli_fukui.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/esd_input_generator.py` (+5 more) |
| `temperature` | own name | `delfin/agent/api_client.py`, `delfin/agent/bench_ab.py`, `delfin/agent/cli.py`, `delfin/agent/run_manifest.py` (+14 more) |
| `maxiter` | own name | `delfin/agent/generate_fact_tasks.py`, `delfin/agent/job_fix.py`, `delfin/api.py`, `delfin/casscf_active_space.py` (+17 more) |
| `qmmm_option` | own name | `delfin/occupier.py`, `delfin/tools/occupier_opi.py`, `delfin/xyz_io.py` |
| `print_MOs` | own name | `delfin/common/orca_blocks.py`, `delfin/recalc_control.py` |
| `print_Loewdin_population_analysis` | own name | `delfin/common/orca_blocks.py`, `delfin/recalc_control.py` |
| `PAL` | own name | `delfin/agent/api_client.py`, `delfin/agent/bench_holdout_glue.py`, `delfin/agent/calc_policy.py`, `delfin/agent/cluster_efficiency.py` (+37 more) |
| `maxcore` | own name | `delfin/agent/api_client.py`, `delfin/agent/calc_budget.py`, `delfin/agent/calc_policy.py`, `delfin/agent/cluster_efficiency.py` (+60 more) |
| `parallel_workflows` | own name | `delfin/cli.py`, `delfin/guppy_sampling.py`, `delfin/hyperpol.py`, `delfin/occupier.py` (+7 more) |
| `pal_jobs` | own name | `delfin/cli.py`, `delfin/guppy_sampling.py`, `delfin/hyperpol.py`, `delfin/occupier.py` (+7 more) |
| `orca_parallel_strategy` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/occupier.py`, `delfin/recalc_control.py` |
| `enable_job_timeouts` | own name | `delfin/recalc_control.py`, `delfin/workflows/scheduling/pool.py` |
| `job_timeout_hours` | own name | `delfin/recalc_control.py`, `delfin/workflows/scheduling/pool.py` |
| `opt_timeout_hours` | own name | `delfin/recalc_control.py`, `delfin/workflows/scheduling/pool.py` |
| `frequency_timeout_hours` | own name | `delfin/recalc_control.py`, `delfin/workflows/scheduling/pool.py` |
| `sp_timeout_hours` | own name | `delfin/recalc_control.py`, `delfin/workflows/scheduling/pool.py` |
| `enable_auto_recovery` | own name | `delfin/occupier_flat_extraction.py`, `delfin/orca.py`, `delfin/recalc_control.py`, `delfin/tools/adapters/orca.py` |
| `max_recovery_attempts` | own name | `delfin/orca.py`, `delfin/recalc_control.py` |
| `enable_adaptive_parallelism` | own name | `delfin/recalc_control.py`, `delfin/workflows/scheduling/pool.py` |
| `enable_performance_metrics` | own name | `delfin/recalc_control.py`, `delfin/workflows/scheduling/pool.py` |
| `MANTA_QUALITY` | own name | `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_SEEDS` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_NUM_CONFS` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_CONSTRUCTION` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_MAX_ISOMERS` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflows/pipeline.py` |
| `MANTA_BINDING_MODES` | own name | `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_HAPTO` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_UFF` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_DETERMINISTIC` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_ENV` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_CLEAN_GATE` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_TOPOLOGY_GATE` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_DEDUP` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_SCREEN` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_SCREEN_ABOVE` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_SCREEN_KEEP` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_OPT` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_OPT_METHOD` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_MULTIPLICITIES` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_REFINE` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_REFINE_TOPK` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_RMSD_CUTOFF` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py`, `delfin/workflows/pipeline.py` |
| `MANTA_ENERGY_WINDOW` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py`, `delfin/workflows/pipeline.py` |
| `MANTA_PARALLEL_JOBS` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `MANTA_TIME_BUDGET` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/common/manta_settings.py` |
| `hyperpol_xTB` | own name | `delfin/cli.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/dashboard/tab_help.py` (+5 more) |
| `hyperpol_xTB_xyz` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/pipeline.py` |
| `hyperpol_xTB_preopt` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/pipeline.py` |
| `hyperpol_xTB_engine` | own name | `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/pipeline.py` |
| `hyperpol_xTB_bfw` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/pipeline.py` |
| `hyperpol_xTB_wavelengths` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/pipeline.py` |
| `hyperpol_xTB_energy_window` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/pipeline.py` |
| `tadf_xTB` | own name | `delfin/cli.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/dashboard/tab_help.py` (+5 more) |
| `tadf_xTB_xyz` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/pipeline.py` |
| `tadf_xTB_preopt` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/contrib/tadf_xtb_workflow.py` (+1 more) |
| `tadf_xTB_excited_method` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/contrib/tadf_xtb_workflow.py` (+1 more) |
| `tadf_xTB_bfw` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/contrib/tadf_xtb_workflow.py` (+1 more) |
| `tadf_xTB_energy_window` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/workflow_connectivity.py`, `delfin/workflows/contrib/tadf_xtb_workflow.py` (+1 more) |
| `tadf_xTB_run_t1_opt` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/dashboard/tab_help.py`, `delfin/workflow_connectivity.py` (+2 more) |
| `thermodynamics` | read as `stability_constant` (alias) | `delfin/cli.py`, `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py` (+4 more) |
| `thermodynamics_mode` | read as `stability_constant_mode` (alias) | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/stability_constant.py` (+2 more) |
| `thermodynamics_reaction` | read as `stability_reaction` (alias) | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/stability_constant.py` (+1 more) |
| `n_explicit_solvent` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/recalc_control.py`, `delfin/stability_constant.py` (+1 more) |
| `logK_exp` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/recalc_control.py`, `delfin/stability_constant.py` (+1 more) |
| `thdy_smiles_converter` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/stability_constant.py` (+1 more) |
| `thdy_preopt` | own name | `delfin/common/control_allowed_values.py`, `delfin/common/control_validator.py`, `delfin/config.py`, `delfin/stability_constant.py` (+1 more) |

## 6. Keys read without a literal occurrence

These are read through a mechanism other than a literal key lookup:

- `keyword:basename` — parsed by shape, not by name: `delfin/common/orca_overrides.py:_OVERRIDE_KEY_RE` matches `keyword:<target>` for any job name, so the key is the prefix `keyword:` plus a target, and no literal `keyword:basename` appears in the code.
- `additions:basename` — parsed by shape, not by name: `delfin/common/orca_overrides.py:_OVERRIDE_KEY_RE` matches `additions:<target>` for any job name, so the key is the prefix `additions:` plus a target, and no literal `additions:basename` appears in the code.

