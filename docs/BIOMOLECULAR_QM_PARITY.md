# Biomolecular parity with DELFIN's QM foundation

This is the implementation checklist for the requested QM-equivalent foundation.
A catalog contract, an implemented adapter, a tested interface and a validated
scientific method are different states. External scientific runs are not possible
in this development session. No package/model installation is performed here.

| Foundation | QM mechanism to reuse | Bio status / remaining acceptance criteria |
| --- | --- | --- |
| Central installation | `installer.Tool`, `install`, `present`, `repair` | In progress: shared bio family, isolated recipe/version prefixes, explicit managed repair and import/version checks. Extend verified recipes to remaining free tools. |
| Cluster discovery | `qm_runtime.ToolSpec`, system/module discovery | Executable profiles and Python module discovery reuse QM mechanisms. Compatible module environments are tried before installation; package checks and workers retain the same module context, including PYTHONPATH. Virtual-environment Python symlinks retain their environment. Real cluster acceptance remains pending; explicit interpreter overrides remain available. |
| Execution provisioning | `qm_health.ensure_tool`, shared auto-install opt-out | In progress: no installation from lists or read-only preflight; provision only after all settings/inputs validate. Complete adapters and provisioning failure tests. |
| Predictor execution | Versioned real process adapters | Native input generation and an ESMFold runner with hashed local assets exist (mocked engine boundary only). Remaining predictors, model acquisition and explicit MSA/template/network policies are open. |
| Resource scheduling | `WorkflowJob`, workflow/global manager, Slurm infrastructure | Serial DAG exists. Integrate CPU/GPU/memory/walltime requests, cluster submission, parallel batches and active-job cancellation/restart. |
| Scientific identity | Typed artifacts, provenance and result checks | Ligand graphs, mapping, units and checksums exist. Checked OpenMM system bundles now preserve topology/atom order, coordinates/units, serialized system and force-field hashes. Trajectory identity and full chain/residue remapping remain open. |
| Preparation | Checked adapters and immutable chemical intent | Canonical PDBFixer subset and Meeko exist. Add explicit protonation ensembles, supported cofactors/metals/modifications and prepared-state validation. |
| Analysis and visualization | Registered analysis blocks and structured reports | Offline views, descriptors, RMSD, PLIP/PoseBusters and assay joins exist. ProLIF canonical-PDB adapter now checks molecule sidecars, atom mappings and contact indices (mocked boundary only). Alignment, ensemble/trajectory processing and comparative reports remain open. |
| Simulation | Real engine adapters plus restart artifacts | Protein-only ff14SB vacuum parameterization now has an engine adapter and checked system bundle (mocked boundary only). Ligand/solvent/periodic parameterization, minimization, equilibration, MD and trajectory analysis remain open. |
| Recovery | Failure classification and chemistry-preserving retry rules | Failure stop and checksummed DAG resume exist. Add bounded engine-specific recovery without changing ligand state, method, force field or target. |
| QM/MM connection | Existing ORCA QM/MM validation and adapter | Implement bio-region selection, boundary/caps, charge/multiplicity and embedding manifests. |
| Agent integration | Tool contracts, examples, guides, safety gates | Registry/planning/execution exist. Cover new environment/recovery/resource capabilities in actual agent-facing guides and tests. |
| Verification | Scientific probes and reference benchmark suites | Local parser/adapter tests and one prior synthetic Vina smoke exist. Provide external acceptance suite and reference-data contracts; record real results only when run externally. |

Implementation stays on Draft PR #746 in its isolated checkout. Nothing in this
checklist marks an unimplemented integration or an unrun external check complete.
