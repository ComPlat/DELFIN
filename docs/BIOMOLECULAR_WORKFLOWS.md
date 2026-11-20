# Biomolecular workflow foundation

The detailed foundation specification, installer coverage and acceptance matrix
are maintained in [BIOMOLECULAR_FOUNDATION/README.md](BIOMOLECULAR_FOUNDATION/README.md).

This layer provides typed planning, native predictor inputs, executable processing
and selected external engine adapters. Execution can provision supported missing tools through the central installer;
read-only planning never installs. Model weights are not downloaded by this layer. Unimplemented backends remain explicit execution blockers.

## Agent entry points

The DELFIN tools MCP server exposes two read-only tools:

- `bio_backends(task=None)` returns backend contracts, semantic input/output
  types, prerequisites, documentation and JSON schemas for plans and requests.
- `bio_plan(template="docking")` returns a planning graph and its diagnostics.
  Alternatively pass `plan={...}` to validate an agent-authored graph.

Both are available during planning. Neither probes engines, accesses the
network, writes files or submits jobs. The tools platform's `get_guide` explains
this route. Regular registered adapters include:

| Adapter | Work performed |
| --- | --- |
| `bio_backend_catalog` | Writes the catalog as JSON; optional `probe=True` searches PATH without running commands. |
| `bio_workflow_plan` | Writes a template/custom plan and diagnostics. Invalid graph: failed result. Structurally valid but unimplemented graph: successful planning result with `execution_ready=false` and explicit blockers. |
| `bio_prediction_input` | Writes native predictor input and a manifest with request, seeds, MSA policy, dialect and SHA-256. Never starts a model. |

Each adapter provides a valid `example_kwargs` in its contract. Structured
requests cannot be inferred from a placeholder such as an empty dictionary.

## Catalog coverage and status

The catalog has 46 backend/operation contracts, of which 28 have execution
bindings. Additional registered adapters handle prediction import and workflow
execution. Catalog coverage is not a claim that every engine runs in DELFIN.

| Domain | Entries |
| --- | --- |
| Prediction | AlphaFold 2, AlphaFold-Multimer, AlphaFold 3, ColabFold, ESMFold, Boltz-2, Chai-1, RoseTTAFold All-Atom |
| Preparation | Existing RDKit ligand preparation, ligand-state enumeration, Molscrub, explicit receptor preparation, Meeko receptor/ligand conversion |
| Docking | AutoDock Vina, GNINA, smina, DiffDock, HADDOCK3, Glide, GOLD |
| Conversion and validation | Meeko pose export, PoseBusters, model selection, complex assembly |
| Analysis/MedChem | ProLIF, PLIP, properties/scaffolds, assay joining |
| Simulation | Parameterization, OpenMM, GROMACS, Amber, NAMD, MDAnalysis |
| QM connection and design | Explicit QM/MM-region construction, RFdiffusion, ProteinMPNN |

`integration=planned` means the execution adapter still needs implementing.
`integration=adapter` identifies an implemented execution binding; its runtime
may still be unavailable or its input subset unsuitable.
A found executable does not verify its version, dependencies, weights, license
terms, GPU support or scientific applicability. Catalog and planning calls
therefore never claim `execution_ready=true`.

## Executable foundation and verification status

| Area | Available implementation | Verification / boundary |
| --- | --- | --- |
| Ligand processing | RDKit preparation, descriptors/scaffolds, tautomer enumeration, explicit record selection, ligand-derived binding box, SVG grid | Real local adapter tests; no pH populations or measured ADMET |
| Assay processing | Exact chemical-state joins with units, censoring, conditions, replicates and unmatched rows | Real local MedChem DAG tests; no unit conversion or averaging |
| Structure processing | Explicit model/chain/water/hetero/altloc selection and atom mapping | Real PDB/mmCIF tests; coordinate selection is not chemical receptor preparation |
| Receptor preparation | PDBFixer/OpenMM worker with exact versions, pH, water and atom-repair policies, original-atom checks | Boundary tests only; packages absent locally; canonical hydrogen-free single-model PDB subset |
| Complex assembly | Explicit pose selection, shared-frame assertion, PDB serial/CONECT remapping and authoritative ligand SDF sidecar | Real local assembly tests; no coordinate alignment or chemistry inference |
| Prediction results | Explicit import, raw confidence assets, SHA-256 validation, model selection | Local tests, including changed confidence and path escape refusal; engine provenance is caller supplied |
| Visualization | Offline vendored 3Dmol viewer, pose/chain/water controls, binding box; PyMOL and ChimeraX scripts | Export tests; desktop scripts syntax checked; actual browser/WebGL and desktop applications not verified |
| Docking | Vina, GNINA, smina adapters with exact version banner, explicit seed/box/CPU count, timeout, checked outputs | Vina 1.2.5 synthetic live I/O smoke passed; GNINA/smina not installed/tested live; no scientific benchmark |
| Preparation/export | Meeko ligand/receptor PDBQT and pose SDF export in a selected Python environment | Implemented optional adapters; Meeko absent locally. Receptor subset: one prepared canonical PDB model, no unresolved altlocs/cofactors/metals |
| Validation | PoseBusters redocking checks retaining raw per-pose CSV values | Mocked boundary tests; package absent locally; failed checks remain visible even when report execution succeeds |
| Simulation setup | OpenMM ff14SB protein-only vacuum or TIP3P/PME parameterization with an explicit constraints/charge policy | Checked system/topology/coordinates/force-field bundle; mocked engine boundary, no real parameterization validated here |
| Prediction | ESMFold v1 single-protein inference, verified local weights, explicit device and runtime versions | Mocked engine boundary; no inference or GPU validation performed here |
| Interactions | ProLIF fingerprints with molecule sidecars, atom mapping and parent contact indices | Canonical hydrogenated PDB subset; mocked engine boundary; chemistry/coordinates/completeness checked |
| Interactions | PLIP XML report for an explicitly assembled PDB complex | Implemented optional adapter; package absent locally |
| Execution | `bio_workflow_preflight`, `bio_workflow_execute` | Actual typed DAG through `run_step`, stop on failure, checksummed outputs, resume refusal after changed inputs/settings/runtime/code |

`bio_workflow_preflight` reads actual files and may run bounded version probes.
It is separate from the read-only planning MCP tools. Exact package versions are
required for optional Python tools; select their environment with `python`.
Docking requires `expected_version` equal to the complete stripped `--version`
banner. External calculations use DELFIN's tracked process runtime. Provisioning is a
separate checked phase, not a fallback hidden inside an engine subprocess. Timeouts terminate the process group. Cancellation callbacks stop
between DAG nodes; active subprocess cancellation uses DELFIN's process manager.
Execution inside a job is sequential. The shared `delfin-step --slurm` CLI can submit the entire Bio DAG with explicit CPU, total-memory, walltime and GPU/partition requests. OpenMM currently accepts CPU/Reference only; the executor does not parallelize DAG nodes or submit individual nodes.

Successful resume verifies recorded files plus plan, input hashes, relevant code,
Python dependency versions and probed executable identity. It is not an atomic
snapshot of all installed libraries, drivers or weights. Only one writer may run
against a workflow directory. Output records retain the original coordinate
frame; no implicit structure alignment or receptor preparation occurs.

Run an offline example with synthetic structures:

```sh
python -m examples.biomolecular_foundation.run_processing --output /tmp/delfin-processing
```

Open `/tmp/delfin-processing/run/view_1/view.html` locally. The directory includes
its JavaScript, input coordinates, desktop viewer scripts and manifest; no web
server is needed. `plan.json` documents all graph edges; `run/run.json` records
actual node outcomes. Invoke custom plans with
`run_step("bio_workflow_execute", plan=plan, work_dir=directory)`; the planning
DAG itself is not a generic DELFIN pipeline spec.

### Remaining foundation work

Predictor execution (AlphaFold family, ColabFold, Boltz, Chai and
RoseTTAFold), automatic model asset acquisition, GPU scheduling, richer receptor chemistry support,
Molscrub, ligand MD parameterization and additional solvent/force-field models and
GROMACS/Amber/NAMD execution remain unimplemented bindings. OpenMM minimization/dynamics and MDAnalysis statistics have bounded adapters for checked vacuum or fixed-volume TIP3P/PME systems and JSONL trajectories. Periodic analysis requires explicit `as_stored` coordinate policy; unwrapping is not implemented.
The same applies to DiffDock, HADDOCK3, licensed docking, protein design and
QM/MM-region construction. Existing native predictor input generation does not
remove these blockers. Scientific reference datasets and installed-engine tests
are needed before calling any complete docking/MedChem/MD workflow validated.

External adapter interfaces follow the official
[Vina tutorial](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html),
[GNINA CLI](https://github.com/gnina/gnina),
[Meeko CLI](https://github.com/forlilab/Meeko/tree/develop/meeko/cli),
[PoseBusters CLI](https://github.com/maabuu/posebusters/blob/main/posebusters/cli.py)
and [PLIP](https://github.com/pharmai/plip). Pin installed versions; documentation
of a CLI is not evidence of a successful local scientific run.

## Explicit wiring

Plans have version 1, named external `inputs` and ordered `nodes`. Each node
names a backend and binds every input to `inputs.NAME` or `NODE.OUTPUT`.
References must point backward; duplicate IDs, unknown ports, mismatched
semantic types and forward/cyclic dependencies are rejected.

The docking template follows:

```text
structure → receptor preparation → Meeko receptor ─┐
ligand → ligand preparation → Meeko ligand ────────┼→ Vina → Meeko export
binding box ─────────────────────────────────────┘             │
                           prepared receptor + reference ─────┼→ pose validation
                           prepared receptor ─────────────────┴→ interactions
```

The vocabulary distinguishes `macromolecule`, `prepared_receptor`,
`receptor_pdbqt`, `ligand`, `ligand_pdbqt`, `poses_pdbqt` and `pose_set`.
Predictions produce a `prediction_bundle`, which requires explicit model
selection. Cofolded ligands require explicit extraction; selection alone does
not provide a prepared docking receptor. Simulation requires a parameterized
system, not merely a coordinate file. An `alphafold3_input` cannot feed Boltz.

The planning validator checks graph structure only. The separate execution
preflight validates bound adapter settings, reads external artifacts and probes
selected runtimes. It does not establish scientific coordinate compatibility.

Templates are provided for `processing`, `docking`, `medchem`, `md` and each predictor ID.
The public JSON schemas and semantic validator complement each other: the
schema describes shape, while the validator checks IDs, references and ports.

## Native prediction inputs

Example request (synthetic sequence, for interface demonstration only):

```json
{
  "schema_version": 1,
  "name": "example",
  "seeds": [7],
  "msa_policy": "single_sequence",
  "entities": [
    {"id": "A", "kind": "protein", "sequence": "MKWVTFISLL"},
    {"id": "L", "kind": "ligand", "smiles": "CCO"}
  ]
}
```

```python
from delfin.tools import run_step

# request is the dictionary above; choose a fresh output directory.
result = run_step("bio_prediction_input", backend="alphafold3",
                  request=request, work_dir="af3_input")
result.require()
print(result.artifacts)  # native input + manifest, no predicted structure
```

| Exporter | Current supported subset |
| --- | --- |
| AlphaFold 3 | Basic version-1 local AF3 JSON, unmodified protein/DNA/RNA and SMILES ligands; explicit seeds. |
| Boltz-2 | Version-1 YAML, emitted using JSON syntax (valid YAML), typed entities and explicit single-sequence MSA fields. |
| Chai-1 | Typed FASTA for unmodified polymers and SMILES ligands; single-sequence mode. |
| AlphaFold 2 / Multimer | Protein FASTA; explicit engine-pipeline MSA policy; multiple chains only for the Multimer exporter. |
| ColabFold / ESMFold | One protein FASTA record per request. Multimer and batch semantics are deliberately not inferred. |

Native input generation is tested against documented shapes, not installed
model parsers or inference. Engine versions must still be pinned and checked by
future execution adapters. FASTA cannot encode all settings: runners must read
the manifest and apply model preset, seeds and MSA policy. Boltz/Chai seeds also
remain runner settings. Single-sequence mode is an explicit computational
choice, not an accuracy recommendation. Custom MSAs, templates, modifications,
CCD components, restraints and covalent bonds need dedicated conversion support;
unknown fields fail rather than being silently dropped.

The old `structure_prediction_engine` remains a limited detection interface.
Its ESMFold example uses `esm-fold`; fabricated generic FASTA commands for
AlphaFold 3 and Boltz have been removed. Use the new catalog and converters.

## Chemical identity and artifact contracts

The SDF preparation path now preserves the original molecular graph through
RDKit, including formal charges, isotopes and stereochemistry. Formal charges
are not Gasteiger partial charges. Prepared SDF properties include
`DELFIN_IDENTITY`: a canonical isomeric-SMILES chemical-state ID, the RDKit
version, formal charge and a mapping from output atoms to this operation's
input atom indices. Added hydrogens have no input index. Removed fragments do
not renumber the original indices. This is a per-operation mapping, not a claim
of globally unique atom labels or full cross-engine identity tracking. Identity
is versioned with the canonicalizer and does not merge protonation states or
tautomers. Pose comparisons use the preserved graph when present.

Adapters can opt into `required_artifacts`. The runner then fails a nominally
successful step if a required file is missing, empty or unreadable. The new
planning/input adapters and ligand preparation use this contract. Legacy
`produces` declarations remain report-only for compatibility. A nonempty file
is not proof of chemical validity: each execution adapter must parse and check
its own outputs before reporting success.

## Extending the scaffold

1. Add/update a `Backend` contract in `delfin/tools/bio/workflows.py`, naming
   semantic ports, prerequisites, limits and the official documentation.
2. Implement a regular `StepAdapter` with declared parameters, input wiring,
   output ports, valid example kwargs and required artifacts. Keep external
   installations isolated from DELFIN's Python environment where needed.
3. Add native input/output converters. Preserve model/chain/atom mapping,
   coordinate frame, units, engine/version, settings, seeds, input hashes and
   model/force-field provenance. Never relabel confidence as affinity.
4. Test malformed inputs, missing assets, failures, cancellations and artifact
   validation; then run an explicitly recorded integration/reference example.
5. Associate the catalog entry with the tested adapter. Implement an explicit
   graph-to-runtime binding before calling a planned chain runnable.

## Official format references

- [AlphaFold 3 input](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md)
- [Boltz prediction inputs](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)
- [Chai inference](https://github.com/chaidiscovery/chai-lab)
- [ESMFold CLI](https://github.com/facebookresearch/esm)
- [ColabFold](https://github.com/sokrypton/ColabFold)
- [Vina/Meeko preparation and pose export](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html)

Other engine documentation links are included in each catalog entry.

## Validation of this update

The focused suite covers real RDKit processing, artifact integrity, DAG execution,
resume rejection, failure/cancellation handling, optional-engine boundaries,
synthetic live Vina I/O and existing biomolecular readers/viewer helpers.
Three pre-existing failures reproduce on unchanged main
`03a02d048f674a972d6c2e2551908eb5ece24b6c`:
`test_declared_required_params_are_enforced` (`orca_th_opi.output_file`),
`test_every_adapter_has_a_contract` (`orca_qmmm_opi`), and
`test_validation_names_upstream_producers_for_missing_input` (hint truncation).

The planning CLI writes 46 catalog entries, six plans and three native input
sets without executing an engine. The separate processing example executes a
real four-node offline DAG. Synthetic spatial SDF fixtures explicitly carry the
3D dimension marker; this does not make them chemically realistic benchmarks.
Native predictor formats have local content tests; real engine-parser/inference
compatibility remains an integration check. The optional `jsonschema` package is
absent, so Draft-2020-12 meta-validation was not performed.

Assay tables use `schema_version: 1` and `records`. Each measurement requires
`measurement_id`, the exact DELFIN `chemical_state_id`, `assay_id`, `endpoint`,
finite numeric `value`, explicit `unit`, `relation`, `source` and a `conditions`
object. Replicates need distinct measurement IDs. The public `assay_table` schema
is included in the backend catalog. Endpoint names and units are retained as
supplied; normalization and quality assessment remain separate operations.

## Receptor preparation and complex assembly

`bio_receptor_prepare` binds `receptor_prepare` to PDBFixer/OpenMM in an explicit
Python environment. Required settings are exact `pdbfixer_version` and
`openmm_version`, `ph` (0..14), `seed`, `waters` (`keep`/`drop`) and
`missing_atoms` (`reject`/`repair`). It uses the CPU Reference platform. Missing supported packages can be provisioned
in a separate managed environment before execution. This first subset accepts one hydrogen-free canonical protein
PDB with resolved alternate locations. Modified residues, cofactors, metals,
nonprotein components and detected missing residues are refused. Existing
protonation requires a separate policy; no residues or loops are built implicitly.
Hydrogen assignment uses default residue rules rather than an electrostatic pKa
calculation. Missing-residue detection depends on sequence records. Output atom
identities and original heavy-atom coordinates are checked; additions/removals
are reported. Original annotations remain in `input.pdb`. No system force-field
parameterization or scientific validation is implied.

The API sequence follows the [PDBFixer manual](https://github.com/openmm/pdbfixer/blob/master/Manual.html)
and [implementation](https://github.com/openmm/pdbfixer/blob/master/pdbfixer/pdbfixer.py).
PDBFixer/OpenMM are absent in this development environment, so adapter boundary
checks use a simulated worker; no successful live receptor repair is claimed.

`bio_assemble_complex` binds `assemble_complex` to native DELFIN processing.
Required settings are a 1-based `record_id` and `shared_coordinate_frame=true`.
The latter is the caller's assertion, not a computed alignment check. An optional
single-character `ligand_chain` defaults to `L`; collisions are refused. Original
receptor coordinate fields and residue identities are preserved, serial numbers
and CONECT references are remapped, and the selected ligand's graph is preserved
in `selected_ligand.sdf`. PDB output alone cannot preserve every chemical detail;
ligand coordinates are rounded to PDB precision. Unsupported anisotropic records,
serial overflow, unresolved alternate locations and multiple receptor models
are refused rather than silently discarded.

## Shared installation and health checks

The central `delfin.installer` now includes a `bio` family. Initial conda-forge
recipes cover Vina, Meeko, PDBFixer, OpenMM, MDAnalysis, ProLIF, PLIP and PoseBusters. Top-level releases
are pinned; transitive dependencies are resolved by conda and their exact
versions/builds/hashes are recorded from `conda-meta`. These recipes have
interface tests, not successful installation evidence in this environment.

Executable discovery reuses the QM resolver, including module profiles, PATH,
overrides and managed binaries. Python selection checks the tool-specific
`DELFIN_<TOOL>_PYTHON`, current/available interpreters, interpreters beside known
entry points and the managed prefix. A caller-selected nondefault Python stays
explicit. Python-only cluster module discovery remains to extend.

Read-only preflight never installs. Actual workflow execution validates all
inputs/settings first, then may provision supported missing runtimes. The same
`DELFIN_AUTO_INSTALL_QM_TOOLS=0` opt-out disables automatic installation. Existing
compatible environments win; managed installs go under `DELFIN_BIO_TOOLS_ROOT`
(default `~/.local/share/delfin/bio_tools`) in separate recipe/version prefixes.
No cluster or DELFIN Python environment is modified. Installation uses an
existing micromamba/mamba/conda; the existing DELFIN bootstrap supplies that
manager when needed. Installations are locked per prefix, version/import-checked
and logged. Failed automatic attempts are not repeated indefinitely.

Use `python -m delfin.installer --install meeko pdbfixer` for the central install
route. Explicit repair can only modify a managed prefix whose ownership and
recipe manifest match; unrelated environments are refused. Repair checks are
technical checks, not scientific reference calculations. Install timeouts now
terminate the installer process group, including children holding log pipes.

Automatic installation is not yet supported for every catalogued backend.
Unsupported engines, model assets, licensed software and incompatible versions
remain explicit blockers. No installation or model download was run during this
implementation session; tests simulate installers and package availability.

### ESMFold and ProLIF execution boundaries

`esmfold_predict` consumes one protein FASTA record and a local `assets_file`.
The asset manifest has `schema_version: 1`, `engine: "esmfold"` and a `files`
array of `{path, sha256}` entries relative to the manifest. Required entries are
`checkpoints/esmfold_3B_v1.pt`, `checkpoints/esm2_t36_3B_UR50D.pt` and
`checkpoints/esm2_t36_3B_UR50D-contact-regression.pt`. All declared files must
exist and match their hashes. The worker disables Torch Hub downloads. Select
an already prepared Python environment with exact `expected_version` (fair-esm)
and `torch_version`, explicit `device` (`cpu` or `cuda:N`) and `seed`.
An automatic fair-esm installation recipe is not yet provided. Returned models
must match the requested sequence and contain a complete backbone; confidence
values must be finite and within their declared ranges. These checks do not
establish prediction accuracy.

`prolif_analyze` requires explicit `interactions`, `expected_version`,
`bond_policy="rdkit_proximity"` and `shared_coordinate_frame=True`. Both partners
must already have explicit hydrogens. The current receptor subset is one
canonical protein PDB without waters, alternate locations or cofactors. RDKit
infers receptor bonds; the adapter does not prepare chemistry or align poses.
Each pose retains a molecule sidecar and source atom mapping. Exports are checked
for unchanged ligand chemistry, mapped bonds and coordinates; receptor atom
identity and coordinates are checked, and contact parent indices must be in
bounds. Interaction geometry is not binding affinity. Tests use synthetic data
and mocked engine boundaries; installed ProLIF execution remains to be validated.

Interfaces follow the [ESMFold implementation](https://github.com/facebookresearch/esm)
and [ProLIF docking documentation](https://prolif.readthedocs.io/en/stable/notebooks/docking.html).

### Python tools supplied by cluster modules

For managed Python recipe families, runtime discovery checks compatible existing
interpreters, then matching cluster modules, before considering automatic
installation. Module lookup reuses DELFIN's QM discovery and honors
`DELFIN_NO_MODULE_PROBE=1`. A bounded, tracked login-shell probe resolves Python
from the loaded module's PATH, including dependencies loaded by that module.
Each candidate must import the recipe's modules and match the requested package
version. Modules exposing the same Python path are checked separately.

The chosen module is recorded as `environment_module` in runtime provenance and
used again for distribution checks, console entry points and worker execution.
Multi-package workers check all dependencies in that same module environment.
Module-backed Python uses `-s -P`: it excludes user-site packages and automatic
script/current-directory imports while retaining the module's PYTHONPATH. This
requires Python 3.11 or newer. Ordinary isolated environments retain `-I`.
Explicit virtual-environment Python paths preserve their symlinks so that Python
can find `pyvenv.cfg` and select the intended environment.

Tests simulate module discovery and launches, and exercise a real temporary
Python venv without installing packages. Actual cluster module stacks still need
acceptance testing. Module identity and exact distribution versions are recorded;
this is not a complete snapshot of modulefiles, shared libraries or GPU drivers.

### OpenMM system foundation

The `parameterize` DAG backend now binds to `bio_openmm_parameterize`. Its first
supported subset is an already hydrogenated canonical single-model protein PDB
with `forcefield="amber14/protein.ff14SB.xml"`, `environment="vacuum"`, explicit
`constraints` (`none` or `HBonds`) and integer `expected_charge` in elementary
charge units. `expected_version` pins OpenMM; the shared environment resolver and
managed installer apply. There is no implicit protonation change, solvation,
ligand parameter assignment or switch from periodic to vacuum conditions.
Nonstandard residues, waters, ligands, alternate locations and CRYST1 boxes are
refused. OpenMM must match all residue templates; unmatched atoms fail the run.

The worker uses OpenMM's `ForceField.createSystem` with `NoCutoff`, stores the
actual packaged ff14SB XML and round-trips the serialized `System`. It exports
indexed topology atoms/bonds and positions with explicit nanometer units.
`simulation_system.json` hashes these files and the unchanged input structure,
and records runtime, expected charge, constraints and boundary-condition policy.
Native artifact validation checks all checksums, atom ordering/identity, coordinate
conversion from PDB angstroms, particle counts/masses, net charge and constraint
pairs. Missing or modified nested files invalidate the bundle. Virtual sites and
force-field includes are outside this first schema. The checks do not establish
scientific accuracy or complete validation of every force parameter.

Parameterization itself is system construction only. Minimization, dynamics and checkpoint continuation are separate explicit `bio_openmm_run` operations.
Protein-only vacuum construction is an explicit technical subset, not a recommended
production protocol for a solvated protein or a protein-ligand complex. Cluster
acceptance must exercise actual OpenMM template matching and compare the resulting
system with a trusted reference before adding dynamics.

References: [OpenMM simulation setup](https://docs.openmm.org/latest/userguide/application/02_running_sims.html),
[system serialization](https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Saving%20Systems%20to%20XML%20Files.html)
and [packaged ff14SB definition](https://github.com/openmm/openmm/blob/8.4.0/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml).

### Minimization, dynamics and trajectory analysis

`bio_openmm_run` / backend `openmm` consumes a checked system bundle. Set exact
OpenMM version, `platform` (`CPU` or `Reference`), `temperature_k`,
`friction_per_ps`, `timestep_ps`, nonzero `seed`, `constraint_tolerance`, `steps`,
`report_interval` and `minimize` explicitly. Reference requires one CPU core.
The integrator is LangevinMiddle; the system's chemical state and boundary
conditions remain those of the bundle. Initial constraints are applied before
optional bounded minimization; initial velocities are sampled at the requested
temperature. `steps=0, minimize=True` performs minimization only. A returned
minimizer with non-increasing finite energy is not an independent convergence
certificate. No equilibrated-production status is inferred.

The trajectory is a streamed JSONL file plus a manifest with nm/ps/kJ mol-1 units,
full nested system identity and exact step/time coverage. Initial and final frames
are always recorded. `max_frames` bounds the requested frame count. Checkpoints
include the engine binary, serialized State and Integrator, protocol, final
step/time and associated trajectory. Validation checks positions/velocities,
serialized integrator parameters and agreement with the final trajectory frame.
Set `checkpoint_file` to the checkpoint manifest and `minimize=False` to continue.
Changing the system, runtime version, platform, CPU count, seed or integrator
settings is refused. OpenMM may reject a binary checkpoint across hardware/build
changes; there is no silent state-file fallback. Checkpoints are currently written
on successful completion, not periodically during interrupted jobs.

`bio_mdanalysis` / backend `mdanalysis` checks trajectory/system identity, then
uses explicit zero-based `atom_indices`, `alignment` (`none` or `fit`),
`weighting` (`uniform` or `mass`), `discard_frames` and `stride`. Fitting uses
MDAnalysis rotation matrices against the first retained frame and requires at
least three non-collinear selected atoms. Statistics include selected-atom RMSD,
per-atom population RMSF and radius of gyration with the selected weighting.
No periodic unwrapping, equilibration detection, affinity inference or convergence
claim is made. MDAnalysis 2.9.0 has an isolated conda-forge recipe; no solver or
installed MDAnalysis engine was run here.

Tests cover mocked OpenMM boundaries, missing/corrupt outputs, protocol-mismatch
restart refusal, streaming numerical analysis with a mocked fitting boundary,
and the three-node parameterize/run/analyze DAG with checksum-verified resume.
Cluster reference runs remain required.

### Submit an explicit Bio plan through the existing Slurm CLI

The shared CLI accepts `--params-file parameters.json`, a JSON object such as
`{"plan": {...}}`. CLI parameters override file values. The MD plan builder
`python -m examples.biomolecular_foundation.build_md_plan --help` lists all
required scientific choices and writes this file without running any engine.
Use absolute input paths on a filesystem visible to compute nodes.

```sh
delfin-step bio_workflow_execute --params-file /shared/md-parameters.json \
  --work-dir /shared/md-review-run --cores 4 --slurm \
  --slurm-time 01:00:00 --slurm-memory-mb 8192
```

For a compatible GPU backend (for example ESMFold configured with an allocated
CUDA device), additionally supply `--slurm-gpus 1 --slurm-partition GPU_PARTITION`.
GPU requests require an explicit partition; this does not add GPU support to
CPU-only adapters. Module/tool version probes and provisioning occur within the
compute job. Script generation snapshots nested parameters as JSON, quotes shell
arguments, preserves the submission working directory and returns the step's
actual exit code. Existing submission files are not overwritten. Slurm job IDs
are printed by the existing CLI; native scheduler status/cancellation remain
available through DELFIN's existing scheduler integration. No real job was
submitted during development.

References: [OpenMM Context/checkpoints](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.Context.html),
[OpenMM minimization](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LocalEnergyMinimizer.html),
[MDAnalysis fitting](https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html).

### Explicit TIP3P solvent and fixed periodic boxes

`bio_openmm_parameterize` also accepts `environment="tip3p_pme"` with the
following `solvation` object (values below illustrate the schema, not a validated
production protocol):

```json
{"box_size_nm": [4.0, 4.0, 4.0], "ionic_strength_m": 0.15,
 "neutralize": true, "positive_ion": "Na+", "negative_ion": "Cl-",
 "seed": 7, "solute_charge": 0, "cutoff_nm": 0.9,
 "ewald_tolerance": 0.0005}
```

The input is still one already hydrogenated canonical protein without an existing
periodic box. The worker uses the packaged Amber14 TIP3P/ion definition and
OpenMM Modeller; Na+ or K+ with Cl- are supported. The final expected charge must
match the explicit neutralization policy. The box must contain the protein's
coordinate extent, and the cutoff must be less than half its shortest side.
Protein atoms keep their identities, ordering and positions. Added water/ions,
ion charges/counts, both force-field files, orthorhombic box vectors and PME
parameters are checked. Ion-count validation follows the pinned OpenMM Modeller
water-concentration convention, including rounding and neutralization ions.

Version 2 system bundles retain both original protein and solvated PDB plus the
actual solvent force-field XML. The worker rereads its written PDB so topology,
coordinates and system consistently use its stored precision. PDB atom, residue
and chain-ID limits are checked; larger systems requiring mmCIF are refused.
Water flexibility follows the explicitly chosen `none`/`HBonds` constraints;
no additional rigid-water or barostat policy is silently introduced.

The dynamics adapter supports these fixed-volume periodic systems and records
box vectors in every frame. Changed/missing boxes and checkpoint disagreement
fail validation. NPT/barostats and varying-box trajectories are not supported.
MDAnalysis rejects periodic input by default; `periodic_policy="as_stored"`
explicitly opts into analysis of stored coordinates without periodic unwrapping.
Select appropriate whole-molecule coordinates for the scientific question.
The plan builder accepts `--environment tip3p_pme --solvation-file policy.json
--periodic-policy as_stored`.

Tests mock the OpenMM boundary and reject source changes, charge/box/cutoff errors
and invalid policies. No real solvation or PME calculation has been run here.
References: [OpenMM Modeller](https://docs.openmm.org/latest/api-python/generated/openmm.app.modeller.Modeller.html)
and [OpenMM 8.4 Modeller implementation](https://github.com/openmm/openmm/blob/8.4.0/wrappers/python/openmm/app/modeller.py).
