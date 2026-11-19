# Biomolecular workflow scaffold

This is the planning and input-contract layer for building workflows with the
DELFIN agent. It does not install or execute the catalogued external engines.
Use the existing tools platform for executable adapters; keep unimplemented
backends as explicit blockers while constructing a workflow.

## Agent entry points

The DELFIN tools MCP server exposes two read-only tools:

- `bio_backends(task=None)` returns backend contracts, semantic input/output
  types, prerequisites, documentation and JSON schemas for plans and requests.
- `bio_plan(template="docking")` returns a planning graph and its diagnostics.
  Alternatively pass `plan={...}` to validate an agent-authored graph.

Both are available during planning. Neither probes engines, accesses the
network, writes files or submits jobs. The tools platform's `get_guide` explains
this route. Three regular registered adapters are also available:

| Adapter | Work performed |
| --- | --- |
| `bio_backend_catalog` | Writes the catalog as JSON; optional `probe=True` searches PATH without running commands. |
| `bio_workflow_plan` | Writes a template/custom plan and diagnostics. Invalid graph: failed result. Structurally valid but unimplemented graph: successful planning result with `execution_ready=false` and explicit blockers. |
| `bio_prediction_input` | Writes native predictor input and a manifest with request, seeds, MSA policy, dialect and SHA-256. Never starts a model. |

Each adapter provides a valid `example_kwargs` in its contract. Structured
requests cannot be inferred from a placeholder such as an empty dictionary.

## Catalog coverage and status

The initial catalog has 38 backend/operation contracts. These are extension
points, not a claim that 38 integrations run in DELFIN.

| Domain | Entries |
| --- | --- |
| Prediction | AlphaFold 2, AlphaFold-Multimer, AlphaFold 3, ColabFold, ESMFold, Boltz-2, Chai-1, RoseTTAFold All-Atom |
| Preparation | Existing RDKit ligand preparation, ligand-state enumeration, Molscrub, explicit receptor preparation, Meeko receptor/ligand conversion |
| Docking | AutoDock Vina, GNINA, smina, DiffDock, HADDOCK3, Glide, GOLD |
| Conversion and validation | Meeko pose export, PoseBusters, model selection, complex assembly |
| Analysis/MedChem | ProLIF, PLIP, properties/scaffolds/similarity, assay joining |
| Simulation | Parameterization, OpenMM, GROMACS, Amber, NAMD, MDAnalysis |
| QM connection and design | Explicit QM/MM-region construction, RFdiffusion, ProteinMPNN |

`integration=planned` means the execution adapter still needs implementing.
`integration=adapter` currently identifies the existing `ligand_prepare` path.
A found executable does not verify its version, dependencies, weights, license
terms, GPU support or scientific applicability. Catalog and planning calls
therefore never claim `execution_ready=true`.

Execution integrations should reuse DELFIN's existing runtime/scheduler,
environment handling and job controls. They must implement version checks,
resource requests, timeout/cancellation and restart semantics. These controls
are not implemented by the planning scaffold. Do not route a planning graph
directly into `run_application` or submit it as a pipeline spec.

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

This validator checks graph structure. It does not yet validate engine-specific
settings, inspect external paths, establish coordinate compatibility or compile
a graph into executable tasks. Those are reported as remaining execution work.

Templates are provided for `docking`, `medchem`, `md` and each predictor ID.
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

The selected Python 3.12 suite completed with 383 passed and three failures.
The two existing contract failures (`test_declared_required_params_are_enforced`
and `test_every_adapter_has_a_contract`) and the upstream-producer hint failure
(`test_validation_names_upstream_producers_for_missing_input`) reproduce on
unchanged main `03a02d048f674a972d6c2e2551908eb5ece24b6c`.

The command-line example writes 38 catalog entries, six plans and three native
input sets. No inference/docking/MD engine was executed. The native formats are
checked by local content tests; real engine-parser compatibility remains an
integration check. The optional `jsonschema` package is absent here, so a
Draft-2020-12 meta-validation run was not performed.
