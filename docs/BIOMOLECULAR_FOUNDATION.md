# Biomolecular foundation — draft scope and next steps

This draft extends DELFIN's existing, composable tool layer to biomolecular
work. The goal is one vocabulary and one execution/reporting interface for
structures, ligands, binding sites and poses. Docking and structure prediction
should become optional adapters in that system, rather than separate scripts
with incompatible inputs and unverifiable outputs.

## What is implemented

| Surface | Present in this draft | Boundary |
| --- | --- | --- |
| Structure readers | PDB, mmCIF, SDF/MOL and MOL2; compressed input; models, chains, residue identity, alternate locations, waters and hetero components | Reading a structure does not prepare a receptor or repair missing chemistry. |
| Shared input types | `macromolecule`, `ligand`, `binding_site`, `complex`, `pose_set` and `sequence` | These describe files or descriptors, not a successful scientific workflow. |
| `biostructure_summary` | Composition and source-reported missing residues and links, with a JSON artifact | An inventory of the file, not a structure-quality certificate. |
| `ligand_prepare` | RDKit-based coordinates, hydrogens and Gasteiger charges; explicit fragment policy and seed | Protonation is taken as supplied; there is no pH-dependent microstate or tautomer enumeration. |
| `pose_rmsd` | Symmetry-aware heavy-atom RMSD in the existing coordinate frame; failed poses retained in the report | It evaluates supplied poses; it does not run docking or estimate an affinity. |
| `structure_confidence` | Per-residue confidence and optional PAE input, with a JSON artifact | File interpretation and confidence provenance still need the checks below. |
| `structure_prediction_engine` | FASTA parsing, engine detection on a supplied PATH and described commands | Detection only. No predictor is executed, installed or downloaded. |
| Chat card payload | Python-derived chain, hetero-component and water facts, plus an optional binding-site descriptor | This is the decision/data layer; the specialised browser view is not yet complete. |

The shared validator accepts a consumed port when the caller supplies its wired
input parameter explicitly. The runner also reports declared products missing
from a successful result. That report does not make the step fail automatically;
callers must not mistake it for proof that every declared artifact exists.

## Work to finish before claiming the foundation is complete

1. **Connect the protein view end to end.** `chat_viewer` attaches `bio` facts,
   but `_show_molecule_card_html` currently calls the generic 3D renderer without
   passing them on. `_mol3d_payload` recognises PDB but not mmCIF. Wire cartoon
   styling, chain selection, hetero components, water visibility and binding-site
   display into the rendering path. Test the final card and browser behaviour,
   not only the JSON payload.
2. **Make confidence provenance explicit.** Equal per-residue B-factors within
   0–100 do not by themselves identify pLDDT. Preserve the producer and metric
   identity, offer an explicit interpretation where needed, and distinguish
   unavailable confidence from an experimental origin. A missing metric does
   not establish how a structure was obtained. Check PAE dimensions against the
   selected model and keep model, chain and residue mapping in the report.
3. **Validate descriptors at the boundary.** Reject non-finite coordinates,
   sizes, scores and confidence values; check declared units and preserve residue
   identity across models and alternate locations. Maintain path checks for
   complex components and plain report filenames on direct and runner routes.
4. **Add a small recorded end-to-end example.** A receptor, bound reference
   ligand and pose set should pass through summary, preparation and pose
   evaluation via ordinary discovery and pipeline wiring. Pin inputs, expected
   artifacts and provenance. An offline recorded run and an operator-run live
   calculation are separate pieces of evidence.

## Next adapters to add, in order

### Receptor and ligand preparation

Expose explicit choices for model, chains, alternate locations, waters, ions,
cofactors and covalently attached ligands. Record every deletion, addition and
protonation assumption. Produce a new file; never silently overwrite the input.
Add pH-dependent ligand states and tautomers as optional procedures with their
method and limitations named. Keep the current as-given path available.

### Docking execution

Start with one versioned, optional docking backend. Define receptor, ligand and
binding-site inputs; deterministic controls where supported; resource requests;
termination and failure states; and a `pose_set` artifact carrying score units,
engine version and run settings. Use DELFIN's execution, job and permission
paths. Detection or input preparation alone must not report a completed docking
run. Validate with redocking against a reference ligand, retaining every failed
pose and the original coordinate frame.

### Structure prediction execution

Keep the current detection step as a preflight. Add separate engine-specific
adapters with version-checked input schemas, model/weight availability checks,
resource declarations and explicit run provenance. Verify the actual command
and accepted input format for each supported engine version; a generic FASTA
command template is not an integration test. Return predicted structures and
their confidence artifacts as ordinary ports for subsequent analysis.

### Reporting and the connection to quantum chemistry

Combine the input structure, preparation choices, binding site, poses, metrics
and warnings into machine-readable and human-readable reports. A docking score
must keep its own meaning and must not be relabelled as a measured binding free
energy. An optional later path can extract a documented pocket or cluster for
QM/QM-MM work, preserving atom mapping, capping choices, charge and multiplicity.

## Integration check on 18 November 2026

The draft was merged locally with main at
`03a02d048f674a972d6c2e2551908eb5ece24b6c`, without conflicts. Revalidation
uses Python 3.12 and the actual merged tree; the older PR description's test
count does not describe the expanded draft.

The review reproduced a report-name path escape in `structure_confidence` and
`structure_prediction_engine`. Both now validate a plain filename before a
runner creates its work directory, and direct adapter execution returns a named
parameter failure for the same input. Regression tests also retain default and
custom report names. This check addresses the filename parameter; it is not a
replacement for host filesystem isolation.

No predictor, docking program, model download or live chemistry job was launched
as part of the integration review. The PR remains a draft.

The selected post-merge suite completed with **390 passed and 2 failed** on
Python 3.12. The two failing checks,
`test_declared_required_params_are_enforced` and
`test_every_adapter_has_a_contract` in `tests/test_tool_contracts.py`, also fail
on the unchanged main revision named above (141 passed, 2 failed when running
the same two tool-layer/contract files). They are not reported as new failures
of this draft and are not fixed as part of this update.

The report-name regression file failed in 24 cases before the fix (the four
default/custom-name controls passed) and all 28 cases pass with the fix.
Thirteen additional fake-model chat integration cases could not start their
local HTTP endpoint because this execution sandbox refuses socket creation;
they do not establish a browser or live-chat verification result.
