"""Typed biomolecular planning contracts; catalog entries are not executors.

No imports of external engines, subprocesses, downloads or filesystem writes.
Plans use explicit named edges so branching never means 'the last ligand'.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import os
import re
import shutil
from typing import Mapping


@dataclass(frozen=True)
class Backend:
    id: str
    task: str
    inputs: tuple[tuple[str, str], ...]
    outputs: tuple[tuple[str, str], ...]
    requirements: tuple[str, ...]
    documentation: str
    scope: str
    commands: tuple[str, ...] = ()
    adapter: str | None = None

    def describe(self, env: Mapping[str, str] | None = None) -> dict:
        row = asdict(self)
        row["inputs"] = dict(self.inputs)
        row["outputs"] = dict(self.outputs)
        row["integration"] = "adapter" if self.adapter else "planned"
        row["execution_ready"] = False
        row["readiness_note"] = "Requires separate environment, version, input and resource validation"
        row["command_paths"] = None if env is None else {
            cmd: shutil.which(cmd, path=env.get("PATH", "")) for cmd in self.commands
        }
        row["version_verified"] = False
        row["requirements_verified"] = False
        return row


def _backend(id, task, inputs, outputs, requirements, url, scope, commands=(), adapter=None):
    from delfin.tools.bio.execution import BINDINGS
    if adapter is None and id in BINDINGS:
        adapter = BINDINGS[id].adapter
    return Backend(id, task, tuple(inputs.items()), tuple(outputs.items()),
                   tuple(requirements), url, scope, tuple(commands), adapter)


_GH = "https://github.com/"
_PREDICTORS = (
    ("alphafold2", "google-deepmind/alphafold", "Protein monomer; FASTA and local databases", ()),
    ("alphafold_multimer", "google-deepmind/alphafold", "Protein complexes; multimer model preset", ()),
    ("alphafold3", "google-deepmind/alphafold3", "Protein/DNA/RNA/ligand complexes; AF3 JSON, not generic FASTA", ()),
    ("colabfold", "sokrypton/ColabFold", "Protein monomer/multimer; explicit MSA server or local MSA policy", ("colabfold_batch",)),
    ("esmfold", "facebookresearch/esm", "Protein sequence inference; no ligand cofolding", ("esm-fold",)),
    ("boltz2", "jwohlwend/boltz", "Biomolecular cofolding; YAML entities; model-specific affinity separate from confidence", ("boltz",)),
    ("chai1", "chaidiscovery/chai-lab", "Biomolecular cofolding; typed FASTA; explicit MSA/download policy", ("chai-lab",)),
    ("rosettafold_all_atom", "baker-laboratory/RoseTTAFold-All-Atom", "All-atom assemblies; engine-specific Hydra configuration", ()),
)

BACKENDS = tuple(
    _backend(id, "prediction", {"request": id + "_input"},
             {"models": "prediction_bundle"},
             ("versioned runner", "model weights", "compute resources", "MSA/template policy", "code and model terms"),
             _GH + repo, scope, commands)
    for id, repo, scope, commands in _PREDICTORS
) + (
    _backend("rdkit_prepare", "preparation", {"ligand": "ligand"}, {"ligand": "ligand"},
             ("rdkit",), "https://www.rdkit.org/docs/", "Coordinates, hydrogens and charges; protonation as supplied",
             adapter="ligand_prepare"),
    _backend("ligand_states", "preparation", {"ligand": "ligand"}, {"states": "ligand_state_set"},
             ("rdkit", "enumeration limits", "stereochemistry policy"),
             "https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html",
             "Executable RDKit tautomer enumeration; pH-dependent states require a separate method"),
    _backend("molscrub", "preparation", {"ligand": "ligand"}, {"states": "ligand_state_set"},
             ("molscrub", "pH range", "enumeration limits"), _GH + "forlilab/molscrub",
             "Ligand protonation/tautomer and coordinate preparation", ("scrub.py",)),
    _backend("receptor_prepare", "preparation", {"structure": "macromolecule"}, {"receptor": "prepared_receptor"},
             ("pdbfixer/openmm versions", "single canonical PDB model without hydrogens/altlocs", "explicit pH/water/repair policy"),
             "https://github.com/openmm/pdbfixer", "Restricted PDBFixer repair and OpenMM hydrogen assignment; nonstandard chemistry and missing residues are refused"),
    _backend("meeko_receptor", "preparation", {"receptor": "prepared_receptor"}, {"receptor": "receptor_pdbqt"},
             ("meeko", "validated receptor chemistry"), "https://meeko.readthedocs.io/",
             "Receptor typing for AutoDock; does not replace biological preparation", ("mk_prepare_receptor.py",)),
    _backend("meeko_ligand", "preparation", {"ligand": "ligand"}, {"ligand": "ligand_pdbqt"},
             ("meeko", "explicit hydrogens", "atom mapping"), "https://meeko.readthedocs.io/",
             "Ligand torsions and PDBQT with recoverable chemical identity", ("mk_prepare_ligand.py",)),
    _backend("vina", "docking", {"receptor": "receptor_pdbqt", "ligand": "ligand_pdbqt", "site": "binding_box"},
             {"poses": "poses_pdbqt"}, ("vina", "seed", "search effort", "CPU allocation"),
             "https://autodock-vina.readthedocs.io/en/latest/docking_basic.html",
             "Rigid-receptor small-molecule docking; Vina scores are not measured affinity", ("vina",)),
    _backend("meeko_export", "conversion", {"poses": "poses_pdbqt"}, {"poses": "pose_set"},
             ("meeko", "embedded chemical mapping"), "https://meeko.readthedocs.io/",
             "Recover SDF chemistry from Meeko/Vina PDBQT", ("mk_export.py",)),
    _backend("gnina", "docking", {"receptor": "prepared_receptor", "ligand": "ligand", "site": "binding_box"},
             {"poses": "pose_set"}, ("gnina", "CNN model", "CPU/GPU policy"), _GH + "gnina/gnina",
             "Docking/CNN rescoring; keep affinity and CNN metrics distinct", ("gnina",)),
    _backend("smina", "docking", {"receptor": "prepared_receptor", "ligand": "ligand", "site": "binding_box"},
             {"poses": "pose_set"}, ("smina", "scoring function"), "https://sourceforge.net/projects/smina/",
             "Small-molecule docking and scoring", ("smina",)),
    _backend("diffdock", "docking", {"receptor": "prepared_receptor", "ligand": "ligand"},
             {"poses": "pose_set"}, ("checkpoint", "GPU environment", "sampling controls"), _GH + "gcorso/DiffDock",
             "Learned pose generation; ranking confidence is not binding affinity"),
    _backend("haddock3", "docking", {"partners": "partner_structures", "restraints": "restraints"},
             {"complexes": "complex_ensemble"}, ("HADDOCK configuration", "partner preparation"),
             _GH + "haddocking/haddock3", "Integrative docking of prepared partners", ("haddock3",)),
    _backend("glide", "docking", {"receptor": "prepared_receptor", "ligand": "ligand", "site": "binding_box"},
             {"poses": "pose_set"}, ("licensed installation", "receptor grid", "suite-specific conversion"),
             "https://www.schrodinger.com/platform/products/glide/", "Commercial backend extension point"),
    _backend("gold", "docking", {"receptor": "prepared_receptor", "ligand": "ligand", "site": "binding_box"},
             {"poses": "pose_set"}, ("licensed installation", "GOLD configuration"),
             "https://www.ccdc.cam.ac.uk/solutions/software/gold/", "Commercial backend extension point"),
    _backend("posebusters", "validation", {"poses": "pose_set", "reference": "ligand", "receptor": "prepared_receptor"},
             {"report": "pose_validation"}, ("posebusters", "validation thresholds"), _GH + "maabuu/posebusters",
             "Chemical and geometric plausibility; retain failed poses", ("bust",)),
    _backend("prolif", "analysis", {"poses": "pose_set", "receptor": "prepared_receptor"},
             {"interactions": "interaction_report"}, ("prolif", "bond orders/hydrogens"),
             "https://prolif.readthedocs.io/", "Protein-ligand interaction fingerprints"),
    _backend("plip", "analysis", {"complex": "complex_structure"}, {"interactions": "interaction_report"},
             ("plip", "complex assembly with atom mapping"), _GH + "pharmai/plip",
             "Interaction analysis of an assembled complex", ("plip",)),
    _backend("medchem_properties", "medchem", {"ligands": "ligand_library"}, {"properties": "property_table"},
             ("rdkit", "compound IDs", "descriptor definitions"), "https://www.rdkit.org/docs/",
             "Computed descriptors and Murcko scaffolds; no automatic ADMET validity claim"),
    _backend("assay_join", "medchem", {"properties": "property_table", "assays": "assay_table"},
             {"dataset": "medchem_dataset"}, ("compound/state identity", "assay units/conditions", "measurement provenance"),
             "https://www.ebi.ac.uk/chembl/", "Join measured endpoints without conflating IC50, Ki and docking scores"),
    _backend("select_model", "conversion", {"models": "prediction_bundle"}, {"structure": "macromolecule"},
             ("explicit model/rank", "chain mapping", "confidence provenance"),
             _GH + "ComPlat/DELFIN", "Choose one predicted structure; cofolded ligands require explicit extraction"),
    _backend("assemble_complex", "conversion", {"receptor": "prepared_receptor", "poses": "pose_set"},
             {"complex": "complex_structure"}, ("explicit pose ID", "common coordinate frame", "atom mapping"),
             _GH + "ComPlat/DELFIN", "Assemble a selected pose and receptor without coordinate realignment"),
    _backend("parameterize", "simulation", {"complex": "complex_structure"}, {"system": "simulation_system"},
             ("OpenMM version", "canonical hydrogenated protein", "explicit ff14SB/vacuum/constraints/charge policies"),
             "https://docs.openmm.org/latest/userguide/", "Checked system bundle; current execution subset is protein-only ff14SB in vacuum, without ligand/solvent parameterization"),
    *tuple(_backend(id, "simulation", {"system": "simulation_system"}, {"trajectory": "trajectory", "checkpoint": "checkpoint"},
                    ("engine-specific system conversion", "integrator/ensemble", "resources", "restart state"), url,
                    "Molecular dynamics; engine-specific parameters and trajectory provenance", commands)
           for id, url, commands in (
               ("openmm", "https://docs.openmm.org/latest/userguide/", ()),
               ("gromacs", "https://manual.gromacs.org/current/", ("gmx",)),
               ("amber", "https://ambermd.org/", ("sander",)),
               ("namd", "https://www.ks.uiuc.edu/Research/namd/", ("namd3",)),
           )),
    _backend("mdanalysis", "analysis", {"trajectory": "trajectory", "system": "simulation_system"},
             {"report": "trajectory_report"}, ("MDAnalysis", "topology/trajectory atom mapping"),
             "https://docs.mdanalysis.org/", "Trajectory analysis with topology, units and periodic boundaries"),
    _backend("qmmm_region", "qm", {"complex": "complex_structure"}, {"region": "qm_region"},
             ("atom mapping", "boundary/caps", "charge/multiplicity", "embedding"),
             _GH + "ComPlat/DELFIN", "Bridge to existing ORCA QM/MM; a protein complex is not an XYZ-ready QM region"),
    _backend("rfdiffusion", "design", {"design": "design_request"}, {"backbones": "backbone_ensemble"},
             ("checkpoint", "contigs/constraints", "GPU environment"), _GH + "RosettaCommons/RFdiffusion",
             "Backbone generation; does not establish binding or experimental function"),
    _backend("proteinmpnn", "design", {"backbone": "macromolecule"}, {"sequences": "sequence_set"},
             ("checkpoint", "fixed/designed positions", "sampling controls"), _GH + "dauparas/ProteinMPNN",
             "Sequence design conditioned on a selected backbone"),
    _backend("structure_select", "processing", {"structure": "macromolecule"}, {"structure": "macromolecule"},
             ("explicit model/chains/altloc policy",), _GH + "ComPlat/DELFIN", "Coordinate selection with atom mapping, not receptor preparation"),
    _backend("ligand_depict", "visualization", {"ligands": "ligand_library"}, {"image": "svg"},
             ("rdkit",), "https://www.rdkit.org/docs/", "Offline 2D structure grid"),
    _backend("visualize", "visualization", {"structure": "macromolecule"}, {"view": "html_view", "manifest": "view_manifest"},
             ("browser WebGL for viewing",), "https://3dmol.csb.pitt.edu/", "Offline 3Dmol viewer plus PyMOL and ChimeraX entry scripts"),
    _backend("visualize_poses", "visualization", {"structure": "macromolecule", "poses": "pose_set", "site": "binding_box"},
             {"view": "html_view", "manifest": "view_manifest"}, ("shared coordinate frame", "browser WebGL"),
             "https://3dmol.csb.pitt.edu/", "Explicit receptor/poses/box wiring into an offline viewer"),
    _backend("ligand_select", "processing", {"ligands": "ligand_library"}, {"ligand": "ligand"},
             ("explicit record_id",), "https://www.rdkit.org/docs/", "Select one ligand without changing chemistry"),
    _backend("pose_select", "processing", {"poses": "pose_set"}, {"ligand": "ligand"},
             ("explicit record_id",), "https://www.rdkit.org/docs/", "Select one pose without alignment or score guessing"),
    _backend("state_select", "processing", {"states": "ligand_state_set"}, {"ligand": "ligand"},
             ("explicit record_id",), "https://www.rdkit.org/docs/", "Select a chemical state; regenerate 3D before docking"),
    _backend("box_from_ligand", "processing", {"ligand": "ligand"}, {"site": "binding_box"},
             ("3D ligand in receptor coordinate frame",), _GH + "ComPlat/DELFIN", "Heavy-atom bounding box plus padding in angstrom"),
)
_INDEX = {backend.id: backend for backend in BACKENDS}
PORT_TYPES = frozenset(kind for b in BACKENDS for _, kind in (*b.inputs, *b.outputs))


def backend_catalog(*, task: str | None = None, probe: bool = False,
                    env: Mapping[str, str] | None = None) -> dict:
    from delfin.tools.bio.schemas import plan_schema, prediction_request_schema, assay_table_schema
    if task is not None and task not in {b.task for b in BACKENDS}:
        raise ValueError(f"unknown biomolecular task: {task}")
    environment = (os.environ if env is None else env) if probe else None
    return {"schema_version": 1, "kind": "biomolecular_backend_catalog",
            "execution_performed": False,
            "schemas": {"plan": plan_schema(), "prediction_request": prediction_request_schema(), "assay_table": assay_table_schema()},
            "backends": [b.describe(environment) for b in BACKENDS if task is None or b.task == task]}


def validate_plan(plan: dict) -> dict:
    """Validate an explicitly ordered DAG. Never convert it to a runnable job.

    Schema: external inputs have semantic types; each node binds every input
    to 'inputs.ID' or 'NODE.PORT'. A producer must precede its consumers.
    Missing implementations are blockers, not structural errors, so agents
    can save and improve a well-formed plan before engines are integrated.
    """
    errors, blockers, steps = [], [], []
    if not isinstance(plan, dict) or type(plan.get("schema_version")) is not int or plan["schema_version"] != 1:
        return {"valid": False, "execution_ready": False, "errors": ["schema_version must be 1"], "blockers": [], "steps": []}
    if set(plan) - {"schema_version", "name", "inputs", "nodes"}:
        errors.append("unknown plan fields")
    inputs, nodes = plan.get("inputs", {}), plan.get("nodes", [])
    if not isinstance(inputs, dict) or not isinstance(nodes, list) or not nodes:
        return {"valid": False, "execution_ready": False, "errors": ["inputs must be an object and nodes a nonempty list"], "blockers": [], "steps": []}
    available = {}
    ident = re.compile(r"[A-Za-z][A-Za-z0-9_]*\Z")
    for name, item in inputs.items():
        if not isinstance(name, str) or not ident.fullmatch(name) or not isinstance(item, dict) or not isinstance(item.get("type"), str) or item["type"] not in PORT_TYPES:
            errors.append(f"invalid external input {name!r}")
            continue
        if set(item) - {"type", "path", "description"}:
            errors.append(f"unknown fields in input {name}")
        available[f"inputs.{name}"] = item["type"]
        if not isinstance(item.get("path"), str) or not item["path"].strip():
            blockers.append(f"inputs.{name}: bind an actual input artifact before execution")
    seen = {"inputs"}
    for index, node in enumerate(nodes):
        if not isinstance(node, dict):
            errors.append(f"node {index} must be an object")
            continue
        name, backend_id = node.get("id"), node.get("backend")
        if not isinstance(name, str) or not ident.fullmatch(name) or name in seen:
            errors.append(f"invalid or duplicate node id {name!r}")
            continue
        seen.add(name)
        if not isinstance(backend_id, str) or backend_id not in _INDEX:
            errors.append(f"{name}: unknown backend {backend_id!r}")
            continue
        backend = _INDEX[backend_id]
        if set(node) - {"id", "backend", "inputs", "settings"}:
            errors.append(f"{name}: unknown node fields")
        bindings = node.get("inputs", {})
        if not isinstance(bindings, dict):
            errors.append(f"{name}: inputs must be an object")
            continue
        expected = dict(backend.inputs)
        for port in sorted(set(bindings) - set(expected)):
            errors.append(f"{name}: unknown input port {port}")
        for port, kind in expected.items():
            ref = bindings.get(port)
            actual = available.get(ref) if isinstance(ref, str) else None
            if actual is None:
                errors.append(f"{name}.{port}: missing/forward reference {ref!r}; expected {kind}")
            elif actual != kind:
                errors.append(f"{name}.{port}: expected {kind}, received {actual}; explicit conversion required")
        settings = node.get("settings", {})
        if not isinstance(settings, dict):
            errors.append(f"{name}: settings must be an object")
        if backend.adapter is None:
            blockers.append(f"{name}: {backend.id} execution adapter not implemented")
        blockers.extend(f"{name}: verify {requirement}" for requirement in backend.requirements)
        steps.append({"id": name, "backend": backend.id, "adapter": backend.adapter,
                      "integration": "adapter" if backend.adapter else "planned",
                      "settings_validation": "not_performed", "outputs": dict(backend.outputs)})
        for port, kind in backend.outputs:
            available[f"{name}.{port}"] = kind
    return {"valid": not errors, "execution_ready": False, "execution_performed": False,
            "errors": errors, "blockers": blockers, "steps": steps,
            "note": "Planning graph only; external paths, settings, resources and engine versions are not execution-validated"}


def workflow_template(name: str) -> dict:
    if name == "processing":
        return {"schema_version": 1, "name": name,
                "inputs": {"structure": {"type": "macromolecule"}, "ligands": {"type": "ligand_library"}},
                "nodes": [
                    {"id": "select", "backend": "structure_select", "inputs": {"structure": "inputs.structure"}},
                    {"id": "view", "backend": "visualize", "inputs": {"structure": "select.structure"}},
                    {"id": "properties", "backend": "medchem_properties", "inputs": {"ligands": "inputs.ligands"}},
                    {"id": "depict", "backend": "ligand_depict", "inputs": {"ligands": "inputs.ligands"}},
                ]}
    if name == "docking":
        return {"schema_version": 1, "name": name,
                "inputs": {"structure": {"type": "macromolecule"}, "ligand": {"type": "ligand"}, "site": {"type": "binding_box"}},
                "nodes": [
                    {"id": "prepare_receptor", "backend": "receptor_prepare", "inputs": {"structure": "inputs.structure"}},
                    {"id": "prepare_ligand", "backend": "rdkit_prepare", "inputs": {"ligand": "inputs.ligand"}},
                    {"id": "receptor", "backend": "meeko_receptor", "inputs": {"receptor": "prepare_receptor.receptor"}},
                    {"id": "ligand", "backend": "meeko_ligand", "inputs": {"ligand": "prepare_ligand.ligand"}},
                    {"id": "dock", "backend": "vina", "inputs": {"receptor": "receptor.receptor", "ligand": "ligand.ligand", "site": "inputs.site"}},
                    {"id": "export", "backend": "meeko_export", "inputs": {"poses": "dock.poses"}},
                    {"id": "validate", "backend": "posebusters", "inputs": {"poses": "export.poses", "reference": "prepare_ligand.ligand", "receptor": "prepare_receptor.receptor"}},
                    {"id": "contacts", "backend": "prolif", "inputs": {"poses": "export.poses", "receptor": "prepare_receptor.receptor"}},
                ]}
    if name in {b.id for b in BACKENDS if b.task == "prediction"}:
        return {"schema_version": 1, "name": name,
                "inputs": {"request": {"type": name + "_input"}},
                "nodes": [{"id": "predict", "backend": name, "inputs": {"request": "inputs.request"}},
                          {"id": "select", "backend": "select_model", "inputs": {"models": "predict.models"}}]}
    if name == "medchem":
        return {"schema_version": 1, "name": name,
                "inputs": {"ligands": {"type": "ligand_library"}, "assays": {"type": "assay_table"}},
                "nodes": [{"id": "properties", "backend": "medchem_properties", "inputs": {"ligands": "inputs.ligands"}},
                          {"id": "join", "backend": "assay_join", "inputs": {"properties": "properties.properties", "assays": "inputs.assays"}}]}
    if name == "md":
        return {"schema_version": 1, "name": name, "inputs": {"complex": {"type": "complex_structure"}},
                "nodes": [{"id": "parameters", "backend": "parameterize", "inputs": {"complex": "inputs.complex"}},
                          {"id": "simulate", "backend": "openmm", "inputs": {"system": "parameters.system"}},
                          {"id": "analyze", "backend": "mdanalysis", "inputs": {"system": "parameters.system", "trajectory": "simulate.trajectory"}}]}
    raise ValueError(f"unknown biomolecular template: {name}")
