"""Content-checked, relocatable artifact manifests for biomolecular jobs."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path):
    def reject(value):
        raise ValueError(f"non-finite JSON number: {value}")
    def finite_float(value):
        number = float(value)
        if not math.isfinite(number):
            reject(value)
        return number
    return json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=reject, parse_float=finite_float)


def write_json(path: Path, value) -> None:
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def child_path(root: Path, relative: str) -> Path:
    if not isinstance(relative, str) or not relative or Path(relative).is_absolute():
        raise ValueError("artifact path must be a nonempty relative path")
    target = (Path(root) / relative).resolve()
    if not target.is_relative_to(Path(root).resolve()):
        raise ValueError("artifact path escapes its bundle")
    return target


def validate_artifact(path: Path, kind: str) -> dict:
    path = Path(path).resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"{kind}: missing or empty file: {path}")
    details = {}
    if kind in {"macromolecule", "prepared_receptor", "complex_structure"}:
        from delfin.tools.bio import read_structure
        structure = read_structure(path)
        if not structure.first_model.atoms:
            raise ValueError("structure has no atoms")
        if any(not all(math.isfinite(v) for v in atom.position)
               for model in structure.models for atom in model.atoms):
            raise ValueError("structure has non-finite coordinates")
        details["models"] = len(structure.models)
    elif kind in {"ligand", "ligand_library", "pose_set", "ligand_state_set"}:
        from delfin.tools.bio.processing import load_ligands
        records = load_ligands(path)
        if kind == "ligand" and len(records) != 1:
            raise ValueError("ligand input must contain exactly one molecule")
        if kind == "pose_set" and any(not mol.GetNumConformers() or not mol.GetConformer().Is3D() for mol in records):
            raise ValueError("pose_set requires existing 3D coordinates for every record")
        details["molecules"] = len(records)
    elif kind in {"binding_box", "binding_site"}:
        from delfin.tools.bio import BindingBox, binding_site_from_dict
        value = read_json(path)
        (BindingBox.from_dict(value) if kind == "binding_box" else binding_site_from_dict(value))
    elif kind == "assay_table":
        from delfin.tools.bio.assays import validate_assays
        details["measurements"] = len(validate_assays(read_json(path)))
    elif kind == "prediction_bundle":
        bundle = read_json(path)
        if bundle.get("schema_version") != 1 or not isinstance(bundle.get("models"), list) or not bundle["models"]:
            raise ValueError("invalid prediction bundle")
        ids = set()
        for model in bundle["models"]:
            if not isinstance(model.get("id"), str) or model["id"] in ids:
                raise ValueError("prediction model IDs must be unique strings")
            ids.add(model["id"])
            artifact = child_path(path.parent, model["path"])
            validate_artifact(artifact, "macromolecule")
            if file_hash(artifact) != model["sha256"]:
                raise ValueError("prediction model checksum mismatch")
            if "confidence" in model:
                confidence = model["confidence"]
                asset = child_path(path.parent, confidence["path"])
                read_json(asset)
                if file_hash(asset) != confidence["sha256"]:
                    raise ValueError("prediction confidence checksum mismatch")
    elif kind in {"receptor_pdbqt", "ligand_pdbqt", "poses_pdbqt"}:
        atoms = 0
        for line in path.read_text().splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                values = [float(line[i:i + 8]) for i in (30, 38, 46)]
                if not all(math.isfinite(v) for v in values):
                    raise ValueError("PDBQT has non-finite coordinates")
                atoms += 1
        if not atoms:
            raise ValueError("PDBQT contains no atoms")
        details["atom_records"] = atoms
    elif path.suffix.lower() == ".json":
        read_json(path)
    return {"type": kind, "path": str(path), "sha256": file_hash(path),
            "bytes": path.stat().st_size, **details}
