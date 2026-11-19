"""Coordinate selection and explicit prediction-result import."""
from __future__ import annotations

from pathlib import Path
import shutil

from delfin.tools.bio import read_structure
from delfin.tools.bio.artifacts import child_path, file_hash, read_json, validate_artifact, write_json


def _token(value) -> str:
    text = str(value)
    if not text:
        return "."
    if any(c in text for c in "\r\n"):
        raise ValueError("multiline coordinate identifiers are unsupported")
    reserved = text.startswith(("_", "$", "[", "]")) or text.lower().startswith(("data_", "save_")) or text.lower() in {"loop_", "stop_", "global_"}
    if not reserved and not any(c.isspace() for c in text) and not any(c in text for c in "'\"#;") and text not in {".", "?"}:
        return text
    if "'" not in text:
        return "'" + text + "'"
    if '"' not in text:
        return '"' + text + '"'
    raise ValueError("coordinate identifier contains both quote types")


def select_structure(source: Path, output: Path, *, model=None, chains=None,
                     waters="keep", hetero="keep", altloc="all") -> dict:
    if waters not in {"keep", "drop"} or hetero not in {"keep", "drop"}:
        raise ValueError("waters and hetero must be keep or drop")
    if altloc not in {"all", "primary"}:
        raise ValueError("altloc must be all or primary")
    if model is not None and type(model) is not int:
        raise ValueError("model must be an integer")
    if chains is not None and (not isinstance(chains, list) or not all(isinstance(c, str) for c in chains)):
        raise ValueError("chains must be a list of chain IDs")
    structure = read_structure(source)
    chosen = structure.first_model if model is None else structure.model(model)
    if chains is not None and set(chains) - {c.chain_id for c in chosen.chains}:
        raise ValueError("selected chain is absent from the chosen model")
    headers = ("group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
               "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z",
               "occupancy", "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num")
    lines = ["data_delfin_selection", "#", "loop_", *["_atom_site." + name for name in headers]]
    selected, removed = [], []
    for residue in chosen.residues:
        reason = None
        if chains is not None and residue.chain_id not in chains:
            reason = "chain"
        elif residue.is_water and waters == "drop":
            reason = "water"
        elif residue.is_hetero and not residue.is_modified and not residue.is_water and hetero == "drop":
            reason = "hetero"
        atoms = residue.primary_atoms if altloc == "primary" else residue.atoms
        chosen_ids = {id(atom) for atom in atoms}
        for atom in residue.atoms:
            identity = {"model": chosen.model_id, "chain": residue.chain_id, "residue": residue.seq_id,
                        "insertion_code": residue.insertion_code, "residue_name": residue.original_name or residue.name,
                        "atom_name": atom.name, "serial": atom.serial, "altloc": atom.altloc}
            why = reason or ("altloc" if id(atom) not in chosen_ids else None)
            if why:
                removed.append({**identity, "reason": why})
                continue
            number = len(selected) + 1
            row = ["HETATM" if atom.is_hetero else "ATOM", number, atom.element, atom.name,
                   atom.altloc, residue.original_name or residue.name, residue.chain_id, 1, residue.seq_id,
                   residue.insertion_code, atom.x, atom.y, atom.z, atom.occupancy, atom.b_factor,
                   residue.seq_id, residue.original_name or residue.name, residue.chain_id, atom.name, chosen.model_id]
            lines.append(" ".join(_token(value) for value in row))
            selected.append({"output_serial": number, **identity})
    if not selected:
        raise ValueError("selection contains no atoms")
    lines.append("#")
    Path(output).write_text("\n".join(lines) + "\n", encoding="utf-8")
    validate_artifact(output, "macromolecule")
    return {"schema_version": 1, "source_sha256": file_hash(source), "source_model": chosen.model_id,
            "selected": selected, "removed": removed,
            "other_models_omitted": [m.model_id for m in structure.models if m.model_id != chosen.model_id],
            "operation": "coordinate_selection_only", "receptor_prepared": False,
            "metadata_note": "Source connection, assembly, missing-residue and modification annotations remain in the original input; selected CIF contains coordinates only"}


def import_predictions(source_dir: Path, models: list, output_dir: Path, *, engine: str, engine_version: str) -> Path:
    if not isinstance(models, list) or not models:
        raise ValueError("models must explicitly list id/path for every model to import")
    if not isinstance(engine, str) or not engine.strip() or not isinstance(engine_version, str) or not engine_version.strip():
        raise ValueError("engine and engine_version are required provenance")
    prepared, ids = [], set()
    for entry in models:
        if not isinstance(entry, dict) or set(entry) - {"id", "path", "confidence"}:
            raise ValueError("model entry accepts id, path and optional confidence JSON path")
        identity = entry.get("id")
        if not isinstance(identity, str) or not identity or identity in ids:
            raise ValueError("model IDs must be nonempty unique strings")
        ids.add(identity)
        path = child_path(source_dir, entry["path"])
        details = validate_artifact(path, "macromolecule")
        confidence = None
        if entry.get("confidence"):
            confidence = child_path(source_dir, entry["confidence"])
            read_json(confidence)
        prepared.append((identity, path, details, confidence))
    rows = []
    for index, (identity, source, details, confidence) in enumerate(prepared, 1):
        suffix = "".join(source.suffixes)
        dest = output_dir / f"model_{index}{suffix}"
        shutil.copyfile(source, dest)
        row = {"id": identity, "path": dest.name, "sha256": details["sha256"], "source_path": str(source)}
        if confidence:
            target = output_dir / f"confidence_{index}.json"
            shutil.copyfile(confidence, target)
            row["confidence"] = {"path": target.name, "sha256": file_hash(target), "interpretation": "raw_engine_output"}
        rows.append(row)
    path = output_dir / "prediction_bundle.json"
    write_json(path, {"schema_version": 1, "engine": engine, "engine_version": engine_version,
                     "provenance_origin": "caller_supplied; engine not executed by import", "models": rows})
    return path


def select_prediction(bundle_path: Path, model_id: str, output_dir: Path) -> Path:
    validate_artifact(bundle_path, "prediction_bundle")
    bundle = read_json(bundle_path)
    matches = [m for m in bundle["models"] if m["id"] == model_id]
    if len(matches) != 1:
        raise ValueError("explicit model_id does not identify exactly one model")
    source = child_path(bundle_path.parent, matches[0]["path"])
    dest = output_dir / ("selected" + "".join(source.suffixes))
    shutil.copyfile(source, dest)
    write_json(output_dir / "selection.json", {"bundle_sha256": file_hash(bundle_path), "model_id": model_id,
                                               "selected_sha256": file_hash(dest), "receptor_prepared": False})
    return dest
