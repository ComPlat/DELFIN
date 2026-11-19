"""Build native prediction input files without executing prediction engines.

The common request intentionally covers unmodified polymers and SMILES ligands
only. Unknown fields are refused rather than silently losing restraints,
modifications or covalent bonds during conversion.
"""
from __future__ import annotations

import json
import re


EXPORTERS = ("alphafold2", "alphafold_multimer", "alphafold3", "colabfold", "esmfold", "boltz2", "chai1")


def validate_request(request: dict) -> dict:
    if not isinstance(request, dict) or type(request.get("schema_version")) is not int or request["schema_version"] != 1:
        raise ValueError("prediction request schema_version must be 1")
    unknown = set(request) - {"schema_version", "name", "entities", "seeds", "msa_policy"}
    if unknown:
        raise ValueError(f"unsupported request fields: {sorted(unknown)}")
    name = request.get("name")
    if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{0,79}", name):
        raise ValueError("name must be a simple identifier (1–80 characters)")
    seeds = request.get("seeds")
    if not isinstance(seeds, list) or not seeds or any(type(s) is not int or not 0 <= s < 2**31 for s in seeds):
        raise ValueError("seeds must be a nonempty list of integers in [0, 2**31)")
    if len(set(seeds)) != len(seeds):
        raise ValueError("seeds must be unique")
    if request.get("msa_policy") not in ("single_sequence", "engine_pipeline"):
        raise ValueError("explicit msa_policy required: single_sequence or engine_pipeline")
    entities = request.get("entities")
    if not isinstance(entities, list) or not entities:
        raise ValueError("entities must be a nonempty list")
    ids = set()
    alphabets = {"protein": set("ACDEFGHIKLMNPQRSTVWYX"), "dna": set("ACGTN"), "rna": set("ACGUN")}
    for item in entities:
        if not isinstance(item, dict) or set(item) - {"id", "kind", "sequence", "smiles"}:
            raise ValueError("unsupported entity fields; modifications/restraints need a dedicated converter")
        id = item.get("id")
        if not isinstance(id, str) or not re.fullmatch(r"[A-Z]+", id) or id in ids:
            raise ValueError("entity IDs must be unique uppercase letters")
        ids.add(id)
        kind = item.get("kind")
        if kind == "ligand":
            value = item.get("smiles")
            if "sequence" in item or not isinstance(value, str) or not value or any(c.isspace() for c in value):
                raise ValueError("ligand requires SMILES, without sequence or whitespace")
            from rdkit import Chem
            mol = Chem.MolFromSmiles(value)
            if mol is None or mol.GetNumAtoms() == 0:
                raise ValueError("invalid ligand SMILES")
        elif isinstance(kind, str) and kind in alphabets:
            value = item.get("sequence")
            if "smiles" in item or not isinstance(value, str) or not value or set(value) - alphabets[kind]:
                raise ValueError(f"invalid {kind} sequence; explicit unmodified uppercase sequence required")
        else:
            raise ValueError(f"unsupported entity kind: {kind!r}")
    return request


def build_prediction_input(backend: str, request: dict) -> dict:
    validate_request(request)
    if backend not in EXPORTERS:
        raise ValueError(f"no input exporter implemented for {backend!r}")
    entities = request["entities"]
    policy = request["msa_policy"]
    notes = ["Input conversion only; engine/version, weights, resources and model terms remain unchecked",
             "Seeds are recorded in the manifest; runners must apply them where not encoded in the native file"]
    if backend == "alphafold3":
        rows = []
        for e in entities:
            value = {"id": e["id"]}
            value["smiles" if e["kind"] == "ligand" else "sequence"] = e.get("smiles", e.get("sequence"))
            if policy == "single_sequence":
                if e["kind"] == "protein":
                    value.update(unpairedMsa="", pairedMsa="", templates=[])
                elif e["kind"] == "rna":
                    value["unpairedMsa"] = ""
            rows.append({e["kind"]: value})
        payload = {"name": request["name"], "modelSeeds": request["seeds"],
                   "sequences": rows, "dialect": "alphafold3", "version": 1}
        filename, content, dialect = "input.json", json.dumps(payload, indent=2), "alphafold3-v1-basic"
    elif backend == "boltz2":
        if policy != "single_sequence":
            raise ValueError("Boltz exporter currently requires explicit single_sequence; local/server MSAs need a dedicated converter")
        rows = []
        for e in entities:
            value = {"id": e["id"], "smiles" if e["kind"] == "ligand" else "sequence": e.get("smiles", e.get("sequence"))}
            if e["kind"] == "protein":
                value["msa"] = "empty"
            rows.append({e["kind"]: value})
        filename = "input.yaml"
        content = json.dumps({"version": 1, "sequences": rows}, indent=2)
        dialect = "boltz-yaml-v1 (JSON subset of YAML)"
    elif backend == "chai1":
        if policy != "single_sequence":
            raise ValueError("Chai exporter currently requires explicit single_sequence; MSAs need a dedicated converter")
        filename = "input.fasta"
        content = "".join(f">{e['kind']}|name={e['id']}\n{e.get('sequence', e.get('smiles'))}\n" for e in entities)
        dialect = "chai-typed-fasta"
    else:
        if any(e["kind"] != "protein" for e in entities):
            raise ValueError(f"{backend} exporter accepts protein entities only")
        if backend != "alphafold_multimer" and len(entities) != 1:
            raise ValueError(f"{backend} exporter accepts one protein; multimer/batch semantics need an explicit converter")
        if backend in ("alphafold2", "alphafold_multimer") and policy != "engine_pipeline":
            raise ValueError("AlphaFold 2 FASTA export requires engine_pipeline MSA policy")
        if backend == "esmfold" and policy != "single_sequence":
            raise ValueError("ESMFold uses single_sequence")
        filename = "input.fasta"
        content = "".join(f">{e['id']}\n{e['sequence']}\n" for e in entities)
        dialect = "protein-fasta"
        notes.append("FASTA does not encode model preset or MSA policy; runner must apply the manifest")
    return {"filename": filename, "content": content.rstrip() + "\n",
            "backend": backend, "artifact_type": backend + "_input", "dialect": dialect,
            "request": request, "execution_ready": False, "execution_performed": False,
            "notes": notes}
