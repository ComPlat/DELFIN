"""Executable structure/ligand processing and visualization adapters."""
from __future__ import annotations

from pathlib import Path
import time

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import ErrorKind, StepStatus
from delfin.tools.bio.artifacts import read_json, write_json


class BioProcessing(StepAdapter):
    produces_geometry = False
    category = "structure"
    data_keys = (DataKeySpec("operation", "str"),)

    def __init__(self, name, description, params, outputs, wires=(), requires_python=(), example=None):
        self.name, self.description, self.params = name, description, tuple(params)
        self.produces = tuple(outputs)
        self.required_artifacts = self.produces
        self.wires = tuple(wires)
        self.consumes = tuple(cap for cap, _ in wires)
        self.requires_python = tuple(requires_python)
        self.example_kwargs = example

    def validate_params(self, **kwargs):
        for param in self.params:
            if param.required and kwargs.get(param.name) is None:
                raise ValueError(f"{param.name} is required")
        if self.name == "bio_ligand_select":
            value = kwargs.get("record_id")
            if type(value) is not int or value < 1:
                raise ValueError("record_id must be a positive 1-based integer")
        if self.name == "bio_assemble_complex":
            if type(kwargs.get("record_id")) is not int or kwargs["record_id"] < 1:
                raise ValueError("record_id must be a positive 1-based integer")
            if kwargs.get("shared_coordinate_frame") is not True:
                raise ValueError("shared_coordinate_frame=True must be explicitly asserted")
        if self.name == "bio_ligand_states":
            value = kwargs.get("max_states", 32)
            if type(value) is not int or not 1 <= value <= 1024:
                raise ValueError("max_states must be in 1..1024")
        if self.name == "bio_structure_select":
            for key in ("waters", "hetero"):
                if kwargs.get(key, "keep") not in {"keep", "drop"}:
                    raise ValueError(f"{key} must be keep or drop")
            if kwargs.get("altloc", "all") not in {"all", "primary"}:
                raise ValueError("altloc must be all or primary")

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start = time.monotonic()
        work_dir = Path(work_dir)
        try:
            self.validate_params(**kwargs)
            work_dir.mkdir(parents=True, exist_ok=True)
            if any(work_dir.iterdir()):
                raise ValueError("processing requires a fresh empty step directory")
            data, artifacts = self._run(work_dir, kwargs)
            data["operation"] = self.name
            report = work_dir / "processing_report.json"
            write_json(report, data)
            artifacts["json_report"] = report
            return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=data, artifacts=artifacts)
        except FileNotFoundError as exc:
            kind = ErrorKind.MISSING_INPUT
            error = str(exc)
        except (ImportError, ModuleNotFoundError) as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, f"Optional dependency unavailable: {exc}"
        except (ValueError, KeyError, TypeError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error=error, error_kind=kind)

    def _run(self, directory, options):
        from delfin.tools.bio import processing, structures, visualization
        if self.name == "bio_assemble_complex":
            from delfin.tools.bio.assembly import assemble_complex
            data = assemble_complex(Path(options["receptor_file"]), Path(options["poses_file"]), directory,
                                    record_id=options["record_id"], shared_coordinate_frame=options["shared_coordinate_frame"],
                                    ligand_chain=options.get("ligand_chain", "L"))
            return data, {"complex_structure": directory / "complex.pdb", "ligand": directory / "selected_ligand.sdf"}
        if self.name == "bio_assay_join":
            from delfin.tools.bio.assays import join_assays
            data = join_assays(Path(options["properties_file"]), Path(options["assays_file"]))
            path = directory / "dataset.json"
            write_json(path, data)
            return data, {"medchem_dataset": path}
        if self.name == "bio_ligand_select":
            path = directory / "selected.sdf"
            data = processing.select_ligand(Path(options["ligands_file"]), path, record_id=options["record_id"])
            return data, {"ligand": path}
        if self.name == "bio_box_from_ligand":
            path = directory / "box.json"
            data = processing.box_from_ligand(Path(options["ligand_file"]), padding=options.get("padding", 4.0))
            write_json(path, data)
            return {"box": data, "coordinate_frame": "input ligand; no alignment", "unit": "angstrom"}, {"binding_box": path}
        if self.name == "bio_ligand_properties":
            data = processing.ligand_properties(Path(options["ligands_file"]))
            path = directory / "properties.json"
            write_json(path, data)
            return data, {"property_table": path}
        if self.name == "bio_ligand_states":
            path = directory / "states.sdf"
            data = processing.enumerate_states(Path(options["ligand_file"]), path, max_states=options.get("max_states", 32))
            return data, {"ligand_state_set": path}
        if self.name == "bio_ligand_depict":
            path = directory / "ligands.svg"
            path.write_text(processing.ligand_grid_svg(Path(options["ligands_file"])), encoding="utf-8")
            return {"method": "RDKit 2D depiction; input coordinates unchanged"}, {"svg": path}
        if self.name == "bio_structure_select":
            path = directory / "selected.cif"
            data = structures.select_structure(Path(options["structure_file"]), path,
                                                **{k: options[k] for k in ("model", "chains", "waters", "hetero", "altloc") if k in options})
            return data, {"macromolecule": path}
        if self.name == "bio_prediction_import":
            path = structures.import_predictions(Path(options["source_dir"]), options["models"], directory,
                                                 engine=options["engine"], engine_version=options["engine_version"])
            return read_json(path), {"prediction_bundle": path}
        if self.name == "bio_model_select":
            path = structures.select_prediction(Path(options["bundle_file"]), options["model_id"], directory)
            return read_json(directory / "selection.json"), {"macromolecule": path}
        if self.name == "bio_visualize":
            data = visualization.export_view(Path(options["structure_file"]), directory,
                                             **{k: options[k] for k in ("poses_file", "site_file", "model") if k in options})
            return data, {"html_view": directory / "view.html", "view_manifest": directory / "view_manifest.json"}
        raise ValueError(f"unknown processing operation {self.name}")


register(BioProcessing("bio_assemble_complex", "Assemble one explicit pose and a single prepared PDB receptor, preserving atom mapping and ligand SDF",
                      [ParamSpec("receptor_file", "path", required=True), ParamSpec("poses_file", "path", required=True),
                       ParamSpec("record_id", "int", required=True), ParamSpec("shared_coordinate_frame", "bool", required=True),
                       ParamSpec("ligand_chain", "str", default="L")],
                      ["complex_structure", "ligand", "json_report"], [("prepared_receptor", "receptor_file"), ("pose_set", "poses_file")],
                      ["rdkit"], example={"receptor_file": "receptor.pdb", "poses_file": "poses.sdf", "record_id": 1, "shared_coordinate_frame": True}))
register(BioProcessing("bio_ligand_properties", "Compute RDKit descriptors and scaffolds for every ligand record",
                      [ParamSpec("ligands_file", "path", required=True)], ["property_table", "json_report"],
                      [("ligand_library", "ligands_file")], ["rdkit"]))
register(BioProcessing("bio_assay_join", "Join measured endpoints by exact chemical state; preserve units, inequalities, conditions, replicates and unmatched rows",
                      [ParamSpec("properties_file", "path", required=True), ParamSpec("assays_file", "path", required=True)],
                      ["medchem_dataset", "json_report"], [("property_table", "properties_file"), ("assay_table", "assays_file")]))
register(BioProcessing("bio_ligand_select", "Select one explicit 1-based ligand record, preserving its chemical state",
                      [ParamSpec("ligands_file", "path", required=True), ParamSpec("record_id", "int", required=True)],
                      ["ligand", "json_report"], [("ligand_library", "ligands_file")], ["rdkit"],
                      example={"ligands_file": "ligands.sdf", "record_id": 1}))
register(BioProcessing("bio_box_from_ligand", "Build an angstrom binding box from one 3D ligand and explicit heavy-atom padding",
                      [ParamSpec("ligand_file", "path", required=True), ParamSpec("padding", "float", default=4.0, unit="angstrom")],
                      ["binding_box", "json_report"], [("ligand", "ligand_file")], ["rdkit"]))
register(BioProcessing("bio_ligand_states", "Enumerate RDKit tautomers with explicit limits; no pH population model",
                      [ParamSpec("ligand_file", "path", required=True), ParamSpec("max_states", "int", default=32)],
                      ["ligand_state_set", "json_report"], [("ligand", "ligand_file")], ["rdkit"]))
register(BioProcessing("bio_ligand_depict", "Export an RDKit SVG grid without changing input coordinates",
                      [ParamSpec("ligands_file", "path", required=True)], ["svg", "json_report"],
                      [("ligand_library", "ligands_file")], ["rdkit"]))
register(BioProcessing("bio_structure_select", "Explicit model/chain/water/hetero/altloc coordinate selection; not receptor preparation",
                      [ParamSpec("structure_file", "path", required=True), ParamSpec("model", "int"), ParamSpec("chains", "list"),
                       ParamSpec("waters", "str", default="keep", enum=("keep", "drop")),
                       ParamSpec("hetero", "str", default="keep", enum=("keep", "drop")),
                       ParamSpec("altloc", "str", default="all", enum=("all", "primary"))],
                      ["macromolecule", "json_report"], [("macromolecule", "structure_file")]))
register(BioProcessing("bio_prediction_import", "Import explicitly selected model files and raw confidence with checksums",
                      [ParamSpec("source_dir", "path", required=True), ParamSpec("models", "list", required=True),
                       ParamSpec("engine", "str", required=True), ParamSpec("engine_version", "str", required=True)],
                      ["prediction_bundle", "json_report"]))
register(BioProcessing("bio_model_select", "Select an explicit model from a checksum-verified prediction bundle",
                      [ParamSpec("bundle_file", "path", required=True), ParamSpec("model_id", "str", required=True)],
                      ["macromolecule", "json_report"], [("prediction_bundle", "bundle_file")]))
register(BioProcessing("bio_visualize", "Export an offline 3Dmol receptor/pose viewer and PyMOL/ChimeraX scripts",
                      [ParamSpec("structure_file", "path", required=True), ParamSpec("poses_file", "path"),
                       ParamSpec("site_file", "path"), ParamSpec("model", "int")],
                      ["html_view", "view_manifest", "json_report"], [("macromolecule", "structure_file")], ["rdkit"]))
