"""Optional Meeko and PLIP integrations with isolated interpreter selection."""
from __future__ import annotations

import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec, DataKeySpec
from delfin.tools._types import ErrorKind, StepStatus
from delfin.tools.bio.artifacts import file_hash, validate_artifact, write_json
from delfin.tools.bio.optional_runtime import run_console


class OptionalProcessing(StepAdapter):
    category = "structure"
    produces_geometry = False
    data_keys = (DataKeySpec("runtime", "dict"),)

    def __init__(self, name, distribution, entry, input_kind, output_kind):
        self.name, self.distribution, self.entry = name, distribution, entry
        self.description = f"Run installed {entry} in an explicit Python environment with exact distribution version and output validation"
        self.consumes = (input_kind,)
        self.produces = (output_kind, "json_report")
        self.required_artifacts = self.produces
        self.wires = ((input_kind, "input_file"),)
        self.params = (ParamSpec("input_file", "path", required=True),
                       ParamSpec("expected_version", "str", required=True),
                       ParamSpec("python", "path", default=sys.executable),
                       ParamSpec("timeout", "float", default=600.0))
        self.example_kwargs = {"input_file": "input.pdb", "expected_version": "0.7.1"}

    def validate_params(self, **kwargs):
        if not kwargs.get("input_file") or not kwargs.get("expected_version"):
            raise ValueError("input_file and expected_version are required")
        value = kwargs.get("timeout", 600.0)
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
            raise ValueError("timeout must be positive and finite")

    def runtime_preflight(self, **kwargs):
        from delfin.tools.bio.optional_runtime import probe_distribution
        return probe_distribution(self.distribution, python=kwargs.get("python", sys.executable),
                                  expected_version=kwargs["expected_version"], work_dir=Path.cwd())

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start = time.monotonic()
        work_dir = Path(work_dir)
        try:
            self.validate_params(**kwargs)
            work_dir.mkdir(parents=True, exist_ok=True)
            if any(work_dir.iterdir()):
                raise ValueError("optional processing requires an empty step directory")
            source = Path(kwargs["input_file"]).resolve()
            validate_artifact(source, self.consumes[0])
            if self.name == "meeko_ligand_prepare":
                from delfin.tools.bio.processing import load_ligands
                molecules = load_ligands(source)
                if len(molecules) != 1 or not molecules[0].GetNumConformers():
                    raise ValueError("Meeko requires one ligand with coordinates")
                if any(atom.GetTotalNumHs(includeNeighbors=False) for atom in molecules[0].GetAtoms()):
                    raise ValueError("prepare explicit ligand hydrogens before Meeko")
                from rdkit import Chem
                with Chem.SDWriter(str(work_dir / "input.sdf")) as writer:
                    writer.write(molecules[0])
                output = work_dir / "ligand.pdbqt"
                argv = ["-i", "input.sdf", "-o", output.name]
            elif self.name == "meeko_receptor_prepare":
                from delfin.tools.bio import read_structure
                structure = read_structure(source)
                if structure.source_format != "pdb" or len(structure.models) != 1:
                    raise ValueError("Meeko receptor integration currently requires one explicitly prepared PDB model")
                standard = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
                if any(r.name not in standard or r.is_modified or r.altlocs for r in structure.first_model.residues):
                    raise ValueError("This Meeko receptor subset requires canonical residues and resolved altlocs; cofactors/metals need an explicit template integration")
                shutil.copyfile(source, work_dir / "input.pdb")
                output = work_dir / "receptor.pdbqt"
                argv = ["--read_pdb", "input.pdb", "-o", "receptor", "-p"]
            elif self.name == "meeko_pose_export":
                shutil.copyfile(source, work_dir / "input.pdbqt")
                output = work_dir / "poses.sdf"
                argv = ["input.pdbqt", "-s", output.name]
            elif self.name == "plip_analyze":
                from delfin.tools.bio import read_structure
                if read_structure(source).source_format != "pdb":
                    raise ValueError("PLIP integration requires an assembled PDB complex")
                shutil.copyfile(source, work_dir / "input.pdb")
                output = work_dir / "plip" / "report.xml"
                argv = ["-f", "input.pdb", "-o", "plip", "-x"]
            else:
                raise ValueError("unknown optional operation")
            runtime = run_console(self.distribution, self.entry, argv, python=kwargs.get("python", sys.executable),
                                  expected_version=kwargs["expected_version"], work_dir=work_dir,
                                  timeout=kwargs.get("timeout", 600.0))
            validate_artifact(output, self.produces[0])
            data = {"runtime": runtime, "source_sha256": file_hash(source), "output_sha256": file_hash(output)}
            if self.name == "plip_analyze":
                tree = ET.parse(output)
                data["binding_sites"] = len(tree.findall(".//bindingsite"))
            report = work_dir / "processing_report.json"
            write_json(report, data)
            return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=data,
                                     artifacts={self.produces[0]: output, "json_report": report, "log": work_dir / "engine.log"})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except FileNotFoundError as exc:
            kind, error = ErrorKind.MISSING_INPUT, str(exc)
        except (ValueError, KeyError, TypeError, ET.ParseError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        except RuntimeError as exc:
            kind, error = ErrorKind.TOOL_FAILED, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error=error, error_kind=kind)


register(OptionalProcessing("meeko_ligand_prepare", "meeko", "mk_prepare_ligand.py", "ligand", "ligand_pdbqt"))
register(OptionalProcessing("meeko_receptor_prepare", "meeko", "mk_prepare_receptor.py", "prepared_receptor", "receptor_pdbqt"))
register(OptionalProcessing("meeko_pose_export", "meeko", "mk_export.py", "poses_pdbqt", "pose_set"))
register(OptionalProcessing("plip_analyze", "plip", "plip", "complex_structure", "interaction_report"))
