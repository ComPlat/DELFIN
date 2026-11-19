"""Version-checked optional docking runners using DELFIN's process runtime."""
from __future__ import annotations

import math
from pathlib import Path
import shutil
import subprocess
import time

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import ErrorKind, StepStatus
from delfin.tools.bio import BindingBox
from delfin.tools.bio.artifacts import file_hash, read_json, validate_artifact, write_json


class DockingAdapter(StepAdapter):
    category = "structure"
    produces_geometry = False
    data_keys = (DataKeySpec("engine_version", "str"), DataKeySpec("scores", "list"))

    def __init__(self, engine):
        self.engine = engine
        self.name = engine + "_dock"
        self.description = f"Run optional {engine} docking with an exact version banner, prepared inputs, explicit box/seed and timeout"
        self.requires_binaries = (engine,)
        self.requires_python = () if engine == "vina" else ("rdkit",)
        self.produces = ("poses_pdbqt" if engine == "vina" else "pose_set", "json_report")
        self.required_artifacts = self.produces
        self.consumes = ("receptor_pdbqt" if engine == "vina" else "prepared_receptor",
                         "ligand_pdbqt" if engine == "vina" else "ligand", "binding_box")
        self.wires = tuple(zip(self.consumes, ("receptor_file", "ligand_file", "site_file")))
        self.params = (
            ParamSpec("receptor_file", "path", required=True), ParamSpec("ligand_file", "path", required=True),
            ParamSpec("site_file", "path", required=True),
            ParamSpec("expected_version", "str", required=True, description="Exact stripped stdout/stderr banner from --version"),
            ParamSpec("seed", "int", required=True), ParamSpec("exhaustiveness", "int", default=8),
            ParamSpec("num_modes", "int", default=9), ParamSpec("timeout", "float", default=3600.0, unit="s"),
        )
        self.example_kwargs = {"receptor_file": "receptor.pdbqt" if engine == "vina" else "receptor.pdb",
                               "ligand_file": "ligand.pdbqt" if engine == "vina" else "ligand.sdf", "site_file": "box.json",
                               "expected_version": "AutoDock Vina v1.2.5" if engine == "vina" else "replace with installed version banner",
                               "seed": 7}

    def validate_params(self, **kwargs):
        for key in ("receptor_file", "ligand_file", "site_file", "expected_version", "seed"):
            if kwargs.get(key) is None or kwargs.get(key) == "":
                raise ValueError(f"{key} is required")
        for key, default, minimum in (("seed", None, 0), ("exhaustiveness", 8, 1), ("num_modes", 9, 1)):
            value = kwargs.get(key, default)
            if type(value) is not int or not minimum <= value < 2**31:
                raise ValueError(f"{key} must be an integer in {minimum}..2**31-1")
        timeout = kwargs.get("timeout", 3600.0)
        if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or timeout <= 0:
            raise ValueError("timeout must be a positive finite number")

    def runtime_preflight(self, **kwargs):
        from delfin import qm_runtime
        resolved = qm_runtime.resolve_tool(self.engine)
        if resolved is None:
            raise ImportError(f"optional engine is not installed: {self.engine}")
        version = qm_runtime.run_tool(resolved.path, ["--version"], capture_output=True,
                                      timeout=30, allow_install=False, track_process=True)
        banner = ((version.stdout or "") + (version.stderr or "")).strip()
        if version.returncode or banner != kwargs["expected_version"].strip():
            raise ValueError(f"version mismatch: expected {kwargs['expected_version']!r}, found {banner!r}")
        return {"engine_version": banner, "executable": resolved.path,
                "executable_sha256": file_hash(Path(resolved.path))}

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        from delfin import qm_runtime
        start = time.monotonic()
        work_dir = Path(work_dir)
        try:
            self.validate_params(**kwargs)
            if type(cores) is not int or cores < 1:
                raise ValueError("cores must be positive")
            work_dir.mkdir(parents=True, exist_ok=True)
            if any(work_dir.iterdir()):
                raise ValueError("docking requires a fresh empty step directory")
            receptor, ligand = Path(kwargs["receptor_file"]).resolve(), Path(kwargs["ligand_file"]).resolve()
            validate_artifact(receptor, self.consumes[0])
            validate_artifact(ligand, self.consumes[1])
            identity = None
            if self.engine != "vina":
                from delfin.tools.bio.processing import load_ligands, chemical_identity
                reference = load_ligands(ligand)
                if len(reference) != 1 or not reference[0].GetNumConformers():
                    raise ValueError("one ligand record with coordinates is required per docking job")
                identity = chemical_identity(reference[0])["chemical_state_id"]
            box = BindingBox.from_dict(read_json(kwargs["site_file"]))
            runtime = self.runtime_preflight(**kwargs)
            banner = runtime["engine_version"]
            rec_name = "receptor" + receptor.suffix
            lig_name = "ligand" + ligand.suffix
            shutil.copyfile(receptor, work_dir / rec_name)
            shutil.copyfile(ligand, work_dir / lig_name)
            output = work_dir / ("poses.pdbqt" if self.engine == "vina" else "poses.sdf")
            argv = ["--receptor", rec_name, "--ligand", lig_name, "--out", output.name,
                    "--seed", str(kwargs["seed"]), "--cpu", str(cores),
                    "--exhaustiveness", str(kwargs.get("exhaustiveness", 8)), "--num_modes", str(kwargs.get("num_modes", 9))]
            for index, axis in enumerate("xyz"):
                argv.extend(["--center_" + axis, str(box.center[index]), "--size_" + axis, str(box.size[index])])
            log = work_dir / "engine.log"
            with log.open("w") as stream:
                result = qm_runtime.run_tool(runtime["executable"], argv, cwd=work_dir, stdout=stream, stderr=subprocess.STDOUT,
                                             timeout=kwargs.get("timeout", 3600.0), allow_install=False, track_process=True)
            if result.returncode:
                return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                         error=f"{self.engine} exited {result.returncode}; see engine.log", error_kind=ErrorKind.TOOL_FAILED,
                                         artifacts={"log": log})
            validate_artifact(output, self.produces[0])
            scores = []
            if self.engine == "vina":
                text = output.read_text()
                for line in text.splitlines():
                    if line.startswith("REMARK VINA RESULT:"):
                        score = float(line.split()[3])
                        if not math.isfinite(score):
                            raise ValueError("non-finite Vina score")
                        scores.append({"pose": len(scores) + 1, "name": "vina_score", "value": score, "unit": "kcal/mol", "kind": "docking_score"})
                starts = sum(line.startswith("MODEL ") for line in text.splitlines())
                ends = sum(line.startswith("ENDMDL") for line in text.splitlines())
                if not scores or starts != ends or starts != len(scores):
                    raise ValueError("incomplete Vina model/score records")
            else:
                from delfin.tools.bio.processing import load_ligands, chemical_identity
                for mol in load_ligands(output):
                    if chemical_identity(mol)["chemical_state_id"] != identity:
                        raise ValueError("docked molecule changed chemical state")
                    scores.append({"pose": len(scores) + 1, "raw_engine_properties": mol.GetPropsAsDict(includePrivate=False, includeComputed=False),
                                   "note": "Engine metrics retain their own names; no conversion to measured affinity"})
            data = {"engine": self.engine, **runtime, "argv": argv,
                    "inputs": {"receptor": file_hash(receptor), "ligand": file_hash(ligand)},
                    "binding_box": box.to_dict(), "seed": kwargs["seed"], "cores": cores,
                    "scores": scores, "scientific_validation": "not_performed"}
            report = work_dir / "docking_report.json"
            write_json(report, data)
            return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start, data=data,
                                     artifacts={self.produces[0]: output, "json_report": report, "log": log})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except FileNotFoundError as exc:
            kind, error = ErrorKind.MISSING_INPUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except (ValueError, KeyError, TypeError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error=error, error_kind=kind)


for _engine in ("vina", "gnina", "smina"):
    register(DockingAdapter(_engine))
