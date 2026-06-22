"""Generic Python function adapter — run ANY Python callable as a pipeline step.

Makes every installed Python module (RDKit, ASE, numpy, pandas, scipy, …)
usable in pipelines without writing a custom adapter.

Usage::

    # Call any module function
    run_step("python_func",
             module="rdkit.Chem",
             function="MolToSmiles",
             args=[mol],
             geometry="mol.xyz")

    # In a pipeline YAML
    steps:
      - step: python_func
        module: rdkit.Chem.Descriptors
        function: MolWt
        args_from_geometry: true    # passes geometry path as first arg

    # Call with keyword arguments
    run_step("python_func",
             module="numpy",
             function="mean",
             kwargs={"a": [1, 2, 3]})

    # Call a method on the geometry (ASE-style)
    run_step("python_func",
             module="ase.io",
             function="read",
             args=["{geometry}"],   # resolved at runtime
             store_as="atoms")

The return value of the function is stored in ``result.data["return_value"]``.
If the function returns a Path-like object ending in ``.xyz``, it is also
set as the output geometry.
"""

from __future__ import annotations

import importlib
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from delfin.tools._registry import register


class PythonFuncAdapter(StepAdapter):
    """Call any Python function as a pipeline step."""

    name = "python_func"
    description = "Call any Python function (RDKit, ASE, numpy, pandas, …)"
    produces_geometry = False  # depends on the function

    def validate_params(self, **kwargs: Any) -> None:
        if "module" not in kwargs:
            raise ValueError("'module' parameter is required (e.g. 'rdkit.Chem')")
        if "function" not in kwargs:
            raise ValueError("'function' parameter is required (e.g. 'MolFromSmiles')")

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        start = time.monotonic()

        module_path = kwargs["module"]
        func_name = kwargs["function"]
        func_args = list(kwargs.get("args", []))
        func_kwargs = dict(kwargs.get("kwargs", {}))
        store_as = kwargs.get("store_as", "return_value")

        # Resolve "{geometry}" placeholder in args
        if geometry:
            func_args = [str(geometry) if a == "{geometry}" else a for a in func_args]
            func_kwargs = {
                k: str(geometry) if v == "{geometry}" else v
                for k, v in func_kwargs.items()
            }

        # If args_from_geometry is set, pass geometry as first argument
        if kwargs.get("args_from_geometry") and geometry:
            func_args.insert(0, str(geometry))

        try:
            mod = importlib.import_module(module_path)
            func = getattr(mod, func_name)
        except (ImportError, AttributeError) as e:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"Cannot import {module_path}.{func_name}: {e}",
            )

        try:
            result = func(*func_args, **func_kwargs)
        except Exception as e:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"{module_path}.{func_name}() raised: {e}",
            )

        # Store result
        data: dict[str, Any] = {store_as: result}

        # Check if result is a geometry path
        geom_out = None
        if isinstance(result, (str, Path)):
            p = Path(result)
            if p.suffix == ".xyz" and p.is_file():
                geom_out = p

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=geom_out,
            data=data,
        )


class PythonScriptAdapter(StepAdapter):
    """Execute a Python script file as a pipeline step.

    The script is run in a subprocess with the work_dir as cwd.
    Geometry path and extra kwargs are passed as environment variables.
    """

    name = "python_script"
    description = "Run a Python script as a pipeline step"
    produces_geometry = False

    def validate_params(self, **kwargs: Any) -> None:
        if "script" not in kwargs:
            raise ValueError("'script' parameter is required (path to .py file)")

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        import subprocess
        import sys
        import json

        start = time.monotonic()
        script = Path(kwargs["script"]).resolve()

        if not script.is_file():
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"Script not found: {script}",
            )

        # Build environment
        import os
        env = os.environ.copy()
        env["DELFIN_WORK_DIR"] = str(work_dir)
        env["DELFIN_CORES"] = str(cores)
        if geometry:
            env["DELFIN_GEOMETRY"] = str(geometry)

        # Pass extra kwargs as JSON
        extra = {k: v for k, v in kwargs.items()
                 if k not in ("script", "_prev_artifacts")}
        if extra:
            env["DELFIN_PARAMS"] = json.dumps(extra)

        proc = subprocess.run(
            [sys.executable, str(script)],
            cwd=work_dir,
            env=env,
            capture_output=True,
            text=True,
            timeout=kwargs.get("timeout"),
        )

        # Check for output geometry
        geom_out = None
        for candidate in (work_dir / "output.xyz", work_dir / "result.xyz"):
            if candidate.is_file():
                geom_out = candidate
                break

        data: dict[str, Any] = {}
        # Try to read JSON result file if script wrote one
        result_json = work_dir / "result.json"
        if result_json.is_file():
            try:
                data = json.loads(result_json.read_text())
            except Exception:
                pass

        if proc.returncode != 0:
            error = proc.stderr.strip() or f"Script exited with code {proc.returncode}"
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                geometry=geom_out,
                data=data,
                error=error,
            )

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=geom_out,
            data=data,
        )


register(PythonFuncAdapter())
register(PythonScriptAdapter())
