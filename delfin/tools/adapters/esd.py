"""ESD adapter — excited-state dynamics via the existing engine.

A thin building-block wrapper around :func:`delfin.esd_module.execute_esd_jobs`,
mirroring how ``imag_fix`` wraps ``run_IMAG`` and ``occupier`` wraps
``run_OCCUPIER``.  It exposes the ESD engine (ISC/RISC/IC rates, fluorescence /
phosphorescence lifetimes) as a registered pipeline step without modifying the
engine or the CONTROL.txt-driven workflow.

ESD is driven by the ``ESD_*`` keys in CONTROL.txt (``ESD_modul=yes``, ``states``,
``ISCs``, ``ICs`` …) and writes its results under ``./ESD/`` relative to the
current working directory (``delfin/esd_module.py``: ``Path("ESD").resolve()``).
Since the engine resolves that directory against the CWD, this adapter runs it in
a **subprocess with ``cwd=work_dir``** — isolating the working directory per
process, so it is thread-safe inside parallel pipeline branches and never calls
``os.chdir`` in the parent (honouring the StepAdapter contract).

The engine takes explicit parameters (charge / solvent / metals / basis sets).
Any of these may be overridden via this adapter's params; otherwise they are
derived from CONTROL.txt, exactly like ``workflows/contrib/esd_workflow.py``.

Precondition: ``work_dir`` must contain a valid ``CONTROL.txt`` (with
``ESD_modul=yes`` and the desired states) and the geometry/workspace ESD needs —
staged beforehand, supplied via ``control_file`` / ``geometry``, or produced by
an upstream step.  ESD manages its own ORCA resources via the ``PAL`` key in
CONTROL.txt; the pipeline ``cores`` argument is not injected.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import ErrorKind, StepResult, StepStatus

# Run the engine in a child process whose working directory is work_dir.  The
# child reads CONTROL.txt (and an optional overrides sidecar) and calls
# execute_esd_jobs exactly as the contrib workflow wrapper does.  delfin.esd_module
# is imported in the CHILD, so this adapter module stays cheap to import.
_CHILD_CODE = r"""
import json, sys
from pathlib import Path
from delfin.config import read_control_file
from delfin.esd_module import execute_esd_jobs

config = read_control_file("CONTROL.txt")
ov = {}
_p = Path("_esd_overrides.json")
if _p.is_file():
    try:
        ov = json.loads(_p.read_text())
    except Exception:
        ov = {}

try:
    _charge = int(ov.get("charge", config.get("charge", 0) or 0))
except (TypeError, ValueError):
    _charge = 0

res = execute_esd_jobs(
    config=config,
    charge=_charge,
    solvent=ov.get("solvent", config.get("solvent", "") or ""),
    metals=ov.get("metals", config.get("metals", []) or []),
    main_basisset=ov.get("main_basisset", config.get("main_basisset", "def2-SVP") or "def2-SVP"),
    metal_basisset=ov.get("metal_basisset", config.get("metal_basisset", "") or ""),
)

out = {
    "completed": len(getattr(res, "completed", []) or []),
    "failed": len(getattr(res, "failed", {}) or {}),
    "skipped": len(getattr(res, "skipped", {}) or {}),
    "success": bool(getattr(res, "success", True)),
}
Path("_esd_result.json").write_text(json.dumps(out))
sys.exit(0 if out["success"] else 1)
"""


class EsdAdapter(StepAdapter):
    name = "esd"
    description = "Excited-state dynamics: ISC/RISC/IC rates, fluorescence, phosphorescence (ESD engine)"
    produces_geometry = False
    category = "dft"
    params = (
        ParamSpec("control_file", "path", default="CONTROL.txt",
                  description="CONTROL.txt to run; copied into work_dir if it lives elsewhere"),
        ParamSpec("charge", "int", description="Override molecular charge (default: from CONTROL.txt)"),
        ParamSpec("solvent", "str", description="Override solvent (default: from CONTROL.txt)"),
        ParamSpec("metals", "list", description="Override metal list (default: from CONTROL.txt)"),
        ParamSpec("main_basisset", "str", description="Override main basis set (default: from CONTROL.txt)"),
        ParamSpec("metal_basisset", "str", description="Override metal basis set (default: from CONTROL.txt)"),
        ParamSpec("timeout", "int", unit="s", description="Optional wall-clock timeout for the ESD run"),
    )
    data_keys = (
        DataKeySpec("completed", "int", "", "Number of completed ESD jobs"),
        DataKeySpec("failed", "int", "", "Number of failed ESD jobs"),
        DataKeySpec("skipped", "int", "", "Number of skipped ESD jobs"),
    )
    requires_binaries = ("orca",)

    def validate_params(self, **kwargs: Any) -> None:
        # No hard requirement: a CONTROL.txt may already sit in work_dir.
        # execute() fails cleanly (MISSING_INPUT) if none is found.
        return

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        start = time.monotonic()

        # --- Stage CONTROL.txt into work_dir ---
        dest_control = work_dir / "CONTROL.txt"
        control_arg = kwargs.get("control_file", "CONTROL.txt")
        if control_arg:
            src_control = Path(control_arg)
            if src_control.is_file() and src_control.resolve() != dest_control.resolve():
                shutil.copy2(src_control, dest_control)
        if not dest_control.is_file():
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error="ESD needs a CONTROL.txt in work_dir "
                      "(with ESD_modul=yes), or a 'control_file' path pointing to one",
                error_kind=ErrorKind.MISSING_INPUT,
            )

        # --- Stage the input geometry referenced by CONTROL.txt, if provided ---
        if geometry is not None:
            input_name = "input.xyz"
            try:
                from delfin.config import read_control_file
                cfg = read_control_file(str(dest_control))
                input_name = cfg.get("input_file") or input_name
            except Exception:
                pass
            dest_geom = work_dir / input_name
            if not dest_geom.exists():
                try:
                    shutil.copy2(geometry, dest_geom)
                except OSError as exc:
                    return self._make_result(
                        self.name, StepStatus.FAILED, work_dir, start,
                        error=f"could not stage geometry into work_dir: {exc}",
                        error_kind=ErrorKind.INTERNAL,
                    )

        # --- Write optional parameter overrides for the child ---
        overrides = {
            k: kwargs[k]
            for k in ("charge", "solvent", "metals", "main_basisset", "metal_basisset")
            if k in kwargs
        }
        sidecar = work_dir / "_esd_overrides.json"
        try:
            if overrides:
                sidecar.write_text(json.dumps(overrides))
            elif sidecar.exists():
                sidecar.unlink()
        except OSError:
            pass

        # --- Run the engine in a cwd-isolated subprocess ---
        timeout = kwargs.get("timeout")
        try:
            proc = subprocess.run(
                [sys.executable, "-c", _CHILD_CODE],
                cwd=str(work_dir),
                env=dict(os.environ),
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"ESD timed out after {timeout}s",
                error_kind=ErrorKind.TIMEOUT,
            )
        except Exception as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"ESD subprocess error: {exc}",
                error_kind=ErrorKind.INTERNAL,
            )

        # Persist captured output for inspection.
        try:
            (work_dir / "esd_stdout.log").write_text(proc.stdout or "")
            if proc.stderr:
                (work_dir / "esd_stderr.log").write_text(proc.stderr)
        except OSError:
            pass

        # --- Collect results ---
        data: dict[str, Any] = {}
        result_file = work_dir / "_esd_result.json"
        if result_file.is_file():
            try:
                data = json.loads(result_file.read_text())
            except (OSError, json.JSONDecodeError):
                data = {}

        esd_dir = work_dir / "ESD"
        artifacts = {"esd_dir": esd_dir} if esd_dir.is_dir() else {}

        if proc.returncode == 0:
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                output_file=esd_dir if esd_dir.is_dir() else None,
                data=data,
                artifacts=artifacts,
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            data=data,
            artifacts=artifacts,
            error=(f"ESD failed (exit {proc.returncode}); "
                   f"{data.get('failed', '?')} job(s) failed, "
                   f"{data.get('skipped', '?')} skipped"),
            error_kind=ErrorKind.TOOL_FAILED,
        )


register(EsdAdapter())
