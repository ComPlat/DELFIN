"""OCCUPIER adapter — adaptive spin-state / redox via the existing engine.

A thin building-block wrapper around :func:`delfin.occupier.run_OCCUPIER`,
mirroring how :mod:`delfin.tools.adapters.imag`'s ``imag_fix`` wraps
``delfin.imag.run_IMAG``.  It exposes the OCCUPIER engine as a registered
pipeline step without touching the engine or the CONTROL.txt-driven workflow.

The OCCUPIER engine is **parameterless and CWD-driven**: it reads ``CONTROL.txt``
from the current working directory and writes ``OCCUPIER.txt`` / ``occupier.log``
there (``delfin/occupier.py``).  Since it cannot be pointed elsewhere without
modifying it (out of scope), this adapter runs it in a **subprocess with
``cwd=work_dir``** — that isolates the working directory per process, so it is
thread-safe inside parallel pipeline branches and never calls ``os.chdir`` in
the parent (honouring the StepAdapter contract).

Precondition: ``work_dir`` must contain a valid ``CONTROL.txt`` (typically with
``method=OCCUPIER``) and the input geometry it references — either staged there
beforehand, supplied via the ``control_file`` / ``geometry`` parameters, or left
in place from an upstream step.  OCCUPIER manages its own ORCA resources via the
``PAL`` key in CONTROL.txt; the pipeline ``cores`` argument is not injected.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec
from delfin.tools._types import ErrorKind, StepResult, StepStatus

# Run the parameterless, CWD-driven engine in a child process whose working
# directory is work_dir.  delfin.occupier is imported in the CHILD, so this
# adapter module stays cheap to import during registry discovery.
_CHILD_CODE = "from delfin.occupier import run_OCCUPIER; run_OCCUPIER()"


class OccupierAdapter(StepAdapter):
    name = "occupier"
    description = "Adaptive spin-state identification + redox potentials (OCCUPIER engine)"
    produces_geometry = False
    category = "dft"
    params = (
        ParamSpec("control_file", "path", default="CONTROL.txt",
                  description="CONTROL.txt to run; copied into work_dir if it lives elsewhere"),
        ParamSpec("timeout", "int", unit="s",
                  description="Optional wall-clock timeout for the OCCUPIER run"),
    )
    # OCCUPIER drives a whole sub-workflow from CONTROL.txt rather than consuming
    # a single pipeline capability port, so it declares no consumes/produces tags.
    requires_binaries = ("orca",)

    def validate_params(self, **kwargs: Any) -> None:
        # No hard requirement here: a CONTROL.txt may already sit in work_dir.
        # The execute() step fails cleanly (MISSING_INPUT) if none is found.
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
                error="OCCUPIER needs a CONTROL.txt in work_dir "
                      "(or a 'control_file' path pointing to one)",
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
                error=f"OCCUPIER timed out after {timeout}s",
                error_kind=ErrorKind.TIMEOUT,
            )
        except Exception as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=f"OCCUPIER subprocess error: {exc}",
                error_kind=ErrorKind.INTERNAL,
            )

        # Persist captured output for inspection.
        try:
            (work_dir / "occupier_stdout.log").write_text(proc.stdout or "")
            if proc.stderr:
                (work_dir / "occupier_stderr.log").write_text(proc.stderr)
        except OSError:
            pass

        # --- Collect results ---
        report = work_dir / "OCCUPIER.txt"
        log = work_dir / "occupier.log"
        artifacts = {}
        if report.is_file():
            artifacts["occupier_report"] = report
        if log.is_file():
            artifacts["occupier_log"] = log

        if proc.returncode == 0 and report.is_file():
            return self._make_result(
                self.name, StepStatus.SUCCESS, work_dir, start,
                output_file=report,
                artifacts=artifacts,
            )
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start,
            output_file=log if log.is_file() else None,
            artifacts=artifacts,
            error=(f"OCCUPIER failed (exit {proc.returncode}); "
                   "OCCUPIER.txt was not produced"),
            error_kind=ErrorKind.TOOL_FAILED,
        )


register(OccupierAdapter())
