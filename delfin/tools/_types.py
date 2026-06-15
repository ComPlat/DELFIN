"""Core types for the DELFIN tools adapter layer."""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional


class StepStatus(enum.Enum):
    SUCCESS = "success"
    FAILED = "failed"
    SKIPPED = "skipped"


class ErrorKind(enum.Enum):
    """Optional, machine-readable failure category.

    The free-text :attr:`StepResult.error` stays the primary message; this lets
    pipeline policies (retry, branch, reactive) react to *why* a step failed
    without substring-matching.  Adapters opt in by setting it on failure; until
    they do it stays :attr:`NONE` and nothing changes.
    """

    NONE = "none"
    MISSING_PARAM = "missing_param"
    MISSING_INPUT = "missing_input"
    BINARY_NOT_FOUND = "binary_not_found"
    CONVERGENCE = "convergence"          # SCF / geometry did not converge
    TIMEOUT = "timeout"
    PARSE = "parse"                      # could not parse tool output
    TOOL_FAILED = "tool_failed"          # nonzero exit / generic failure
    INTERNAL = "internal"                # adapter raised an exception


# Hartree ↔ electron-volt (CODATA).
_EH_TO_EV = 27.211386245988


@dataclass
class StepResult:
    """Uniform result returned by every ``run_step()`` call."""

    step_name: str
    status: StepStatus

    # Absolute path to output geometry (XYZ).  None for analysis-only steps.
    geometry: Optional[Path] = None

    # Primary output / log file for inspection.
    output_file: Optional[Path] = None

    # Isolated working directory where the step ran.
    work_dir: Optional[Path] = None

    # Structured data extracted by the adapter (energies, descriptors, …).
    data: Dict[str, Any] = field(default_factory=dict)

    # Auxiliary output files (GBW, hessian, conformer ensemble, …).
    artifacts: Dict[str, Path] = field(default_factory=dict)

    # Human-readable error message on failure.
    error: Optional[str] = None

    # Wall-clock runtime in seconds.
    elapsed_seconds: float = 0.0

    # Optional machine-readable failure category (appended; default keeps the
    # historical StepResult(step_name=, status=) construction identical).
    error_kind: "ErrorKind" = ErrorKind.NONE

    @property
    def ok(self) -> bool:
        return self.status == StepStatus.SUCCESS

    def require(self) -> "StepResult":
        """Return *self* if successful, raise :class:`StepError` otherwise."""
        if not self.ok:
            raise StepError(self.step_name, self.error or "step failed")
        return self

    # --- Canonical energy accessors -------------------------------------
    # The free `data` dict is left untouched; these normalise the common
    # energy keys/units so cross-adapter aggregation (redox, screening) does
    # not have to guess between energy_Eh and energy_eV.

    def energy_eh(self) -> Optional[float]:
        """Best available energy in Hartree, or ``None`` if not present."""
        d = self.data
        if "energy_Eh" in d:
            return d["energy_Eh"]
        if "gibbs_Eh" in d:
            return d["gibbs_Eh"]
        if "energy_eV" in d and d["energy_eV"] is not None:
            return d["energy_eV"] / _EH_TO_EV
        return None

    def energy_ev(self) -> Optional[float]:
        """Best available energy in electron-volts, or ``None`` if not present."""
        if "energy_eV" in self.data and self.data["energy_eV"] is not None:
            return self.data["energy_eV"]
        eh = self.energy_eh()
        return None if eh is None else eh * _EH_TO_EV


class StepError(RuntimeError):
    """Raised when a step fails and the caller used ``.require()``."""

    def __init__(self, step_name: str, message: str) -> None:
        self.step_name = step_name
        super().__init__(f"[{step_name}] {message}")
