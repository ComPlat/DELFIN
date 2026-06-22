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

    @property
    def ok(self) -> bool:
        return self.status == StepStatus.SUCCESS

    def require(self) -> "StepResult":
        """Return *self* if successful, raise :class:`StepError` otherwise."""
        if not self.ok:
            raise StepError(self.step_name, self.error or "step failed")
        return self


class StepError(RuntimeError):
    """Raised when a step fails and the caller used ``.require()``."""

    def __init__(self, step_name: str, message: str) -> None:
        self.step_name = step_name
        super().__init__(f"[{step_name}] {message}")
