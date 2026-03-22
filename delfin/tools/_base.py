"""Abstract base class for tool adapters."""

from __future__ import annotations

import abc
import shutil
import time
from pathlib import Path
from typing import Any, Optional

from delfin.tools._types import StepResult, StepStatus


class StepAdapter(abc.ABC):
    """Base class every tool adapter inherits from.

    Subclasses must set :attr:`name` and :attr:`description` and implement
    :meth:`execute`.  They must **never** call ``os.chdir()``; all file I/O
    goes through the *work_dir* parameter.
    """

    name: str = ""
    description: str = ""
    produces_geometry: bool = True

    @abc.abstractmethod
    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        """Run the tool inside *work_dir* and return a :class:`StepResult`."""
        ...

    def validate_params(self, **kwargs: Any) -> None:
        """Optional pre-flight check.  Raise ``ValueError`` for bad params."""

    # ------------------------------------------------------------------
    # Helpers available to all adapters
    # ------------------------------------------------------------------

    @staticmethod
    def _copy_geometry_to_workdir(
        geometry: Path,
        work_dir: Path,
        dest_name: str = "input.xyz",
    ) -> Path:
        dest = work_dir / dest_name
        shutil.copy2(geometry, dest)
        return dest

    @staticmethod
    def _make_result(
        name: str,
        status: StepStatus,
        work_dir: Path,
        start_time: float,
        **kwargs: Any,
    ) -> StepResult:
        return StepResult(
            step_name=name,
            status=status,
            work_dir=work_dir,
            elapsed_seconds=time.monotonic() - start_time,
            **kwargs,
        )
