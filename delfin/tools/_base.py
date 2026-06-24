"""Abstract base class for tool adapters."""

from __future__ import annotations

import abc
import shutil
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional, Tuple

from delfin.tools._types import StepResult, StepStatus

if TYPE_CHECKING:  # pragma: no cover
    from delfin.tools._spec import DataKeySpec, ParamSpec, StepContract


class StepAdapter(abc.ABC):
    """Base class every tool adapter inherits from.

    Subclasses must set :attr:`name` and :attr:`description` and implement
    :meth:`execute`.  They must **never** call ``os.chdir()``; all file I/O
    goes through the *work_dir* parameter.

    Subclasses *may* additionally declare a machine-readable **contract** —
    :attr:`category`, :attr:`params`, :attr:`produces`, :attr:`consumes`,
    :attr:`data_keys`, :attr:`requires_binaries`, :attr:`requires_python`.
    These are all optional and default empty, so existing adapters that
    declare nothing keep working unchanged.  See :mod:`delfin.tools._spec`
    and :meth:`contract`.
    """

    name: str = ""
    description: str = ""
    produces_geometry: bool = True

    # --- Optional declarative contract (all default empty / opt-in) ---------
    # See delfin.tools._spec for the dataclasses and the capability-tag
    # vocabulary.  Consumers assemble these into a StepContract via contract().
    category: str = ""
    params: Tuple["ParamSpec", ...] = ()
    produces: Tuple[str, ...] = ()          # capability tags emitted
    consumes: Tuple[str, ...] = ()          # capability tags required upstream
    data_keys: Tuple["DataKeySpec", ...] = ()
    requires_binaries: Tuple[str, ...] = ()
    requires_python: Tuple[str, ...] = ()
    # opt-in input auto-wiring: {capability: kwarg} or ((capability, kwarg), ...)
    wires: Any = ()

    def contract(self) -> "StepContract":
        """Assemble this adapter's declarative :class:`StepContract`.

        The primitive :attr:`produces_geometry` flag is folded into the
        ``produces`` capability set as the ``geometry`` tag, so callers have a
        single uniform view of every port the adapter exposes.
        """
        from delfin.tools._spec import StepContract

        produced = set(self.produces)
        if self.produces_geometry:
            produced.add("geometry")
        wires = self.wires.items() if isinstance(self.wires, dict) else self.wires
        return StepContract(
            name=self.name,
            description=self.description,
            category=self.category,
            produces_geometry=self.produces_geometry,
            params=tuple(self.params),
            produces=frozenset(produced),
            consumes=frozenset(self.consumes),
            data_keys=tuple(self.data_keys),
            requires_binaries=frozenset(self.requires_binaries),
            requires_python=frozenset(self.requires_python),
            wires=tuple((str(c), str(k)) for c, k in wires),
        )

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
