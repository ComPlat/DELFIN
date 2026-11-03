"""OPI-backed ORCA geometry / frequency steps that chain into a workflow.

Registers ``orca_opt_opi``, ``orca_freq_opi`` and ``orca_opt_freq_opi`` and
make the single point ``orca_sp_opi`` compose: their ``StepResult`` carries the
optimized geometry (OPT), the ``.gbw`` artifact (offered to the next ``orca_*``
step's ``moread`` start orbitals via the engine's wiring rule) and the
frequency step's thermal corrections (fed to an energy profile).  A step ORCA
did not finish returns ``FAILED`` with its reason, which halts the chain.

All input building / running / reading lives in ``delfin.tools.opi_chain``,
which reuses ``orca_sp_opi``'s validated helpers by import — nothing in that
module is edited.
"""

from __future__ import annotations

import time

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec
from delfin.tools._types import ErrorKind, StepStatus
from delfin.tools.opi_chain import (
    read_geometry_xyz,
    run_opi_job,
    step_result_from,
)

_GEOMETRY_PARAMS = (
    ParamSpec("functional", "str", default="r2SCAN-3c", description="DFT functional"),
    ParamSpec("basis", "str", default="def2-SVP", description="Basis set"),
    ParamSpec("ri", "str", default="", description="RI approximation keyword"),
    ParamSpec("aux_basis", "str", default="", description="Auxiliary basis"),
    ParamSpec("dispersion", "str", default="D4", description="Dispersion correction"),
)


def _geometry_kwargs(geometry, kwargs) -> tuple[list[str], list[list[float]]]:
    """The elements/coordinates to build the structure: the previous step's
    geometry wins; otherwise the caller's explicit ``elements``/``coordinates``."""
    if geometry is not None:
        elements, coordinates = read_geometry_xyz(geometry)
        return elements, coordinates
    try:
        elements = kwargs["elements"]
        coordinates = kwargs["coordinates"]
    except KeyError as exc:
        raise ValueError(f"missing required {exc} (or a chained geometry)") from exc
    return list(elements), list(coordinates)


class OrcaOptOpiAdapter(StepAdapter):
    """ORCA geometry optimization driven by OPI; the optimized geometry and
    the .gbw flow to the next step."""

    name = "orca_opt_opi"
    description = "OPI-backed ORCA geometry optimization"
    category = "dft"
    produces_geometry = True
    consumes = ("geometry",)
    produces = ("gbw", "qc_output", "property_json")
    requires_python = ("opi",)
    params = _GEOMETRY_PARAMS + (
        ParamSpec("maxcore", "int", default=1000, unit="MB",
                  description="Memory per core"),
    )

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        kwargs.pop("_prev_artifacts", None)
        start = time.monotonic()
        try:
            elements, coordinates = _geometry_kwargs(geometry, kwargs)
        except ValueError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.MISSING_PARAM,
            )
        res = run_opi_job(
            job_kind="OPT",
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            functional=kwargs.get("functional", "r2SCAN-3c"),
            basis=kwargs.get("basis", "def2-SVP"),
            ri=kwargs.get("ri", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            dispersion=kwargs.get("dispersion", "D4"),
            extra_keywords=kwargs.get("extra_keywords"),
            scf_options=kwargs.get("scf_options"),
            moread=kwargs.get("moread"),
            maxcore=kwargs.get("maxcore", 1000),
            ncores=cores,
            base_name=kwargs.get("base_name", "calc"),
            work_dir=work_dir,
            timeout=kwargs.get("timeout"),
        )
        return step_result_from(self, res, work_dir, start)


class OrcaFreqOpiAdapter(StepAdapter):
    """OPI-backed ORCA frequency + thermochemistry; the thermal corrections
    (free energy / enthalpy / ZPE) are exposed as data for an energy profile,
    and the .gbw / .hess flow as artifacts."""

    name = "orca_freq_opi"
    description = "OPI-backed ORCA frequency / thermochemistry"
    category = "dft"
    produces_geometry = False
    consumes = ("geometry",)
    produces = ("gbw", "qc_output", "hessian", "property_json")
    requires_python = ("opi",)
    params = _GEOMETRY_PARAMS + (
        ParamSpec("freq_temp", "float", default=298.15, unit="K",
                  description="Thermochemistry temperature"),
        ParamSpec("freq_pressure", "float", default=1.0, unit="atm",
                  description="Thermochemistry pressure"),
    )

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        kwargs.pop("_prev_artifacts", None)
        start = time.monotonic()
        try:
            elements, coordinates = _geometry_kwargs(geometry, kwargs)
        except ValueError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.MISSING_PARAM,
            )
        freq_options = {}
        if kwargs.get("freq_temp") is not None:
            freq_options["temp"] = float(kwargs["freq_temp"])
        if kwargs.get("freq_pressure") is not None:
            freq_options["pressure"] = float(kwargs["freq_pressure"])
        res = run_opi_job(
            job_kind="FREQ",
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            functional=kwargs.get("functional", "r2SCAN-3c"),
            basis=kwargs.get("basis", "def2-SVP"),
            ri=kwargs.get("ri", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            dispersion=kwargs.get("dispersion", "D4"),
            extra_keywords=kwargs.get("extra_keywords"),
            scf_options=kwargs.get("scf_options"),
            freq_options=freq_options or None,
            moread=kwargs.get("moread"),
            maxcore=kwargs.get("maxcore", 1000),
            ncores=cores,
            base_name=kwargs.get("base_name", "calc"),
            work_dir=work_dir,
            timeout=kwargs.get("timeout"),
        )
        return step_result_from(self, res, work_dir, start)


class OrcaOptFreqOpiAdapter(StepAdapter):
    """OPI-backed ORCA combined optimization + frequency; produces the
    optimized geometry, the .hessian and the thermochemistry."""

    name = "orca_opt_freq_opi"
    description = "OPI-backed ORCA optimization + frequency"
    category = "dft"
    produces_geometry = True
    consumes = ("geometry",)
    produces = ("gbw", "qc_output", "hessian", "property_json")
    requires_python = ("opi",)
    params = _GEOMETRY_PARAMS + (
        ParamSpec("freq_temp", "float", default=298.15, unit="K",
                  description="Thermochemistry temperature"),
        ParamSpec("freq_pressure", "float", default=1.0, unit="atm",
                  description="Thermochemistry pressure"),
    )

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        kwargs.pop("_prev_artifacts", None)
        start = time.monotonic()
        try:
            elements, coordinates = _geometry_kwargs(geometry, kwargs)
        except ValueError as exc:
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error=str(exc), error_kind=ErrorKind.MISSING_PARAM,
            )
        freq_options = {}
        if kwargs.get("freq_temp") is not None:
            freq_options["temp"] = float(kwargs["freq_temp"])
        if kwargs.get("freq_pressure") is not None:
            freq_options["pressure"] = float(kwargs["freq_pressure"])
        res = run_opi_job(
            job_kind="OPT_FREQ",
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            functional=kwargs.get("functional", "r2SCAN-3c"),
            basis=kwargs.get("basis", "def2-SVP"),
            ri=kwargs.get("ri", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            dispersion=kwargs.get("dispersion", "D4"),
            extra_keywords=kwargs.get("extra_keywords"),
            scf_options=kwargs.get("scf_options"),
            freq_options=freq_options or None,
            moread=kwargs.get("moread"),
            maxcore=kwargs.get("maxcore", 1000),
            ncores=cores,
            base_name=kwargs.get("base_name", "calc"),
            work_dir=work_dir,
            timeout=kwargs.get("timeout"),
        )
        return step_result_from(self, res, work_dir, start)


register(OrcaOptOpiAdapter())
register(OrcaFreqOpiAdapter())
register(OrcaOptFreqOpiAdapter())
