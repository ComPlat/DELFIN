"""OPI building block: geometry optimisation with an MLIP through ORCA's ExtOpt.

The block (``orca_mlip_opt_opi``) drives an ORCA ``! ExtOpt`` optimisation
whose energy / gradients come from an MLIP through the orca-mlips command,
reusing the shared input builder
:func:`delfin.mlp_tools.orca_mlips.build_orca_mlips_input` (the same route the
GOAT MLIP block uses).  The ``%method`` ``ProgExt`` / ``Ext_Params`` block is
written with OPI's typed objects inside that shared builder; this block adds
the chemistry guard-rails and the result reader.

Manual ground truth (ORCA 6.1):
* ``optimizations.txt`` 4.1.10 "Optimizing with External Methods" -- ``!ExtOpt``
  "must be combined with other keywords like !Opt or !GOAT to select the
  specific run type"; the program that supplies energies / gradients is named
  by ``ProgExt`` and its arguments by ``Ext_Params`` in a ``%method`` block.
  ExtOpt discards all electronic-structure keywords that are not consumed by
  the external program, so the chemistry is whatever the MLIP (and its
  element / charge coverage) actually computes -- never silently changed here.
* The MLIP reports the "current total energy in Eh", so the ``FINAL SINGLE
  POINT ENERGY`` line in the ``.out`` is read in Hartree (unit ``Eh``).

Refusals (each raises :class:`ValueError` with the reason; never "fixed"):
* ``run_type`` other than ``opt`` / ``opt_freq`` (this block is Opt-scoped;
  GOAT, NEB, IRC and scan are separate blocks).
* A backend that is not an orca-mlips plugin (``uma``, ``orb``, ``mace``,
  ``aimnet2``).
* An MLIP that is not installed -- the refusal says how to install it.
* Elements the chosen model does not cover (only when the registry model entry
  enumerates its element set; a universal / model-dependent model is not
  guessed at).
* A GPU ``device`` -- GPU execution goes only through the cluster submission
  path that already exists, never locally on this host.
* A multiplicity below the physical minimum (``mult < 1``).

The reader returns typed values with units through :func:`delfin.orca_read.read`
(final energy in Eh), the optimisation geometry (the ``*xyzfile``-named file),
and the geometry-convergence and normal-termination flags.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional, Sequence

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import ErrorKind, StepResult, StepStatus

log = get_logger(__name__)

# The two job kinds this Opt block drives: a plain ``! Opt ExtOpt`` and a
# combined ``! Opt Freq ExtOpt``.  Other ExtOpt run types (GOAT, NEB, IRC,
# scan) belong to their own blocks.
_OPT_RUN_TYPES = ("opt", "opt_freq")

#: device values that would run the MLIP on a GPU.
_GPU_DEVICES = ("gpu", "cuda", "cuda0", "cuda1", "auto")

_MANUAL_EXT_OPT = "ORCA 6.1 optimizations.txt 4.1.10"


def _normalize_element(elem: str) -> str:
    """A registry element symbol "cl" -> "Cl"; registry sets use proper case."""
    el = elem.strip()
    if len(el) == 1:
        return el.upper()
    return el[0].upper() + el[1:].lower()


def _registry_backend_spec(backend: str):
    """The :mod:`~delfin.mlp_tools.registry` model spec behind an orca-mlips
    plugin, or ``None`` when the backend is not an orca-mlips plugin."""
    import delfin.mlp_tools.orca_mlips as om
    from delfin.mlp_tools import registry

    name = backend.strip().lower()
    reg_name = om.MLIPS_REGISTRY_BACKEND.get(name)
    if reg_name is None:
        return None
    return registry.get_backend(reg_name)


def validate_mlip_opt(
    *,
    run_type: str = "opt",
    backend: str = "uma",
    device: str = "cpu",
    elements: Sequence[str],
    charge: int = 0,
    mult: int = 1,
    check_installed: bool = True,
) -> list[str]:
    """The Opt/ExtOpt-MLIP guard-rails; returns a list of error messages.

    Empty list = the inputs are chemistry-safe to hand to the MLIP via
    ``! Opt ExtOpt``.  Raises nothing; the caller decides.

    ``check_installed`` gates whether the orca-mlips program (ExtOpt's
    ``ProgExt``) must be present.  Building an input never needs the program,
    so callers that only write the ``.inp`` pass ``check_installed=False``;
    callers that would actually launch ORCA (the step adapter's ``execute``)
    keep the default and get the refusal with an install hint.
    """
    import delfin.mlp_tools.orca_mlips as om

    errors: list[str] = []
    rt = str(run_type).strip().lower()
    if rt not in _OPT_RUN_TYPES:
        errors.append(
            f"run_type {run_type!r} not supported by the Opt block: use one of "
            f"{', '.join(_OPT_RUN_TYPES)}"
        )
    name = str(backend).strip().lower()
    if name not in om._MLIPS_BACKENDS:
        errors.append(
            f"backend {backend!r} is not an orca-mlips plugin: use one of "
            f"{', '.join(om._MLIPS_BACKENDS)}"
        )
    else:
        if check_installed:
            install = om.detect_orca_mlips()
            if not install.get("available"):
                errors.append(
                    "the orca-mlips command line program is not installed; "
                    f"install it with e.g. {install.get('install_hint') or 'pip install orca-mlips'}"
                )
        spec = _registry_backend_spec(name)
        if spec is not None and getattr(spec, "elements", None) is not None:
            covered = {_normalize_element(e) for e in spec.elements}
            for el in elements:
                if _normalize_element(el) not in covered:
                    errors.append(
                        f"element {el!r} is not covered by the {backend!r} model; "
                        f"choose a backend that covers it or split the structure"
                    )
    if str(device).strip().lower() in _GPU_DEVICES:
        errors.append(
            f"device {device!r} runs the MLIP on a GPU, which is only supported "
            "through the cluster submission path, never locally on this host"
        )
    if mult is not None and int(mult) < 1:
        errors.append(f"multiplicity {mult} is below the physical minimum of 1")
    return errors


def build_mlip_opt_input(
    work_dir: Path | str,
    *,
    run_type: str = "opt",
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int = 0,
    mult: int = 1,
    backend: str = "uma",
    model: Optional[str] = None,
    device: str = "cpu",
    hessian: Optional[str] = None,
    base_name: str = "extopt",
    ncores: int = 1,
    extra_keywords: Sequence[str] = (),
    **kwargs: Any,
) -> Path:
    """Write an ``! ExtOpt`` Opt input driven by an MLIP; return its ``.inp``.

    Refusals (see module docstring) raise :class:`ValueError` before any file
    is written.  Otherwise the input is delegated to
    :func:`~delfin.mlp_tools.orca_mlips.build_orca_mlips_input`, which writes
    the ``%method ProgExt / Ext_Params`` block with OPI's typed objects and
    the ``! Opt ExtOpt`` (or ``! Opt Freq ExtOpt``) simple-keyword line.

    ``run_type="opt_freq"`` builds the combined ``! Opt Freq ExtOpt`` run; the
    frequency block's Hessian is read back from the MLIP's dumped Hessian
    exactly as the shared builder's ``freq`` (``two_step``) route does
    (``%geom inhess Read``), driven by ``hessian="two_step"``.
    """
    work_dir = Path(work_dir)
    errors = validate_mlip_opt(
        run_type=run_type,
        backend=backend,
        device=device,
        elements=elements,
        charge=charge,
        mult=mult,
        check_installed=False,  # building an input never needs the MLIP present
    )
    if errors:
        raise ValueError("; ".join(errors))

    from delfin.mlp_tools.orca_mlips import build_orca_mlips_input

    rt = str(run_type).strip().lower()
    if rt == "opt_freq":
        keywords = ("Freq",) if "freq" not in {k.lower() for k in extra_keywords} else ()
        extra_keywords = tuple(extra_keywords) + keywords
    return build_orca_mlips_input(
        "opt",
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        backend=backend,
        model=model,
        device=device,
        hessian=hessian,
        base_name=base_name,
        work_dir=work_dir,
        extra_keywords=extra_keywords,
        ncores=ncores,
        **{k: v for k, v in kwargs.items() if k in _BUILDER_KWARGS},
    )


_BUILDER_KWARGS = frozenset(
    {
        "solvent",
        "solvent_model",
        "precision",
        "wrapper_path",
        "goat_options",
        "goat_variant",
        "neb_end_coordinates",
        "scan",
        "xyz_name",
    }
)


def _xyz_filename_from_inp(work_dir: Path, base_name: str) -> Optional[str]:
    """The geometry file named by the ``*xyzfile`` line of the input."""
    inp = work_dir / f"{base_name}.inp"
    if not inp.exists():
        return None
    try:
        text = inp.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    for line in text.splitlines():
        stripped = line.strip()
        # e.g. "* xyzfile 0 1 extopt.xyz" or "*xyzfile 0 1 extopt.xyz"
        if stripped.startswith("* xyzfile") or stripped.startswith("*xyzfile"):
            parts = stripped.split()
            if len(parts) >= 4:
                return parts[-1]
    return None


def read_mlip_opt(work_dir: Path | str, *, base_name: str = "extopt") -> dict:
    """Read an ``! Opt ExtOpt`` MLIP run: energy (Eh), geometry, convergence.

    Returns a dict with ``energy_Eh`` (float, Hartree = how the MLIP reports
    its ExtOpt energy, ORCA 6.1 optimizations.txt 4.1.10), ``energy_unit``
    (``"Eh"``), ``geometry`` (the ``*xyzfile``-named ``.xyz`` path), and the
    ``geometry_converged`` / ``terminated_normally`` flags.  ``error`` is
    ``None`` on success, otherwise a short reason.
    """
    from delfin.orca_read import Quantity, read

    work_dir = Path(work_dir)
    out_path = work_dir / f"{base_name}.out"
    out_text = ""
    if out_path.exists():
        try:
            out_text = out_path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            out_text = ""

    terminated = "ORCA TERMINATED NORMALLY" in out_text
    geom_converged = "THE OPTIMIZATION HAS CONVERGED" in out_text

    energy_Eh: Optional[float] = None
    energy_unit = "Eh"
    if terminated:
        reading = read(work_dir, base_name, Quantity.FINAL_SINGLE_POINT_ENERGY)
        value = getattr(reading, "value", None)
        if isinstance(value, (int, float)):
            energy_Eh = float(value)
            unit = getattr(reading, "unit", "")
            if unit:
                energy_unit = unit

    geometry = work_dir / (
        _xyz_filename_from_inp(work_dir, base_name)
        or f"{base_name}.xyz"
    )

    if not terminated:
        error = "the ORCA run did not terminate normally"
        if not out_text.strip():
            error = "no .out file to read"
    else:
        error = None
    return {
        "energy_Eh": energy_Eh,
        "energy_unit": energy_unit,
        "geometry": geometry,
        "geometry_converged": geom_converged,
        "terminated_normally": terminated,
        "error": error,
    }


class OrcaMlipOptOpiAdapter(StepAdapter):
    """Opt (and Opt+Freq) with an MLIP via ExtOpt, OPI-backed.

    Mirrors the other OPI blocks' registry position but drives an ``! ExtOpt``
    optimisation whose energy / gradients come from an MLIP (orca-mlips).
    Refusals (MLIP not installed, elements not covered, GPU locally,
    impossible multiplicity, unknown run type / backend) stop the step with a
    named reason; chemistry (method, charge, multiplicity) is never changed.
    """

    name = "orca_mlip_opt_opi"
    description = (
        "OPI-backed ORCA geometry optimisation with an MLIP via ExtOpt "
        "(energy/gradients from orca-mlips)"
    )
    category = "optimisation"
    produces_geometry = True
    consumes = ("geometry",)
    produces = ("qc_output", "geometry")
    requires_python = ("opi",)

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        kwargs.pop("_prev_artifacts", None)
        start = time.monotonic()
        try:
            elements = kwargs.get("elements")
            coordinates = kwargs.get("coordinates")
            if elements is None or coordinates is None:
                raise KeyError("elements/coordinates")
        except KeyError as exc:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=f"missing required {exc}",
                error_kind=ErrorKind.MISSING_PARAM,
            )

        run_type = kwargs.get("run_type", "opt")
        backend = kwargs.get("backend", "uma")
        device = kwargs.get("device", "cpu")
        base_name = kwargs.get("base_name", "extopt")

        try:
            errors = validate_mlip_opt(
                run_type=run_type,
                backend=backend,
                device=device,
                elements=elements,
                charge=kwargs.get("charge", 0),
                mult=kwargs.get("mult", 1),
            )
            if errors:
                raise ValueError("; ".join(errors))
            inp = build_mlip_opt_input(
                work_dir,
                run_type=run_type,
                elements=elements,
                coordinates=coordinates,
                charge=kwargs.get("charge", 0),
                mult=kwargs.get("mult", 1),
                backend=backend,
                model=kwargs.get("model"),
                device=device,
                base_name=base_name,
                ncores=cores,
            )
        except ValueError as exc:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=str(exc),
                error_kind=ErrorKind.INVALID_INPUT,
            )

        from delfin.orca import run_orca

        out_path = work_dir / f"{base_name}.out"
        ok = run_orca(
            str(inp),
            str(out_path),
            timeout=kwargs.get("timeout"),
            working_dir=work_dir,
        )
        res = read_mlip_opt(work_dir, base_name=base_name)

        return self._make_result(
            self.name,
            StepStatus.SUCCESS if ok and res.get("terminated_normally")
            else StepStatus.FAILED,
            work_dir,
            start,
            output_file=out_path,
            data={
                "meaningful": bool(res.get("terminated_normally")),
                "energy_Eh": res.get("energy_Eh"),
                "energy_unit": res.get("energy_unit"),
                "geometry": str(res.get("geometry")),
                "geometry_converged": res.get("geometry_converged"),
                "terminated_normally": res.get("terminated_normally"),
            },
            error=res.get("error"),
            error_kind=(ErrorKind.NONE if res.get("terminated_normally")
                        else ErrorKind.CONVERGENCE if res.get("geometry_converged") is False
                        else ErrorKind.TOOL_FAILED),
        )


register(OrcaMlipOptOpiAdapter())
