"""Registered step adapter for the OPI-backed CO2 coordinator input builder.

This adapter makes ``delfin.tools.co2_coordinator_opi`` reachable through
the step registry (``run_step`` / ``list_adapters``): it writes the same
ORCA ``.inp`` files the classic CO2 Coordinator writes, but with OPI's
input classes.  The classic coordinator remains the workflow's driver —
this step is the parity-proof reference for its inputs.

The step builds inputs only; it does not launch ORCA.  It requires OPI
(``requires_python=("opi",)``), imported lazily so the module imports
without OPI present.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec
from delfin.tools._types import StepResult, StepStatus

logger = get_logger(__name__)


class Co2CoordinatorOpiAdapter(StepAdapter):
    """Build the CO2-coordinator ORCA inputs with OPI's input classes.

    Mirrors ``Co2Coordinator6``'s two input writers (orientation
    single-points and the relaxed distance scan) into a registered step.
    The inputs written are equivalent to the classic writer's — the parity
    test ``tests/test_opi_co2_coordinator_parity.py`` proves it.
    """

    name = "co2_coordinator_opi"
    description = "OPI-backed ORCA input builder for the CO2 coordinator"
    category = "dft"
    produces_geometry = False
    consumes = ("geometry",)
    produces = ("orca_input",)
    requires_python = ("opi",)
    #: ``job`` is declared so the choice is DISCOVERABLE: without it
    #: ``delfin-step co2_coordinator_opi --describe`` printed no params at all,
    #: and a caller could not learn that the relaxed-scan writer exists.  It is
    #: optional with the historical default, so ``_validate``'s required-param
    #: check is untouched; its enum is advisory (a mismatched value raises a
    #: WARNING, and the adapter's own gate is the hard refusal).
    params = (
        ParamSpec(
            "job", "str", required=False, default="sp",
            enum=("sp", "scan"),
            description=("Which classic CO2 input to write: the orientation "
                         "single points (sp) or the relaxed distance scan"
                         " (scan)")),
        # Declared so the alternative the adapter names actually works from
        # the CLI: a caller that names both indices cannot also pass an
        # ase.Atoms through argv, and without this the only route left was a
        # keyword no `--describe` mentioned.
        ParamSpec(
            "start_distance", "float", required=False, default=None,
            unit="A",
            description=("Where the relaxed scan starts (metal-C, Angstrom). "
                         "Measured from the atoms object when omitted, which "
                         "needs an ase.Atoms rather than a symbol")),
    )

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        **kwargs: Any,
    ) -> StepResult:
        del geometry  # classic writers take atoms + paths as keyword args
        start = time.monotonic()
        # Which of the classic coordinator's two inputs to write.  The default
        # is the single point, so a caller that names no job gets exactly what
        # it got before ``scan`` was mounted: same writer, same kwargs, same
        # file name.  The gate runs BEFORE the scan writer is imported, so
        # refusing an unknown job costs nothing and needs neither OPI nor the
        # classic module.
        job = kwargs.pop("job", "sp")
        if job not in ("sp", "scan"):
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=f"unknown job {job!r}; use one of ['scan', 'sp']",
            )
        if job == "scan":
            return self._execute_scan(work_dir, start, kwargs)
        try:
            from delfin.tools.co2_coordinator_opi import (
                write_co2_sp_opi_input,
            )

            outdir = kwargs.get("outdir", work_dir / "calc")
            tag = kwargs.get("tag", "calc")
            inp_path = write_co2_sp_opi_input(
                atoms=kwargs["atoms"],
                xyz_path=str(kwargs.get("xyz_path", outdir / "input.xyz")),
                outdir=str(outdir),
                orca_keywords=kwargs.get("orca_keywords", "PBE0 def2-SVP SP"),
                broken_sym=kwargs.get("broken_sym", ""),
                charge=kwargs.get("charge", 0),
                multiplicity=kwargs.get("mult", 1),
                PAL=kwargs.get("PAL", 4),
                maxcore=kwargs.get("maxcore", 1000),
                tag=tag,
                metal_symbol=kwargs.get("metal_symbol"),
                metal_basis=kwargs.get("metal_basis"),
                control_args=kwargs.get("control_args"),
                qmmm_range=kwargs.get("qmmm_range"),
                co2_indices=kwargs.get("co2_indices"),
            )
        except KeyError as exc:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=f"missing required argument {exc}",
            )
        return self._make_result(
            self.name,
            StepStatus.SUCCESS,
            work_dir,
            start,
            output_file=inp_path,
            data={"input_file": str(inp_path)},
        )


    def _execute_scan(
        self, work_dir: Path, start: float, kwargs: dict
    ) -> StepResult:
        """Write the relaxed CO2 distance-scan input (``job="scan"``).

        The classic coordinator does not pick the scan axis from a caller's
        index list: after choosing the best orientation it finds the metal with
        ``detect_metal_index`` and takes the metal-C distance of THAT structure
        as the scan's start (``CO2_Coordinator6.py:1496-1498``), with the CO2
        carbon located by symbol (``_co2_axis_center_indices``,
        ``CO2_Coordinator6.py:605-614``).  An explicit ``metal_index`` /
        ``co2_c_index`` still wins, so a caller that knows its own axis is not
        second-guessed.

        The indices must not be inferred by POSITION.  DELFIN's own ``co2.xyz``
        is written O, C, O (``CO2_Coordinator6.py:137-141``,
        ``chain_setup.py:201-207``), so ``co2_indices[0]`` is an OXYGEN; taking
        it as the carbon would scan the M-O distance and label the result a
        C scan.  The carbon is therefore the single carbon among
        ``co2_indices``, and an ambiguous substrate is refused rather than
        guessed.
        """
        atoms = kwargs.get("atoms")
        if atoms is None:
            return self._refuse_scan(work_dir, start, "missing required "
                                     "argument 'atoms'")
        co2_indices = kwargs.get("co2_indices")
        try:
            metal_index = kwargs.get("metal_index")
            co2_c_index = kwargs.get("co2_c_index")
            if metal_index is None or co2_c_index is None:
                # Only the INFERENCE path reads atoms, so only that path
                # requires an atoms object.  ``delfin-step --atoms=Ir`` arrives
                # as a string (argv cannot carry an ase.Atoms); a caller that
                # names both indices and the start distance never needs one.
                _require_atoms_object(atoms)
            if metal_index is None:
                from delfin.co2.CO2_Coordinator6 import detect_metal_index

                metal_index = detect_metal_index(atoms)
            co2_c_index = kwargs.get("co2_c_index")
            if co2_c_index is None:
                co2_c_index = _co2_carbon_index(atoms, co2_indices)
            start_distance = kwargs.get("start_distance")
            if start_distance is None:
                # The scan's start is the metal-C distance of the structure,
                # which only an atoms object can give.  A caller naming both
                # indices has been told it "never needs" one (:169-173) --
                # true only when it also names start_distance, which the
                # message below says.  Without this check the same caller got
                # `AttributeError: 'str' object has no attribute 'position'`
                # out of _distance, the one failure shape a caller cannot act
                # on and the only one that does not name its missing input.
                if not hasattr(atoms, "positions"):
                    raise ValueError(
                        "a scan needs a start distance: name start_distance, "
                        "or pass an ase.Atoms through run_step (the step "
                        f"registry) so it can be measured -- 'atoms' was "
                        f"{type(atoms).__name__}")
                start_distance = _distance(atoms, metal_index, co2_c_index)
        except ValueError as exc:
            return self._refuse_scan(work_dir, start, str(exc))

        try:
            from delfin.tools.co2_coordinator_opi import (
                write_co2_scan_opi_input,
            )

            scan_dir = kwargs.get("scan_dir", work_dir / "relaxed_surface_scan")
            inp_path = write_co2_scan_opi_input(
                atoms,
                str(kwargs.get("xyz_path", work_dir / "input.xyz")),
                int(metal_index),
                int(co2_c_index),
                start_distance=float(start_distance),
                end_distance=kwargs.get("end_distance", 1.7),
                steps=kwargs.get("steps", 5),
                orca_keywords=kwargs.get(
                    "orca_keywords", "GFN2-XTB OPT ALPB(DMF)"),
                broken_sym=kwargs.get("broken_sym", ""),
                charge=kwargs.get("charge", -2),
                multiplicity=kwargs.get("mult", 3),
                PAL=kwargs.get("PAL", 4),
                maxcore=kwargs.get("maxcore", 2000),
                metal_symbol=kwargs.get("metal_symbol"),
                metal_basis=kwargs.get("metal_basis"),
                control_args=kwargs.get("control_args"),
                qmmm_range=kwargs.get("qmmm_range"),
                co2_indices=co2_indices,
                scan_dir=str(scan_dir),
            )
        except KeyError as exc:
            return self._refuse_scan(work_dir, start,
                                     f"missing required argument {exc}")
        return self._make_result(
            self.name,
            StepStatus.SUCCESS,
            work_dir,
            start,
            output_file=inp_path,
            data={"input_file": str(inp_path)},
        )

    def _refuse_scan(self, work_dir: Path, start: float, error: str):
        return self._make_result(
            self.name, StepStatus.FAILED, work_dir, start, error=error)


def _require_atoms_object(atoms) -> None:
    """Refuse a non-atoms ``atoms`` with a message that says what to pass."""
    if not hasattr(atoms, "symbols") or not hasattr(atoms, "positions"):
        raise ValueError(
            f"inferring the scan axis needs an ase.Atoms for 'atoms', got "
            f"{type(atoms).__name__}; pass an atoms object through run_step "
            f"(the step registry), or name metal_index, co2_c_index and "
            f"start_distance explicitly")


def _co2_carbon_index(atoms, co2_indices) -> int:
    """The CO2 carbon's index, found by SYMBOL inside ``co2_indices``.

    By symbol, not by position: DELFIN's own ``co2.xyz`` is O, C, O, so the
    first index of the CO2 group is an oxygen.  A group without exactly one
    carbon is ambiguous and is refused rather than guessed.
    """
    if not co2_indices:
        raise ValueError(
            "a scan needs the CO2 carbon: pass co2_indices (or co2_c_index)")
    group = [int(i) for i in co2_indices]
    carbons = [i for i in group if atoms[i].symbol == "C"]
    if len(carbons) != 1:
        raise ValueError(
            f"a scan needs exactly one carbon in co2_indices {group}; found "
            f"{len(carbons)} -- pass co2_c_index to name it")
    return carbons[0]


def _distance(atoms, i: int, j: int) -> float:
    """``|r_i - r_j|`` in Angstrom, the classic's scan start
    (``CO2_Coordinator6.py:1498``)."""
    from math import dist

    return dist(tuple(atoms[int(i)].position), tuple(atoms[int(j)].position))


register(Co2CoordinatorOpiAdapter())
