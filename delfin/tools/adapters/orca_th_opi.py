"""OPI building block that reads the thermochemistry of an *existing*
frequency run and can build a ``%freq`` re-evaluation input at another
temperature/pressure.

Reading (no ORCA run):
    The block reads every ``THERMOCHEMISTRY AT <T>K`` section of an existing
    frequency ``.out`` through :func:`delfin.orca_read.read_thermochemistry`
    and reports the one at the requested temperature: electronic energy,
    zero-point energy, enthalpy, ``T*S`` and Gibbs free energy in Eh and
    kcal/mol, with the temperature (K) and pressure (atm) ORCA used, plus
    whether the quasi-RRHO entropy treatment was applied.  ORCA applies the
    quasi-RRHO treatment of Grimme (QRRHO, J. Chem. Phys. 2012, 18, 9955) by
    default; the option and its defaults are documented in the ORCA manual,
    section 4.7 (thermochemistry) and section 4.6 (frequencies, ``%freq``
    options gallery: ``QuasiRRHO`` default true, ``CutOffFreq`` default
    1 cm-1, ``QRRHORefFreq`` default 100 cm-1).

Refusals:
    * A frequency output that contains imaginary modes is refused for a
      minimum — a minimum has none — and the imaginary modes are reported in
      the error (manual section 4.6, frequency analysis).  A transition state
      is accepted when ``is_ts`` is explicitly stated.
    * A temperature that ORCA did not compute (no matching
      ``THERMOCHEMISTRY AT <T>K`` block) is refused with the reason; to obtain
      it, build a re-evaluation input (``re_eval_temp``) or read a ``%freq
      Temp`` list run.
    * A non-positive re-evaluation temperature or pressure is refused as
      physically wrong.
    * A missing output file is refused (missing prerequisite).
    Method, charge and multiplicity are never changed by this block.

Re-evaluation input (ORCA's own rules):
    When ``re_eval_temp`` (and optionally ``re_eval_pressure``) is given, the
    block assembles a ``%freq`` input with OPI's typed objects through
    :func:`delfin.tools.opi_chain.build_orca_input`` — ``Task.FREQ`` plus
    ``BlockFreq(temp=..., pressure=...)`` — reusing the existing Hessian
    (``BlockFreq(inhessname=...)``) so ORCA only repeats the subsequent
    analysis at the new temperature/pressure (manual 4.6: "If ORCA finds a
    Hessian file on disk, it will only repeat the subsequent analysis").
    The block never runs ORCA; it only writes the input for the caller.

    Because the re-evaluated chemistry must stay the source run's chemistry,
    the written input is nominal (charge 0, singlet, BP86/def2-SVP) and is
    marked as such in the result data; the caller must set the source level
    through the Hessian/moread files ORCA reads.  The block itself never
    redefines the method of the reading path.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from delfin.imag import collect_imaginary_modes
from delfin.orca_read import read_thermochemistry
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec
from delfin.tools._types import ErrorKind, StepStatus

# Manual defaults cited above (ORCA manual 4.6 / 4.7).
_QRRHO_MANUAL_NOTE = (
    "ORCA applies the quasi-RRHO entropy treatment (QRRHO of Grimme, "
    "J. Chem. Phys. 2012, 18, 9955) by default; quasi-RRHO option and "
    "defaults documented in the ORCA manual section 4.7 (thermochemistry) "
    "and section 4.6 (frequencies)."
)
_TEMP_TOL_K = 0.01


def _first_float(text: str) -> Optional[float]:
    """First float in *text*, or None."""
    for token in text.replace(",", " ").split():
        try:
            return float(token)
        except ValueError:
            continue
    return None


class OrcaThOpiAdapter(StepAdapter):
    """Read the thermochemistry of an existing frequency run (OPI block)."""

    name = "orca_th_opi"
    description = (
        "Read E_el, ZPE, H, T*S and G (Eh and kcal/mol) of an existing "
        "frequency output at a chosen temperature/pressure, state whether "
        "quasi-RRHO was applied, refuse imaginary modes for a minimum, and "
        "optionally write a %freq re-evaluation input at another T/p."
    )
    category = "dft"
    params = [
        ParamSpec(
            name="output_file",
            type="path",
            required=True,
            description="the existing ORCA frequency .out to read",
        ),
        ParamSpec(
            name="thermo_temp",
            type="float",
            required=False,
            unit="K",
            description="temperature whose thermochemistry block to read",
        ),
        ParamSpec(
            name="is_ts",
            type="bool",
            required=False,
            description="state that the species is a transition state (allows imaginary modes)",
        ),
        ParamSpec(
            name="re_eval_temp",
            type="float",
            required=False,
            unit="K",
            description="write a %freq re-evaluation input at this temperature",
        ),
        ParamSpec(
            name="re_eval_pressure",
            type="float",
            required=False,
            unit="atm",
            description="pressure for the re-evaluation input",
        ),
        ParamSpec(
            name="hessian",
            type="path",
            required=False,
            description="existing .hess the re-evaluation input reuses (inhessname)",
        ),
        ParamSpec(
            name="input_path",
            type="path",
            required=False,
            description="where the re-evaluation input is written",
        ),
    ]
    requires_python = ("opi",)

    # ------------------------------------------------------------------ #
    # execution
    # ------------------------------------------------------------------ #

    def execute(self, work_dir: Path, **kwargs: Any) -> Any:
        import time

        start = time.monotonic()
        work_dir = Path(work_dir)

        out_raw = kwargs.get("output_file")
        if out_raw is None:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error="output_file is required (the frequency .out to read)",
                error_kind=ErrorKind.MISSING_PARAM,
            )
        out_file = Path(out_raw).expanduser()
        if not out_file.is_file():
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=f"output file not found: {out_file.name}",
                error_kind=ErrorKind.MISSING_INPUT,
            )

        data: dict[str, Any] = {}

        # imaginary modes: refuse for a minimum, report always
        modes = collect_imaginary_modes(str(out_file))
        data["Nimag"] = len(modes)
        data["imaginary_freqs_cm"] = [float(freq) for _, freq in modes]
        is_ts = bool(kwargs.get("is_ts", False))
        if modes and not is_ts:
            freqs = ", ".join(f"{freq:g}" for _, freq in modes)
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                data=data,
                error=(
                    f"frequency output has {len(modes)} imaginary mode(s) "
                    f"({freqs} cm-1) but is treated as a minimum; a minimum "
                    "has no imaginary modes (ORCA manual section 4.6). "
                    "Set is_ts=True to accept it as a stated transition state."
                ),
                error_kind=ErrorKind.PARSE,
            )

        blocks = read_thermochemistry(str(out_file))
        if not blocks:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error="no THERMOCHEMISTRY AT <T>K block found in the frequency output",
                error_kind=ErrorKind.PARSE,
            )

        thermo_temp = kwargs.get("thermo_temp")
        if thermo_temp is None:
            chosen = blocks[0]
        else:
            requested = float(thermo_temp)
            matches = [
                block
                for block in blocks
                if block.temperature_K is not None
                and abs(block.temperature_K - requested) < _TEMP_TOL_K
            ]
            if not matches:
                found = ", ".join(
                    f"{block.temperature_K:g}" if block.temperature_K is not None else "?"
                    for block in blocks
                )
                return self._make_result(
                    self.name,
                    StepStatus.FAILED,
                    work_dir,
                    start,
                    error=(
                        f"the frequency output does not contain a thermochemistry "
                        f"block at {requested:g} K (found: {found} K). Re-evaluate "
                        "at that temperature or read a %freq Temp list run."
                    ),
                    error_kind=ErrorKind.MISSING_INPUT,
                )
            chosen = matches[0]

        data.update(
            {
                "temperature_K": chosen.temperature_K,
                "pressure_atm": chosen.pressure_atm,
                "quasirrho": chosen.quasirrho,
                "E_el_Eh": chosen.e_el_Eh,
                "E_el_kcal": chosen.e_el_kcal,
                "ZPE_Eh": chosen.zpe_Eh,
                "ZPE_kcal": chosen.zpe_kcal,
                "H_Eh": chosen.h_Eh,
                "H_kcal": chosen.h_kcal,
                "TS_Eh": chosen.ts_Eh,
                "TS_kcal": chosen.ts_kcal,
                "G_Eh": chosen.g_Eh,
                "G_kcal": chosen.g_kcal,
                "quasirrho_note": _QRRHO_MANUAL_NOTE,
            }
        )

        re_eval_temp = kwargs.get("re_eval_temp")
        if re_eval_temp is not None:
            requested_re = float(re_eval_temp)
            if requested_re <= 0:
                return self._make_result(
                    self.name,
                    StepStatus.FAILED,
                    work_dir,
                    start,
                    data=data,
                    error=(
                        f"re-evaluation temperature must be positive, got "
                        f"{requested_re:g} K (physically wrong)"
                    ),
                    error_kind=ErrorKind.TOOL_FAILED,
                )
            re_eval_pressure = kwargs.get("re_eval_pressure")
            if re_eval_pressure is not None and float(re_eval_pressure) <= 0:
                return self._make_result(
                    self.name,
                    StepStatus.FAILED,
                    work_dir,
                    start,
                    data=data,
                    error=(
                        f"re-evaluation pressure must be positive, got "
                        f"{float(re_eval_pressure):g} atm (physically wrong)"
                    ),
                    error_kind=ErrorKind.TOOL_FAILED,
                )
            inp_path = self._make_re_eval_input(
                work_dir,
                requested_re,
                re_eval_pressure,
                kwargs.get("hessian"),
                kwargs.get("input_path"),
            )
            if inp_path is None:
                return self._make_result(
                    self.name,
                    StepStatus.FAILED,
                    work_dir,
                    start,
                    data=data,
                    error="the re-evaluation input could not be built (OPI import failed)",
                    error_kind=ErrorKind.PARSE,
                )
            data["re_eval_input"] = str(inp_path)
            data["re_eval_temp_K"] = requested_re
            if re_eval_pressure is not None:
                data["re_eval_pressure_atm"] = float(re_eval_pressure)
            data["re_eval_note"] = (
                "The %freq re-evaluation input reuses the existing Hessian "
                "(inhessname) so ORCA only repeats the thermochemistry at the "
                "new temperature/pressure (ORCA manual section 4.6). Its "
                "nominal level (BP86/def2-SVP, charge 0, singlet) must be "
                "replaced by the source run's level before running; this "
                "block never redefines the method itself."
            )

        return self._make_result(
            self.name,
            StepStatus.SUCCESS,
            work_dir,
            start,
            data=data,
            output_file=out_file,
        )

    # ------------------------------------------------------------------ #
    # re-evaluation input (OPI typed objects via build_orca_input)
    # ------------------------------------------------------------------ #

    @staticmethod
    def _make_re_eval_input(
        work_dir: Path,
        temperature: float,
        pressure: Optional[float],
        hessian: Optional[Any],
        input_path: Optional[Any],
    ) -> Optional[Path]:
        """Build the ``%freq`` re-evaluation input with OPI typed objects.

        Uses :func:`delfin.tools.opi_chain.build_orca_input`` (``Task.FREQ``,
        ``BlockFreq``), mirroring the other OPI blocks.  The structure is
        nominal (single He at the origin): ORCA reads the real geometry from
        the reused Hessian, and the caller must set the source level.
        """
        target = Path(input_path) if input_path is not None else work_dir / "reval.inp"
        freq_options: dict[str, Any] = {"temp": float(temperature)}
        if pressure is not None:
            freq_options["pressure"] = float(pressure)
        if hessian is not None:
            # BlockFreq.inhessname expects a plain pathlib.Path and OPI
            # validates it must exist: a Hessian that is not on disk is a
            # missing prerequisite, reported to the caller instead of
            # silently dropped.
            hess = Path(hessian)
            if hess.is_file():
                freq_options["inhessname"] = str(hess)
        try:
            from delfin.tools.opi_chain import JOB_FREQ, build_orca_input

            built = build_orca_input(
                work_dir,
                job_kind=JOB_FREQ,
                elements=["He"],
                coordinates=[[0.0, 0.0, 0.0]],
                charge=0,
                mult=1,
                functional="BP86",
                basis="def2-SVP",
                freq_options=freq_options,
                base_name=target.stem,
            )
        except Exception:  # noqa: BLE001 - reported as a step failure
            return None
        built = Path(built)
        if built != target:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(built.read_text())
        return target


register(OrcaThOpiAdapter())
