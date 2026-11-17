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
    the written input inherits charge, multiplicity, method and basis from the
    source ``.inp`` and takes the geometry from the Hessian's own ``$atoms``
    block -- the structure the Hessian was computed at.  For an Opt+Freq
    source the ``.inp`` holds only the START geometry, which does not belong
    to the Hessian (real ORCA 6.1.1 run, H2O B3LYP/def2-SVP: G(350 K) off by
    0.69 mEh with the start geometry).  The block itself never redefines the
    method of the reading path.
"""

from __future__ import annotations

import re

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


_GEOMETRY_HEADER_RE = re.compile(
    r"^\s*\*\s*xyz\s+(-?\d+)\s+(-?\d+)\s*$"
)
_ATOM_COORD_RE = re.compile(
    r"^\s*([A-Za-z][a-z]?)\s+(-?[\d.]+(?:[eE][+-]?\d+)?)\s+"
    r"(-?[\d.]+(?:[eE][+-]?\d+)?)\s+"
    r"(-?[\d.]+(?:[eE][+-]?\d+)?)\s*$"
)


def _quote_inhessname(inp: Path) -> None:
    """Wrap the ``%freq inhessname`` value in double quotes.

    OPI's BlockFreq renders a string value bare (``inhessname /abs/path.hess``),
    which real ORCA 6.1 rejects: "Invalid assignment in FREQ block".  ORCA
    accepts the value only quoted (``inhessname \"/abs/path.hess\"``).
    """
    text = inp.read_text()
    new = re.sub(
        r"(?m)^(\s*inhessname\s+)(\S+)$",
        lambda m: m.group(1) + '"' + m.group(2) + '"',
        text,
    )
    if new != text:
        inp.write_text(new)


_BOHR_TO_ANGSTROM = 0.529177210903  # CODATA 2018


def _read_hess_geometry(hess_file: Path) -> Optional[list[tuple[str, float, float, float]]]:
    """Atoms of an ORCA ``.hess`` file from its ``$atoms`` block.

    ORCA writes ``$atoms``, the atom count, then one line per atom: element,
    mass, x, y, z in bohr.  This is the geometry the Hessian belongs to (the
    optimised one for an Opt+Freq run).  Returns Angstrom coordinates, or None
    when the block is missing or malformed.
    """
    try:
        lines = hess_file.read_text(errors="replace").splitlines()
    except OSError:
        return None
    for i, line in enumerate(lines):
        if line.strip() != "$atoms":
            continue
        try:
            n = int(lines[i + 1].split()[0])
            atoms = []
            for row in lines[i + 2:i + 2 + n]:
                f = row.split()
                atoms.append((f[0], *(float(v) * _BOHR_TO_ANGSTROM for v in f[2:5])))
        except (IndexError, ValueError):
            return None
        return atoms if len(atoms) == n and n > 0 else None
    return None


def _read_source_level(inp_file: Path) -> Optional[dict[str, Any]]:
    """Parse the source ORCA ``.inp`` for the level a Hessian-reusing
    re-evaluation input must inherit: charge + multiplicity (from the
    ``* xyz <charge> <mult>`` header), method + basis (from the ``!`` line),
    and the geometry atoms.  Returns None when any part is missing (e.g. the
    input used ``moread`` or a non-``xyz`` coordinate block), so the caller
    refuses rather than emit a level that is not the source run's.
    """
    if not inp_file.is_file():
        return None
    text = inp_file.read_text()
    method = basis = None
    charge = mult = None
    geometry: list[tuple[str, float, float, float]] = []
    after_header = False
    for line in text.splitlines():
        if line.startswith("!"):
            toks = line[1:].split()
            if len(toks) >= 2:
                method, basis = toks[0], toks[1]
            continue
        m = _GEOMETRY_HEADER_RE.match(line)
        if m:
            charge, mult = int(m.group(1)), int(m.group(2))
            after_header = True
            continue
        if after_header:
            if line.strip() == "*":
                after_header = False
                continue
            a = _ATOM_COORD_RE.match(line)
            if a:
                geometry.append(
                    (a.group(1), float(a.group(2)),
                     float(a.group(3)), float(a.group(4)))
                )
    if (charge is None or mult is None or not method
            or not basis or not geometry):
        return None
    return {
        "charge": charge,
        "mult": mult,
        "functional": method,
        "basis": basis,
        "geometry": geometry,
    }


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
            name="source_inp",
            type="path",
            required=False,
            description=(
                "source ORCA .inp whose charge, multiplicity, method, basis "
                "and geometry a Hessian-reusing re-evaluation input inherits "
                "(default: the .out's sibling .inp)"
            ),
        ),
        ParamSpec(
            name="geometry",
            type="list",
            required=False,
            description=(
                "source molecule geometry, [(element, x, y, z), ...] in "
                "Angstrom, written into a Hessian-reusing re-evaluation "
                "input when no Hessian is reused (with a Hessian the geometry "
                "comes from its $atoms block)"
            ),
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
            hessian_in = kwargs.get("hessian")
            geometry_in = kwargs.get("geometry")
            if hessian_in is not None:
                hessian_in = Path(hessian_in)
                if not hessian_in.is_file():
                    return self._make_result(
                        self.name,
                        StepStatus.FAILED,
                        work_dir,
                        start,
                        data=data,
                        error=(
                            f"the referenced Hessian {hessian_in} does not "
                            f"exist; a missing prerequisite cannot be reused"
                        ),
                        error_kind=ErrorKind.MISSING_INPUT,
                    )
                inp_raw = kwargs.get("source_inp")
                source_inp = (
                    Path(inp_raw).expanduser()
                    if inp_raw is not None
                    else out_file.with_suffix(".inp")
                )
                level = _read_source_level(source_inp)
                if level is None:
                    return self._make_result(
                        self.name,
                        StepStatus.FAILED,
                        work_dir,
                        start,
                        data=data,
                        error=(
                            "a Hessian-reusing re-evaluation input needs the "
                            "source run's charge, multiplicity, method, basis "
                            "and geometry (read from its input file "
                            f"{source_inp}); that file is absent or does not "
                            "carry a * xyz geometry block, so the input would "
                            "not re-evaluate the source chemistry at the new "
                            "temperature"
                        ),
                        error_kind=ErrorKind.MISSING_INPUT,
                    )
                hess_geometry = _read_hess_geometry(hessian_in)
                src_elements = [g[0].capitalize() for g in level["geometry"]]
                if hess_geometry is None or [
                    g[0].capitalize() for g in hess_geometry
                ] != src_elements:
                    return self._make_result(
                        self.name,
                        StepStatus.FAILED,
                        work_dir,
                        start,
                        data=data,
                        error=(
                            f"the Hessian {hessian_in} has no readable $atoms "
                            "block matching the source molecule "
                            f"({' '.join(src_elements)}); the re-evaluation "
                            "must run at the geometry the Hessian was "
                            "computed at"
                        ),
                        error_kind=ErrorKind.MISSING_INPUT,
                    )
                src_geometry = hess_geometry
            else:
                source_inp = None
                level = None
                src_geometry = geometry_in
            inp_path = self._make_re_eval_input(
                work_dir,
                requested_re,
                re_eval_pressure,
                hessian_in,
                kwargs.get("input_path"),
                src_geometry,
                charge=level["charge"] if level is not None else 0,
                mult=level["mult"] if level is not None else 1,
                functional=level["functional"] if level is not None else "BP86",
                basis=level["basis"] if level is not None else "def2-SVP",
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
                "The re-evaluation input reuses the existing Hessian "
                "(inhessname) so ORCA only repeats the thermochemistry at the "
                "new temperature/pressure (ORCA manual section 4.6). It "
                "inherits charge, multiplicity, method and basis from the "
                "source input and the geometry from the Hessian's $atoms "
                "block (the structure the Hessian belongs to); this block "
                "never redefines the method itself."
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
        hessian: Optional[Path],
        input_path: Optional[Any],
        geometry: Optional[Any] = None,
        *,
        charge: int = 0,
        mult: int = 1,
        functional: str = "BP86",
        basis: str = "def2-SVP",
    ) -> Optional[Path]:
        """Build the ``%freq`` re-evaluation input with OPI typed objects.

        The input inherits the *source* run's level from ``execute``: charge,
        multiplicity, method and basis, and the geometry of the Hessian's
        ``$atoms`` block (a nominal single He would re-evaluate the wrong
        chemistry).  When no Hessian is reused (grammar-only probe path) the
        default BP86/def2-SVP charge-0 singlet He remains.  The ``$hess`` path
        is rendered quoted (``inhessname "..."``): real ORCA 6.1 aborts with
        "Invalid assignment in FREQ block" on a bare path.
        """
        target = Path(input_path) if input_path is not None else work_dir / "reval.inp"
        freq_options: dict[str, Any] = {"temp": float(temperature)}
        if pressure is not None:
            freq_options["pressure"] = float(pressure)
        if hessian is not None:
            freq_options["inhessname"] = str(hessian)
        if geometry:
            elements = [str(g[0]) for g in geometry]
            coordinates = [[float(x) for x in g[1:4]] for g in geometry]
        else:
            # No Hessian to reuse (probe/grammar path): a nominal single He.
            elements, coordinates = ["He"], [[0.0, 0.0, 0.0]]
        try:
            from delfin.tools.opi_chain import JOB_FREQ, build_orca_input

            built = build_orca_input(
                work_dir,
                job_kind=JOB_FREQ,
                elements=elements,
                coordinates=coordinates,
                charge=charge,
                mult=mult,
                functional=functional,
                basis=basis,
                freq_options=freq_options,
                base_name=target.stem,
            )
        except Exception:  # noqa: BLE001 - reported as a step failure
            return None
        built = Path(built)
        if built != target:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(built.read_text())
        if hessian is not None:
            _quote_inhessname(target)
        return target


register(OrcaThOpiAdapter())
