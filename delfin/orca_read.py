"""One place where DELFIN reads a quantity from a finished ORCA run.

Phase 1 (tests/test_opi_equivalence_q04.py and tests/test_opi_equivalence_q05.py)
measured ORCA's Python interface OPI (``pip install "delfin[orca]"``) against
DELFIN's own text readers on ORCA 6.1.1 runs: identical values for the final
single-point energy, the Gibbs free energy, ZPE, total enthalpy, inner energy,
the <S**2> deviation, multiplicity, charge and normal termination. This module
is the single reader that uses OPI for those quantities and DELFIN's
established text parsers for everything else (or wherever OPI is not
available). Callers keep their own parsers until they switch here, one PR each.

Rules, in order:

1. OPI is used only when: it is importable, the ``.out`` header says
   ORCA >= 6.1.1, the property JSON OPI reads is present, and the run
   terminated normally with a converged SCF (both checked on the text,
   the same way OPI checks them). Otherwise the text parser answers.
2. Shadow check: when both readers give a value, they must agree within
   the tolerance of the phase-1 tests (1e-6 Eh for energies, 1e-6 for
   <S**2> -- the text prints 8 resp. 6 decimals). On a mismatch the TEXT
   value is returned (the established one) with a warning naming the
   quantity, both values and the run's base name.
3. Setting ``orca.reader`` in ``~/.delfin_settings.json``: ``"auto"``
   (default, rules above) or ``"text"`` (never OPI).
4. A missing or partial run never raises: value ``None``, source
   ``"none"``, and the note says why.

Quantities deliberately NOT routed through OPI (measured in phase 1):

* failure reason -- OPI says THAT a run failed, never WHY;
* geometry convergence -- DELFIN answers ``None`` for a run without an
  optimisation, OPI answers ``False`` (absent "HURRAY" marker), which
  would turn "not applicable" into "failed";
* entropy -- both sides give T*S in Eh, not the molar entropy; the
  naming decision is pending, so the text reader stays;
* number of electrons -- DELFIN derives it from the input geometry and
  charge; OPI reads it from the property JSON, which a silently
  converted older run may not carry.

The optional helper :func:`ensure_property_json` may create the property
JSON with ``orca_2json`` -- but only inside the run's own directory, and
only when ``orca_2json`` is on PATH. It never writes anywhere else.
"""
from __future__ import annotations

import json
import logging
import re
import shutil
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Optional, Union

logger = logging.getLogger(__name__)

#: The .out prints <S**2> to 6 and energies to 8 decimals; OPI's JSON
#: carries full doubles. The phase-1 equivalence tests compared at the
#: text's printed precision; the shadow check uses the same.
SHADOW_TOL = 1e-6

#: Geometry shadow-check tolerance, Angstrom (see _geometries_agree).
GEOMETRY_SHADOW_TOL = 1e-6

#: Minimum ORCA version OPI supports (its own contract).
_MIN_ORCA_VERSION = (6, 1, 1)

_VERSION_RE = re.compile(r"Program\s+Version\s+(\d+)\.(\d+)(?:\.(\d+))?")


class Quantity(Enum):
    """A quantity DELFIN reads from a finished ORCA run."""

    #: final single point (electronic) energy, Eh
    FINAL_SINGLE_POINT_ENERGY = "final_single_point_energy"
    #: final Gibbs free energy G, Eh
    GIBBS = "gibbs"
    #: zero-point energy (correction), Eh
    ZPE = "zpe"
    #: total enthalpy H, Eh
    ENTHALPY = "enthalpy"
    #: inner energy U (total thermal energy), Eh
    INNER_ENERGY = "inner_energy"
    #: entropy term as printed by ORCA: T*S, Eh (NOT the molar entropy)
    ENTROPY = "entropy"
    #: <S**2> deviation (expected - ideal); 0.0 when closed shell
    S2_DEVIATION = "s2_deviation"
    #: spin multiplicity
    MULTIPLICITY = "multiplicity"
    #: total charge
    CHARGE = "charge"
    #: total number of electrons
    NUMBER_OF_ELECTRONS = "number_of_electrons"
    #: True when the .out ends with ORCA TERMINATED NORMALLY
    TERMINATED_NORMALLY = "terminated_normally"
    #: why the run failed (outcome phrase); None when it did not fail
    FAILURE_REASON = "failure_reason"
    #: geometry-optimisation converged; None when the run did not optimise
    GEOMETRY_CONVERGED = "geometry_converged"
    #: the run's final geometry, [(element, x, y, z), ...] in Angstrom
    FINAL_GEOMETRY = "final_geometry"
    #: IR modes of the IR SPECTRUM block (list of delfin.ir_spectrum.IRMode)
    IR_MODES = "ir_modes"
    #: (mode index, frequency) imaginary modes printed in the .out,
    #: most negative first
    IMAGINARY_MODES = "imaginary_modes"
    #: all $vibrational_frequencies of a .hess file (imaginary = negative)
    VIBRATIONAL_FREQUENCIES = "vibrational_frequencies"


@dataclass(frozen=True)
class Reading:
    """One quantity read from a finished ORCA run.

    ``source`` is ``"opi"`` (ORCA's structured output via OPI),
    ``"text"`` (DELFIN's own text parsers) or ``"none"`` (nothing was
    read -- the note says why; never an error).
    """

    value: Union[float, int, bool, str, None]
    unit: str
    source: str
    note: str = ""


# ---------------------------------------------------------------------------
# settings
# ---------------------------------------------------------------------------

def _reader_mode() -> str:
    """``orca.reader`` from the user settings: "auto" (default) or "text"."""
    try:
        from delfin.user_settings import load_settings

        settings = load_settings()
        mode = str((settings.get("orca") or {}).get("reader") or "auto").strip().lower()
    except Exception:
        return "auto"
    return mode if mode in ("auto", "text") else "auto"


# ---------------------------------------------------------------------------
# OPI gate: importability, version, property files, run health
# ---------------------------------------------------------------------------

def _opi_importable() -> bool:
    """Whether OPI is installed (imported lazily; absence is normal)."""
    try:
        import opi.output.core  # noqa: F401
    except Exception:
        return False
    return True


def _out_header_version(text: str) -> Optional[tuple]:
    """ORCA version from the .out header, e.g. (6, 1, 1); None if absent."""
    m = _VERSION_RE.search(text)
    if not m:
        return None
    return (int(m.group(1)), int(m.group(2)), int(m.group(3) or 0))


def _version_at_least(version: Optional[tuple], minimum: tuple) -> bool:
    return version is not None and version >= minimum


def _head_text(path: Path, size: int = 65536) -> str:
    """The first *size* bytes of a text file ('' when unreadable)."""
    try:
        with path.open(encoding="utf-8", errors="replace") as f:
            return f.read(size)
    except OSError:
        return ""


#: How much of the output's tail a single-quantity text read touches. The
#: termination marker, the final single-point energy, the thermochemistry
#: block and the convergence/outcome phrases all live in the last megabyte
#: of an ORCA run (the thermo block of the h2o_optfreq fixture sits about
#: 140 kB from the end; a 500-step optimization adds only a few kB per
#: cycle), so 1 MiB carries them at ~2.2 MB tracemalloc peak instead of
#: the whole 100-500 MB file.
_TAIL_BYTES = 1 * 1024 * 1024


def _tail_text(path: Path, size: int = _TAIL_BYTES) -> str:
    """The last *size* bytes of a text file ('' when unreadable).

    Binary seek+read like doc_server.calc_indexer._tail_text, not a text
    mode f.read(): the text-mode read of a 2 MB tail allocated ~8.5 MB
    (decoder buffering), the binary one stays near the tail's own size.
    """
    try:
        with path.open("rb") as fh:
            fh.seek(0, 2)
            n = fh.tell()
            fh.seek(max(0, n - size))
            return fh.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def _terminated_normally_text(out_path: Path) -> bool:
    """DELFIN's text read of normal termination (delfin.smart_recalc)."""
    from delfin.smart_recalc import has_ok_marker

    return has_ok_marker(str(out_path))


def _scf_converged_text(text: str) -> bool:
    """DELFIN's text read of SCF convergence (delfin.api)."""
    return "SCF CONVERGED" in text


def _opi_ready(work_dir: Path, base_name: str,
               out_path: Path, text: str) -> tuple:
    """(opi_usable, reason) -- the gate of rule 1, each step cheap."""
    if not _opi_importable():
        return False, "OPI not installed"
    version = _out_header_version(_head_text(out_path))
    if not _version_at_least(version, _MIN_ORCA_VERSION):
        return False, (f"ORCA version {version} in output header "
                       f"predates the minimum {'.'.join(map(str, _MIN_ORCA_VERSION))}")
    if not (work_dir / f"{base_name}.property.json").is_file():
        return False, "property JSON missing"
    if not (work_dir / f"{base_name}.property.txt").is_file():
        return False, "property txt missing"
    if not _terminated_normally_text(out_path):
        return False, "run did not terminate normally"
    if not _scf_converged_text(text):
        return False, "SCF not converged"
    return True, ""


# ---------------------------------------------------------------------------
# OPI readers
# ---------------------------------------------------------------------------

_OPI_GETTERS = {
    "final_energy": "get_final_energy",
    "free_energy": "get_free_energy",
    "zpe": "get_zpe",
    "enthalpy": "get_enthalpy",
    "inner_energy": "get_inner_energy",
    "entropy": "get_entropy",
}


def _opi_values(work_dir: Path, base_name: str) -> dict:
    """All OPI reads for one run in a single parse; ``None`` when absent.

    ``{"error": ...}`` when OPI could not parse the run at all; every
    OPI failure degrades to the text parser, never to an exception.
    """
    from opi.output.core import Output

    result: dict = {}
    try:
        out = Output(base_name, working_dir=work_dir, version_check=False)
        out.parse()
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}

    def get(getter_name, kind=float):
        try:
            v = getattr(out, getter_name)()
        except Exception:
            return None
        return None if v is None else kind(v)

    for key, getter in _OPI_GETTERS.items():
        result[key] = get(getter)
    result["charge"] = get("get_charge", int)
    result["mult"] = get("get_mult", int)

    try:
        s2 = out.get_s2()
    except Exception:
        s2 = None
    result["s2"] = s2

    result["geometry"] = None
    try:
        structure = out.get_structure()
        if structure is not None and structure.atoms:
            result["geometry"] = [
                (a.element, float(a.coordinates.x), float(a.coordinates.y),
                 float(a.coordinates.z)) for a in structure.atoms]
    except Exception:
        pass
    return result


def _opi_final_single_point_energy(work_dir: Path, base_name: str) -> Optional[float]:
    """One OPI final-energy read; kept as a narrow seam for callers that
    only need this one quantity without parsing the whole run."""
    return _opi_values(work_dir, base_name).get("final_energy")


# ---------------------------------------------------------------------------
# text readers (DELFIN's established parsers, called -- never copied)
# ---------------------------------------------------------------------------

def _text_values(out_path: Path, inp_path: Path,
                 quantities: Optional[set] = None,
                 tail: Optional[str] = None) -> dict:
    """Text reads for one run, each quantity from only the region it needs.

    Every parser here reads either the bounded tail (``_TAIL_BYTES``), a
    pre-read bounded text, or the small ``.inp`` -- never the whole output,
    which used to cost ~200 MB peak on a 12 MB file. *quantities* limits
    the work to the keys the caller wants (all keys when omitted); *tail*
    passes in an already-read bounded tail so one ``read()`` call scans
    the output once. The values are identical to the whole-file readers
    because the markers they match live near the end of an ORCA run.
    """
    wanted = quantities or {
        "final_energy", "free_energy", "zpe", "s2", "charge", "mult",
        "term", "scf", "geometry_converged", "failure_reason",
        "enthalpy", "inner_energy", "entropy", "nelectrons",
    }
    values: dict = {}
    # frequency quantities read by their own leaf parsers, not the tail
    if {"ir_modes", "imaginary_modes"} & wanted:
        if "ir_modes" in wanted:
            from delfin.ir_spectrum import _parse_ir_spectrum_text
            values["ir_modes"] = _parse_ir_spectrum_text(out_path) or None
        if "imaginary_modes" in wanted:
            from delfin.imag import _collect_imaginary_modes_text
            values["imaginary_modes"] = \
                _collect_imaginary_modes_text(str(out_path)) or None
        return values
    if "term" in wanted:
        values["term"] = _terminated_normally_text(out_path)
    if "geometry" in wanted:
        # the run's final geometry comes from the ORCA-written xyz next to
        # the output -- the established chaining-reader semantics, not the
        # Cartesian blocks of the .out (those never carried the chaining
        # hand-off). The xyz read is small (a few kB per run).
        values["geometry"] = _final_geometry_text(out_path)
    if {"scf", "final_energy", "free_energy", "zpe", "enthalpy",
        "inner_energy", "entropy", "geometry_converged",
        "failure_reason", "s2"} & wanted:
        if tail is None:
            tail = _tail_text(out_path)
    tail = tail or ""
    if "scf" in wanted:
        values["scf"] = _scf_converged_text(tail)
    if {"final_energy", "free_energy", "zpe", "enthalpy", "inner_energy",
        "entropy", "geometry_converged", "failure_reason", "s2"} & wanted:
        _fill_out_text_values(values, out_path, tail, wanted)
    if {"charge", "mult", "nelectrons"} & wanted:
        header = _coord_header(inp_path)
        values["charge"] = int(header.group(1)) if header else None
        values["mult"] = int(header.group(2)) if header else None
        if "nelectrons" in wanted:
            values["nelectrons"] = _text_nelectrons(inp_path,
                                                    values["charge"])
    return values


def _fill_out_text_values(values: dict, out_path: Path, tail: str,
                          wanted: set) -> None:
    """Output-side quantities, each read from the bounded tail only."""
    from io import StringIO

    from delfin.parser import extract_last_uhf_deviation

    if {"final_energy", "free_energy", "zpe", "enthalpy", "entropy",
        "inner_energy"} & wanted:
        _fill_thermo_and_energy(values, out_path, tail, wanted)
    if "s2" in wanted:
        # extract_last_uhf_deviation reads its argument in full; the
        # argument is the bounded tail, so the read stays bounded.
        from io import StringIO

        try:
            values["s2"] = extract_last_uhf_deviation(StringIO(tail))
        except Exception:
            values["s2"] = None
    if "geometry_converged" in wanted:
        values["geometry_converged"] = _geometry_converged_text(tail)
    if "failure_reason" in wanted:
        values["failure_reason"] = _failure_reason_text(values, out_path)


def _final_energy_text(values: dict, out_path: Path, tail: str) -> None:
    """The final single point energy with the old reader's job selection.

    On a multi-job output (several ``JOB NUMBER`` sections) the bounded
    tail can cut across job boundaries and the last energy in the tail
    belongs to a foreign job (11 of 200 archived runs picked a different
    job than the whole-file reader).  The energy is therefore always read
    from the whole file through :func:`find_electronic_energy`, the
    pre-switch reader, so the job selection is identical by construction;
    ``tail`` is unused here.
    """
    from delfin.energies import find_electronic_energy
    values["final_energy"] = find_electronic_energy(str(out_path))


def _fill_thermo_and_energy(values: dict, out_path: Path, tail: str,
                            wanted: set) -> None:
    """Energy and thermochemistry quantities from the bounded tail."""
    from delfin.energies import gibbs_energy_in_text

    # energies first: the text functions accept pre-read text, so the
    # bounded tail is reused instead of re-reading the file.
    if "final_energy" in wanted:
        _final_energy_text(values, out_path, tail)
    if "free_energy" in wanted:
        values["free_energy"] = gibbs_energy_in_text(tail, str(out_path))
    if "inner_energy" in wanted:
        values["inner_energy"] = _last_float_after(
            tail, r"Total thermal energy\s+(-?\d+\.\d+)\s*Eh")
    if "zpe" in wanted:
        values["zpe"] = _last_float_in_text(tail, _ZPE_PATTERNS, str(out_path))
    if "enthalpy" in wanted:
        values["enthalpy"] = _last_float_after(
            tail, r"Total Enthalpy\s+\.\.\.\s+(-?\d+\.\d+)\s*Eh")
    if "entropy" in wanted:
        values["entropy"] = _last_float_after(
            tail, r"Final entropy term\s+\.\.\.\s+(-?\d+\.\d+)\s*Eh")


_ZPE_PATTERNS = (
    r"Zero[\s-]?point\s+energy.{0,200}?([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)"
    r"\s*(?:E[hH]|a\.u\.)?",
    r"Zero[\s-]?point\s+correction.{0,200}?([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)"
    r"\s*(?:E[hH]|a\.u\.)?",
)


def _last_float_in_text(text: str, patterns, path: str) -> Optional[float]:
    """delfin.energies._last_float_in_text on a pre-read bounded text (the
    original takes the same arguments; importing it here keeps one copy)."""
    from delfin.energies import _last_float_in_text as energies_last_float

    return energies_last_float(text, patterns, path)


_GEOMETRY_DONE_PATTERNS = (
    "OPTIMIZATION RUN DONE",
    "THE OPTIMIZATION HAS CONVERGED",
)
_GEOMETRY_FAILED_PATTERNS = (
    "OPTIMIZATION DID NOT CONVERGE",
    "FAILED TO CONVERGE THE GEOMETRY OPTIMIZATION",
)


def _geometry_converged_text(text: str) -> Optional[bool]:
    """Optimization convergence from a bounded tail (delfin.api
    parse_orca_output semantics: True/False markers, None when absent)."""
    if any(p in text for p in _GEOMETRY_DONE_PATTERNS):
        return True
    if any(p in text for p in _GEOMETRY_FAILED_PATTERNS):
        return False
    return None


def _failure_reason_text(values: dict, out_path: Path) -> Optional[str]:
    """The folder outcome phrase, or None. Same value as before: the old
    text path took it from parse_orca_output's ``outcome`` field, which is
    ``outcome_of_folder``'s phrase; the folder scan reads only exit-code
    and log side files plus a bounded output tail, never the whole .out."""
    if values.get("scf"):
        return None
    from delfin.doc_server.calc_indexer import outcome_of_folder

    return outcome_of_folder(out_path.parent) or None


def _final_geometry_text(out_path: Path) -> Optional[list]:
    """The run's final geometry from the ORCA-written ``<base>.xyz``.

    Exact copy of the chaining readers' established semantics
    (delfin.imag._read_xyz: count line, blank, then ``el x y z`` rows,
    floats; a malformed count line raises). Returns ``None`` when the xyz
    is missing; None is a legitimate answer for a run without one.
    """
    xyz_path = out_path.with_suffix(".xyz")
    if not xyz_path.is_file():
        return None
    lines = xyz_path.read_text(encoding="utf-8",
                               errors="ignore").splitlines()
    count = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + count]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]),
                      float(parts[3])))
    return atoms


def _coord_header(inp_path: Path):
    """``* xyz <charge> <mult>`` from an ORCA input; None when absent."""
    try:
        text = inp_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    return re.search(r"\*\s*(?:xyz|xyzfile)\s+(-?\d+)\s+(\d+)", text)


def _text_inner_energy(out_path: Path) -> Optional[float]:
    """U = Total thermal energy (electronic + ZPE + thermal corrections), Eh."""
    return _last_float_after(_tail_text(out_path),
                             r"Total thermal energy\s+(-?\d+\.\d+)\s*Eh")


def _last_float_after(text: str, pattern: str) -> Optional[float]:
    matches = re.findall(pattern, text, flags=re.IGNORECASE)
    if not matches:
        return None
    try:
        return float(matches[-1])
    except ValueError:
        return None


_ATOMIC_NUMBERS = {
    "H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "NE": 10, "NA": 11, "MG": 12, "AL": 13, "SI": 14, "P": 15,
    "S": 16, "CL": 17, "AR": 18, "K": 19, "CA": 20, "SC": 21, "TI": 22,
    "V": 23, "CR": 24, "MN": 25, "FE": 26, "CO": 27, "NI": 28, "CU": 29,
    "ZN": 30, "GA": 31, "GE": 32, "AS": 33, "SE": 34, "BR": 35, "KR": 36,
    "RB": 37, "SR": 38, "Y": 39, "ZR": 40, "NB": 41, "MO": 42, "TC": 43,
    "RU": 44, "RH": 45, "PD": 46, "AG": 47, "CD": 48, "IN": 49, "SN": 50,
    "SB": 51, "TE": 52, "I": 53, "XE": 54, "CS": 55, "BA": 56, "LA": 57,
    "CE": 58, "PR": 59, "ND": 60, "PM": 61, "SM": 62, "EU": 63, "GD": 64,
    "TB": 65, "DY": 66, "HO": 67, "ER": 68, "TM": 69, "YB": 70, "LU": 71,
    "HF": 72, "TA": 73, "W": 74, "RE": 75, "OS": 76, "IR": 77, "PT": 78,
    "AU": 79, "HG": 80, "TL": 81, "PB": 82, "BI": 83, "PO": 84, "AT": 85,
    "RN": 86, "FR": 87, "RA": 88, "AC": 89, "TH": 90, "PA": 91, "U": 92,
    "NP": 93, "PU": 94, "AM": 95, "CM": 96, "BK": 97, "CF": 98, "ES": 99,
    "FM": 100, "MD": 101, "NO": 102, "LR": 103,
}


def _text_nelectrons(inp_path: Path, charge: Optional[int]) -> Optional[int]:
    """Total electrons = sum of atomic numbers in the input minus charge.

    This is the physical electron count. With an ECP the electrons ORCA
    treats explicitly are fewer (the core electrons are replaced), so this
    value is not ORCA's "Number of Electrons" for such runs.

    Same derivation as OCCUPIER's ``calculate_total_electrons`` minus its
    CONTROL indirection: the run's own ``.inp`` carries the geometry.
    """
    if charge is None:
        return None
    try:
        text = inp_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    m = re.search(r"\*\s*xyz\s+-?\d+\s+\d+\s*\n(.*?)(?:\n\s*\*\s*|\Z)",
                  text, re.DOTALL)
    if not m:
        return None
    total = 0
    for line in m.group(1).splitlines():
        parts = line.split()
        if not parts or parts[0].startswith(("*", "$", "%")):
            continue
        z = _ATOMIC_NUMBERS.get(parts[0].upper())
        if z is None:
            return None
        total += z
    return total - charge if total else None


# ---------------------------------------------------------------------------
# the reader
# ---------------------------------------------------------------------------

#: quantity -> (OPI getter key, text key, unit)
_ENERGY_QUANTITIES = {
    Quantity.FINAL_SINGLE_POINT_ENERGY: ("final_energy", "final_energy", "Eh"),
    Quantity.GIBBS: ("free_energy", "free_energy", "Eh"),
    Quantity.ZPE: ("zpe", "zpe", "Eh"),
    Quantity.ENTHALPY: ("enthalpy", "enthalpy", "Eh"),
    Quantity.INNER_ENERGY: ("inner_energy", "inner_energy", "Eh"),
}

#: The text-read key(s) each Quantity needs -- the only work ``read()`` does
#: on the text path. Thermochemistry quantities share the tail scan but each
#: matches only its own pattern.
_TEXT_KEYS = {
    Quantity.TERMINATED_NORMALLY: {"term"},
    Quantity.GEOMETRY_CONVERGED: {"scf", "geometry_converged",
                                  "failure_reason"},
    Quantity.FAILURE_REASON: {"scf", "geometry_converged",
                              "failure_reason"},
    Quantity.S2_DEVIATION: {"s2"},
    Quantity.CHARGE: {"charge", "mult", "nelectrons"},
    Quantity.MULTIPLICITY: {"charge", "mult", "nelectrons"},
    Quantity.NUMBER_OF_ELECTRONS: {"charge", "mult", "nelectrons"},
    Quantity.FINAL_GEOMETRY: {"geometry"},
}
for _q, (_opi_key, _text_key, _unit) in _ENERGY_QUANTITIES.items():
    _TEXT_KEYS[_q] = {_text_key}
# GIBBS/ZPE/ENTHALPY/INNER_ENERGY/ENTROPY sit in the thermochemistry block;
# ENTROPY is not in _ENERGY_QUANTITIES but uses the same tail scan.
_TEXT_KEYS[Quantity.ENTROPY] = {"entropy"}
#: The frequency quantities read by their own leaf parsers (not the tail
#: scan): IR modes and imaginary modes from the .out, vibrational
#: frequencies from the .hess. Their text readers are the established
#: parsers delfin.ir_spectrum.parse_ir_spectrum and
#: delfin.imag.collect_imaginary_modes / hessian_frequencies.
_TEXT_KEYS[Quantity.IR_MODES] = {"ir_modes"}
_TEXT_KEYS[Quantity.IMAGINARY_MODES] = {"imaginary_modes"}
_TEXT_KEYS[Quantity.VIBRATIONAL_FREQUENCIES] = {"vib_freqs"}

#: quantities whose OPI comparison showed identical values AND where OPI
#: answers the same question -- everything else goes to the text parser.
_OPI_FIRST = (set(_ENERGY_QUANTITIES) | {Quantity.S2_DEVIATION}
              | {Quantity.FINAL_GEOMETRY})


def read(work_dir: Union[str, Path], base_name: str, quantity: Quantity) -> Reading:
    """Read one *quantity* from the finished ORCA run in *work_dir*.

    Prefers OPI where the phase-1 comparison measured identical values,
    falls back to DELFIN's established text parsers otherwise (see the
    module docstring for the exact rules). Never raises for a missing or
    partial run: the Reading then has value ``None`` and source ``"none"``.

    Args:
        work_dir: directory holding ``<base_name>.out`` and, for OPI,
            ``<base_name>.property.json``.
        base_name: the run's base name (its ``.out`` without extension).
        quantity: one of :class:`Quantity`.

    Returns:
        Reading: ``value`` (Eh for energies, counts for charge/multipli-
        city/electrons, bool for terminated/converged, str for the failure
        reason), ``unit`` ("" for dimensionless), ``source`` ("opi",
        "text" or "none") and a ``note`` when a fallback or mismatch
        happened.
    """
    work_dir = Path(work_dir)
    out_path = work_dir / f"{base_name}.out"
    inp_path = work_dir / f"{base_name}.inp"

    def none(note: str) -> Reading:
        return Reading(None, "", "none", note)

    if not out_path.is_file():
        return none(f"no {base_name}.out in the run directory")

    text = ""

    # quantities that stay on the text reader by design
    if quantity not in _OPI_FIRST:
        return _read_text_only(work_dir, base_name, quantity, inp_path)

    text = _tail_text(out_path)
    text_values = _text_values(out_path, inp_path,
                               _TEXT_KEYS.get(quantity), tail=text)
    if quantity is Quantity.FINAL_GEOMETRY:
        text_value = text_values.get("geometry")
    else:
        text_value = text_values.get(
            quantity == Quantity.S2_DEVIATION and "s2"
            or _ENERGY_QUANTITIES[quantity][1])

    if _reader_mode() == "text":
        return _reading_from_text(quantity, text_value, text_values,
                                  note="orca.reader=text")

    usable, reason = _opi_ready(work_dir, base_name, out_path, text)
    if not usable:
        return _reading_from_text(quantity, text_value, text_values, note=reason)

    opi_values = _opi_values(work_dir, base_name)
    if "error" in opi_values:
        return _reading_from_text(quantity, text_value, text_values,
                                  note=f"OPI parse failed ({opi_values['error']})")

    if quantity is Quantity.S2_DEVIATION:
        opi_value = _opi_s2_deviation(opi_values.get("s2"))
        if opi_value is None:
            # closed shell: OPI has no UHF block; the deviation is 0.0
            opi_value = 0.0
    elif quantity is Quantity.FINAL_GEOMETRY:
        opi_value = opi_values.get("geometry")
    else:
        opi_value = opi_values.get(_ENERGY_QUANTITIES[quantity][0])

    if opi_value is None:
        return _reading_from_text(quantity, text_value, text_values,
                                  note="OPI carries no value for this quantity")

    # shadow check: both readers answered; the text value is the record
    if text_value is None:
        if quantity is Quantity.FINAL_GEOMETRY:
            return Reading(opi_value, "Angstrom", "opi",
                           "text parser has no value; OPI value unshadowed")
        return Reading(opi_value, _ENERGY_QUANTITIES.get(quantity, (None, None, "Eh"))[2]
                       if quantity in _ENERGY_QUANTITIES else "", "opi",
                       "text parser has no value; OPI value unshadowed")
    if quantity is Quantity.FINAL_GEOMETRY:
        if _geometries_agree(text_value, opi_value):
            return Reading(text_value, "Angstrom", "opi",
                           "confirmed by text shadow check")
        logger.warning(
            "orca_read mismatch for %s on run '%s': OPI geometry differs "
            "from the output text xyz -- returning the text value",
            quantity.name, base_name)
        return Reading(text_value, "Angstrom", "text",
                       "OPI disagreed; text value returned")
    if _values_agree(text_value, opi_value):
        return Reading(text_value,
                       _ENERGY_QUANTITIES[quantity][2]
                       if quantity in _ENERGY_QUANTITIES else "",
                       "opi", "confirmed by text shadow check")
    logger.warning(
        "orca_read mismatch for %s on run '%s': OPI says %r, the output text "
        "%r -- returning the text value", quantity.name, base_name,
        opi_value, text_value)
    return Reading(text_value,
                   _ENERGY_QUANTITIES[quantity][2]
                   if quantity in _ENERGY_QUANTITIES else "",
                   "text", "OPI disagreed; text value returned")


def _values_agree(a: float, b: float) -> bool:
    return abs(float(a) - float(b)) <= SHADOW_TOL


def _geometries_agree(text_geom, opi_geom) -> bool:
    """The xyz-derived and the OPI-derived geometry are the same structure.

    The xyz prints 14 decimals in Angstrom; OPI converts the property
    JSON's Bohr doubles with AU_TO_ANGST, so the two agree to ~1e-12 A
    with any rounding of the unit conversion; 1e-6 A leaves that margin
    while still telling structures apart (any real geometry change moves
    atoms by far more).
    """
    if not text_geom or not opi_geom or len(text_geom) != len(opi_geom):
        return False
    for (te, tx, ty, tz), (oe, ox, oy, oz) in zip(text_geom, opi_geom):
        if str(te).capitalize() != str(oe).capitalize():
            return False
        if any(abs(float(a) - float(b)) > GEOMETRY_SHADOW_TOL
               for a, b in ((tx, ox), (ty, oy), (tz, oz))):
            return False
    return True


def _opi_s2_deviation(s2) -> Optional[float]:
    """OPI's (expected, ideal) tuple as DELFIN's deviation value."""
    if not s2:
        return None
    expected, ideal = s2
    return float(expected) - float(ideal)


def _reading_from_text(quantity: Quantity, text_value, text_values: dict,
                       note: str) -> Reading:
    """Build the answer from DELFIN's text readers (the established one)."""
    if quantity is Quantity.FINAL_GEOMETRY:
        if text_value is None:
            return Reading(None, "", "none",
                           f"text parser found nothing ({note})")
        return Reading(text_value, "Angstrom", "text", note)
    if quantity is Quantity.IR_MODES:
        value = text_values.get("ir_modes") if text_value is None else text_value
        if value is None:
            return Reading(None, "", "none",
                           f"text parser found nothing ({note})")
        return Reading(value, "cm-1/km-mol", "text", note)
    if quantity is Quantity.IMAGINARY_MODES:
        value = text_values.get("imaginary_modes") if text_value is None else text_value
        if value is None:
            return Reading(None, "", "none",
                           f"text parser found nothing ({note})")
        return Reading(value, "cm-1", "text", note)
    if quantity in _ENERGY_QUANTITIES or quantity is Quantity.S2_DEVIATION:
        if text_value is None:
            return Reading(None, "", "none", f"text parser found nothing ({note})")
        unit = _ENERGY_QUANTITIES[quantity][2] if quantity in _ENERGY_QUANTITIES else ""
        return Reading(text_value, unit, "text", note)
    if quantity is Quantity.TERMINATED_NORMALLY:
        return Reading(bool(text_values.get("term")), "", "text", note)
    if quantity is Quantity.GEOMETRY_CONVERGED:
        value = text_values.get("geometry_converged")
        return Reading(value, "", "none" if value is None else "text", note)
    if quantity is Quantity.FAILURE_REASON:
        value = text_values.get("failure_reason")
        return Reading(value, "", "none" if value is None else "text", note)
    if quantity is Quantity.MULTIPLICITY:
        value = text_values.get("mult")
        return Reading(value, "", "none" if value is None else "text", note)
    if quantity is Quantity.CHARGE:
        value = text_values.get("charge")
        return Reading(value, "", "none" if value is None else "text", note)
    if quantity is Quantity.NUMBER_OF_ELECTRONS:
        value = text_values.get("nelectrons")
        return Reading(value, "", "none" if value is None else "text", note)
    if quantity is Quantity.ENTROPY:
        value = text_values.get("entropy")
        if value is None:
            return Reading(None, "", "none",
                           f"text parser found nothing ({note})")
        return Reading(value, "Eh (T*S)", "text",
                       "ORCA prints T*S, not the molar entropy; naming "
                       "decision pending" + (f" ({note})" if note else ""))
    return Reading(None, "", "none", f"unsupported quantity {quantity.name}")


def _read_text_only(work_dir: Path, base_name: str, quantity: Quantity,
                    inp_path: Path) -> Reading:
    """Quantities that never route through OPI (see module docstring)."""
    out_path = work_dir / f"{base_name}.out"
    if quantity is Quantity.VIBRATIONAL_FREQUENCIES:
        return _reading_vib_freqs(work_dir, base_name, note="")
    text_values = _text_values(out_path, inp_path, _TEXT_KEYS.get(quantity))
    return _reading_from_text(quantity, None, text_values, note="")


def _reading_vib_freqs(work_dir: Path, base_name: str, note: str) -> Reading:
    """$vibrational_frequencies of <base>.hess via imag.hessian_frequencies."""
    from delfin.imag import _hessian_frequencies_text

    hess_path = work_dir / f"{base_name}.hess"
    if not hess_path.is_file():
        return Reading(None, "", "none", f"no {base_name}.hess in the run dir"
                       + (f" ({note})" if note else ""))
    freqs = _hessian_frequencies_text(hess_path)
    if not freqs:
        return Reading(None, "", "none", f"no $vibrational_frequencies in "
                       f"{base_name}.hess" + (f" ({note})" if note else ""))
    return Reading(freqs, "cm-1", "text", note)


def read_hess(hess_path: Union[str, Path], quantity: Quantity) -> Reading:
    """Read one *quantity* from a single .hess file (no run directory).

    Only vibrational-frequencies quantities are defined on a lone .hess;
    anything else is answered with a ``none`` Reading. Never raises for
    a missing file. The text reader is delfin.imag.hessian_frequencies
    -- the established parser; the .hess carries no OPI counterpart.
    """
    if quantity is not Quantity.VIBRATIONAL_FREQUENCIES:
        return Reading(None, "", "none",
                       f"{quantity.name} is not read from a .hess file")
    hess_path = Path(hess_path)
    if not hess_path.is_file():
        return Reading(None, "", "none", f"{hess_path.name} does not exist")
    from delfin.imag import _hessian_frequencies_text
    freqs = _hessian_frequencies_text(hess_path)
    if not freqs:
        return Reading(None, "", "none",
                       f"no $vibrational_frequencies in {hess_path.name}")
    return Reading(freqs, "cm-1", "text", "")


# ---------------------------------------------------------------------------
# optional helper: create OPI's property JSON in a DELFIN work dir
# ---------------------------------------------------------------------------

def ensure_property_json(work_dir: Union[str, Path], base_name: str) -> Optional[Path]:
    """Create ``<base>.property.json`` with ``orca_2json``, when allowed.

    Only inside the run's own directory (a DELFIN work dir, never a path
    outside the run), and only when ``orca_2json`` is on PATH. Returns
    the created file's path, or ``None`` when it was not created
    (already present, no ``orca_2json``, or no ``.property.txt`` to
    convert). Refuses to run anywhere but *work_dir* itself.
    """
    work_dir = Path(work_dir).resolve()
    target = work_dir / f"{base_name}.property.json"
    if target.is_file():
        return target
    if not (work_dir / f"{base_name}.property.txt").is_file():
        return None
    if shutil.which("orca_2json") is None:
        logger.info("orca_2json not on PATH; property JSON not created for '%s'",
                    base_name)
        return None
    from opi.execution.core import Runner

    runner = Runner(working_dir=work_dir)
    runner.create_property_json(base_name)
    return target if target.is_file() else None


# ---------------------------------------------------------------------------
# relativistic / spin-property readers (ZFS via %eprnmr, SOC-TDDFT)
# ---------------------------------------------------------------------------

ZFS_CM1 = "cm-1"


@dataclass
class ZfsParameters:
    """Zero-field splitting read from ORCA's ``D_Tensor`` property record.

    ORCA writes the axial and rhombic ZFS parameters into the run's
    ``.property.json`` under a ``D_Tensor`` record whose ``D`` (axial,
    cm-1) and ``E`` (rhombic, cm-1) fields are documented in the manual
    (property-file documentation, ``D_Tensor`` record).  ``E_over_D`` is
    dimensionless and derived as ``E / D`` when ``D != 0``, else ``None``.
    """

    d_cm1: Optional[float] = None      # axial ZFS parameter D
    e_cm1: Optional[float] = None      # rhombic ZFS parameter E
    e_over_d: Optional[float] = None   # derived, dimensionless (E/D)
    source: str = "none"               # "opi" | "json" | "none"

    @property
    def unit(self) -> tuple[str, str, str]:
        return (ZFS_CM1, ZFS_CM1, "")


def _json_find_d_tensor(node: Any, depth: int = 0) -> Optional[dict]:
    """Return the first nested dict that carries ORCA's ``D_Tensor`` record.

    ORCA's property JSON nests each property record as a flat dict of the
    documented keys (the ``D_Tensor`` record's ``D``/``E``/``d_raw``/
    ``d_eigenvalues``/``d_eigenvectors`` plus ``method``/``level``/...
    -- cf. the manual's property-file documentation).  The nesting depends
    on the ORCA version, so we locate the record by its keys instead of
    guessing a fixed path: a dict that carries an axial ``D`` and rhombic
    ``E`` scalar together with the ``d_raw`` tensor must *be* the record.
    """
    if depth > 6:
        return None
    if isinstance(node, dict):
        lo = {k.lower(): v for k, v in node.items()}
        if "d_tensor" in lo and isinstance(lo["d_tensor"], dict):
            return lo["d_tensor"]
        if ("d" in lo and "e" in lo
                and _coerce_float(lo["d"]) is not None
                and _coerce_float(lo["e"]) is not None):
            return node
        for value in node.values():
            found = _json_find_d_tensor(value, depth + 1)
            if found is not None:
                return found
    elif isinstance(node, list):
        for value in node:
            found = _json_find_d_tensor(value, depth + 1)
            if found is not None:
                return found
    return None


def _json_loads_d_tensor(json_text: str) -> Optional[ZfsParameters]:
    """ZFS D/E from a *D_Tensor* property-JSON record given as text.

    ``json_text`` is ORCA's ``<base>.property.json`` content (synthetic in
    tests).  Returns ``None`` when no ``D_Tensor`` record is found.  Never
    raises on malformed JSON.
    """
    try:
        data = json.loads(json_text)
    except (ValueError, TypeError):
        return None
    node = _json_find_d_tensor(data)
    if node is None:
        return None
    d_f = _coerce_float(node.get("D"))
    e_f = _coerce_float(node.get("E"))
    if d_f is None and e_f is None:
        return None
    return ZfsParameters(
        d_cm1=d_f,
        e_cm1=e_f,
        e_over_d=(e_f / d_f) if (d_f and d_f != 0.0) else None,
        source="json",
    )


def _coerce_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_zfs_parameters(work_dir: Union[str, Path],
                        base_name: str) -> Optional[ZfsParameters]:
    """ZFS D, E (cm-1) and E/D of a finished ``%eprnmr DTensor`` run.

    ORCA holds the ZFS parameters in the run's ``.property.json`` under a
    ``D_Tensor`` record (documented fields ``D`` axial, ``E`` rhombic);
    that record is the authoritative source because the manual documents
    no stable plain-text ``.out`` table for the two scalars.  Returns
    ``None`` when the property JSON is absent, unreadable or carries no
    ``D_Tensor`` record.  Never raises.
    """
    json_path = Path(work_dir) / f"{base_name}.property.json"
    if not json_path.is_file():
        return None
    try:
        text = json_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    return _json_loads_d_tensor(text)


# section markers: the SOC-TDDFT spectrum mirrors ROCIS (manual: tddft.txt
# "Spin-orbit coupling" -- "similarly to what is done in ROCIS"; rocis.txt
# documents the .out table verbatim).
_SOC_SPECTRUM_RE = re.compile(
    r"SOC\s+CORRECTED\s+ABSORPTION\s+SPECTRUM\s+VIA\s+TRANSITION"
    r"\s+ELECTRIC\s+DIPOLE\s+MOMENTS",
    re.IGNORECASE,
)
_SOC_STATE_RE = re.compile(
    r"^\s*State\s+(\d+):\s+([-+]?\d+\.\d+)\s+cm\*\*-1\s+"
    r"([-+]?\d+\.\d+)\s+eV",
    re.IGNORECASE,
)


@dataclass
class SocState:
    """One SOC-corrected TD-DFT excited state (transition)."""

    states: int                     # the state list index (0-based)
    e_cm1: float                    # transition energy, cm**-1
    e_eV: float                     # transition energy, eV

    @property
    def unit(self) -> tuple[str, str]:
        return (ZFS_CM1 + "(cm**-1)", "eV")


def _parse_soc_states(text: str) -> list[SocState]:
    """SOC transitions from the SOC-TDDFT ``.out`` text section."""
    states: list[SocState] = []
    started = False
    for line in text.splitlines():
        if not started:
            if _SOC_SPECTRUM_RE.search(line):
                started = True
            continue
        m = _SOC_STATE_RE.match(line)
        if m:
            states.append(SocState(int(m.group(1)), float(m.group(2)),
                                   float(m.group(3))))
        elif states and line.strip() and "State" not in line:
            break  # section ended
    return states


def read_soc_states(work_dir: Union[str, Path],
                    base_name: str) -> list[SocState]:
    """SOC-corrected TD-DFT excited states read from the ``.out`` text.

    Anchors on the ``SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION
    ELECTRIC DIPOLE MOMENTS`` section (mirrors ROCIS, manual rocis.txt;
    SOC-TDDFT reuses it per tddft.txt) and reads ``State N:  <E> cm**-1
    <E> eV`` rows from it.  Returns ``[]`` when the section is absent.
    Never raises.
    """
    out_path = Path(work_dir) / f"{base_name}.out"
    if not out_path.is_file():
        return []
    text = _tail_text(out_path)
    return _parse_soc_states(text)


# ---------------------------------------------------------------------------
# thermochemistry reading from an existing frequency .out (text only)
# ---------------------------------------------------------------------------

#: 1 Eh in kcal/mol (CODATA 2018; DELFIN's standard constant, see pka.py).
_HARTREE_TO_KCAL = 627.5094740631

#: The temperature is carried in the block header that ORCA prints over each
#: thermochemistry section (manual, section 4.7). A ``%freq Temp`` list run
#: (manual, section 4.6.1) prints one such block per temperature.
_THERMO_BLOCK_RE = re.compile(r"THERMOCHEMISTRY AT\s+(-?\d+(?:\.\d+)?)\s*K\b",
                              re.IGNORECASE)
_THERMO_QUANTITY_RE = re.compile(
    r"{label}\s+\.\.\.\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*{unit}",
    re.IGNORECASE)
#: ``<label> ... <eh> Eh  <kcal> kcal/mol`` -- ORCA prints ZPE and the
#: entropy term in both units on one line (manual, section 4.7).
_THERMO_BOTH_RE = re.compile(
    r"{label}\s+\.\.\.\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*Eh\s+"
    r"(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*kcal/mol",
    re.IGNORECASE)
#: The line ORCA prints when the entropy is evaluated with the quasi-RRHO
#: treatment of S. Grimme (manual, section 4.7). When this exact line is
#: absent but a vibrational-entropy line is present, plain RRHO is used.
_QUASI_RRHO_LINE = "QRRHO of S. Grimme"


@dataclass(frozen=True)
class ThermoReading:
    """The thermochemistry of one temperature block of a frequency .out.

    All energy quantities carry both hartree and kcal/mol. ORCA prints the
    electronic energy, the total enthalpy and the Gibbs free energy only in
    hartree (manual, section 4.7); their kcal/mol values are computed here
    with :data:`_HARTREE_TO_KCAL`. ``zpe`` and ``ts`` are the values ORCA
    itself prints in both units.

    ``quasirrho`` is ``True`` when ORCA used the quasi-RRHO entropy
    treatment (manual 4.7, ``%freq QuasiRRHO`` default true), ``False`` for
    plain RRHO, ``None`` when no vibrational-entropy line was found.
    """

    #: temperature of this block, K
    temperature_K: Optional[float]
    #: pressure of this block, atm
    pressure_atm: Optional[float]
    #: electronic (single-point) energy, Eh and kcal/mol
    e_el_Eh: Optional[float]
    e_el_kcal: Optional[float]
    #: zero-point energy correction, Eh and kcal/mol
    zpe_Eh: Optional[float]
    zpe_kcal: Optional[float]
    #: total enthalpy H, Eh and kcal/mol
    h_Eh: Optional[float]
    h_kcal: Optional[float]
    #: entropy term as ORCA prints it: T*S, Eh and kcal/mol
    ts_Eh: Optional[float]
    ts_kcal: Optional[float]
    #: Gibbs free energy G, Eh and kcal/mol
    g_Eh: Optional[float]
    g_kcal: Optional[float]
    #: True for quasi-RRHO, False for RRHO, None when undetermined
    quasirrho: Optional[bool]
    #: reading source ("text" always here)
    source: str = "text"
    #: human note (empty on success)
    note: str = ""

    @property
    def e_el(self) -> tuple[Optional[float], Optional[float]]:
        """(Eh, kcal/mol) electronic energy."""
        return (self.e_el_Eh, self.e_el_kcal)

    @property
    def zpe(self) -> tuple[Optional[float], Optional[float]]:
        """(Eh, kcal/mol) zero-point energy."""
        return (self.zpe_Eh, self.zpe_kcal)

    @property
    def enthalpy(self) -> tuple[Optional[float], Optional[float]]:
        """(Eh, kcal/mol) total enthalpy H."""
        return (self.h_Eh, self.h_kcal)

    @property
    def ts(self) -> tuple[Optional[float], Optional[float]]:
        """(Eh, kcal/mol) entropy term T*S."""
        return (self.ts_Eh, self.ts_kcal)

    @property
    def gibbs(self) -> tuple[Optional[float], Optional[float]]:
        """(Eh, kcal/mol) Gibbs free energy G."""
        return (self.g_Eh, self.g_kcal)


def _thermo_quantity(block: str, label: str, unit: str) -> Optional[float]:
    """First ``<label> ... <value> <unit>`` on one line within *block*."""
    pattern = _THERMO_QUANTITY_RE.pattern.format(label=re.escape(label),
                                                 unit=re.escape(unit))
    m = re.search(pattern, block)
    if m is None:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None


def _thermo_both(block: str, label: str) -> tuple[Optional[float], Optional[float]]:
    """``(Eh, kcal/mol)`` of the ``<label> ... <eh> Eh <kcal> kcal/mol`` line."""
    pattern = _THERMO_BOTH_RE.pattern.format(label=re.escape(label))
    m = re.search(pattern, block)
    if m is None:
        return (None, None)
    try:
        return (float(m.group(1)), float(m.group(2)))
    except ValueError:
        return (None, None)


def read_thermochemistry(out_path: Union[str, Path]) -> list[ThermoReading]:
    """Read every thermochemistry block of a frequency run's ``.out``.

    ORCA prints one ``THERMOCHEMISTRY AT <T>K`` block per temperature
    (manual, section 4.7). A run requested at several temperatures with
    ``%freq Temp 290,295,...`` (section 4.6.1) carries one block per
    temperature; this returns one :class:`ThermoReading` per block, in file
    order. This is a pure text reader and never needs OPI.

    A missing file (or a file with no thermochemistry block) yields an empty
    list. It never raises for a missing or unparseable output.
    """
    try:
        text = Path(out_path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    starts = [(m.start(), float(m.group(1)))
              for m in _THERMO_BLOCK_RE.finditer(text)]
    readings: list[ThermoReading] = []
    for i, (start, temp) in enumerate(starts):
        end = starts[i + 1][0] if i + 1 < len(starts) else len(text)
        block = text[start:end]

        e_el = _thermo_quantity(block, "Electronic energy", "Eh")
        zpe_eh, zpe_kcal = _thermo_both(block, "Zero point energy")
        h = _thermo_quantity(block, "Total Enthalpy", "Eh")
        ts_eh, ts_kcal = _thermo_both(block, "Final entropy term")
        if ts_eh is None:
            # ORCA also prints the entropy as "Total entropy correction"
            # (manual 4.7), but with the opposite sign: -T*S, the term already
            # applied in G = H - T*S, whereas "Final entropy term" prints +T*S
            # (verified: h2o_optfreq real output lines 3437 vs 3469). Negate so
            # the reported entropy term stays +T*S, sign-stable whichever of
            # the two labels survived a truncated-but-terminated output.
            corr_eh, corr_kcal = _thermo_both(block, "Total entropy correction")
            if corr_eh is not None:
                ts_eh, ts_kcal = -corr_eh, -corr_kcal
        g = _thermo_quantity(block, "Final Gibbs free energy", "Eh")
        pressure = _thermo_quantity(block, "Pressure", "atm")

        if "QRRHO of S. Grimme" in block:
            quasirrho = True
        elif "RRHO" in block:
            quasirrho = False
        else:
            quasirrho = None

        def _kcal(eh: Optional[float]) -> Optional[float]:
            return None if eh is None else eh * _HARTREE_TO_KCAL

        readings.append(ThermoReading(
            temperature_K=temp,
            pressure_atm=pressure,
            e_el_Eh=e_el, e_el_kcal=_kcal(e_el),
            zpe_Eh=zpe_eh, zpe_kcal=zpe_kcal,
            h_Eh=h, h_kcal=_kcal(h),
            ts_Eh=ts_eh, ts_kcal=ts_kcal,
            g_Eh=g, g_kcal=_kcal(g),
            quasirrho=quasirrho,
        ))
    return readings
