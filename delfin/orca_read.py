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

import logging
import re
import shutil
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional, Union

logger = logging.getLogger(__name__)

#: The .out prints <S**2> to 6 and energies to 8 decimals; OPI's JSON
#: carries full doubles. The phase-1 equivalence tests compared at the
#: text's printed precision; the shadow check uses the same.
SHADOW_TOL = 1e-6

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
    if "term" in wanted:
        values["term"] = _terminated_normally_text(out_path)
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


def _fill_thermo_and_energy(values: dict, out_path: Path, tail: str,
                            wanted: set) -> None:
    """Energy and thermochemistry quantities from the bounded tail."""
    from delfin.energies import (electronic_energy_in_text,
                                 gibbs_energy_in_text)

    # energies first: the text functions accept pre-read text, so the
    # bounded tail is reused instead of re-reading the file.
    if "final_energy" in wanted:
        values["final_energy"] = electronic_energy_in_text(tail, str(out_path))
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
}
for _q, (_opi_key, _text_key, _unit) in _ENERGY_QUANTITIES.items():
    _TEXT_KEYS[_q] = {_text_key}
# GIBBS/ZPE/ENTHALPY/INNER_ENERGY/ENTROPY sit in the thermochemistry block;
# ENTROPY is not in _ENERGY_QUANTITIES but uses the same tail scan.
_TEXT_KEYS[Quantity.ENTROPY] = {"entropy"}

#: quantities whose OPI comparison showed identical values AND where OPI
#: answers the same question -- everything else goes to the text parser.
_OPI_FIRST = set(_ENERGY_QUANTITIES) | {Quantity.S2_DEVIATION}


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
    text_value = text_values.get(quantity == Quantity.S2_DEVIATION
                                 and "s2" or _ENERGY_QUANTITIES[quantity][1])

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
    else:
        opi_value = opi_values.get(_ENERGY_QUANTITIES[quantity][0])

    if opi_value is None:
        return _reading_from_text(quantity, text_value, text_values,
                                  note="OPI carries no value for this quantity")

    # shadow check: both readers answered; the text value is the record
    if text_value is None:
        return Reading(opi_value, _ENERGY_QUANTITIES.get(quantity, (None, None, "Eh"))[2]
                       if quantity in _ENERGY_QUANTITIES else "", "opi",
                       "text parser has no value; OPI value unshadowed")
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


def _opi_s2_deviation(s2) -> Optional[float]:
    """OPI's (expected, ideal) tuple as DELFIN's deviation value."""
    if not s2:
        return None
    expected, ideal = s2
    return float(expected) - float(ideal)


def _reading_from_text(quantity: Quantity, text_value, text_values: dict,
                       note: str) -> Reading:
    """Build the answer from DELFIN's text readers (the established one)."""
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
    text_values = _text_values(out_path, inp_path, _TEXT_KEYS.get(quantity))
    return _reading_from_text(quantity, None, text_values, note="")


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
