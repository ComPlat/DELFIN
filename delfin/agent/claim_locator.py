"""Locate and verify a reported number+unit in the real output file of a run.

VR1 (verify-before-report): a number with a unit in an answer must be traced to
the real output FILE of the run -- file path, line number, parsed value and
unit -- converted to a common unit and compared with an explicit per-quantity
tolerance, before the number is reported.

Division of labour, stated honestly:

* Energy claims reuse DELFIN's own parsers (``delfin.energies``:
  ``find_electronic_energy`` / ``find_gibbs_energy`` / ``find_ZPE``) for the
  run's reported values in Eh.  The value's ``file:line`` is found by locating
  the labelled energy line in the file text.
* Frequency, absolute-coordinate (distance) and dipole values have no in-repo
  value parser in ``delfin/`` for a ``.out`` file (DELFIN's full-mode, geometry
  and dipole extraction lives in the ops server's typed tools, not in this
  repo).  For them the claim value is located directly on its output line.
  That is locating a claimed number, not re-parsing ORCA semantics.

Verdicts
  * ``verified``  -- a value on a known output line matches the claim (value
    converted to the quantity's canonical unit, within that quantity's
    tolerance).
  * ``mismatch``  -- the output computes the quantity, but the claimed value
    matches none of the output's values within tolerance.
  * ``no-source`` -- the output does not report the quantity at all.
  * ``estimate``  -- the unit maps to no physical quantity (a count, a year, a
    setting); such numbers are marked, never flagged as wrong.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, asdict
from pathlib import Path

from delfin import energies

# ---------------------------------------------------------------------------
# Unit registry
# ---------------------------------------------------------------------------

# Conversion factors, claim unit -> canonical unit of the quantity.
# 1 Hartree (Eh) = 627.509 kcal/mol = 2625.500 kJ/mol = 27.2114 eV.
_EHR = {
    "Eh": 1.0,
    "Hartree": 1.0,
    "Ha": 1.0,
    "eV": 1.0 / 27.2114,
    "meV": 1.0 / 27211.4,
    "kcal/mol": 1.0 / 627.509,
    "kcal mol-1": 1.0 / 627.509,
    "kJ/mol": 1.0 / 2625.500,
    "kJ mol-1": 1.0 / 2625.500,
}
_FREQ = {"cm-1": 1.0, "cm**-1": 1.0, "cm^-1": 1.0}
_DIST = {"Angstrom": 1.0, "A": 1.0, "bohr": 0.529177}
_DIPOLE = {"Debye": 1.0, "D": 1.0}

# quantity -> (canonical unit, tolerance in canonical unit)
QUANTITY_TOLERANCE = {
    "energy": ("Eh", 1.0 / 627.509),       # 1.0 kcal/mol, in Eh
    "frequency": ("cm-1", 1.0),            # 1.0 cm-1
    "distance": ("Angstrom", 0.01),        # 0.01 Angstrom
    "dipole": ("Debye", 0.1),              # 0.1 Debye
}

_UNIT_QUANTITY = {}
for _factor_map, _q in (
    (_EHR, "energy"),
    (_FREQ, "frequency"),
    (_DIST, "distance"),
    (_DIPOLE, "dipole"),
):
    for _u, _f in _factor_map.items():
        _UNIT_QUANTITY[_u] = (_q, _f)

_FLOAT_RE = r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?"
_FLOAT = re.compile(_FLOAT_RE)

# A prose number followed (optional single space) by a known unit token.
# Units are tried longest-first so "kcal/mol" wins over "mol" and "Angstrom"
# over "A".  The trailing look-ahead rejects a unit that is a prefix of a
# longer token ("eV" inside "eVs", "A" inside "Angstrom" as the second unit).
_UNIT_KEYS = sorted(_UNIT_QUANTITY, key=len, reverse=True)
_CLAIM_TOKEN = re.compile(
    r"(" + _FLOAT_RE + r")\s*("
    + "|".join(re.escape(u) for u in _UNIT_KEYS)
    + r")(?![A-Za-z0-9/])"
)

_ANCHORS = {
    # quantity -> head-of-block or labelled-line markers, in file order.
    "frequency": re.compile(r"VIBRATIONAL FREQUENCIES"),
    "distance": re.compile(r"CARTESIAN COORDINATES"),
    "dipole": re.compile(r"DIPOLE|Dipole"),
}
_ENERGY_LABELS = {
    "FINAL SINGLE POINT ENERGY": energies.find_electronic_energy,
    "Final Gibbs free energy": energies.find_gibbs_energy,
    "Zero point energy": energies.find_ZPE,
}


@dataclass
class VerificationRecord:
    """One claim's verification result, fully machine-readable."""

    verdict: str
    value: float
    unit: str
    quantity: str
    file: str | None = None
    line: int | None = None
    parsed_value: float | None = None
    parsed_unit: str | None = None
    tolerance: float | None = None

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _largest_out(calc_folder: Path) -> Path | None:
    """Largest ``*.out`` in the folder (same convention as result_critic)."""
    outs = sorted(calc_folder.glob("*.out"), key=lambda p: p.stat().st_size)
    return outs[-1] if outs else None


def _first_float(text: str) -> float | None:
    m = _FLOAT.search(text)
    return float(m.group()) if m else None


def _label_line(lines: list[str], label: str) -> int:
    """1-based line number of the line containing ``label`` and a float."""
    for i, line in enumerate(lines, start=1):
        if label in line and _FLOAT.search(line):
            return i
    return 0


def _energy_candidates(
    out: Path,
) -> list[tuple[float, str, int, str]]:
    """(value_in_Eh, label, line, canonical_unit) for each energy DELFIN reads."""
    lines = out.read_text(errors="replace").splitlines()
    result: list[tuple[float, str, int, str]] = []
    for label, finder in _ENERGY_LABELS.items():
        value = finder(str(out))
        if value is None:
            continue
        line = _label_line(lines, label)
        result.append((value, label, line, "Eh"))
    return result


def _block_values(out: Path, marker: re.Pattern, unit: str) -> list[tuple[float, int]]:
    """(value, line) for every real float token after a block marker's lines.

    A token that is a row counter (a float immediately followed by ``:``, e.g.
    the ``0:`` in ``   0:     3441.3``) is an index, not a value, and is
    skipped; every other float on a data line is a candidate value.
    """
    lines = out.read_text(errors="replace").splitlines()
    values: list[tuple[float, int]] = []
    in_block = False
    for i, line in enumerate(lines, start=1):
        if marker.search(line):
            in_block = True
            continue
        if in_block and line.strip() == "":
            break
        if in_block:
            for m in _FLOAT.finditer(line):
                after = line[m.end():m.end() + 1]
                if after == ":":
                    continue  # row counter, not a value
                values.append((float(m.group()), i))
    return values


_ATOM_LINE = re.compile(
    r"^\s*([A-Z][a-z]?)\s+"
    r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)\s*$"
)


def _geometry_atoms(out: Path) -> list[tuple[str, float, float, float]]:
    """Per-atom ``(element, x, y, z)`` from the CARTESIAN COORDINATES block.

    Coordinates are in Angstrom (the ``(ANGSTROEM)`` header).  A line that is
    not a clean ``ELEM  x  y  z`` atom row is skipped; the block ends at the
    first blank line, the same boundary ``_block_values`` uses.
    """
    marker = _ANCHORS["distance"]
    lines = out.read_text(errors="replace").splitlines()
    atoms: list[tuple[str, float, float, float]] = []
    in_block = False
    for line in lines:
        if marker.search(line):
            in_block = True
            continue
        if in_block and line.strip() == "":
            break
        if not in_block:
            continue
        m = _ATOM_LINE.match(line)
        if m:
            atoms.append(
                (m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4)))
            )
    return atoms


def _interatomic_distances(out: Path) -> list[tuple[float, int]]:
    """All unique pairwise interatomic distances in Angstrom, in block order.

    A distance needs TWO atoms: the candidates are the Euclidean separations
    between every pair of atoms read from the coordinate block.  A single
    atom's coordinate (x, y or z) is a position, never a distance (VR1-1).
    The line reported is the first of the two atom rows.
    """
    atoms = _geometry_atoms(out)
    distances: list[tuple[float, int]] = []
    for a, (_, ax, ay, az) in enumerate(atoms):
        line_i = _geometry_atom_line(out, a)
        for b in range(a + 1, len(atoms)):
            _, bx, by, bz = atoms[b]
            d = math.sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)
            distances.append((d, line_i))
    return distances


def _geometry_atom_line(out: Path, index: int) -> int:
    """1-based line of the ``index``-th atom row in the coordinate block."""
    marker = _ANCHORS["distance"]
    lines = out.read_text(errors="replace").splitlines()
    in_block = False
    seen = 0
    for i, line in enumerate(lines, start=1):
        if marker.search(line):
            in_block = True
            continue
        if in_block and line.strip() == "":
            break
        if in_block and _ATOM_LINE.match(line):
            if seen == index:
                return i
            seen += 1
    return 0


def _candidate_values(out: Path, quantity: str) -> list[tuple[float, int, str]]:
    """Output's candidate (value, line, canonical_unit) for a quantity."""
    if quantity == "energy":
        return [(v, ln, unit) for v, _label, ln, unit in _energy_candidates(out)]
    canonical = QUANTITY_TOLERANCE[quantity][0]
    if quantity == "distance":
        # Interatomic distances are computed from the coordinate block (two
        # atoms each).  A bare coordinate is never a distance (VR1-1).
        return [(d, ln, canonical) for d, ln in _interatomic_distances(out)]
    marker = _ANCHORS[quantity]
    raw = _block_values(out, marker, canonical)
    return [(v, ln, canonical) for v, ln in raw]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def verify_claim(value: float, unit: str, calc_folder: str | Path) -> VerificationRecord:
    """Locate ``value``+``unit`` in the run's real output and classify it.

    ``calc_folder`` is the calculation folder; its largest ``*.out`` is the
    evidence file (the convention result_critic uses).  Returns a
    fully-populated :class:`VerificationRecord`.
    """
    root = Path(calc_folder)
    unit_q = _UNIT_QUANTITY.get(unit)
    if unit_q is None:
        return VerificationRecord("estimate", value, unit, "unknown")
    quantity, to_canonical = unit_q
    canonical, tolerance = QUANTITY_TOLERANCE[quantity]
    claim_canonical = value * to_canonical

    out = _largest_out(root)
    if out is None:
        return VerificationRecord(
            "no-source", value, unit, quantity, file=str(calc_folder)
        )

    candidates = _candidate_values(out, quantity)
    if not candidates:
        return VerificationRecord(
            "no-source",
            value,
            unit,
            quantity,
            file=str(out),
            parsed_unit=canonical,
            tolerance=tolerance,
        )

    best = None
    for cand_value, cand_line, cand_unit in candidates:
        if abs(claim_canonical - cand_value) <= tolerance:
            best = (cand_value, cand_line, cand_unit)
            break
    if best is None:
        return VerificationRecord(
            "mismatch",
            value,
            unit,
            quantity,
            file=str(out),
            parsed_unit=canonical,
            tolerance=tolerance,
        )
    cand_value, cand_line, cand_unit = best
    return VerificationRecord(
        "verified",
        value,
        unit,
        quantity,
        file=str(out),
        line=cand_line,
        parsed_value=cand_value,
        parsed_unit=cand_unit,
        tolerance=tolerance,
    )


def locate_claims(text: str) -> list[VerificationRecord]:
    """Every number+unit token in prose text, as a pre-verification record.

    A ``(value, unit)`` pair whose unit maps to a known physical quantity
    (energy / frequency / distance / dipole) becomes one
    :class:`VerificationRecord` with verdict ``"located"`` (no file/line yet):
    the claim to be verified against the real run output by the caller via
    :func:`verify_claim`.

    Prose numbers without a unit (counts, years, indices, settings) and
    numbers whose unit maps to no physical quantity are NOT claims and are
    skipped.  They can therefore never produce a false ``mismatch`` /
    ``no-source`` alarm -- the "prose that isn't a result" control.
    """
    claims: list[VerificationRecord] = []
    for m in _CLAIM_TOKEN.finditer(text):
        value = float(m.group(1))
        raw_unit = m.group(2)
        map_key = _UNIT_QUANTITY.get(raw_unit)
        if map_key is None:
            continue
        quantity, _factor = map_key
        claims.append(VerificationRecord("located", value, raw_unit, quantity))
    return claims

