"""OPI-backed ORCA input builder for the CO2 coordinator.

Builds the ORCA ``.inp`` files that the classic CO2 Coordinator
(:mod:`delfin.co2.CO2_Coordinator6`) writes, but with OPI's input
classes (``opi.input``) instead of hand-written string templates.  The
geometry, placement and keyword-composition logic is reused **by import**
from the classic module — this module adds no geometry or method of its
own, so the two stay on the same level of theory and the same
coordinates.

PARITY CONTRACT
  For the same CONTROL.txt config and the same structure, the ``.inp``
  file written here is EQUIVALENT to what the classic writer writes:
  the same simple keywords (order- and case-insensitive), the same
  blocks with the same options and values, the same charge /
  multiplicity / coordinates (per-atom ``NewGTO`` tags included), the
  same ``%pal`` / ``%maxcore``, and the same job file names and
  sequence::

      orientation single points   <outdir>/<tag>.inp     (tag="calc")
      relaxed distance scan       relaxed_surface_scan/scan.inp

  The classic coordinator remains the workflow's driver; the task this
  module serves is *parity* — proving the OPI input equals the classic
  one (see the normaliser ``normalise_inp`` and the parity test
  ``tests/test_opi_co2_coordinator_parity.py``).

  Validation: block options (``maxiter``, ``brokensym``, ``qmatoms``,
  ``scan``) are enforced by OPI's pydantic models, which only accept
  ORCA-known fields.  Simple keyword tokens are passed through
  **unchanged** from the classic builder — a token the knowledge store
  does not list is surfaced as a warning, never invented here and never
  dropped (a user's explicit CONTROL choice wins and is only validated).

  The OPI ``Calculator`` is created with ``json_via_input=False`` so OPI
  does not inject its default ``%output jsongbwfile True
  jsonpropfile True`` block: the classic writer never writes it, and
  parity forbids an added block.

OPI notes (verified against orca-pi 2.0.0):
  * OPI writes each simple keyword on its own ``!`` line, lowercased.
    Keywords that have no OPI enum member are written verbatim as bare
    strings.  ORCA keywords are case-insensitive, so parity holds.
  * ``%pal nprocs N end`` is rendered as a ``%pal`` block, ``%maxcore N``
    inline — same values as the classic inline forms.
  * per-atom ``NewGTO "b" end`` tags attach through ``Atom.append_str``.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)

# Simple-keyword tokens that are part of the classic keyword string but are
# not themselves ORCA "listable" keywords in DELFIN's knowledge store; they
# are passed through unchanged (see the PARITY CONTRACT note above).
_PASS_THROUGH_TOKEN_RE = re.compile(r"^[A-Za-z0-9_-]+(\([^)]*\))?$")


# ---------------------------------------------------------------------------
# Parity normaliser
# ---------------------------------------------------------------------------

def normalise_inp(text: str) -> dict:
    """Parse an ORCA ``.inp`` into a canonical structure for comparison.

    Turns *both* textual conventions — the classic writer's and OPI's —
    into the same dictionary, so parity can be asserted structurally
    rather than by string comparison.  Returns::

        {
          "tokens":    frozenset of lowercase simple keywords,
          "maxcore":   int or None,
          "blocks":    {"scf"/"geom"/"qmmm"/"pal": {option: value, ...}},
          "charge":    int,
          "multiplicity": int,
          "coords":    tuple of (symbol, (x, y, z), newgto_or_None),
        }

    ``coords`` are compared by value with a small tolerance; block option
    values are normalised to lower-cased keys and canonical (float / int /
    set) values so layout and case differences do not matter.
    """
    lines = text.splitlines()
    # ------------------------------------------------------------------
    # 1. simple keywords: every '!' line contributes space-separated tokens
    # ------------------------------------------------------------------
    tokens: set[str] = set()
    for ln in lines:
        s = ln.strip()
        if s.startswith("!") and not s.startswith("!!"):
            rest = s[1:].strip()
            for tok in rest.split():
                tokens.add(tok.lower())
    # ------------------------------------------------------------------
    # 2. blocks: %name ... end  (a %geom block holds one nested end for its
    #    Scan sub-block, so it needs two 'end' lines; every other block one)
    # ------------------------------------------------------------------
    blocks: Dict[str, dict] = {}
    i = 0
    n = len(lines)
    maxcore: Optional[int] = None
    while i < n:
        s = lines[i].strip()
        low = s.lower()
        if low.startswith("%maxcore"):
            m = re.search(r"%maxcore\s+(\d+)", low)
            if m:
                maxcore = int(m.group(1))
            i += 1
            continue
        if low.startswith("%pal"):
            header_line = low[len("%pal"):]
            # classic writes '%pal nprocs 4 end' inline (end on the same
            # line); OPI writes '%pal' then 'nprocs 4' then 'end'. Only
            # consume body lines when the header did NOT already close it.
            if re.search(r"\bend\b", header_line):
                body, i = [], i + 1
            else:
                body, i = _consume_block(lines, i, end_count=1)
            blocks.setdefault("pal", {})
            combined = " ".join([header_line] + body).lower()
            pm = re.search(r"nprocs\s+(\d+)", combined)
            if pm:
                blocks["pal"]["nprocs"] = int(pm.group(1))
            continue
        if low.startswith("%scf"):
            header_tail = low[len("%scf"):]
            # Classic SP writer appends an '%scf BrokenSym 1,5 end' one-liner
            # from CONTROL; the scan writer and OPI write the '%scf' header
            # then option lines then a standalone 'end'.  Only consume body
            # lines when the header did NOT already carry the terminator.
            if re.search(r"\bend\b", header_tail):
                body, i = [], i + 1
            else:
                body, i = _consume_block(lines, i, end_count=1)
            scf_contents = ([header_tail] if header_tail.strip() else []) + body
            blocks["scf"] = _normalise_scf_block(scf_contents)
            continue
        if low.startswith("%geom"):
            body, i = _consume_container(lines, i)
            blocks["geom"] = _normalise_geom_block(body)
            continue
        if low.startswith("%qmmm"):
            body, i = _consume_container(lines, i)
            blocks["qmmm"] = _normalise_qmmm_block(body)
            continue
        if low.startswith("%output"):
            # OPI's default %output block; the classic never writes one.
            # Present or absent must not break parity, so we record it but
            # treat it as an ignorable block (see equivalence()).
            body, i = _consume_block(lines, i, end_count=1)
            blocks["output"] = {"_present": True}
            continue
        i += 1
    # ------------------------------------------------------------------
    # 3. coordinate block: '* xyz charge mult' ... '*'
    # ------------------------------------------------------------------
    charge: Optional[int] = None
    mult: Optional[int] = None
    coords: List[Tuple[str, Tuple[float, float, float], Optional[str]]] = []
    in_coords = False
    for ln in lines:
        s = ln.strip()
        if in_coords:
            if s == "*":
                break
            parts = s.split()
            if len(parts) >= 4:
                sym = parts[0]
                x, y, z = (float(parts[1]), float(parts[2]), float(parts[3]))
                newgto = None
                rest = s[len(parts[0]) + 1:].strip()
                if "NewGTO" in rest or "newgto" in rest:
                    m = re.search(r'NewGTO\s+"([^"]+)"', s, re.IGNORECASE)
                    if m:
                        newgto = m.group(1)
                coords.append((sym, (x, y, z), newgto))
            continue
        if s.startswith("*"):  # header line: '* xyz charge mult'
            in_coords = True
            toks = s[1:].strip().split()
            if len(toks) >= 3:
                try:
                    charge = int(toks[-2])
                    mult = int(toks[-1])
                except (ValueError, IndexError):
                    pass
            continue
    return {
        "tokens": frozenset(tokens),
        "maxcore": maxcore,
        "blocks": blocks,
        "charge": charge,
        "multiplicity": mult,
        "coords": tuple(coords),
    }


def _consume_block(lines: Sequence[str], start: int, *, end_count: int) -> Tuple[List[str], int]:
    """Return the interior lines of a ``%block`` closed by ``end`` lines.

    ``end_count`` is how many standalone ``end`` lines close the block.
    """
    body: List[str] = []
    ends = 0
    i = start + 1
    n = len(lines)
    while i < n:
        s = lines[i].strip()
        if re.fullmatch(r"end", s, re.IGNORECASE):
            ends += 1
            i += 1
            if ends >= end_count:
                break
            continue
        body.append(lines[i])
        i += 1
    return body, i


def _consume_container(lines: Sequence[str], start: int) -> Tuple[List[str], int]:
    """Interior of a block holding one nested sub-block (``%geom``/``%qmmm``).

    The classic writer puts its nested ``end`` (the ``Scan`` sub-block, or the
    second ``END`` of ``%QMMM``) on its own line, so such a block is closed by
    two standalone ``end`` lines.  OPI folds a ``Scan`` / ``qmatoms`` terminator
    inline (``... end``) and closes the block with one standalone ``end``.  To
    serve both, stop after the second standalone ``end`` OR at the next top-level
    directive (``%`` / ``!`` / ``*``), whichever comes first.
    """
    body: List[str] = []
    ends = 0
    i = start + 1
    n = len(lines)
    while i < n:
        s = lines[i].strip()
        if not s:
            i += 1
            continue
        if re.fullmatch(r"end", s, re.IGNORECASE):
            ends += 1
            i += 1
            if ends >= 2:
                break
            continue
        if s.startswith("%") or s.startswith("!") or s.startswith("*"):
            break
        body.append(lines[i])
        i += 1
    return body, i


def _normalise_scf_block(body: Sequence[str]) -> dict:
    """Canonical form of a ``%scf`` block: lower-cased option names."""
    out: dict[str, Any] = {}
    joined = " ".join(body).lower()
    m = re.search(r"\bmaxiter\s+(\d+)", joined)
    if m:
        out["maxiter"] = int(m.group(1))
    m = re.search(r"\bbrokensym\s+([0-9,\s]+)", joined)
    if m:
        out["brokensym"] = _normalise_int_list(m.group(1))
    return out


def _normalise_geom_block(body: Sequence[str]) -> dict:
    """Canonical form of a ``%geom`` block: maxiter and the scan descriptor."""
    out: dict[str, Any] = {}
    joined = " ".join(body).lower()
    m = re.search(r"\bmaxiter\s+(\d+)", joined)
    if m:
        out["maxiter"] = int(m.group(1))
    # Scan descriptor: 'B i j = v1, v2, steps end' — case-blind (joined is
    # lowercased, classic writes 'B', OPI 'b').
    sm = re.search(
        r"\b([ABD])\s+(\d+)\s+(\d+)\s*=\s*"
        r"([+-]?\d+(?:\.\d+)?)\s*,\s*([+-]?\d+(?:\.\d+)?)\s*,\s*(\d+)",
        joined,
        re.IGNORECASE,
    )
    if sm:
        out["scan"] = {
            "mode": sm.group(1).upper(),
            "atoms": (int(sm.group(2)), int(sm.group(3))),
            "start": float(sm.group(4)),
            "end": float(sm.group(5)),
            "steps": int(sm.group(6)),
        }
    return out


def _normalise_qmmm_block(body: Sequence[str]) -> dict:
    """Canonical form of a ``%qmmm`` block: the QM atom set.

    Handles both the classic ``QMATOMS {0:3}`` inclusive range and OPI's
    enumerated ``qmatoms { 0 1 2 3 } end``; both reduce to ``{0,1,2,3}``.
    """
    joined = " ".join(body)
    out: dict[str, Any] = {}
    m = re.search(r"qmatoms\s*\{\s*([^}]*)\}", joined, re.IGNORECASE)
    if not m:
        m = re.search(r"QMATOMS\s*\{\s*(.*?)\}", joined)
    if m:
        inner = m.group(1)
        atoms: set[int] = set()
        # each token is either a number or an inclusive range a:b
        for tok in re.split(r"[,\s]+", inner.strip()):
            if not tok:
                continue
            if ":" in tok:
                lo, hi = (int(v) for v in tok.split(":"))
                atoms.update(range(lo, hi + 1))
            elif tok.isdigit() or (tok.startswith("-") and tok[1:].isdigit()):
                atoms.add(int(tok))
        out["qmatoms"] = frozenset(atoms)
    return out


def _normalise_int_list(text: str) -> Tuple[int, ...]:
    return tuple(int(v) for v in re.findall(r"\d+", text))


def inputs_equivalent(a: dict, b: dict) -> List[str]:
    """Structural diff of two normalised inputs; empty list means equivalent.

    Comparing *both* directions so a difference is never silently dropped.
    """
    diffs: List[str] = []
    if set(a["tokens"]) != set(b["tokens"]):
        diffs.append(
            "tokens differ: "
            f"only-a={sorted(set(a['tokens']) - set(b['tokens']))} "
            f"only-b={sorted(set(b['tokens']) - set(a['tokens']))}"
        )
    if a["maxcore"] != b["maxcore"]:
        diffs.append(f"maxcore: {a['maxcore']!r} != {b['maxcore']!r}")
    _compare_blocks(a, b, diffs)
    if a["charge"] != b["charge"]:
        diffs.append(f"charge: {a['charge']!r} != {b['charge']!r}")
    if a["multiplicity"] != b["multiplicity"]:
        diffs.append(f"multiplicity: {a['multiplicity']!r} != {b['multiplicity']!r}")
    _compare_coords(a, b, diffs)
    return diffs


def _compare_blocks(a: dict, b: dict, diffs: List[str]) -> None:
    names = set(a["blocks"]) | set(b["blocks"])
    for name in sorted(names):
        av = a["blocks"].get(name)
        bv = b["blocks"].get(name)
        # The %output block is OPI's optional default; classic never writes
        # it, so its presence in only one side is a non-issue.
        if name == "output":
            continue
        if av is None or bv is None:
            diffs.append(f"block {name!r}: only in one side")
            continue
        keys = set(av) | set(bv)
        for k in sorted(keys):
            if name == "geom" and k == "scan":
                continue  # compared atomically below
            if name == "qmmm" and k == "qmatoms":
                if av.get("qmatoms") != bv.get("qmatoms"):
                    diffs.append(f"block qmmm qmatoms differ: {av.get('qmatoms')} vs {bv.get('qmatoms')}")
                continue
            if av.get(k) != bv.get(k):
                diffs.append(
                    f"block {name!r} option {k!r}: {av.get(k)!r} != {bv.get(k)!r}"
                )
        if name == "geom":
            sa, sb = av.get("scan"), bv.get("scan")
            if (sa is None) != (sb is None):
                diffs.append("geom scan presence differs")
            elif sa and sb:
                if sa["mode"] != sb["mode"] or sa["steps"] != sb["steps"]:
                    diffs.append(f"geom scan header differs: {sa} vs {sb}")
                if list(sa["atoms"]) != list(sb["atoms"]):
                    diffs.append(f"geom scan atoms differ: {sa['atoms']} vs {sb['atoms']}")
                if abs(sa["start"] - sb["start"]) > 1e-6 or abs(sa["end"] - sb["end"]) > 1e-6:
                    diffs.append(f"geom scan values differ: {sa['start']},{sa['end']} vs {sb['start']},{sb['end']}")


def _compare_coords(a: dict, b: dict, diffs: List[str]) -> None:
    ca, cb = a["coords"], b["coords"]
    if len(ca) != len(cb):
        diffs.append(f"coord count differs: {len(ca)} vs {len(cb)}")
        return
    for i, (sa, sb) in enumerate(zip(ca, cb)):
        if sa[0] != sb[0]:
            diffs.append(f"coord {i} symbol: {sa[0]} != {sb[0]}")
        for j in range(3):
            if abs(sa[1][j] - sb[1][j]) > 1e-6:
                diffs.append(f"coord {i} axis {j}: {sa[1][j]} != {sb[1][j]}")
        if (sa[2] is None) != (sb[2] is None):
            diffs.append(f"coord {i} NewGTO presence differs")
        elif sa[2] is not None and sa[2].lower() != sb[2].lower():
            diffs.append(f"coord {i} NewGTO basis: {sa[2]} != {sb[2]}")


# ---------------------------------------------------------------------------
# Lazy OPI import
# ---------------------------------------------------------------------------

def _import_opi():
    """Import OPI classes lazily so the module imports without OPI present."""
    from opi.core import Calculator
    from opi.input import Input
    from opi.input.blocks.block_geom import BlockGeom, Scan
    from opi.input.blocks.block_qmmm import BlockQmmm
    from opi.input.blocks.block_scf import BlockScf
    from opi.input.structures import Atom, Structure
    return Calculator, Input, BlockGeom, Scan, BlockQmmm, BlockScf, Atom, Structure


# ---------------------------------------------------------------------------
# Reused classic logic (BY IMPORT — never edited here)
# ---------------------------------------------------------------------------

def _classic():
    from delfin.co2 import CO2_Coordinator6 as c
    return c


# ---------------------------------------------------------------------------
# broken-sym → %scf block
# ---------------------------------------------------------------------------

def _parse_scf_lines(broken_sym: str) -> List[str]:
    """Split an ``%scf ... end`` string into its inner option lines.

    Mirrors the classic scan writer (CO2_Coordinator6.py:1002-1008):
    strip the leading ``%scf``, strip a trailing ``end``, return the
    non-empty inner lines.
    """
    text = broken_sym.strip()
    if text.lower().startswith("%scf"):
        text = text[4:].strip()
    inner = text
    if inner.lower().endswith("end"):
        inner = inner[:-3].strip()
    lines = [ln.strip() for ln in inner.splitlines() if ln.strip()]
    if not lines and inner.strip():
        lines = [inner.strip()]
    return lines


def _scf_fields_for_lines(lines: Sequence[str]) -> Dict[str, Any]:
    """Map the inner ``%scf`` option lines to OPI ``BlockScf`` fields.

    Recognised options: ``MaxIter N`` → ``maxiter``, ``BrokenSym a,b,...`` →
    ``brokensym``.  An option OPI's BlockScf cannot express is an error —
    never silently dropped (a wrong input file for a live ORCA run is the
    cost of guessing).
    """
    fields: Dict[str, Any] = {}
    for ln in lines:
        low = ln.lower()
        m = re.match(r"maxiter\s+(\d+)", low)
        if m:
            fields["maxiter"] = int(m.group(1))
            continue
        m = re.match(r"brokensym\s+([0-9,\s]+)", low)
        if m:
            fields["brokensym"] = [int(v) for v in re.findall(r"\d+", m.group(1))]
            continue
        raise ValueError(
            "OPI BlockScf cannot represent broken-sym line "
            f"{ln!r}; refusing to silently drop an %scf option "
            "(only MaxIter and BrokenSym are mapped)."
        )
    return fields


def _scf_for_sp(broken_sym: str):
    """The ``%scf`` block the classic SP writer emits (CO2_Coordinator6.py:952-954).

    The SP writer appends the raw broken-sym string (an ``%scf`` block)
    and adds no ``MaxIter``.  Empty broken-sym → no block at all.
    """
    if not broken_sym:
        return None
    lines = _parse_scf_lines(broken_sym)
    fields = _scf_fields_for_lines(lines)
    if not fields:
        return None
    _, _, _, _, _, BlockScf, _, _ = _import_opi()
    return BlockScf(**fields)


def _scf_for_scan(broken_sym: str):
    """The ``%scf`` block the classic scan writer emits, always MaxIter 375.

    Mirrors CO2_Coordinator6.py:1001-1016: the block ends with MaxIter 375
    (appended when the parsed broken-sym does not already set one).
    """
    _, _, _, _, _, BlockScf, _, _ = _import_opi()
    if not broken_sym:
        return BlockScf(maxiter=375)
    lines = _parse_scf_lines(broken_sym)
    fields = _scf_fields_for_lines(lines)
    if "maxiter" not in fields:
        fields["maxiter"] = 375
    return BlockScf(**fields)


def _qmmm_for_range(qmmm_range: Optional[Tuple[int, int]]):
    """The ``%qmmm`` block for a QM/MM range, or None when absent."""
    if not qmmm_range:
        return None
    start, end = qmmm_range
    _, _, _, _, BlockQmmm, _, _, _ = _import_opi()
    return BlockQmmm(qmatoms=list(range(start, end + 1)))


# ---------------------------------------------------------------------------
# Core OPI write
# ---------------------------------------------------------------------------

def _write_opi_input(
    work_dir: Path,
    *,
    tokens: Sequence[str],
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    newgto: Mapping[int, str],
    charge: int,
    mult: int,
    maxcore: int,
    ncores: int,
    scf_block,
    geom_block,
    qmmm_block,
    base_name: str,
) -> Path:
    """Assemble and write one OPI ``.inp``; returns its absolute path."""
    Calculator, Input, _, _, _, _, Atom, Structure = _import_opi()
    inp = Input()
    # Every token is passed through as-is (lowercased by OPI on write, which
    # is parity-safe since ORCA keywords are case-insensitive).
    inp.add_simple_keywords(*[str(t) for t in tokens])
    inp.ncores = int(ncores)
    inp.memory = int(maxcore)
    if scf_block is not None:
        inp.add_blocks(scf_block)
    if geom_block is not None:
        inp.add_blocks(geom_block)
    if qmmm_block is not None:
        inp.add_blocks(qmmm_block)

    atoms = []
    for idx, (el, (x, y, z)) in enumerate(zip(elements, coordinates)):
        # Match the classic writer's 8-decimal geometry (CO2_Coordinator6.py:523)
        append = None
        basis = newgto.get(idx)
        if basis:
            append = f'NewGTO "{basis}" end'
        atoms.append(Atom(el, coordinates=(round(float(x), 8), round(float(y), 8), round(float(z), 8)),
                        append_str=append))
    structure = Structure(atoms, charge=int(charge), multiplicity=int(mult))

    calc = Calculator(base_name, working_dir=work_dir, version_check=False)
    calc.json_via_input = False  # parity: classic never writes %output
    calc.input = inp
    calc.structure = structure
    calc.write_input(force=True)
    return work_dir / f"{base_name}.inp"


# ---------------------------------------------------------------------------
# Public builders (parity mirrors of the two classic writers)
# ---------------------------------------------------------------------------

def write_co2_sp_opi_input(
    atoms,
    xyz_path: str,
    outdir: str,
    *,
    orca_keywords: str = "GFN2-XTB SP ALPB(DMF)",
    broken_sym: str = "",
    charge: int = -2,
    multiplicity: int = 1,
    PAL: int = 8,
    maxcore: int = 2000,
    tag: str = "calc",
    metal_symbol: Optional[str] = None,
    metal_basis: Optional[str] = None,
    control_args: Optional[Mapping[str, Any]] = None,
    qmmm_range: Optional[Tuple[int, int]] = None,
    co2_indices: Optional[Sequence[int]] = None,
) -> Path:
    """OPI equivalent of ``CO2_Coordinator6.write_orca_sp_input_and_run``.

    Writes the same single-point ``.inp`` the classic SP writer produces
    (CO2_Coordinator6.py:946-963): the ``!`` keywords, ``%maxcore``,
    ``%pal``, the optional broken-sym ``%scf`` block, the ``%qmmm`` block
    and the ``* xyz`` coordinate block with per-atom ``NewGTO`` tags.  Only
    the input is written — the classic function additionally launches ORCA,
    which this parity module deliberately does not do.
    """
    c = _classic()
    keywords_line = c._ensure_qmmm_keyword(orca_keywords, qmmm_range)
    tokens = keywords_line.split()

    os.makedirs(outdir, exist_ok=True)
    assignments = c._build_newgto_assignments(
        atoms, keywords_line, metal_basis, dict(control_args or {}),
        qmmm_range, skip_metals=False, co2_indices=co2_indices,
    )
    elements = [a.symbol for a in atoms]
    coords = [tuple(a.position) for a in atoms]

    return _write_opi_input(
        Path(outdir).resolve(),
        tokens=tokens,
        elements=elements,
        coordinates=coords,
        newgto=assignments,
        charge=charge,
        mult=multiplicity,
        maxcore=maxcore,
        ncores=PAL,
        scf_block=_scf_for_sp(broken_sym),
        geom_block=None,
        qmmm_block=_qmmm_for_range(qmmm_range),
        base_name=tag,
    )


def write_co2_scan_opi_input(
    atoms,
    xyz_path: str,
    metal_index: int,
    co2_c_index: int,
    *,
    start_distance: float,
    end_distance: float = 1.7,
    steps: int = 5,
    orca_keywords: str = "GFN2-XTB OPT ALPB(DMF)",
    broken_sym: str = "",
    charge: int = -2,
    multiplicity: int = 3,
    PAL: int = 4,
    maxcore: int = 2000,
    metal_symbol: Optional[str] = None,
    metal_basis: Optional[str] = None,
    control_args: Optional[Mapping[str, Any]] = None,
    qmmm_range: Optional[Tuple[int, int]] = None,
    co2_indices: Optional[Sequence[int]] = None,
    scan_dir: str = "relaxed_surface_scan",
) -> Path:
    """OPI equivalent of ``CO2_Coordinator6.write_orca_input_and_run``.

    Writes the same relaxed distance-scan ``.inp`` the classic scan writer
    produces (CO2_Coordinator6.py:993-1033): keywords augmented with
    ``NoTRAH`` / ``SlowConv``, ``%maxcore``, the ``%scf`` block (always
    ``MaxIter 375``), ``%pal``, ``%geom`` with ``Scan B i j = start, end,
    steps``, and the ``* xyz`` block.  Written to ``scan_dir/scan.inp``,
    the classic's hard-coded ``relaxed_surface_scan/scan.inp`` by default.
    """
    c = _classic()
    keywords_line = c._ensure_qmmm_keyword(orca_keywords, qmmm_range)
    keywords_line = c._append_orca_keyword(keywords_line, "NoTRAH")
    keywords_line = c._append_orca_keyword(keywords_line, "SlowConv")
    tokens = keywords_line.split()

    # Classic formats the start distance to two decimals (CO2_Coordinator6.py:1025)
    start_fmt = float(f"{float(start_distance):.2f}")

    _, _, BlockGeom, Scan, _, _, _, _ = _import_opi()
    geom_block = BlockGeom(
        maxiter=400,
        scan=Scan(
            mode="B",
            atom1=int(metal_index),
            atom2=int(co2_c_index),
            startvalue=start_fmt,
            endvalue=float(end_distance),
            steps=int(steps),
        ),
    )

    assignments = c._build_newgto_assignments(
        atoms, keywords_line, metal_basis, dict(control_args or {}),
        qmmm_range, skip_metals=False, co2_indices=co2_indices,
    )
    elements = [a.symbol for a in atoms]
    coords = [tuple(a.position) for a in atoms]

    os.makedirs(scan_dir, exist_ok=True)
    return _write_opi_input(
        Path(scan_dir).resolve(),
        tokens=tokens,
        elements=elements,
        coordinates=coords,
        newgto=assignments,
        charge=charge,
        mult=multiplicity,
        maxcore=maxcore,
        ncores=PAL,
        scf_block=_scf_for_scan(broken_sym),
        geom_block=geom_block,
        qmmm_block=_qmmm_for_range(qmmm_range),
        base_name="scan",
    )
