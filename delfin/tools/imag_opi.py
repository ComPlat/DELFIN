"""IMAG expressed on OPI's building blocks, with input parity.

The classic IMAG writer (``delfin.imag``) is a transformer of the step's own
ORCA input: the single-point candidate is ``delfin.imag._candidate_job``
(imag.py:268-285, called at imag.py:410/419) and the re-optimisation input is
``_with_base(_with_geometry(job, atoms)) + rest`` (imag.py:551-554).  Both
re-emit only what the step input already contains — method, per-atom basis,
``%scf``, ``%tddft``, appended ``$new_job`` jobs — with OPT/FREQ tokens taken
off the candidate's '!' line and the coordinates displaced.

``delfin.tools.imag_opi`` runs the SAME classic transform functions BY IMPORT
(not copied) and then re-emits the transformed text with OPI 2.0.0 input
classes.  The acceptance criterion is parity: for the same source input and
displacement the OPI file must be structurally equivalent to the classic one
— same simple keywords (order/case-insensitive), same blocks with the same
options and values, same charge/multiplicity/coordinates (per-atom ``NewGTO``
tails included), same ``%pal``/``%maxcore``/``%base``, same file names.  The
parity test ``tests/test_opi_imag_parity.py`` proves it with
``delfin.tools.occupier_opi.normalise_inp`` / ``inputs_equivalent`` — the
same normalisation the OCCUPIER parity test uses, reused so the contract
cannot drift between the writers.

Emission mapping (every decision quotes the classic line it reproduces):
- '!' tokens and ``%maxcore``/``%pal nprocs`` → OPI simple keywords / memory /
  ncores (imag.py keeps the step input's own lines untouched).
- ``%scf`` → OPI ``BlockScf`` via ``delfin.tools.occupier_opi.
  _scf_fields_for_lines``: MaxIter/BrokenSym/APMethod map; any other option
  raises — an ``%scf`` option is never silently dropped (the OCCUPIER
  precedent, tests/test_opi_occupier_parity.py::test_unknown_scf_option…).
- ``%moinp "path"`` → arbitrary string, as the OCCUPIER OPI writer does
  (OPI's ``Input.moinp`` requires a real subpath).
- Any other ``%`` block (``%tddft``, ``%cpcm``, ``%geom``, …) → verbatim
  arbitrary string before the coordinates: OPI 2.0.0 has no block class for
  these options and the classic writer passes the step input's text through
  unchanged, so a re-encoding could only lose information.
- ``%base "<name>"`` → arbitrary string at the top; the classic writer writes
  the base line inside the file (imag.py:279-284) and OPI's ``Calculator``
  only derives the file name from its base argument.
- coordinates → OPI ``Atom``/``Structure`` with the classic row tail kept as
  ``append_str`` (imag.py:239-265 keeps everything after x y z on each line,
  "so a metal's per-atom NewGTO survives").
- everything from the first ``$new_job`` → verbatim arbitrary string at the
  bottom (the rerun's appended jobs are carried unchanged, imag.py:551-554).

The classic module stays the workflow's driver; this module only writes
inputs (no ORCA launch).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from delfin.common.logging import get_logger
from delfin.imag import (
    _candidate_job,
    _split_first_job,
    _with_base,
    _with_geometry,
)

logger = get_logger(__name__)

__all__ = [
    "ImagOpiError",
    "write_imag_candidate_opi_input",
    "write_imag_rerun_opi_input",
]


class ImagOpiError(ValueError):
    """The OPI writer cannot reproduce the classic IMAG input faithfully."""


def _import_opi():
    """OPI's input classes (lazy; the module imports without OPI present)."""
    try:
        from opi.core import Calculator
        from opi.input import Input
        from opi.input.arbitrary_string import ArbitraryString, ArbitraryStringPos
        from opi.input.blocks.block_scf import BlockScf
        from opi.input.structures import Atom, Structure
    except Exception as exc:  # pragma: no cover - exercised via ImportError paths
        raise ImagOpiError(
            "OPI is required to write IMAG inputs on OPI's building "
            f"blocks, but 'import opi' failed: {exc}"
        ) from exc
    return {
        "Calculator": Calculator,
        "Input": Input,
        "ArbitraryString": ArbitraryString,
        "ArbitraryStringPos": ArbitraryStringPos,
        "BlockScf": BlockScf,
        "Atom": Atom,
        "Structure": Structure,
    }


# ---------------------------------------------------------------------------
# Parsing the classic transform's output (the ground truth being re-emitted)
# ---------------------------------------------------------------------------

_SCF_UNMAPPED = (
    "OPI BlockScf cannot represent %scf line {line!r}; refusing to "
    "silently drop an %scf option (only MaxIter, BrokenSym and APMethod "
    "are mapped)."
)


def _consume_scf(lines: List[str], start: int) -> Tuple[List[str], int]:
    """Interior lines of a ``%scf`` block (standalone ``end`` closes it)."""
    body: List[str] = []
    i = start + 1
    while i < len(lines):
        s = lines[i].strip()
        if re.fullmatch(r"end", s, re.IGNORECASE):
            return body, i + 1
        body.append(lines[i])
        i += 1
    raise ImagOpiError("%scf block without end")


def _consume_passthrough(lines: List[str], start: int) -> Tuple[List[str], int]:
    """Verbatim lines of a ``%`` block OPI has no class for.

    Nesting by directive: a ``%name`` line opens one level, a standalone
    ``end`` closes one.  Inline ``%pal nprocs 4 end`` never reaches here
    (handled separately); nested containers (``%geom`` constraints) close
    with their own ``end`` first.
    """
    body: List[str] = [lines[start]]
    depth = 1
    i = start + 1
    while i < len(lines):
        s = lines[i].strip()
        if re.match(r"%(?!\w*(?:end)\b)", s) and re.match(r"^\s*%\w+", s):
            depth += 1
        elif re.fullmatch(r"end", s, re.IGNORECASE):
            depth -= 1
            body.append(lines[i])
            i += 1
            if depth <= 0:
                return body, i
            continue
        body.append(lines[i])
        i += 1
    raise ImagOpiError(f"{lines[start].strip()!r} block without end")


def _parse_first_job(job: str) -> Dict[str, Any]:
    """The first job of a classic IMAG transform, into canonical parts."""
    lines = job.splitlines()
    tokens: List[str] = []
    maxcore: Optional[int] = None
    nprocs: Optional[int] = None
    scf_lines: List[str] = []
    passthrough: List[str] = []
    moinp: Optional[str] = None
    base_line: Optional[str] = None
    charge: Optional[int] = None
    multiplicity: Optional[int] = None
    coords: List[Tuple[str, Tuple[float, float, float], str]] = []
    rest: List[str] = []
    in_coords = False
    i = 0
    n = len(lines)
    while i < n:
        s = lines[i].strip()
        if in_coords:
            if s == "*":
                in_coords = False
                i += 1
                continue
            parts = s.split()
            if len(parts) >= 4:
                sym = parts[0]
                xyz = (float(parts[1]), float(parts[2]), float(parts[3]))
                tail = " ".join(parts[4:])
                coords.append((sym, xyz, tail))
            i += 1
            continue
        if s.startswith("$"):
            # appended jobs from here on: carried verbatim to the bottom
            rest = lines[i:]
            break
        low = s.lower()
        if s.startswith("!") and not s.startswith("!!"):
            tokens.extend(s[1:].split())
            i += 1
            continue
        if low.startswith("%maxcore"):
            m = re.search(r"%maxcore\s+(\d+)", low)
            if m:
                maxcore = int(m.group(1))
            i += 1
            continue
        if low.startswith("%pal"):
            if re.search(r"\bend\b", low[len("%pal"):]):
                combined = low
                i += 1
            else:
                body, i = _consume_scf(lines, i)  # same standalone-end shape
                combined = low + " " + " ".join(body).lower()
            pm = re.search(r"nprocs\s+(\d+)", combined)
            if pm:
                nprocs = int(pm.group(1))
            continue
        if low.startswith("%scf"):
            body, i = _consume_scf(lines, i)
            scf_lines.extend(ln.strip() for ln in body if ln.strip())
            continue
        if low.startswith("%moinp"):
            mm = re.search(r'%moinp\s+"([^"]+)"', s, re.IGNORECASE)
            if mm:
                moinp = mm.group(1)
            i += 1
            continue
        if low.startswith("%base"):
            base_line = s
            i += 1
            continue
        if s.startswith("*"):
            # geometry header: '* xyz charge multiplicity [unit]'
            toks = s[1:].strip().split()
            if len(toks) >= 3 and toks[0].lower() == "xyz":
                try:
                    charge = int(toks[1])
                    multiplicity = int(toks[2])
                except (ValueError, IndexError) as exc:
                    raise ImagOpiError(
                        f"cannot read charge/multiplicity from {s!r}"
                    ) from exc
                in_coords = True
                i += 1
                continue
            raise ImagOpiError(
                f"geometry header {s!r} is not '* xyz c m'; the IMAG writer "
                "only re-emits inline-coordinate jobs"
            )
        if s.startswith("%"):
            body, i = _consume_passthrough(lines, i)
            passthrough.append("\n".join(body))
            continue
        i += 1
    if charge is None or multiplicity is None:
        raise ImagOpiError("no '* xyz charge multiplicity' header in the job")
    if not coords:
        raise ImagOpiError("no coordinates in the job")
    return {
        "tokens": tokens,
        "maxcore": maxcore,
        "nprocs": nprocs,
        "scf_lines": scf_lines,
        "passthrough": passthrough,
        "moinp": moinp,
        "base_line": base_line,
        "charge": charge,
        "multiplicity": multiplicity,
        "coords": coords,
        "rest": "\n".join(rest),
    }


def _scf_block_from_lines(scf_lines: Sequence[str], opi: Mapping[str, Any]):
    """``BlockScf`` for the classic job's ``%scf`` lines, or None.

    Reuses ``delfin.tools.occupier_opi._scf_fields_for_lines`` so the option
    mapping — and its refusal to drop an unmappable option — is the OCCUPIER
    writer's, not a second implementation.
    """
    from delfin.tools.occupier_opi import _scf_fields_for_lines

    if not scf_lines:
        return None
    fields = _scf_fields_for_lines(scf_lines)
    if "brokensym" in fields and opi["BlockScf"] is not None:
        from opi.input.blocks.util import NumList

        fields["brokensym"] = NumList(fields["brokensym"])
    return opi["BlockScf"](**fields)


# ---------------------------------------------------------------------------
# Core OPI write
# ---------------------------------------------------------------------------

def _write_imag_opi_input(
    work_dir,
    *,
    classic_text: str,
    base: str,
) -> Path:
    """Re-emit a classic IMAG transform's text with OPI's input classes.

    ``classic_text`` is what ``delfin.imag`` itself writes (the candidate job
    or the rerun input); the OPI file must be structurally equivalent to it.
    """
    opi = _import_opi()
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    job, _rest_from_split = _split_first_job(classic_text)
    parsed = _parse_first_job(job)
    rest = parsed["rest"] or _rest_from_split

    Input = opi["Input"]
    Atom = opi["Atom"]
    Structure = opi["Structure"]
    Calculator = opi["Calculator"]

    inp = Input()
    inp.add_simple_keywords(*parsed["tokens"])
    if parsed["maxcore"] is not None:
        inp.memory = int(parsed["maxcore"])
    if parsed["nprocs"] is not None:
        inp.ncores = int(parsed["nprocs"])

    # %base: the classic transform writes it inside the file
    # (imag.py:279-284); OPI's Calculator only names the file from its base
    # argument, so the line is carried as a top arbitrary string.  When the
    # classic text has none the writer's base argument is used.
    base_line = parsed["base_line"] or f'%base "{base}"'
    inp.add_arbitrary_string(base_line, pos=opi["ArbitraryStringPos"].TOP)

    # unmapped % blocks verbatim, before the coordinates (see module docstring)
    for block_text in parsed["passthrough"]:
        inp.add_arbitrary_string(block_text)

    # %moinp as a raw line: OPI's Input.moinp requires a real subpath, the
    # classic writer writes the line unconditionally (occupier_opi precedent)
    if parsed["moinp"]:
        inp.add_arbitrary_string(f'%moinp "{parsed["moinp"]}"')

    scf_block = _scf_block_from_lines(parsed["scf_lines"], opi)
    if scf_block is not None:
        inp.add_blocks(scf_block)

    atoms_opi = []
    for sym, (x, y, z), tail in parsed["coords"]:
        atoms_opi.append(
            Atom(
                str(sym),
                coordinates=(round(float(x), 8), round(float(y), 8),
                             round(float(z), 8)),
                append_str=tail or None,
            )
        )
    structure = Structure(
        atoms_opi,
        charge=int(parsed["charge"]),
        multiplicity=int(parsed["multiplicity"]),
    )

    calc = Calculator(base, working_dir=work_dir, version_check=False)
    calc.json_via_input = False  # parity: the classic writer never writes %output
    calc.input = inp
    calc.structure = structure
    if rest.strip():
        # appended $new_job jobs, verbatim at the bottom (imag.py:551-554)
        inp.add_arbitrary_string(rest, pos=opi["ArbitraryStringPos"].BOTTOM)
    calc.write_input(force=True)
    return work_dir / f"{base}.inp"


def write_imag_candidate_opi_input(
    work_dir,
    *,
    source_text: str,
    atoms: Sequence[Tuple[str, float, float, float]],
    base: str,
    optimise: bool = False,
) -> Path:
    """OPI equivalent of ``delfin.imag._candidate_job`` + write.

    The classic candidate text is produced by the classic function itself
    (imported, not copied) and then re-emitted with OPI's input classes.
    """
    classic_text = _candidate_job(source_text, list(atoms), base, optimise=optimise)
    return _write_imag_opi_input(work_dir, classic_text=classic_text, base=base)


def write_imag_rerun_opi_input(
    work_dir,
    *,
    source_text: str,
    atoms: Sequence[Tuple[str, float, float, float]],
    base: str,
) -> Path:
    """OPI equivalent of the re-optimisation input (imag.py:551-554):
    ``_with_base(_with_geometry(job, atoms)) + rest``, written with OPI."""
    job, rest = _split_first_job(source_text)
    classic_text = _with_base(_with_geometry(job, list(atoms)), base) + rest
    return _write_imag_opi_input(work_dir, classic_text=classic_text, base=base)
