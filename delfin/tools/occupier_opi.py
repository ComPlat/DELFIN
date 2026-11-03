"""OCCUPIER expressed on OPI's building blocks, with input parity.

``delfin.tools.occupier_opi`` writes the ORCA ``.inp`` files the classic
OCCUPIER writer (``delfin.occupier.read_and_modify_file_OCCUPIER``,
delfin/occupier.py:414-612) writes, but builds them with OPI's input
classes (``opi.input``).  The acceptance criterion is parity: for the
same CONTROL config, structure and broken-sym string, the OPI file must
be structurally equivalent to the classic one — same simple keywords
(order/case-insensitive), same blocks with the same options and values,
same charge/multiplicity/coordinates (per-atom ``NewGTO`` included),
same ``%pal``/``%maxcore``, same file names.  The parity test
``tests/test_opi_occupier_parity.py`` proves it over the repository's
own CONTROL settings.

The classic module remains the workflow's driver; this module only
writes inputs (no ORCA launch).  Every decision below quotes the classic
line it reproduces.

Reused BY IMPORT (not copied): ``delfin.occupier._is_composite_3c``,
``delfin.utils.resolve_level_of_theory``,
``delfin.common.solvation.build_solvation_keyword``,
``delfin.common.control_validator.resolve_occupier_compare`` — the exact
helpers the classic writer itself calls (delfin/occupier.py:484-498).
Selection logic (which spin states / occupations run, in which order)
stays in the classic module; this module writes one FoB input at a time
with the same signature shape (occupier.py:414-416).
"""

from __future__ import annotations

import re

from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)

__all__ = (
    "OccupierOpiError",
    "normalise_inp",
    "inputs_equivalent",
    "write_occupier_fob_opi_input",
)


class OccupierOpiError(ValueError):
    """A rule violation the caller must resolve; the rule is named."""


# ---------------------------------------------------------------------------
# Structural normaliser (same contract as the CO2 parity module)
# ---------------------------------------------------------------------------

def normalise_inp(text: str) -> dict:
    """Parse an ORCA ``.inp`` into a canonical structure for comparison.

    The block/coordinate parsing helpers are imported from
    ``delfin.tools.co2_coordinator_opi`` (the first parity module) so the
    structural comparison cannot drift between the two OPI writers; the
    OCCUPIER-specific extensions (``%moinp`` line, ``%freq Temp`` block)
    live here because the CO2 writer never emits those.

    Same contract as ``delfin.tools.co2_coordinator_opi.normalise_inp``
    (structurally equal dictionaries must come from the classic writer's
    and OPI's textual conventions).  Returns::

        {
          "tokens":    frozenset of lowercase simple keywords,
          "maxcore":   int or None,
          "blocks":    {"scf"/"freq"/"qmmm"/"pal"/...: {option: value}},
          "moinp":     lowercase path string or None,
          "charge":    int,
          "multiplicity": int,
          "coords":    tuple of (symbol, (x, y, z), newgto_or_None),
        }
    """
    import re

    lines = text.splitlines()
    # 1. simple keywords: every '!' line contributes tokens
    tokens: set[str] = set()
    for ln in lines:
        s = ln.strip()
        if s.startswith("!") and not s.startswith("!!"):
            for tok in s[1:].strip().split():
                tokens.add(tok.lower())
    # 2. blocks: %name ... end
    blocks: Dict[str, dict] = {}
    moinp: Optional[str] = None
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
            # classic: '%pal nprocs N end' inline; OPI: separate lines.
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
            if re.search(r"\bend\b", header_tail):
                body, i = [], i + 1
            else:
                body, i = _consume_block(lines, i, end_count=1)
            scf_contents = ([header_tail] if header_tail.strip() else []) + body
            # The classic OCCUPIER writer emits TWO %scf blocks (the
            # CONTROL maxiter at occupier.py:1068/:549 and the raw
            # broken-sym string at :1072-1073/:553-554); merge them so
            # the later block does not erase the earlier options.
            parsed = _normalise_scf_block(scf_contents)
            if "scf" in blocks:
                merged = dict(blocks["scf"])
                merged.update(parsed)
                blocks["scf"] = merged
            else:
                blocks["scf"] = parsed
            continue
        if low.startswith("%freq"):
            body, i = _consume_block(lines, i, end_count=1)
            joined = " ".join(body).lower()
            fm = re.search(r"\btemp\s+([0-9.eE+-]+)", joined)
            blocks["freq"] = {"temp": float(fm.group(1))} if fm else {}
            continue
        if low.startswith("%qmmm"):
            body, i = _consume_container(lines, i)
            blocks["qmmm"] = _normalise_qmmm_block(body)
            continue
        if low.startswith("%moinp"):
            # classic: '%moinp "path"'; OPI renders the same shape.
            mm = re.search(r'%moinp\s+"([^"]+)"', s, re.IGNORECASE)
            if mm:
                moinp = mm.group(1).lower()
            i += 1
            continue
        i += 1
    # 3. coordinate block
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
                try:
                    sym = parts[0]
                    x, y, z = (float(parts[1]), float(parts[2]), float(parts[3]))
                except ValueError:
                    continue
                newgto = None
                if "newgto" in s.lower():
                    m = re.search(r'NewGTO\s+"([^"]+)"', s, re.IGNORECASE)
                    if m:
                        newgto = m.group(1)
                coords.append((sym, (x, y, z), newgto))
            continue
        if s.startswith("*"):
            in_coords = True
            toks = s[1:].strip().split()
            if len(toks) >= 3 and toks[0].lower() == "xyz":
                try:
                    charge = int(toks[1])
                    mult = int(toks[2])
                except (ValueError, IndexError):
                    pass
            continue
    return {
        "tokens": frozenset(tokens),
        "maxcore": maxcore,
        "blocks": blocks,
        "moinp": moinp,
        "charge": charge,
        "multiplicity": mult,
        "coords": tuple(coords),
    }


# _consume_block / _consume_container are imported from the CO2 parity
# module (identical parsing contract); the OCCUPIER %scf normaliser is
# local because it must keep APMethod (see normalise_inp).
from delfin.tools.co2_coordinator_opi import (  # noqa: E402
    _consume_block,
    _consume_container,
)


def _normalise_scf_block(body: Sequence[str]) -> dict:
    """Canonical form of a ``%scf`` block, extended with APMethod.

    Deliberately NOT the CO2 module's helper: the OCCUPIER broken-sym
    string carries ``APMethod N`` (occupier.py:1413), which the CO2 scan
    writer never emits. Dropping it here would make parity silently
    accept a missing APMethod.
    """
    out: dict[str, Any] = {}
    joined = " ".join(body).lower()
    m = re.search(r"\bmaxiter\s+(\d+)", joined)
    if m:
        out["maxiter"] = int(m.group(1))
    m = re.search(r"\bbrokensym\s+([0-9,\s]+)", joined)
    if m:
        out["brokensym"] = tuple(int(v) for v in re.findall(r"\d+", m.group(1)))
    m = re.search(r"\bapmethod\s+(\d+)", joined)
    if m:
        out["apmethod"] = int(m.group(1))
    return out


def _normalise_qmmm_block(body: Sequence[str]) -> dict:
    joined = " ".join(body)
    out: dict[str, Any] = {}
    m = re.search(r"qmatoms\s*\{\s*([^}]*)\}", joined, re.IGNORECASE)
    if not m:
        m = re.search(r"QMATOMS\s*\{\s*(.*?)\}", joined)
    if m:
        inner = m.group(1)
        atoms: set[int] = set()
        for tok in re.split(r"[,\s:]+", inner.strip()):
            if not tok:
                continue
            if tok.isdigit() or (tok.startswith("-") and tok[1:].isdigit()):
                atoms.add(int(tok))
        out["qmatoms"] = frozenset(atoms)
    return out


def inputs_equivalent(a: dict, b: dict) -> List[str]:
    """Structural diff of two normalised inputs; empty list = equivalent."""
    diffs: List[str] = []
    if set(a["tokens"]) != set(b["tokens"]):
        diffs.append(
            "tokens differ: "
            f"only-classic={sorted(set(a['tokens']) - set(b['tokens']))} "
            f"only-opi={sorted(set(b['tokens']) - set(a['tokens']))}"
        )
    if a["maxcore"] != b["maxcore"]:
        diffs.append(f"maxcore: {a['maxcore']!r} != {b['maxcore']!r}")
    if a["moinp"] != b["moinp"]:
        diffs.append(f"moinp: {a['moinp']!r} != {b['moinp']!r}")
    names = set(a["blocks"]) | set(b["blocks"])
    for name in sorted(names):
        av = a["blocks"].get(name)
        bv = b["blocks"].get(name)
        if av is None or bv is None:
            diffs.append(f"block {name!r}: only in one side")
            continue
        keys = set(av) | set(bv)
        for k in sorted(keys):
            if av.get(k) != bv.get(k):
                diffs.append(
                    f"block {name!r} option {k!r}: {av.get(k)!r} != {bv.get(k)!r}"
                )
    if a["charge"] != b["charge"]:
        diffs.append(f"charge: {a['charge']!r} != {b['charge']!r}")
    if a["multiplicity"] != b["multiplicity"]:
        diffs.append(f"multiplicity: {a['multiplicity']!r} != {b['multiplicity']!r}")
    ca, cb = a["coords"], b["coords"]
    if len(ca) != len(cb):
        diffs.append(f"coord count differs: {len(ca)} vs {len(cb)}")
        return diffs
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
    return diffs


# ---------------------------------------------------------------------------
# OPI import (lazy; the module imports without OPI present)
# ---------------------------------------------------------------------------

def _import_opi():
    try:
        from opi.core import Calculator
        from opi.input import Input
        from opi.input.arbitrary_string import ArbitraryString
        from opi.input.blocks.block_freq import BlockFreq
        from opi.input.blocks.block_qmmm import BlockQmmm
        from opi.input.blocks.block_scf import BlockScf
        from opi.input.blocks.util import NumList
        from opi.input.structures import Atom, Structure
    except ImportError as exc:  # pragma: no cover - environment guard
        raise OccupierOpiError(
            "OPI is required to write OCCUPIER inputs on OPI's building "
            f"blocks, but 'import opi' failed: {exc}"
        ) from exc
    return {
        "Calculator": Calculator,
        "Input": Input,
        "ArbitraryString": ArbitraryString,
        "BlockFreq": BlockFreq,
        "BlockQmmm": BlockQmmm,
        "BlockScf": BlockScf,
        "NumList": NumList,
        "Atom": Atom,
        "Structure": Structure,
    }


# ---------------------------------------------------------------------------
# broken-sym string -> OPI BlockScf (the classic writer appends the raw
# CONTROL string verbatim; occupier.py:1072-1073 / :553-554)
# ---------------------------------------------------------------------------

def _parse_scf_lines(broken_sym: str) -> List[str]:
    """Split an ``%scf ... end`` string into its inner option lines."""
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

    Recognised options: ``MaxIter N`` → maxiter, ``BrokenSym a,b`` →
    brokensym, ``APMethod N`` → apmethod (ORCA 6.1 manual,
    contents/spectroscopyproperties/magnx.txt: "The AP method is
    requested via the tag APMethod in the %scf block").  An option OPI's
    BlockScf cannot represent is an error — never silently dropped.
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
        m = re.match(r"apmethod\s+(\d+)", low)
        if m:
            fields["apmethod"] = int(m.group(1))
            continue
        raise ValueError(
            "OPI BlockScf cannot represent broken-sym line "
            f"{ln!r}; refusing to silently drop an %scf option "
            "(only MaxIter, BrokenSym and APMethod are mapped)."
        )
    return fields


def _scf_for_fob(broken_sym: str):
    """The ``%scf`` block(s) the classic writer emits for one FoB.

    The classic writer writes TWO %scf occurrences for a BS FoB:
    ``%scf maxiter N end`` from CONTROL (occupier.py:1068/:549) and the
    raw broken-sym string (occupier.py:1072-1073/:553-554).  OPI keeps
    one block per type, so both collapse into ONE BlockScf with the same
    option set — ORCA's parser reads both spellings identically and the
    normaliser reduces them to the same structure.
    """
    opi = _import_opi()
    BlockScf = opi["BlockScf"]
    NumList = opi["NumList"]
    fields: Dict[str, Any] = {}
    # the raw broken-sym string first; the classic writer appends it
    # verbatim, so an option there wins in the classic file too
    if broken_sym and broken_sym.strip():
        fields = _scf_fields_for_lines(_parse_scf_lines(broken_sym))
    if "brokensym" in fields:
        fields["brokensym"] = NumList(fields["brokensym"])
    if not fields:
        return None
    return BlockScf(**fields)


def _freq_block(config: Mapping[str, Any]):
    """The ``%freq Temp`` block when OCCUPIER_compare resolves to G.

    Classic: FREQ keyword + ``%freq\\n  Temp {temperature}\\nend`` from
    ``delfin.xyz_io._build_freq_block`` (xyz_io.py:497-503), appended when
    freq_flag is set (occupier.py:1076-1079 / :556-560).
    """
    from delfin.common.control_validator import resolve_occupier_compare

    if resolve_occupier_compare(dict(config)) != "G":
        return None
    opi = _import_opi()
    temperature = config.get("temperature", "298.15")
    return opi["BlockFreq"](temp=float(temperature))


def _qmmm_block(qmmm_range):
    """The ``%qmmm`` block for a QM/XTB split, or None.

    Classic: ``delfin.xyz_io.build_qmmm_block`` (xyz_io.py:310-327) writes
    ``%QMMM / QMAtoms {s:e} end / Charge_Medium 0 / Mult_Medium 1 / end``.
    """
    if not qmmm_range:
        return None
    start, end = qmmm_range
    if end < start:
        return None
    opi = _import_opi()
    BlockQmmm = opi["BlockQmmm"]
    return BlockQmmm(
        qmatoms=list(range(start, end + 1)),
        charge_medium=0,
        mult_medium=1,
    )


# ---------------------------------------------------------------------------
# The '!' line — mirrored token assembly (occupier.py:508-540)
# ---------------------------------------------------------------------------

def _bang_tokens(
    config: Mapping[str, Any],
    *,
    qmmm_range: Optional[Tuple[int, int]],
    broken_sym: str,
    pass_wavefunction: str,
) -> List[str]:
    """The '!' line tokens, in the classic writer's assembly order.

    Reproduces delfin/occupier.py:500-545 (and the job-level variant at
    :1014-1064, same rules): 3c composites carry their own basis and
    dispersion (ORCA 6.1.1 manual section 3.6.3, quoted at
    delfin/occupier.py:502-507); the generic branch emits
    functional, rel token, main basis, dispersion, RI, aux-JK, solvation,
    geom_opt, FREQ, initial guess.  ``MORead`` is appended when the
    broken-sym string mentions moinp (occupier.py:543-544).
    """
    from delfin.common.solvation import build_solvation_keyword
    from delfin.common.control_validator import resolve_occupier_compare
    from delfin.occupier import _is_composite_3c
    from delfin.utils import resolve_level_of_theory

    found_metals = list(config.get("_found_metals") or [])
    main_basis, _metal_basis, rel_token, aux_jk_token = resolve_level_of_theory(
        found_metals, dict(config),
        config.get("main_basisset"), config.get("metal_basisset"),
    )
    implicit = build_solvation_keyword(
        config.get("implicit_solvation_model"), config.get("solvent")
    )
    initial_guess = (str(config.get("initial_guess", "")).split() or [""])[0]
    freq_flag = "FREQ" if resolve_occupier_compare(dict(config)) == "G" else ""
    functional_str = str(config["functional"])

    tokens: List[str] = []
    if qmmm_range:
        tokens.append(str(config.get("qmmm_option", "QM/XTB")).strip())

    if _is_composite_3c(functional_str):
        tokens.append(functional_str)
        for extra in (
            implicit,
            str(config.get("geom_opt_OCCUPIER", config.get("geom_opt", ""))),
            str(config.get("ri_jkx", "")),
        ):
            if extra:
                tokens.append(extra)
    else:
        tokens.extend([
            functional_str,
            rel_token,
            main_basis,
            str(config.get("disp_corr", "")),
            str(config.get("ri_jkx", "")),
            aux_jk_token,
            implicit,
            str(config.get("geom_opt_OCCUPIER", config.get("geom_opt", ""))),
        ])
    if freq_flag:
        tokens.append(freq_flag)
    tokens.append(initial_guess)
    if broken_sym and "moinp" in broken_sym.lower():
        tokens.append("MORead")
    if str(pass_wavefunction or "no").strip().lower() in (
        "yes", "true", "1", "on", "y"
    ):
        tokens.append("MORead")
    return [t for t in tokens if t]


# ---------------------------------------------------------------------------
# Per-atom NewGTO tags — the classic geometry loop (occupier.py:589-604)
# ---------------------------------------------------------------------------

def _metal_and_first_sphere(
    atoms: Sequence[Tuple[str, Tuple[float, float, float]]],
    found_metals: Sequence[str],
    metal_basis: Optional[str],
    enable_first: bool,
    sphere_scale: float,
    radii_map,
):
    """Reproduce the classic metal/first-sphere tagging decision.

    Mirrors occupier.py:572-587: metals found by symbol (capitalised),
    first-sphere by the same distance rule as the classic helper — the
    classic's own ``_first_coordination_sphere_indices`` is imported and
    used so the rule cannot drift.
    """
    from delfin.occupier import _first_coordination_sphere_indices

    metal_syms = {m.strip().capitalize() for m in (found_metals or [])}
    metal_indices = [i for i, (el, _c) in enumerate(atoms)
                     if el.capitalize() in metal_syms]
    atom_dicts = [
        {"line_idx": i, "raw": el, "elem": el, "x": c[0], "y": c[1], "z": c[2]}
        for i, (el, c) in enumerate(atoms)
    ]
    first_sphere: set[int] = set()
    if enable_first and metal_indices and metal_basis:
        first_sphere = _first_coordination_sphere_indices(
            atom_dicts, metal_indices, sphere_scale, radii_map
        )
    newgto: Dict[int, str] = {}
    if metal_basis:
        for i in metal_indices:
            newgto[i] = metal_basis
        if enable_first:
            for i in sorted(first_sphere):
                newgto[i] = metal_basis
    return newgto


# ---------------------------------------------------------------------------
# Core OPI write
# ---------------------------------------------------------------------------

def write_occupier_fob_opi_input(
    work_dir,
    *,
    config: Mapping[str, Any],
    charge: int,
    multiplicity: int,
    PAL: int,
    maxcore: int,
    solvent: Optional[str],
    found_metals: Sequence[str],
    broken_sym: str = "",
    from_index: int = 1,
    base_name: str = "input",
    atoms: Optional[Sequence[Tuple[str, Tuple[float, float, float]]]] = None,
    xyz_path=None,
    qmmm_range: Optional[Tuple[int, int]] = None,
    pass_wavefunction: str = "",
    geom_path_ready=None,
) -> Path:
    """Write one OCCUPIER FoB ``.inp`` with OPI's input classes.

    OPI equivalent of ``delfin.occupier.read_and_modify_file_OCCUPIER``
    (delfin/occupier.py:414-612) for one fragment-of-basis job.  Only the
    input is written — the classic flow additionally launches ORCA, which
    this parity module deliberately does not do.
    """
    if atoms is None and xyz_path is None:
        raise OccupierOpiError(
            "write_occupier_fob_opi_input needs either atoms= or xyz_path="
        )
    opi = _import_opi()
    from delfin.utils import resolve_level_of_theory

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # ---- geometry source: classic reads input{n}.xyz / input.xyz ------
    # (occupier.py:425: "input.xyz" if from_index == 1 else f"input{from_index}.xyz")
    if atoms is None:
        xyz_file = "input.xyz" if from_index == 1 else f"input{from_index}.xyz"
        xyz_file = Path(xyz_path) if xyz_path else work_dir / xyz_file
        lines = xyz_file.read_text().splitlines()
        body = lines[2:]
        atoms = []
        for line in body:
            ls = line.strip()
            if not ls or ls == "*":
                break
            parts = ls.split()
            if len(parts) < 4:
                continue
            atoms.append((parts[0], (float(parts[1]), float(parts[2]), float(parts[3]))))
    if not atoms:
        raise OccupierOpiError("no atoms to write")

    # ---- level of theory / solvation, the classic writer's own helpers -
    cfg = dict(config)
    cfg["_found_metals"] = list(found_metals)
    main_basis, metal_basis, _rel, _aux = resolve_level_of_theory(
        list(found_metals), cfg,
        cfg.get("main_basisset"), cfg.get("metal_basisset"),
    )

    # 3c composite carries its own basis: no per-atom NewGTO override
    # (occupier.py:566-571, quoting ORCA 6.1.1 manual section 3.6.3)
    from delfin.occupier import _is_composite_3c
    if _is_composite_3c(str(cfg["functional"])):
        metal_basis = None

    enable_first = str(cfg.get(
        "first_coordination_sphere_metal_basisset", "no")).lower() in (
        "yes", "true", "1", "on")
    sphere_scale_raw = str(cfg.get("first_coordination_sphere_scale", "")).strip()
    if enable_first and sphere_scale_raw:
        sphere_scale = float(sphere_scale_raw)
        radii_all = None
    elif enable_first:
        from delfin.occupier import load_covalent_radii
        radii_all = load_covalent_radii(
            source=str(cfg.get("covalent_radii_source", "pyykko2009")))
        sphere_scale = 1.20
    else:
        radii_all = None
        sphere_scale = float(sphere_scale_raw or 1.20)

    newgto = _metal_and_first_sphere(
        list(atoms), found_metals, metal_basis, enable_first,
        sphere_scale, radii_all,
    )

    tokens = _bang_tokens(
        cfg,
        qmmm_range=qmmm_range,
        broken_sym=broken_sym,
        pass_wavefunction=pass_wavefunction,
    )

    # ---- assemble the OPI Input --------------------------------------
    Input = opi["Input"]
    Atom = opi["Atom"]
    Structure = opi["Structure"]
    Calculator = opi["Calculator"]

    inp = Input()
    inp.add_simple_keywords(*[str(t) for t in tokens])
    inp.memory = int(maxcore)
    inp.ncores = int(PAL)

    scf_block = _scf_for_fob(broken_sym)
    # %scf maxiter from CONTROL merges into the same block (classic writes
    # it as its own %scf, occupier.py:1068/:549; one block per type in OPI)
    maxiter = cfg.get("maxiter_occupier")
    if maxiter is not None:
        if scf_block is None:
            scf_block = _scf_for_fob("")
        if scf_block is None:
            scf_block = opi["BlockScf"](maxiter=int(maxiter))
        elif getattr(scf_block, "maxiter", None) is None:
            scf_block.maxiter = int(maxiter)
    if scf_block is not None:
        inp.add_blocks(scf_block)

    freq_block = _freq_block(cfg)
    if freq_block is not None:
        inp.add_blocks(freq_block)

    qmmm_block = _qmmm_block(qmmm_range)
    if qmmm_block is not None:
        inp.add_blocks(qmmm_block)

    # %moinp: the classic job loop emits it as a raw line
    # '%moinp "<gbw>"' (occupier.py:1402) when pass_wavefunction is on and
    # the source's .gbw exists. OPI's Input.moinp renders exactly
    # '%moinp "<path>"' (opi/input/core.py:~216). When the file does not
    # exist (parity tests), an arbitrary string is used instead — the
    # classic writer writes the same line unconditionally on an existing
    # gbw, and OPI's Input.moinp requires a real subpath.
    if str(pass_wavefunction or "no").strip().lower() in (
        "yes", "true", "1", "on", "y"
    ):
        gbw_candidate = (
            "input.gbw" if from_index == 1 else f"input{from_index}.gbw"
        )
        inp.add_arbitrary_string(
            f'%moinp "{gbw_candidate}"'
        )

    atoms_opi = []
    for idx, (el, (x, y, z)) in enumerate(atoms):
        append = None
        basis = newgto.get(idx)
        if basis:
            append = f'NewGTO "{basis}" end'
        atoms_opi.append(
            Atom(str(el).capitalize(),
                 coordinates=(round(float(x), 8), round(float(y), 8),
                              round(float(z), 8)),
                 append_str=append)
        )
    structure = Structure(atoms_opi, charge=int(charge),
                          multiplicity=int(multiplicity))

    base = base_name if from_index == 1 else f"{base_name}{from_index}"
    calc = Calculator(base, working_dir=work_dir, version_check=False)
    calc.json_via_input = False  # parity: classic never writes %output
    calc.input = inp
    calc.structure = structure
    calc.write_input(force=True)
    return work_dir / f"{base}.inp"
