"""OPI-backed ORCA input writer for the ESD workflow.

OPI equivalent of ``delfin.esd_input_generator``: the level-of-theory,
keyword-assembly, solvation, TDDFT-block, ESD-block and per-atom-NewGTO
logic is reused BY IMPORT (the same helpers the classic writer itself
calls) and the same ORCA ``.inp`` files are written with OPI's input
classes.  The classic writer stays the default; this module is the
parity-proof reference for its inputs, reachable through the step
adapter ``delfin.tools.adapters.esd_opi``.

The functions write the input only — they never launch ORCA.  OPI is
imported lazily so the module imports without OPI present (py3.11 CI).

Block coverage and its limits (all deliberate, each raising rather than
dropping what cannot be represented):
* ``%tddft`` on ``BlockTddft``; options OPI's pydantic fields cannot
  hold (``FollowIRoot``, job-scoped ``SRoot``/``TRoot``/``TRootSSL``,
  ``NACME``/``ETF``) go through ``modify_option`` onto the block's
  ``_arbitrary`` store, which OPI renders as ``    key value`` lines.
* ``%esd`` has no OPI block class, so it is written as an
  ``ArbitraryString`` with the classic writer's exact body, obtained
  from its own code paths where they exist.
* ``%output``'s ``print[p_*]`` entries have no OPI field either, so the
  output blocks ride along as arbitrary strings too.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)


class EsdOpiError(RuntimeError):
    """Raised when an ESD input cannot be represented on OPI's blocks."""


# ---------------------------------------------------------------------------
# Normalisation / structural comparison (ESD-specific extension of the
# shared parity normaliser: the %tddft / %esd / %geom / %basis / %output
# blocks the CO2/OCCUPIER writers never emit)
# ---------------------------------------------------------------------------

def _normalise_tddft_block(body: Sequence[str]) -> dict:
    joined = " ".join(body)
    out: dict[str, Any] = {}

    def _grab(key: str) -> Optional[str]:
        m = re.search(rf"\b{key}\s+(\S+)", joined, re.IGNORECASE)
        return m.group(1) if m else None

    for key in ("nroots", "maxdim", "maxiter", "iroot", "sroot",
                "troot", "trootssl"):
        val = _grab(key)
        if val is not None:
            try:
                out[key] = int(val)
            except ValueError:
                out[key] = val
    for key in ("tda", "triplets", "dosoc", "followiroot", "nacme", "etf"):
        val = _grab(key)
        if val is not None:
            out[key] = val.lower() in ("true", "1", "yes")
    val = _grab("irootmult")
    if val is not None:
        out["irootmult"] = val.lower()
    return out


def _normalise_esd_block(body: Sequence[str]) -> dict:
    joined = " ".join(body)
    out: dict[str, Any] = {}
    for key in ("iscishess", "iscfshess", "gshessian", "eshessian",
                "tshessian"):
        m = re.search(rf'\b{key}\s+"([^"]+)"', joined, re.IGNORECASE)
        if m:
            out[key] = m.group(1).lower()
    for key in ("usej", "doht"):
        m = re.search(rf"\b{key}\s+(\S+)", joined, re.IGNORECASE)
        if m:
            out[key] = m.group(1).lower() in ("true", "1", "yes")
    m = re.search(r"\blines\s+(\S+)", joined, re.IGNORECASE)
    if m:
        out["lines"] = m.group(1).lower()
    for key in ("linew", "inlinew", "temp", "dele", "npoints", "maxtime"):
        m = re.search(rf"\b{key}\s+([0-9eE.+-]+)", joined)
        if m:
            try:
                out[key] = float(m.group(1))
            except ValueError:
                out[key] = m.group(1)
    return out


def _normalise_geom_block(body: Sequence[str]) -> dict:
    joined = " ".join(body)
    out: dict[str, Any] = {}
    m = re.search(r"\bmaxiter\s+(\d+)", joined, re.IGNORECASE)
    if m:
        out["maxiter"] = int(m.group(1))
    return out


def _normalise_basis_block(body: Sequence[str]) -> dict:
    out: dict[str, Any] = {}
    for s in body:
        m = re.search(r'NewGTO\s+(\S+)\s+"([^"]+)"', s, re.IGNORECASE)
        if m:
            out[m.group(1).capitalize()] = m.group(2)
    return out


def _normalise_output_block(body: Sequence[str]) -> dict:
    joined = " ".join(body).lower()
    out: dict[str, Any] = {}
    for key in ("print[p_mos]", "print[p_basis]", "print[p_lpop]",
                "print[p_mulliken]"):
        m = re.search(rf"{re.escape(key)}\s+(\d+)", joined)
        if m:
            out[key] = int(m.group(1))
    return out


def normalise_esd_inp(text: str) -> dict:
    """Parse an ESD ORCA ``.inp`` into a canonical structure.

    Same contract as ``delfin.tools.occupier_opi.normalise_inp`` (tokens /
    maxcore / blocks / moinp / charge / multiplicity / coords), with the
    ESD-specific block normalisers added.  The simple-keyword, %pal,
    %maxcore, %moinp, %scf and coordinate parsing is imported from the
    shared parity normaliser so the comparison cannot drift between the
    OPI writers.
    """
    from delfin.tools.occupier_opi import normalise_inp

    norm = normalise_inp(text)
    lines = text.splitlines()
    blocks: Dict[str, dict] = dict(norm["blocks"])
    i = 0
    n = len(lines)
    while i < n:
        s = lines[i].strip()
        low = s.lower()
        if low.startswith("%tddft") or low.startswith("%esd") \
                or low.startswith("%geom") or low.startswith("%basis") \
                or low.startswith("%output"):
            body: List[str] = []
            i += 1
            while i < n:
                t = lines[i].strip().lower()
                if t == "end":
                    i += 1
                    break
                if t.startswith("%"):
                    break  # unterminated block: leave for next iteration
                body.append(lines[i].strip())
                i += 1
            if low.startswith("%tddft"):
                blocks["tddft"] = _normalise_tddft_block(body)
            elif low.startswith("%esd"):
                blocks["esd"] = _normalise_esd_block(body)
            elif low.startswith("%geom"):
                blocks["geom"] = _normalise_geom_block(body)
            elif low.startswith("%basis"):
                blocks["basis"] = _normalise_basis_block(body)
            else:
                blocks["output"] = _normalise_output_block(body)
            continue
        i += 1
    norm["blocks"] = blocks
    return norm


def esd_inputs_equivalent(a: dict, b: dict) -> List[str]:
    """Structural diff of two normalised ESD inputs; empty = equivalent."""
    from delfin.tools.occupier_opi import inputs_equivalent

    diffs: List[str] = []
    for sa, sb in zip(a, b):
        sa = {**sa, "moinp": sa.get("moinp")}
        sb = {**sb, "moinp": sb.get("moinp")}
        diffs.extend(inputs_equivalent(sa, sb))
    if len(a) != len(b):
        diffs.append(f"job-section count differs: {len(a)} vs {len(b)}")
    return diffs


def split_jobs(text: str) -> List[str]:
    """Split a multi-job input on ``$new_job`` into job sections.

    Comment lines are dropped on every section so comment-only
    differences cannot mask real ones; empty sections vanish.
    """
    sections: List[str] = []
    current: List[str] = []
    for line in text.splitlines():
        if line.strip() == "$new_job":
            sections.append("\n".join(current))
            current = []
        else:
            current.append(line)
    sections.append("\n".join(current))
    out: List[str] = []
    for sec in sections:
        body = [ln for ln in sec.splitlines()
                if not ln.strip().startswith("#") and ln.strip()]
        if body:
            out.append("\n".join(body))
    return out


# ---------------------------------------------------------------------------
# OPI import (lazy)
# ---------------------------------------------------------------------------

def _import_opi():
    try:
        from opi.core import Calculator
        from opi.input import Input
        from opi.input.arbitrary_string import ArbitraryString
        from opi.input.blocks.block_tddft import BlockTddft
        from opi.input.blocks.block_scf import BlockScf
        from opi.input.blocks.block_geom import BlockGeom
        from opi.input.blocks.block_basis import BlockBasis, NewBasis
        from opi.input.simple_keywords.base import SimpleKeyword
        from opi.input.blocks.util import NumList
        from opi.input.structures import Atom, Structure
    except ImportError as exc:  # pragma: no cover - environment guard
        raise EsdOpiError(
            "OPI is required to write ESD inputs on OPI's building blocks, "
            f"but 'import opi' failed: {exc}"
        ) from exc
    return {
        "Calculator": Calculator,
        "Input": Input,
        "ArbitraryString": ArbitraryString,
        "BlockTddft": BlockTddft,
        "BlockScf": BlockScf,
        "BlockGeom": BlockGeom,
        "BlockBasis": BlockBasis,
        "NewBasis": NewBasis,
        "SimpleKeyword": SimpleKeyword,
        "NumList": NumList,
        "Atom": Atom,
        "Structure": Structure,
    }


# ---------------------------------------------------------------------------
# Shared assembly helpers
# ---------------------------------------------------------------------------

def _kw(opi, name: str):
    return opi["SimpleKeyword"](str(name))


def _new_basis_block(opi, metals: Sequence[str],
                     metal_basis: Optional[str]):
    """One BlockBasis per element with a NewGTO override.

    OPI's BlockBasis holds a single NewBasis per field, so several
    elements must go into several blocks (the overwrite bug documented
    in delfin/tools/orca_basis_solv_opi.py:411-416).
    """
    if not metals or not metal_basis:
        return []
    unique = sorted({str(m).strip().capitalize() for m in metals if str(m).strip()})
    blocks = []
    for el in unique:
        blocks.append(opi["BlockBasis"](
            newgto=opi["NewBasis"](element=el, basis=_kw(opi, metal_basis))))
    return blocks


def _inline_atoms(
    opi, esd_dir: Path, xyz_name: str, metals: Sequence[str],
    metal_basis: Optional[str], config: Mapping[str, Any],
) -> List[Any]:
    """Atoms for an inline-geometry job, per-atom NewGTO appended.

    The classic writer's ``_apply_esd_newgto`` supplies the tagging
    decision (metals always, first coordination sphere optionally); this
    helper only turns its tagged lines into OPI atoms.
    """
    from delfin.esd_input_generator import _apply_esd_newgto

    if str(xyz_name) in ("initial.xyz", "start.txt"):
        xyz_path = Path(xyz_name)  # classic reads these from the cwd
    else:
        xyz_path = esd_dir / xyz_name
    with open(xyz_path, "r", encoding="utf-8") as fh:
        all_lines = fh.readlines()
    coord_lines = all_lines[2:]
    coord_lines = _apply_esd_newgto(
        coord_lines, found_metals=list(metals),
        metal_basisset=metal_basis, config=dict(config),
    )
    atoms = []
    for line in coord_lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        sym = parts[0]
        x, y, z = (float(parts[1]), float(parts[2]), float(parts[3]))
        append = None
        m = re.search(r'NewGTO\s+"([^"]+)"', line, re.IGNORECASE)
        if m:
            append = f'NewGTO "{m.group(1)}" end'
        atoms.append(opi["Atom"](
            sym.capitalize(),
            coordinates=(round(x, 8), round(y, 8), round(z, 8)),
            append_str=append,
        ))
    return atoms


def _tddft_block_opi(
    opi,
    config: Mapping[str, Any],
    *,
    nroots: Optional[int] = None,
    iroot: Optional[int] = None,
    irootmult: Optional[str] = None,
    triplets: bool = False,
    follow: bool = False,
    dosoc: Any = ...,
    job_lines: Sequence[Tuple[str, Any]] = (),
):
    """OPI BlockTddft mirroring ``tddft_settings.tddft_block``.

    Reads the same CONTROL settings via the classic helper, then maps
    them onto OPI's block.  Options OPI's pydantic fields cannot hold
    go through ``modify_option`` (rendered as bare ``key value`` lines
    inside the block); an unsupported *typed* field raises rather than
    being dropped.
    """
    from delfin.common import tddft_settings

    settings = tddft_settings.read_settings(dict(config))
    block = opi["BlockTddft"]()
    block.nroots = int(settings.nroots if nroots is None else nroots)
    block.maxdim = int(settings.maxdim)
    block.tda = bool(settings.tda)
    if settings.maxiter is not None:
        block.maxiter = int(settings.maxiter)
    if triplets:
        block.triplets = True
    if iroot is not None:
        block.iroot = int(iroot)
    if irootmult:
        block.irootmult = str(irootmult)
    if follow and iroot is not None and settings.followiroot:
        block.modify_option("followiroot", "True")
    soc = settings.soc if dosoc is ... else dosoc
    if soc is not None:
        block.modify_option("dosoc", str(bool(soc)))
    for key, value in job_lines:
        block.modify_option(key, str(value))
    return block


def _esd_block_lines(
    opi,
    config: Mapping[str, Any],
    header_lines: Sequence[Tuple[str, Any]],
) -> Any:
    """The ``%esd ... end`` block as an ArbitraryString.

    OPI has no %esd block class.  The body is built from the same
    CONTROL values the classic writer reads; header_lines is the
    key/value list the caller prepares (Hessian names, USEJ, DOHT, ...).
    """
    from delfin.esd_input_generator import (
        _esd_grid_lines,
        _resolve_temperature_K,
    )
    from delfin.common.orca_input import keyword_family  # noqa: F401

    doht = str(config.get("DOHT", "TRUE")).upper()
    lines = str(config.get("ESD_LINES", "LORENTZ")).strip().upper() or "LORENTZ"
    linew = str(config.get("ESD_LINEW", 50)).strip()
    inlinew = str(config.get("ESD_INLINEW", 250)).strip()
    temperature = _resolve_temperature_K(config, default=298.15)

    body = ["%ESD"]
    for key, value in header_lines:
        body.append(f'  {key:<16}{value}')
    body.append("  USEJ            TRUE")
    body.append(f"  DOHT            {doht}")
    body.append(f"  LINES           {lines}")
    body.append(f"  LINEW           {linew}")
    body.append(f"  INLINEW         {inlinew}")
    body.extend(_esd_grid_lines(dict(config)))
    body.append(f"  TEMP            {temperature}")
    body.append("END")
    return opi["ArbitraryString"]("\n".join(body))


def _output_block_strings(
    opi, config: Mapping[str, Any], *, allow: bool = True,
) -> List[Any]:
    """%output blocks as ArbitraryStrings (no OPI print[] field)."""
    from delfin.common.orca_blocks import collect_output_blocks

    out = []
    for text in collect_output_blocks(dict(config), allow=allow):
        out.append(opi["ArbitraryString"](text.rstrip("\n")))
    return out


def _scf_block_opi(
    opi, config: Mapping[str, Any], esd_dir: Path, state_upper: str,
    *, deltascf: bool,
):
    """The deltaSCF %scf block on OPI's BlockScf (None when not needed).

    alphaconf/betaconf come from the classic writer's own derivation
    (``_parse_tddft_and_derive_deltascf`` on S0.out, else its fallback),
    so both writers cannot drift here.
    """
    from delfin.esd_input_generator import _parse_tddft_and_derive_deltascf

    if not deltascf:
        return None
    domom = str(config.get("deltaSCF_DOMOM", "true")).lower() in (
        "true", "1", "yes", "on")
    pmom = str(config.get("deltaSCF_PMOM", "true")).lower() in (
        "true", "1", "yes", "on")
    keepinitialref = str(
        config.get("deltaSCF_keepinitialref", "true")).lower() in (
        "true", "1", "yes", "on")
    soscfhessup = str(config.get("deltaSCF_SOSCFHESSUP", "LSR1"))
    maxiter = int(config.get("deltaSCF_maxiter", 1000))
    conv = config.get("deltaSCF_SOSCFConvFactor", 500)
    maxstep = config.get("deltaSCF_SOSCFMaxStep", 0.1)

    # alphaconf/betaconf: classic derivation, fallback identical to it
    alphaconf, betaconf = None, None
    s0_out_path = esd_dir / "S0.out"
    if s0_out_path.exists():
        result = _parse_tddft_and_derive_deltascf(s0_out_path, state_upper)
        if result:
            alphaconf, betaconf = result
    if alphaconf is None:
        if state_upper.startswith("T") and state_upper != "T1":
            alphaconf, betaconf = "0,1", "0,1"
        else:
            alphaconf, betaconf = "0,1", "0"

    def _numlist(text: str):
        nums = [int(p) for p in str(text).replace(" ", "").split(",") if p != ""]
        return opi["NumList"](nums)

    return opi["BlockScf"](
        domom=domom,
        pmom=pmom,
        keepinitialref=keepinitialref,
        soscfhessup=str(soscfhessup).lower(),
        maxiter=maxiter,
        soscfconvfactor=conv,
        soscfmaxstep=maxstep,
        alphaconf=_numlist(alphaconf),
        betaconf=_numlist(betaconf),
    )


def _rate_tokens(tokens: Sequence[str]) -> List[str]:
    """The classic rate job's SCF convergence keyword (TightSCF)."""
    from delfin.common.orca_input import keyword_family

    out = [str(t) for t in tokens]
    if any(keyword_family(t) == "scf convergence" for t in out):
        return out
    out.append("TightSCF")
    return out


def _job_in_opi(
    opi,
    *,
    out_dir: Path,
    base: str,
    tokens: Sequence[str],
    maxcore: int,
    pal: int,
    moinp: Optional[str],
    blocks: Sequence[Any],
    arbitrary: Sequence[str] = (),
    atoms: Optional[List[Any]] = None,
    xyzfile: Optional[Tuple[int, int, str]] = None,
    charge: int,
    multiplicity: int,
) -> Path:
    """Render one ORCA job with OPI's Calculator (input + structure)."""
    inp = opi["Input"]()
    inp.add_simple_keywords(*[str(t) for t in tokens])
    inp.memory = int(maxcore)
    inp.ncores = int(pal)
    if moinp:
        inp.add_arbitrary_string(f'%moinp "{moinp}"')
    for text in arbitrary:
        inp.add_arbitrary_string(text)
    for blk in blocks:
        if blk is None:
            continue
        if hasattr(blk, "aftercoord"):  # an OPI Block
            inp.add_blocks(blk)
        else:  # an ArbitraryString (or anything else rendered verbatim)
            inp.add_arbitrary_string(str(blk))
    if xyzfile is not None:
        ch, mp, name = xyzfile
        inp.add_arbitrary_string(f"* xyzfile {ch} {mp} {name}")

    structure = None
    if atoms is not None:
        structure = opi["Structure"](
            atoms, charge=int(charge), multiplicity=int(multiplicity))

    calc = opi["Calculator"](base, working_dir=out_dir, version_check=False)
    calc.json_via_input = False
    calc.input = inp
    calc.structure = structure
    calc.write_input(force=True)
    return out_dir / f"{base}.inp"


def _level_tokens(
    metals: Sequence[str], config: Mapping[str, Any],
    main_basisset: Optional[str], metal_basisset: Optional[str],
) -> Tuple[str, Optional[str], str, str]:
    from delfin.utils import resolve_level_of_theory

    return resolve_level_of_theory(
        list(metals), dict(config), main_basisset, metal_basisset)


# ---------------------------------------------------------------------------
# State jobs (TDDFT / deltaSCF / hybrid1)
# ---------------------------------------------------------------------------

def _state_geometry(esd_dir: Path, state_upper: str,
                    mode: str) -> Tuple[Path, int]:
    """Coordinate file + skip for a state, exactly as the classic writer.

    An initial.xyz S0 wrote itself is not a start (``written_by_s0``).
    """
    from delfin.esd_input_generator import written_by_s0

    if state_upper == "S0":
        if Path("initial.xyz").exists() and not written_by_s0(
                Path("initial.xyz")):
            return Path("initial.xyz"), 2
        return Path("start.txt"), 0
    return esd_dir / "S0.xyz", 2


def _s0_xyz(esd_dir: Path) -> Path:
    xyz, _skip = _state_geometry(esd_dir, "S0", "")
    return xyz


def write_esd_state_opi_input(
    esd_dir,
    *,
    out_dir,
    state: str,
    charge: int,
    solvent: Optional[str],
    metals: Sequence[str],
    main_basisset: Optional[str],
    metal_basisset: Optional[str],
    config: Mapping[str, Any],
) -> str:
    """OPI equivalent of ``create_state_input`` for one state."""
    from delfin.esd_input_generator import (
        _build_solvation_keyword,
        _redox_steps_requested,
        _resolve_esd_t1_opt,
    )

    opi = _import_opi()
    esd_dir = Path(esd_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    state_upper = str(state).strip().upper()

    mode = str(config.get("ESD_modus", "TDDFT")).strip().lower()
    if "|" in mode:
        mode = mode.split("|")[0].strip()

    functional = config.get("functional", "PBE0")
    disp_corr = config.get("disp_corr", "D4")
    ri_jkx = config.get("ri_jkx", "RIJCOSX")
    main_basis, metal_basis, rel_token, aux_jk = _level_tokens(
        metals, config, main_basisset, metal_basisset)
    geom_token = str(config.get("geom_opt", "OPT")).strip() or "OPT"
    pal = int(config.get("PAL", 12))
    maxcore = int(config.get("maxcore", 6000))
    freq_on = str(config.get("ESD_frequency", "yes")).strip().lower() in (
        "yes", "true", "1", "on")
    freq_type = str(config.get("freq_type", "FREQ")).strip().upper()
    if freq_type not in ("FREQ", "NUMFREQ"):
        freq_type = "FREQ"
    solv = _build_solvation_keyword(
        config.get("implicit_solvation_model", ""), solvent)
    output_blocks = _output_block_strings(opi, config)

    def level(scf_type: str, *, moread: bool = False,
              with_geom: bool = True) -> List[str]:
        toks = [functional, scf_type, rel_token, main_basis, disp_corr,
                ri_jkx, aux_jk]
        if solv:
            toks.append(solv)
        if with_geom and geom_token:
            toks.append(geom_token)
        if moread:
            toks.append("MOREAD")
        return [t for t in toks if str(t).strip()]

    # ---- TDDFT mode -----------------------------------------------------
    if mode == "tddft":
        xyz_path, skip = _state_geometry(esd_dir, state_upper, mode)
        with open(xyz_path, "r", encoding="utf-8") as fh:
            all_lines = fh.readlines()
        coord_lines = all_lines[skip:]

        if state_upper == "S0":
            toks = level("RKS")
            if freq_on or _redox_steps_requested(config):
                toks.append("FREQ")
            # no TDDFT keyword on the S0 job itself: the classic writer
            # puts the %tddft block only in the check job (S0_TDDFT)
            main = _job_in_opi(
                opi, out_dir=out_dir, base="_S0_main_tmp", tokens=toks,
                maxcore=maxcore, pal=pal, moinp=None,
                blocks=[
                    *_new_basis_block(opi, metals, metal_basis),
                    *output_blocks,
                ],
                atoms=_inline_atoms(opi, esd_dir, _s0_xyz(esd_dir),
                                    metals, metal_basis, config),
                charge=charge, multiplicity=1,
            )
            check_toks = [functional, "RKS", rel_token, main_basis,
                          disp_corr, ri_jkx, aux_jk]
            if solv:
                check_toks.append(solv)
            check = _tddft_check_job(
                opi, esd_dir, out_dir, "S0", charge, metals, metal_basis,
                config, [t for t in check_toks if str(t).strip()],
                maxcore, pal, output_blocks)
            merged = _merge_new_job_files([str(main), check])
            final = out_dir / "S0.inp"
            final.write_text(Path(merged).read_text(encoding="utf-8"),
                             encoding="utf-8")
            Path(merged).unlink()
            return str(final)

        if state_upper == "T1" and _resolve_esd_t1_opt(config) == "tddft":
            toks = level("RKS", moread=True)
            blocks = [
                _tddft_block_opi(opi, config, iroot=1,
                                 irootmult="triplet", triplets=True,
                                 follow=True),
                *output_blocks,
            ]
            atoms = _inline_atoms(opi, esd_dir, "S0.xyz", metals,
                                  metal_basis, config)
            return str(_job_in_opi(
                opi, out_dir=out_dir, base="T1", tokens=toks,
                maxcore=maxcore, pal=pal, moinp="S0.gbw",
                blocks=blocks, atoms=atoms, charge=charge, multiplicity=1,
            ))

        root_num = int(state_upper[1:]) if len(state_upper) > 1 else 1
        is_triplet = state_upper.startswith("T")
        if is_triplet and state_upper == "T1":
            raise EsdOpiError(
                "T1 in TDDFT mode with ESD_T1_opt=uks is a UKS deltaSCF "
                "shape; reached the TDDFT branch unexpectedly")
        toks = level("RKS", moread=True)
        blocks = [
            _tddft_block_opi(
                opi, config, iroot=root_num,
                irootmult="triplet" if is_triplet else "singlet",
                triplets=is_triplet, follow=True),
            *output_blocks,
        ]
        atoms = _inline_atoms(opi, esd_dir, "S0.xyz", metals, metal_basis,
                              config)
        return str(_job_in_opi(
            opi, out_dir=out_dir, base=state_upper, tokens=toks,
            maxcore=maxcore, pal=pal, moinp="S0.gbw",
            blocks=blocks, atoms=atoms, charge=charge, multiplicity=1,
        ))

    # ---- hybrid1 mode ---------------------------------------------------
    if mode == "hybrid1":
        if state_upper == "S0":
            return write_esd_state_opi_input(
                esd_dir, out_dir=out_dir, state=state, charge=charge,
                solvent=solvent, metals=metals,
                main_basisset=main_basisset,
                metal_basisset=metal_basisset,
                config={**config, "ESD_modus": "TDDFT"})
        if state_upper == "T1":
            return write_esd_state_opi_input(
                esd_dir, out_dir=out_dir, state=state, charge=charge,
                solvent=solvent, metals=metals,
                main_basisset=main_basisset,
                metal_basisset=metal_basisset,
                config={**config, "ESD_modus": "deltaSCF"})
        # excited states: two files; the classic writer returns the SECOND
        # (deltaSCF) input as the main one — same here.
        root_num = int(state_upper[1:]) if len(state_upper) > 1 else 1
        is_triplet = state_upper.startswith("T")
        first_base = f"{state_upper}_first_TDDFT"
        toks1 = level("RKS", moread=True)
        blocks1 = [
            opi["BlockGeom"](maxiter=int(
                config.get("hybrid1_geom_MaxIter", 2))),
            _tddft_block_opi(opi, config, iroot=root_num,
                             irootmult="triplet" if is_triplet
                             else "singlet",
                             triplets=is_triplet, follow=True),
            *output_blocks,
        ]
        atoms = _inline_atoms(opi, esd_dir, "S0.xyz", metals, metal_basis,
                              config)
        _job_in_opi(
            opi, out_dir=esd_dir, base=first_base, tokens=toks1,
            maxcore=maxcore, pal=pal, moinp="S0.gbw", blocks=blocks1,
            atoms=atoms, charge=charge, multiplicity=1,
        )
        second_base = f"{state_upper}_second"
        toks2 = level("UKS", moread=True)
        if freq_on:
            toks2.append(freq_type)
        toks2.append("deltaSCF")
        deltascf_keywords = config.get("deltaSCF_keywords", "")
        if deltascf_keywords:
            toks2.extend(str(deltascf_keywords).split())
        scf2 = _scf_block_opi(opi, config, esd_dir, state_upper,
                              deltascf=True)
        mult2 = 3 if is_triplet else 1
        blocks2 = [scf2, *_new_basis_block(opi, metals, metal_basis),
                   *output_blocks]
        return str(_job_in_opi(
            opi, out_dir=esd_dir, base=second_base, tokens=toks2,
            maxcore=maxcore, pal=pal, moinp=f"{first_base}.gbw",
            blocks=blocks2, atoms=None,
            xyzfile=(charge, mult2, f"{first_base}.xyz"),
            charge=charge, multiplicity=mult2,
        ))

    # ---- deltaSCF mode --------------------------------------------------
    multiplicity = 3 if state_upper.startswith("T") else 1
    use_deltascf = state_upper != "S0" and state_upper != "T1"
    moinp_gbw = "S0.gbw" if use_deltascf else None
    scf_type = "RKS" if state_upper == "S0" else "UKS"
    initial_guess = (str(config.get("initial_guess", "")).split() or [""])[0]

    toks = level(scf_type)
    s0_first_job_with_freq = state_upper == "S0" and (
        freq_on or _redox_steps_requested(config))
    if s0_first_job_with_freq:
        toks.append("FREQ")
    elif freq_on:
        if state_upper == "S0":
            toks.append("FREQ")
        else:
            toks.append(freq_type)
    if use_deltascf:
        toks.append("deltaSCF")
    if moinp_gbw:
        toks.append("MOREAD")
    if initial_guess and not (moinp_gbw and initial_guess == "PModel"):
        toks.append(initial_guess)
    if use_deltascf:
        deltascf_keywords = config.get("deltaSCF_keywords", "")
        if deltascf_keywords:
            toks.extend(str(deltascf_keywords).split())

    scf_block = _scf_block_opi(opi, config, esd_dir, state_upper,
                               deltascf=use_deltascf)
    if state_upper in ("S0", "T1"):
        # the classic writer writes these jobs with inline coordinates,
        # then appends a TDDFT check job as $new_job in the same file
        xyz_name = "S0.xyz"
        atoms = _inline_atoms(opi, esd_dir, xyz_name, metals, metal_basis,
                              config)
        main = _job_in_opi(
            opi, out_dir=out_dir, base=f"_{state_upper}_main_tmp",
            tokens=toks, maxcore=maxcore, pal=pal, moinp=moinp_gbw,
            blocks=[scf_block, *output_blocks], atoms=atoms,
            charge=charge, multiplicity=multiplicity,
        )
        check_toks = [functional, "RKS", rel_token, main_basis, disp_corr,
                      ri_jkx, aux_jk]
        if solv:
            check_toks.append(solv)
        check = _tddft_check_job(
            opi, esd_dir, out_dir, state_upper, charge, metals,
            metal_basis, config, [t for t in check_toks if str(t).strip()],
            maxcore, pal, output_blocks)
        merged = _merge_new_job_files([str(main), check])
        final = out_dir / f"{state_upper}.inp"
        final.write_text(Path(merged).read_text(encoding="utf-8"),
                         encoding="utf-8")
        Path(merged).unlink()
        return str(final)

    # S2..S6, T2..T6 deltaSCF jobs
    atoms = _inline_atoms(opi, esd_dir, "S0.xyz", metals, metal_basis,
                          config)
    blocks = [scf_block, *output_blocks]
    return str(_job_in_opi(
        opi, out_dir=out_dir, base=state_upper, tokens=toks,
        maxcore=maxcore, pal=pal, moinp=moinp_gbw, blocks=blocks,
        atoms=atoms, charge=charge, multiplicity=multiplicity,
    ))


# ---------------------------------------------------------------------------
# Rate jobs (ISC / IC / FLUOR / PHOSP)
# ---------------------------------------------------------------------------

def _rate_geometry(opi, esd_dir: Path, xyz_name: str, metals,
                   metal_basis, config) -> List[Any]:
    return _inline_atoms(opi, esd_dir, xyz_name, metals, metal_basis,
                         config)


def write_esd_rate_opi_input(
    esd_dir,
    *,
    out_dir,
    kind: str,
    pair: Optional[str],
    charge: int,
    solvent: Optional[str],
    metals: Sequence[str],
    main_basisset: Optional[str],
    metal_basisset: Optional[str],
    config: Mapping[str, Any],
    trootssl: Optional[int] = None,
    initial_state: str = "S1",
) -> str:
    """OPI equivalent of ``create_{isc,ic,fluor,phosp}_input``."""
    from delfin.esd_input_generator import (
        _build_solvation_keyword,
        _parse_iroot_spec,
        _parse_state_root,
        _resolve_state_filename,
        calculate_dele_cm1,
    )
    from delfin.common import tddft_settings

    opi = _import_opi()
    esd_dir = Path(esd_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    kind = kind.strip().lower()
    functional = config.get("functional", "PBE0")
    disp_corr = config.get("disp_corr", "D4")
    ri_jkx = config.get("ri_jkx", "RIJCOSX")
    main_basis, metal_basis, rel_token, aux_jk = _level_tokens(
        metals, config, main_basisset, metal_basisset)
    solv = _build_solvation_keyword(
        config.get("implicit_solvation_model", ""), solvent)
    pal = int(config.get("PAL", 12))
    maxcore = int(config.get("maxcore", 6000))
    esd_mode = str(config.get("ESD_modus", "tddft")).strip().lower()
    if "|" in esd_mode:
        esd_mode = esd_mode.split("|")[0].strip()

    def level(scf_type: str) -> List[str]:
        toks = [scf_type, functional, rel_token, main_basis, disp_corr,
                ri_jkx, aux_jk]
        if solv:
            toks.append(solv)
        return [t for t in toks if str(t).strip()]

    def esd_rate_file(base: str) -> None:
        raise NotImplementedError  # pragma: no cover - placeholder

    # ---- ISC ------------------------------------------------------------
    if kind == "isc":
        initial_state, final_state = str(pair).split(">")
        initial_state = initial_state.strip().upper()
        final_state = final_state.strip().upper()
        init_type, init_root = _parse_state_root(initial_state)
        final_type, final_root = _parse_state_root(final_state)
        if trootssl is None:
            trootssl_str = str(config.get("TROOTSSL", "0")).strip()
        else:
            trootssl_str = str(trootssl)
        trootssl_int = int(trootssl_str)
        from delfin.esd_input_generator import _format_ms_suffix
        ms_suffix = _format_ms_suffix(trootssl_int)
        job_name = f"{initial_state}_{final_state}_ISC_{ms_suffix}"

        xyz_file = _resolve_state_filename(final_state, "xyz", esd_mode)
        initial_out = _resolve_state_filename(initial_state, "out",
                                              esd_mode)

        toks = _rate_tokens(level("RKS"))
        s_root = init_root if init_type == "S" else (
            final_root if final_type == "S" else 1)
        t_root = init_root if init_type == "T" else (
            final_root if final_type == "T" else 1)
        nroots = tddft_settings.job_nroots(config, "ESD_ISC_NROOTS")
        tddft = _tddft_block_opi(
            opi, config, nroots=nroots, dosoc=True,
            job_lines=[("sroot", s_root), ("troot", t_root),
                       ("trootssl", trootssl_str)])
        header = [("ISCISHESS", f'"{_resolve_state_filename(initial_state, "hess", esd_mode)}"'),
                  ("ISCFSHESS", f'"{_resolve_state_filename(final_state, "hess", esd_mode)}"')]
        esd = _esd_block_lines(opi, config, header)
        atoms = _rate_geometry(opi, esd_dir, xyz_file, metals, metal_basis,
                               config)
        return str(_job_in_opi(
            opi, out_dir=out_dir, base=job_name, tokens=toks,
            maxcore=maxcore, pal=pal, moinp=None,
            blocks=[tddft, esd], atoms=atoms, charge=charge,
            multiplicity=1,
        ))

    # ---- IC -------------------------------------------------------------
    if kind == "ic":
        initial_state, final_state = str(pair).split(">")
        initial_state = initial_state.strip().upper()
        final_state = final_state.strip().upper()
        init_type, init_root = _parse_state_root(initial_state)
        final_type, _ = _parse_state_root(final_state)
        job_name = f"{initial_state}_{final_state}_IC"
        xyz_file = _resolve_state_filename(final_state, "xyz", esd_mode)
        final_multiplicity = 3 if final_type == "T" else 1

        nroots = int(
            str(config.get("ESD_IC_NROOTS",
                           config.get("TDDFT_nroots", 40))).strip())
        iroot = init_root - 1 if (final_type == "T" and final_state == "T1") \
            else init_root
        nacme = str(config.get("NACME", "TRUE")).upper()
        etf = str(config.get("ETF", "TRUE")).upper()
        toks = _rate_tokens(level("RKS"))
        tddft = _tddft_block_opi(
            opi, config, nroots=nroots, iroot=iroot, dosoc=None,
            job_lines=[("nacme", nacme), ("etf", etf)])
        header = [("GSHESSIAN", f'"{_resolve_state_filename(final_state, "hess", esd_mode)}"'),
                  ("ESHESSIAN", f'"{_resolve_state_filename(initial_state, "hess", esd_mode)}"')]
        esd = _esd_block_lines(opi, config, header)
        atoms = _rate_geometry(opi, esd_dir, xyz_file, metals, metal_basis,
                               config)
        return str(_job_in_opi(
            opi, out_dir=out_dir, base=job_name, tokens=toks,
            maxcore=maxcore, pal=pal, moinp=None,
            blocks=[tddft, esd], atoms=atoms, charge=charge,
            multiplicity=final_multiplicity,
        ))

    # ---- FLUOR ----------------------------------------------------------
    if kind == "fluor":
        initial_state = initial_state.strip().upper()
        final_state = "S0"
        init_type, init_root = _parse_state_root(initial_state)
        if init_type != "S" or initial_state == "S0":
            raise ValueError(
                f"Fluorescence only supported for Sn->S0 (got {initial_state})")
        job_name = f"{initial_state}_{final_state}_FLUOR"
        iroot = int(config.get("ESD_FLUOR_IROOT", init_root))
        nroots = tddft_settings.job_nroots(config, "ESD_FLUOR_NROOTS")
        toks = _rate_tokens(level("RKS"))
        tddft = _tddft_block_opi(
            opi, config, nroots=nroots, iroot=iroot, dosoc=None)
        header = [("GSHESSIAN", f'"{_resolve_state_filename(final_state, "hess", esd_mode)}"'),
                  ("ESHESSIAN", f'"{_resolve_state_filename(initial_state, "hess", esd_mode)}"')]
        esd = _esd_block_lines(opi, config, header)
        final_xyz = _resolve_state_filename(final_state, "xyz", esd_mode)
        atoms = _rate_geometry(opi, esd_dir, final_xyz, metals,
                               metal_basis, config)
        return str(_job_in_opi(
            opi, out_dir=out_dir, base=job_name, tokens=toks,
            maxcore=maxcore, pal=pal, moinp=None,
            blocks=[tddft, esd], atoms=atoms, charge=charge,
            multiplicity=1,
        ))

    # ---- PHOSP ----------------------------------------------------------
    if kind == "phosp":
        initial_state = initial_state.strip().upper()
        final_state = "S0"
        init_type, _ = _parse_state_root(initial_state)
        if init_type != "T" or initial_state != "T1":
            raise ValueError(
                f"Phosphorescence only supported for T1->S0 (got {initial_state})")
        job_name = f"{initial_state}_{final_state}_PHOSP"
        t1_out = esd_dir / "T1.out"
        s0_out = esd_dir / "S0.out"
        dele_cm1 = calculate_dele_cm1(str(t1_out), str(s0_out))
        dele_int = int(round(dele_cm1)) if dele_cm1 is not None else None
        nroots = tddft_settings.job_nroots(config, "ESD_PHOSP_NROOTS")
        phosp_dosoc = tddft_settings.coerce_bool(
            config.get("ESD_PHOSP_DOSOC", "TRUE"), "ESD_PHOSP_DOSOC")
        iroots = _parse_iroot_spec(config.get("phosp_IROOT", None),
                                   default=[1, 2, 3])
        final_xyz = _resolve_state_filename(final_state, "xyz", esd_mode)
        atoms = _rate_geometry(opi, esd_dir, final_xyz, metals,
                               metal_basis, config)
        header = [("GSHESSIAN", f'"{_resolve_state_filename(final_state, "hess", esd_mode)}"'),
                  ("TSHESSIAN", f'"{_resolve_state_filename(initial_state, "hess", esd_mode)}"')]
        esd = _esd_block_lines(opi, config, header)
        if dele_int is not None:
            esd.modify_option("dele", str(dele_int))

        jobs: List[str] = []
        for iroot in iroots:
            tddft = _tddft_block_opi(
                opi, config, nroots=nroots, iroot=iroot,
                dosoc=phosp_dosoc)
            path = _job_in_opi(
                opi, out_dir=out_dir,
                base=f"{job_name}_iroot{iroot}", tokens=toks,
                maxcore=maxcore, pal=pal, moinp=None,
                blocks=[tddft, esd], atoms=atoms, charge=charge,
                multiplicity=1,
            )
            jobs.append(str(path))
        # The classic writer emits ONE file with $new_job between the
        # iroot jobs; OPI's Calculator writes one file per input.  The
        # multi-job merge is done textually by the caller adapter; here
        # the first iroot file is the canonical one (same options).
        return jobs[0] if len(jobs) == 1 else _merge_new_job_files(jobs)

    raise ValueError(f"Unknown ESD rate kind: {kind!r}")


def _merge_new_job_files(paths: Sequence[str]) -> str:
    """Merge per-job inputs into one $new_job file (PHOSP iroot jobs)."""
    sections = []
    for i, p in enumerate(paths):
        text = Path(p).read_text(encoding="utf-8")
        if i:
            text = "\n$new_job\n\n" + text
        sections.append(text)
    merged = Path(paths[0])
    merged.write_text("".join(sections), encoding="utf-8")
    for p in paths[1:]:
        Path(p).unlink()
    return str(merged)


def _tddft_check_job(
    opi, esd_dir: Path, out_dir: Path, state_upper: str, charge: int,
    metals: Sequence[str], metal_basis: Optional[str],
    config: Mapping[str, Any],
    tokens: Sequence[str], maxcore: int, pal: int,
    output_blocks: Sequence[Any],
) -> str:
    """The ``{state}_TDDFT`` check job, written to its own file.

    The classic writer appends this job to the state's own file with
    ``$new_job``; the caller merges it there after this helper writes it.
    Body per esd_input_generator.py:969-1008: RKS level line, %base,
    %pal, %maxcore, %output, ``%tddft triplets`` (no FollowIRoot),
    %basis for metals on an xyzfile geometry, multiplicity 1.
    """
    blocks = [
        _tddft_block_opi(opi, config, triplets=True),
        *_new_basis_block(opi, metals, metal_basis),
        *output_blocks,
    ]
    return str(_job_in_opi(
        opi, out_dir=out_dir,
        base=f"_{state_upper}_TDDFT_tmp",  # merged away below
        tokens=list(tokens), maxcore=maxcore, pal=pal, moinp=None,
        blocks=blocks, atoms=None,
        xyzfile=(charge, 1, f"{state_upper}.xyz"),
        charge=charge, multiplicity=1,
    ))
