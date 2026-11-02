"""DELFIN's ORCA geometry-optimisation / frequency building block on OPI.

This module builds ORCA geometry-optimisation, frequency, relaxed surface
scan, NEB and IRC inputs with the ORCA Python Interface (OPI, orca-pi 2.0.0)
and reads the results back through OPI.  It is the geometry-aware sibling of
``delfin/tools/adapters/orca_sp_opi.py`` (a single-point DFT block): the
method keywords, ``%pal``/``%maxcore`` and geometry come from OPI's classes
(``opi.input.Input``, ``opi.input.structures``, the typed blocks ``BlockGeom``
/ ``BlockFreq`` / ``BlockNeb`` / ``BlockIrc``), the binary is still launched
by DELFIN's existing runner (``delfin.orca.run_orca``), and every keyword and
block option is checked against DELFIN's ORCA knowledge store before the
input is written.  Everything here is a wrapper, not a replacement.

Two rules bind the reader, as in the single-point block:

* A value from OPI is only meaningful once the matching convergence flag says
  so.  The final energy and the converged geometry are reported only after
  ``terminated_normally()``; thermochemistry and frequencies only after the
  frequency calculation finished (``get_ir`` / ``get_zpe`` are only meaningful
  once the frequency job ran to completion); an optimisation that did not
  converge keeps ``opt_converged=False``.
* A minimum with an imaginary mode is reported as such, never silently
  accepted: when a converged frequency run has one or more imaginary modes,
  ``character`` is ``"transition-state-like"`` and the imaginary
  wavenumbers are listed.

An unknown keyword or block option is an error returned to the caller, never
dropped; OPI's typed blocks cannot express a multi-coordinate ``%geom Scan``
(their ``BlockGeom.scan`` holds a single ``Scan``), so asking for more than
one scan coordinate is refused rather than spilled as a hand-written string
template.

The module adds two adapters (``orca_opt_opi``) to the existing step registry;
nothing in the existing modules is changed.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import ErrorKind, StepResult, StepStatus
from delfin.tools.adapters.orca_sp_opi import (
    _failure_lines,
    _validated_keyword,
    validate_method,
)

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Task composition: task name -> (list of OPI SimpleKeyword enum members)
# ---------------------------------------------------------------------------

# The opt level is a simple keyword its own right (manual optimizations.txt:
# "!TightOpt", "!VeryTightOpt" choose a set of convergence thresholds).  The
# job-type keywords Task.OPT/FREQ/NUMFREQ and Neb.NEB/NEB_TS/IRC live on the
# Task / Neb enum boxes.  Mapping task-name to the enum members reproduces the
# ORCA input line exactly; every name is validated against known_keyword below.
_TASK_TO_ENUMS: dict[str, tuple[str, ...]] = {
    "opt": ("opt",),
    "tightopt": ("tightopt",),
    "verytightopt": ("verytightopt",),
    "optts": ("optts",),
    "freq": ("freq",),
    "numfreq": ("numfreq",),
    "opt+freq": ("opt", "freq"),
    "tightopt+freq": ("tightopt", "freq"),
    "verytightopt+freq": ("verytightopt", "freq"),
    "neb": ("neb",),
    "neb-ts": ("neb-ts",),
    "irc": ("irc",),
}


@dataclass
class OrcaOptResult:
    """Typed result of one OPI-backed ORCA opt / freq / scan / neb / irc run.

    Numeric and structural fields sourced from OPI are populated only when the
    matching convergence flags make them meaningful (see module docstring):
    :attr:`meaningful` summarises that state, and every numeric field is
    ``None`` when it would not be trustworthy.
    """

    ok: bool
    work_dir: Path
    output_file: Optional[Path] = None

    # Convergence flags (always read; they gate every numeric field).
    terminated_normally: Optional[bool] = None
    scf_converged: Optional[bool] = None
    opt_converged: Optional[bool] = None
    freq_converged: Optional[bool] = None
    meaningful: bool = False

    # Geometry + energy; reported only when the run terminated normally.
    geometry: Any = None            # opi Structure (converged geometry)
    final_energy: Optional[float] = None
    el_energy: Optional[float] = None

    # Frequencies + thermochemistry; only after the frequency block completed.
    n_modes: Optional[int] = None
    n_imag: Optional[int] = None
    imag_freqs: tuple[float, ...] = ()
    lowest_freq: Optional[float] = None
    zpe: Optional[float] = None
    inner_energy: Optional[float] = None
    enthalpy: Optional[float] = None
    entropy: Optional[float] = None
    free_energy: Optional[float] = None
    character: Optional[str] = None

    # Method knowledge findings (warnings) surfaced, never silently dropped.
    findings: tuple[dict, ...] = ()
    # Lines of the .out text that explain a failure (OPI does not give them).
    failure_lines: tuple[str, ...] = ()
    error: Optional[str] = None


def validate_opt_method(
    *,
    functional: str,
    basis: str,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    task: str = "opt",
    geom_options: Optional[Mapping[str, Any]] = None,
    freq_options: Optional[Mapping[str, Any]] = None,
    neb_options: Optional[Mapping[str, Any]] = None,
    irc_options: Optional[Mapping[str, Any]] = None,
    scan: Optional[Sequence[tuple]] = None,
    cards: Optional[Sequence[dict]] = None,
) -> tuple[list[str], list[dict]]:
    """Check every keyword and block option against DELFIN's ORCA knowledge.

    Reuses ``orca_sp_opi.validate_method`` for the method keywords and then
    validates the job-type keyword (``task``), every ``%geom``/``%freq``/
    ``%neb``/``%irc`` block option and the scan coordinate against
    ``delfin.orca_knowledge``.  Returns ``(errors, warnings)``; a non-empty
    ``errors`` must stop the step (the run wrapper refuses to launch).
    """
    errors, warnings = validate_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        cards=cards,
    )

    from delfin.orca_knowledge import known_block_option

    task_key = (task or "opt").lower()
    if task_key not in _TASK_TO_ENUMS:
        errors.append(
            f"task: keyword {task!r} is not a supported opt/freq/neb/irc "
            "job type (OPI block model); refused."
        )
    else:
        _validated_keyword(errors, task_key.split("+")[0], "task")

    for block, opts in (
        ("geom", geom_options),
        ("freq", freq_options),
        ("neb", neb_options),
        ("irc", irc_options),
    ):
        if opts:
            for option in opts:
                if not known_block_option(block, option):
                    errors.append(
                        f"{block}: block option {option!r} is not one ORCA "
                        f"knows (OPI block model); refused."
                    )

    if scan:
        if len(scan) > 1:
            errors.append(
                "scan: OPI's BlockGeom.scan holds one coordinate; a relaxed "
                "surface scan with multiple coordinates cannot be expressed "
                "through OPI and is refused (manual optimizations_scans.txt "
                "allows up to 3; no string template is emitted)."
            )
        else:
            for entry in scan:
                _validate_scan_entry(errors, entry)
    return errors, warnings


def _validate_scan_entry(errors: list[str], entry: tuple) -> None:
    """Validate one ``(coordinate, start, end, steps)`` scan entry.

    ``coordinate`` is a tuple ``(mode, atom1[, atom2[, atom3[, atom4]]])``
    with mode in {"A", "B", "D"}, matching ORCA's ``Scan`` block (manual
    optimizations_scans.txt, table 4.5: B = two atoms bond, A = three atoms
    angle, D = four atoms dihedral).  Atoms are 0-indexed in OPI's atom model
    (the ``Input`` consistency, matching the ``Structure``/``Atom`` order);
    a negative index is refused up front.
    """
    if not isinstance(entry, (tuple, list)) or len(entry) < 3:
        errors.append(
            f"scan: entry {entry!r} must be (coordinate, start, end[, steps])."
        )
        return
    coord, start, end = entry[0], entry[1], entry[2]
    steps = entry[3] if len(entry) > 3 else None
    mode = coord[0] if isinstance(coord, (tuple, list)) and coord else ""
    modeU = str(mode).upper()
    if modeU not in {"A", "B", "D"}:
        errors.append(f"scan: coordinate mode {mode!r} not one of A/B/D.")
        return
    if not isinstance(start, (int, float)) or not isinstance(end, (int, float)):
        errors.append(f"scan: start/end {start!r}, {end!r} must be numbers.")
    if steps is not None and (not isinstance(steps, int) or steps < 0):
        errors.append(f"scan: steps {steps!r} must be a non-negative integer.")
    for atom in coord[1:]:
        if not isinstance(atom, int) or atom < 0:
            errors.append(f"scan: atom index {atom!r} must be a non-negative int.")


def _enum_for_keyword(keyword: str) -> tuple[Any, str]:
    """Map one validated simple-keyword name onto its OPI enum member and the
    enum box name (e.g. 'opt' -> (Opt.OPT, 'Opt'), 'freq' -> (Task.FREQ,
    'Task')), so the input is built from OPI's typed enums, never a string.
    """
    from opi.input.simple_keywords.opt import Opt
    from opi.input.simple_keywords.task import Task
    from opi.input.simple_keywords.neb import Neb

    kw = keyword.lower()
    box_candidates = (Opt, Task, Neb)
    lookup = {str(m).lower(): m for b in box_candidates for m in b.__dict__.values()}
    attr = lookup.get(kw)
    if attr is None:
        # Fall back to a bare validated string (the caller already validated
        # the name against the knowledge store).
        return str(keyword), "bare"
    return attr, type(attr).__name__


def _build_opi_input(
    work_dir: Path,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str,
    basis: str,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    task: str = "opt",
    geom_options: Optional[Mapping[str, Any]] = None,
    freq_options: Optional[Mapping[str, Any]] = None,
    neb_options: Optional[Mapping[str, Any]] = None,
    irc_options: Optional[Mapping[str, Any]] = None,
    scan: Optional[Sequence[tuple]] = None,
    scan_as_single: bool = False,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "calc",
) -> Path:
    """Assemble the ORCA input with OPI's classes and write it with OPI.

    Returns the absolute path of the written ``<base_name>.inp``.  ``%pal``/
    ``%maxcore`` map onto ``Input.ncores``/``Input.memory``; the job-type
    simple keyword and the ``%geom`` / ``%freq`` / ``%neb`` / ``%irc`` blocks
    come from OPI's typed pydantic classes (``BlockGeom``, ``BlockFreq``,
    ``BlockNeb``, ``BlockIrc``, ``Scan``, ``Constraints``).  The file is
    written by OPI's ``Calculator.write_input()`` — no hand-written string
    template.  ``opts`` (a block-name -> option dict) must already be
    validated (run wrappers call ``validate_opt_method`` first).
    """
    # Structural rules no OPI block can express are refused before building.
    if scan and len(scan) > 1:
        raise ValueError(
            "scan: OPI's BlockGeom.scan holds one coordinate; a relaxed "
            "surface scan with multiple coordinates cannot be expressed "
            "through OPI and is refused."
        )
    if geom_options and str(geom_options.get("inhess") or "").lower() == "read":
        if not (geom_options.get("inhessname") or geom_options.get("inhessname2")):
            raise ValueError(
                "inhess=read needs an inhessname/inhessname2 (the Hessian file "
                "to read); the manual requires the file when reusing a Hessian."
            )
    # --- job-type simple keywords, from OPI's typed enums -----------------
    writer = _build_common_keywords(
        work_dir,
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task=task,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    writer.add_blocks_from(
        task=task,
        geom_options=geom_options,
        freq_options=freq_options,
        neb_options=neb_options,
        irc_options=irc_options,
        scan=scan,
    )
    return writer.write()


class _Writer:
    """Small mutable holder for assembling an OPI input (opaque to callers)."""

    def __init__(
        self,
        work_dir: Path,
        elements: Sequence[str],
        coordinates: Sequence[Sequence[float]],
        charge: int,
        mult: int,
        base_name: str,
    ) -> None:
        from opi.input import Input
        from opi.input.structures import Atom, Structure

        self.work_dir = Path(work_dir)
        self.base_name = base_name
        self.inp = Input()
        atoms = [
            Atom(str(el), coordinates=(float(x), float(y), float(z)))
            for el, (x, y, z) in zip(elements, coordinates)
        ]
        self.structure = Structure(atoms, charge=charge, multiplicity=mult)

    def set_keywords(self, keywords: Sequence[Any]) -> None:
        self.inp.add_simple_keywords(*keywords)

    def set_pal(self, ncores: int, maxcore: int) -> None:
        self.inp.ncores = ncores
        self.inp.memory = maxcore

    def add_blocks_from(
        self,
        *,
        task: str,
        geom_options: Optional[Mapping[str, Any]],
        freq_options: Optional[Mapping[str, Any]],
        neb_options: Optional[Mapping[str, Any]],
        irc_options: Optional[Mapping[str, Any]],
        scan: Optional[Sequence[tuple]],
    ) -> None:
        from opi.input.blocks.block_geom import (
            BlockGeom,
            Scan,
            TSMode,
        )
        from opi.input.blocks.block_irc import BlockIrc
        from opi.input.blocks.block_freq import BlockFreq
        from opi.input.blocks.block_neb import BlockNeb

        blocks: list[Any] = []

        geom = dict(geom_options or {})
        # A TSMode string for OptTS is accepted by OPI's validator; map the
        # option name to the OPI block field.  Constraints: BlockGeom's
        # ``constraints`` field validator (block_geom.Constraint.constraint_init)
        # already coerces a list of strings ("{ B 0 1 C }", ...) into OPI
        # Constraint objects, so pass the list straight through -- freezing
        # atoms ({ C N C }), bonds ({ B N1 N2 C }) or angles ({ A 1 2 3 C })
        # per manual optimisation section 4.5.
        if "ts_mode" in geom:
            # TSMode.from_string accepts only the bracketed string form
            # ('{ M 1 }' per manual section 4.5.1); map an int x (follow mode
            # x uphill) or a list ['M', 1] onto it -- anything else is a user
            # error and stays one.
            raw = geom.pop("ts_mode")
            if isinstance(raw, int):
                raw = f"{{ M {raw} }}"
            elif isinstance(raw, (list, tuple)) and raw and str(raw[0]).upper() \
                    in {"A", "B", "D", "M"}:
                raw = "{ " + " ".join(str(part) for part in raw) + " }"
            geom["ts_mode"] = TSMode.from_string(raw)
        # keep constraints as a list of strings for OPI's coercion
        # Single-coordinate scan -> OPI Scan object (or its string).
        if scan:
            if len(scan) > 1:
                raise ValueError(
                    "scan: OPI's BlockGeom.scan holds one coordinate; a "
                    "relaxed surface scan with multiple coordinates cannot be "
                    "expressed through OPI and is refused."
                )
            if "scan" in geom and geom["scan"] is not None:
                # explicit scan block option wins; it must be a Scan already
                pass
            else:
                entry = scan[0]
                coord, start, end = entry[0], entry[1], entry[2]
                steps = entry[3] if len(entry) > 3 else None
                mode = str(coord[0]).upper()
                kwargs: dict[str, Any] = {"mode": mode}
                for i, atom in enumerate(coord[1:], start=1):
                    kwargs[f"atom{i}"] = atom
                if steps is not None:
                    kwargs.update(startvalue=float(start), endvalue=float(end),
                                  steps=int(steps))
                else:
                    kwargs.update(valuelist=[float(start), float(end)])
                geom["scan"] = Scan(**kwargs)
        if geom:
            blocks.append(BlockGeom(**geom))
        if freq_options:
            blocks.append(BlockFreq(**dict(freq_options)))
        if neb_options:
            blocks.append(BlockNeb(**dict(neb_options)))
        if irc_options:
            blocks.append(BlockIrc(**dict(irc_options)))
        if blocks:
            self.inp.add_blocks(*blocks)

    def write(self) -> Path:
        from opi.core import Calculator

        calc = Calculator(self.base_name, working_dir=self.work_dir,
                          version_check=False)
        calc.input = self.inp
        calc.structure = self.structure
        calc.write_input(force=True)
        return self.work_dir / f"{self.base_name}.inp"


def _build_common_keywords(
    work_dir: Path,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str,
    basis: str,
    ri: str,
    aux_basis: str,
    dispersion: str,
    extra_keywords: Optional[Sequence[str]],
    task: str,
    maxcore: int,
    ncores: int,
    base_name: str,
) -> "_Writer":
    """Shared method-keyword assembly (functional/basis/ri/dispersion/...)."""
    from opi.input.simple_keywords.approximation import Approximation
    from opi.input.simple_keywords.base import SimpleKeyword
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dft import Dft
    from opi.input.simple_keywords.dispersion_correction import (
        DispersionCorrection,
    )

    def _lookup(enum_cls: type, keyword: str) -> Any:
        target = str(keyword).lower()
        for attr in vars(enum_cls).values():
            if isinstance(attr, SimpleKeyword) and str(attr).lower() == target:
                return attr
        return str(keyword)

    writer = _Writer(work_dir, elements, coordinates, charge, mult, base_name)
    keywords: list[Any] = []
    if functional:
        keywords.append(_lookup(Dft, functional))
    keywords.append(_lookup(BasisSet, basis))
    if ri:
        keywords.append(_lookup(Approximation, ri))
    if aux_basis:
        keywords.append(_lookup(BasisSet, aux_basis))
    if dispersion:
        keywords.append(_lookup(DispersionCorrection, dispersion))
    # The job-type keyword(s): 'opt+freq' -> Opt.OPT and Task.FREQ.
    for name in _TASK_TO_ENUMS[(task or "opt").lower()]:
        enum_member, _box = _enum_for_keyword(name)
        keywords.append(enum_member)
    for kw in extra_keywords or ():
        keywords.append(str(kw))
    writer.set_keywords(keywords)
    writer.set_pal(ncores=ncores, maxcore=maxcore)
    return writer


def _read_results(
    work_dir: Path,
    base_name: str,
    *,
    has_freq: bool = False,
    has_geom: bool = True,
) -> OrcaOptResult:
    """Read a finished run with OPI, honouring convergence-gated reporting.

    ``has_freq`` says the run carried a frequency block (so frequencies and
    thermochemistry are expected and gated on the frequency job finishing);
    ``has_geom`` says it carried a geometry step (converged geometry /
    ``opt_converged`` gate).  The matching convergence flag(s) have to hold
    before the corresponding numbers are reported.
    """
    from opi.output.core import Output

    out_path = work_dir / f"{base_name}.out"
    try:
        out = Output(basename=base_name, working_dir=work_dir)
        out.parse()
    except Exception as exc:  # OPI could not read the run at all
        return OrcaOptResult(
            ok=False,
            work_dir=work_dir,
            output_file=out_path if out_path.is_file() else None,
            meaningful=False,
            failure_lines=_failure_lines(out_path),
            error=f"OPI could not parse the run: {exc}",
        )

    terminated = out.terminated_normally()
    scf_ok = out.scf_converged()
    opt_ok = out.geometry_optimization_converged() if has_geom else True
    result = OrcaOptResult(
        ok=False,
        work_dir=work_dir,
        output_file=out_path if out_path.is_file() else None,
        terminated_normally=terminated,
        scf_converged=scf_ok,
        opt_converged=opt_ok if has_geom else None,
        meaningful=False,
    )
    # The frequency block is meaningful only when the job terminated and the
    # opt (if any) converged; OPI's get_ir reads the final IR block, which only
    # exists when the frequency job ran.
    freq_ok = bool(terminated and (not has_geom or opt_ok))
    result.freq_converged = freq_ok if has_freq else None
    result.meaningful = bool(terminated and scf_ok and opt_ok)
    result.ok = result.meaningful

    if not result.meaningful:
        result.failure_lines = _failure_lines(out_path)
        result.error = (
            "run did not finish normally: "
            f"terminated_normally={terminated}, scf_converged={scf_ok}, "
            f"opt_converged={opt_ok if has_geom else 'n/a'}"
        )
        return result

    # --- geometry + energy (gated on terminated_normally + opt conv.) -----
    try:
        result.final_energy = out.get_final_energy()
        result.el_energy = out.get_el_energy()
        result.geometry = out.get_structure(index=-1)
    except Exception:
        logger.debug("OPI geometry/energy read failed for %s", work_dir,
                     exc_info=True)
    if result.geometry is None:
        result.geometry = out.get_structure()  # fall back to the single frame

    # --- frequencies + thermochemistry (gated on the freq block finishing) --
    if has_freq and result.freq_converged:
        try:
            ir = out.get_ir() or {}
            freqs = [m.wavenumber for m in ir.values()]
            result.n_modes = len(freqs)
            result.n_imag = sum(1 for f in freqs if f < 0)
            result.imag_freqs = tuple(sorted(f for f in freqs if f < 0))
            positives = [f for f in freqs if f >= 0]
            result.lowest_freq = min(positives) if positives else None
            result.zpe = out.get_zpe()
            result.inner_energy = out.get_inner_energy()
            result.enthalpy = out.get_enthalpy()
            result.entropy = out.get_entropy()
            result.free_energy = out.get_free_energy()
            # A converged minimum with an imaginary mode is reported as such.
            if result.n_imag and result.n_imag > 0:
                result.character = (
                    f"transition-state-like ({result.n_imag} imaginary mode(s))"
                )
            else:
                result.character = "minimum (no imaginary modes)"
        except Exception:
            logger.debug("OPI freq/thermo read failed for %s", work_dir,
                         exc_info=True)
    return result


def _run(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str,
    basis: str,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    task: str = "opt",
    geom_options: Optional[Mapping[str, Any]] = None,
    freq_options: Optional[Mapping[str, Any]] = None,
    neb_options: Optional[Mapping[str, Any]] = None,
    irc_options: Optional[Mapping[str, Any]] = None,
    scan: Optional[Sequence[tuple]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "calc",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> OrcaOptResult:
    """Common validate -> build -> run -> read pipeline."""
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else _make_workdir()
    wdir.mkdir(parents=True, exist_ok=True)

    errors, findings = validate_opt_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task=task,
        geom_options=geom_options,
        freq_options=freq_options,
        neb_options=neb_options,
        irc_options=irc_options,
        scan=scan,
    )
    if errors:
        return OrcaOptResult(
            ok=False,
            work_dir=wdir,
            meaningful=False,
            findings=tuple(findings),
            error="; ".join(errors),
        )

    inp_path = _build_opi_input(
        wdir,
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task=task,
        geom_options=geom_options,
        freq_options=freq_options,
        neb_options=neb_options,
        irc_options=irc_options,
        scan=scan,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    out_path = wdir / f"{base_name}.out"
    success = run_orca(str(inp_path), str(out_path), timeout=timeout,
                       working_dir=wdir)
    has_freq = bool(task and ("freq" in task.lower() or "numfreq" in task.lower()))
    has_geom = bool(task and not task.lower().startswith("freq") and "neb" not in task.lower() and "irc" not in task.lower())
    result = _read_results(wdir, base_name, has_freq=has_freq, has_geom=has_geom)
    result.findings = tuple(findings)
    if not success and result.ok:
        result.ok = False
        result.meaningful = False
        result.error = result.error or "ORCA subprocess reported failure"
        result.failure_lines = result.failure_lines or _failure_lines(out_path)
    return result


def run_opt_freq(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str,
    basis: str,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    opt_level: str = "opt",
    geom_options: Optional[Mapping[str, Any]] = None,
    freq_options: Optional[Mapping[str, Any]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "optfreq",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> OrcaOptResult:
    """Geometry-optimise and compute frequencies in one ORCA job.

    Opt + Freq in a single job (``! Opt Freq``): optimise with the chosen
    convergence level (``opt``/``tightopt``/``verytightopt`` for the ``! Opt``
    / ``!TightOpt`` / ``!VeryTightOpt`` thresholds, manual optimisation
    Table 4.1) and read back the converged geometry, final energy,
    frequencies and thermochemistry through OPI, gated on the convergence
    flags; a minimum with an imaginary mode is reported as such.
    """
    level = (opt_level or "opt").lower()
    if level not in {"opt", "tightopt", "verytightopt"}:
        return OrcaOptResult(
            ok=False,
            work_dir=Path(work_dir) if work_dir else Path.cwd(),
            meaningful=False,
            error=f"opt_level {opt_level!r}: must be opt/tightopt/verytightopt",
        )
    task = f"{level}+freq"
    return _run(
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task=task,
        geom_options=geom_options,
        freq_options=freq_options,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
        work_dir=work_dir,
        timeout=timeout,
    )


def run_opt_only(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str,
    basis: str,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    opt_level: str = "opt",
    geom_options: Optional[Mapping[str, Any]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "opt",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> OrcaOptResult:
    """Geometry-optimise only (no frequency block)."""
    level = (opt_level or "opt").lower()
    if level not in {"opt", "tightopt", "verytightopt", "optts"}:
        return OrcaOptResult(
            ok=False,
            work_dir=Path(work_dir) if work_dir else Path.cwd(),
            meaningful=False,
            error=f"opt_level {opt_level!r}: must be opt/tightopt/"
                  "verytightopt/optts",
        )
    return _run(
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task=level,
        geom_options=geom_options,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
        work_dir=work_dir,
        timeout=timeout,
    )


def run_scan(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str,
    basis: str,
    scan: Sequence[tuple],
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    geom_options: Optional[Mapping[str, Any]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "scan",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> OrcaOptResult:
    """Relaxed surface scan of one coordinate (OPI supports a single Scan).

    OPI's ``BlockGeom.scan`` holds one ``Scan``, so exactly one coordinate is
    accepted; a request for several is refused (manual optimizations_scans
    allows up to 3, but no string template is emitted).
    """
    return _run(
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task="opt",
        geom_options=geom_options,
        scan=scan,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
        work_dir=work_dir,
        timeout=timeout,
    )


def run_neb(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    product_elements: Sequence[str],
    product_coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str,
    basis: str,
    nimages: int,
    product_file: str = "product.xyz",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    neb_options: Optional[Mapping[str, Any]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "neb",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> OrcaOptResult:
    """NEB set-up: %neb NImages, product structure file, PreOpt.

    The replica path between the reactant and the product (a separate XYZ
    file) is set through ``%neb`` options; ``product_file`` is written into
    the working directory by the caller before this runs.  Only set-up +
    run + read-back is wired here; per-image energies live in the ORCA
    trajectory file (not surfaced by OPI 2.0.0).
    """
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else _make_workdir()
    wdir.mkdir(parents=True, exist_ok=True)

    errors, findings = validate_opt_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task="neb",
        neb_options=neb_options,
    )
    if errors:
        return OrcaOptResult(
            ok=False,
            work_dir=wdir,
            meaningful=False,
            findings=tuple(findings),
            error="; ".join(errors),
        )

    # Write the product XYZ into the run dir so %neb's product_pdbfile field
    # (an OPI InputFilePath) can point at it.
    _write_xyz(wdir / product_file, product_elements, product_coordinates)

    opts = dict(neb_options or {})
    opts.setdefault("nimages", nimages)
    opts.setdefault("preopt", True)
    opts.setdefault("neb_end_xyzfile", product_file)
    inp_path = _build_opi_input(
        wdir,
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task="neb",
        neb_options=opts,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    out_path = wdir / f"{base_name}.out"
    success = run_orca(str(inp_path), str(out_path), timeout=timeout,
                       working_dir=wdir)
    result = _read_results(wdir, base_name, has_freq=False, has_geom=False)
    result.findings = tuple(findings)
    if not success and result.ok:
        result.ok = False
        result.meaningful = False
        result.error = result.error or "ORCA subprocess reported failure"
        result.failure_lines = result.failure_lines or _failure_lines(out_path)
    return result


def run_irc(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str,
    basis: str,
    irc_options: Optional[Mapping[str, Any]] = None,
    hess_filename: str = "ts.hess",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "irc",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> OrcaOptResult:
    """IRC set-up: %irc MaxIter, Direction (and Hessian reuse)."""
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else _make_workdir()
    wdir.mkdir(parents=True, exist_ok=True)

    errors, findings = validate_opt_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task="irc",
        irc_options=irc_options,
    )
    if errors:
        return OrcaOptResult(
            ok=False,
            work_dir=wdir,
            meaningful=False,
            findings=tuple(findings),
            error="; ".join(errors),
        )
    opts = dict(irc_options or {})
    opts.setdefault("maxiter", 100)
    opts.setdefault("hess_filename", hess_filename)
    inp_path = _build_opi_input(
        wdir,
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        task="irc",
        irc_options=opts,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    out_path = wdir / f"{base_name}.out"
    success = run_orca(str(inp_path), str(out_path), timeout=timeout,
                       working_dir=wdir)
    result = _read_results(wdir, base_name, has_freq=False, has_geom=False)
    result.findings = tuple(findings)
    if not success and result.ok:
        result.ok = False
        result.meaningful = False
        result.error = result.error or "ORCA subprocess reported failure"
        result.failure_lines = result.failure_lines or _failure_lines(out_path)
    return result


def _write_xyz(path: Path, elements: Sequence[str],
               coordinates: Sequence[Sequence[float]]) -> None:
    """Write a plain XYZ file (no header guards; caller's own structure)."""
    lines = [str(len(elements)), "generated by orca_opt_opi"]
    for el, (x, y, z) in zip(elements, coordinates):
        lines.append(f"{el} {x:.8f} {y:.8f} {z:.8f}")
    path.write_text("\n".join(lines) + "\n")


def _make_workdir() -> Path:
    import tempfile

    return Path(tempfile.mkdtemp(prefix="orca_opt_opi_", dir=Path.cwd()))


# --------------------------------------------------------------------------
# Adapter for delfin.tools (registered as ``orca_opt_opi``)
# --------------------------------------------------------------------------


class OrcaOptOpiAdapter(StepAdapter):
    """An ORCA geometry-optimisation + frequency step driven by OPI.

    Builds and reads through OPI (single-point sibling of ``orca_sp_opi``),
    launching through DELFIN's own runner.  Values are reported only once the
    convergence flags make them meaningful; an unknown keyword / block option
    stops the step; a minimum with an imaginary mode is reported as such.
    """

    name = "orca_opt_opi"
    description = "OPI-backed ORCA geometry optimisation + frequencies"
    category = "dft"
    produces_geometry = True
    consumes = ("geometry",)
    produces = ("gbw", "qc_output", "property_json")
    requires_python = ("opi",)

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        kwargs.pop("_prev_artifacts", None)
        start = time.monotonic()
        try:
            elements = kwargs["elements"]
            coordinates = kwargs["coordinates"]
        except KeyError as exc:
            return self._make_result(
                self.name,
                StepStatus.FAILED,
                work_dir,
                start,
                error=f"missing required {exc}",
                error_kind=ErrorKind.MISSING_PARAM,
            )
        res = run_opt_freq(
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            functional=kwargs.get("functional", "BP86"),
            basis=kwargs.get("basis", "def2-SVP"),
            ri=kwargs.get("ri", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            dispersion=kwargs.get("dispersion", ""),
            extra_keywords=kwargs.get("extra_keywords"),
            opt_level=kwargs.get("opt_level", "opt"),
            geom_options=kwargs.get("geom_options"),
            freq_options=kwargs.get("freq_options"),
            maxcore=kwargs.get("maxcore", 1000),
            ncores=cores,
            base_name=kwargs.get("base_name", "optfreq"),
            work_dir=work_dir,
            timeout=kwargs.get("timeout"),
        )
        error_kind = ErrorKind.NONE
        if not res.ok:
            if res.scf_converged is False or res.opt_converged is False:
                error_kind = ErrorKind.CONVERGENCE
            else:
                error_kind = ErrorKind.TOOL_FAILED
        data: dict[str, Any] = {
            "meaningful": res.meaningful,
            "final_energy": res.final_energy,
            "opt_converged": res.opt_converged,
            "freq_converged": res.freq_converged,
            "n_modes": res.n_modes,
            "n_imag": res.n_imag,
            "imag_freqs": list(res.imag_freqs),
            "lowest_freq": res.lowest_freq,
            "zpe": res.zpe,
            "enthalpy": res.enthalpy,
            "free_energy": res.free_energy,
            "character": res.character,
            "findings": list(res.findings),
            "failure_lines": list(res.failure_lines),
        }
        return self._make_result(
            self.name,
            StepStatus.SUCCESS if res.ok else StepStatus.FAILED,
            work_dir,
            start,
            output_file=res.output_file,
            data=data,
            error=res.error,
            error_kind=error_kind,
        )


register(OrcaOptOpiAdapter())
