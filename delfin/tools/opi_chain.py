"""Composability machinery for chaining OPI-backed ORCA steps.

Builds on the validated helpers in ``delfin.tools.adapters.orca_sp_opi``
(``validate_method`` / ``_build_opi_input`` / ``read_orca_output``) WITHOUT
editing that module.  The chain semantics are the Pipeline engine's own:

- **Opt geometry feeds the next step** — an OPT adapter that returns
  ``StepResult.geometry`` (the ``<base>.xyz`` ORCA writes) makes the next step
  inherit it, because ``Pipeline.run`` propagates ``result.geometry`` into the
  next step's ``geometry`` kwarg.
- **the .gbw is offered as MORead start orbitals** — an adapter that returns
  ``result.artifacts["gbw"]`` lets the wiring rule ``gbw -> moread`` inject the
  path into the next ``orca_*`` step's ``moread`` kwarg; the builder then
  writes ``! moread`` + ``%moinp "<gbw>"`` so ORCA starts from those orbitals.
- **a frequency step's thermal correction is available to an energy profile**
  — the freq adapter exposes ``data["free_energy_Eh"]`` / ``["enthalpy_Eh"]`` /
  ``["zpe_Eh"]`` (OPI ``get_free_energy`` etc.; the thermochemistry the manual
  describes in section 4.7) that a reporting / energy-profile step can read.
- **a failed or unconverged step stops the chain with its reason** — every
  adapter here returns ``StepStatus.FAILED`` with ``error`` / ``error_kind``
  when ORCA did not finish, so the pipeline's ``stop_on_failure`` halts.

The MORead block is written the way the ORCA 6.1 manual documents it —
``! moread`` plus a top-level ``%moinp "<file>.gbw"`` (section 2.20.5,
``essentialelements/initialguess.txt``) — and DELFIN's production adapter
emits the same form.  OPI's ``ArbitraryString`` provides the sanctioned escape
hatch for block syntax the typed enums express differently.
"""

from __future__ import annotations

import logging
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from delfin.tools.adapters.orca_sp_opi import (  # noqa: F401  (re-exported)
    _build_opi_input,
    _make_workdir,
    read_orca_output,
    validate_method,
)
from delfin.tools._types import StepResult

logger = logging.getLogger(__name__)

JOB_OPT = "OPT"
JOB_FREQ = "FREQ"
JOB_OPT_FREQ = "OPT_FREQ"
_JOB_KINDS = frozenset({JOB_OPT, JOB_FREQ, JOB_OPT_FREQ})

#: ``%freq`` block options the chain steps may carry, matched against OPI's
#: BlockFreq fields and validated against DELFIN's knowledge store.
_FREQ_OPTION_KEYS = ("temp", "pressure", "numfreq", "anfreq", "quasirrho",
                     "cutofffreq", "qrrhoreffreq", "scalfreq")


@dataclass
class OpiChainResult:
    """Outcome of one OPI-backed ORCA step in a chain, convergence-gated.

    Numbers are filled only when OPI reports a normal termination *and* a
    converged SCF; a failed run keeps them ``None`` and carries the failure.
    """

    ok: bool
    work_dir: Path
    meaningful: bool = False
    terminated_normally: bool = False
    scf_converged: bool = False
    energy_Eh: Optional[float] = None
    thermal: dict = field(default_factory=dict)
    geometry_path: Optional[Path] = None
    artifacts: dict = field(default_factory=dict)
    failure_lines: tuple = ()
    error: Optional[str] = None
    output_file: Optional[Path] = None
    geometry_optimization_converged: bool = False


def validate_chain_method(
    *,
    functional: str,
    basis: str,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    freq_options: Optional[Mapping[str, Any]] = None,
    cards: Optional[Sequence[dict]] = None,
) -> tuple[list[str], list[dict]]:
    """Validate a chain step's keywords / block options against the knowledge
    store, reusing ``validate_method`` for the DFT surface and additionally
    checking each ``%freq`` block option is one ORCA knows."""
    errors, warnings = validate_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        scf_options=scf_options,
        cards=cards,
    )
    for key in (freq_options or {}):
        if key not in _FREQ_OPTION_KEYS:
            from delfin.orca_knowledge import known_block_option

            if not known_block_option("freq", key):
                errors.append(
                    f"freq: block option {key!r} is not one ORCA knows; refused."
                )
    return errors, warnings


def build_orca_input(
    work_dir: Path,
    *,
    job_kind: str,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str = "",
    basis: str,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    method_keywords: Optional[Sequence[str]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    freq_options: Optional[Mapping[str, Any]] = None,
    moread: Optional[str] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "calc",
) -> Path:
    """Assemble and write the ORCA input with OPI's classes for an *OPT*,
    *FREQ* or *OPT+FREQ* job, mirroring ``_build_opi_input``'s writer.

    ``! Opt`` / ``! Freq`` are added as OPI simple keywords (``Opt.OPT``,
    ``Task.FREQ``), the ``%freq`` block via ``BlockFreq``.  An optional
    MOREAD *gbw* writes ``! moread`` + ``%moinp "<gbw>"`` (the manual's form,
    section 2.20.5), the gbw being copied into the working directory first so
    ORCA finds it by name.  The file is written by OPI's
    ``Calculator.write_input``.
    """
    from opi.core import Calculator
    from opi.input import ArbitraryString, ArbitraryStringPos, Input
    from opi.input.blocks.block_freq import BlockFreq
    from opi.input.simple_keywords.approximation import Approximation
    from opi.input.simple_keywords.base import SimpleKeyword
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dft import Dft
    from opi.input.simple_keywords.dispersion_correction import (  # noqa: E501
        DispersionCorrection,
    )
    from opi.input.simple_keywords.opt import Opt
    from opi.input.simple_keywords.scf import Scf
    from opi.input.simple_keywords.task import Task
    from opi.input.structures import Atom, Structure

    def _lookup(enum_cls: type, keyword: str) -> Any:
        target = str(keyword).lower()
        for attr in vars(enum_cls).values():
            if isinstance(attr, SimpleKeyword) and str(attr).lower() == target:
                return attr
        return str(keyword)

    inp = Input()
    keywords: list[Any] = []
    if functional:
        keywords.append(_lookup(Dft, functional))
    if basis:
        keywords.append(_lookup(BasisSet, basis))
    if ri:
        keywords.append(_lookup(Approximation, ri))
    if aux_basis:
        keywords.append(_lookup(BasisSet, aux_basis))
    if dispersion:
        keywords.append(_lookup(DispersionCorrection, dispersion))
    for m in method_keywords or ():
        keywords.append(str(m))
    if job_kind in (JOB_OPT, JOB_OPT_FREQ):
        keywords.append(Opt.OPT)
    if job_kind in (JOB_FREQ, JOB_OPT_FREQ):
        keywords.append(Task.FREQ)
    for kw in extra_keywords or ():
        keywords.append(str(kw))
    if moread:
        keywords.append(Scf.MOREAD)
    inp.add_simple_keywords(*keywords)
    inp.ncores = ncores
    inp.memory = maxcore

    blocks = []
    if scf_options:
        from opi.input.blocks.block_scf import BlockScf

        blocks.append(BlockScf(**scf_options))
    if job_kind in (JOB_FREQ, JOB_OPT_FREQ) and freq_options:
        blocks.append(BlockFreq(**freq_options))
    unwritten: list[Any] = []
    if moread:
        # ORCA reads the guess orbitals from a .gbw named in %moinp (manual
        # 2.20.5); copy it into the working dir so ORCA finds it by name.
        moread_src = Path(moread)
        src = moread_src if moread_src.is_file() else moread_src.with_suffix(".gbw")
        if src.is_file():
            shutil.copy2(src, work_dir / f"{base_name}.gbw")
            unwritten.append(ArbitraryString(
                f'%moinp "{base_name}.gbw"', ArbitraryStringPos.BEFORE_COORDS))
        else:
            raise FileNotFoundError(
                f"MOREAD gbw {src} does not exist (refused to write a dangling "
                "MORead reference)")
    if blocks:
        inp.add_blocks(*blocks)
    if unwritten:
        for s in unwritten:
            inp.add_arbitrary_string(s.string, pos=s.pos)

    atoms = [
        Atom(str(el), coordinates=(float(x), float(y), float(z)))
        for el, (x, y, z) in zip(elements, coordinates)
    ]
    structure = Structure(atoms, charge=charge, multiplicity=mult)

    calc = Calculator(base_name, working_dir=work_dir, version_check=False)
    calc.input = inp
    calc.structure = structure
    calc.write_input(force=True)
    return work_dir / f"{base_name}.inp"


def collect_artifacts(work_dir: Path, base_name: str) -> dict:
    """The artifacts ORCA leaves behind in *work_dir* for a finished job
    (gbw / hess / densities), mirroring the classic ``orca.py`` collector."""
    artifacts: dict = {}
    for ext, key in ((".gbw", "gbw"), (".hess", "hess"),
                     (".densitiesinfo", "densities")):
        f = work_dir / f"{base_name}{ext}"
        if f.is_file():
            artifacts[key] = f
    return artifacts


def read_chain_output(
    work_dir: Path,
    base_name: str,
    *,
    job_kind: str,
) -> OpiChainResult:
    """Read a finished OPT/FREQ chain step with OPI, convergence-gated.

    Reuses ``read_orca_output`` for the core (terminated / scf / energy) and
    additionally pulls the thermal corrections (freq) and the optimized
    geometry file (OPT), honouring the same "no number before convergence"
    rule.
    """
    core = read_orca_output(work_dir, base_name)
    out_path = work_dir / f"{base_name}.out"
    result = OpiChainResult(
        ok=core.ok,
        work_dir=work_dir,
        meaningful=core.meaningful,
        terminated_normally=bool(core.terminated_normally),
        scf_converged=bool(core.scf_converged),
        energy_Eh=core.energy_Eh,
        failure_lines=tuple(core.failure_lines),
        error=core.error,
        output_file=out_path if out_path.is_file() else None,
        artifacts=collect_artifacts(work_dir, base_name),
    )
    # The optimized geometry ORCA reports for an OPT is the <base>.xyz file it
    # writes; this is what the classic orca.py adapter feeds the next step.
    if job_kind in (JOB_OPT, JOB_OPT_FREQ):
        xyz = work_dir / f"{base_name}.xyz"
        if xyz.is_file():
            result.geometry_path = xyz
    if core.meaningful:
        try:
            from opi.output.core import Output

            if job_kind in (JOB_OPT, JOB_OPT_FREQ):
                _out = Output(basename=base_name, working_dir=work_dir)
                _out.parse()
                result.geometry_optimization_converged = bool(
                    _out.geometry_optimization_converged())
            if job_kind in (JOB_FREQ, JOB_OPT_FREQ):
                # Thermal corrections the freq step carries (thermochemistry,
                # manual 4.7) and an energy profile can read downstream.
                try:
                    from opi.output.core import Output as _O

                    _t = Output(basename=base_name, working_dir=work_dir)
                    _t.parse()
                    for key, fn in (
                        ("free_energy_Eh", _t.get_free_energy),
                        ("enthalpy_Eh", _t.get_enthalpy),
                        ("zpe_Eh", _t.get_zpe),
                        ("entropy", _t.get_entropy),
                    ):
                        try:
                            val = fn()
                            if val is not None:
                                result.thermal[key] = float(val)
                        except Exception:  # OPI read quirk; keep the rest
                            logger.debug("thermal read failed for %s", key,
                                         exc_info=True)
                except Exception:
                    logger.debug("OPI thermal read failed for %s", work_dir,
                                 exc_info=True)
        except Exception:
            logger.debug("OPI geometry read failed for %s", work_dir,
                         exc_info=True)
    return result


def run_opi_job(
    *,
    job_kind: str,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str = "",
    basis: str,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    method_keywords: Optional[Sequence[str]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    freq_options: Optional[Mapping[str, Any]] = None,
    moread: Optional[str] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "calc",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> OpiChainResult:
    """Run one OPT / FREQ / OPT+FREQ ORCA job on OPI and return its result.

    Steps, in order: validate against the knowledge store (abort before any
    file is written), build the input with OPI, launch through
    ``delfin.orca.run_orca``, read back with OPI gating every number on the
    convergence flags (``read_chain_output``).
    """
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else _make_workdir()
    wdir.mkdir(parents=True, exist_ok=True)

    errors, _findings = validate_chain_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        scf_options=scf_options,
        freq_options=freq_options,
    )
    if errors:
        return OpiChainResult(
            ok=False,
            work_dir=wdir,
            meaningful=False,
            error="; ".join(errors),
        )

    inp_path = build_orca_input(
        wdir,
        job_kind=job_kind,
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        method_keywords=method_keywords,
        extra_keywords=extra_keywords,
        scf_options=scf_options,
        freq_options=freq_options,
        moread=moread,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    out_path = wdir / f"{base_name}.out"
    subprocess_ok = run_orca(
        str(inp_path), str(out_path), timeout=timeout, working_dir=wdir)
    result = read_chain_output(wdir, base_name, job_kind=job_kind)
    if not subprocess_ok and result.ok:
        # ORCA subprocess returned failure even though OPI read meaningful
        # values; defer to the subprocess verdict, as run_sp_dft does.
        result.ok = False
        result.meaningful = False
        result.error = result.error or "ORCA subprocess reported failure"
        result.failure_lines = result.failure_lines or ()
    return result


def step_result_from(
    adapter: Any,
    res: OpiChainResult,
    work_dir: Path,
    start: float,
    *,
    data_extra: Optional[dict] = None,
) -> StepResult:
    """Assemble a chain ``StepResult`` from an :class:`OpiChainResult`.

    Puts the geometry / artifacts / thermal data onto the result so the
    pipeline can propagate them, and marks a failed or unconverged step
    ``FAILED`` with its reason (which halts ``stop_on_failure``).
    """
    from delfin.tools._types import ErrorKind, StepStatus

    extra = dict(data_extra or {})
    data: dict = {
        "meaningful": res.meaningful,
        "energy_Eh": res.energy_Eh,
        "thermal": dict(res.thermal),
        "geometry_optimization_converged": res.geometry_optimization_converged,
        "failure_lines": list(res.failure_lines),
    }
    data.update(extra)
    if res.error:
        data["error"] = res.error
    ok = res.ok and (not res.scf_converged or res.meaningful)
    error_kind = (ErrorKind.NONE if ok
                  else ErrorKind.CONVERGENCE if res.scf_converged is False
                  else ErrorKind.TOOL_FAILED)
    return StepResult(
        step_name=getattr(adapter, "name", "?"),
        status=(StepStatus.SUCCESS if ok else StepStatus.FAILED),
        work_dir=work_dir,
        elapsed_seconds=time.monotonic() - start,
        geometry=res.geometry_path,
        artifacts=dict(res.artifacts),
        data=data,
        output_file=res.output_file,
        error=res.error,
        error_kind=error_kind,
    )


def read_geometry_xyz(path) -> tuple[list[str], list[list[float]]]:
    """Parse an ORCA/XYZ ``.xyz`` file into ``(elements, coordinates)``.

    The chain passes the previous step's geometry as a ``.xyz`` kwarg (ORCA
    writes ``<base>.xyz`` for an OPT).  Format: line 1 = atom count, line 2 =
    title, then ``symbol x y z`` per atom.  Returns the atom list a chained
    executor needs to rebuild the OPI structure.
    """
    # Drop only *leading* blank lines; the title line (line 2 of the XYZ
    # format) is allowed to be blank and must keep its position so the atom
    # block is indexed as line 2 .. 2+natoms, not shifted by the filter.
    lines = Path(path).read_text().splitlines()
    while lines and not lines[0].strip():
        lines.pop(0)
    if len(lines) < 3:
        raise ValueError(f"{path} is not a three-line-plus XYZ file")
    try:
        natoms = int(lines[0].split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path}: first line must be the atom count") from exc
    # After the count+title (which may be blank), take the next `natoms`
    # non-blank lines as the atom block.
    atom_candidates = [ln for ln in lines[2:] if ln.strip()]
    atoms = atom_candidates[:natoms]
    if len(atoms) < natoms:
        raise ValueError(
            f"{path}: expected {natoms} atoms after the title, found {len(atoms)}")
    elements: list[str] = []
    coordinates: list[list[float]] = []
    for ln in atoms:
        parts = ln.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: malformed atom line {ln!r}")
        elements.append(parts[0])
        coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elements, coordinates
