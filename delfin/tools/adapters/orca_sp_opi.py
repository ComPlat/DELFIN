"""DELFIN's first ORCA building block on OPI: a single-point SCF/DFT step.

This module builds a single-point DFT ORCA calculation entirely with the ORCA
Python Interface (OPI, orca-pi 2.0.0): the ``!`` simple keywords and the
``%scf`` block come from OPI's input classes (``opi.input``), the geometry
from ``opi.input.structures``, and the ``.inp`` file is written by OPI's
``Calculator.write_input()``.  It is a wrapper, not a replacement: the binary
is still launched by DELFIN's existing runner (``delfin.orca.run_orca`` — the
same local-now / Slurm-ready launcher the classic ``orca_sp`` adapter uses)
and every keyword is checked against DELFIN's ORCA knowledge store
(``delfin.orca_knowledge``) before the input is written.

Two rules bind the reader:

* A value from OPI is only meaningful once the convergence flags say so.  A
  failed run still yields a HOMO-LUMO gap (from unconverged orbitals), so
  ``hl_gap_eV`` (and every other number) is returned only when
  ``terminated_normally()`` and ``scf_converged()`` are both true.
* OPI says *that* a run failed, not *why*.  The reason lives in the ``.out``
  text, so on a failed run we keep DELFIN's text inspection of ``.out`` to
  surface the lines that explain the failure.

The module adds a new adapter ``orca_sp_opi`` to the existing step registry
(``delfin.tools``); nothing in the existing modules is changed.
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

logger = get_logger(__name__)

# Simple-input keywords a plain single-point run always understands, so a
# call that omits them does not trip the knowledge-store check.
_JOB_NEUTRAL_KEYWORDS = ()


@dataclass
class OrcaSpResult:
    """Typed result of one OPI-backed single-point DFT run.

    Fields sourced from OPI are populated only after the convergence flags
    (``terminated_normally`` and ``scf_converged``) make them meaningful:
    :attr:`meaningful` summarises that state, and every numeric field is
    ``None`` when it would not be trustworthy.
    """

    ok: bool
    work_dir: Path
    output_file: Optional[Path] = None

    # Convergence flags (always read; they gate every numeric field).
    terminated_normally: Optional[bool] = None
    scf_converged: Optional[bool] = None
    # True exactly when terminated_normally and scf_converged are both true,
    # i.e. the numeric fields below are meaningful.
    meaningful: bool = False

    # Numbers returned only when meaningful.
    energy_Eh: Optional[float] = None
    s2: Optional[tuple] = None                # (expectation value, exact)
    homo_energy_eh: Optional[float] = None
    lumo_energy_eh: Optional[float] = None
    hl_gap_eV: Optional[float] = None

    # Method knowledge findings (warnings) surfaced, never silently dropped.
    findings: tuple[dict, ...] = ()
    # Lines of the .out text that explain a failure (OPI does not give them).
    failure_lines: tuple[str, ...] = ()
    error: Optional[str] = None


def _validated_keyword(error_list: list[str], word: str, field: str) -> None:
    """Append a knowledge error to ``error_list`` when ORCA does not know
    ``word`` (checked against delfin.orca_knowledge.known_keyword)."""
    from delfin.orca_knowledge import known_keyword

    if not word:
        return
    if not known_keyword(word):
        error_list.append(
            f"{field}: keyword {word!r} is not one ORCA knows "
            "(OPI catalog / ORCA parser); refused."
        )


def validate_method(
    *,
    functional: str,
    basis: str,
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    cards: Optional[Sequence[dict]] = None,
) -> tuple[list[str], list[dict]]:
    """Check every keyword and block option against DELFIN's ORCA knowledge
    store.

    Returns ``(errors, warnings)``.  *errors* enumerate keywords / block
    options ORCA does not know; a non-empty ``errors`` must stop the step
    (``run_sp_dft`` refuses to launch).  *warnings* are the findings of the
    applicable verified knowledge cards (``delfin.orca_knowledge.findings``),
    reported and never silently dropped.
    """
    errors: list[str] = []
    _validated_keyword(errors, functional, "functional")
    _validated_keyword(errors, basis, "basis")
    _validated_keyword(errors, ri, "ri")
    _validated_keyword(errors, aux_basis, "aux_basis")
    _validated_keyword(errors, dispersion, "dispersion")
    for kw in extra_keywords or ():
        _validated_keyword(errors, kw, "keyword")

    if scf_options:
        from delfin.orca_knowledge import known_block_option

        for option in scf_options:
            if not known_block_option("scf", option):
                errors.append(
                    f"scf: block option {'scf'} option {option!r} is not one ORCA "
                    "knows (OPI block model); refused."
                )

    # Findings of the applicable verified knowledge cards; empty unless a
    # card with a check applies (the shipped store has none as of this phase,
    # but the mechanism stands and any future card's warning lands here).
    calc = {
        "keywords": [
            w for w in (functional, basis, ri, aux_basis, dispersion)
            if w
        ] + [str(k) for k in (extra_keywords or ())],
        "blocks": {"scf": dict(scf_options or {})},
        "method_family": "dft",
        "job": "sp",
    }
    from delfin.orca_knowledge import findings as knowledge_findings

    warnings = knowledge_findings(calc, cards=cards)
    return errors, warnings


def _build_opi_input(
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
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    maxcore: int,
    ncores: int,
    base_name: str,
) -> Path:
    """Assemble the ORCA input with OPI's classes and write it with OPI.

    Returns the absolute path of the written ``<base_name>.inp``.  Every
    keyword is added through ``opi.input`` (``!`` simple keywords from the
    typed enums), ``%pal``/``%maxcore`` map onto ``Input.ncores``/``Input.
    memory``, and the ``%scf`` block via the ``BlockScf`` pydantic block.  The
    file is written by OPI's ``Calculator.write_input()`` — no hand-written
    string template.
    """
    def _lookup(enum_cls: type, keyword: str) -> Any:
        """Return the OPI SimpleKeyword enum member whose ORCA keyword matches
        ``keyword`` (case-blind), or the bare validated string when OPI has no
        enum member for it (e.g. some basis sets).  The keyword is never
        rewritten or dropped.
        """
        target = str(keyword).lower()
        for attr in vars(enum_cls).values():
            if isinstance(attr, SimpleKeyword) and str(attr).lower() == target:
                return attr
        return str(keyword)

    from opi.core import Calculator
    from opi.input import Input
    from opi.input.blocks.block_scf import BlockScf
    from opi.input.simple_keywords.approximation import Approximation
    from opi.input.simple_keywords.base import SimpleKeyword
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dft import Dft
    from opi.input.simple_keywords.dispersion_correction import (  # noqa: E501
        DispersionCorrection,
    )
    from opi.input.structures import Atom, Structure

    inp = Input()
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
    for kw in extra_keywords or ():
        keywords.append(str(kw))
    inp.add_simple_keywords(*keywords)
    inp.ncores = ncores
    inp.memory = maxcore
    if scf_options:
        inp.add_blocks(BlockScf(**scf_options))

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


def _failure_lines(out_path: Path) -> tuple[str, ...]:
    """The .out lines that explain a failure, so a survivor of OPI's binary
    ``terminated_normally`` verdict has the *reason* to act on."""
    if not out_path.is_file():
        return ()
    try:
        text = out_path.read_text(errors="replace")
    except OSError:
        return ()
    marks = (
        "SCF NOT CONVERGED",
        "ERROR",
        "CONVERGENCE FAILURE",
        "aborted",
        "out-of-memory",
        "not converged",
        "FATAL",
    )
    out: list[str] = []
    for line in text.splitlines():
        if any(m in line for m in marks):
            if len(out) < 40:
                out.append(line.strip())
    return tuple(out)


def read_orca_output(work_dir: Path, base_name: str) -> OrcaSpResult:
    """Read a finished run with OPI, honouring convergence-gated reporting.

    The numeric fields are filled in only when ``terminated_normally()`` and
    ``scf_converged()`` both return true; a failed run keeps them ``None`` and
    instead fills :attr:`failure_lines` from the ``.out`` text.
    """
    from opi.output.core import Output

    out_path = work_dir / f"{base_name}.out"
    try:
        out = Output(basename=base_name, working_dir=work_dir)
        out.parse()
    except Exception as exc:  # OPI could not read the run at all
        return OrcaSpResult(
            ok=False,
            work_dir=work_dir,
            output_file=out_path if out_path.is_file() else None,
            meaningful=False,
            failure_lines=_failure_lines(out_path),
            error=f"OPI could not parse the run: {exc}",
        )

    terminated = out.terminated_normally()
    scf_ok = out.scf_converged()
    meaningful = bool(terminated and scf_ok)

    result = OrcaSpResult(
        ok=meaningful,
        work_dir=work_dir,
        output_file=out_path if out_path.is_file() else None,
        terminated_normally=terminated,
        scf_converged=scf_ok,
        meaningful=meaningful,
    )
    if meaningful:
        try:
            result.energy_Eh = out.get_final_energy()
            result.s2 = tuple(out.get_s2()) if out.get_s2() is not None else None
        except Exception:
            logger.debug("OPI property read failed for %s", work_dir, exc_info=True)
        try:
            result.hl_gap_eV = out.get_hl_gap()
            homo = out.get_homo()
            lumo = out.get_lumo()
            if homo is not None:
                result.homo_energy_eh = homo.orbitalenergy
            if lumo is not None:
                result.lumo_energy_eh = lumo.orbitalenergy
        except Exception:
            logger.debug("OPI orbital read failed for %s", work_dir, exc_info=True)
    else:
        result.failure_lines = _failure_lines(out_path)
        result.error = (
            "run did not finish normally: "
            f"terminated_normally={terminated}, scf_converged={scf_ok}"
        )
    return result


def run_sp_dft(
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
    scf_options: Optional[Mapping[str, Any]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "calc",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> OrcaSpResult:
    """Run one single-point DFT calculation on OPI and read the result.

    Steps, in order:

    1. validate every keyword / block option against ``delfin.orca_knowledge``
       (errors abort before any file is written; warnings are carried into the
       result),
    2. build and write the input with OPI's classes,
    3. launch through DELFIN's ``delfin.orca.run_orca`` (not a new runner),
    4. read back with OPI, gating every number on the convergence flags.

    ``work_dir`` must already exist; when it is ``None`` a fresh directory is
    created under the current directory (mirroring ``delfin.tools.run_step``).
    """
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else _make_workdir()
    wdir.mkdir(parents=True, exist_ok=True)

    errors, findings = validate_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        scf_options=scf_options,
    )
    if errors:
        return OrcaSpResult(
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
        scf_options=scf_options,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    out_path = wdir / f"{base_name}.out"
    success = run_orca(str(inp_path), str(out_path), timeout=timeout, working_dir=wdir)
    result = read_orca_output(wdir, base_name)
    result.findings = tuple(findings)
    if not success and result.ok:
        # ORCA subprocess returned failure even though OPI read meaningful
        # values; defer to the subprocess verdict and surface the lines.
        result.ok = False
        result.meaningful = False
        result.error = result.error or "ORCA subprocess reported failure"
        result.failure_lines = result.failure_lines or _failure_lines(out_path)
    return result


def _make_workdir() -> Path:
    import tempfile

    return Path(tempfile.mkdtemp(prefix="orca_sp_opi_", dir=Path.cwd()))


# --------------------------------------------------------------------------
# Adapter for delfin.tools (registered as ``orca_sp_opi``)
# --------------------------------------------------------------------------


class OrcaSpOpiAdapter(StepAdapter):
    """A single-point DFT step driven by OPI, launchable via ``run_step``.

    Mirrors ``OrcaSpAdapter``'s position in the registry but reads and builds
    through OPI.  Values are reported only once the convergence flags make
    them meaningful; an unknown keyword / block option stops the step.
    """

    name = "orca_sp_opi"
    description = "OPI-backed ORCA single-point DFT calculation"
    category = "dft"
    produces_geometry = False
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
        res = run_sp_dft(
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            functional=kwargs.get("functional", "B3LYP"),
            basis=kwargs.get("basis", "def2-SVP"),
            ri=kwargs.get("ri", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            dispersion=kwargs.get("dispersion", ""),
            extra_keywords=kwargs.get("extra_keywords"),
            scf_options=kwargs.get("scf_options"),
            maxcore=kwargs.get("maxcore", 1000),
            ncores=cores,
            base_name=kwargs.get("base_name", "calc"),
            work_dir=work_dir,
            timeout=kwargs.get("timeout"),
        )
        return self._make_result(
            self.name,
            StepStatus.SUCCESS if res.ok else StepStatus.FAILED,
            work_dir,
            start,
            output_file=res.output_file,
            data={
                "meaningful": res.meaningful,
                "energy_Eh": res.energy_Eh,
                "s2": res.s2,
                "hl_gap_eV": res.hl_gap_eV,
                "findings": list(res.findings),
                "failure_lines": list(res.failure_lines),
            },
            error=res.error,
            # No kind on success; an unconverged SCF is a convergence failure,
            # anything else (abort, missing output) a tool failure.
            error_kind=(ErrorKind.NONE if res.ok
                        else ErrorKind.CONVERGENCE if res.scf_converged is False
                        else ErrorKind.TOOL_FAILED),
        )


register(OrcaSpOpiAdapter())
