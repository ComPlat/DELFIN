"""DELFIN's DLPNO-CCSD(T) single-point building block on OPI.

This module builds a local-correlation single point entirely with the ORCA
Python Interface (OPI, orca-pi 2.0.0): the method keyword, the orbital and
auxiliary basis, the PNO level and the ``%mdci`` block come from OPI's
typed input classes, the input is written by OPI's
``Calculator.write_input()``, the binary is launched by DELFIN's existing
runner (``delfin.orca.run_orca``) and every keyword is checked against
DELFIN's ORCA knowledge store before anything is written.
``orca_sp_opi.py`` helpers are imported, never edited.

Manual ground truth (ORCA 6.1 manual, modelchemistries/mdci.txt):

* Section 3.10.10.1 "Basic Usage": the example
  ``! cc-pVTZ cc-pVTZ/C DLPNO-CCSD(T) TightPNO TightSCF`` -- the auxiliary
  basis is the ``/C`` fitting basis of the *same family* as the orbital
  basis (cc-pVnZ -> cc-pVnZ/C, def2-XZVPP -> def2-XZVPP/C).  RIJCOSX or
  RIJK may be given for the SCF reference (section 3.10.10: scaling
  behaviour of DLPNO-CCSD(T) "using RIJCOSX"; the /C fitting basis is the
  pair-correlation one, not the MP2/JK auxiliary).
* PNO level: "Sensible defaults have been determined from extensive
  benchmark calculations and are accessible via LoosePNO, NormalPNO and
  TightPNO keywords in the simple input line.  These keywords represent
  the recommended way to control the accuracy of DLPNO calculations."
  The TightPNO example is annotated: "# Tight settings for increased
  accuracy, e.g. when investigating weak interactions or conformational
  equilibria" -- hence this module's *suggested* default
  DLPNO-CCSD(T)/TightPNO.  NormalPNO is ORCA's default and need not be
  spelled out.
* Two-point CBS extrapolation: "Using the keyword
  ``! Extrapolate(X/Y,basis)``, where X and Y are the corresponding
  successive cardinal numbers and basis is the type of basis set
  requested (= cc, aug-cc, cc-core, ano, saug-ano, aug-ano, def2) ORCA
  will calculate the SCF and optionally the MP2 or MDCI energies with two
  basis sets and separately extrapolate."  Printed exponents (cc family):
  "Alpha(2/3): 4.420, Beta(2/3): 2.460; Alpha(3/4): 5.460,
  Beta(3/4): 3.050"; def2-family table: alpha_23 = 10.39, beta_23 = 2.40,
  alpha_34 = 7.88, beta_34 = 2.97.  "Since the beta values for 2/3 are
  close to 2.4, we always take this value.  Likewise, all 3/4 and higher
  extrapolations are done with beta=3."  "The extrapolation method should
  only be used with very tight SCF convergence criteria"; "Auxiliary
  basis sets are automatically chosen and can not be changed"; "only
  implemented for single points".  The generalization
  ``ExtrapolateEP2(n/m,bas,[method,method-details])`` runs DLPNO-CCSD(T)
  as the cheap method (default) with LoosePNO/NormalPNO/TightPNO as
  method-details.
* DLPNO-CCSD(T1): "an improved perturbative triples correction for
  DLPNO-CCSD ... the triples amplitudes are computed iteratively, which
  can reproduce more accurately the canonical (T) energies" (available
  since ORCA 4.1).  Nomenclature: the ORCA keyword DLPNO-CCSD(T) is the
  semi-canonical correction, DLPNO-CCSD(T1) the iterative one.
* T1 diagnostic readback: "A rule of thumb says, that for a value of the
  diagnostic of larger than 0.02 the results are not to be trusted."
  (The manual also argues the T1 diagnostic measures orbital relaxation,
  not multi-reference character, and offers the S-diagnostic for
  canonical CCSD via ``%mdci DoSDiagnostic true``.)

Reader rules, inherited from the orca_sp_opi contract: numbers are
returned only when the run ``terminated_normally()`` and the SCF
converged; a failed run keeps them ``None`` and carries the ``.out``
failure lines instead.  The DLPNO-CCSD(T) energy is read through OPI
(``get_energies()`` key "MDCI(SD(T))"); the T1 diagnostic is parsed from
the ``.out`` text, where ORCA prints it (T1(DLPNO) line).
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Sequence

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import ErrorKind, StepResult, StepStatus

logger = get_logger(__name__)

#: families for the /C auxiliary basis: orbital-basis regex -> aux format
_FAMILY_PATTERNS: tuple[tuple[str, str], ...] = (
    (r"^cc-pV([DTQ])Z$", r"cc-pV\1Z/C"),
    (r"^aug-cc-pV([DTQ])Z$", r"aug-cc-pV\1Z/C"),
    (r"^def2-(S|TZ|Q)VPP$", r"def2-\1VPP/C"),
    (r"^def2-(S|TZ|Q)V[PZ]?$", r"def2-\1V/C"),
)

#: CBS pairs -> the manual's Extrapolate(X/Y,basis) family keyword
_CBS_FAMILIES: dict[tuple[str, str], str] = {
    ("cc-pVTZ", "cc-pVQZ"): "Extrapolate(3/4,cc)",
    ("cc-pVDZ", "cc-pVTZ"): "Extrapolate(2/3,cc)",
    ("def2-TZVPP", "def2-QZVPP"): "Extrapolate(3/4,def2)",
}

#: cardinal number of the +VnZ orbital bases, for the successive check
_CARDINALS = {"D": 2, "T": 3, "Q": 4}

#: Simple keywords a DLPNO single point always understands
_PNO_LEVELS = ("LoosePNO", "NormalPNO", "TightPNO")
_DLPNO_METHODS = ("DLPNO-CCSD(T)", "DLPNO-CCSD(T1)", "DLPNO-CCSD")
#: SCF reference approximations that may accompany DLPNO (manual 3.10.10)
_DLPNO_RI = ("RIJCOSX", "RIJK")


@dataclass
class DlpnoResult:
    """Typed result of one OPI-backed DLPNO-CCSD(T) single point.

    Numeric fields are populated only after the run
    ``terminated_normally()`` and ``scf_converged()`` are both true.
    """

    ok: bool
    work_dir: Path
    output_file: Optional[Path] = None

    terminated_normally: Optional[bool] = None
    scf_converged: Optional[bool] = None
    meaningful: bool = False

    # Numbers returned only when meaningful.
    scf_energy_Eh: Optional[float] = None
    ccsd_energy_Eh: Optional[float] = None     # MDCI(SD)
    ccsd_t_energy_Eh: Optional[float] = None   # MDCI(SD(T))
    energy_Eh: Optional[float] = None          # final (== the CCSD(T) energy)
    t1_diagnostic: Optional[float] = None
    triples_correction: str = ""               # "(T)", "(T1)" or ""
    aux_basis_used: str = ""

    findings: tuple[dict, ...] = ()
    failure_lines: tuple[str, ...] = ()
    error: Optional[str] = None


def _aux_for(basis: str) -> Optional[str]:
    """The /C auxiliary basis of the orbital basis's family, or None when
    the family is unknown (manual 3.10.10.1: cc-pVTZ -> cc-pVTZ/C)."""
    b = basis.strip()
    for pattern, fmt in _FAMILY_PATTERNS:
        m = re.match(pattern, b)
        if m:
            return re.sub(r"\\1", m.group(1), fmt)
    return None


def suggested_default() -> dict:
    """The module's *suggested* default, with the manual's reason.

    DLPNO-CCSD(T)/TightPNO: the manual's TightPNO example is annotated
    "Tight settings for increased accuracy, e.g. when investigating weak
    interactions or conformational equilibria" (mdci.txt, section
    3.10.10.1).  This is a suggestion the user confirms; an explicit user
    choice is only validated, never overridden.
    """
    return {
        "method": "DLPNO-CCSD(T)",
        "pno": "TightPNO",
        "scf": "TightSCF",
        "reason": (
            "Manual 6.1, mdci.txt section 3.10.10.1: 'Tight settings for "
            "increased accuracy, e.g. when investigating weak interactions "
            "or conformational equilibria'; 'These keywords represent the "
            "recommended way to control the accuracy of DLPNO calculations.'"
        ),
    }


def validate_dlpno(
    *,
    method: str = "DLPNO-CCSD(T)",
    pno: str = "TightPNO",
    scf: str = "TightSCF",
    basis: str = "",
    aux_basis: str = "",
    ri: str = "",
    cbs: Optional[tuple[str, str]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    mdci_options: Optional[dict] = None,
    cards: Optional[Sequence[dict]] = None,
) -> tuple[list[str], list[dict]]:
    """Check every keyword and option against DELFIN's ORCA knowledge store.

    ``cbs`` requests a two-point CBS extrapolation; the pair must be two
    successive cardinal numbers of one basis family (cc-pVTZ/cc-pVQZ or
    def2-TZVPP/def2-QZVPP), because the manual's Extrapolate(X/Y,basis)
    takes "the corresponding successive cardinal numbers" of one basis
    type.  When ``cbs`` is given, ``basis`` must be empty (the automatic
    procedure chooses both basis sets; "All other basis set choices ...
    will be ignored").

    Returns ``(errors, warnings)``; a non-empty ``errors`` must stop the
    step.
    """
    from delfin.tools.adapters.orca_sp_opi import _validated_keyword

    errors: list[str] = []

    if method not in _DLPNO_METHODS:
        errors.append(
            f"method: {method!r} is not one of {', '.join(_DLPNO_METHODS)} "
            "(mdci.txt 3.10.10); refused."
        )
    if pno and pno not in _PNO_LEVELS:
        errors.append(
            f"pno: {pno!r} is not one of {', '.join(_PNO_LEVELS)} "
            "(mdci.txt 3.10.10.1); refused."
        )
    _validated_keyword(errors, scf, "scf")
    _validated_keyword(errors, ri, "ri")
    if ri and ri not in _DLPNO_RI:
        errors.append(
            f"ri: {ri!r} is not one of {', '.join(_DLPNO_RI)} for a DLPNO "
            "reference (mdci.txt 3.10.10); refused."
        )
    for kw in extra_keywords or ():
        _validated_keyword(errors, kw, "keyword")

    if cbs:
        lo, hi = (str(x).strip() for x in cbs)
        family = _CBS_FAMILIES.get((lo, hi))
        if family is None:
            errors.append(
                f"cbs: ({lo!r}, {hi!r}) is not a supported two-point pair; "
                "the manual's Extrapolate(X/Y,basis) needs successive "
                "cardinal numbers of one family, e.g. (cc-pVTZ, cc-pVQZ) "
                "or (def2-TZVPP, def2-QZVPP) (mdci.txt 3.10.8); refused."
            )
        if basis:
            errors.append(
                "cbs: give no orbital basis with a CBS extrapolation -- the "
                "manual: 'All other basis set choices, externally defined "
                "bases etc. will be ignored -- the automatic procedure only "
                "works with the default basis sets!' (mdci.txt 3.10.8); "
                "refused."
            )
        if aux_basis:
            errors.append(
                "cbs: give no auxiliary basis with a CBS extrapolation -- "
                "the manual: 'Auxiliary basis sets are automatically chosen "
                "and can not be changed' (mdci.txt 3.10.8); refused."
            )
    else:
        if not basis:
            errors.append("basis: an orbital basis is required for a DLPNO "
                          "single point; refused.")
        elif not aux_basis:
            aux = _aux_for(basis)
            if aux is None:
                errors.append(
                    f"basis: {basis!r} has no known /C auxiliary-basis "
                    "family; give aux_basis explicitly (mdci.txt "
                    "3.10.10.1); refused."
                )
        else:
            want = _aux_for(basis)
            if want is not None and aux_basis != want:
                errors.append(
                    f"aux_basis: {basis!r} pairs with {want!r} (same "
                    "family), not {0!r} (mdci.txt 3.10.10.1); refused."
                    .format(aux_basis)
                )

    if mdci_options:
        from delfin.orca_knowledge import known_block_option

        for option in mdci_options:
            if not known_block_option("mdci", option):
                errors.append(
                    f"mdci: block option {option!r} is not one ORCA knows "
                    "(OPI block model / ORCA parser); refused."
                )

    keywords = [w for w in (method, pno, scf, ri, basis, aux_basis) if w]
    calc = {
        "keywords": keywords + [str(k) for k in (extra_keywords or ())]
        + ([family] if cbs and (family := _CBS_FAMILIES.get(
            (str(cbs[0]).strip(), str(cbs[1]).strip()))) else []),
        "blocks": {"mdci": dict(mdci_options or {})},
        "method_family": "dlpno",
        "job": "sp",
    }
    from delfin.orca_knowledge import findings as knowledge_findings

    warnings = knowledge_findings(calc, cards=cards)
    return errors, warnings


def build_dlpno_input(
    work_dir: Path,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int = 0,
    mult: int = 1,
    method: str = "DLPNO-CCSD(T)",
    pno: str = "TightPNO",
    scf: str = "TightSCF",
    basis: str = "",
    aux_basis: str = "",
    ri: str = "",
    cbs: Optional[tuple[str, str]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    mdci_options: Optional[dict] = None,
    maxcore: int = 2000,
    ncores: int = 2,
    base_name: str = "dlpno",
) -> Path:
    """Assemble and write the DLPNO ORCA input with OPI's classes.

    Validation runs first; an unknown keyword or option raises
    ``ValueError`` before any file is written.  Returns the path of the
    written ``<base_name>.inp``.
    """
    errors, _ = validate_dlpno(
        method=method,
        pno=pno,
        scf=scf,
        basis=basis,
        aux_basis=aux_basis,
        ri=ri,
        cbs=cbs,
        extra_keywords=extra_keywords,
        mdci_options=mdci_options,
    )
    if errors:
        raise ValueError("; ".join(errors))

    from opi.core import Calculator
    from opi.input import Input
    from opi.input.blocks.block_mdci import BlockMdci
    from opi.input.simple_keywords.approximation import Approximation
    from opi.input.simple_keywords.base import SimpleKeyword
    from opi.input.simple_keywords.aux_basis_set import AuxBasisSet
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dlpno import Dlpno
    from opi.input.simple_keywords.wft import Wft
    from opi.input.structures import Atom, Structure

    def _lookup(enum_cls: type, keyword: str) -> Any:
        target = str(keyword).lower()
        for attr in vars(enum_cls).values():
            if isinstance(attr, SimpleKeyword) and str(attr).lower() == target:
                return attr
        return str(keyword)

    inp = Input()
    if cbs:
        lo, hi = (str(x).strip() for x in cbs)
        inp.add_simple_keywords(_CBS_FAMILIES[(lo, hi)])
    else:
        if not aux_basis:
            # The /C auxiliary basis of the same family (mdci.txt
            # 3.10.10.1); validate_dlpno already refused an unknown family.
            aux_basis = _aux_for(basis) or ""
        inp.add_simple_keywords(_lookup(BasisSet, basis))
        inp.add_simple_keywords(_lookup(AuxBasisSet, aux_basis))
    inp.add_simple_keywords(_lookup(Wft, method))
    if pno and pno != "NormalPNO":
        # NormalPNO is ORCA's default ("no need to give NormalPNO
        # explicitly", mdci.txt 3.10.10.1); the others are spelled out.
        inp.add_simple_keywords(_lookup(Dlpno, pno))
    inp.add_simple_keywords(str(scf))
    if ri:
        inp.add_simple_keywords(_lookup(Approximation, ri))
    for kw in extra_keywords or ():
        inp.add_simple_keywords(str(kw))
    inp.ncores = ncores
    inp.memory = maxcore
    if mdci_options:
        inp.add_blocks(BlockMdci(**dict(mdci_options)))

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


#: where ORCA prints the T1 diagnostic in a DLPNO run
_T1_RE = re.compile(r"^T1\(DLPNO\)\s*:\s*([0-9.eE+-]+)\s*$", re.MULTILINE)


def _read_t1_diagnostic(out_path: Path) -> Optional[float]:
    """Parse the T1 diagnostic from the .out text (manual: 'a value of the
    diagnostic of larger than 0.02' -- the number ORCA prints)."""
    if not out_path.is_file():
        return None
    try:
        text = out_path.read_text(errors="replace")
    except OSError:
        return None
    match = _T1_RE.search(text)
    return float(match.group(1)) if match else None


def _failure_lines(out_path: Path) -> tuple[str, ...]:
    """The .out lines that explain a failure (same marks as orca_sp_opi)."""
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
    return tuple(
        line.strip()
        for line in text.splitlines()
        if any(m in line for m in marks)
    )[:40]


def read_dlpno_output(work_dir: Path, base_name: str) -> DlpnoResult:
    """Read a finished DLPNO run with OPI, gating every number on the
    convergence flags."""
    from opi.output.core import Output

    out_path = work_dir / f"{base_name}.out"
    try:
        out = Output(basename=base_name, working_dir=work_dir)
        out.parse()
    except Exception as exc:
        return DlpnoResult(
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

    result = DlpnoResult(
        ok=meaningful,
        work_dir=work_dir,
        output_file=out_path if out_path.is_file() else None,
        terminated_normally=terminated,
        scf_converged=scf_ok,
        meaningful=meaningful,
    )
    if not meaningful:
        result.failure_lines = _failure_lines(out_path)
        result.error = (
            "run did not finish normally: "
            f"terminated_normally={terminated}, scf_converged={scf_ok}"
        )
        return result

    energies = out.get_energies() or {}
    for key, energy in energies.items():
        base_key = key.split("_")[0]
        # OPI 2.0.0's Energy model stores the number in ``totalenergy``
        # (list[list[float]]); take the last entry ( opi/output/models/json/
        # property/properties/energy.py, class Energy).
        raw = getattr(energy, "totalenergy", None) or getattr(energy, "value", None)
        if not raw:
            continue
        value = float(raw[-1][-1]) if isinstance(raw[0], (list, tuple)) else float(raw)
        if value is None:
            continue
        if base_key == "SCF":
            result.scf_energy_Eh = float(value)
        elif base_key == "MDCI(SD)":
            result.ccsd_energy_Eh = float(value)
        elif base_key == "MDCI(SD(T))":
            result.ccsd_t_energy_Eh = float(value)
    result.energy_Eh = (
        float(out.get_final_energy()) if out.get_final_energy() is not None else None
    )
    result.t1_diagnostic = _read_t1_diagnostic(out_path)
    # Which triples correction ran: from the printed T1 marker (DLPNO runs
    # always print T1(DLPNO)) plus the input file's method keyword.
    try:
        inp_text = (work_dir / f"{base_name}.inp").read_text(errors="replace")
    except OSError:
        inp_text = ""
    lowered = inp_text.lower()
    if "dlpno-ccsd(t1)" in lowered:
        result.triples_correction = "(T1)"
    elif "dlpno-ccsd(t)" in lowered:
        result.triples_correction = "(T)"
    if result.ccsd_t_energy_Eh is None and result.ccsd_energy_Eh is not None:
        # A plain DLPNO-CCSD run has no (T) energy; the SD energy is the
        # correlation result.
        result.energy_Eh = result.energy_Eh or result.ccsd_energy_Eh
    result.aux_basis_used = ""
    match = re.search(r"NewAuxJGTO\s+(\S+)", inp_text) or re.search(
        r"([A-Za-z0-9-]+/[Cc])\b", inp_text
    )
    if match:
        result.aux_basis_used = match.group(1).strip('"').replace("/c", "/C")
        # normalize OPI's lowercase keyword output back to the canonical
        # basis name (cc-pvtz/C -> cc-pVTZ/C); kept as a small map because
        # ORCA itself is case-insensitive for basis names
        canonical = {"cc-pvtz/c": "cc-pVTZ/C", "cc-pvqz/c": "cc-pVQZ/C",
                     "def2-tzvpp/c": "def2-TZVPP/C", "def2-qzvpp/c": "def2-QZVPP/C"}
        result.aux_basis_used = canonical.get(
            result.aux_basis_used.lower(), result.aux_basis_used
        )
    return result


def run_dlpno(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int = 0,
    mult: int = 1,
    method: str = "DLPNO-CCSD(T)",
    pno: str = "TightPNO",
    scf: str = "TightSCF",
    basis: str = "",
    aux_basis: str = "",
    ri: str = "",
    cbs: Optional[tuple[str, str]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    mdci_options: Optional[dict] = None,
    maxcore: int = 2000,
    ncores: int = 2,
    base_name: str = "dlpno",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> DlpnoResult:
    """Run one DLPNO-CCSD(T) single point and read the result back.

    Steps, in order: validate against the knowledge store (errors abort
    before any file is written), build and write the input with OPI,
    launch through DELFIN's ``delfin.orca.run_orca``, read back with OPI
    with every number gated on the convergence flags.
    """
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else Path.cwd() / "dlpno_run"
    wdir.mkdir(parents=True, exist_ok=True)

    errors, findings = validate_dlpno(
        method=method,
        pno=pno,
        scf=scf,
        basis=basis,
        aux_basis=aux_basis,
        ri=ri,
        cbs=cbs,
        extra_keywords=extra_keywords,
        mdci_options=mdci_options,
    )
    if errors:
        return DlpnoResult(
            ok=False,
            work_dir=wdir,
            meaningful=False,
            findings=tuple(findings),
            error="; ".join(errors),
        )

    inp_path = build_dlpno_input(
        wdir,
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        method=method,
        pno=pno,
        scf=scf,
        basis=basis,
        aux_basis=aux_basis,
        ri=ri,
        cbs=cbs,
        extra_keywords=extra_keywords,
        mdci_options=mdci_options,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    out_path = wdir / f"{base_name}.out"
    success = run_orca(str(inp_path), str(out_path), timeout=timeout, working_dir=wdir)
    result = read_dlpno_output(wdir, base_name)
    result.findings = tuple(findings)
    if not success and result.ok:
        result.ok = False
        result.meaningful = False
        result.error = result.error or "ORCA subprocess reported failure"
        result.failure_lines = result.failure_lines or _failure_lines(out_path)
    return result


# --------------------------------------------------------------------------
# Adapter for delfin.tools (registered as ``orca_dlpno_opi``)
# --------------------------------------------------------------------------


class OrcaDlpnoOpiAdapter(StepAdapter):
    """A DLPNO-CCSD(T) single-point step driven by OPI, launchable via
    ``run_step``.  Mirrors ``OrcaSpOpiAdapter``'s position in the registry.

    Suggested default: DLPNO-CCSD(T)/TightPNO + TightSCF -- see
    :func:`suggested_default` for the manual's reason.  The caller's
    explicit choice wins and is only validated.
    """

    name = "orca_dlpno_opi"
    description = "OPI-backed ORCA DLPNO-CCSD(T)/(T1) single point"
    category = "wft"
    produces_geometry = False
    consumes = ("geometry",)
    produces = ("gbw", "qc_output")
    requires_python = ("opi",)

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 2,
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
        res = run_dlpno(
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            method=kwargs.get("method", "DLPNO-CCSD(T)"),
            pno=kwargs.get("pno", "TightPNO"),
            scf=kwargs.get("scf", "TightSCF"),
            basis=kwargs.get("basis", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            ri=kwargs.get("ri", ""),
            cbs=kwargs.get("cbs"),
            extra_keywords=kwargs.get("extra_keywords"),
            mdci_options=kwargs.get("mdci_options"),
            maxcore=kwargs.get("maxcore", 2000),
            ncores=cores,
            base_name=kwargs.get("base_name", "dlpno"),
            work_dir=work_dir,
            timeout=kwargs.get("timeout"),
        )
        data: dict[str, Any] = {
            "meaningful": res.meaningful,
            "energy_Eh": res.energy_Eh,
            "ccsd_energy_Eh": res.ccsd_energy_Eh,
            "ccsd_t_energy_Eh": res.ccsd_t_energy_Eh,
            "t1_diagnostic": res.t1_diagnostic,
            "triples_correction": res.triples_correction,
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
            error_kind=(ErrorKind.NONE if res.ok
                        else ErrorKind.CONVERGENCE if res.scf_converged is False
                        else ErrorKind.TOOL_FAILED),
        )


register(OrcaDlpnoOpiAdapter())
