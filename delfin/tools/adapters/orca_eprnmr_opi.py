"""DELFIN's ORCA EPR/NMR building block on OPI.

This module builds an ``%eprnmr`` block — EPR g-tensor, hyperfine A-tensors,
NMR shieldings and spin-spin J couplings on chosen nuclei — entirely with
OPI's typed input classes (orca-pi 2.0.0): ``gtensor`` / ``nmrshielding`` /
``nuclei`` come from ``opi.input.blocks.block_eprnmr`` (``BlockEprnmr``,
``Nuclei``, ``NucleiFlag``), the ``!`` simple keywords from OPI's typed
enums, and the ``.inp`` file is written by OPI's ``Calculator.write_input()``.
The binary is still launched by DELFIN's existing runner (``delfin.orca.
run_orca``) and every keyword is checked against DELFIN's ORCA knowledge
store (``delfin.orca_knowledge``) before the input is written.

Manual ground truth (ORCA 6.1 manual, delfin-knowledge/orca/manual-6.1-text/
contents/spectroscopyproperties/) and the defaults this block chooses:

* ``gtensor`` / ``nmrshielding`` : the ``%eprnmr`` top-level flags that
  request the EPR g-tensor and the NMR shieldings — epr.txt section 5.22.1.2
  ("Basic usage": ``gtensor`` in the ``%eprnmr`` block with the ``SOMF`` SOC
  operator) and nmr.txt section 5.21.11 (keyword list: ``gtensor true``,
  ``nmrshielding 1``).  The g-tensor needs an open-shell reference; the
  default SOC treatment with the SOMF operator is left to ORCA rather than
  silently changed.
* ``nuclei`` : selects which nuclei get which property — epr.txt 5.22.2.2
  ("Counting of atom numbers starts from 1") for the index convention, with
  the typed property flags:
  ``aiso``/``adip`` (hyperfine), ``shift`` (NMR shielding tensor), and
  ``ssall`` (all spin-spin coupling terms, in Hz) — nmr.txt 5.21.11 and
  5.21.2.
* Basis advice (advisory only, never a silent switch): for hyperfine,
  epr.txt 5.22.2.3 "A note on basis sets" recommends core-property bases
  (EPR-II/EPR-III, IGLO-II/III, CP, uncontracted Partridge, the pcH bases by
  Jensen); for shielding, nmr.txt recommends pcSseg-n or TZVPP/QZVP.  The
  caller's explicit basis always wins.

Refusals (the reason is returned, nothing is written):

* EPR properties (g-tensor or any hyperfine request) on a closed-shell
  singlet (``mult == 1``): there is no unpaired spin for the EPR transition,
  so the numbers would be meaningless.  The caller's method, charge and
  multiplicity are never "fixed".
* an explicit nucleus index outside the molecule (a 1-based index greater
  than the number of atoms, per epr.txt 5.22.2.2).

Results come back as typed values with units via
:class:`EprNmrResult`: the g-tensor principal values and ``g_iso``
(dimensionless), one :class:`A tensor <ATensor>` per nucleus (MHz), the
isotropic shieldings (ppm) and the J couplings (Hz).  Shieldings and J
couplings are read through DELFIN's existing NMR reader
(``delfin.nmr_spectrum.parse_nmr_orca``), which is a hand-written ORCA
``.out``-text parser; the g-tensor and A-tensors are read by this module's
own hand-written parsers from the ``.out`` text (never real ORCA output is
committed).
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


@dataclass(frozen=True)
class NucleusSelection:
    """One ``Nuclei`` selection: which nuclei of one element get properties.

    ``indices`` is ``None`` for "all" of ``element`` (the ``Nuclei = all X``
    form), otherwise the explicit 1-based atom indices (manual epr.txt
    5.22.2.2: "Counting of atom numbers starts from 1").  The boolean flags
    map onto OPI's typed ``NucleiFlag`` fields.  Only hydrogen is usually
    selected with ``isotope`` (the ``ist`` override).
    """

    element: str
    indices: Optional[tuple[int, ...]] = None
    aiso: bool = False
    adip: bool = False
    shift: bool = False
    j_coupling: bool = False
    isotope: Optional[int] = None


@dataclass(frozen=True)
class Gtensor:
    """The EPR g-tensor: not a matrix but the principal values ORCA prints.

    ``principal`` are the diagonalised values ORCA prints on the ``g(tot)``
    line of the ``ELECTRONIC G-MATRIX`` block (epr.txt 5.22.1.2);
    dimensionless.  ``iso`` is the isotropic ``g_iso``.
    """

    iso: float
    principal: tuple[float, float, float]


@dataclass(frozen=True)
class ATensor:
    """One hyperfine A-tensor, as ORCA prints it (MHz, SAI convention).

    ``atom_index`` is the ORCA-printed nucleus number (1-based, epr.txt
    5.22.2.2).  ``iso_MHz`` is ``A(iso)``; ``principal_MHz`` the three
    principal values from the ``A(Tot)`` row.
    """

    atom_index: int
    element: str
    iso_MHz: float
    principal_MHz: tuple[float, float, float]


@dataclass(frozen=True)
class IsotropicShielding:
    """One nucleus' isotropic shielding (ppm), from the shielding summary.

    ``atom_index`` is the ORCA-printed nucleus index in the shielding summary
    (0-based there; ORCA writes both 0-based in the summary and 1-based in
    the tensor blocks).  ``iso_ppm`` is the absolute isotropic shielding.
    """

    atom_index: int
    element: str
    iso_ppm: float


@dataclass(frozen=True)
class JCoupling:
    """One spin-spin coupling constant (Hz), between two nuclei.

    ``nuclei`` is the unordered pair of the ORCA indices (0-based, as in the
    coupling summary); ``j_hz`` is the isotropic J in Hz (nmr.txt 5.21.2).
    """

    nuclei: tuple[int, int]
    j_hz: float


@dataclass
class EprNmrResult:
    """What a finished EPR/NMR run means, with units on every value.

    ``meaningful`` is False unless the run terminated normally and the SCF
    converged; only then are the property values populated (an unconverged
    run still prints an SCF energy).  ``findings`` carry the knowledge-card
    warnings and the basis advice; ``failure_lines`` the ``.out`` lines that
    explain a failure.
    """

    ok: bool = False
    work_dir: Optional[Path] = None
    output_file: Optional[Path] = None
    terminated_normally: bool = False
    scf_converged: bool = False
    meaningful: bool = False
    scf_energy_Eh: Optional[float] = None

    g_tensor: Optional[Gtensor] = None
    a_tensors: tuple[ATensor, ...] = ()
    shieldings: tuple[IsotropicShielding, ...] = ()
    j_couplings: tuple[JCoupling, ...] = ()

    findings: tuple[dict, ...] = ()
    failure_lines: tuple[str, ...] = ()
    error: Optional[str] = None


class EprNmrInputError(ValueError):
    """The input was refused; nothing was written."""


#: top-level ``%eprnmr`` block options this block accepts and validates
_EPRNMR_OPTIONS = ("gtensor", "nmrshielding")

#: ORCA's printed markers, from the manual
_G_TOT_RE = re.compile(
    r"^\s*g\(tot\)\s+"
    r"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+iso=\s*([-\d.eE+]+)",
    re.MULTILINE,
)
_NUCLEUS_HEADER_RE = re.compile(
    r"^\s*Nucleus\s+(\d+)\s*([A-Za-z]+)\s*:", re.MULTILINE
)
_A_TOT_RE = re.compile(
    r"^\s*A\(Tot\)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+A\(iso\)=([-\d.eE+]+)",
    re.MULTILINE,
)


def _validated_keyword(error_list: list[str], word: str, field: str) -> None:
    """Append a knowledge error when ORCA does not know ``word``."""
    from delfin.orca_knowledge import known_keyword

    if not word:
        return
    if not known_keyword(word):
        error_list.append(
            f"{field}: keyword {word!r} is not one ORCA knows "
            "(OPI catalog / ORCA parser); refused."
        )


def _nucleus_flag_names(sel: NucleusSelection) -> tuple[str, ...]:
    """The NucleiFlag names the selection turns on."""
    names = []
    if sel.aiso:
        names.append("aiso")
    if sel.adip:
        names.append("adip")
    if sel.shift:
        names.append("shift")
    if sel.j_coupling:
        names.append("ssall")
    return tuple(names)


def validate_eprnmr(
    *,
    gtensor: bool = False,
    shielding: bool = False,
    j_coupling: bool = False,
    mult: int = 1,
    elements: Sequence[str],
    nuclei: Optional[NucleusSelection] = None,
    functional: str = "",
    basis: str = "",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    nmrshielding_level: int = 1,
    cards: Optional[Sequence[dict]] = None,
) -> tuple[list[str], list[dict]]:
    """Check an EPR/NMR request against DELFIN's ORCA knowledge store.

    Returns ``(errors, warnings)``; a non-empty ``errors`` must stop the step
    (nothing is written).  *warnings* carry the knowledge-card findings and
    the basis-set advice (never a silent change).
    """
    errors: list[str] = []
    _validated_keyword(errors, functional, "functional")
    _validated_keyword(errors, basis, "basis")
    _validated_keyword(errors, ri, "ri")
    _validated_keyword(errors, aux_basis, "aux_basis")
    _validated_keyword(errors, dispersion, "dispersion")
    for kw in extra_keywords or ():
        _validated_keyword(errors, kw, "keyword")

    from delfin.orca_knowledge import known_block_option

    for opt in _EPRNMR_OPTIONS:
        if not known_block_option("eprnmr", opt):
            errors.append(
                f"eprnmr: block option {opt!r} is not one ORCA knows "
                "(OPI block model); refused."
            )

    # -- physical refusals --------------------------------------------------
    epr_requested = bool(gtensor) or (nuclei is not None and (nuclei.aiso or nuclei.adip))
    if epr_requested and mult == 1:
        errors.append(
            "EPR properties (g-tensor / hyperfine A) have no meaning on a "
            f"closed-shell singlet (mult==1, found {mult}); refused.  The "
            "g-tensor and hyperfine coupling need an unpaired spin; the "
            "multiplicity is not corrected silently."
        )

    if nuclei is not None and nuclei.indices:
        n_atoms = len(elements)
        for idx in nuclei.indices:
            if idx < 1 or idx > n_atoms:
                errors.append(
                    f"nuclei: index {idx} is outside the molecule "
                    f"({n_atoms} atoms); ORCA numbers atoms from 1 "
                    "(epr.txt 5.22.2.2); refused."
                )
                continue
            # A 1-based index list, in range, must still be refused: ORCA's
            # EPRNMR grammar (epr.txt 5.22.2.4, Table 5.14) takes the element
            # only after "all" -- ``nuclei = i1,i2,...`` carries no element --
            # but OPI's typed ``Nuclei`` always emits ``<element>`` onto the
            # line (e.g. ``nuclei = 2,3 H { aiso }``), which ORCA's parser
            # rejects ("Unrecognized symbol ... Last token: H").  DELFIN will
            # not hand-write around OPI's typed model, so an explicit index
            # list is refused and ``all <element>`` is the supported path.
            errors.append(
                f"nuclei: an explicit index list ({idx!r}) is refused — OPI's "
                "typed Nuclei writes the element on every selection, which "
                "ORCA only accepts after 'all' (Table 5.14, epr.txt 5.22.2.4); "
                "select the element with 'all <element>' instead."
            )

    # -- basis advice (advisory, never a silent switch) ---------------------
    warnings: list[dict] = []
    _basis_advice(warnings, epr_requested=epr_requested,
                  shielding=bool(shielding) or (nuclei is not None and nuclei.shift),
                  basis=basis)

    calc = {
        "keywords": list(filter(None, (functional, basis, ri, aux_basis, dispersion))),
        "blocks": {
            "eprnmr": [o for o in _EPRNMR_OPTIONS
                       if (o == "gtensor" and gtensor)
                       or (o == "nmrshielding" and shielding)],
            "nuclei": list(_nucleus_flag_names(nuclei)) if nuclei else [],
        },
        "method_family": "dft" if functional else "",
        "job": "eprnmr",
    }
    try:
        from delfin.orca_knowledge import findings as knowledge_findings

        warnings += list(knowledge_findings(calc, cards=cards) or ())
    except Exception:  # knowledge cards are advisory
        pass
    return errors, warnings


def _basis_advice(
    warnings: list[dict], *, epr_requested: bool, shielding: bool, basis: str
) -> None:
    """Append basis-set advice as warnings (never a silent basis switch).

    Cited from the manual so the advice is checked, not guessed:

    * hyperfine: epr.txt 5.22.2.3 "A note on basis sets" recommends the
      core-property bases EPR-II/EPR-III, IGLO-II/III, CP, uncontracted
      Partridge, or the pcH sets by Jensen.
    * shielding: nmr.txt requires "basis sets including tight exponents"; it
      recommends "pcSseg-n", with "TZVPP or QZVP ... sufficient in most
      cases".
    """
    if epr_requested:
        warnings.append(
            {
                "topic": "basis_advice",
                "severity": "info",
                "statement": (
                    "Hyperfine A-tensors: epr.txt 5.22.2.3 recommends a "
                    "core-property basis (EPR-II/III, IGLO-II/III, CP, "
                    "uncontracted Partridge, pcH) for a hyperfine "
                    f"calculation; the supplied basis {basis or '(default)'!r} "
                    "is used unchanged (advisory, not a silent switch)."
                ),
            }
        )
    if shielding:
        warnings.append(
            {
                "topic": "basis_advice",
                "severity": "info",
                "statement": (
                    "NMR shieldings: nmr.txt recommends pcSseg-n or, when "
                    "absent, TZVPP/QZVP ('basis sets including tight "
                    "exponents'); the supplied basis "
                    f"{basis or '(default)'!r} is used unchanged (advisory, "
                    "not a silent switch)."
                ),
            }
        )


def build_eprnmr_input(
    work_dir: Path,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    functional: str = "PBE0",
    basis: str = "def2-SVP",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    gtensor: bool = False,
    shielding: bool = False,
    j_coupling: bool = False,
    nmrshielding_level: int = 1,
    nuclei: Optional[NucleusSelection] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "eprnmr",
) -> Path:
    """Assemble the ORCA input with OPI's classes and write it with OPI.

    ``gtensor``/``shielding``/``j_coupling`` (the top-level ``%eprnmr``
    flags) and the per-nucleus ``NucleiFlag`` flags are all typed: the
    ``%eprnmr`` block is a ``BlockEprnmr`` and the nucleus selection a
    ``Nuclei`` with a ``NucleiFlag``.  Raises :class:`EprNmrInputError` on a
    refused input (after validation finds nothing — nothing is written).
    """
    errors, warnings = validate_eprnmr(
        gtensor=gtensor,
        shielding=shielding,
        j_coupling=j_coupling,
        mult=mult,
        elements=elements,
        nuclei=nuclei,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        nmrshielding_level=nmrshielding_level,
    )
    if errors:
        raise EprNmrInputError("; ".join(errors))

    from opi.core import Calculator
    from opi.input import Input
    from opi.input.blocks.block_eprnmr import BlockEprnmr, Nuclei, NucleiFlag
    from opi.input.simple_keywords.approximation import Approximation
    from opi.input.simple_keywords.base import SimpleKeyword
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dft import Dft
    from opi.input.simple_keywords.dispersion_correction import (  # noqa: E501
        DispersionCorrection,
    )
    from opi.input.structures import Atom, Structure

    def _lookup(enum_cls: type, keyword: str) -> Any:
        """OPI SimpleKeyword enum member for ``keyword``, or the bare
        validated string when OPI has none (e.g. some basis sets)."""
        target = str(keyword).lower()
        for attr in vars(enum_cls).values():
            if isinstance(attr, SimpleKeyword) and str(attr).lower() == target:
                return attr
        return str(keyword)

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

    # -- the %eprnmr block, fully typed --------------------------------
    block_kwargs: dict[str, Any] = {}
    if gtensor:
        block_kwargs["gtensor"] = 1
    if shielding:
        block_kwargs["nmrshielding"] = max(1, int(nmrshielding_level))
    if nuclei is not None:
        flag_kwargs: dict[str, Any] = {}
        flag_kwargs["aiso"] = nuclei.aiso
        flag_kwargs["adip"] = nuclei.adip
        flag_kwargs["shift"] = nuclei.shift
        flag_kwargs["ssall"] = nuclei.j_coupling
        if nuclei.isotope is not None:
            flag_kwargs["ist"] = nuclei.isotope
        # validate_eprnmr refuses an explicit index list (OPI's Nuclei writes
        # the element on every selection, which ORCA only accepts after
        # "all"); the only supported mode is therefore "all <element>".
        block_kwargs["nuclei"] = Nuclei(
            mode="all",
            atom=nuclei.element,
            flags=NucleiFlag(**flag_kwargs),
        )
    inp.add_blocks(BlockEprnmr(**block_kwargs))

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
    """The ``.out`` lines that explain a failure (same marks as the other
    OPI blocks)."""
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


def _read_g_tensor(text: str) -> Optional[Gtensor]:
    """The g-tensor principal values from the ``ELECTRONIC G-MATRIX`` block.

    Parse the ``g(tot)`` row, which ORCA prints as the three principal values
    and ``iso=`` (epr.txt 5.22.1.2).
    """
    m = _G_TOT_RE.search(text)
    if not m:
        return None
    principal = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    return Gtensor(iso=float(m.group(4)), principal=principal)


def _read_a_tensors(text: str) -> tuple[ATensor, ...]:
    """Every hyperfine A-tensor ORCA printed, in the order printed.

    ORCA prints one ``Raw HFC matrix (all values in MHz)`` block per
    selected nucleus, preceded by its ``Nucleus  nX :`` header; the
    principal values and ``A(iso)`` come from the ``A(Tot)`` row (epr.txt
    5.22.2.2).
    """
    tensors: list[ATensor] = []
    for header in _NUCLEUS_HEADER_RE.finditer(text):
        idx = int(header.group(1))
        element = header.group(2)
        block_head = header.end()
        block = text[block_head:]
        m = _A_TOT_RE.search(block)
        if not m:
            continue
        tensors.append(
            ATensor(
                atom_index=idx,
                element=element,
                iso_MHz=float(m.group(4)),
                principal_MHz=(float(m.group(1)), float(m.group(2)),
                               float(m.group(3))),
            )
        )
    return tuple(tensors)


def _read_shieldings_and_j(out_path: Path) -> tuple[tuple[IsotropicShielding, ...],
                                                    tuple[JCoupling, ...]]:
    """Isotropic shieldings (ppm) and J couplings (Hz) from the ``.out`` text.

    Delegated to DELFIN's existing NMR reader (``delfin.nmr_spectrum.
    parse_nmr_orca``), which parses the ``CHEMICAL SHIELDING SUMMARY`` and
    the ``SUMMARY OF ISOTROPIC COUPLING CONSTANTS J (Hz)`` tables — the
    same ORCA text this block launches (nmr.txt 5.21.11, 5.21.2).
    """
    from delfin.nmr_spectrum import parse_nmr_orca

    if not out_path.is_file():
        return (), ()
    try:
        parsed = parse_nmr_orca(out_path)
    except Exception as exc:
        logger.warning("could not parse NMR text from %s: %s", out_path, exc)
        return (), ()
    shieldings = tuple(
        IsotropicShielding(atom_index=s.atom_index, element=s.element,
                           iso_ppm=s.isotropic_ppm)
        for s in parsed.shieldings
    )
    couplings = tuple(
        JCoupling(nuclei=(min(i, j), max(i, j)), j_hz=value)
        for (i, j), value in parsed.j_couplings.items()
    )
    return shieldings, couplings


def read_eprnmr_output(
    work_dir: Path,
    base_name: str,
    *,
    elements: Optional[Sequence[str]] = None,
) -> EprNmrResult:
    """Read a finished EPR/NMR run from its ``.out`` text.

    The shieldings and J couplings come through DELFIN's existing NMR reader;
    the g-tensor and A-tensors through this module's hand-written parsers of
    the same text.  No value is reported unless the run terminated normally
    and the SCF converged (an unconverged run still prints an SCF energy).
    """
    from delfin.orca_read import _scf_converged_text, _terminated_normally_text

    out_path = work_dir / f"{base_name}.out"
    text = ""
    if out_path.is_file():
        try:
            text = out_path.read_text(errors="replace")
        except OSError:
            text = ""
    terminated = _terminated_normally_text(out_path)
    scf_ok = _scf_converged_text(text)
    meaningful = bool(terminated and scf_ok)

    result = EprNmrResult(
        ok=meaningful,
        work_dir=Path(work_dir),
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

    # SCF energy (the reference), if ORCA printed it
    m = re.search(r"FINAL SINGLE POINT ENERGY\s+([-\d.eE+]+)", text)
    if m:
        result.scf_energy_Eh = float(m.group(1))

    result.g_tensor = _read_g_tensor(text)
    result.a_tensors = _read_a_tensors(text)
    result.shieldings, result.j_couplings = _read_shieldings_and_j(out_path)
    return result


def run_eprnmr(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int = 0,
    mult: int = 1,
    functional: str = "PBE0",
    basis: str = "def2-SVP",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    gtensor: bool = False,
    shielding: bool = False,
    j_coupling: bool = False,
    nmrshielding_level: int = 1,
    nuclei: Optional[NucleusSelection] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "eprnmr",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> EprNmrResult:
    """Run one EPR/NMR calculation on OPI and read the result.

    Steps, in order: validate against the knowledge store (errors abort
    before any file is written), build and write the input with OPI, launch
    through DELFIN's ``delfin.orca.run_orca``, read back with the module's
    reader.  ``work_dir`` must already exist; the result's ``findings`` carry
    the basis advice.
    """
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else Path.cwd() / "eprnmr_run"
    wdir.mkdir(parents=True, exist_ok=True)

    errors, findings = validate_eprnmr(
        gtensor=gtensor,
        shielding=shielding,
        j_coupling=j_coupling,
        mult=mult,
        elements=elements,
        nuclei=nuclei,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        nmrshielding_level=nmrshielding_level,
    )
    if errors:
        return EprNmrResult(
            ok=False,
            work_dir=wdir,
            meaningful=False,
            findings=tuple(findings),
            error="; ".join(errors),
        )

    inp_path = build_eprnmr_input(
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
        gtensor=gtensor,
        shielding=shielding,
        j_coupling=j_coupling,
        nmrshielding_level=nmrshielding_level,
        nuclei=nuclei,
        extra_keywords=extra_keywords,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    out_path = wdir / f"{base_name}.out"
    success = run_orca(str(inp_path), str(out_path), timeout=timeout,
                       working_dir=wdir)
    result = read_eprnmr_output(wdir, base_name, elements=elements)
    result.findings = tuple(findings)
    if not success and result.ok:
        result.ok = False
        result.meaningful = False
        result.error = result.error or "ORCA subprocess reported failure"
        result.failure_lines = result.failure_lines or _failure_lines(out_path)
    return result


# --------------------------------------------------------------------------
# Adapter for delfin.tools (registered as ``orca_eprnmr_opi``)
# --------------------------------------------------------------------------


class OrcaEprNmrOpiAdapter(StepAdapter):
    """An EPR/NMR ``%eprnmr`` step driven by OPI, launchable via ``run_step``.

    Builds the g-tensor / hyperfine / shielding / J-coupling request with
    OPI's typed classes and reads the result with units.  Refuses EPR
    properties on a closed-shell singlet and a nucleus index outside the
    molecule; basis advice is advisory only.
    """

    name = "orca_eprnmr_opi"
    description = "OPI-backed ORCA EPR/NMR properties (g-tensor, hyperfine, shielding, J)"
    category = "properties"
    produces_geometry = False
    consumes = ("geometry",)
    produces = ("qc_output", "property_json")
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
        nuclei = kwargs.get("nuclei")
        if isinstance(nuclei, dict):
            nuclei = NucleusSelection(**nuclei)
        res = run_eprnmr(
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            functional=kwargs.get("functional", "PBE0"),
            basis=kwargs.get("basis", "def2-SVP"),
            ri=kwargs.get("ri", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            dispersion=kwargs.get("dispersion", ""),
            gtensor=kwargs.get("gtensor", False),
            shielding=kwargs.get("shielding", False),
            j_coupling=kwargs.get("j_coupling", False),
            nmrshielding_level=kwargs.get("nmrshielding_level", 1),
            nuclei=nuclei,
            extra_keywords=kwargs.get("extra_keywords"),
            maxcore=kwargs.get("maxcore", 1000),
            ncores=cores,
            base_name=kwargs.get("base_name", "eprnmr"),
            work_dir=work_dir,
            timeout=kwargs.get("timeout"),
        )
        data: dict[str, Any] = {
            "meaningful": res.meaningful,
            "scf_energy_Eh": res.scf_energy_Eh,
            "g_tensor": {
                "iso": res.g_tensor.iso,
                "principal": list(res.g_tensor.principal),
            } if res.g_tensor else None,
            "a_tensors": [
                {
                    "atom_index": a.atom_index,
                    "element": a.element,
                    "iso_MHz": a.iso_MHz,
                    "principal_MHz": list(a.principal_MHz),
                }
                for a in res.a_tensors
            ],
            "shieldings": [
                {
                    "atom_index": s.atom_index,
                    "element": s.element,
                    "iso_ppm": s.iso_ppm,
                }
                for s in res.shieldings
            ],
            "j_couplings": [
                {"nuclei": list(c.nuclei), "j_hz": c.j_hz}
                for c in res.j_couplings
            ],
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


register(OrcaEprNmrOpiAdapter())
