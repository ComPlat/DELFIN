"""DELFIN's ORCA building block for relativistic and spin properties on OPI.

Three related single-job blocks live here, all built entirely with the ORCA
Python Interface (OPI, orca-pi 2.0.0) and checked against DELFIN's ORCA
knowledge store before the input is written:

* ``rel_sp``   — scalar-relativistic single point (ZORA / DKH1 / DKH2) with
  the matching recontracted basis; refuses a non-relativistic basis with a
  relativistic Hamiltonian unless :paramref:`allow_nonrel_basis` is set.
* ``zfs``      — zero-field splitting via the ``%eprnmr`` ``DTensor`` block
  (axial ``D`` and rhombic ``E``, cm^-1, plus their ratio ``E/D``).
* ``soc_tddft``— SOC-corrected excited states via ``TDDFT`` with the
  ``%tddft`` ``soc`` flag; reads the SOC-corrected transition list.

The geometry, the ``!`` simple keywords and the ``%rel`` / ``%eprnmr`` /
``%tddft`` blocks come from OPI's typed objects (``opi.input``); ``%pal`` /
``%maxcore`` map onto ``Input.ncores``/``Input.memory``; the ``.inp`` is
written by OPI's ``Calculator.write_input()`` — no hand-written string
template.

Manual sections each default is cited against (ORCA 6.1 text manual:

* ``relativity.txt`` 2.12 "Relativistic effects" — DKH/ZORA variants,
  ``%rel`` block, relativistic (recontracted) basis requirement.
* ``epr.txt`` — the ``%eprnmr`` block and the ZFS ``D``/``E`` values (cm^-1).
* ``tddft.txt`` 2.5 "TD-DFT" — ``%tddft`` ``soc`` option for SOC-TDDFT.

Reading goes through ``delfin.orca_read``: ZFS D/E from ORCA's ``D_Tensor``
property record, SOC states from the ``.out`` SOC-corrected spectrum section.
"""
from __future__ import annotations

from dataclasses import dataclass
import time
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import ErrorKind, StepResult, StepStatus

logger = get_logger(__name__)

# Relativistic Hamiltonian values (case-blind) that require a matching
# recontracted / decontracted basis.  ORCA manual relativity.txt 2.12.6:
# "the basis set must be of the appropriate (recontracted) type".
_JOBS = ("rel_sp", "zfs", "soc_tddft")

_RELATIVISTIC_HAMILTONIANS = (
    "ZORA", "DKH", "DKH1", "DKH2", "X2C", "IORA", "IORAMM", "RESC", "REX2",
)
# Basis markers that satisfy the relativistic-basis requirement.
_RELATIVISTIC_BASIS_MARKERS = (
    "ZORA", "DKH", "SARC", "RELES", "DYALL", "UNCONTRACTED", "X2C",
)
# Products of a complex Dirac-Kohn-Sham step; the defect (dE_rad in atomic
# units, ppm) is read back as the primary number.
_DK_PRODUCTS = ("DK",)


@dataclass
class OrcaRelSpinResult:
    """Typed result of one OPI-backed relativistic/spin ORCA step.

    Fields sourced from OPI are populated only after the run gate and the
    per-job convergence check make them meaningful; :attr:`meaningful`
    summarises that state and every numeric field is ``None`` otherwise.
    """

    ok: bool
    work_dir: Path
    output_file: Optional[Path] = None
    job: str = "rel_sp"

    # Convergence gate (always read; gates every numeric field).
    terminated_normally: Optional[bool] = None
    # True exactly when terminated_normally and the job-specific check
    # (e.g. SCF convergence) both hold.
    meaningful: bool = False

    # Scalar-relativistic single point.
    energy_Eh: Optional[float] = None

    # Zero-field splitting (cm^-1 and dimensionless ratio).
    zfs_D_cm1: Optional[float] = None
    zfs_E_cm1: Optional[float] = None
    zfs_E_over_D: Optional[float] = None

    # SOC-corrected excited states (TDDFT).
    soc_states: Optional[tuple] = None  # tuple[delfin.orca_read.SocState, ...]

    findings: tuple = ()
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


def _validate_block_values(block: str, options: Mapping[str, Any]) -> list[str]:
    """Check the *values* of block options against OPI's typed block model.

    ``known_block_option`` only knows the option *names*; the value check
    runs OPI's own pydantic validation: the block class is constructed with
    the options and a ``ValidationError`` comes back as a refusal reason.
    When OPI is not importable (plain CI interpreter) the check is skipped
    and the caller relies on the name check above.
    """
    errors: list[str] = []
    try:
        import opi  # noqa: F401
    except ImportError:
        return errors
    if not options:
        return errors
    module = {
        "rel": "opi.input.blocks.block_rel",
        "eprnmr": "opi.input.blocks.block_eprnmr",
        "tddft": "opi.input.blocks.block_tddft",
    }.get(block)
    cls_name = {"rel": "BlockRel", "eprnmr": "BlockEprnmr",
                "tddft": "BlockTddft"}.get(block)
    if module is None or cls_name is None:
        return errors
    import importlib

    try:
        block_cls = getattr(importlib.import_module(module), cls_name)
    except ImportError:
        return errors
    try:
        block_cls(**dict(options))
    except Exception as exc:  # pydantic ValidationError and friends
        errors.append(
            f"{block}: block option value rejected by OPI's typed model "
            f"({exc}); refusing."
        )
    return errors


def validate_method(
    *,
    job: str,
    functional: str = "",
    hamiltonian: str = "",
    basis: str,
    extra_keywords: Optional[Sequence[str]] = None,
    rel_options: Optional[Mapping[str, Any]] = None,
    epr_options: Optional[Mapping[str, Any]] = None,
    tddft_options: Optional[Mapping[str, Any]] = None,
) -> tuple[list[str], list[dict]]:
    """Check every keyword and block option against DELFIN's ORCA knowledge
    store.  Returns ``(errors, warnings)``; a non-empty ``errors`` stops the
    step.  *warnings* are findings of the applicable verified knowledge cards,
    reported and never silently dropped.

    :paramref:`basis` is checked against the relativistic-basis requirement
    whenever a relativistic Hamiltonian is requested and
    ``allow_nonrel_basis`` is not set (handled by the caller before any
    launch, so a mismatched basis is refused with a reason here).
    """
    errors: list[str] = []
    if not functional:
        errors.append(
            "functional is required: ORCA otherwise runs the pure SCF "
            "Hartree-Fock silently ('Ab initio Hamiltonian Method .... "
            "Hartree-Fock'), which changes the chemistry; refusing."
        )
    else:
        _validated_keyword(errors, functional, "functional")
    if hamiltonian:
        _validated_keyword(errors, hamiltonian, "hamiltonian")
    _validated_keyword(errors, basis, "basis")
    for kw in extra_keywords or ():
        _validated_keyword(errors, kw, "keyword")

    # Block options must be in the ORCA knowledge store's block model.
    from delfin.orca_knowledge import known_block_option

    for block, opts in (
        ("rel", rel_options),
        ("eprnmr", epr_options),
        ("tddft", tddft_options),
    ):
        for option in (opts or {}):
            if not known_block_option(block, option):
                errors.append(
                    f"{block}: block option {option!r} is not one ORCA knows "
                    "(OPI block model); refused."
                )
        errors.extend(_validate_block_values(block, opts or {}))

    calc = {
        "keywords": [k for k in (functional, hamiltonian, basis) if k]
        + [str(k) for k in (extra_keywords or ())],
        "blocks": {
            "rel": dict(rel_options or {}),
            "eprnmr": dict(epr_options or {}),
            "tddft": dict(tddft_options or {}),
        },
        "method_family": "dft" if job != "soc_tddft" else "tddft",
        "job": job,
    }
    from delfin.orca_knowledge import findings as knowledge_findings

    return errors, knowledge_findings(calc)


def validate_relspin(
    *,
    job: str,
    functional: str = "",
    basis: str,
    hamiltonian: str = "",
    elements: Optional[Sequence[str]] = None,
    charge: int = 0,
    mult: int = 1,
    nroots: int = 8,
    epr_options: Optional[Mapping[str, Any]] = None,
    rel_options: Optional[Mapping[str, Any]] = None,
    tddft_options: Optional[Mapping[str, Any]] = None,
    allow_nonrel_basis: bool = False,
) -> tuple[list[str], list[dict]]:
    """Validate one ``relspin`` job end to end, refusing what is physically
    wrong before any input is written.

    Checks, in order: the job name is one of the three modes; every keyword
    and block option is in DELFIN's ORCA knowledge store
    (:func:`validate_method`); the multiplicity is consistent with the
    electron count (:func:`_check_mult_triplet`, when a geometry is given);
    ZFS requires an open shell; SOC-TDDFT requires ``nroots >= 1``; and a
    relativistic Hamiltonian requires a matching recontracted basis unless
    ``allow_nonrel_basis`` lifts it.  Returns ``(errors, warnings)`` with
    the refusal reason inside each error string.
    """
    errors: list[str] = []
    warnings: list[dict] = []

    if job not in _JOBS:
        errors.append(
            f"unknown job {job!r}: relspin offers "
            f"{', '.join(sorted(_JOBS))}; refusing."
        )
        return errors, warnings

    kw_errors, kw_warnings = validate_method(
        job=job,
        functional=functional,
        hamiltonian=hamiltonian,
        basis=basis,
        rel_options=rel_options,
        epr_options=epr_options,
        tddft_options=tddft_options,
    )
    errors.extend(kw_errors)
    warnings.extend(kw_warnings)

    spin_error = _check_mult_triplet(
        charge=charge, mult=mult, elements=elements)
    if spin_error:
        errors.append(spin_error)

    if job == "zfs" and mult == 1:
        errors.append(
            "zero-field splitting needs an open shell: multiplicity 1 is a "
            "singlet, which has no D/E parameters (manual epr.txt); "
            "refusing rather than changing the spin state."
        )
    if job == "soc_tddft" and nroots < 1:
        errors.append(
            "SOC-TDDFT needs nroots >= 1; got 0. Refusing an empty "
            "excited-state list rather than silently running a ground state."
        )
    rel = _relativity_requested((hamiltonian,) if hamiltonian else ())
    if rel:
        basis_error = _basis_check_error(
            hamiltonian=rel,
            basis=basis,
            allow_nonrel_basis=allow_nonrel_basis,
        )
        if basis_error:
            errors.append(basis_error)
    return errors, warnings


def _relativity_requested(keywords: Sequence[str]) -> Optional[str]:
    """Return the relativistic Hamiltonian keyword among ``keywords``
    (case-blind) when one is present, else ``None``."""
    for kw in keywords:
        up = str(kw).upper()
        for h in _RELATIVISTIC_HAMILTONIANS:
            if h in up:
                return h
    return None


def _basis_is_relativistic(basis: str) -> bool:
    b = basis.upper()
    return any(m in b for m in _RELATIVISTIC_BASIS_MARKERS)


def _n_electrons(elements: Sequence[str], charge: int) -> int:
    """Total electron count of the system: sum of the atomic numbers minus
    the charge.  Atomic numbers come from RDKit's periodic table (the same
    source atom_mapping.py uses)."""
    from delfin.atom_mapping import _periodic_table

    pt = _periodic_table()
    return sum(pt.GetAtomicNumber(el.strip().capitalize()) for el in elements) \
        - charge


def _check_mult_triplet(*, charge: int, mult: int,
                        elements: Optional[Sequence[str]] = None
                        ) -> Optional[str]:
    """Refuse a physically wrong spin state.

    For a fixed electron count a non-integer total spin cannot be a pure
    spin-multiplicity: an even number of electrons must give an odd
    multiplicity.  When the geometry is known the check runs on the real
    electron count (sum of atomic numbers minus charge); without a geometry
    only the trivial ``mult >= 1`` guard runs.
    Returns a refusal reason or ``None``.
    """
    if mult < 1:
        return f"multiplicity must be >= 1, got {mult} (physically impossible)"
    if elements is None:
        return None
    n_electrons = _n_electrons(elements, charge)
    if n_electrons < 1:
        return (
            "the electron count is not positive: "
            f"{n_electrons} electrons (elements {list(elements)}, charge "
            f"{charge}); refusing rather than launching a meaningless job."
        )
    # Electron count parity: N_even -> S integer -> 2S+1 odd.
    if n_electrons % 2 == 0 and mult % 2 == 0:
        return (
            "multiplicity is inconsistent with the electron count: "
            f"{n_electrons} electrons force an odd multiplicity, got {mult}. "
            "Refusing rather than silently changing the spin state."
        )
    return None


def _basis_check_error(
    *, hamiltonian: str, basis: str, allow_nonrel_basis: bool
) -> Optional[str]:
    """Refuse a relativistic Hamiltonian with a non-relativistic basis unless
    ``allow_nonrel_basis`` explicitly lifts the requirement."""
    if not _basis_is_relativistic(basis):
        if not allow_nonrel_basis:
            return (
                f"relativistic Hamiltonian {hamiltonian!r} needs a matching "
                f"(recontracted) basis, got {basis!r}.  Set "
                "allow_nonrel_basis=True only if the atomic set is "
                "deliberately non-relativistic."
            )
    return None


def build_rel_sp_input(
    wdir: Path,
    *,
    hamiltonian: str,
    functional: str = "",
    basis: str,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    rel_options: Optional[Mapping[str, Any]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    maxcore: int,
    ncores: int,
    base_name: str,
    allow_nonrel_basis: bool = False,
) -> Path:
    """Assemble a scalar-relativistic single-point ORCA input with OPI and
    validate every keyword/block option first.

    Returns the absolute path of the written ``<base_name>.inp``.  Raises
    ``ValueError`` with a plain reason on any knowledge / physics refusal.
    """
    errors, _ = validate_method(
        job="rel_sp",
        functional=functional,
        hamiltonian=hamiltonian,
        basis=basis,
        extra_keywords=extra_keywords,
        rel_options=rel_options,
    )
    if errors:
        raise ValueError("; ".join(errors))

    re = _check_mult_triplet(
        charge=charge, mult=mult, elements=elements)
    if re:
        raise ValueError(re)
    h = _relativity_requested((hamiltonian,))
    if h:
        be = _basis_check_error(
            hamiltonian=h, basis=basis, allow_nonrel_basis=allow_nonrel_basis
        )
        if be:
            raise ValueError(be)

    from opi.input.blocks.block_rel import BlockRel
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dft import Dft
    from opi.input.simple_keywords.relativistic_correction import (
        RelativisticCorrection,
    )
    from opi.input import Input

    inp = Input()
    keywords: list[Any] = [
        _lookup(Dft, functional), _lookup(BasisSet, basis)
    ]
    h_kw = _lookup(RelativisticCorrection, hamiltonian)
    keywords.append(h_kw)
    for kw in extra_keywords or ():
        keywords.append(str(kw))
    inp.add_simple_keywords(*keywords)
    inp.ncores = ncores
    inp.memory = maxcore

    opts = dict(rel_options or {})
    # Default the scalar-relativistic %rel method to match the simple
    # keyword: ZORA -> zora, DKH/DKH1/DKH2 -> dkh (manual sells %rel method;
    # the DKH order is carried by the DKH1/DKH2 keyword itself).
    if "method" not in opts:
        method = "zora" if "ZORA" in hamiltonian.upper() else "dkh"
        opts["method"] = method
    inp.add_blocks(BlockRel(**opts))

    return _structure_and_write(
        inp, elements, coordinates, charge, mult,
        wdir=wdir, base_name=base_name,
    )


def _lookup(enum_cls: type, keyword: str) -> Any:
    """Return the OPI SimpleKeyword enum member whose ORCA keyword matches
    ``keyword`` (case-blind) or the bare validated string when OPI has no
    member for it."""
    from opi.input.simple_keywords.base import SimpleKeyword

    for member in vars(enum_cls).values():
        if isinstance(member, SimpleKeyword):
            if str(member).upper() == str(keyword).upper():
                return member
    return str(keyword)


def _write_input(inp: Any, *, wdir: Path, base_name: str,
                 structure: Any = None) -> Path:
    """Write ``inp`` to ``<base_name>.inp`` under ``wdir`` with OPI.

    OPI renders the coordinate block from ``Calculator.structure``, so the
    geometry is attached there (as the sp adapter does); the structure is
    optional so a block-only input can still be written.
    """
    from opi.core import Calculator

    wdir.mkdir(parents=True, exist_ok=True)
    calc = Calculator(base_name, working_dir=wdir, version_check=False)
    calc.input = inp
    if structure is not None:
        calc.structure = structure
    calc.write_input(force=True)
    path = wdir / f"{base_name}.inp"
    if not path.exists():
        raise RuntimeError(f"OPI did not produce {path}")
    return path


def _read_output_with_reason(wdir: Path, base_name: str) -> tuple:
    """``(out, reason)``: OPI's parsed output handle for the finished run, or
    ``None`` plus the reason OPI could not parse it.

    OPI's ``Output`` takes the run's base name and working directory and
    reads ``<base_name>.out`` (plus the property JSON when present) from
    there."""
    from opi.output.core import Output

    try:
        out = Output(basename=base_name, working_dir=wdir)
        out.parse()
    except Exception as exc:  # noqa: BLE001 - reported, never swallowed
        first = str(exc).strip().splitlines()[0] if str(exc).strip() else ""
        return None, f"OPI could not parse the output ({type(exc).__name__}: {first})"
    return out, None


def _read_output(wdir: Path, base_name: str) -> Any:
    """OPI's parsed output handle, or ``None`` when OPI cannot parse it."""
    return _read_output_with_reason(wdir, base_name)[0]


def _terminated_normally(wdir: Path, base_name: str) -> bool:
    """DELFIN's own text check for ORCA's normal-termination marker.

    The ZFS and SOC-TDDFT readers do not need OPI's parsed model (D/E come
    from the property JSON, the SOC states from the ``.out`` table), and
    OPI 2.0.0's property model rejects the property JSON of a real ORCA
    6.1.1 ``%tddft soc true`` run (SOC states carry multiplicity 0), so
    gating on OPI's parse would fail a run that ended normally."""
    from delfin.orca_read import _terminated_normally_text

    return _terminated_normally_text(wdir / f"{base_name}.out")


def _scf_converged(wdir: Path, base_name: str) -> bool:
    """DELFIN's own SCF-convergence text check (delfin.orca_read)."""
    from delfin.orca_read import _scf_converged_text

    try:
        text = (wdir / f"{base_name}.out").read_text(errors="replace")
    except OSError:
        return False
    return _scf_converged_text(text)


# ---------------------------------------------------------------------------
# run_* entry points (one per job mode).  Each validates, writes the OPI
# input, launches ORCA through DELFIN's existing launcher, and reads the
# typed result.
# ---------------------------------------------------------------------------

def read_rel_sp_result(wdir: Path, base_name: str) -> OrcaRelSpinResult:
    """Parse a finished scalar-relativistic single-point run into
    :class:`OrcaRelSpinResult` (final energy via OPI's parsed output)."""
    out, reason = _read_output_with_reason(wdir, base_name)
    ok = terminated = out is not None and getattr(
        out, "terminated_normally", lambda: False
    )()
    energy = None
    if ok and not _scf_converged(wdir, base_name):
        # An unconverged SCF still prints an energy; never report it.
        ok = False
        reason = ("the SCF did not converge (no 'SCF CONVERGED' in the "
                  "output); the energy is not reported")
    elif ok:
        energy = out.get_final_energy()
    elif reason is None:
        reason = "ORCA did not terminate normally"
    return OrcaRelSpinResult(
        ok=ok,
        work_dir=wdir,
        output_file=wdir / f"{base_name}.out",
        job="rel_sp",
        terminated_normally=terminated,
        meaningful=ok and energy is not None,
        energy_Eh=energy if ok else None,
        error=reason,
    )


def build_zfs_input(
    wdir: Path,
    *,
    functional: str = "",
    basis: str,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    epr_options: Optional[Mapping[str, Any]] = None,
    hamiltonian: str = "",
    rel_options: Optional[Mapping[str, Any]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    maxcore: int,
    ncores: int,
    base_name: str,
    allow_nonrel_basis: bool = False,
) -> Path:
    """Assemble a ZFS (``%eprnmr``) ORCA input with OPI, validated first.

    The single-triplet requirement (an unpaired-electron system) is enforced
    here; a singlet (``mult == 1``) is refused because a singlet has no
    zero-field splitting.
    """
    errors, _ = validate_method(
        job="zfs",
        functional=functional,
        hamiltonian=hamiltonian,
        basis=basis,
        extra_keywords=extra_keywords,
        rel_options=rel_options,
        epr_options=epr_options,
    )
    if errors:
        raise ValueError("; ".join(errors))
    if mult < 2:
        raise ValueError(
            f"ZFS needs an open shell (multiplicity >= 2); got mult={mult}."
            "  A singlet has no zero-field splitting; refusing rather than "
            "changing the spin state."
        )
    re = _check_mult_triplet(
        charge=charge, mult=mult, elements=elements)
    if re:
        raise ValueError(re)
    h = _relativity_requested((hamiltonian,))
    if h:
        be = _basis_check_error(
            hamiltonian=h, basis=basis, allow_nonrel_basis=allow_nonrel_basis
        )
        if be:
            raise ValueError(be)

    from opi.input.blocks.block_eprnmr import BlockEprnmr
    from opi.input.blocks.block_rel import BlockRel
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dft import Dft
    from opi.input.simple_keywords.relativistic_correction import (
        RelativisticCorrection,
    )
    from opi.input import Input

    inp = Input()
    keywords: list[Any] = [
        _lookup(Dft, functional), _lookup(BasisSet, basis)
    ]
    if hamiltonian:
        keywords.append(_lookup(RelativisticCorrection, hamiltonian))
    for kw in extra_keywords or ():
        keywords.append(str(kw))
    inp.add_simple_keywords(*keywords)
    inp.ncores = ncores
    inp.memory = maxcore

    # ZFS tensor: spin-spin and spin-orbit (the total D and E) unless the
    # caller names a narrower DSS/DSOC mix.
    epr = dict(epr_options or {})
    if "dtensor" not in epr:
        # epr.txt (DTensor keywords): DTensor SSAndSo makes ORCA compute both
        # the spin-spin and the spin-orbit contribution, so D and E are the
        # total zero-field-splitting parameters; DTensor SS would leave the
        # spin-orbit part out of the read-back D and E.
        epr["dtensor"] = "ssandso"
    blocks: list[Any] = [BlockEprnmr(**epr)]
    if hamiltonian or rel_options:
        ropts = dict(rel_options or {})
        if "method" not in ropts:
            method = "zora" if "ZORA" in hamiltonian.upper() else "dkh"
            ropts["method"] = method
        blocks.append(BlockRel(**ropts))
    inp.add_blocks(*blocks)

    return _structure_and_write(
        inp, elements, coordinates, charge, mult,
        wdir=wdir, base_name=base_name,
    )


def read_zfs_result(wdir: Path, base_name: str) -> OrcaRelSpinResult:
    """Parse a finished ZFS run into :class:`OrcaRelSpinResult` (D, E, E/D in
    cm^-1 via delfin.orca_read)."""
    from delfin.orca_read import read_zfs_parameters

    if not _terminated_normally(wdir, base_name):
        return OrcaRelSpinResult(ok=False, work_dir=wdir,
                                 output_file=wdir / f"{base_name}.out",
                                 job="zfs", terminated_normally=False,
                                 error="ORCA did not terminate normally")
    z = read_zfs_parameters(wdir, base_name)
    if z is None:
        return OrcaRelSpinResult(ok=True, work_dir=wdir,
                                 output_file=wdir / f"{base_name}.out",
                                 job="zfs", terminated_normally=True,
                                 meaningful=False, error="no D_Tensor record")
    return OrcaRelSpinResult(
        ok=True,
        work_dir=wdir,
        output_file=wdir / f"{base_name}.out",
        job="zfs",
        terminated_normally=True,
        meaningful=True,
        zfs_D_cm1=z.d_cm1,
        zfs_E_cm1=z.e_cm1,
        zfs_E_over_D=z.e_over_d,
    )


def build_soc_tddft_input(
    wdir: Path,
    *,
    functional: str = "",
    basis: str,
    nroots: int,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    tddft_options: Optional[Mapping[str, Any]] = None,
    hamiltonian: str = "",
    rel_options: Optional[Mapping[str, Any]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    maxcore: int,
    ncores: int,
    base_name: str,
    allow_nonrel_basis: bool = False,
) -> Path:
    """Assemble a SOC-corrected TDDFT ORCA input with OPI, validated first.

    Requires at least one excited root (``nroots >= 1``); with ``nroots < 1``
    the request is physically empty and refused.
    """
    errors, _ = validate_method(
        job="soc_tddft",
        functional=functional,
        hamiltonian=hamiltonian,
        basis=basis,
        extra_keywords=extra_keywords,
        rel_options=rel_options,
        tddft_options=tddft_options,
    )
    if errors:
        raise ValueError("; ".join(errors))
    if nroots < 1:
        raise ValueError(
            f"SOC-TDDFT needs nroots >= 1; got {nroots}.  Refusing an empty "
            "excited-state list rather than silently running a ground state."
        )
    re = _check_mult_triplet(
        charge=charge, mult=mult, elements=elements)
    if re:
        raise ValueError(re)
    h = _relativity_requested((hamiltonian,))
    if h:
        be = _basis_check_error(
            hamiltonian=h, basis=basis, allow_nonrel_basis=allow_nonrel_basis
        )
        if be:
            raise ValueError(be)

    from opi.input.blocks.block_tddft import BlockTddft
    from opi.input.blocks.block_rel import BlockRel
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dft import Dft
    from opi.input.simple_keywords.relativistic_correction import (
        RelativisticCorrection,
    )
    from opi.input import Input

    inp = Input()
    keywords: list[Any] = [
        _lookup(Dft, functional), _lookup(BasisSet, basis)
    ]
    if hamiltonian:
        keywords.append(_lookup(RelativisticCorrection, hamiltonian))
    for kw in extra_keywords or ():
        keywords.append(str(kw))
    inp.add_simple_keywords(*keywords)
    inp.ncores = ncores
    inp.memory = maxcore

    tddft = dict(tddft_options or {})
    tddft["nroots"] = nroots
    if "soc" not in tddft:
        tddft["soc"] = True
    blocks: list[Any] = [BlockTddft(**tddft)]
    if hamiltonian or rel_options:
        ropts = dict(rel_options or {})
        if "method" not in ropts:
            method = "zora" if "ZORA" in hamiltonian.upper() else "dkh"
            ropts["method"] = method
        blocks.append(BlockRel(**ropts))
    inp.add_blocks(*blocks)

    return _structure_and_write(
        inp, elements, coordinates, charge, mult,
        wdir=wdir, base_name=base_name,
    )


def read_soc_tddft_result(wdir: Path, base_name: str) -> OrcaRelSpinResult:
    """Parse a finished SOC-TDDFT run into :class:`OrcaRelSpinResult`
    (``soc_states`` via delfin.orca_read)."""
    from delfin.orca_read import read_soc_states

    if not _terminated_normally(wdir, base_name):
        return OrcaRelSpinResult(ok=False, work_dir=wdir,
                                 output_file=wdir / f"{base_name}.out",
                                 job="soc_tddft", terminated_normally=False,
                                 error="ORCA did not terminate normally")
    states = read_soc_states(wdir, base_name)
    if not states:
        # Fail closed: a SOC-TDDFT step without its SOC-corrected states
        # has no result, however normally ORCA ended.
        return OrcaRelSpinResult(
            ok=False, work_dir=wdir, output_file=wdir / f"{base_name}.out",
            job="soc_tddft", terminated_normally=True, meaningful=False,
            error=("no SOC-corrected states: the 'SOC CORRECTED ABSORPTION "
                   "SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' table "
                   "is missing or empty in the output"),
        )
    return OrcaRelSpinResult(
        ok=True,
        work_dir=wdir,
        output_file=wdir / f"{base_name}.out",
        job="soc_tddft",
        terminated_normally=True,
        meaningful=True,
        soc_states=tuple(states),
    )


def _structure_and_write(
    inp: Any, elements: Sequence[str], coordinates: Sequence[Sequence[float]],
    charge: int, mult: int, *, wdir: Path, base_name: str,
) -> Path:
    """Attach the geometry to ``inp`` and write it with OPI."""
    from opi.input.structures import Atom, Structure

    atoms = [
        Atom(el, coordinates=(float(x), float(y), float(z)))
        for el, (x, y, z) in zip(elements, coordinates)
    ]
    structure = Structure(atoms, charge=charge, multiplicity=mult)
    return _write_input(inp, wdir=wdir, base_name=base_name,
                        structure=structure)


def run_rel_sp(
    work_dir: Path,
    *,
    hamiltonian: str,
    basis: str,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "rel_sp",
    rel_options: Optional[Mapping[str, Any]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    allow_nonrel_basis: bool = False,
    functional: str = "",
) -> OrcaRelSpinResult:
    """Build, launch and read a scalar-relativistic single point."""
    try:
        build_rel_sp_input(
            work_dir,
            functional=functional,
            hamiltonian=hamiltonian,
            basis=basis,
            elements=elements,
            coordinates=coordinates,
            charge=charge,
            mult=mult,
            rel_options=rel_options,
            extra_keywords=extra_keywords,
            maxcore=maxcore,
            ncores=ncores,
            base_name=base_name,
            allow_nonrel_basis=allow_nonrel_basis,
        )
    except ValueError as exc:
        return OrcaRelSpinResult(ok=False, work_dir=work_dir, job="rel_sp",
                                 error=str(exc))
    launch_error = _launch(work_dir, base_name)
    if launch_error:
        return OrcaRelSpinResult(ok=False, work_dir=work_dir, job="rel_sp",
                                 output_file=work_dir / f"{base_name}.out",
                                 error=launch_error)
    return read_rel_sp_result(work_dir, base_name)


def run_zfs(
    work_dir: Path,
    *,
    basis: str,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "zfs",
    epr_options: Optional[Mapping[str, Any]] = None,
    hamiltonian: str = "",
    rel_options: Optional[Mapping[str, Any]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    allow_nonrel_basis: bool = False,
    functional: str = "",
) -> OrcaRelSpinResult:
    """Build, launch and read a zero-field-splitting run."""
    try:
        build_zfs_input(
            work_dir,
            functional=functional,
            basis=basis,
            elements=elements,
            coordinates=coordinates,
            charge=charge,
            mult=mult,
            epr_options=epr_options,
            hamiltonian=hamiltonian,
            rel_options=rel_options,
            extra_keywords=extra_keywords,
            maxcore=maxcore,
            ncores=ncores,
            base_name=base_name,
            allow_nonrel_basis=allow_nonrel_basis,
        )
    except ValueError as exc:
        return OrcaRelSpinResult(ok=False, work_dir=work_dir, job="zfs",
                                 error=str(exc))
    launch_error = _launch(work_dir, base_name)
    if launch_error:
        return OrcaRelSpinResult(ok=False, work_dir=work_dir, job="zfs",
                                 output_file=work_dir / f"{base_name}.out",
                                 error=launch_error)
    return read_zfs_result(work_dir, base_name)


def run_soc_tddft(
    work_dir: Path,
    *,
    basis: str,
    nroots: int,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "soc_tddft",
    tddft_options: Optional[Mapping[str, Any]] = None,
    hamiltonian: str = "",
    rel_options: Optional[Mapping[str, Any]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    allow_nonrel_basis: bool = False,
    functional: str = "",
) -> OrcaRelSpinResult:
    """Build, launch and read a SOC-corrected TDDFT run."""
    try:
        build_soc_tddft_input(
            work_dir,
            functional=functional,
            basis=basis,
            nroots=nroots,
            elements=elements,
            coordinates=coordinates,
            charge=charge,
            mult=mult,
            tddft_options=tddft_options,
            hamiltonian=hamiltonian,
            rel_options=rel_options,
            extra_keywords=extra_keywords,
            maxcore=maxcore,
            ncores=ncores,
            base_name=base_name,
            allow_nonrel_basis=allow_nonrel_basis,
        )
    except ValueError as exc:
        return OrcaRelSpinResult(ok=False, work_dir=work_dir, job="soc_tddft",
                                 error=str(exc))
    launch_error = _launch(work_dir, base_name)
    if launch_error:
        return OrcaRelSpinResult(ok=False, work_dir=work_dir,
                                 job="soc_tddft",
                                 output_file=work_dir / f"{base_name}.out",
                                 error=launch_error)
    return read_soc_tddft_result(work_dir, base_name)


def _launch(wdir: Path, base_name: str) -> Optional[str]:
    """Launch ORCA on the written input through DELFIN's existing runner.

    Returns an error string on failure to launch, else ``None``.  The real
    calculation is never run in CI; the caller's launcher is stubbed there."""
    from delfin.orca import run_orca

    inp = wdir / f"{base_name}.inp"
    out = wdir / f"{base_name}.out"
    try:
        run_orca(str(inp), str(out), working_dir=wdir)
    except Exception as exc:  # noqa: BLE001 - report the launch failure
        return f"failed to launch ORCA: {exc}"
    # The runner's return value is ignored on purpose: a %eprnmr input goes
    # through run_orca's automatic NMR-report path, which returns False when
    # the run has no H shieldings to plot -- that says nothing about ORCA's
    # own termination.  The readers gate on terminated_normally() from the
    # parsed .out instead.
    return None


class OrcaRelSpinOpiAdapter(StepAdapter):
    """DELFIN step wrapper for the relativistic/spin OPI block."""

    name = "orca_relspin_opi"
    description = (
        "Relativistic and spin ORCA blocks on OPI: ZORA/DKH scalar "
        "relativity, zero-field splitting (%eprnmr D/E), SOC-TDDFT excited "
        "states."
    )
    category = "relativistic"
    produces_geometry = False
    consumes = ()
    produces = ("result.json",)
    requires_python = ()

    def execute(self, work_dir: Path, *, cores: int = 1,
                **kwargs: Any) -> StepResult:
        """Run one relspin job (``job``: rel_sp | zfs | soc_tddft).

        Standard step contract: ``cores`` becomes the input's ``%pal``;
        the other kwargs go to the matching ``run_*`` function.  The numbers
        are returned in ``StepResult.data``; a FAILED result always says why.
        """
        start = time.monotonic()
        wdir = Path(work_dir)
        kwargs.pop("_prev_artifacts", None)
        kwargs.pop("geometry", None)
        job = kwargs.pop("job", "rel_sp")
        runners = {"rel_sp": run_rel_sp, "zfs": run_zfs,
                   "soc_tddft": run_soc_tddft}
        if job not in runners:
            return self._make_result(
                self.name, StepStatus.FAILED, wdir, start,
                error=f"unknown job {job!r}; use one of {sorted(runners)}",
                error_kind=ErrorKind.MISSING_PARAM,
            )
        kwargs.setdefault("base_name", job)
        try:
            result = runners[job](wdir, ncores=cores, **kwargs)
        except TypeError as exc:  # a missing or unknown parameter
            return self._make_result(
                self.name, StepStatus.FAILED, wdir, start,
                error=f"{job}: {exc}", error_kind=ErrorKind.MISSING_PARAM,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("relspin %s failed: %s", job, exc)
            return self._make_result(
                self.name, StepStatus.FAILED, wdir, start,
                error=f"{job}: {exc}", error_kind=ErrorKind.INTERNAL,
            )
        data = self._result_data(result)
        if not result.ok:
            return self._make_result(
                self.name, StepStatus.FAILED, wdir, start, data=data,
                output_file=result.output_file,
                error=result.error or f"{job} run failed",
                error_kind=ErrorKind.TOOL_FAILED,
            )
        return self._make_result(
            self.name, StepStatus.SUCCESS, wdir, start, data=data,
            output_file=result.output_file,
        )

    @staticmethod
    def _result_data(result: OrcaRelSpinResult) -> dict:
        return {
            "job": result.job,
            "terminated_normally": result.terminated_normally,
            "meaningful": result.meaningful,
            "energy_Eh": result.energy_Eh,
            "zfs_D_cm1": result.zfs_D_cm1,
            "zfs_E_cm1": result.zfs_E_cm1,
            "zfs_E_over_D": result.zfs_E_over_D,
            "soc_states": [
                {"states": s.states, "e_cm1": s.e_cm1, "e_eV": s.e_eV,
                 "wavelength_nm": s.wavelength_nm, "fosc": s.fosc,
                 "mult": s.mult}
                for s in (result.soc_states or ())
            ],
            "work_dir": str(result.work_dir),
            "output_file": str(result.output_file) if result.output_file else None,
            "error": result.error,
        }


register(OrcaRelSpinOpiAdapter())
