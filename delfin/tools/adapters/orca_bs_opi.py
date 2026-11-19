"""DELFIN's broken-symmetry (BS) single-point building block on OPI.

This module builds a broken-symmetry DFT calculation entirely with the ORCA
Python Interface (OPI, orca-pi 2.0.0) and reads its result with DELFIN's
established reader.  It follows the same shape as the other ``*_opi``
adapter modules in this package: validate every keyword / block option
against DELFIN's ORCA knowledge store (``delfin.orca_knowledge``), build the
input with OPI's typed classes, write it with OPI's
``Calculator.write_input()``, launch through ``delfin.orca.run_orca`` and
read the answer back without hand-writing strings wherever OPI has a type.

Manual ground truth (ORCA 6.1 manual, ``spectroscopyproperties/magnx.txt``,
section 5.27 "Broken-Symmetry Wavefunctions and Exchange Couplings"):

* The phenomenological Hamiltonian is the Heisenberg-Dirac-van Vleck form
  (eq. 5.167)
  ::

      H(HDvV) = -2 J_AB (S_A . S_B)

  with ``J_AB < 0`` the antiferromagnetic (BS lowest) case (lines 9-38).
  For ``N_A`` unpaired electrons localised on site A and ``N_B`` on site B,
  the high-spin state has ``S = (N_A + N_B)/2`` and the broken-symmetry
  determinant has ``M_S = (N_A - N_B)/2`` — the manual describes the method
  as "two separate spin-unrestricted SCF calculations: (a) the calculation
  for the high-spin state ... and (b) the 'broken symmetry' calculation"
  (lines 56-66).  "Make sure that in your input coordinates 'site A' is the
  site that contains the larger number of unpaired electrons" (lines 127-129)
  is reported as a warning, not a refusal: the manual only cautions.
* The keyword is ``%scf BrokenSym NA,NB end`` (lines 114-124).  The
  alternative ``FlipSpin`` mechanism exchanges the alpha/beta spin blocks of
  the density on listed atoms *after* converging the high-spin wavefunction
  and needs ``FinalMs``; "This value MUST be given since the program cannot
  determine it by itself" (lines 156-186).
* Several formulae extract ``J_AB`` from ``E_HS``/``E_BS`` (eqs. 5.168-5.170):
  J(1) = -(E_HS - E_BS)/S_max^2 (Noodleman/Ginsberg, weak-coupling limit),
  J(2) = -(E_HS - E_BS)/(S_max(S_max+1)) (Bencini/Gatteschi, strong-coupling
  limit) and J(3), Yamaguchi's approximate spin projection.  All three are
  reported; the preferred one is Yamaguchi's approximate spin projection
  (section 5.27.1, eq. 5.170), which ORCA itself prefers: "We prefer the
  last definition because it is approximately valid over the whole coupling
  strength regime while the first equation implies the weak coupling limit
  and the second the strong coupling limit" (lines 107-112):
  ::

      J_AB = -(E_HS - E_BS) / ( <S^2>_HS - <S^2>_BS )

  The method is attributed by ORCA to Yamaguchi and co-workers (section
  5.27.1 "Approximate Spin Projection Method", citations [728, 730, 731];
  the human-facing form of the formula is the standard Yamaguchi 1986
  result, e.g. Yamaguchi, Int. J. Quantum Chem. 1986, 29, 521).  When the
  BS run carries ``APMethod 3`` (``3 = Yamaguchi``, manual lines 284-288)
  ORCA itself prints a "Spin-Hamiltonian Analysis based on
  H(HDvV)= -2J*SA*SB" block whose ``J(3)`` is this same Yamaguchi value;
  the module surfaces it as a cross-check alongside its own computation.

This module runs **one** ORCA job: the BS input makes ORCA converge the
high-spin state first and then the broken-symmetry determinant (the manual's
"two separate SCF calculations" happen inside that run).  E(High-Spin),
E(BrokenSym), <S**2> of both and E_HS - E_BS in cm**-1 are read from the
run's own summary block.  The FINAL SINGLE POINT ENERGY of an ``APMethod``
run is NOT E(BrokenSym) (real ORCA 6.1.1, stretched H2: FINAL SP -0.999153,
E(BrokenSym) -0.996998 Eh) and must not enter J.  ``J`` is reported in
``cm**-1``.  Nothing in the existing modules is changed.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec
from delfin.tools._types import ErrorKind, StepStatus

#: 1 Eh in reciprocal centimetres (fundamental constant; hartree as an
#: energy is equivalent to 219474.6313702 cm-1 when h, c and k are exact).
_EH_TO_CM1 = 219474.6313702


@dataclass
class BsOutput:
    """Everything a broken-symmetry pair yields, as typed values with units.

    ``ok`` is True only when both the high-spin and the broken-symmetry runs
    finished normally (``meaningful``); a failed either-side keeps the
    numeric fields ``None`` and surfaces ``failure_lines`` / ``error``.
    ``j_yamaguchi_cm`` is the exchange coupling from Eq. 5.170 and
    ``orca_j3_cm`` the same value as ORCA's own Spin-Hamiltonian Analysis
    J(3) block (present only when the BS run could report it).
    """

    ok: bool
    work_dir: Optional[Path] = None
    output_file: Optional[Path] = None
    meaningful: bool = False
    #: high-spin reference energy, Eh
    hs_energy_Eh: Optional[float] = None
    #: broken-symmetry determinant energy, Eh
    bs_energy_Eh: Optional[float] = None
    #: high-spin <S^2> (dimensionless)
    hs_s2: Optional[float] = None
    #: broken-symmetry <S^2> (dimensionless)
    bs_s2: Optional[float] = None
    #: E_HS - E_BS as ORCA prints it, cm**-1
    delta_e_cm: Optional[float] = None
    #: S of the high-spin state (S_max)
    s_max: Optional[float] = None
    #: J(1), eq. 5.168, weak-coupling limit, cm**-1
    j_noodleman_cm: Optional[float] = None
    #: J(2), eq. 5.169, strong-coupling limit, cm**-1
    j_bencini_cm: Optional[float] = None
    #: J(3), eq. 5.170, Yamaguchi (preferred over the whole regime), cm**-1
    j_yamaguchi_cm: Optional[float] = None
    #: J(1)/J(2)/J(3) as ORCA's own Spin-Hamiltonian Analysis prints them
    orca_j1_cm: Optional[float] = None
    orca_j2_cm: Optional[float] = None
    orca_j3_cm: Optional[float] = None
    findings: tuple[dict, ...] = ()
    failure_lines: tuple[str, ...] = ()
    error: Optional[str] = None


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


def _validate_scf_options(errors: list[str], scf_options: Optional[Mapping[str, Any]]) -> None:
    """Every ``%scf`` block option must be one the knowledge store knows.

    ``brokensym``, ``flipspin``, ``finalms`` and ``apmethod`` are OPI
    ``BlockScf`` fields verified present in the pinned OPI model (the task
    checked ``BlockScf.brokensym`` / ``.flipspin`` / ``.finalms`` /
    ``.apmethod``); the knowledge store knows them as ``scf`` block options.
    """
    from delfin.orca_knowledge import known_block_option

    for option in scf_options or {}:
        if not known_block_option("scf", option):
            errors.append(
                f"scf: block option {option!r} is not one ORCA knows "
                "(OPI block model); refused."
            )


def validate_bs(
    *,
    spin_a: int,
    spin_b: int,
    mult: int,
    flip_atoms: Optional[Sequence[int]] = None,
    final_ms: Optional[int] = None,
    method: str = "B3LYP",
    functional: str = "B3LYP",
    basis: str = "def2-SVP",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    cards: Optional[Sequence[dict]] = None,
) -> tuple[list[str], list[dict]]:
    """Check a broken-symmetry request against physics and the knowledge store.

    Refused (each with its reason) when:

    * the reference is closed-shell (``mult <= 1``) — a broken-symmetry
      determinant needs an open-shell high-spin reference (manual 5.27);
    * ``spin_a + spin_b`` is not the reference's unpaired-electron count.
      The high-spin reference of multiplicity ``mult = 2S+1`` has
      ``2S = mult-1`` unpaired electrons, so ``spin_a + spin_b == mult-1``
      must hold (manual: N_A+N_B unpaired electrons, S=(N_A+N_B)/2);
    * either spin count is negative;
    * ``FlipSpin`` mode (``flip_atoms``) is used without ``final_ms`` —
      "This value MUST be given" (manual lines 184-185);
    * both the ``BrokenSym`` and the ``FlipSpin`` mechanism are requested at
      once (two conflicting SCF mechanisms);
    * a keyword / block option ORCA does not know is requested.

    Returns ``(errors, warnings)``; a non-empty ``errors`` must stop the step
    before anything is written.  Warnings are the findings of the applicable
    knowledge cards, reported and never silently dropped.
    """
    errors: list[str] = []

    if mult <= 1:
        errors.append(
            "broken symmetry needs an open-shell (multiplicity > 1) "
            f"high-spin reference, got mult={mult} (closed-shell reference "
            "refused: there is no spin to break)."
        )
    if spin_a < 0 or spin_b < 0:
        errors.append(
            f"spin counts must be non-negative; got M={spin_a}, N={spin_b}."
        )
    if errors:
        # Broken symmetry is meaningless on a closed shell; short-circuit so
        # the electron-count check (which assumes an open-shell reference) is
        # not asked to make sense of mult=1.
        pass

    if mult > 1 and spin_a + spin_b != mult - 1:
        errors.append(
            "M+N inconsistent with the electron count: a high-spin reference "
            f"of multiplicity {mult} has {mult - 1} unpaired electrons, but "
            f"BrokenSym {spin_a},{spin_b} splits only {spin_a + spin_b}; "
            "refused (would silently change the spin state)."
        )

    if flip_atoms is not None:
        if not flip_atoms:
            errors.append("FlipSpin was requested with an empty atom list.")
        if final_ms is None:
            errors.append(
                "FlipSpin needs FinalMs: ORCA 'This value MUST be given since "
                "the program cannot determine it by itself' (manual 5.27)."
            )
        if scf_options and "brokensym" in scf_options:
            errors.append(
                "both BrokenSym and FlipSpin were requested; the two SCF "
                "mechanisms conflict, choose one (refused)."
            )
    if scf_options and "flipspin" in scf_options and final_ms is None:
        errors.append(
            "FlipSpin (in %scf) needs FinalMs; ORCA 'This value MUST be "
            "given' (manual 5.27)."
        )

    # Keyword / block-option knowledge checks come after the physics so a
    # physically-wrong request is reported as such even if its spelling is
    # also unknown.
    _validated_keyword(errors, functional, "functional")
    _validated_keyword(errors, basis, "basis")
    _validated_keyword(errors, ri, "ri")
    _validated_keyword(errors, aux_basis, "aux_basis")
    _validated_keyword(errors, dispersion, "dispersion")
    if method:
        _validated_keyword(errors, method, "method")
    for kw in extra_keywords or ():
        _validated_keyword(errors, kw, "keyword")
    _validate_scf_options(errors, scf_options)

    # Findings of the applicable verified knowledge cards.
    calc = {
        "keywords": [
            w for w in (method, functional, basis, ri, aux_basis, dispersion)
            if w
        ] + [str(k) for k in (extra_keywords or ())],
        "blocks": {"scf": dict(scf_options or {})},
        "method_family": "dft",
        "job": "bs",
    }
    from delfin.orca_knowledge import findings as knowledge_findings

    warnings = knowledge_findings(calc, cards=cards)
    return errors, warnings


def _lookup(enum_cls: type, keyword: str) -> Any:
    """The OPI SimpleKeyword enum member whose ORCA keyword is ``keyword``
    (case-blind), or the bare validated string when OPI has no member for it.
    """
    from opi.input.simple_keywords.base import SimpleKeyword

    target = str(keyword).lower()
    for attr in vars(enum_cls).values():
        if isinstance(attr, SimpleKeyword) and str(attr).lower() == target:
            return attr
    return str(keyword)


def build_bs_input(
    work_dir: Path,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    spin_a: int,
    spin_b: int,
    flip_atoms: Optional[Sequence[int]] = None,
    final_ms: Optional[int] = None,
    functional: str = "B3LYP",
    basis: str = "def2-SVP",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "calc_bs",
    hs_reference: bool = False,
) -> Path:
    """Assemble and write a broken-symmetry ORCA input with OPI's classes.

    The ``!`` simple keywords are added through ``opi.input`` (typed enums),
    ``%pal`` / ``%maxcore`` map onto ``Input.ncores`` / ``Input.memory`` and
    the ``%scf`` block (``BrokenSym NA,NB`` or ``FlipSpin`` + ``FinalMs``,
    plus ``APMethod 3`` for the Yamaguchi AP) goes through the
    ``BlockScf`` pydantic block, so every typed field OPI knows is carried
    by OPI — no hand-written ``%scf`` string.  The file is written by OPI's
    ``Calculator.write_input()``.

    Raises ``ValueError`` (before writing anything) when the request is
    refused by :func:`validate_bs`.  With ``hs_reference=True`` the input is
    the identical high-spin *reference* single point (same method and
    multiplicity, no ``%scf`` broken-symmetry block) that a BrokenSym pair
    is built on top of, so the BS electron-count check is not applied to it.
    """
    from opi.core import Calculator
    from opi.input import Input
    from opi.input.blocks.block_scf import BlockScf
    from opi.input.simple_keywords.approximation import Approximation
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dft import Dft
    from opi.input.simple_keywords.dispersion_correction import (
        DispersionCorrection,
    )
    from opi.input.structures import Atom, Structure

    if not hs_reference:
        errors, _ = validate_bs(
            spin_a=spin_a, spin_b=spin_b, mult=mult,
            flip_atoms=flip_atoms, final_ms=final_ms,
            functional=functional, basis=basis, ri=ri,
            aux_basis=aux_basis, dispersion=dispersion,
            extra_keywords=extra_keywords, scf_options=scf_options,
        )
        if errors:
            raise ValueError("; ".join(errors))

    scf = dict(scf_options or {})
    if not hs_reference:
        if flip_atoms is not None:
            scf["flipspin"] = list(flip_atoms)
            scf["finalms"] = final_ms
        else:
            scf["brokensym"] = [spin_a, spin_b]
            # APMethod 3 = Yamaguchi (manual 5.27): ORCA then prints its own
            # Spin-Hamiltonian Analysis J(3) block for the cross-check.
            scf.setdefault("apmethod", 3)

    inp = Input()
    keywords: list[Any] = []
    for enum, kw in ((Dft, functional), (BasisSet, basis),
                     (Approximation, ri), (BasisSet, aux_basis),
                     (DispersionCorrection, dispersion)):
        if kw:
            keywords.append(_lookup(enum, kw))
    for kw in extra_keywords or ():
        keywords.append(str(kw))
    inp.add_simple_keywords(*keywords)
    inp.ncores = ncores
    inp.memory = maxcore
    inp.add_blocks(BlockScf(**scf))

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


_NUM = r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"


def _bs_summary(text: str) -> dict:
    """The last broken-symmetry summary block of an ORCA ``.out``.

    ORCA prints, after the BS SCF::

        S(High-Spin)      =   1.0
        <S**2>(High-Spin) =   2.0000
        <S**2>(BrokenSym) =   0.9845
        E(High-Spin)      =    -0.994776 Eh
        E(BrokenSym)      =    -0.996998 Eh
        E(High-Spin)-E(BrokenSym)= 0.0605 eV    487.723 cm**-1 (...)

    Returns the values found (keys s_max, hs_s2, bs_s2, hs_e, bs_e,
    delta_cm); a missing line leaves its key out.
    """
    import re

    pats = {
        "s_max": r"S\(High-Spin\)\s*=\s*" + _NUM,
        "hs_s2": r"<S\*\*2>\(High-Spin\)\s*=\s*" + _NUM,
        "bs_s2": r"<S\*\*2>\(BrokenSym\)\s*=\s*" + _NUM,
        "hs_e": r"E\(High-Spin\)\s*=\s*" + _NUM + r"\s*Eh",
        "bs_e": r"E\(BrokenSym\)\s*=\s*" + _NUM + r"\s*Eh",
        "delta_cm": r"E\(High-Spin\)-E\(BrokenSym\)\s*=.*?" + _NUM + r"\s*cm\*\*-1",
    }
    out = {}
    for key, pat in pats.items():
        vals = re.findall(pat, text)
        if vals:
            out[key] = float(vals[-1])
    return out


def read_bs_output(work_dir: Path, base_name: str) -> BsOutput:
    """Read a finished broken-symmetry run (one ORCA job, ``APMethod 3``).

    E(High-Spin), E(BrokenSym), both <S**2>, S_max and E_HS - E_BS (cm**-1)
    come from ORCA's own summary block; J(1)/J(2)/J(3) (eqs. 5.168-5.170)
    are computed from the cm**-1 difference and ORCA's printed J values are
    surfaced next to them as a cross-check.
    """
    from delfin.parser import extract_last_J_block

    out_path = work_dir / f"{base_name}.out"
    text = out_path.read_text(errors="replace") if out_path.is_file() else ""
    terminated = "ORCA TERMINATED NORMALLY" in text
    summ = _bs_summary(text)
    failure_lines = _failure_lines(out_path)
    needed = ("s_max", "hs_s2", "bs_s2", "hs_e", "bs_e", "delta_cm")
    if not terminated or any(k not in summ for k in needed):
        return BsOutput(
            ok=False, work_dir=work_dir, output_file=out_path,
            meaningful=False, failure_lines=failure_lines,
            error=(
                "broken-symmetry run not readable: "
                + ("ORCA did not terminate normally" if not terminated else
                   "summary block lacks " + ", ".join(
                       k for k in needed if k not in summ))
            ),
        )
    d_cm = summ["delta_cm"]
    # the printed cm**-1 value is |E_HS - E_BS|; the sign comes from the
    # energies themselves
    if summ["hs_e"] < summ["bs_e"]:
        d_cm = -d_cm
    s_max = summ["s_max"]
    ds2 = summ["hs_s2"] - summ["bs_s2"]
    if s_max <= 0 or ds2 == 0:
        return BsOutput(
            ok=False, work_dir=work_dir, output_file=out_path,
            meaningful=False, failure_lines=failure_lines,
            error=f"J undefined: S_max={s_max}, <S**2>HS-<S**2>BS={ds2}",
        )
    orca = {}
    try:
        orca = extract_last_J_block(str(out_path), raise_on_missing=False) or {}
    except Exception:
        orca = {}
    return BsOutput(
        ok=True, work_dir=work_dir, output_file=out_path, meaningful=True,
        hs_energy_Eh=summ["hs_e"], bs_energy_Eh=summ["bs_e"],
        hs_s2=summ["hs_s2"], bs_s2=summ["bs_s2"],
        delta_e_cm=d_cm, s_max=s_max,
        j_noodleman_cm=-d_cm / s_max ** 2,
        j_bencini_cm=-d_cm / (s_max * (s_max + 1)),
        j_yamaguchi_cm=-d_cm / ds2,
        orca_j1_cm=orca.get("J1"), orca_j2_cm=orca.get("J2"),
        orca_j3_cm=orca.get("J3"),
        failure_lines=failure_lines,
    )


def _failure_lines(out_path: Path) -> tuple[str, ...]:
    """The ``.out`` lines that explain a failure (OPI/ORCA does not give
    them through the reader)."""
    if not out_path.is_file():
        return ()
    marks = ("SCF NOT CONVERGED", "ERROR", "CONVERGENCE FAILURE",
             "aborted", "out-of-memory", "not converged", "FATAL")
    return tuple(
        ln.strip() for ln in out_path.read_text(errors="replace").splitlines()
        if any(m in ln for m in marks)
    )[:40]


def run_bs(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    spin_a: int,
    spin_b: int,
    flip_atoms: Optional[Sequence[int]] = None,
    final_ms: Optional[int] = None,
    functional: str = "B3LYP",
    basis: str = "def2-SVP",
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
) -> BsOutput:
    """Run the high-spin reference and the broken-symmetry single point.

    Steps, in order:

    1. validate the physics + every keyword / block option (abort before any
       file is written when refused),
    2. build and write the input with OPI's classes,
    3. launch the BS job through ``delfin.orca.run_orca`` (ORCA converges
       the high-spin state and then the BS determinant in that one job),
    4. read ORCA's own summary block, compute J(1)/J(2)/J(3) and surface
       ORCA's printed J values as a cross-check.  ``work_dir`` must already exist; when it is ``None`` a fresh
    directory is created under the current one.
    """
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else _make_workdir()
    wdir.mkdir(parents=True, exist_ok=True)

    errors, findings = validate_bs(
        spin_a=spin_a, spin_b=spin_b, mult=mult, flip_atoms=flip_atoms,
        final_ms=final_ms, functional=functional, basis=basis, ri=ri,
        aux_basis=aux_basis, dispersion=dispersion,
        extra_keywords=extra_keywords, scf_options=scf_options,
    )
    if errors:
        return BsOutput(
            ok=False, work_dir=wdir, meaningful=False,
            findings=tuple(findings), error="; ".join(errors),
        )

    common = dict(
        elements=elements, coordinates=coordinates, charge=charge, mult=mult,
        functional=functional, basis=basis, ri=ri, aux_basis=aux_basis,
        dispersion=dispersion, extra_keywords=extra_keywords,
        maxcore=maxcore, ncores=ncores,
    )

    # One ORCA job: the BS input converges the high-spin state first and
    # then the broken-symmetry determinant (APMethod 3), and prints both.
    bs_name = f"{base_name}_bs"
    build_bs_input(wdir, spin_a=spin_a, spin_b=spin_b,
                   flip_atoms=flip_atoms, final_ms=final_ms,
                   scf_options=scf_options, base_name=bs_name, **common)
    bs_out = wdir / f"{bs_name}.out"
    run_orca(str(wdir / f"{bs_name}.inp"), str(bs_out), timeout=timeout,
             working_dir=wdir)

    result = read_bs_output(wdir, bs_name)
    result.findings = tuple(findings)
    result.work_dir = wdir
    return result


def _make_workdir() -> Path:
    import tempfile

    return Path(tempfile.mkdtemp(prefix="orca_bs_opi_", dir=Path.cwd()))


# ---------------------------------------------------------------------------
# Adapter for delfin.tools (registered as ``orca_bs_opi``)
# ---------------------------------------------------------------------------


class OrcaBsOpiAdapter(StepAdapter):
    """OPI-backed broken-symmetry single point; computes the exchange
    coupling J by the Yamaguchi formula on top of a high-spin reference."""

    name = "orca_bs_opi"
    description = "OPI-backed ORCA broken-symmetry single point (Yamaguchi J)"
    produces_geometry = False
    category = "dft"
    consumes = ("geometry",)
    produces = ("gbw", "qc_output", "property_json")
    requires_python = ("opi",)
    params = (
        ParamSpec("spin_a", "int", required=True,
                  description="M: unpaired electrons on site A for BrokenSym"),
        ParamSpec("spin_b", "int", required=True,
                  description="N: unpaired electrons on site B for BrokenSym"),
        ParamSpec("charge", "int", required=True, description="Molecular charge"),
        ParamSpec("mult", "int", default=3, description="Spin multiplicity (2S+1) of the high-spin reference"),
        ParamSpec("functional", "str", default="B3LYP", description="Functional"),
        ParamSpec("basis", "str", default="def2-SVP", description="Basis set"),
        ParamSpec("final_ms", "int", default=None, description="FinalMs for FlipSpin mode"),
        ParamSpec("flip_atoms", "list", default=None, description="FlipSpin atom indices (0-based)"),
        ParamSpec("maxcore", "int", default=1000, unit="MB", description="Memory per core"),
        ParamSpec("base_name", "str", default="calc", description="Base name for input/output files"),
    )
    data_keys = (
        DataKeySpec("hs_energy_Eh", "float", "Eh", "E(High-Spin) from ORCA's BS summary"),
        DataKeySpec("bs_energy_Eh", "float", "Eh", "E(BrokenSym) from ORCA's BS summary"),
        DataKeySpec("hs_s2", "float", "", "<S**2>(High-Spin)"),
        DataKeySpec("bs_s2", "float", "", "<S**2>(BrokenSym)"),
        DataKeySpec("delta_e_cm", "float", "cm^-1", "E_HS - E_BS"),
        DataKeySpec("j_noodleman_cm", "float", "cm^-1", "J(1), eq. 5.168, weak-coupling limit"),
        DataKeySpec("j_bencini_cm", "float", "cm^-1", "J(2), eq. 5.169, strong-coupling limit"),
        DataKeySpec("j_yamaguchi_cm", "float", "cm^-1", "J(3), eq. 5.170, Yamaguchi (preferred)"),
        DataKeySpec("orca_j3_cm", "float", "cm^-1", "J(3) as ORCA prints it (cross-check)"),
    )

    def validate_params(self, **kwargs: Any) -> None:
        # Presence checks only; physics (spin/electron-count consistency) is
        # validated in validate_bs and reported as a FAILED step, so a
        # physically-wrong request is refused with its reason, never "fixed".
        for req in ("spin_a", "spin_b", "charge"):
            if req not in kwargs:
                raise ValueError(f"'{req}' parameter is required")

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        self.validate_params(**kwargs)
        kwargs.pop("_prev_artifacts", None)
        start = time.monotonic()
        res = run_bs(
            elements=kwargs["elements"],
            coordinates=kwargs["coordinates"],
            charge=kwargs["charge"],
            mult=kwargs.get("mult", 3),
            spin_a=kwargs.get("spin_a", 1),
            spin_b=kwargs.get("spin_b", 1),
            flip_atoms=kwargs.get("flip_atoms"),
            final_ms=kwargs.get("final_ms"),
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
                "hs_energy_Eh": res.hs_energy_Eh,
                "bs_energy_Eh": res.bs_energy_Eh,
                "hs_s2": res.hs_s2,
                "bs_s2": res.bs_s2,
                "delta_e_cm": res.delta_e_cm,
                "s_max": res.s_max,
                "j_noodleman_cm": res.j_noodleman_cm,
                "j_bencini_cm": res.j_bencini_cm,
                "j_yamaguchi_cm": res.j_yamaguchi_cm,
                "orca_j1_cm": res.orca_j1_cm,
                "orca_j2_cm": res.orca_j2_cm,
                "orca_j3_cm": res.orca_j3_cm,
                "hamiltonian": "H = -2J S_A.S_B",
                "findings": list(res.findings),
                "failure_lines": list(res.failure_lines),
            },
            error=res.error,
            error_kind=(ErrorKind.NONE if res.ok
                        else ErrorKind.CONVERGENCE if res.failure_lines
                        else ErrorKind.TOOL_FAILED),
        )


register(OrcaBsOpiAdapter())
