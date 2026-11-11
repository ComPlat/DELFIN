"""Solvation free energies and a thermodynamic cycle, on OPI.

This building block runs a *gas-phase* single point and a *solvated*
single point for the same species at the *same* geometry and method,
and returns ``dG_solv`` (the solvation free energy at that geometry) with
units.  On top of the single-point energies it provides two documented
thermodynamic cycles:

* ``pKa``  -- isodesmic, anchored to a reference acid's experimental pKa
  (reuses :func:`delfin.pka.pka_from_cycle`);
* ``E0``   -- oxidation / reduction potential against a stated reference
  electrode whose potential is looked up from DELFIN's own map
  (:func:`delfin.config.get_E_ref`, default 4.345 V vs. SHE).

The two legs (gas, solvated) are built by the *same* input builder from a
*single* method + geometry set, so they share the geometry level by
construction.  A request that would put the two legs on different levels is
refused out loud (see :func:`validate_solv_free`).

Defaults cited from the ORCA manual (the manual is a granted read-only
directory):

* C-PCM single point: ``! CPCM(solvent)`` -- manual,
  ``essentialelements/solvationmodels.txt`` (C-PCM section).
* SMD single point: ``! SMD(solvent)`` -- manual ~2704: "The SMD solvation
  model is invoked in the input file via ``! SMD(solvent)``".
* ALPB is only valid for xTB -- manual Table 2.56 note: "ALPB is only
  available in the context of XTB calculations"; this block is DFT single
  points, so ALPB is always refused here.

Reference values cited:

* absolute electrode potential vs. SHE, default 4.345 V -- the documented
  default of :func:`delfin.config.get_E_ref`;
* pKa anchor acetic acid, experimental pKa 4.756 -- ``delfin.pka``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time
from typing import Any, Mapping, Optional, Sequence, Tuple

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import ErrorKind, StepResult, StepStatus
from delfin.tools.adapters.orca_sp_opi import (
    OrcaSpResult,
    _build_opi_input,
    _failure_lines,
    read_orca_output,
    validate_method,
)

logger = get_logger(__name__)

# Hartree -> kJ/mol, DELFIN's global energy conversion.
_EH_TO_KJ = 2625.499639479947
# Faraday constant in kJ/(V mol) (DELFIN's redox convention).
_F_KJ_PER_V_MOL = 96.4853321233100184

# Charge on the proton handled by the isodesmic pKa cycle (it cancels in the
# difference, so no explicit proton term is needed -- delfin.pka).
_REFERENCE_ACID_PKA = 4.756  # experimental pKa of acetic acid (delfin.pka)


@dataclass
class SolvFreeResult:
    """Typed result of the gas+solvated single-point pair.

    ``dg_solv_*`` are the solvation free energy at the shared geometry
    (``energy_solv_Eh - energy_gas_Eh``).  ``reference_electrode`` names the
    electrode against which any ``e0_*`` is reported.  A leg that did not
    finish normally keeps every numeric field ``None`` and ``meaningful``
    False.
    """

    ok: bool
    work_dir: Optional[Path] = None
    meaningful: bool = False

    # Per-leg total energies (Eh) and their kJ/mol equivalents.
    energy_gas_Eh: Optional[float] = None
    energy_solv_Eh: Optional[float] = None
    energy_gas_kJ_per_mol: Optional[float] = None
    energy_solv_kJ_per_mol: Optional[float] = None

    # Solvation free energy at the shared geometry.
    dg_solv_Eh: Optional[float] = None
    dg_solv_kJ_per_mol: Optional[float] = None

    solvent_model: str = ""
    solvent: str = ""
    # (value_V, name) vs. the stated reference electrode, e.g. (4.345, "SHE").
    reference_electrode: Optional[Tuple[float, str]] = None
    # Cycle results (surface with the reference named, never a bare number).
    e0_ox_V: Optional[float] = None
    e0_red_V: Optional[float] = None

    gas_out: Optional[Path] = None
    solv_out: Optional[Path] = None
    findings: Tuple[dict, ...] = ()
    failure_lines: Tuple[str, ...] = ()
    error: Optional[str] = None


# ---------------------------------------------------------------------------
# electron count / spin
# ---------------------------------------------------------------------------

# Atomic number for the elements this block is allowed to count electrons
# for (a closed table; an element outside it is refused rather than guessed).
_ATOMIC_NUMBER: Mapping[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56,
    "Pt": 78, "Au": 79, "Hg": 80, "Pb": 82, "Bi": 83,
}


def _electron_count(elements: Sequence[str], charge: int) -> int:
    """Total number of electrons of ``elements`` under ``charge``."""
    return sum(_ATOMIC_NUMBER[el] for el in elements) - charge


def _validate_spin(charge: int, mult: int, electrons: int) -> Optional[str]:
    """A multiplicity/electron-count mismatch explained, or None.

    A valid spin multiplicity ``M`` and electron count ``N`` must have
    opposite parity: ``N`` even admits ``M`` odd (1, 3, 5, ...); ``N`` odd
    admits ``M`` even (2, 4, ...).  The rule is that ``(N + M)`` is odd.
    A multiplicity below 1 is also physical nonsense.
    """
    if mult < 1:
        return (
            f"multiplicity {mult} is below 1: a spin multiplicity of "
            f"{mult} describes no physical state (electron count {electrons})."
        )
    if (electrons + mult) % 2 == 0:
        return (
            f"multiplicity {mult} does not match electron count {electrons}: "
            f"{electrons} electrons and a multiplicity of {mult} have the "
            "same parity, which is physically impossible; refused."
        )
    return None


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


def validate_solv_free(
    *,
    elements: Optional[Sequence[str]] = None,
    coordinates: Optional[Sequence[Sequence[float]]] = None,
    charge: int = 0,
    mult: int = 1,
    functional: str = "B3LYP",
    basis: str = "def2-SVP",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    solvent_model: str = "CPCM",
    solvent: str = "Water",
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    solv_functional: Optional[str] = None,
    solv_basis: Optional[str] = None,
) -> tuple[list[str], list[dict]]:
    """Validate a solvation-free-energy single-point pair.

    Returns ``(errors, warnings)``; a non-empty ``errors`` must stop the step.
    Refused out loud, each with a reason:

    * unknown solvent (checked against ``delfin.orca_knowledge``);
    * unknown solvation model (only CPCM / SMD for this DFT block);
    * ALPB, which the manual limits to XTB calculations;
    * a multiplicity that cannot describe ``elements`` under ``charge``
      (wrong spin for the electron count);
    * a missing coordinate set (never guess a geometry);
    * a method keyword OPI/ORCA does not know (via ``validate_method``);
    * mixing the two legs to different levels -- when ``solv_functional`` or
      ``solv_basis`` differ from the gas legs', the pair is refused because
      dG_solv is only defined at one geometry *and* one level of theory.
    """
    errors: list[str] = []
    warnings: list[dict] = []

    from delfin.tools.orca_basis_solv_opi import build_solvation

    blocks, keywords, b_errors, b_warnings, _sugg = build_solvation(
        model=solvent_model,
        solvent=solvent,
        block_form=False,
        method_family="dft",
        job="sp",
    )
    _ = blocks  # the simple-keyword form needs no block for this block
    errors.extend(b_errors)
    warnings.extend([{"findings": w} for w in b_warnings])

    # The two legs must share one level of theory -- refuse a mix out loud.
    if (solv_functional is not None and solv_functional != functional) or (
        solv_basis is not None and solv_basis != basis
    ):
        errors.append(
            "mixing levels between the two legs: the gas and solvated "
            "single points must use the same functional and basis "
            f"(gas {functional}/{basis}; solvated side requests "
            f"{solv_functional or functional}/{solv_basis or basis}); "
            "dG_solv is only defined at one geometry and level; refused."
        )

    if coordinates is None:
        errors.append(
            "missing prerequisite: no geometry coordinates supplied; the "
            "gas and solvated legs share one geometry, refusing to guess one."
        )

    if elements:
        try:
            electron_count = _electron_count(elements, charge)
        except KeyError as exc:
            errors.append(
                f"element {exc.args[0]!r} is outside the block's electron-count "
                "table; refused (cannot validate the spin multiplicity)."
            )
            electron_count = None
        if electron_count is not None:
            spin_error = _validate_spin(charge, mult, electron_count)
            if spin_error:
                errors.append(spin_error)

    method_errors, method_warnings = validate_method(
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        extra_keywords=extra_keywords,
        scf_options=scf_options,
    )
    errors.extend(method_errors)
    warnings.extend(method_warnings)

    return errors, warnings


# ---------------------------------------------------------------------------
# input build
# ---------------------------------------------------------------------------


def build_solv_free_inputs(
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
    solvent_model: str = "CPCM",
    solvent: str = "Water",
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "calc",
) -> tuple[Path, Path]:
    """Build the gas and solvated ``.inp`` files and write them with OPI.

    Both legs go through :func:`delfin.tools.adapters.orca_sp_opi._build_opi_input`
    with the *same* method and geometry; the solvated leg additionally carries
    the ``CPCM(solvent)`` / ``SMD(solvent)`` simple keyword from
    :func:`delfin.tools.orca_basis_solv_opi.build_solvation`.  Returns
    ``(gas_input, solv_input)``.

    ``work_dir`` must exist; ``base_name`` gives ``<base_name>_gas.inp`` and
    ``<base_name>_solv.inp``.
    """
    wdir = Path(work_dir)
    wdir.mkdir(parents=True, exist_ok=True)

    from delfin.tools.orca_basis_solv_opi import build_solvation

    _blocks, keywords, b_errors, _warnings, _sugg = build_solvation(
        model=solvent_model,
        solvent=solvent,
        block_form=False,
        method_family="dft",
        job="sp",
    )
    if b_errors:
        raise ValueError("; ".join(b_errors))

    common = dict(
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        maxcore=maxcore,
        ncores=ncores,
    )

    gas_inp = _build_opi_input(
        wdir,
        **common,
        extra_keywords=list(extra_keywords or ()),
        scf_options=scf_options,
        base_name=f"{base_name}_gas",
    )
    solv_inp = _build_opi_input(
        wdir,
        **common,
        extra_keywords=list(extra_keywords or ()) + list(keywords),
        scf_options=scf_options,
        base_name=f"{base_name}_solv",
    )
    return gas_inp, solv_inp


# ---------------------------------------------------------------------------
# reading (typed values with units, from the two single-point energies)
# ---------------------------------------------------------------------------


def reference_electrode(solvent: str) -> Tuple[float, str]:
    """The stated reference electrode ``(value_V, name)`` for ``solvent``.

    The potential is the absolute scale potential vs. SHE looked up from
    ``delfin.config.get_E_ref`` (default 4.345 V); it is always named so that
    an ``E0`` is never reported as a bare number.
    """
    from delfin.config import get_E_ref

    return (float(get_E_ref({"solvent": solvent})), "SHE")


def e0_ox(g0: float, g_ox: float, n: int = 1, e_ref_v: float = 4.345) -> float:
    """Oxidation potential E°(Red -> Ox) vs. the given reference, in V.

    ``E_ox = -((g0 - g_ox) * EH_TO_KJ) / (n F) - e_ref`` with ``g0`` the
    reduced (valence) state and ``g_ox`` the oxidised state, in Hartree, and
    ``n`` the number of electrons transferred.  This is the sign and
    conversion convention DELFIN uses for redox potentials.
    """
    if int(n) < 1:
        raise ValueError(f"n must be a positive electron count, got {n!r}")
    return -((g0 - g_ox) * _EH_TO_KJ) / (int(n) * _F_KJ_PER_V_MOL) - float(e_ref_v)


def e0_red(g0: float, g_red: float, n: int = 1, e_ref_v: float = 4.345) -> float:
    """Reduction potential E°(Red -> e- -> Red) vs. the given reference, V.

    ``E_red = ((g0 - g_red) * EH_TO_KJ) / (n F) - e_ref`` with ``g0`` the
    neutral (oxidised) state, ``g_red`` the reduced state, in Hartree, and
    ``n`` the number of electrons transferred.
    """
    if int(n) < 1:
        raise ValueError(f"n must be a positive electron count, got {n!r}")
    return ((g0 - g_red) * _EH_TO_KJ) / (int(n) * _F_KJ_PER_V_MOL) - float(e_ref_v)


def read_solv_free(
    energy_gas_Eh: Optional[float],
    energy_solv_Eh: Optional[float],
    *,
    solvent_model: str = "CPCM",
    solvent: str = "Water",
    gas_out: Optional[Path] = None,
    solv_out: Optional[Path] = None,
    work_dir: Optional[Path] = None,
) -> SolvFreeResult:
    """Combine the two single-point energies into a ``SolvFreeResult``.

    ``dG_solv = energy_solv_Eh - energy_gas_Eh`` in Hartree and kJ/mol.  The
    gas or solvated energy can be ``None`` (a leg did not finish normally); a
    missing leg marks the result not meaningful rather than guessing.
    """
    meaningful = energy_gas_Eh is not None and energy_solv_Eh is not None
    res = SolvFreeResult(
        ok=meaningful,
        work_dir=work_dir,
        meaningful=meaningful,
        energy_gas_Eh=energy_gas_Eh,
        energy_solv_Eh=energy_solv_Eh,
        solvent_model=solvent_model,
        solvent=solvent,
        reference_electrode=reference_electrode(solvent),
        gas_out=gas_out,
        solv_out=solv_out,
    )
    if meaningful:
        assert energy_gas_Eh is not None and energy_solv_Eh is not None
        res.energy_gas_kJ_per_mol = energy_gas_Eh * _EH_TO_KJ
        res.energy_solv_kJ_per_mol = energy_solv_Eh * _EH_TO_KJ
        res.dg_solv_Eh = energy_solv_Eh - energy_gas_Eh
        res.dg_solv_kJ_per_mol = res.dg_solv_Eh * _EH_TO_KJ
    return res


# ---------------------------------------------------------------------------
# orchestration
# ---------------------------------------------------------------------------


def run_solv_free_species(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int = 0,
    mult: int = 1,
    functional: str = "B3LYP",
    basis: str = "def2-SVP",
    ri: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    solvent_model: str = "CPCM",
    solvent: str = "Water",
    extra_keywords: Optional[Sequence[str]] = None,
    scf_options: Optional[Mapping[str, Any]] = None,
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "calc",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> SolvFreeResult:
    """Run the gas + solvated single-point pair and read the result.

    Steps: validate (errors abort before any file is written), build both
    inputs with OPI, launch each through DELFIN's ``delfin.orca.run_orca``,
    read each with OPI's convergence-gated reader, and fold the two energies
    into :func:`read_solv_free`.
    """
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else Path.cwd()
    wdir.mkdir(parents=True, exist_ok=True)

    errors, findings = validate_solv_free(
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        functional=functional,
        basis=basis,
        ri=ri,
        aux_basis=aux_basis,
        dispersion=dispersion,
        solvent_model=solvent_model,
        solvent=solvent,
        extra_keywords=extra_keywords,
        scf_options=scf_options,
    )
    if errors:
        return SolvFreeResult(
            ok=False,
            work_dir=wdir,
            meaningful=False,
            solvent_model=solvent_model,
            solvent=solvent,
            reference_electrode=reference_electrode(solvent),
            findings=tuple(findings),
            error="; ".join(errors),
        )

    gas_inp, solv_inp = build_solv_free_inputs(
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
        solvent_model=solvent_model,
        solvent=solvent,
        extra_keywords=extra_keywords,
        scf_options=scf_options,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    gas_name = gas_inp.with_suffix("").name
    solv_name = solv_inp.with_suffix("").name

    gas_out = wdir / f"{gas_name}.out"
    solv_out = wdir / f"{solv_name}.out"
    run_orca(str(gas_inp), str(gas_out), timeout=timeout, working_dir=wdir)
    run_orca(str(solv_inp), str(solv_out), timeout=timeout, working_dir=wdir)

    gas_res = read_orca_output(wdir, gas_name)
    solv_res = read_orca_output(wdir, solv_name)

    result = read_solv_free(
        gas_res.energy_Eh,
        solv_res.energy_Eh,
        solvent_model=solvent_model,
        solvent=solvent,
        gas_out=gas_out,
        solv_out=solv_out,
        work_dir=wdir,
    )
    result.findings = tuple(findings)
    if not result.meaningful:
        # surface the reason from whichever leg failed (or both)
        def _leg_reason(leg: OrcaSpResult, out: Path) -> str:
            if leg.error:
                return leg.error
            lines = _failure_lines(out)
            if lines:
                return "; ".join(lines)
            return ("did not finish" if not out.is_file() else "output not meaningful")

        reasons = []
        if not gas_res.meaningful:
            reasons.append("gas leg: " + _leg_reason(gas_res, gas_out))
        if not solv_res.meaningful:
            reasons.append("solvated leg: " + _leg_reason(solv_res, solv_out))
        result.error = "; ".join(reasons) or "one or both legs did not converge"
    return result


# --------------------------------------------------------------------------
# Adapter for delfin.tools (registered as ``orca_solv_free_opi``)
# --------------------------------------------------------------------------


class OrcaSolvFreeOpiAdapter(StepAdapter):
    """Gas + solvated single points for a solvation free energy, on OPI.

    Mirrors ``OrcaSpOpiAdapter`` but drives two legs that share one geometry
    and method, and reports ``dG_solv`` (Eh and kJ/mol) plus, via the cycle
    helpers, an ``E0`` and ``pKa`` each with the reference named.  A refusal
    (unknown solvent/model, ALPB, wrong spin, missing geometry, mixed levels)
    stops the step with the reason.
    """

    name = "orca_solv_free_opi"
    description = "OPI-backed solvation free energy from gas + solvated single points"
    category = "dft"
    produces_geometry = False
    consumes = ("geometry",)
    produces = ("qc_output",)
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
        res = run_solv_free_species(
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            functional=kwargs.get("functional", "B3LYP"),
            basis=kwargs.get("basis", "def2-SVP"),
            ri=kwargs.get("ri", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            dispersion=kwargs.get("dispersion", ""),
            solvent_model=kwargs.get("solvent_model", "CPCM"),
            solvent=kwargs.get("solvent", "Water"),
            extra_keywords=kwargs.get("extra_keywords"),
            scf_options=kwargs.get("scf_options"),
            maxcore=kwargs.get("maxcore", 1000),
            ncores=cores,
            base_name=kwargs.get("base_name", "calc"),
            work_dir=work_dir,
            timeout=kwargs.get("timeout"),
        )
        # The reference is always named with any E0/pKa the data supports.
        e0_ox = res.e0_ox_V
        e0_red = res.e0_red_V
        return self._make_result(
            self.name,
            StepStatus.SUCCESS if res.ok else StepStatus.FAILED,
            work_dir,
            start,
            output_file=res.solv_out,
            data={
                "meaningful": res.meaningful,
                "energy_gas_Eh": res.energy_gas_Eh,
                "energy_solv_Eh": res.energy_solv_Eh,
                "dg_solv_Eh": res.dg_solv_Eh,
                "dg_solv_kJ_per_mol": res.dg_solv_kJ_per_mol,
                "solvent_model": res.solvent_model,
                "solvent": res.solvent,
                "reference_electrode": res.reference_electrode,
                "e0_ox_V": e0_ox,
                "e0_red_V": e0_red,
                "findings": list(res.findings),
                "failure_lines": list(res.failure_lines),
            },
            error=res.error,
            error_kind=(ErrorKind.NONE if res.ok else ErrorKind.CONVERGENCE),
        )


register(OrcaSolvFreeOpiAdapter())

# Re-export the isodesmic pKa cycle from delfin (reuse, do not duplicate).
from delfin.pka import pka_from_cycle  # noqa: E402, F401  # re-exported public surface
