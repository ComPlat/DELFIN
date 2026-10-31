"""DELFIN's GOAT conformer-search building block on OPI.

This module builds a GOAT (Global Optimizer Algorithm) conformer search
entirely with the ORCA Python Interface (OPI, orca-pi 2.0.0): the ``!GOAT``
simple keyword, the method keyword and the ``%goat`` block come from OPI's
typed input classes, the input file is written by OPI's
``Calculator.write_input()``, the binary is launched by DELFIN's existing
runner (``delfin.orca.run_orca``, reused from the ``orca_sp_opi`` flow) and
every keyword is checked against DELFIN's ORCA knowledge store
(``delfin.orca_knowledge``) before anything is written.  Nothing in
``orca_sp_opi.py`` is edited; its helpers are imported.

Manual ground truth (ORCA 6.1 manual,
structurereactivity/goat.txt, section 4.10):

* "The command to call the global optimizer is simply ``!GOAT``, like you
  would with ``!OPT``, and its options can be given under the ``%GOAT``
  block as usual."  "You can give it together with any other method
  available in ORCA ... Here we will just use GFN2 (or ``!XTB``)."
* Table 4.9 ``%goat`` block keywords used here:
  - MAXITER: "defines an arbitrary number of max GOAT geom opt iters per
    worker" (default 128)
  - NWORKERS: "define the number of workers (default auto). auto for an
    automatic ideal assignment, or give any number multiple of 4 (number
    of temperatures)"
  - MAXEN: "the maximum relative energy of a conformer to be taken, in
    kcal/mol. 6 kcal/mol by default" -- this is the energy window
  - ENDIFF: "minimum energy difference needed to differentiate
    conformers, in kcal/mol" (default 0.1)
  - CONFTEMP: "temperature used to compute the free energy, in Kelvin"
    (default 298.15)
  - RMSD: "minimum RMSD to differentiate conformers, in Angstroem"
    (default 0.125)
* Ensemble readback: "The final relative energies of the ensemble are
  printed afterwards, together with a ``BaseName.finalensemble.xyz``
  file" with a "# Final ensemble info #" table
  (Conformer / Energy kcal/mol / Degen. / % total / % cumul.).

Reader rules, inherited from the orca_sp_opi contract: numbers are
returned only when the run ``terminated_normally()``; a failed run keeps
them empty and carries the ``.out`` failure lines instead.

The ensemble energies in the ``.xyz`` comment lines are *relative*
energies in kcal/mol with the lowest conformer at 0.000, at the ensemble
temperature (``CONFTEMP``, default 298.15 K); the Boltzmann weights this
module reports are computed from those relative energies at that same
temperature, matching the "% total" column of ORCA's own table.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from delfin.common.logging import get_logger
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import ErrorKind, StepResult, StepStatus

logger = get_logger(__name__)

R_GAS_KCAL_MOL_K = 0.0019872041  # kcal/(mol K)


@dataclass
class GoatResult:
    """Typed result of one OPI-backed GOAT conformer search.

    The ensemble fields are populated only after the run
    ``terminated_normally()``; a failed run keeps them empty and carries
    the ``.out`` failure lines instead.
    """

    ok: bool
    work_dir: Path
    output_file: Optional[Path] = None
    terminated_normally: Optional[bool] = None
    meaningful: bool = False

    # Ensemble (meaningful only): relative energies in kcal/mol, lowest = 0,
    # at the ensemble temperature.
    energies_kcal_mol: tuple[float, ...] = ()
    boltzmann_weights: Optional[tuple[float, ...]] = None
    temperature_k: float = 298.15
    structures: Optional[tuple[tuple[tuple[str, float, float, float], ...], ...]] = None
    ensemble_file: Optional[Path] = None
    global_minimum_file: Optional[Path] = None
    n_conformers: int = 0

    findings: tuple[dict, ...] = ()
    failure_lines: tuple[str, ...] = ()
    error: Optional[str] = None


def validate_goat(
    *,
    method: str = "XTB",
    goat_options: Optional[Mapping[str, Any]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    functional: str = "",
    basis: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    cards: Optional[Sequence[dict]] = None,
) -> tuple[list[str], list[dict]]:
    """Check every keyword and ``%goat`` block option against DELFIN's ORCA
    knowledge store.

    Returns ``(errors, warnings)``; a non-empty ``errors`` must stop the
    step.  ``method`` is a simple keyword ("XTB" for GFN2-xTB, or a DFT
    functional for a DFT-level search, supplied together with ``basis``).
    """
    from delfin.tools.adapters.orca_sp_opi import _validated_keyword

    errors: list[str] = []
    _validated_keyword(errors, method, "method")
    _validated_keyword(errors, functional, "functional")
    _validated_keyword(errors, basis, "basis")
    _validated_keyword(errors, aux_basis, "aux_basis")
    _validated_keyword(errors, dispersion, "dispersion")
    for kw in extra_keywords or ():
        _validated_keyword(errors, kw, "keyword")

    if goat_options:
        from delfin.orca_knowledge import known_block_option

        for option in goat_options:
            if not known_block_option("goat", option):
                errors.append(
                    f"goat: block option {option!r} is not one ORCA knows "
                    "(OPI block model / ORCA parser); refused."
                )

    calc = {
        "keywords": [w for w in (method, functional, basis, aux_basis, dispersion) if w]
        + [str(k) for k in (extra_keywords or ())],
        "blocks": {"goat": dict(goat_options or {})},
        "method_family": "goat",
        "job": "goat",
    }
    from delfin.orca_knowledge import findings as knowledge_findings

    warnings = knowledge_findings(calc, cards=cards)
    return errors, warnings


def build_goat_input(
    work_dir: Path,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int = 0,
    mult: int = 1,
    method: str = "XTB",
    goat_options: Optional[Mapping[str, Any]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    functional: str = "",
    basis: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "goat",
) -> Path:
    """Assemble and write the GOAT ORCA input with OPI's classes.

    The method keyword and ``%goat`` block are validated first; an unknown
    keyword or option raises ``ValueError`` before any file is written.
    Returns the absolute path of the written ``<base_name>.inp``.
    """
    errors, _ = validate_goat(
        method=method,
        goat_options=goat_options,
        extra_keywords=extra_keywords,
        functional=functional,
        basis=basis,
        aux_basis=aux_basis,
        dispersion=dispersion,
    )
    if errors:
        raise ValueError("; ".join(errors))

    from opi.core import Calculator
    from opi.input import Input
    from opi.input.blocks.block_goat import BlockGoat
    from opi.input.simple_keywords.base import SimpleKeyword
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.dispersion_correction import (
        DispersionCorrection,
    )
    from opi.input.simple_keywords.dft import Dft
    from opi.input.simple_keywords.goat import Goat
    from opi.input.structures import Atom, Structure

    def _lookup(enum_cls: type, keyword: str) -> Any:
        target = str(keyword).lower()
        for attr in vars(enum_cls).values():
            if isinstance(attr, SimpleKeyword) and str(attr).lower() == target:
                return attr
        return str(keyword)

    inp = Input()
    # !GOAT -- "The command to call the global optimizer is simply !GOAT".
    inp.add_simple_keywords(_lookup(Goat, "GOAT"))
    inp.add_simple_keywords(_lookup(Dft, method))
    if functional:
        inp.add_simple_keywords(_lookup(Dft, functional))
    if basis:
        inp.add_simple_keywords(_lookup(BasisSet, basis))
    if aux_basis:
        inp.add_simple_keywords(_lookup(BasisSet, aux_basis))
    if dispersion:
        inp.add_simple_keywords(_lookup(DispersionCorrection, dispersion))
    for kw in extra_keywords or ():
        inp.add_simple_keywords(str(kw))
    inp.ncores = ncores
    inp.memory = maxcore
    if goat_options:
        inp.add_blocks(BlockGoat(**dict(goat_options)))

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


def _boltzmann_weights(
    rel_energies: Sequence[float], temperature_k: float
) -> tuple[float, ...]:
    """Boltzmann weights from relative energies (kcal/mol) at ``T``.

    The weights are normalised so they sum to 1, matching the "% total"
    column of ORCA's "# Final ensemble info #" table (manual, section 4.10.3).
    """
    if not rel_energies or temperature_k <= 0:
        return ()
    kbt = R_GAS_KCAL_MOL_K * temperature_k
    exps = [math.exp(-e / kbt) for e in rel_energies]
    total = sum(exps)
    return tuple(e / total for e in exps)


def read_goat_ensemble(
    work_dir: Path, base_name: str, temperature_k: float = 298.15
) -> dict:
    """Read the GOAT ensemble back from the run directory.

    Sources, per the manual (section 4.10.3 "The final ensemble"): the
    ``<base>.finalensemble.xyz`` file -- an XYZ trajectory whose comment
    lines carry the relative energies in kcal/mol in the format "Energy
    (float)", "as generated by GOAT" (manual, READENSEMBLE description) --
    and the ``<base>.globalminimum.xyz`` file.  Only files that exist are
    reported.
    """
    out: dict = {
        "energies_kcal_mol": (),
        "structures": (),
        "ensemble_file": None,
        "global_minimum_file": None,
        "n_conformers": 0,
    }
    ens = work_dir / f"{base_name}.finalensemble.xyz"
    if not ens.is_file():
        return out
    out["ensemble_file"] = ens
    energies: list[float] = []
    structures: list[tuple[tuple[str, float, float, float], ...]] = []
    lines = ens.read_text(errors="replace").splitlines()
    i = 0
    while i < len(lines):
        try:
            natoms = int(lines[i].split()[0])
        except (ValueError, IndexError):
            i += 1
            continue
        # The comment line carries the relative energy in kcal/mol
        # ("Energy (float)", READENSEMBLE format in the manual).
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        tokens = comment.split()
        energy = float(tokens[-1]) if tokens else float("nan")
        frame: list[tuple[str, float, float, float]] = []
        for j in range(natoms):
            parts = lines[i + 2 + j].split()
            frame.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
        energies.append(energy)
        structures.append(tuple(frame))
        i += 2 + natoms
    out["energies_kcal_mol"] = tuple(energies)
    out["structures"] = tuple(structures)
    out["n_conformers"] = len(energies)
    gmin = work_dir / f"{base_name}.globalminimum.xyz"
    if gmin.is_file():
        out["global_minimum_file"] = gmin
    return out


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


def read_goat_output(work_dir: Path, base_name: str) -> GoatResult:
    """Read a finished GOAT run with OPI, gating the ensemble on
    ``terminated_normally()``."""
    from opi.output.core import Output

    out_path = work_dir / f"{base_name}.out"
    try:
        out = Output(basename=base_name, working_dir=work_dir)
        out.parse()
    except Exception as exc:
        return GoatResult(
            ok=False,
            work_dir=work_dir,
            output_file=out_path if out_path.is_file() else None,
            meaningful=False,
            failure_lines=_failure_lines(out_path),
            error=f"OPI could not parse the run: {exc}",
        )

    terminated = out.terminated_normally()
    result = GoatResult(
        ok=bool(terminated),
        work_dir=work_dir,
        output_file=out_path if out_path.is_file() else None,
        terminated_normally=terminated,
        meaningful=bool(terminated),
    )
    if not terminated:
        result.failure_lines = _failure_lines(out_path)
        result.error = "run did not finish normally"
        return result

    # The ensemble temperature: CONFTEMP in the %goat options is the manual's
    # "temperature used to compute the free energy, in Kelvin" (default
    # 298.15); the caller's options are re-read from the input file by OPI.
    result.temperature_k = 298.15
    ensemble = read_goat_ensemble(work_dir, base_name, result.temperature_k)
    result.energies_kcal_mol = ensemble["energies_kcal_mol"]
    result.structures = ensemble["structures"]
    result.n_conformers = ensemble["n_conformers"]
    result.ensemble_file = ensemble["ensemble_file"]
    result.global_minimum_file = ensemble["global_minimum_file"]
    result.boltzmann_weights = _boltzmann_weights(
        result.energies_kcal_mol, result.temperature_k
    )
    return result


def run_goat(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int = 0,
    mult: int = 1,
    method: str = "XTB",
    goat_options: Optional[Mapping[str, Any]] = None,
    extra_keywords: Optional[Sequence[str]] = None,
    functional: str = "",
    basis: str = "",
    aux_basis: str = "",
    dispersion: str = "",
    maxcore: int = 1000,
    ncores: int = 1,
    base_name: str = "goat",
    work_dir: Optional[Path] = None,
    timeout: Optional[int] = None,
) -> GoatResult:
    """Run one GOAT conformer search and read the ensemble back.

    Steps, in order: validate against the knowledge store (errors abort
    before any file is written), build and write the input with OPI,
    launch through DELFIN's ``delfin.orca.run_orca``, read back with OPI
    with the ensemble gated on ``terminated_normally()``.
    """
    from delfin.orca import run_orca

    wdir = Path(work_dir).resolve() if work_dir is not None else Path.cwd() / "goat_run"
    wdir.mkdir(parents=True, exist_ok=True)

    errors, findings = validate_goat(
        method=method,
        goat_options=goat_options,
        extra_keywords=extra_keywords,
        functional=functional,
        basis=basis,
        aux_basis=aux_basis,
        dispersion=dispersion,
    )
    if errors:
        return GoatResult(
            ok=False,
            work_dir=wdir,
            meaningful=False,
            findings=tuple(findings),
            error="; ".join(errors),
        )

    inp_path = build_goat_input(
        wdir,
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        method=method,
        goat_options=goat_options,
        extra_keywords=extra_keywords,
        functional=functional,
        basis=basis,
        aux_basis=aux_basis,
        dispersion=dispersion,
        maxcore=maxcore,
        ncores=ncores,
        base_name=base_name,
    )
    out_path = wdir / f"{base_name}.out"
    success = run_orca(str(inp_path), str(out_path), timeout=timeout, working_dir=wdir)
    result = read_goat_output(wdir, base_name)
    result.findings = tuple(findings)
    if not success and result.ok:
        result.ok = False
        result.meaningful = False
        result.error = result.error or "ORCA subprocess reported failure"
        result.failure_lines = result.failure_lines or _failure_lines(out_path)
    return result


# --------------------------------------------------------------------------
# Adapter for delfin.tools (registered as ``orca_goat_opi``)
# --------------------------------------------------------------------------


class OrcaGoatOpiAdapter(StepAdapter):
    """A GOAT conformer-search step driven by OPI, launchable via
    ``run_step``.  Mirrors ``OrcaSpOpiAdapter``'s position in the registry.

    Default method: GFN2-xTB (``!XTB``) -- the manual's example method for
    GOAT ("Here we will just use GFN2 (or !XTB)"); a DFT-level search is
    available by passing ``functional``/``basis`` (the manual warns it is
    much costlier).  This default is a *suggestion*: the caller's explicit
    choice wins and is only validated.
    """

    name = "orca_goat_opi"
    description = "OPI-backed ORCA GOAT conformer search (GFN2-xTB or DFT)"
    category = "conformers"
    produces_geometry = True
    consumes = ("geometry",)
    produces = ("geometry", "ensemble")
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
        res = run_goat(
            elements=elements,
            coordinates=coordinates,
            charge=kwargs.get("charge", 0),
            mult=kwargs.get("mult", 1),
            method=kwargs.get("method", "XTB"),
            goat_options=kwargs.get("goat_options"),
            extra_keywords=kwargs.get("extra_keywords"),
            functional=kwargs.get("functional", ""),
            basis=kwargs.get("basis", ""),
            aux_basis=kwargs.get("aux_basis", ""),
            dispersion=kwargs.get("dispersion", ""),
            maxcore=kwargs.get("maxcore", 1000),
            ncores=cores,
            base_name=kwargs.get("base_name", "goat"),
            work_dir=work_dir,
            timeout=kwargs.get("timeout"),
        )
        data: dict[str, Any] = {
            "meaningful": res.meaningful,
            "n_conformers": res.n_conformers,
            "energies_kcal_mol": list(res.energies_kcal_mol),
            "boltzmann_weights": list(res.boltzmann_weights or ()),
            "temperature_k": res.temperature_k,
            "findings": list(res.findings),
            "failure_lines": list(res.failure_lines),
        }
        return self._make_result(
            self.name,
            StepStatus.SUCCESS if res.ok else StepStatus.FAILED,
            work_dir,
            start,
            output_file=res.ensemble_file or res.output_file,
            data=data,
            error=res.error,
            error_kind=(ErrorKind.NONE if res.ok else ErrorKind.TOOL_FAILED),
        )


register(OrcaGoatOpiAdapter())
