"""Unified calculator factory for DELFIN.

Single entry point for ALL computational backends — ML potentials,
DFT, semi-empirical, force fields, and periodic codes.  Every backend
returns a standard ASE Calculator so downstream code never needs to
know which engine is running.

Usage::

    from delfin.calculators import create_calculator, available_backends

    # ML potentials (pre-trained, fast)
    calc = create_calculator("ani2x")
    calc = create_calculator("mace_off", device="cuda")

    # DFT / ab-initio (needs binary in PATH)
    calc = create_calculator("orca", method="B3LYP", basis="def2-SVP")
    calc = create_calculator("gaussian", method="B3LYP", basis="6-31G*")
    calc = create_calculator("vasp", xc="PBE", kpts=[4,4,4])

    # Semi-empirical
    calc = create_calculator("xtb", method="GFN2-xTB")
    calc = create_calculator("mopac", method="PM7")

    # Use like any ASE calculator
    atoms.calc = calc
    energy = atoms.get_potential_energy()   # eV
    forces = atoms.get_forces()            # eV/Å
"""

from __future__ import annotations

import shutil
from typing import Optional

from ase.calculators.calculator import Calculator

from delfin.common.logging import get_logger

logger = get_logger(__name__)


# ═══════════════════════════════════════════════════════════════════════
#  Backend registry
# ═══════════════════════════════════════════════════════════════════════

# alias → (canonical_key, category)
_REGISTRY: dict[str, tuple[str, str]] = {}

def _register(canonical: str, category: str, *aliases: str):
    _REGISTRY[canonical] = (canonical, category)
    for a in aliases:
        _REGISTRY[a] = (canonical, category)

# ── ML Potentials ─────────────────────────────────────────────────────
_register("ani2x",      "mlp", "ani-2x", "torchani")
_register("aimnet2",    "mlp", "aimnet")
_register("mace_off",   "mlp", "mace-off", "mace")
_register("chgnet",     "mlp")
_register("m3gnet",     "mlp", "matgl", "megnet")
_register("schnetpack", "mlp", "schnet", "painn")
_register("nequip",     "mlp", "allegro")
_register("alignn",     "mlp")

# ── DFT / ab-initio ──────────────────────────────────────────────────
_register("orca",       "qm")
_register("gaussian",   "qm", "g16", "g09")
_register("turbomole",  "qm")
_register("nwchem",     "qm")
_register("psi4",       "qm")
_register("qchem",      "qm", "q-chem")
_register("cp2k",       "qm")
_register("espresso",   "qm", "quantum-espresso", "qe", "pw")
_register("aims",       "qm", "fhi-aims")
_register("siesta",     "qm")
_register("gpaw",       "qm")
_register("gamess",     "qm")
_register("dalton",     "qm")
_register("molpro",     "qm")
_register("abinit",     "qm")
_register("crystal",    "qm")
_register("vasp",       "qm")
_register("fleur",      "qm")
_register("openmolcas",  "qm", "molcas")
_register("pyscf",      "qm")

# ── Semi-empirical ───────────────────────────────────────────────────
_register("xtb",        "semiempirical", "gfn2-xtb", "gfn-ff", "gfnff")
_register("dftb",       "semiempirical", "dftb+", "dftbplus")
_register("mopac",      "semiempirical", "pm7", "pm6")

# ── MD engines ───────────────────────────────────────────────────────
_register("lammps",     "md", "lmp")
_register("openmm",     "md")
_register("gromacs",    "md", "gmx")


# ═══════════════════════════════════════════════════════════════════════
#  Factory functions — MLP backends (delegate to mlp_tools)
# ═══════════════════════════════════════════════════════════════════════

def _create_mlp(key: str, **kwargs) -> Calculator:
    from delfin.mlp_tools.calculators import create_calculator as _mlp_create
    return _mlp_create(key, **kwargs)


# ═══════════════════════════════════════════════════════════════════════
#  Factory functions — QM backends (via ASE built-in calculators)
# ═══════════════════════════════════════════════════════════════════════

def _create_orca(method: str = "B3LYP", basis: str = "def2-SVP",
                 charge: int = 0, mult: int = 1, pal: int = 1,
                 extra_input: str = "", **kwargs) -> Calculator:
    from ase.calculators.orca import ORCA, OrcaProfile
    orca_path = shutil.which("orca")
    if not orca_path:
        raise FileNotFoundError("ORCA executable not found in PATH")
    simple_input = f"{method} {basis} {extra_input}".strip()
    profile = OrcaProfile(command=orca_path)
    calc = ORCA(
        profile=profile,
        orcasimpleinput=simple_input,
        orcablocks=f"%pal nprocs {pal} end",
        charge=charge,
        mult=mult,
        **kwargs,
    )
    logger.info("Created ORCA calculator (%s/%s, charge=%d, mult=%d)", method, basis, charge, mult)
    return calc


def _create_gaussian(method: str = "B3LYP", basis: str = "6-31G*",
                     charge: int = 0, mult: int = 1, nprocs: int = 1,
                     mem: str = "4GB", **kwargs) -> Calculator:
    from ase.calculators.gaussian import Gaussian
    calc = Gaussian(
        method=method,
        basis=basis,
        charge=charge,
        mult=mult,
        nprocshared=nprocs,
        mem=mem,
        **kwargs,
    )
    logger.info("Created Gaussian calculator (%s/%s)", method, basis)
    return calc


def _create_turbomole(**kwargs) -> Calculator:
    from ase.calculators.turbomole import Turbomole
    calc = Turbomole(**kwargs)
    logger.info("Created Turbomole calculator")
    return calc


def _create_nwchem(method: str = "dft", xc: str = "B3LYP",
                   basis: str = "6-31G*", **kwargs) -> Calculator:
    from ase.calculators.nwchem import NWChem
    calc = NWChem(
        dft={"xc": xc},
        basis=basis,
        **kwargs,
    )
    logger.info("Created NWChem calculator (%s/%s)", xc, basis)
    return calc


def _create_psi4(method: str = "b3lyp", basis: str = "6-31G*",
                 charge: int = 0, mult: int = 1, **kwargs) -> Calculator:
    from ase.calculators.psi4 import Psi4
    calc = Psi4(
        method=method,
        basis=basis,
        charge=charge,
        multiplicity=mult,
        **kwargs,
    )
    logger.info("Created Psi4 calculator (%s/%s)", method, basis)
    return calc


def _create_qchem(method: str = "B3LYP", basis: str = "6-31G*",
                  charge: int = 0, mult: int = 1, **kwargs) -> Calculator:
    from ase.calculators.qchem import QChem
    calc = QChem(
        method=method,
        basis=basis,
        charge=charge,
        multiplicity=mult,
        **kwargs,
    )
    logger.info("Created Q-Chem calculator (%s/%s)", method, basis)
    return calc


def _create_cp2k(xc: str = "PBE", basis: str = "DZVP-MOLOPT-SR-GTH",
                 **kwargs) -> Calculator:
    from ase.calculators.cp2k import CP2K
    calc = CP2K(
        xc=xc,
        basis_set=basis,
        **kwargs,
    )
    logger.info("Created CP2K calculator (%s/%s)", xc, basis)
    return calc


def _create_espresso(pseudopotentials: dict | None = None,
                     ecutwfc: float = 40.0, kpts: tuple = (1, 1, 1),
                     xc: str = "PBE", **kwargs) -> Calculator:
    from ase.calculators.espresso import Espresso
    input_data = {"system": {"ecutwfc": ecutwfc}}
    calc = Espresso(
        pseudopotentials=pseudopotentials or {},
        input_data=input_data,
        kpts=kpts,
        **kwargs,
    )
    logger.info("Created Quantum ESPRESSO calculator (ecutwfc=%.1f)", ecutwfc)
    return calc


def _create_aims(xc: str = "PBE", **kwargs) -> Calculator:
    from ase.calculators.aims import Aims
    calc = Aims(xc=xc, **kwargs)
    logger.info("Created FHI-aims calculator (%s)", xc)
    return calc


def _create_siesta(xc: str = "PBE", **kwargs) -> Calculator:
    from ase.calculators.siesta import Siesta
    calc = Siesta(xc=xc, **kwargs)
    logger.info("Created SIESTA calculator (%s)", xc)
    return calc


def _create_gpaw(xc: str = "PBE", mode: str = "fd", **kwargs) -> Calculator:
    from gpaw import GPAW
    calc = GPAW(xc=xc, mode=mode, **kwargs)
    logger.info("Created GPAW calculator (%s, mode=%s)", xc, mode)
    return calc


def _create_vasp(xc: str = "PBE", kpts: tuple = (1, 1, 1),
                 encut: float = 400, **kwargs) -> Calculator:
    from ase.calculators.vasp import Vasp
    calc = Vasp(xc=xc, kpts=kpts, encut=encut, **kwargs)
    logger.info("Created VASP calculator (%s, encut=%.0f)", xc, encut)
    return calc


def _create_gamess(**kwargs) -> Calculator:
    from ase.calculators.gamess_us import GAMESSUS
    calc = GAMESSUS(**kwargs)
    logger.info("Created GAMESS calculator")
    return calc


def _create_dalton(**kwargs) -> Calculator:
    from ase.calculators.dalton import Dalton
    calc = Dalton(**kwargs)
    logger.info("Created Dalton calculator")
    return calc


def _create_abinit(**kwargs) -> Calculator:
    from ase.calculators.abinit import Abinit
    calc = Abinit(**kwargs)
    logger.info("Created ABINIT calculator")
    return calc


def _create_crystal(**kwargs) -> Calculator:
    from ase.calculators.crystal import CRYSTAL
    calc = CRYSTAL(**kwargs)
    logger.info("Created CRYSTAL calculator")
    return calc


def _create_fleur(**kwargs) -> Calculator:
    from ase.calculators.fleur import FLEUR
    calc = FLEUR(**kwargs)
    logger.info("Created FLEUR calculator")
    return calc


def _create_openmolcas(**kwargs) -> Calculator:
    from ase.calculators.openmx import OpenMX
    # OpenMolcas doesn't have a native ASE calculator;
    # use generic command-line calculator pattern
    from ase.calculators.genericfileio import GenericFileIOCalculator
    raise NotImplementedError(
        "OpenMolcas ASE calculator is not yet available. "
        "Use pymolcas CLI wrapper instead."
    )


def _create_molpro(**kwargs) -> Calculator:
    # Molpro has no built-in ASE calculator; community implementations exist
    raise NotImplementedError(
        "Molpro ASE calculator requires the molpro-python interface. "
        "Install via: pip install pymolpro"
    )


def _create_pyscf(method: str = "b3lyp", basis: str = "6-31G*",
                  charge: int = 0, mult: int = 1, **kwargs) -> Calculator:
    from ase.calculators.pyscf import PySCF
    calc = PySCF(
        method=method,
        basis=basis,
        charge=charge,
        spin=max(0, mult - 1),
        **kwargs,
    )
    logger.info("Created PySCF calculator (%s/%s)", method, basis)
    return calc


# ── Semi-empirical ───────────────────────────────────────────────────

def _create_xtb(method: str = "GFN2-xTB", charge: int = 0,
                mult: int = 1, **kwargs) -> Calculator:
    from ase.calculators.orca import ORCA
    # xTB through ORCA interface (most reliable)
    orca_path = shutil.which("orca")
    if orca_path:
        calc = ORCA(
            orcasimpleinput=f"{method}",
            charge=charge,
            mult=mult,
            **kwargs,
        )
        logger.info("Created xTB calculator via ORCA (%s)", method)
        return calc
    # Fallback: standalone xTB
    from ase.calculators.xtb import XTB  # type: ignore[import-not-found]
    calc = XTB(method=method, charge=charge, uhf=max(0, mult - 1), **kwargs)
    logger.info("Created standalone xTB calculator (%s)", method)
    return calc


def _create_dftb(**kwargs) -> Calculator:
    from ase.calculators.dftb import Dftb
    calc = Dftb(**kwargs)
    logger.info("Created DFTB+ calculator")
    return calc


def _create_mopac(method: str = "PM7", charge: int = 0,
                  mult: int = 1, **kwargs) -> Calculator:
    from ase.calculators.mopac import MOPAC
    task = method
    if mult > 1:
        task += f" UHFSINGLET" if mult == 1 else f" UHF"
    calc = MOPAC(method=task, charge=charge, **kwargs)
    logger.info("Created MOPAC calculator (%s)", method)
    return calc


# ── MD engines ───────────────────────────────────────────────────────

def _create_lammps(**kwargs) -> Calculator:
    from ase.calculators.lammpsrun import LAMMPS
    calc = LAMMPS(**kwargs)
    logger.info("Created LAMMPS calculator")
    return calc


def _create_openmm(forcefield: str = "amber14-all.xml", **kwargs) -> Calculator:
    try:
        from openmm.app import ForceField
        from openmmml import MLPotential
    except ImportError:
        raise ImportError("OpenMM not installed. Install via: conda install -c conda-forge openmm")
    raise NotImplementedError(
        "OpenMM ASE calculator requires setup via openmmml. "
        "See: https://github.com/openmm/openmm-ml"
    )


def _create_gromacs(**kwargs) -> Calculator:
    from ase.calculators.gromacs import Gromacs
    calc = Gromacs(**kwargs)
    logger.info("Created GROMACS calculator")
    return calc


# ═══════════════════════════════════════════════════════════════════════
#  Factory dispatch table
# ═══════════════════════════════════════════════════════════════════════

_FACTORIES = {
    # MLP — delegated to mlp_tools
    "ani2x":      lambda **kw: _create_mlp("ani2x", **kw),
    "aimnet2":    lambda **kw: _create_mlp("aimnet2", **kw),
    "mace_off":   lambda **kw: _create_mlp("mace_off", **kw),
    "chgnet":     lambda **kw: _create_mlp("chgnet", **kw),
    "m3gnet":     lambda **kw: _create_mlp("m3gnet", **kw),
    "schnetpack": lambda **kw: _create_mlp("schnetpack", **kw),
    "nequip":     lambda **kw: _create_mlp("nequip", **kw),
    "alignn":     lambda **kw: _create_mlp("alignn", **kw),
    # QM
    "orca":       _create_orca,
    "gaussian":   _create_gaussian,
    "turbomole":  _create_turbomole,
    "nwchem":     _create_nwchem,
    "psi4":       _create_psi4,
    "qchem":      _create_qchem,
    "cp2k":       _create_cp2k,
    "espresso":   _create_espresso,
    "aims":       _create_aims,
    "siesta":     _create_siesta,
    "gpaw":       _create_gpaw,
    "vasp":       _create_vasp,
    "gamess":     _create_gamess,
    "dalton":     _create_dalton,
    "abinit":     _create_abinit,
    "crystal":    _create_crystal,
    "fleur":      _create_fleur,
    "openmolcas": _create_openmolcas,
    "molpro":     _create_molpro,
    "pyscf":      _create_pyscf,
    # Semi-empirical
    "xtb":        _create_xtb,
    "dftb":       _create_dftb,
    "mopac":      _create_mopac,
    # MD
    "lammps":     _create_lammps,
    "openmm":     _create_openmm,
    "gromacs":    _create_gromacs,
}


# ═══════════════════════════════════════════════════════════════════════
#  Public API
# ═══════════════════════════════════════════════════════════════════════

def create_calculator(backend: str, **kwargs) -> Calculator:
    """Create an ASE calculator for any supported backend.

    Parameters
    ----------
    backend : str
        Backend name (e.g. "orca", "ani2x", "gaussian", "vasp", "xtb").
        Case-insensitive, aliases accepted (e.g. "g16" → "gaussian").
    **kwargs
        Backend-specific parameters (method, basis, charge, mult, device, ...).

    Returns
    -------
    ASE Calculator ready for ``atoms.calc = calc``.

    Examples
    --------
    >>> calc = create_calculator("orca", method="B3LYP", basis="def2-SVP")
    >>> calc = create_calculator("ani2x", device="cuda")
    >>> calc = create_calculator("vasp", xc="PBE", kpts=[4,4,4], encut=500)
    """
    entry = _REGISTRY.get(backend.lower().strip())
    if entry is None:
        all_names = sorted(set(k for k, (_, cat) in _REGISTRY.items() if k == _REGISTRY[k][0]))
        raise ValueError(
            f"Unknown backend '{backend}'.\n"
            f"Available: {', '.join(all_names)}"
        )
    key, category = entry
    factory = _FACTORIES.get(key)
    if factory is None:
        raise NotImplementedError(f"Backend '{key}' is registered but has no factory yet.")
    return factory(**kwargs)


def available_backends(category: Optional[str] = None) -> list[str]:
    """Return canonical names of all registered backends.

    Parameters
    ----------
    category : optional filter — "mlp", "qm", "semiempirical", "md", or None for all.
    """
    seen = set()
    result = []
    for alias, (canonical, cat) in _REGISTRY.items():
        if canonical in seen:
            continue
        if category and cat != category:
            continue
        seen.add(canonical)
        result.append(canonical)
    return sorted(result)


def backend_info(backend: str) -> dict:
    """Return metadata about a backend."""
    entry = _REGISTRY.get(backend.lower().strip())
    if entry is None:
        return {"name": backend, "registered": False}
    key, category = entry
    aliases = [a for a, (k, _) in _REGISTRY.items() if k == key and a != key]
    return {
        "name": key,
        "category": category,
        "aliases": aliases,
        "has_factory": key in _FACTORIES,
        "registered": True,
    }
