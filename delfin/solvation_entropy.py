"""Solution-phase entropy corrections for routine DELFIN workflows.

This module implements the solution-entropy model described by Ariai and
Gellrich in a clean-room Python implementation for DELFIN:

    J. Ariai and U. Gellrich, Phys. Chem. Chem. Phys. 2023, 25, 14005.
    DOI: 10.1039/d3cp00970j

The underlying solution-phase entropy formalism follows Garza:

    A. J. Garza, J. Chem. Theory Comput. 2019, 15, 3204.

The legacy shell scripts by Jama Ariai in the development workspace were used
only as scientific and behavioral reference material. Their code is not copied
here. DELFIN's implementation is licensed under DELFIN's LGPL-3.0-or-later
license. vdW radii are a compact Bondi/Rowland-Truhlar style table for this
calculation path; see the README method note for attribution.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


R_GAS_CAL = 8.314462618 / 4.184
AU_TO_KCAL = 627.5095
CM_TO_K = 1.43877696


@dataclass(frozen=True)
class SolventParams:
    """Physical solvent parameters used by the solution-entropy model."""

    name: str
    permittivity: float
    molar_mass_g_mol: float
    density_g_ml: float
    vdw_volume_A3: float
    source: str = ""


@dataclass(frozen=True)
class MoleculeGeometry:
    """Simple XYZ geometry container."""

    symbols: tuple[str, ...]
    coords_A: tuple[tuple[float, float, float], ...]

    @property
    def n_atoms(self) -> int:
        return len(self.symbols)


@dataclass(frozen=True)
class Thermochemistry:
    """Thermochemistry parsed from a QC output."""

    program: str
    temperature_K: float | None = None
    rotational_entropy_cal_mol_K: float | None = None
    frequencies_cm1: tuple[float, ...] = ()
    point_group: str | None = None
    symmetry_number: int | None = None


@dataclass(frozen=True)
class EntropyComponents:
    """Solution entropy components in cal mol-1 K-1 unless stated otherwise."""

    solvent: str
    temperature_K: float
    concentration_standard: str
    vdw_volume_A3: float
    solvent_vdw_volume_A3: float
    molecular_weight_g_mol: float
    radius_gyration_A: float
    cavity_volume_A3: float
    accessible_cavities: float
    S_trans: float
    S_trans_diff: float
    S_rot_diff: float
    S_cav: float
    S_conc: float
    S_solv: float
    S_rot: float | None = None
    S_vib_qrrho: float | None = None
    S_soln: float | None = None
    S_rot_gas: float | None = None
    n_imaginary_ignored: int = 0
    imaginary_frequencies_cm1: tuple[float, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class ReactionEntropyCorrection:
    """Stoichiometric entropy correction for a reaction or barrier."""

    delta_S_solv_cal_mol_K: float
    delta_G_entropy_corr_kcal_mol: float
    delta_G_corrected_kcal_mol: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


VDW_RADII_A = {
    "H": 1.10, "He": 1.40,
    "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55,
    "O": 1.52, "F": 1.47, "Ne": 1.54,
    "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80,
    "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Sc": 2.30, "Ti": 2.15, "V": 2.05,
    "Cr": 2.05, "Mn": 2.05, "Fe": 2.00, "Co": 2.00, "Ni": 1.97,
    "Cu": 1.96, "Zn": 2.01, "Ga": 1.87, "Ge": 2.11, "As": 1.85,
    "Se": 1.90, "Br": 1.85, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49, "Y": 2.40, "Zr": 2.30, "Nb": 2.15,
    "Mo": 2.10, "Tc": 2.05, "Ru": 2.05, "Rh": 2.00, "Pd": 2.05,
    "Ag": 2.10, "Cd": 2.18, "In": 1.93, "Sn": 2.17, "Sb": 2.06,
    "Te": 2.06, "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "La": 2.50, "Hf": 2.25, "Ta": 2.20,
    "W": 2.10, "Re": 2.05, "Os": 2.00, "Ir": 2.00, "Pt": 2.05,
    "Au": 2.10, "Hg": 2.05, "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
    "Po": 1.97, "At": 2.02, "Rn": 2.20,
    "Fr": 3.48, "Ra": 2.83,
}

ATOMIC_WEIGHTS = {
    "H": 1.008, "He": 4.002602, "Li": 6.94, "Be": 9.0121831,
    "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998403163,
    "Ne": 20.1797, "Na": 22.98976928, "Mg": 24.305, "Al": 26.9815385,
    "Si": 28.085, "P": 30.973761998, "S": 32.06, "Cl": 35.45,
    "Ar": 39.948, "K": 39.0983, "Ca": 40.078, "Sc": 44.955908,
    "Ti": 47.867, "V": 50.9415, "Cr": 51.9961, "Mn": 54.938044,
    "Fe": 55.845,
    "Co": 58.933194, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.921595, "Se": 78.971,
    "Br": 79.904, "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62,
    "Y": 88.90584, "Zr": 91.224, "Nb": 92.90637, "Mo": 95.95,
    "Tc": 98.0, "Ru": 101.07, "Rh": 102.90550, "Pd": 106.42,
    "Ag": 107.8682, "Cd": 112.414, "In": 114.818, "Sn": 118.710,
    "Sb": 121.760, "Te": 127.60, "I": 126.90447, "Xe": 131.293,
    "Cs": 132.90545196, "Ba": 137.327, "La": 138.90547,
    "Hf": 178.49, "Ta": 180.94788, "W": 183.84, "Re": 186.207,
    "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
    "Au": 196.966569, "Hg": 200.592, "Tl": 204.38, "Pb": 207.2,
    "Bi": 208.98040, "Po": 209.0, "At": 210.0, "Rn": 222.0,
}

ATOMIC_NUMBER_SYMBOLS = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)

_SOLVENTS = {
    "water": SolventParams("water", 78.3553, 18.015, 0.9970, 20.69, "curated from Ariai/Gellrich-style solvent table"),
    "acetonitrile": SolventParams("acetonitrile", 36.64, 41.052, 0.7857, 43.40, "curated from Ariai/Gellrich-style solvent table"),
    "dmso": SolventParams("dmso", 46.826, 78.133, 1.0955, 71.30, "curated from Ariai/Gellrich-style solvent table"),
    "dmf": SolventParams("dmf", 37.219, 73.095, 0.9445, 76.55, "curated from Ariai/Gellrich-style solvent table"),
    "dimethylformamide": SolventParams("dmf", 37.219, 73.095, 0.9445, 76.55, "curated from Ariai/Gellrich-style solvent table"),
    "methanol": SolventParams("methanol", 32.613, 32.042, 0.7914, 35.32, "curated from Ariai/Gellrich-style solvent table"),
    "ethanol": SolventParams("ethanol", 24.852, 46.069, 0.7893, 52.22, "curated from Ariai/Gellrich-style solvent table"),
    "thf": SolventParams("thf", 7.4257, 72.106, 0.8892, 73.54, "curated from Ariai/Gellrich-style solvent table"),
    "tetrahydrofuran": SolventParams("thf", 7.4257, 72.106, 0.8892, 73.54, "curated from Ariai/Gellrich-style solvent table"),
    "dichloromethane": SolventParams("dichloromethane", 8.93, 84.933, 1.3266, 61.28, "curated from Ariai/Gellrich-style solvent table"),
    "dcm": SolventParams("dichloromethane", 8.93, 84.933, 1.3266, 61.28, "curated from Ariai/Gellrich-style solvent table"),
    "chloroform": SolventParams("chloroform", 4.8069, 119.378, 1.4788, 69.55, "curated from Ariai/Gellrich-style solvent table"),
    "chcl3": SolventParams("chloroform", 4.8069, 119.378, 1.4788, 69.55, "curated from Ariai/Gellrich-style solvent table"),
    "toluene": SolventParams("toluene", 2.3741, 92.139, 0.8669, 96.31, "curated from Ariai/Gellrich-style solvent table"),
    "benzene": SolventParams("benzene", 2.2825, 78.112, 0.8765, 78.86, "curated from Ariai/Gellrich-style solvent table"),
}


def solvent_params(
    solvent: str,
    *,
    density_g_ml: float | None = None,
    molar_mass_g_mol: float | None = None,
    solvent_vdw_volume_A3: float | None = None,
    permittivity: float | None = None,
) -> SolventParams:
    """Resolve solvent parameters from DELFIN's curated table or overrides."""

    key = str(solvent or "").strip().lower()
    base = _SOLVENTS.get(key)
    if base is None:
        if None in (density_g_ml, molar_mass_g_mol, solvent_vdw_volume_A3, permittivity):
            raise ValueError(
                f"unknown solvent {solvent!r}; provide density_g_ml, "
                "molar_mass_g_mol, solvent_vdw_volume_A3 and permittivity"
            )
        return SolventParams(
            key or "custom",
            float(permittivity),
            float(molar_mass_g_mol),
            float(density_g_ml),
            float(solvent_vdw_volume_A3),
            "user-provided",
        )
    return SolventParams(
        base.name,
        float(permittivity if permittivity is not None else base.permittivity),
        float(molar_mass_g_mol if molar_mass_g_mol is not None else base.molar_mass_g_mol),
        float(density_g_ml if density_g_ml is not None else base.density_g_ml),
        float(solvent_vdw_volume_A3 if solvent_vdw_volume_A3 is not None else base.vdw_volume_A3),
        base.source,
    )


def read_xyz(path: str | Path) -> MoleculeGeometry:
    """Read a simple XYZ file."""

    lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
    if not lines:
        raise ValueError(f"empty XYZ file: {path}")
    try:
        n_atoms = int(lines[0].strip())
        coord_lines = lines[2:2 + n_atoms]
    except ValueError:
        coord_lines = lines
        n_atoms = len(coord_lines)
    symbols: list[str] = []
    coords: list[tuple[float, float, float]] = []
    for line in coord_lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        symbols.append(_normalise_symbol(parts[0]))
        coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
    if len(symbols) != n_atoms:
        raise ValueError(f"invalid XYZ coordinate count in {path}")
    return MoleculeGeometry(tuple(symbols), tuple(coords))


def write_report(path: str | Path, components: EntropyComponents) -> Path:
    """Write a JSON report for downstream DELFIN tools."""

    p = Path(path)
    p.write_text(json.dumps(components.to_dict(), indent=2, sort_keys=True), encoding="utf-8")
    return p


def read_report(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def parse_orca_thermochemistry(path: str | Path) -> Thermochemistry:
    """Parse the final ORCA thermochemistry section used by the entropy model."""

    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    marker = re.search(r"THERMOCHEMISTRY", text, flags=re.IGNORECASE)
    section = text[marker.start():] if marker else text
    temp = _last_float(section, r"Temperature\s+\.\.\.\s+([-+]?\d+(?:\.\d+)?)")
    if temp is None:
        temp = _last_float(section, r"Temperature\s+([-+]?\d+(?:\.\d+)?)")
    s_rot_au = _last_float(section, r"Rotational entropy\s+([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)")
    s_rot = None
    if s_rot_au is not None and temp:
        s_rot = s_rot_au * AU_TO_KCAL * 1000.0 / temp
    pg = _last_group(section, r"Point Group\s+([A-Za-z0-9_]+)")
    sym = _last_float(section, r"Symmetry Number\s+([0-9]+)")
    freqs = tuple(float(x) for x in re.findall(r"\bfreq\.\s+([-+]?\d+(?:\.\d+)?)", section))
    return Thermochemistry(
        program="orca",
        temperature_K=temp,
        rotational_entropy_cal_mol_K=s_rot,
        frequencies_cm1=freqs,
        point_group=pg,
        symmetry_number=int(sym) if sym is not None else None,
    )


def parse_gaussian_thermochemistry(path: str | Path) -> Thermochemistry:
    """Parse Gaussian thermochemistry values needed by the entropy model."""

    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    temp = _last_float(text, r"Temperature\s+([-+]?\d+(?:\.\d+)?)")
    s_rot = _last_float(text, r"Rotational\s+([-+]?\d+(?:\.\d+)?)")
    freqs: list[float] = []
    for line in re.findall(r"Frequencies\s+--\s+([^\n]+)", text):
        freqs.extend(float(x) for x in line.split())
    return Thermochemistry(
        program="gaussian",
        temperature_K=temp,
        rotational_entropy_cal_mol_K=s_rot,
        frequencies_cm1=tuple(freqs),
    )


def parse_thermochemistry(path: str | Path) -> Thermochemistry:
    """Detect ORCA/Gaussian and parse thermochemistry."""

    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    if "* O   R   C   A *" in text or "O   R   C   A" in text:
        return parse_orca_thermochemistry(path)
    if "Entering Gaussian System" in text:
        return parse_gaussian_thermochemistry(path)
    raise ValueError(f"could not identify QC output format: {path}")


def calculate_solution_entropy(
    geometry: MoleculeGeometry | str | Path,
    *,
    solvent: str = "benzene",
    temperature_K: float = 298.15,
    thermochemistry: Thermochemistry | None = None,
    concentration_standard: str = "1M",
    symmetry_number: int | None = None,
    density_g_ml: float | None = None,
    molar_mass_g_mol: float | None = None,
    solvent_vdw_volume_A3: float | None = None,
    permittivity: float | None = None,
    volume_samples: int = 60_000,
) -> EntropyComponents:
    """Calculate solution entropy components for one species."""

    geom = read_xyz(geometry) if isinstance(geometry, (str, Path)) else geometry
    solv = solvent_params(
        solvent,
        density_g_ml=density_g_ml,
        molar_mass_g_mol=molar_mass_g_mol,
        solvent_vdw_volume_A3=solvent_vdw_volume_A3,
        permittivity=permittivity,
    )
    mass = molecular_weight(geom)
    vol_m = vdw_volume(geom, samples=volume_samples)
    r_gyr = radius_of_gyration(geom)
    vol_free = (solv.molar_mass_g_mol * 1.0e24) / (6.02214076e23 * solv.density_g_ml) - solv.vdw_volume_A3
    vol_cav = (vol_m ** (1.0 / 3.0) + vol_free ** (1.0 / 3.0)) ** 3
    if vol_free > vol_m:
        x = (
            (vol_free ** (2.0 / 3.0) - vol_m ** (2.0 / 3.0))
            / (vol_free ** (2.0 / 3.0) + solv.vdw_volume_A3 ** (2.0 / 3.0))
        )
    else:
        x = 0.0
    n_cav = 1.0 + (
        (4.0 * vol_cav ** (2.0 / 3.0))
        / (vol_free ** (2.0 / 3.0) + solv.vdw_volume_A3 ** (2.0 / 3.0))
    ) * ((1.0 / (1.0 - x)) - 1.0)
    vol = n_cav * vol_cav

    s_trans_diff = R_GAS_CAL * math.log(vol / temperature_K) - 9.792392
    s_trans = R_GAS_CAL * ((1.5 * math.log(mass * temperature_K)) + math.log(vol) - 6.079431401)

    r_cav = ((3.0 * vol_cav) / (4.0 * math.pi)) ** (1.0 / 3.0)
    if r_cav > r_gyr:
        s_rot_diff = 3.0 * R_GAS_CAL * math.log((r_cav - r_gyr) / r_cav)
    else:
        r_free = ((3.0 * vol_free) / (4.0 * math.pi)) ** (1.0 / 3.0)
        theta0_arg = r_gyr / math.sqrt(r_gyr**2 + r_free**2)
        theta0 = 2.0 * math.acos(max(-1.0, min(1.0, theta0_arg)))
        s_rot_diff = (
            R_GAS_CAL * (2.0 * math.log(theta0 / math.pi) + temperature_K)
            + 3.0 * R_GAS_CAL * math.log(r_free / r_cav)
        )

    y = (3.0 / (4.0 * math.pi)) * ((solv.permittivity - 1.0) / (solv.permittivity + 2.0))
    ratio = (vol_m / solv.vdw_volume_A3) ** (1.0 / 3.0)
    s_cav = -R_GAS_CAL * (
        -math.log(1.0 - y)
        + (3.0 * y * ratio) / (1.0 - y)
        + (((3.0 * y) / (1.0 - y)) + 4.5 * ((y / (1.0 - y)) ** 2)) * (ratio**2)
    )

    standard = concentration_standard.strip().lower()
    s_conc = -R_GAS_CAL * math.log(temperature_K) + 4.942522175
    if standard in {"liquid", "liq", "pure_liquid"}:
        s_conc += R_GAS_CAL * math.log(solv.molar_mass_g_mol / (1000.0 * solv.density_g_ml))
        ref = "liq"
    else:
        ref = "1M"

    s_solv = s_trans_diff + s_cav + s_rot_diff + s_conc

    s_rot_gas = thermochemistry.rotational_entropy_cal_mol_K if thermochemistry else None
    if s_rot_gas is not None and symmetry_number and symmetry_number > 0:
        old = thermochemistry.symmetry_number if thermochemistry else None
        old = old if old and old > 0 else 1
        s_rot_gas += -R_GAS_CAL * math.log(symmetry_number) + R_GAS_CAL * math.log(old)
    s_rot = s_rot_gas + s_rot_diff if s_rot_gas is not None else None
    s_vib, imaginary = _qrrho_entropy(
        thermochemistry.frequencies_cm1 if thermochemistry else (),
        temperature_K,
    )
    if not thermochemistry or not thermochemistry.frequencies_cm1:
        s_vib = None
    s_soln = None
    if s_rot is not None and s_vib is not None:
        s_soln = s_trans + s_rot + s_vib + s_cav + s_conc

    return EntropyComponents(
        solvent=solv.name,
        temperature_K=float(temperature_K),
        concentration_standard=ref,
        vdw_volume_A3=float(vol_m),
        solvent_vdw_volume_A3=float(solv.vdw_volume_A3),
        molecular_weight_g_mol=float(mass),
        radius_gyration_A=float(r_gyr),
        cavity_volume_A3=float(vol_cav),
        accessible_cavities=float(n_cav),
        S_trans=float(s_trans),
        S_trans_diff=float(s_trans_diff),
        S_rot_diff=float(s_rot_diff),
        S_cav=float(s_cav),
        S_conc=float(s_conc),
        S_solv=float(s_solv),
        S_rot=None if s_rot is None else float(s_rot),
        S_vib_qrrho=None if s_vib is None else float(s_vib),
        S_soln=None if s_soln is None else float(s_soln),
        S_rot_gas=None if s_rot_gas is None else float(s_rot_gas),
        n_imaginary_ignored=len(imaginary),
        imaginary_frequencies_cm1=tuple(imaginary),
    )


def calculate_reaction_entropy_correction(
    species: Mapping[str, EntropyComponents | Mapping[str, Any] | float],
    stoichiometry: Mapping[str, float],
    *,
    temperature_K: float = 298.15,
    uncorrected_delta_g_kcal_mol: float | None = None,
) -> ReactionEntropyCorrection:
    """Calculate a stoichiometric solution-entropy correction.

    Stoichiometry convention: products positive, reactants negative.
    """

    delta_s = 0.0
    for name, coeff in stoichiometry.items():
        value = species[name]
        if isinstance(value, EntropyComponents):
            s_solv = value.S_solv
        elif isinstance(value, Mapping):
            s_solv = float(value["S_solv"])
        else:
            s_solv = float(value)
        delta_s += float(coeff) * s_solv
    delta_g = -float(temperature_K) * delta_s / 1000.0
    corrected = (
        None if uncorrected_delta_g_kcal_mol is None
        else float(uncorrected_delta_g_kcal_mol) + delta_g
    )
    return ReactionEntropyCorrection(delta_s, delta_g, corrected)


def molecular_weight(geometry: MoleculeGeometry) -> float:
    total = 0.0
    for sym in geometry.symbols:
        if sym not in ATOMIC_WEIGHTS:
            raise ValueError(f"missing atomic weight for element {sym}")
        total += ATOMIC_WEIGHTS[sym]
    return total


def radius_of_gyration(geometry: MoleculeGeometry) -> float:
    coords = np.asarray(geometry.coords_A, dtype=float)
    masses = np.asarray([ATOMIC_WEIGHTS[s] for s in geometry.symbols], dtype=float)
    center = np.average(coords, axis=0, weights=masses)
    return float(np.sqrt(np.mean(np.sum((coords - center) ** 2, axis=1))))


def vdw_volume(geometry: MoleculeGeometry, *, samples: int = 60_000, seed: int = 112358) -> float:
    """Estimate union-of-vdW-spheres volume in A^3.

    Single atoms use the exact sphere volume. Multi-atom systems use a
    deterministic Monte-Carlo integration over the bounding box. The stochastic
    integration is intentionally isolated here so accuracy/performance can be
    upgraded without changing the public entropy API.
    """

    coords = np.asarray(geometry.coords_A, dtype=float)
    radii = np.asarray([VDW_RADII_A.get(s, 2.0) for s in geometry.symbols], dtype=float)
    if len(radii) == 1:
        return float((4.0 / 3.0) * math.pi * radii[0] ** 3)
    lower = np.min(coords - radii[:, None], axis=0)
    upper = np.max(coords + radii[:, None], axis=0)
    extent = upper - lower
    rng = np.random.default_rng(seed)
    pts = lower + rng.random((int(samples), 3)) * extent
    inside = np.zeros(len(pts), dtype=bool)
    for xyz, r in zip(coords, radii):
        inside |= np.sum((pts - xyz) ** 2, axis=1) <= r * r
    return float(np.prod(extent) * np.mean(inside))


def _qrrho_entropy(frequencies_cm1: Sequence[float], temperature_K: float) -> tuple[float, list[float]]:
    s_vib = 0.0
    imaginary: list[float] = []
    mol_inertia_avg = 1.0
    freq_0 = 100.0
    damp = 4.0
    for freq in frequencies_cm1:
        f = float(freq)
        if f < 0:
            imaginary.append(f)
            continue
        if f == 0:
            continue
        theta_vib = CM_TO_K * f
        coeff = theta_vib / temperature_K
        s_pure = R_GAS_CAL * ((coeff / math.expm1(coeff)) - math.log1p(-math.exp(-coeff)))
        s_free = 0.5 * R_GAS_CAL * (
            5.356746219 + math.log(temperature_K) - math.log(mol_inertia_avg + (35.72353633 * f))
        )
        weight = 1.0 / (1.0 + ((freq_0 / f) ** damp))
        s_vib += weight * s_pure + (1.0 - weight) * s_free
    return float(s_vib), imaginary


def _normalise_symbol(value: str) -> str:
    raw = str(value).strip()
    if raw.isdigit():
        idx = int(raw)
        if not 1 <= idx <= len(ATOMIC_NUMBER_SYMBOLS):
            raise ValueError(f"invalid atomic number: {raw}")
        return ATOMIC_NUMBER_SYMBOLS[idx - 1]
    return raw[:1].upper() + raw[1:].lower()


def _last_float(text: str, pattern: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.IGNORECASE)
    if not matches:
        return None
    value = matches[-1]
    if isinstance(value, tuple):
        value = value[-1]
    return float(value)


def _last_group(text: str, pattern: str) -> str | None:
    matches = re.findall(pattern, text, flags=re.IGNORECASE)
    return str(matches[-1]) if matches else None


__all__ = [
    "EntropyComponents",
    "MoleculeGeometry",
    "ReactionEntropyCorrection",
    "SolventParams",
    "Thermochemistry",
    "calculate_reaction_entropy_correction",
    "calculate_solution_entropy",
    "parse_gaussian_thermochemistry",
    "parse_orca_thermochemistry",
    "parse_thermochemistry",
    "read_report",
    "read_xyz",
    "solvent_params",
    "vdw_volume",
    "write_report",
]
