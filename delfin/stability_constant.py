"""Stability constant (log K) calculation via Born-Haber thermodynamic cycle.

Thermodynamic cycle:
  [M(Solv)_m]^x+  +  Σ L_i  →  [ML_1...L_k(Solv)_{m-n}]^y+  +  n·Solv

  ΔG = G(complex) + n·G(solvent) - G(solv_complex) - Σ(count_i · G(ligand_i))
  log K = -ΔG / (2.303 · R · T)
"""

from __future__ import annotations

import copy
import math
import os
import re
import shutil
import threading
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

from delfin.common.logging import get_logger
from delfin.common.orca_blocks import OrcaInputBuilder, collect_output_blocks, resolve_maxiter

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

HARTREE_TO_KCAL = 627.5095  # 1 Eh = 627.5095 kcal/mol
HARTREE_TO_KJ = 2625.5     # 1 Eh = 2625.5 kJ/mol
R_KCAL = 0.001987204       # gas constant in kcal/(mol·K)

# ---------------------------------------------------------------------------
# Solvent database: name → {smiles, donor atom, display name}
# ---------------------------------------------------------------------------

SOLVENT_DB: Dict[str, Dict[str, str]] = {
    "water":            {"smiles": "O",           "donor": "O",  "name": "Water"},
    "dmso":             {"smiles": "CS(=O)C",     "donor": "O",  "name": "DMSO"},
    "dmf":              {"smiles": "O=CN(C)C",    "donor": "O",  "name": "DMF"},
    "acetonitrile":     {"smiles": "CC#N",        "donor": "N",  "name": "MeCN"},
    "thf":              {"smiles": "C1CCOC1",     "donor": "O",  "name": "THF"},
    "methanol":         {"smiles": "CO",          "donor": "O",  "name": "MeOH"},
    "ethanol":          {"smiles": "CCO",         "donor": "O",  "name": "EtOH"},
    "pyridine":         {"smiles": "c1ccncc1",    "donor": "N",  "name": "Pyridine"},
    "dichloromethane":  {"smiles": "ClCCl",       "donor": "Cl", "name": "DCM"},
    "toluene":          {"smiles": "Cc1ccccc1",   "donor": "C",  "name": "Toluene"},
    "acetone":          {"smiles": "CC(=O)C",     "donor": "O",  "name": "Acetone"},
    "chloroform":       {"smiles": "ClC(Cl)Cl",   "donor": "Cl", "name": "CHCl3"},
    "diethylether":     {"smiles": "CCOCC",       "donor": "O",  "name": "Et2O"},
    "1,4-dioxane":      {"smiles": "C1COCCO1",    "donor": "O",  "name": "Dioxane"},
    "nitromethane":     {"smiles": "C[N+](=O)[O-]", "donor": "O", "name": "MeNO2"},
    "formamide":        {"smiles": "O=CN",        "donor": "O",  "name": "Formamide"},
    "nmp":              {"smiles": "CN1CCCC1=O",  "donor": "O",  "name": "NMP"},
    "dmac":             {"smiles": "CN(C)C(=O)C", "donor": "O",  "name": "DMAc"},
}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class LigandInfo:
    """Information about a single ligand extracted from the complex SMILES."""
    smiles: str
    charge: int
    denticity: int
    label: str  # e.g. "bipy", "acac"


@dataclass
class StabilityAnalysis:
    """Complete analysis of the stability constant reaction."""
    complex_smiles: str
    metals: List[Tuple[str, int]]           # [(symbol, charge), ...]
    ligands: List[LigandInfo]               # all ligands (may have duplicates)
    unique_ligands: List[LigandInfo]         # unique ligands (for computation)
    ligand_counts: Dict[str, int]           # smiles → count
    solvent_name: str
    solvent_smiles: str
    solvent_donor: str
    n_explicit_solvent: int
    n_displaced: int
    solv_complex_smiles: str
    metal_charge: int                       # total metal charge


@dataclass
class SpeciesEnergy:
    """Energy data extracted from an ORCA output file."""
    label: str
    folder: str
    charge: int
    multiplicity: int
    e_el: Optional[float] = None       # Hartree
    g_rrho: Optional[float] = None     # Hartree (thermal free energy correction)
    g_total: Optional[float] = None    # Hartree (E_el + G_RRHO)
    occupier_mult: Optional[int] = None
    occupier_bs: Optional[str] = None
    occupier_idx: Optional[int] = None


@dataclass
class StabilityResult:
    """Final stability constant result."""
    delta_g_hartree: float
    delta_g_kcal: float
    delta_g_kj: float
    log_k: float
    temperature: float
    species: List[SpeciesEnergy]
    analysis: StabilityAnalysis
    logK_exp: Optional[float] = None


@dataclass
class StabilityWorkflowPlan:
    """Prepared SC workflow jobs plus the parsed reaction analysis."""
    analysis: StabilityAnalysis
    sc_dir: Path
    jobs: List["WorkflowJob"]


# ---------------------------------------------------------------------------
# 1. Analysis: extract ligands with denticity from complex SMILES
# ---------------------------------------------------------------------------

def analyze_complex(smiles: str, solvent: str, n_explicit_solvent: int) -> StabilityAnalysis:
    """Parse complex SMILES and build the full reaction analysis."""
    from delfin.build_up_complex import extract_ligands_from_complex, get_ligand_charge
    from delfin.smiles_converter import _METALS, RDKIT_AVAILABLE

    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is required for stability constant analysis")

    from rdkit import Chem

    # Extract metals and ligands
    metals, ligand_smiles_list, error = extract_ligands_from_complex(smiles)
    if error:
        raise ValueError(f"Cannot extract ligands from SMILES '{smiles}': {error}")
    if not metals:
        raise ValueError(f"No metal found in SMILES '{smiles}'")
    if not ligand_smiles_list:
        raise ValueError(f"No ligands found in SMILES '{smiles}'")

    # Determine denticity per ligand by counting coordinating atoms
    denticities = _determine_denticities(smiles, ligand_smiles_list)

    # Build ligand info list
    ligands = []
    for i, (lig_smi, dent) in enumerate(zip(ligand_smiles_list, denticities)):
        charge = get_ligand_charge(lig_smi)
        label = _make_ligand_label(lig_smi, i)
        ligands.append(LigandInfo(smiles=lig_smi, charge=charge, denticity=dent, label=label))

    # Identify unique ligands
    counts: Dict[str, int] = Counter(lig.smiles for lig in ligands)
    seen: Set[str] = set()
    unique_ligands: List[LigandInfo] = []
    for lig in ligands:
        if lig.smiles not in seen:
            seen.add(lig.smiles)
            unique_ligands.append(lig)

    # Total displaced solvent molecules = sum of denticities across ALL ligands
    n_displaced = sum(lig.denticity for lig in ligands)
    if n_displaced > n_explicit_solvent:
        logger.warning(
            "n_displaced (%d) > n_explicit_solvent (%d) — "
            "complex has more donor atoms than explicit solvent molecules. "
            "Consider increasing n_explicit_solvent.",
            n_displaced, n_explicit_solvent,
        )

    # Resolve solvent
    solvent_key = solvent.strip().lower()
    if solvent_key not in SOLVENT_DB:
        raise ValueError(
            f"Solvent '{solvent}' not found in SOLVENT_DB. "
            f"Available: {', '.join(sorted(SOLVENT_DB.keys()))}"
        )
    solv_info = SOLVENT_DB[solvent_key]

    # Total metal charge
    metal_charge = sum(chg for _, chg in metals)

    # Build solvation complex SMILES
    solv_complex_smiles = build_solvation_complex_smiles(
        metals, solv_info["smiles"], solv_info["donor"], n_explicit_solvent,
    )

    return StabilityAnalysis(
        complex_smiles=smiles,
        metals=metals,
        ligands=ligands,
        unique_ligands=unique_ligands,
        ligand_counts=dict(counts),
        solvent_name=solv_info["name"],
        solvent_smiles=solv_info["smiles"],
        solvent_donor=solv_info["donor"],
        n_explicit_solvent=n_explicit_solvent,
        n_displaced=n_displaced,
        solv_complex_smiles=solv_complex_smiles,
        metal_charge=metal_charge,
    )


def _determine_denticities(complex_smiles: str, ligand_smiles_list: List[str]) -> List[int]:
    """Count coordinating atoms per ligand fragment by examining bonds to metals."""
    from rdkit import Chem
    from delfin.smiles_converter import _METALS

    mol = Chem.MolFromSmiles(complex_smiles, sanitize=False)
    if mol is None:
        return [1] * len(ligand_smiles_list)

    # Find metal atoms
    metal_idxs = set()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in _METALS:
            metal_idxs.add(atom.GetIdx())

    # Find coordinating atoms per metal
    coord_atoms: Set[int] = set()
    for bond in mol.GetBonds():
        b, e = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if b in metal_idxs:
            coord_atoms.add(e)
        elif e in metal_idxs:
            coord_atoms.add(b)

    # Now remove metals and get fragments, tracking which coord atoms belong where
    edit = Chem.RWMol(mol)
    # Remove metal bonds
    bonds_to_remove = []
    for bond in edit.GetBonds():
        b, e = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if b in metal_idxs or e in metal_idxs:
            bonds_to_remove.append((b, e))
    for b, e in bonds_to_remove:
        edit.RemoveBond(b, e)

    # Remove metal atoms (reverse order)
    # Build mapping old_idx → new_idx
    removed = sorted(metal_idxs, reverse=True)
    idx_map = {}
    shift = 0
    for old_idx in range(mol.GetNumAtoms()):
        if old_idx in metal_idxs:
            shift += 1
        else:
            idx_map[old_idx] = old_idx - shift

    for m_idx in removed:
        edit.RemoveAtom(m_idx)

    try:
        frag_mol = edit.GetMol()
        frags = Chem.GetMolFrags(frag_mol, asMols=False)
    except Exception:
        return [1] * len(ligand_smiles_list)

    # Count coord atoms in each fragment
    denticities = []
    mapped_coord = {idx_map[c] for c in coord_atoms if c in idx_map}
    for frag_atom_ids in frags:
        dent = sum(1 for a in frag_atom_ids if a in mapped_coord)
        denticities.append(max(1, dent))

    # Pad if mismatch (shouldn't happen but safety)
    while len(denticities) < len(ligand_smiles_list):
        denticities.append(1)

    return denticities[:len(ligand_smiles_list)]


def _make_ligand_label(smiles: str, index: int) -> str:
    """Generate a short human-readable label for a ligand."""
    # Common ligand patterns
    known = {
        "c1ccc(-c2ccccn2)nc1": "bipy",
        "c1ccnc(-c2ccccn2)c1": "bipy",
        "CC(=O)[O-]": "acac",
        "[O-]C(C)=O": "acac",
        "[Cl-]": "Cl",
        "[Br-]": "Br",
        "[I-]": "I",
        "[O-]": "OH",
        "N": "NH3",
        "O": "H2O",
        "[NH3]": "NH3",
    }
    canonical = smiles.strip()
    for pattern, label in known.items():
        if canonical == pattern:
            return label
    # Fallback: use first 10 chars
    short = canonical[:10].replace("[", "").replace("]", "").replace("(", "").replace(")", "")
    return short or f"L{index + 1}"


# ---------------------------------------------------------------------------
# 2. Build solvation complex SMILES
# ---------------------------------------------------------------------------

def build_solvation_complex_smiles(
    metals: List[Tuple[str, int]],
    solvent_smiles: str,
    donor_atom: str,
    n_solvent: int,
) -> str:
    """Build SMILES for [M(Solv)_n]^x+ solvation complex.

    Example: Cu2+ + 6 DMSO → [Cu+2](OS(C)C)(OS(C)C)(OS(C)C)(OS(C)C)(OS(C)C)(OS(C)C)
    """
    parts = []
    for symbol, charge in metals:
        if charge > 0:
            charge_str = f"+{charge}" if charge > 1 else "+"
        elif charge < 0:
            charge_str = f"{charge}" if charge < -1 else "-"
        else:
            charge_str = ""
        parts.append(f"[{symbol}{charge_str}]")

    metal_part = "".join(parts)

    # Attach solvent molecules via donor atom
    solvent_bonds = []
    for _ in range(n_solvent):
        solvent_bonds.append(f"({solvent_smiles})")

    return metal_part + "".join(solvent_bonds)


# ---------------------------------------------------------------------------
# 3. SMILES → XYZ conversion for sub-jobs
# ---------------------------------------------------------------------------

def _convert_smiles_and_write(
    smiles: str,
    output_dir: Path,
    converter: str = "NORMAL",
    config: Optional[Dict[str, Any]] = None,
) -> Path:
    """Convert SMILES to XYZ and write start.txt in output_dir."""
    from delfin.pipeline import _run_guppy_for_smiles
    from delfin.smiles_converter import (
        smiles_to_xyz,
        smiles_to_xyz_architector,
        smiles_to_xyz_quick,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    start_path = output_dir / "start.txt"

    if converter == "QUICK":
        xyz_content, error = smiles_to_xyz_quick(smiles)
    elif converter == "ARCHITECTOR":
        xyz_content, error = smiles_to_xyz_architector(smiles)
    elif converter == "GUPPY":
        run_config = dict(config or {})
        _run_guppy_for_smiles(smiles, start_path, run_config)
        xyz_content = start_path.read_text(encoding="utf-8")
        error = None
    else:
        # NORMAL: try full converter, fallback to quick
        xyz_content, error = smiles_to_xyz(smiles)
        if error:
            logger.warning("NORMAL conversion failed for '%s', trying QUICK: %s", smiles, error)
            xyz_content, error = smiles_to_xyz_quick(smiles)

    if error:
        raise ValueError(f"SMILES conversion failed for '{smiles}': {error}")

    start_path.write_text(xyz_content, encoding="utf-8")
    logger.info("[SC] Wrote %s (%d atoms)", start_path, len([l for l in xyz_content.splitlines() if l.strip()]))
    return start_path


# ---------------------------------------------------------------------------
# 4. Preoptimization (xTB / GOAT / CREST)
# ---------------------------------------------------------------------------

def _is_yes_token(value: Any) -> bool:
    return str(value).strip().lower() in {"yes", "true", "1", "on"}


def _normalized_sc_preopt(config: Dict[str, Any]) -> str:
    preopt = str(config.get("sc_preopt", "no")).strip().upper()
    return "" if preopt == "NO" else preopt


def _resolve_main_workflow_converter(config: Dict[str, Any]) -> str:
    from delfin.pipeline import _resolve_smiles_converter

    return _resolve_smiles_converter(config)


def _resolve_main_preopt_steps(config: Dict[str, Any]) -> List[str]:
    steps: List[str] = []
    if _is_yes_token(config.get("XTB_OPT", "no")):
        steps.append("XTB")
    if _is_yes_token(config.get("XTB_GOAT", "no")):
        steps.append("GOAT")
    if _is_yes_token(config.get("CREST", "no")):
        steps.append("CREST")
    return steps


def _run_preopt(workdir: Path, config: Dict[str, Any], preopt: str,
                charge: int, multiplicity: int, solvent: str) -> None:
    """Run preoptimization in workdir using the specified method."""
    if preopt == "no" or not preopt:
        return

    original_cwd = os.getcwd()
    try:
        os.chdir(workdir)

        if preopt.upper() == "XTB":
            from delfin.xtb_crest import XTB
            XTB(multiplicity, charge, config)

        elif preopt.upper() == "GOAT":
            from delfin.xtb_crest import XTB_GOAT
            XTB_GOAT(multiplicity, charge, config)

        elif preopt.upper() == "CREST":
            from delfin.xtb_crest import run_crest_workflow
            run_crest_workflow(
                config.get("PAL", 1),
                solvent,
                charge,
                multiplicity,
                input_file="start.txt",
                crest_dir="CREST",
            )
    finally:
        os.chdir(original_cwd)


def _run_preopt_sequence(
    workdir: Path,
    config: Dict[str, Any],
    steps: List[str],
    charge: int,
    multiplicity: int,
    solvent: str,
    *,
    label: str,
) -> None:
    from delfin.pipeline import _skip_xtb_goat_after_guppy

    for step in steps:
        if step == "GOAT" and _skip_xtb_goat_after_guppy(config):
            logger.info("[SC] Skipping XTB_GOAT for %s: GUPPY already delivered a GOAT-refined geometry.", label)
            continue
        _run_preopt(workdir, config, step, charge, multiplicity, solvent)


# ---------------------------------------------------------------------------
# 5. ORCA input generation (mirrors read_and_modify_file from xyz_io.py)
# ---------------------------------------------------------------------------

def build_orca_input(
    config: Dict[str, Any],
    coord_file: Path,
    output_file: Path,
    charge: int,
    multiplicity: int,
    found_metals: List[str],
    solvent: str,
    broken_sym: str = "",
    include_freq: bool = True,
) -> None:
    """Generate ORCA input file for a stability constant sub-job.

    Uses the same keyword construction as the main DELFIN workflow
    (read_and_modify_file in xyz_io.py) to ensure maximum error cancellation.
    """
    from delfin.utils import resolve_level_of_theory
    from delfin.xyz_io import (
        _build_bang_line,
        _build_freq_block,
        _implicit_token,
        _apply_per_atom_newgto,
        _load_covalent_radii,
        split_qmmm_sections,
    )

    # Read coordinate lines
    with coord_file.open("r") as f:
        coord_lines = [ln for ln in f.readlines() if ln.strip() and ln.strip() != "*"]

    geom_lines, qmmm_range, _ = split_qmmm_sections(coord_lines, coord_file)

    # Resolve level of theory
    main_basisset = str(config.get("main_basisset", "")).strip() or None
    metal_basisset = str(config.get("metal_basisset", "")).strip() or None
    main, metal, rel_token, aux_jk = resolve_level_of_theory(
        found_metals, config, main_basisset, metal_basisset,
    )
    implicit = _implicit_token(config, solvent)

    # Build ! line
    bang = _build_bang_line(
        config,
        rel_token,
        main,
        aux_jk,
        implicit,
        include_freq=include_freq,
        geom_key="geom_opt",
    )

    # Build input
    output_blocks = collect_output_blocks(config, allow=True)
    builder = OrcaInputBuilder(bang)
    builder.add_resources(config.get("maxcore", 6000), config.get("PAL", 1), resolve_maxiter(config))
    builder.add_broken_sym(broken_sym)
    if include_freq:
        builder.add_block(_build_freq_block(config))
    builder.add_blocks(output_blocks)

    lines = builder.lines

    # Geometry block
    lines.append(f"* xyz {charge} {multiplicity}\n")

    # Per-atom NewGTO for metals
    enable_first = str(config.get('first_coordination_sphere_metal_basisset', 'no')).lower() in ('yes', 'true', '1', 'on')
    radii_all = None
    if enable_first:
        radii_all = _load_covalent_radii(config.get("covalent_radii_source", "pyykko2009"))
    geom = _apply_per_atom_newgto(geom_lines, found_metals, metal, config, radii_all)
    lines.extend(geom)
    lines.append("*\n")

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w") as f:
        f.writelines(lines)
    logger.info("[SC] Wrote ORCA input: %s", output_file)


# ---------------------------------------------------------------------------
# 6. Energy extraction from ORCA .out
# ---------------------------------------------------------------------------

def extract_free_energy(out_path: Path) -> SpeciesEnergy:
    """Extract electronic energy and thermal correction from ORCA output."""
    if not out_path.exists():
        raise FileNotFoundError(f"ORCA output not found: {out_path}")

    text = out_path.read_text(encoding="utf-8", errors="replace")

    # Check for successful termination
    if "ORCA TERMINATED NORMALLY" not in text:
        raise RuntimeError(f"ORCA did not terminate normally: {out_path}")

    # Extract FINAL SINGLE POINT ENERGY
    e_el = None
    for match in re.finditer(r"FINAL SINGLE POINT ENERGY\s+([-+]?\d+\.\d+)", text):
        e_el = float(match.group(1))  # take the last one

    if e_el is None:
        raise ValueError(f"Could not extract electronic energy from {out_path}")

    # Extract thermal free energy correction (from FREQ)
    g_rrho = None
    match = re.search(r"Total thermal Free Energy correction\s+([-+]?\d+\.\d+)\s+Eh", text)
    if match:
        g_rrho = float(match.group(1))

    # Also try "Final Gibbs free energy" as alternative
    if g_rrho is None:
        match = re.search(r"Final Gibbs free energy\s+\.\.\.\s+([-+]?\d+\.\d+)\s+Eh", text)
        if match:
            # This is already E_el + G_corr, so compute G_corr
            g_total_direct = float(match.group(1))
            g_rrho = g_total_direct - e_el

    g_total = e_el + g_rrho if g_rrho is not None else None

    return SpeciesEnergy(
        label="",
        folder=str(out_path.parent),
        charge=0,
        multiplicity=1,
        e_el=e_el,
        g_rrho=g_rrho,
        g_total=g_total,
    )


# ---------------------------------------------------------------------------
# 7. Thermodynamic calculation
# ---------------------------------------------------------------------------

def compute_stability_constant(
    analysis: StabilityAnalysis,
    g_complex: SpeciesEnergy,
    g_solv_complex: SpeciesEnergy,
    g_ligands: Dict[str, SpeciesEnergy],  # smiles → energy
    g_solvent: SpeciesEnergy,
    temperature: float = 298.15,
    logK_exp: Optional[float] = None,
) -> StabilityResult:
    """Compute ΔG and log K from free energies of all species."""

    # Check all species have G_total
    all_species = [g_complex, g_solv_complex, g_solvent] + list(g_ligands.values())
    for sp in all_species:
        if sp.g_total is None:
            raise ValueError(
                f"Species '{sp.label}' ({sp.folder}) is missing G_total. "
                "Ensure FREQ was computed for all species."
            )

    n = analysis.n_displaced

    # ΔG = G(complex) + n·G(solvent) - G(solv_complex) - Σ(count_i · G(ligand_i))
    delta_g = g_complex.g_total + n * g_solvent.g_total - g_solv_complex.g_total
    for lig in analysis.unique_ligands:
        count = analysis.ligand_counts[lig.smiles]
        delta_g -= count * g_ligands[lig.smiles].g_total

    delta_g_kcal = delta_g * HARTREE_TO_KCAL
    delta_g_kj = delta_g * HARTREE_TO_KJ
    log_k = -delta_g_kcal / (2.303 * R_KCAL * temperature)

    species = [g_complex, g_solv_complex, g_solvent] + [
        g_ligands[lig.smiles] for lig in analysis.unique_ligands
    ]

    return StabilityResult(
        delta_g_hartree=delta_g,
        delta_g_kcal=delta_g_kcal,
        delta_g_kj=delta_g_kj,
        log_k=log_k,
        temperature=temperature,
        species=species,
        analysis=analysis,
        logK_exp=logK_exp,
    )


# ---------------------------------------------------------------------------
# 8. Report generation
# ---------------------------------------------------------------------------

def write_stability_report(report_path: Path, result: StabilityResult) -> None:
    """Write STABILITY.txt report."""
    a = result.analysis
    lines = []

    sep = "=" * 70
    thin = "-" * 70

    lines.append(sep)
    lines.append("              DELFIN Stability Constant Report")
    lines.append(sep)
    lines.append("")
    lines.append(f"  Complex SMILES:    {a.complex_smiles}")
    metal_str = ", ".join(f"{sym}{'+' * chg if chg > 0 else ''}{chg if chg != 0 else ''}" for sym, chg in a.metals)
    lines.append(f"  Metal:             {metal_str}")
    lines.append(f"  Solvent:           {a.solvent_name}  (SMILES: {a.solvent_smiles}, donor: {a.solvent_donor})")
    lines.append(f"  Explicit Solvent:  {a.n_explicit_solvent}")
    lines.append(f"  Temperature:       {result.temperature:.2f} K")
    lines.append("")

    # Reaction equation
    lines.append(sep)
    lines.append("  Thermodynamic Cycle (Born-Haber):")
    lines.append(sep)
    lines.append("")
    lines.append("               DeltaG_gas")
    lines.append("  Reactants --------------------------------> Products")
    lines.append("      |                                          |")
    lines.append("      | DeltaG_solv(react)          DeltaG_solv(prod)")
    lines.append("      v                                          v")
    lines.append("  Reactants(solv) --------------------> Products(solv)")
    lines.append("                      DeltaG_aq")
    lines.append("")

    # Build reaction string
    reactants = f"[M(Solv)_{a.n_explicit_solvent}]^{a.metal_charge}+"
    for lig in a.unique_ligands:
        count = a.ligand_counts[lig.smiles]
        charge_str = f"^{lig.charge}" if lig.charge != 0 else ""
        count_str = f"{count} " if count > 1 else ""
        reactants += f" + {count_str}{lig.label}{charge_str}"

    n_remaining = a.n_explicit_solvent - a.n_displaced
    complex_charge = a.metal_charge + sum(lig.charge * a.ligand_counts[lig.smiles] for lig in a.unique_ligands)
    products = f"[Complex(Solv)_{n_remaining}]^{complex_charge}"
    if a.n_displaced > 0:
        products += f" + {a.n_displaced} Solv"

    lines.append(f"  {reactants}")
    lines.append(f"    --> {products}")
    lines.append("")

    # Species details
    lines.append(sep)
    lines.append("  Species Details:")
    lines.append(sep)
    lines.append("")

    for sp in result.species:
        lines.append(f"  {sp.label}")
        lines.append(f"    Folder:       {sp.folder}")
        lines.append(f"    Charge:       {sp.charge:+d}")
        lines.append(f"    Multiplicity: {sp.multiplicity}")
        if sp.occupier_mult is not None:
            bs_str = sp.occupier_bs or "none"
            lines.append(f"    OCCUPIER:     Mult = {sp.occupier_mult}, BS = {bs_str}, Idx = {sp.occupier_idx}")
        lines.append(f"    E(el):        {sp.e_el:.10f} Eh" if sp.e_el is not None else "    E(el):        N/A")
        lines.append(f"    G(RRHO):      {sp.g_rrho:+.10f} Eh" if sp.g_rrho is not None else "    G(RRHO):      N/A")
        lines.append(f"    G(total):     {sp.g_total:.10f} Eh" if sp.g_total is not None else "    G(total):     N/A")
        lines.append("")

    # Energy balance
    lines.append(sep)
    lines.append("  Energy Balance:")
    lines.append(sep)
    lines.append("")

    # Build formula
    terms_plus = [f"G(complex)"]
    if a.n_displaced > 0:
        terms_plus.append(f"{a.n_displaced}*G({a.solvent_name})")
    terms_minus = [f"G([M(Solv)_{a.n_explicit_solvent}])"]
    for lig in a.unique_ligands:
        count = a.ligand_counts[lig.smiles]
        if count > 1:
            terms_minus.append(f"{count}*G({lig.label})")
        else:
            terms_minus.append(f"G({lig.label})")

    lines.append(f"  DeltaG = {' + '.join(terms_plus)} - {' - '.join(terms_minus)}")
    lines.append("")
    lines.append(f"  DeltaG = {result.delta_g_hartree:+.10f} Eh")
    lines.append(f"         = {result.delta_g_kcal:+.4f} kcal/mol")
    lines.append(f"         = {result.delta_g_kj:+.4f} kJ/mol")
    lines.append("")

    # Result
    lines.append(thin)
    lines.append("  Result:")
    lines.append(thin)
    lines.append("")
    lines.append(f"  DeltaG(aq)    = {result.delta_g_kcal:+.4f} kcal/mol")
    lines.append(f"  log K(calc)   = {result.log_k:.2f}")
    if result.logK_exp is not None:
        lines.append(f"  log K(exp)    = {result.logK_exp:.2f}  (from CONTROL.txt)")
        lines.append(f"  Delta log K   = {result.log_k - result.logK_exp:+.2f}")
    lines.append("")
    lines.append(sep)

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info("[SC] Report written: %s", report_path)


# ---------------------------------------------------------------------------
# 9. Main phase function
# ---------------------------------------------------------------------------

def run_stability_constant_phase(ctx: "PipelineContext") -> bool:
    """Execute the stability constant calculation phase.

    This is called from cli.py after all other phases complete.
    """
    from delfin.global_scheduler import GlobalOrcaScheduler

    config = ctx.config
    if str(config.get("stability_constant", "no")).strip().lower() != "yes":
        return True

    logger.info("=" * 60)
    logger.info("  Stability Constant Phase")
    logger.info("=" * 60)

    smiles = str(config.get("SMILES", "")).strip()
    if not smiles:
        logger.error("[SC] No SMILES found in CONTROL.txt")
        return False

    solvent = str(config.get("solvent", "")).strip()
    if not solvent:
        logger.error("[SC] No solvent specified in CONTROL.txt")
        return False

    initial_out = _find_initial_output(ctx.control_file_path.parent.resolve(), config)
    if initial_out is None:
        logger.error("[SC] Cannot find initial.out — ensure calc_initial=yes")
        return False
    logger.info("[SC] Using complex energy from: %s", initial_out)

    try:
        plan = build_stability_constant_plan(ctx, resolved_initial_out=initial_out)
    except (ValueError, RuntimeError) as exc:
        logger.error("[SC] Analysis failed: %s", exc)
        return False

    scheduler = GlobalOrcaScheduler(config, label="stability_constant")
    scheduler.add_jobs(plan.jobs)
    result = scheduler.run()

    if "sc_postprocess" in result.completed:
        logger.info("[SC] Stability constant calculation completed successfully")
        return True

    # Log failures
    for job_id, reason in result.failed.items():
        logger.error("[SC] Job %s failed: %s", job_id, reason)
    for job_id, deps in result.skipped.items():
        logger.warning("[SC] Job %s skipped (missing: %s)", job_id, ", ".join(deps))
    return False


def build_stability_constant_plan(
    ctx: "PipelineContext",
    *,
    initial_completion_dependency: Optional[str] = None,
    resolved_initial_out: Optional[Path] = None,
) -> StabilityWorkflowPlan:
    from delfin.parallel_classic_manually import WorkflowJob

    config = ctx.config
    smiles = str(config.get("SMILES", "")).strip()
    if not smiles:
        raise ValueError("No SMILES found in CONTROL.txt")

    solvent = str(config.get("solvent", "")).strip()
    if not solvent:
        raise ValueError("No solvent specified in CONTROL.txt")

    n_explicit = int(str(config.get("n_explicit_solvent", 6)).strip())
    temperature = float(str(config.get("temperature", "298.15")).strip())
    logK_exp_raw = str(config.get("logK_exp", "")).strip()
    logK_exp = float(logK_exp_raw) if logK_exp_raw else None
    sc_converter = str(config.get("sc_smiles_converter", "NORMAL")).strip().upper()
    sc_preopt = _normalized_sc_preopt(config)
    solv_converter = _resolve_main_workflow_converter(config)
    solv_preopt_steps = _resolve_main_preopt_steps(config)

    original_cwd = ctx.control_file_path.parent.resolve()
    sc_dir = original_cwd / "stability_constant"
    sc_dir.mkdir(parents=True, exist_ok=True)

    logger.info("[SC] Analyzing complex SMILES: %s", smiles)
    analysis = analyze_complex(smiles, solvent, n_explicit)
    logger.info("[SC] Metal: %s", analysis.metals)
    logger.info("[SC] Ligands: %s", [(l.label, l.charge, l.denticity) for l in analysis.ligands])
    logger.info("[SC] Unique ligands: %d", len(analysis.unique_ligands))
    logger.info("[SC] Solv complex SMILES: %s", analysis.solv_complex_smiles)
    logger.info("[SC] Displaced solvent: %d", analysis.n_displaced)
    logger.info(
        "[SC] Solv complex workflow mirrors main system: smiles_converter=%s, preopt=%s",
        solv_converter,
        " + ".join(solv_preopt_steps) if solv_preopt_steps else "none",
    )

    total_cores = max(1, int(str(config.get("PAL", 1)).strip()))
    cwd_lock = threading.RLock()
    jobs: List[WorkflowJob] = []

    solv_dir = sc_dir / "solv_complex"

    def _work_solv_occupier(cores: int) -> None:
        local_config = dict(config)
        _run_solv_complex_occupier(
            solv_dir=solv_dir,
            analysis=analysis,
            config=local_config,
            converter=solv_converter,
            preopt_steps=solv_preopt_steps,
            solvent=solvent,
            cwd_lock=cwd_lock,
        )

    jobs.append(WorkflowJob(
        job_id="sc_solv_occ",
        work=_work_solv_occupier,
        description=f"SC: OCCUPIER for [M(Solv)_{n_explicit}]",
        dependencies=set(),
        cores_min=max(1, min(total_cores, 2)),
        cores_optimal=max(2, min(total_cores, total_cores // 2)),
        cores_max=total_cores,
        working_dir=solv_dir,
    ))

    def _work_solv_freq(cores: int) -> None:
        _run_solv_complex_freq(
            solv_dir=solv_dir,
            analysis=analysis,
            config=config,
            solvent=solvent,
            cores=cores,
            cwd_lock=cwd_lock,
        )

    jobs.append(WorkflowJob(
        job_id="sc_solv_freq",
        work=_work_solv_freq,
        description=f"SC: OPT+FREQ [M(Solv)_{n_explicit}]",
        dependencies={"sc_solv_occ"},
        cores_min=max(1, min(total_cores, 2)),
        cores_optimal=total_cores,
        cores_max=total_cores,
        working_dir=solv_dir,
    ))

    for i, lig in enumerate(analysis.unique_ligands):
        safe_label = re.sub(r"[^a-zA-Z0-9_]", "", lig.label) or f"L{i + 1}"
        lig_dir = sc_dir / f"ligand_{i + 1}_{safe_label}"

        def _work_ligand(cores: int, _lig=lig, _dir=lig_dir) -> None:
            _run_closed_shell_species(
                workdir=_dir,
                smiles=_lig.smiles,
                charge=_lig.charge,
                label=f"Ligand {_lig.label}",
                config=config,
                sc_converter=sc_converter,
                sc_preopt=sc_preopt,
                solvent=solvent,
                cores=cores,
                cwd_lock=cwd_lock,
            )

        jobs.append(WorkflowJob(
            job_id=f"sc_ligand_{i + 1}",
            work=_work_ligand,
            description=f"SC: OPT+FREQ ligand {lig.label}",
            dependencies=set(),
            cores_min=max(1, min(total_cores, 2)),
            cores_optimal=max(2, total_cores // 2),
            cores_max=total_cores,
            working_dir=lig_dir,
        ))

    solvent_label = re.sub(r"[^a-zA-Z0-9_]", "", analysis.solvent_name) or "solvent"
    solvent_dir = sc_dir / f"solvent_{solvent_label}"

    def _work_solvent(cores: int) -> None:
        _run_closed_shell_species(
            workdir=solvent_dir,
            smiles=analysis.solvent_smiles,
            charge=0,
            label=f"Solvent {analysis.solvent_name}",
            config=config,
            sc_converter=sc_converter,
            sc_preopt=sc_preopt,
            solvent=solvent,
            cores=cores,
            cwd_lock=cwd_lock,
        )

    jobs.append(WorkflowJob(
        job_id="sc_solvent",
        work=_work_solvent,
        description=f"SC: OPT+FREQ solvent {analysis.solvent_name}",
        dependencies=set(),
        cores_min=1,
        cores_optimal=max(2, total_cores // 4),
        cores_max=total_cores,
        working_dir=solvent_dir,
    ))

    all_dep_ids = {"sc_solv_freq", "sc_solvent"} | {f"sc_ligand_{i + 1}" for i in range(len(analysis.unique_ligands))}
    if initial_completion_dependency:
        all_dep_ids.add(initial_completion_dependency)

    def _work_postprocess(cores: int) -> None:
        _run_postprocessing(
            sc_dir=sc_dir,
            initial_out=resolved_initial_out,
            analysis=analysis,
            config=config,
            temperature=temperature,
            logK_exp=logK_exp,
        )

    jobs.append(WorkflowJob(
        job_id="sc_postprocess",
        work=_work_postprocess,
        description="SC: Compute DeltaG and log K",
        dependencies=all_dep_ids,
        cores_min=1,
        cores_optimal=1,
        cores_max=1,
        inline=True,
        working_dir=sc_dir,
    ))

    return StabilityWorkflowPlan(
        analysis=analysis,
        sc_dir=sc_dir,
        jobs=jobs,
    )


# ---------------------------------------------------------------------------
# Internal worker functions
# ---------------------------------------------------------------------------

def _find_initial_output(cwd: Path, config: Dict[str, Any]) -> Optional[Path]:
    """Locate the initial.out file from the main DELFIN workflow."""
    method = str(config.get("method", "OCCUPIER")).strip()

    if method == "OCCUPIER":
        candidates = [
            cwd / "initial_OCCUPIER" / "initial.out",
            cwd / "initial_OCCUPIER" / "OCCUPIER_initial.out",
        ]
        # Also check for any .out in initial_OCCUPIER
        occ_dir = cwd / "initial_OCCUPIER"
        if occ_dir.is_dir():
            for f in sorted(occ_dir.glob("*.out")):
                if "ORCA TERMINATED NORMALLY" in f.read_text(encoding="utf-8", errors="replace"):
                    candidates.insert(0, f)
                    break
    else:
        candidates = [cwd / "initial.out"]

    for c in candidates:
        if c.exists():
            try:
                text = c.read_text(encoding="utf-8", errors="replace")
                if "ORCA TERMINATED NORMALLY" in text:
                    return c
            except Exception:
                continue
    return None


def _run_solv_complex_occupier(
    solv_dir: Path,
    analysis: StabilityAnalysis,
    config: Dict[str, Any],
    converter: str,
    preopt_steps: List[str],
    solvent: str,
    cwd_lock: threading.RLock,
) -> None:
    """Set up solvation complex folder and run OCCUPIER (like initial_OCCUPIER)."""
    from delfin.occupier import run_OCCUPIER

    solv_dir.mkdir(parents=True, exist_ok=True)
    occ_dir = solv_dir / "solv_OCCUPIER"
    occ_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: SMILES → XYZ
    start_path = solv_dir / "start.txt"
    if not start_path.exists():
        _convert_smiles_and_write(
            analysis.solv_complex_smiles, solv_dir, converter=converter, config=config,
        )

    # Step 2: Mirror the main-system preopt sequence before OCCUPIER.
    if preopt_steps:
        mult_guess = 1  # initial guess, OCCUPIER will determine the right one
        _run_preopt_sequence(
            solv_dir,
            config,
            preopt_steps,
            analysis.metal_charge,
            mult_guess,
            solvent,
            label="solvation complex",
        )

    # Step 3: Prepare OCCUPIER folder (similar to prepare_occ_folder_only_setup)
    # Copy start.txt → input.xyz in OCCUPIER folder
    input_xyz = occ_dir / "input.xyz"
    coord_text = (solv_dir / "start.txt").read_text(encoding="utf-8")
    coord_lines = [l for l in coord_text.splitlines() if l.strip()]
    n_atoms = len(coord_lines)
    with input_xyz.open("w", encoding="utf-8") as f:
        f.write(f"{n_atoms}\n\n")
        f.write("\n".join(coord_lines) + "\n")

    # Copy input.xyz to input0.xyz (backup)
    shutil.copy(input_xyz, occ_dir / "input0.xyz")

    # Create CONTROL.txt for OCCUPIER folder with solv complex charge
    _write_occupier_control(occ_dir, config, charge=analysis.metal_charge)

    # Step 4: Run OCCUPIER
    with cwd_lock:
        prev_cwd = os.getcwd()
    try:
        os.chdir(occ_dir)
        logger.info("[SC] Running OCCUPIER in %s", occ_dir)
        run_OCCUPIER()
    finally:
        os.chdir(prev_cwd)

    logger.info("[SC] OCCUPIER completed for solvation complex")


def _write_occupier_control(occ_dir: Path, config: Dict[str, Any], charge: int) -> None:
    """Write a CONTROL.txt for an OCCUPIER sub-folder, inheriting main settings."""
    # Find parent CONTROL.txt
    parent_control = Path.cwd() / "CONTROL.txt"
    if not parent_control.exists():
        # Try two levels up (from within stability_constant/solv_complex/solv_OCCUPIER/)
        for p in [occ_dir.parent.parent.parent, occ_dir.parent.parent]:
            candidate = p / "CONTROL.txt"
            if candidate.exists():
                parent_control = candidate
                break

    if parent_control.exists():
        shutil.copy(parent_control, occ_dir / "CONTROL.txt")
    else:
        raise FileNotFoundError("Cannot find parent CONTROL.txt for OCCUPIER setup")

    # Update charge and input_file in the copy
    control_path = occ_dir / "CONTROL.txt"
    text = control_path.read_text(encoding="utf-8")
    text = re.sub(r"charge=[+-]?\d+", f"charge={charge}", text)
    text = re.sub(r"input_file=\S+", "input_file=input.xyz", text)
    control_path.write_text(text, encoding="utf-8")


def _run_solv_complex_freq(
    solv_dir: Path,
    analysis: StabilityAnalysis,
    config: Dict[str, Any],
    solvent: str,
    cores: int,
    cwd_lock: threading.RLock,
) -> None:
    """Run OPT+FREQ for solvation complex using OCCUPIER-determined multiplicity."""
    from delfin.copy_helpers import read_occupier_file
    from delfin.orca import run_orca
    from delfin.utils import search_transition_metals
    from delfin import smart_recalc

    occ_dir = solv_dir / "solv_OCCUPIER"

    # Read OCCUPIER result
    mult, broken_sym, preferred_idx, gbw_path = read_occupier_file(
        str(occ_dir), "OCCUPIER.txt",
        multiplicity=1, broken_sym="",
        min_fspe_index=None,
        config=config,
    )

    if mult is None:
        raise RuntimeError(f"Could not read OCCUPIER results from {occ_dir / 'OCCUPIER.txt'}")

    logger.info("[SC] Solvation complex OCCUPIER result: mult=%d, BS='%s', idx=%s", mult, broken_sym, preferred_idx)

    # Find coordinate file (OCCUPIER may have produced optimized geometry)
    coord_file = solv_dir / f"input_solv_OCCUPIER.xyz"
    if not coord_file.exists():
        coord_file = solv_dir / "start.txt"
    if not coord_file.exists():
        coord_file = occ_dir / f"input{preferred_idx}.xyz" if preferred_idx else occ_dir / "input.xyz"

    # For start.txt format (no header), use directly; for .xyz format, strip header
    if coord_file.suffix == ".xyz":
        lines = coord_file.read_text(encoding="utf-8").splitlines()
        # Strip XYZ header (natoms + comment line)
        coord_lines = [l for l in lines[2:] if l.strip()]
        coord_path = solv_dir / "solv_complex_coords.txt"
        coord_path.write_text("\n".join(coord_lines) + "\n", encoding="utf-8")
    else:
        coord_path = coord_file

    # Detect metals in coordinate file
    try:
        metals = search_transition_metals(str(coord_path))
    except Exception:
        metals = [sym for sym, _ in analysis.metals]

    # Build ORCA input
    inp_path = solv_dir / "solv_complex.inp"
    out_path = solv_dir / "solv_complex.out"

    if smart_recalc.should_skip(inp_path, out_path):
        logger.info("[SC] Skipping solv_complex OPT+FREQ (recalc, inputs unchanged)")
        return

    # Override PAL with allocated cores
    sc_config = dict(config)
    sc_config["PAL"] = cores

    build_orca_input(
        config=sc_config,
        coord_file=coord_path,
        output_file=inp_path,
        charge=analysis.metal_charge,
        multiplicity=mult,
        found_metals=metals,
        solvent=solvent,
        broken_sym=broken_sym or "",
        include_freq=True,
    )

    # Add GBW reference if available
    if gbw_path and Path(gbw_path).exists():
        text = inp_path.read_text(encoding="utf-8")
        gbw_line = f'%moinp "{gbw_path}"\n'
        text = text.replace("\n* xyz", f"\n{gbw_line}* xyz", 1)
        inp_path.write_text(text, encoding="utf-8")

    # Run ORCA
    run_orca(str(inp_path), cores)
    smart_recalc.store_fingerprint(inp_path)

    # Verify success
    if not out_path.exists() or "ORCA TERMINATED NORMALLY" not in out_path.read_text(encoding="utf-8", errors="replace"):
        raise RuntimeError(f"ORCA failed for solvation complex: {out_path}")

    logger.info("[SC] Solvation complex OPT+FREQ completed: %s", out_path)


def _run_closed_shell_species(
    workdir: Path,
    smiles: str,
    charge: int,
    label: str,
    config: Dict[str, Any],
    sc_converter: str,
    sc_preopt: str,
    solvent: str,
    cores: int,
    cwd_lock: threading.RLock,
) -> None:
    """Run OPT+FREQ for a closed-shell species (ligand or solvent)."""
    from delfin.orca import run_orca
    from delfin.utils import search_transition_metals
    from delfin.smiles_converter import contains_metal
    from delfin import smart_recalc

    workdir.mkdir(parents=True, exist_ok=True)

    inp_path = workdir / "calc.inp"
    out_path = workdir / "calc.out"

    if smart_recalc.should_skip(inp_path, out_path):
        logger.info("[SC] Skipping %s (recalc, inputs unchanged)", label)
        return

    # SMILES → XYZ
    start_path = workdir / "start.txt"
    if not start_path.exists():
        _convert_smiles_and_write(smiles, workdir, converter=sc_converter)

    # Preopt
    if sc_preopt:
        _run_preopt(workdir, config, sc_preopt, charge, 1, solvent)

    # Detect metals in coordinates (for basis set selection)
    has_metal = contains_metal(smiles)
    if has_metal:
        try:
            metals = search_transition_metals(str(start_path))
        except Exception:
            metals = []
    else:
        metals = []

    # Build ORCA input (closed-shell, mult=1)
    sc_config = dict(config)
    sc_config["PAL"] = cores

    build_orca_input(
        config=sc_config,
        coord_file=start_path,
        output_file=inp_path,
        charge=charge,
        multiplicity=1,
        found_metals=metals,
        solvent=solvent,
        broken_sym="",
        include_freq=True,
    )

    # Run ORCA
    run_orca(str(inp_path), cores)
    smart_recalc.store_fingerprint(inp_path)

    # Verify
    if not out_path.exists() or "ORCA TERMINATED NORMALLY" not in out_path.read_text(encoding="utf-8", errors="replace"):
        raise RuntimeError(f"ORCA failed for {label}: {out_path}")

    logger.info("[SC] %s OPT+FREQ completed: %s", label, out_path)


def _run_postprocessing(
    sc_dir: Path,
    initial_out: Optional[Path],
    analysis: StabilityAnalysis,
    config: Dict[str, Any],
    temperature: float,
    logK_exp: Optional[float],
) -> None:
    """Collect all energies and compute stability constant."""
    from delfin.copy_helpers import read_occupier_file

    if initial_out is None:
        initial_out = _find_initial_output(sc_dir.parent, config)
    if initial_out is None:
        raise FileNotFoundError("Could not locate initial.out for stability constant post-processing")

    # 1. Complex energy (from main workflow initial.out)
    logger.info("[SC] Extracting complex energy from %s", initial_out)
    g_complex = extract_free_energy(initial_out)
    g_complex.label = f"[Complex] {analysis.complex_smiles}"
    g_complex.folder = str(initial_out.parent.relative_to(sc_dir.parent))
    g_complex.charge = analysis.metal_charge + sum(
        lig.charge * analysis.ligand_counts[lig.smiles] for lig in analysis.unique_ligands
    )

    # Try to get OCCUPIER info for complex
    occ_dir = sc_dir.parent / "initial_OCCUPIER"
    if occ_dir.is_dir() and (occ_dir / "OCCUPIER.txt").exists():
        mult, bs, idx, _ = read_occupier_file(
            str(occ_dir), "OCCUPIER.txt", 1, "", None, config, verbose=False,
        )
        if mult is not None:
            g_complex.multiplicity = mult
            g_complex.occupier_mult = mult
            g_complex.occupier_bs = bs
            g_complex.occupier_idx = idx

    # 2. Solvation complex energy
    solv_out = sc_dir / "solv_complex" / "solv_complex.out"
    logger.info("[SC] Extracting solvation complex energy from %s", solv_out)
    g_solv = extract_free_energy(solv_out)
    g_solv.label = f"[M(Solv)_{analysis.n_explicit_solvent}]^{analysis.metal_charge}+"
    g_solv.folder = "stability_constant/solv_complex"
    g_solv.charge = analysis.metal_charge

    # Get OCCUPIER info for solv complex
    solv_occ_dir = sc_dir / "solv_complex" / "solv_OCCUPIER"
    if solv_occ_dir.is_dir() and (solv_occ_dir / "OCCUPIER.txt").exists():
        mult, bs, idx, _ = read_occupier_file(
            str(solv_occ_dir), "OCCUPIER.txt", 1, "", None, config, verbose=False,
        )
        if mult is not None:
            g_solv.multiplicity = mult
            g_solv.occupier_mult = mult
            g_solv.occupier_bs = bs
            g_solv.occupier_idx = idx

    # 3. Ligand energies
    g_ligands: Dict[str, SpeciesEnergy] = {}
    for i, lig in enumerate(analysis.unique_ligands):
        safe_label = re.sub(r"[^a-zA-Z0-9_]", "", lig.label) or f"L{i + 1}"
        lig_out = sc_dir / f"ligand_{i + 1}_{safe_label}" / "calc.out"
        logger.info("[SC] Extracting ligand %s energy from %s", lig.label, lig_out)
        g_lig = extract_free_energy(lig_out)
        g_lig.label = f"{lig.label} (charge: {lig.charge}, dent: {lig.denticity}, count: {analysis.ligand_counts[lig.smiles]})"
        g_lig.folder = f"stability_constant/ligand_{i + 1}_{safe_label}"
        g_lig.charge = lig.charge
        g_lig.multiplicity = 1
        g_ligands[lig.smiles] = g_lig

    # 4. Solvent energy
    solvent_label = re.sub(r"[^a-zA-Z0-9_]", "", analysis.solvent_name) or "solvent"
    solv_mol_out = sc_dir / f"solvent_{solvent_label}" / "calc.out"
    logger.info("[SC] Extracting solvent energy from %s", solv_mol_out)
    g_solvent = extract_free_energy(solv_mol_out)
    g_solvent.label = f"{analysis.solvent_name} (Solvent, displaced: {analysis.n_displaced})"
    g_solvent.folder = f"stability_constant/solvent_{solvent_label}"
    g_solvent.charge = 0
    g_solvent.multiplicity = 1

    # 5. Compute stability constant
    logger.info("[SC] Computing stability constant...")
    result = compute_stability_constant(
        analysis=analysis,
        g_complex=g_complex,
        g_solv_complex=g_solv,
        g_ligands=g_ligands,
        g_solvent=g_solvent,
        temperature=temperature,
        logK_exp=logK_exp,
    )

    # 6. Write report
    report_path = sc_dir / "STABILITY.txt"
    write_stability_report(report_path, result)

    # 7. Print summary
    logger.info("[SC] " + "=" * 50)
    logger.info("[SC]  DeltaG = %.4f kcal/mol", result.delta_g_kcal)
    logger.info("[SC]  log K  = %.2f", result.log_k)
    if logK_exp is not None:
        logger.info("[SC]  log K(exp) = %.2f  (Delta = %+.2f)", logK_exp, result.log_k - logK_exp)
    logger.info("[SC] " + "=" * 50)

    print(f"\n  Stability Constant Result:")
    print(f"  DeltaG = {result.delta_g_kcal:+.4f} kcal/mol")
    print(f"  log K  = {result.log_k:.2f}")
    if logK_exp is not None:
        print(f"  log K(exp) = {logK_exp:.2f}  (Delta = {result.log_k - logK_exp:+.2f})")
    print(f"  Report: {report_path}\n")
