"""Build up metal complex step by step using ORCA XTB DOCKER.

This module provides functionality to:
1. Parse a metal complex SMILES and extract metal + ligands
2. Convert ligands to XYZ files
3. Create ORCA input files for stepwise complex build-up using XTB DOCKER
4. Run the ORCA jobs sequentially

Workflow:
- Metal starts at origin (0, 0, 0)
- Each ligand is docked onto the growing complex using ORCA's XTB DOCKER method
- After each step, the optimized structure becomes the host for the next ligand
"""

from __future__ import annotations

import argparse
import math
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import Optional, Tuple, List, Dict

from delfin.common.logging import configure_logging, get_logger
from delfin.orca import run_orca, find_orca_executable
from delfin.smiles_converter import RDKIT_AVAILABLE, contains_metal, _METALS, _METAL_SET
from delfin.xyz_io import _load_covalent_radii, _COVALENT_RADII_FALLBACK

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem

logger = get_logger(__name__)



def _find_repo_root(start: Path) -> Path | None:
    cur = start.resolve()
    for _ in range(8):
        if (cur / ".git").exists():
            return cur
        if cur == cur.parent:
            break
        cur = cur.parent
    return None


def _get_git_commit() -> str:
    repo = _find_repo_root(Path(__file__).parent)
    if not repo:
        return "unknown"
    try:
        import subprocess
        out = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo)
        return out.decode().strip()
    except Exception:
        return "unknown"


def _write_build_info(builder_dir: Path, cmdline: str | None) -> None:
    try:
        from datetime import datetime
        import socket
        from delfin import __version__

        info_path = builder_dir / "infos.txt"
        lines = []
        lines.append(f"timestamp: {datetime.now().isoformat(timespec='seconds')}")
        lines.append(f"host: {socket.gethostname()}")
        lines.append(f"delfin_version: {__version__}")
        lines.append(f"git_commit: {_get_git_commit()}")
        if cmdline:
            lines.append(f"command: {cmdline}")
        info_path.write_text("\n".join(lines) + "\n")
    except Exception:
        # Best-effort only
        pass


def extract_ligands_from_complex(smiles: str) -> Tuple[
    Optional[List[Tuple[str, int]]],  # metals: list of (symbol, charge)
    Optional[List[str]],  # ligand_smiles_list
    Optional[str],  # error_message
]:
    """Extract metals and ligands from a metal complex SMILES.

    Based on ChemDarwin's extract_ligands_from_complex function.
    Supports single and multi-metal complexes.

    Args:
        smiles: SMILES string of the metal complex

    Returns:
        Tuple of (metals, ligand_smiles_list, error_message)
        - metals: List of (symbol, charge) tuples for each metal
        - ligand_smiles_list: List of SMILES strings for each ligand
        - error_message: Error description if failed, None on success
    """
    if not RDKIT_AVAILABLE:
        return None, None, "RDKit not available"

    # Parse without sanitization for metal complexes
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None, None, "Could not parse SMILES"

    # Find metal atom(s)
    metal_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in _METALS:
            metal_atoms.append(atom.GetIdx())

    if not metal_atoms:
        return None, None, "No metal found in SMILES"

    # Collect metal info
    metals = []
    for metal_idx in metal_atoms:
        metal_atom = mol.GetAtomWithIdx(metal_idx)
        metals.append((metal_atom.GetSymbol(), metal_atom.GetFormalCharge()))

    # Capture original atom properties FIRST
    orig_props = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        orig_props[idx] = {
            'formal_charge': atom.GetFormalCharge(),
            'explicit_h': atom.GetNumExplicitHs(),
            'no_implicit': atom.GetNoImplicit(),
            'rad_e': atom.GetNumRadicalElectrons(),
            'symbol': atom.GetSymbol(),
        }

    # Create editable mol
    edit_mol = Chem.RWMol(mol)

    # Find coordinating atoms (atoms bonded to metals) BEFORE removing bonds
    coordinating_atoms = set()
    bonds_to_remove = []
    for bond in edit_mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in metal_atoms:
            coordinating_atoms.add(end_idx)
            bonds_to_remove.append((begin_idx, end_idx))
        elif end_idx in metal_atoms:
            coordinating_atoms.add(begin_idx)
            bonds_to_remove.append((begin_idx, end_idx))

    # Remove bonds to ALL metals
    for b in bonds_to_remove:
        edit_mol.RemoveBond(b[0], b[1])

    # Remove metal atoms (in reverse order to preserve indices)
    for metal_idx in sorted(metal_atoms, reverse=True):
        edit_mol.RemoveAtom(metal_idx)

    # Get fragments (ligands)
    try:
        mol_no_metal = edit_mol.GetMol()

        # Build index mapping: old_idx -> new_idx after metal removal
        removed_count = 0
        idx_map = {}
        for old_idx in range(len(orig_props)):
            if old_idx in metal_atoms:
                removed_count += 1
            else:
                idx_map[old_idx] = old_idx - removed_count

        # Apply original atom properties
        # For coordinating atoms: we need to handle the lost metal bond
        for old_idx, props in orig_props.items():
            if old_idx in metal_atoms:
                continue
            new_idx = idx_map.get(old_idx)
            if new_idx is None:
                continue
            try:
                a = mol_no_metal.GetAtomWithIdx(new_idx)
                a.SetNumRadicalElectrons(props['rad_e'])

                if old_idx in coordinating_atoms:
                    # Coordinating atom lost a bond to metal
                    # Set NoImplicit=True and explicit Hs to 0 to prevent H addition
                    a.SetFormalCharge(props['formal_charge'])
                    a.SetNumExplicitHs(0)
                    a.SetNoImplicit(True)
                else:
                    a.SetFormalCharge(props['formal_charge'])
                    a.SetNumExplicitHs(props['explicit_h'])
                    a.SetNoImplicit(props['no_implicit'])
            except Exception:
                pass

        try:
            mol_no_metal.UpdatePropertyCache(strict=False)
        except Exception:
            pass

        frags = Chem.GetMolFrags(mol_no_metal, asMols=True, sanitizeFrags=False)

        ligand_smiles = []
        for frag_mol in frags:
            # Remove atom maps for clean SMILES
            for atom in frag_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            smi = Chem.MolToSmiles(frag_mol, canonical=True)
            ligand_smiles.append(smi)

        return metals, ligand_smiles, None

    except Exception as e:
        return None, None, str(e)


def get_ligand_charge(smiles: str) -> int:
    """Calculate the total charge of a ligand from its SMILES.

    Args:
        smiles: SMILES string of the ligand

    Returns:
        Total formal charge of the ligand
    """
    if not RDKIT_AVAILABLE:
        return 0

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return 0

    total_charge = 0
    for atom in mol.GetAtoms():
        total_charge += atom.GetFormalCharge()

    return total_charge


def get_ligand_size(smiles: str) -> int:
    """Calculate the size of a ligand (number of atoms including hydrogens).

    Args:
        smiles: SMILES string of the ligand

    Returns:
        Number of atoms (including implicit hydrogens)
    """
    if not RDKIT_AVAILABLE:
        return 0

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return 0

        # Add hydrogens to get total atom count
        try:
            mol = Chem.AddHs(mol)
        except Exception:
            pass

        return mol.GetNumAtoms()
    except Exception:
        return 0


def _is_halide_or_pseudohalide(smiles: str) -> bool:
    """Heuristic to flag halides and common pseudohalides for late docking."""
    smi = smiles.strip()

    # Halides (monoatomic anions)
    halides = {"[F-]", "[Cl-]", "[Br-]", "[I-]"}
    if smi in halides:
        return True

    # Common pseudohalides: CN-, SCN-, OCN-, N3- (various SMILES orderings)
    pseudohalides = {
        "[N-]#C", "[C-]#N",          # cyanide
        "N#C[S-]", "[S-]C#N",        # thiocyanate
        "N#C[O-]", "[O-]C#N",        # cyanate
        "[N-]=[N+]=N", "[N-][N+]#N", "N=[N+]=[N-]",  # azide
    }
    return smi in pseudohalides


def ligand_to_xyz(smiles: str, output_path: Path) -> Tuple[bool, Optional[str]]:
    """Convert a ligand SMILES to XYZ file.

    Handles coordinating atoms (from metal complexes) that should not receive H atoms.

    Args:
        smiles: SMILES string of the ligand
        output_path: Path to write the XYZ file

    Returns:
        Tuple of (success, error_message)
    """
    if not RDKIT_AVAILABLE:
        return False, "RDKit not available"

    try:
        # First, try normal parsing (sanitized) - works for most ligands
        mol = Chem.MolFromSmiles(smiles)
        atoms_no_h = []

        if mol is None:
            # Try parsing without sanitization (for complex ligands with unusual valences)
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                # Find atoms that should NOT get H added (NoImplicit=True)
                for atom in mol.GetAtoms():
                    if atom.GetNoImplicit():
                        atoms_no_h.append(atom.GetIdx())

                # Try to sanitize (may partially fail, but that's ok)
                try:
                    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                                     Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                except Exception:
                    pass

                try:
                    mol.UpdatePropertyCache(strict=False)
                except Exception:
                    pass

        if mol is None:
            return False, f"Could not parse SMILES: {smiles}"

        # Add hydrogens only to atoms that are not marked as NoImplicit
        try:
            if atoms_no_h:
                # Get list of atoms that CAN receive H
                all_atoms = set(range(mol.GetNumAtoms()))
                atoms_with_h = list(all_atoms - set(atoms_no_h))
                if atoms_with_h:
                    mol = Chem.AddHs(mol, onlyOnAtoms=atoms_with_h)
            else:
                mol = Chem.AddHs(mol)
        except Exception:
            pass

        # Generate 3D coordinates
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)

        if result != 0:
            # Try with random coordinates
            result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)

        if result != 0:
            return False, f"Could not generate 3D coordinates for: {smiles}"

        # Try UFF optimization (may fail for some molecules)
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            # Try MMFF as fallback
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                pass  # Skip optimization if both fail

        # Write XYZ file
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()

        lines = [str(num_atoms), f"Ligand from SMILES: {smiles}"]
        for i in range(num_atoms):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            symbol = atom.GetSymbol()
            lines.append(f"{symbol:4s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")

        output_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
        return True, None

    except Exception as e:
        return False, str(e)


def create_metal_xyz(metal_symbol: str, output_path: Path) -> None:
    """Create XYZ file with just the metal at origin.

    Args:
        metal_symbol: Element symbol of the metal
        output_path: Path to write the XYZ file
    """
    content = f"1\nMetal center\n{metal_symbol:4s}     0.000000     0.000000     0.000000\n"
    output_path.write_text(content, encoding='utf-8')


def create_docker_input(
    host_xyz: str,
    guest_xyz: str,
    guest_charge: int,
    host_charge: int,
    multiplicity: int,
    output_path: Path,
    pal: int = 32,
    maxcore: int = 1000,
) -> None:
    """Create ORCA input file for XTB DOCKER calculation.

    Args:
        host_xyz: Filename of the host structure (e.g., "step_0.xyz")
        guest_xyz: Filename of the guest ligand (e.g., "ligand1.xyz")
        guest_charge: Formal charge of the guest ligand
        host_charge: Formal charge of the host
        multiplicity: Spin multiplicity
        output_path: Path to write the ORCA input file
        pal: Number of parallel processes (default: 32)
        maxcore: Memory per core in MB (default: 1000)
    """
    content = f"""! XTB2 ALPB(DMF)

%PAL NPROCS {pal} END
%MAXCORE {maxcore}

%DOCKER
    DOCKLEVEL      COMPLETE   
    GUEST          "{guest_xyz}"
    GuestCharge    {guest_charge}
END

* XYZFile {host_charge} {multiplicity} {host_xyz}
"""
    output_path.write_text(content, encoding='utf-8')


def create_goat_input(
    xyz_file: str,
    charge: int,
    multiplicity: int,
    output_path: Path,
    pal: int = 32,
    maxcore: int = 1000,
    include_goat_block: bool = True,
    uphill_atoms: Optional[List[int]] = None,
    topobreak_atoms: Optional[List[int]] = None,
    restrictive: bool = False,
    skip_initial_opt: bool = False,
) -> None:
    """Create ORCA input file for XTB2 GOAT global optimization.

    Args:
        xyz_file: Filename of the structure to optimize
        charge: Total charge
        multiplicity: Spin multiplicity
        output_path: Path to write the ORCA input file
        pal: Number of parallel processes
        maxcore: Memory per core in MB
    """
    content = """! XTB2 GOAT ALPB(DMF)

%PAL NPROCS {pal} END
%MAXCORE {maxcore}

* XYZFile {charge} {multiplicity} {xyz_file}
"""
    if include_goat_block:
        goat_lines = ["%GOAT\n"]
        if skip_initial_opt:
            goat_lines.append("  SKIPINITIALOPT TRUE\n")
        if not restrictive:
            if uphill_atoms:
                uphill_str = _format_atom_indices(uphill_atoms)
                # Escape braces to avoid str.format() interpreting them
                goat_lines.append(f"  UPHILLATOMS {{{{{uphill_str}}}}} END\n")
            goat_lines.append("  FREEFRAGMENTS FALSE\n")
            goat_lines.append("  FREEZEBONDS FALSE\n")
            goat_lines.append("  FREEZEANGLES FALSE\n")
            goat_lines.append("  FREEHETEROATOMS FALSE\n")
            goat_lines.append("  MAXTOPODIFF 2\n")
            if topobreak_atoms:
                topobreak_str = _format_atom_indices(topobreak_atoms)
                goat_lines.append(f"  TOPOBREAK {{{{{topobreak_str}}}}} END\n")
        goat_lines.append("END\n\n")
        content = content.replace(
            "* XYZFile {charge} {multiplicity} {xyz_file}\n",
            "".join(goat_lines) + "* XYZFile {charge} {multiplicity} {xyz_file}\n",
        )
    content = content.format(
        pal=pal,
        maxcore=maxcore,
        charge=charge,
        multiplicity=multiplicity,
        xyz_file=xyz_file,
    )
    output_path.write_text(content, encoding='utf-8')


def _orca_terminated_normally(out_path: Path) -> bool:
    if not out_path.exists():
        return False
    try:
        with out_path.open("rb") as f:
            f.seek(0, 2)
            size = f.tell()
            # Read last 8KB (or full file if smaller)
            f.seek(max(size - 8192, 0))
            tail = f.read().decode("utf-8", errors="ignore")
        return "****ORCA TERMINATED NORMALLY****" in tail
    except Exception:
        return False


def run_goat_optimization(
    builder_dir: Path,
    xyz_file: str,
    charge: int,
    multiplicity: int,
    pal: int,
    maxcore: int,
    prefix: str = "goat",
    include_goat_block: bool = True,
    uphill_atoms: Optional[List[int]] = None,
    uphill_from_xyz: bool = False,
    uphill_include_metals: bool = True,
    restrictive: bool = False,
    skip_initial_opt: bool = False,
) -> Tuple[bool, Optional[str]]:
    """Run GOAT global optimization on a structure.

    Args:
        builder_dir: Working directory
        xyz_file: Input XYZ filename (relative to builder_dir)
        charge: Total charge
        multiplicity: Spin multiplicity
        pal: Number of parallel processes
        maxcore: Memory per core in MB
        prefix: Prefix for output files

    Returns:
        Tuple of (success, optimized_xyz_filename or None)
    """
    goat_input = builder_dir / f"{prefix}.inp"
    goat_output = builder_dir / f"{prefix}.out"
    goat_result = builder_dir / f"{prefix}.globalminimum.xyz"

    topobreak_atoms: Optional[List[int]] = None
    if include_goat_block and uphill_from_xyz:
        uphill_atoms = _compute_uphill_atoms(
            builder_dir / xyz_file,
            include_metals=uphill_include_metals,
        )
        # Always set TOPOBREAK to UPHILLATOMS including metals.
        topobreak_atoms = _compute_uphill_atoms(
            builder_dir / xyz_file,
            include_metals=True,
        )
    create_goat_input(
        xyz_file=xyz_file,
        charge=charge,
        multiplicity=multiplicity,
        output_path=goat_input,
        pal=pal,
        maxcore=maxcore,
        include_goat_block=include_goat_block,
        uphill_atoms=uphill_atoms,
        topobreak_atoms=topobreak_atoms,
        restrictive=restrictive,
        skip_initial_opt=skip_initial_opt,
    )

    # Recalc logic: skip GOAT if previous output terminated normally and result exists
    if _orca_terminated_normally(goat_output) and goat_result.exists():
        logger.info(f"Skipping GOAT: found normal termination in {goat_output.name}")
        shutil.copy(goat_result, builder_dir / xyz_file)
        return True, xyz_file

    logger.info(f"Running GOAT optimization on {xyz_file}...")
    success = run_orca(
        str(goat_input),
        str(goat_output),
        working_dir=builder_dir,
    )

    if not success:
        logger.error(f"GOAT optimization failed for {xyz_file}")
        return False, None

    if not goat_result.exists():
        logger.error(f"GOAT result file not found: {goat_result}")
        return False, None

    # Copy result back to original xyz file
    shutil.copy(goat_result, builder_dir / xyz_file)
    logger.info(f"GOAT optimization completed, updated {xyz_file}")

    return True, xyz_file


def parse_xyz_to_content(xyz_path: Path) -> Tuple[int, str, List[str]]:
    """Parse XYZ file and return atom count, comment, and coordinate lines.

    Args:
        xyz_path: Path to XYZ file

    Returns:
        Tuple of (atom_count, comment, coordinate_lines)
    """
    lines = xyz_path.read_text(encoding='utf-8').strip().split('\n')
    atom_count = int(lines[0].strip())
    comment = lines[1].strip() if len(lines) > 1 else ""
    coords = lines[2:2 + atom_count] if len(lines) > 2 else []
    return atom_count, comment, coords


def _extract_lowest_energy_xyz_from_trajectory(traj_xyz: Path, out_xyz: Path) -> bool:
    """Extract the lowest-energy structure from a multi-XYZ trajectory file."""
    if not traj_xyz.exists():
        return False
    try:
        lines = traj_xyz.read_text(encoding="utf-8", errors="ignore").splitlines()
        if not lines:
            return False
        i = 0
        best_energy = None
        best_block = None
        energy_re = re.compile(r"Epreopt=([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)")
        while i < len(lines):
            try:
                atom_count = int(lines[i].strip())
            except Exception:
                break
            if atom_count <= 0:
                break
            start = i
            end = i + 2 + atom_count
            if end > len(lines):
                break
            comment = lines[i + 1] if i + 1 < len(lines) else ""
            energy = None
            m = energy_re.search(comment)
            if m:
                try:
                    energy = float(m.group(1))
                except Exception:
                    energy = None
            if energy is not None:
                if best_energy is None or energy < best_energy:
                    best_energy = energy
                    best_block = lines[start:end]
            elif best_block is None:
                # Fallback: keep first block if no energies found yet
                best_block = lines[start:end]
            i = end
        if not best_block:
            return False
        out_xyz.write_text("\n".join(best_block) + "\n", encoding="utf-8")
        return True
    except Exception:
        return False


def _format_atom_indices(indices: List[int]) -> str:
    """Format 1-based atom indices as ORCA range string (e.g., '1:3 5 7:9')."""
    if not indices:
        return ""
    idx = sorted(set(int(i) for i in indices))
    parts: List[str] = []
    start = prev = idx[0]
    for cur in idx[1:]:
        if cur == prev + 1:
            prev = cur
            continue
        parts.append(f"{start}:{prev}" if start != prev else f"{start}")
        start = prev = cur
    parts.append(f"{start}:{prev}" if start != prev else f"{start}")
    return " ".join(parts)


def _compute_uphill_atoms(
    xyz_path: Path,
    scale: float = 1.3,
    sphere_depth: int = 2,
    include_metals: bool = True,
    exclude_h: bool = True,
) -> List[int]:
    """Return 1-based indices of metals + coordination spheres up to sphere_depth."""
    _, _, coord_lines = parse_xyz_to_content(xyz_path)
    atoms = []
    for line in coord_lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        elem = re.match(r"([A-Za-z]{1,2})", parts[0])
        if not elem:
            continue
        sym = elem.group(1)
        try:
            x, y, z = map(float, parts[1:4])
        except ValueError:
            continue
        atoms.append({"elem": sym, "x": x, "y": y, "z": z})

    if not atoms:
        return []

    radii = _load_covalent_radii("pyykko2009") or _COVALENT_RADII_FALLBACK

    def _rcov(sym: str) -> float:
        return float(radii.get(sym, 1.20))

    def _dist(a, b) -> float:
        return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)

    metal_indices = [i for i, a in enumerate(atoms) if a["elem"] in _METAL_SET]
    if not metal_indices:
        return []

    # iterative sphere expansion
    selected = set(metal_indices if include_metals else [])
    frontier = set(metal_indices)
    depth = max(1, int(sphere_depth))
    for _ in range(depth):
        next_frontier = set()
        for im in frontier:
            m = atoms[im]
            r_m = _rcov(m["elem"])
            for i, a in enumerate(atoms):
                if i == im:
                    continue
                r_a = _rcov(a["elem"])
                if _dist(m, a) <= scale * (r_m + r_a):
                    if i not in selected:
                        next_frontier.add(i)
        selected.update(next_frontier)
        frontier = next_frontier
        if not frontier:
            break

    # convert to 1-based indices for ORCA, optionally drop metals
    if not include_metals:
        selected = {i for i in selected if i not in metal_indices}
    if exclude_h:
        selected = {i for i in selected if atoms[i]["elem"] != "H"}
    uphill = sorted(i + 1 for i in selected)
    return uphill


def extract_optimized_xyz(orca_out: Path, output_xyz: Path) -> bool:
    """Extract final optimized geometry from ORCA output.

    The XTB DOCKER method writes the docked structure. We look for
    the final coordinates in the output or the generated .xyz file.

    Args:
        orca_out: Path to ORCA output file
        output_xyz: Path to write extracted XYZ

    Returns:
        True if extraction was successful
    """
    # ORCA XTB DOCKER creates a .xyz file with the result
    # The output basename + .xyz should contain the docked structure
    base = orca_out.stem
    parent = orca_out.parent

    # Try various possible output files
    possible_xyz = [
        parent / f"{base}.docker.xyz",  # XTB DOCKER output
        parent / f"{base}.xyz",
        parent / f"{base}_trj.xyz",
        parent / f"{base}_dock.xyz",
    ]

    for xyz_file in possible_xyz:
        if xyz_file.exists():
            shutil.copy(xyz_file, output_xyz)
            return True

    # If no .xyz file found, try to parse from output
    try:
        content = orca_out.read_text(encoding='utf-8', errors='ignore')

        # Look for CARTESIAN COORDINATES (ANGSTROEM) section
        pattern = r'CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)\n\n'
        matches = re.findall(pattern, content, re.DOTALL)

        if matches:
            # Take the last match (final geometry)
            coord_block = matches[-1].strip()
            lines = [l.strip() for l in coord_block.split('\n') if l.strip()]
            atom_count = len(lines)

            xyz_content = f"{atom_count}\nExtracted from {orca_out.name}\n"
            xyz_content += '\n'.join(lines) + '\n'
            output_xyz.write_text(xyz_content, encoding='utf-8')
            return True

    except Exception as e:
        logger.warning(f"Could not extract geometry from {orca_out}: {e}")

    return False


def _run_docker_step(
    builder_dir: Path,
    step: int,
    host_xyz: str,
    guest_xyz: str,
    host_charge: int,
    guest_charge: int,
    multiplicity: int,
    description: str,
    dry_run: bool,
    pal: int,
    maxcore: int,
    use_goat: bool = False,
    uphill_include_metals: bool = False,
    step_goat_preoptimized: bool = False,
    step_goat_swarm: bool = False,
) -> Tuple[bool, str, int]:
    """Run a single DOCKER step, optionally followed by GOAT optimization.

    Args:
        use_goat: If True, run GOAT global optimization after docking

    Returns:
        Tuple of (success, output_xyz_name, new_charge)
    """
    step_input = builder_dir / f"step_{step}.inp"
    step_output = builder_dir / f"step_{step}.out"
    step_xyz = builder_dir / f"step_{step}.xyz"

    logger.info(f"\n=== Step {step}: {description} ===")

    # Create ORCA input
    create_docker_input(
        host_xyz=host_xyz,
        guest_xyz=guest_xyz,
        guest_charge=guest_charge,
        host_charge=host_charge,
        multiplicity=multiplicity,
        output_path=step_input,
        pal=pal,
        maxcore=maxcore,
    )
    logger.info(f"Created {step_input.name}")
    logger.info(f"  Host: {host_xyz} (charge: {host_charge:+d})")
    logger.info(f"  Guest: {guest_xyz} (charge: {guest_charge:+d})")

    new_charge = host_charge + guest_charge

    if dry_run:
        logger.info(f"[DRY RUN] Would run ORCA for step {step}")
        if use_goat:
            logger.info(f"[DRY RUN] Would run GOAT optimization after docking")
        shutil.copy(builder_dir / host_xyz, step_xyz)
    else:
        # Recalc logic: skip ORCA if previous output terminated normally
        if _orca_terminated_normally(step_output):
            logger.info(f"Skipping ORCA for step {step}: found normal termination in {step_output.name}")
            if not step_xyz.exists():
                if not extract_optimized_xyz(step_output, step_xyz):
                    logger.warning(f"Could not extract geometry from {step_output.name}; rerunning ORCA")
                else:
                    logger.info(f"Extracted geometry to {step_xyz.name}")
        if not step_xyz.exists():
            logger.info(f"Running ORCA for step {step}...")
            success = run_orca(
                str(step_input),
                str(step_output),
                working_dir=builder_dir,
            )

            if not success:
                logger.error(f"ORCA failed for step {step}")
                return False, "", 0

            if not extract_optimized_xyz(step_output, step_xyz):
                logger.error(f"Could not extract geometry from step {step}")
                return False, "", 0

            logger.info(f"Step {step} docking completed successfully")

        # Run GOAT optimization if requested
        if use_goat:
            goat_result = builder_dir / f"step_{step}_goat.globalminimum.xyz"
            if goat_result.exists():
                logger.info(f"Skipping GOAT for step {step}: found {goat_result.name}")
            else:
                goat_xyz = f"step_{step}.xyz"
                if step_goat_preoptimized or step_goat_swarm:
                    if step_goat_swarm:
                        traj_name = f"step_{step}.docker.struc1.all.swarm.xyz"
                        out_name = f"step_{step}.swarm_best.xyz"
                        src_label = "swarm"
                    else:
                        traj_name = f"step_{step}.docker.struc1.all.preoptimized.xyz"
                        out_name = f"step_{step}.preopt_first.xyz"
                        src_label = "preoptimized"
                    traj_path = builder_dir / traj_name
                    out_path = builder_dir / out_name
                    if _extract_lowest_energy_xyz_from_trajectory(traj_path, out_path):
                        goat_xyz = out_path.name
                        logger.info(f"Using lowest-energy {src_label} structure for GOAT: {traj_path.name}")
                    else:
                        logger.warning(
                            f"{src_label.capitalize()} trajectory not found or invalid: {traj_path.name}; "
                            f"falling back to {goat_xyz}"
                        )
                goat_success, _ = run_goat_optimization(
                    builder_dir=builder_dir,
                    xyz_file=goat_xyz,
                    charge=new_charge,
                    multiplicity=multiplicity,
                    pal=pal,
                    maxcore=maxcore,
                    prefix=f"step_{step}_goat",
                    uphill_from_xyz=True,
                    uphill_include_metals=uphill_include_metals,
                    restrictive=step_goat_preoptimized or step_goat_swarm,
                    skip_initial_opt=step_goat_preoptimized or step_goat_swarm,
                )
                if not goat_success:
                    logger.warning(f"GOAT optimization failed for step {step}, continuing with docked structure")

    return True, f"step_{step}.xyz", new_charge


def run_build_up(
    smiles: str,
    builder_dir: Path,
    multiplicity: int = 1,
    dry_run: bool = False,
    pal: int = 32,
    maxcore: int = 1000,
    use_goat: bool = False,
    skip_ligand_goat: bool = False,
    uphill_include_metals: bool = False,
    step_goat_preoptimized: bool = False,
    step_goat_swarm: bool = False,
    cmdline: str | None = None,
) -> bool:
    """Run the full complex build-up workflow.

    Workflow depends on number of metals:
    - Mono-metal: Metal as host, ligands docked by size (largest first)
    - Multi-metal: Largest ligand as host, then metals, then remaining ligands

    Args:
        smiles: SMILES string of the metal complex
        builder_dir: Directory to create for the build-up process
        multiplicity: Spin multiplicity (default: 1)
        dry_run: If True, only create input files but don't run ORCA
        pal: Number of parallel processes (default: 32)
        maxcore: Memory per core in MB (default: 1000)
        use_goat: If True, run GOAT global optimization after each docking step
        skip_ligand_goat: If True, skip GOAT for initial ligand structures
        uphill_include_metals: If True, include metal atoms in UPHILLATOMS
        step_goat_preoptimized: If True, use lowest-energy preoptimized trajectory structure for GOAT steps
        step_goat_swarm: If True, use lowest-energy swarm trajectory structure for GOAT steps

    Returns:
        True if successful
    """
    # Check prerequisites
    if not RDKIT_AVAILABLE:
        logger.error("RDKit is required for SMILES parsing")
        return False

    if not dry_run:
        orca_path = find_orca_executable()
        if not orca_path:
            logger.error("ORCA executable not found")
            return False

    # Extract metals and ligands
    logger.info(f"Parsing complex SMILES: {smiles}")
    metals, ligands, error = extract_ligands_from_complex(smiles)

    if error:
        logger.error(f"Failed to extract ligands: {error}")
        return False

    if not metals:
        logger.error("No metals found in complex")
        return False

    if not ligands:
        logger.error("No ligands found in complex")
        return False

    # Log metals
    for i, (metal_symbol, metal_charge) in enumerate(metals, 1):
        logger.info(f"Found metal {i}: {metal_symbol} (charge: {metal_charge:+d})")
    logger.info(f"Found {len(ligands)} ligand(s)")

    # Create builder directory
    builder_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created builder directory: {builder_dir}")
    _write_build_info(builder_dir, cmdline)

    # Calculate ligand charges and sizes. Sorting prefers:
    # 1) Non-halide/pseudohalide ligands
    # 2) Very large ligands first if there's a strong size disparity
    # 3) Charged first, neutral last; more negative first
    # 4) Within same charge, largest first
    ligand_info = []
    for lig_smiles in ligands:
        charge = get_ligand_charge(lig_smiles)
        size = get_ligand_size(lig_smiles)
        ligand_info.append({
            'smiles': lig_smiles,
            'charge': charge,
            'size': size,
        })

    # Detect strong size disparity -> prioritize very large ligands first
    sizes = sorted(info['size'] for info in ligand_info)
    if sizes:
        mid = len(sizes) // 2
        if len(sizes) % 2 == 0:
            median_size = 0.5 * (sizes[mid - 1] + sizes[mid])
        else:
            median_size = float(sizes[mid])
    else:
        median_size = 0.0
    large_cutoff = median_size * 1.5 if median_size > 0 else 0.0
    for info in ligand_info:
        info['is_large'] = info['size'] >= large_cutoff

    ligand_info.sort(
        key=lambda x: (
            _is_halide_or_pseudohalide(x['smiles']),
            not x.get('is_large', False),
            x['charge'] == 0,
            x['charge'],
            -x['size'],
        )
    )

    for i, info in enumerate(ligand_info, 1):
        logger.info(f"  Ligand {i}: {info['smiles']} (charge: {info['charge']:+d}, size: {info['size']} atoms)")

    # Generate ligand XYZ files (numbered by sorted order)
    # Track already optimized SMILES to avoid duplicate GOAT runs
    goat_cache: Dict[str, str] = {}  # smiles -> optimized xyz filename
    run_ligand_goat = use_goat and not skip_ligand_goat

    for i, info in enumerate(ligand_info, 1):
        lig_xyz = builder_dir / f"ligand{i}.xyz"
        info['xyz'] = f"ligand{i}.xyz"
        smiles_key = info['smiles']

        # Check if this SMILES was already GOAT-optimized
        if run_ligand_goat and smiles_key in goat_cache:
            # Copy from already optimized ligand
            src_xyz = builder_dir / goat_cache[smiles_key]
            shutil.copy(src_xyz, lig_xyz)
            logger.info(f"Created {lig_xyz.name} (copied from {goat_cache[smiles_key]}, same SMILES)")
        else:
            # Generate new XYZ
            success, err = ligand_to_xyz(info['smiles'], lig_xyz)
            if not success:
                logger.error(f"Failed to convert ligand {i} to XYZ: {err}")
                return False
            logger.info(f"Created {lig_xyz.name}")

            # Optionally run GOAT on ligands (only for first occurrence of each SMILES)
            if run_ligand_goat and not dry_run:
                goat_success, _ = run_goat_optimization(
                    builder_dir=builder_dir,
                    xyz_file=info['xyz'],
                    charge=info['charge'],
                    multiplicity=1,  # Ligands typically singlet
                    pal=pal,
                    maxcore=maxcore,
                    prefix=f"ligand{i}_goat",
                    include_goat_block=False,
                )
                if goat_success:
                    logger.info(f"GOAT optimized {info['xyz']}")
                    goat_cache[smiles_key] = info['xyz']  # Cache for duplicates
                else:
                    logger.warning(f"GOAT failed for {info['xyz']}, using RDKit geometry")

    # Create metal XYZ files
    for i, (metal_symbol, metal_charge) in enumerate(metals, 1):
        metal_xyz = builder_dir / f"metal{i}.xyz"
        create_metal_xyz(metal_symbol, metal_xyz)
        logger.info(f"Created metal{i}.xyz ({metal_symbol})")

    # === Build-up workflow ===
    step = 0
    is_multi_metal = len(metals) > 1

    if is_multi_metal:
        # Multi-metal: Prefer grouping identical metal-ligand fragments to avoid redundant docking
        logger.info("\n=== Multi-metal workflow ===")

        unique_metals = {(m[0], m[1]) for m in metals}
        single_metal_type = len(unique_metals) == 1

        # Build SMILES -> entries mapping (order preserved from ligand_info)
        smiles_to_entries: Dict[str, List[Dict]] = {}
        smiles_order: List[str] = []
        for info in ligand_info:
            smi = info['smiles']
            if smi not in smiles_to_entries:
                smiles_to_entries[smi] = []
                smiles_order.append(smi)
            smiles_to_entries[smi].append(info)

        if single_metal_type:
            n_metals = len(metals)
            per_metal: Dict[str, int] = {
                smi: len(entries) // n_metals for smi, entries in smiles_to_entries.items()
            }
            has_groupable = any(v > 0 for v in per_metal.values())
        else:
            has_groupable = False

        if single_metal_type and has_groupable:
            # Build one representative metal-ligand unit, then reuse it for identical metals
            metal_symbol, metal_charge = metals[0]
            logger.info(
                f"Grouping identical metals: {n_metals} x {metal_symbol} (charge: {metal_charge:+d})"
            )

            used_idx: Dict[str, int] = {smi: 0 for smi in smiles_to_entries}

            # 1) Build base unit: M + (per_metal ligands for each identical SMILES)
            current_host = "metal1.xyz"
            current_charge = metal_charge
            logger.info("Building base unit from shared ligands")
            for smi in smiles_order:
                repeat = per_metal.get(smi, 0)
                for _ in range(repeat):
                    info = smiles_to_entries[smi][used_idx[smi]]
                    used_idx[smi] += 1
                    step += 1
                    success, new_host, current_charge = _run_docker_step(
                        builder_dir=builder_dir,
                        step=step,
                        host_xyz=current_host,
                        guest_xyz=info['xyz'],
                        host_charge=current_charge,
                        guest_charge=info['charge'],
                        multiplicity=multiplicity,
                        description=f"Docking ligand for unit ({smi})",
                        dry_run=dry_run,
                        pal=pal,
                        maxcore=maxcore,
                        use_goat=use_goat,
                        uphill_include_metals=uphill_include_metals,
                        step_goat_preoptimized=step_goat_preoptimized,
                        step_goat_swarm=step_goat_swarm,
                    )
                    if not success:
                        return False
                    current_host = new_host

            unit_xyz = current_host
            unit_charge = current_charge
            logger.info(f"Base unit ready: {unit_xyz} (charge: {unit_charge:+d})")

            # 2) Dock remaining ligands (those not distributed equally across metals)
            logger.info("Docking remaining ligands onto base unit")
            for smi in smiles_order:
                total = len(smiles_to_entries[smi])
                remaining = total - per_metal.get(smi, 0) * n_metals
                for _ in range(remaining):
                    info = smiles_to_entries[smi][used_idx[smi]]
                    used_idx[smi] += 1
                    step += 1
                    success, new_host, current_charge = _run_docker_step(
                        builder_dir=builder_dir,
                        step=step,
                        host_xyz=current_host,
                        guest_xyz=info['xyz'],
                        host_charge=current_charge,
                        guest_charge=info['charge'],
                        multiplicity=multiplicity,
                        description=f"Docking remaining ligand ({smi})",
                        dry_run=dry_run,
                        pal=pal,
                        maxcore=maxcore,
                        use_goat=use_goat,
                        uphill_include_metals=uphill_include_metals,
                        step_goat_preoptimized=step_goat_preoptimized,
                        step_goat_swarm=step_goat_swarm,
                    )
                    if not success:
                        return False
                    current_host = new_host

            # 3) Attach the remaining identical units
            logger.info("Docking remaining identical metal-ligand units")
            for i in range(2, n_metals + 1):
                unit_copy = f"unit_{i}.xyz"
                shutil.copy(builder_dir / unit_xyz, builder_dir / unit_copy)
                step += 1
                success, new_host, current_charge = _run_docker_step(
                    builder_dir=builder_dir,
                    step=step,
                    host_xyz=current_host,
                    guest_xyz=unit_copy,
                    host_charge=current_charge,
                    guest_charge=unit_charge,
                    multiplicity=multiplicity,
                    description=f"Docking unit {i} (identical fragment)",
                    dry_run=dry_run,
                    pal=pal,
                    maxcore=maxcore,
                    use_goat=use_goat,
                    uphill_include_metals=uphill_include_metals,
                    step_goat_preoptimized=step_goat_preoptimized,
                    step_goat_swarm=step_goat_swarm,
                )
                if not success:
                    return False
                current_host = new_host

        else:
            # Fallback: First ligand by charge/size as host, then metals, then remaining ligands
            logger.info("Fallback multi-metal workflow (no grouping possible)")

            # Start with first ligand in sorted order as host
            first_ligand = ligand_info[0]
            current_host = first_ligand['xyz']
            current_charge = first_ligand['charge']
            logger.info(
                f"Starting with first ligand (charge-priority): {first_ligand['smiles']} "
                f"(charge: {first_ligand['charge']:+d}, size: {first_ligand['size']} atoms)"
            )

            # Dock metals onto the growing complex
            for i, (metal_symbol, metal_charge) in enumerate(metals, 1):
                step += 1
                success, new_host, current_charge = _run_docker_step(
                    builder_dir=builder_dir,
                    step=step,
                    host_xyz=current_host,
                    guest_xyz=f"metal{i}.xyz",
                    host_charge=current_charge,
                    guest_charge=metal_charge,
                    multiplicity=multiplicity,
                    description=f"Docking metal {i} ({metal_symbol})",
                    dry_run=dry_run,
                    pal=pal,
                    maxcore=maxcore,
                    use_goat=use_goat,
                    uphill_include_metals=uphill_include_metals,
                    step_goat_preoptimized=step_goat_preoptimized,
                    step_goat_swarm=step_goat_swarm,
                )
                if not success:
                    return False
                current_host = new_host

            # Dock remaining ligands (skip first, already used as host)
            for i, info in enumerate(ligand_info[1:], 2):
                step += 1
                success, new_host, current_charge = _run_docker_step(
                    builder_dir=builder_dir,
                    step=step,
                    host_xyz=current_host,
                    guest_xyz=info['xyz'],
                    host_charge=current_charge,
                    guest_charge=info['charge'],
                    multiplicity=multiplicity,
                    description=f"Docking ligand {i} ({info['size']} atoms)",
                    dry_run=dry_run,
                    pal=pal,
                    maxcore=maxcore,
                    use_goat=use_goat,
                    uphill_include_metals=uphill_include_metals,
                    step_goat_preoptimized=step_goat_preoptimized,
                    step_goat_swarm=step_goat_swarm,
                )
                if not success:
                    return False
                current_host = new_host

    else:
        # Mono-metal: Metal as host, ligands docked by charge/size order
        logger.info("\n=== Mono-metal workflow ===")

        metal_symbol, metal_charge = metals[0]
        current_host = "metal1.xyz"
        current_charge = metal_charge

        # Dock ligands onto metal (already sorted by charge/size)
        for i, info in enumerate(ligand_info, 1):
            step += 1
            success, new_host, current_charge = _run_docker_step(
                builder_dir=builder_dir,
                step=step,
                host_xyz=current_host,
                guest_xyz=info['xyz'],
                host_charge=current_charge,
                guest_charge=info['charge'],
                multiplicity=multiplicity,
                description=f"Docking ligand {i} ({info['size']} atoms)",
                dry_run=dry_run,
                pal=pal,
                maxcore=maxcore,
                use_goat=use_goat,
                uphill_include_metals=uphill_include_metals,
                step_goat_preoptimized=step_goat_preoptimized,
                step_goat_swarm=step_goat_swarm,
            )
            if not success:
                return False
            current_host = new_host

    logger.info(f"\n=== Build-up complete ===")
    logger.info(f"Final structure: {builder_dir / current_host}")
    logger.info(f"Final charge: {current_charge:+d}")

    # Export final structure to parent directory
    try:
        parent_dir = builder_dir.parent
        if use_goat and not dry_run and step > 0:
            final_xyz = builder_dir / f"step_{step}_goat.globalminimum.xyz"
        else:
            final_xyz = builder_dir / current_host

        if final_xyz.exists():
            out_xyz = parent_dir / "build_complex.xyz"
            shutil.copy(final_xyz, out_xyz)

            # Write start.txt with XYZ content excluding the first two lines
            lines = final_xyz.read_text().splitlines()
            coord_lines = lines[2:] if len(lines) >= 2 else []
            (parent_dir / "start.txt").write_text("\n".join(coord_lines).strip() + ("\n" if coord_lines else ""))

            logger.info(f"Exported final structure to {out_xyz}")
            logger.info(f"Exported coordinates to {parent_dir / 'start.txt'}")
        else:
            logger.warning(f"Final structure not found: {final_xyz}")
    except Exception as exc:
        logger.warning(f"Failed to export final structure: {exc}")

    return True


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point for build_up_complex.

    Args:
        argv: Command line arguments (uses sys.argv if None)

    Returns:
        Exit code (0 for success, non-zero for failure)
    """
    parser = argparse.ArgumentParser(
        prog="delfin-build",
        description="Build up metal complex step by step using ORCA XTB DOCKER",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default="input.txt",
        help="Input file containing SMILES string (default: input.txt)",
    )
    parser.add_argument(
        "-d", "--directory",
        default="builder",
        help="Output directory for build-up files (default: builder)",
    )
    parser.add_argument(
        "-m", "--multiplicity",
        type=int,
        default=1,
        help="Spin multiplicity (default: 1)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only create input files, don't run ORCA",
    )
    parser.add_argument(
        "-p", "--pal",
        type=int,
        default=32,
        help="Number of parallel processes (default: 32)",
    )
    parser.add_argument(
        "--maxcore",
        type=int,
        default=1000,
        help="Memory per core in MB (default: 1000)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output",
    )
    parser.add_argument(
        "--goat",
        action="store_true",
        help="Run GOAT global optimization on ligands and after each docking step",
    )
    parser.add_argument(
        "--no-ligand-goat",
        action="store_true",
        help="Skip GOAT optimization for initial ligand structures (use with --goat)",
    )
    parser.add_argument(
        "--uphill-include-metals",
        action="store_true",
        help="Include metal atoms in UPHILLATOMS for GOAT (default: false)",
    )
    parser.add_argument(
        "--step-goat-preoptimized",
        action="store_true",
        help=(
            "Use lowest-energy structure from preoptimized trajectory "
            "(step_N.docker.struc1.all.preoptimized.xyz) for GOAT steps and run "
            "restrictive GOAT with SKIPINITIALOPT"
        ),
    )
    parser.add_argument(
        "--step-goat-swarm",
        action="store_true",
        help=(
            "Use lowest-energy structure from swarm trajectory "
            "(step_N.docker.struc1.all.swarm.xyz) for GOAT steps and run "
            "restrictive GOAT with SKIPINITIALOPT"
        ),
    )

    args = parser.parse_args(argv)

    # Configure logging
    configure_logging(level=logging.DEBUG if args.verbose else logging.INFO)

    # Read SMILES from input file
    input_path = Path(args.input_file)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        return 1

    try:
        content = input_path.read_text(encoding='utf-8').strip()
    except Exception as e:
        logger.error(f"Could not read input file: {e}")
        return 1

    # Get first non-empty, non-comment line as SMILES
    smiles = None
    for line in content.split('\n'):
        line = line.strip()
        if line and not line.startswith('#'):
            smiles = line
            break

    if not smiles:
        logger.error("No SMILES string found in input file")
        return 1

    if not contains_metal(smiles):
        logger.error("SMILES does not appear to contain a metal complex")
        return 1

    # Run build-up
    builder_dir = Path(args.directory)
    success = run_build_up(
        smiles=smiles,
        builder_dir=builder_dir,
        multiplicity=args.multiplicity,
        dry_run=args.dry_run,
        pal=args.pal,
        maxcore=args.maxcore,
        use_goat=args.goat,
        skip_ligand_goat=args.no_ligand_goat,
        uphill_include_metals=args.uphill_include_metals,
        step_goat_preoptimized=args.step_goat_preoptimized,
        step_goat_swarm=args.step_goat_swarm,
        cmdline=" ".join(sys.argv),
    )

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
