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
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import Optional, Tuple, List, Dict

from delfin.common.logging import configure_logging, get_logger
from delfin.orca import run_orca, find_orca_executable
from delfin.smiles_converter import RDKIT_AVAILABLE, contains_metal, _METALS

if RDKIT_AVAILABLE:
    from rdkit import Chem
    from rdkit.Chem import AllChem

logger = get_logger(__name__)


def extract_ligands_from_complex(smiles: str) -> Tuple[
    Optional[str],  # metal_symbol
    Optional[int],  # metal_charge
    Optional[List[str]],  # ligand_smiles_list
    Optional[str],  # error_message
]:
    """Extract metal and ligands from a metal complex SMILES.

    Based on ChemDarwin's extract_ligands_from_complex function.

    Args:
        smiles: SMILES string of the metal complex

    Returns:
        Tuple of (metal_symbol, metal_charge, ligand_smiles_list, error_message)
        - metal_symbol: Element symbol of the metal (e.g., "Mn")
        - metal_charge: Formal charge of the metal
        - ligand_smiles_list: List of SMILES strings for each ligand
        - error_message: Error description if failed, None on success
    """
    if not RDKIT_AVAILABLE:
        return None, None, None, "RDKit not available"

    # Parse without sanitization for metal complexes
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None, None, None, "Could not parse SMILES"

    # Find metal atom(s)
    metal_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in _METALS:
            metal_atoms.append(atom.GetIdx())

    if not metal_atoms:
        return None, None, None, "No metal found in SMILES"

    # For now, handle single metal center
    metal_idx = metal_atoms[0]
    metal_atom = mol.GetAtomWithIdx(metal_idx)
    metal_symbol = metal_atom.GetSymbol()
    metal_charge = metal_atom.GetFormalCharge()

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

    # Remove bonds to metal first
    bonds_to_remove = []
    for bond in edit_mol.GetBonds():
        if bond.GetBeginAtomIdx() == metal_idx or bond.GetEndAtomIdx() == metal_idx:
            bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    for b in bonds_to_remove:
        edit_mol.RemoveBond(b[0], b[1])

    # Remove metal atom (this shifts indices!)
    edit_mol.RemoveAtom(metal_idx)

    # Get fragments (ligands)
    try:
        mol_no_metal = edit_mol.GetMol()

        # Apply original atom properties
        for old_idx, props in orig_props.items():
            if old_idx == metal_idx:
                continue
            new_idx = old_idx - 1 if old_idx > metal_idx else old_idx
            try:
                a = mol_no_metal.GetAtomWithIdx(new_idx)
                a.SetFormalCharge(props['formal_charge'])
                a.SetNumExplicitHs(props['explicit_h'])
                a.SetNoImplicit(props['no_implicit'])
                a.SetNumRadicalElectrons(props['rad_e'])
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

        return metal_symbol, metal_charge, ligand_smiles, None

    except Exception as e:
        return None, None, None, str(e)


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


def ligand_to_xyz(smiles: str, output_path: Path) -> Tuple[bool, Optional[str]]:
    """Convert a ligand SMILES to XYZ file.

    Args:
        smiles: SMILES string of the ligand
        output_path: Path to write the XYZ file

    Returns:
        Tuple of (success, error_message)
    """
    if not RDKIT_AVAILABLE:
        return False, "RDKit not available"

    try:
        # Try normal parsing first
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Try without sanitization
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                try:
                    mol.UpdatePropertyCache(strict=False)
                except Exception:
                    pass

        if mol is None:
            return False, f"Could not parse SMILES: {smiles}"

        # Add hydrogens
        try:
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
            pass

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
    maxcore: int = 4000,
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
        maxcore: Memory per core in MB (default: 4000)
    """
    content = f"""! XTB2

%PAL NPROCS {pal} END
%MAXCORE {maxcore}

%DOCKER
    GUEST "{guest_xyz}"
    GuestCharge    {guest_charge}
END

* XYZFile {host_charge} {multiplicity} {host_xyz}
"""
    output_path.write_text(content, encoding='utf-8')


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


def run_build_up(
    smiles: str,
    builder_dir: Path,
    multiplicity: int = 1,
    dry_run: bool = False,
    pal: int = 32,
    maxcore: int = 4000,
) -> bool:
    """Run the full complex build-up workflow.

    Args:
        smiles: SMILES string of the metal complex
        builder_dir: Directory to create for the build-up process
        multiplicity: Spin multiplicity (default: 1)
        dry_run: If True, only create input files but don't run ORCA
        pal: Number of parallel processes (default: 32)
        maxcore: Memory per core in MB (default: 4000)

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

    # Extract metal and ligands
    logger.info(f"Parsing complex SMILES: {smiles}")
    metal, metal_charge, ligands, error = extract_ligands_from_complex(smiles)

    if error:
        logger.error(f"Failed to extract ligands: {error}")
        return False

    if not ligands:
        logger.error("No ligands found in complex")
        return False

    logger.info(f"Found metal: {metal} (charge: {metal_charge:+d})")
    logger.info(f"Found {len(ligands)} ligand(s)")

    # Create builder directory
    builder_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created builder directory: {builder_dir}")

    # Calculate ligand charges
    ligand_charges = []
    for i, lig_smiles in enumerate(ligands, 1):
        charge = get_ligand_charge(lig_smiles)
        ligand_charges.append(charge)
        logger.info(f"  Ligand {i}: {lig_smiles} (charge: {charge:+d})")

    # Generate ligand XYZ files
    for i, lig_smiles in enumerate(ligands, 1):
        lig_xyz = builder_dir / f"ligand{i}.xyz"
        success, error = ligand_to_xyz(lig_smiles, lig_xyz)
        if not success:
            logger.error(f"Failed to convert ligand {i} to XYZ: {error}")
            return False
        logger.info(f"Created {lig_xyz.name}")

    # Create initial metal structure
    metal_xyz = builder_dir / "metal.xyz"
    create_metal_xyz(metal, metal_xyz)
    logger.info(f"Created {metal_xyz.name}")

    # Build up complex step by step
    current_host = "metal.xyz"
    current_charge = metal_charge if metal_charge else 0

    for step, (lig_smiles, lig_charge) in enumerate(zip(ligands, ligand_charges), 1):
        logger.info(f"\n=== Step {step}: Adding ligand {step} ===")

        guest_xyz = f"ligand{step}.xyz"
        step_input = builder_dir / f"step_{step}.inp"
        step_output = builder_dir / f"step_{step}.out"
        step_xyz = builder_dir / f"step_{step}.xyz"

        # Create ORCA input
        create_docker_input(
            host_xyz=current_host,
            guest_xyz=guest_xyz,
            guest_charge=lig_charge,
            host_charge=current_charge,
            multiplicity=multiplicity,
            output_path=step_input,
            pal=pal,
            maxcore=maxcore,
        )
        logger.info(f"Created {step_input.name}")
        logger.info(f"  Host: {current_host} (charge: {current_charge:+d})")
        logger.info(f"  Guest: {guest_xyz} (charge: {lig_charge:+d})")

        if dry_run:
            logger.info(f"[DRY RUN] Would run ORCA for step {step}")
            # Create a dummy output for the next step
            shutil.copy(builder_dir / current_host, step_xyz)
        else:
            # Run ORCA
            logger.info(f"Running ORCA for step {step}...")
            success = run_orca(
                str(step_input),
                str(step_output),
                working_dir=builder_dir,
            )

            if not success:
                logger.error(f"ORCA failed for step {step}")
                return False

            # Extract optimized geometry
            if not extract_optimized_xyz(step_output, step_xyz):
                logger.error(f"Could not extract geometry from step {step}")
                return False

            logger.info(f"Step {step} completed successfully")

        # Update for next step
        current_host = f"step_{step}.xyz"
        current_charge += lig_charge

    logger.info(f"\n=== Build-up complete ===")
    logger.info(f"Final structure: {builder_dir / current_host}")
    logger.info(f"Final charge: {current_charge:+d}")

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
        default=4000,
        help="Memory per core in MB (default: 4000)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output",
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
    )

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
