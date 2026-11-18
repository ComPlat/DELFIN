#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Proton coordinator - donor templates, placement, orientation + relaxed scans.

Analogue of CO2_Coordinator6.py for the competing hydrogen-evolution
reaction (HER). A proton donor MOLECULE (never a bare proton) is placed
next to the metal centre with its donor H pointing at the metal. An
orientation scan at fixed M...H distance selects the best rotation, then a
relaxed M...H distance scan maps the proton/hydrogen-transfer path.

Donor chemistry contract (chain_setup.bookkeeping):
- protic donors (TFE, HOAc, phenol, H3O(H2O)n): proton transfer,
  products [M-H] + conjugate base.
- HAT donors (BIH): hydrogen atom transfer, products [M-H] + BI radical.

Folder naming follows the donor: each donor gets its own
``<DonorName>_coordinator/`` directory (H3O+_coordinator, TFE_coordinator,
BIH_coordinator, phenol_coordinator, HOAc_coordinator), driven from a
single CONTROL.txt via ``proton_donor=``.
"""

import os
# --- Headless plot backend, before matplotlib is imported ---
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")

import shutil
import subprocess
import numpy as np
from typing import Dict, List, Optional, Sequence, Tuple

from ase import Atoms
from ase.io import read, write

import matplotlib.pyplot as plt

# === Donor templates =========================================================
# Geometry: ASE-readable xyz text. donor_h_index is 0-based within the donor.
# product: charge/multiplicity of the DONOR-derived transfer product
# (conjugate base / radical) AFTER the H (or H+) has left.
# The combined complex+donor charge/multiplicity are computed by
# chain_setup.bookkeeping().
DONOR_TEMPLATES: Dict[str, Dict] = {}

DONOR_TEMPLATES["tfe"] = {
    # Trifluoroethanol: F3C-CH2-O-H; donor H = hydroxyl H (index 6).
    "xyz": """7

C      0.000000    0.000000    0.000000
C      1.250000    0.350000    0.950000
F      0.050000    1.250000   -0.150000
F     -0.950000   -0.650000   -0.550000
F      0.950000   -0.950000   -0.650000
O      2.300000    1.050000    0.150000
H      2.250000    1.900000   -0.300000
""",
    "n_atoms": 7,
    "donor_h_index": 6,
    "donor_heavy_index": 5,
    "transfer": "protic",          # PT: [M-H] + TFE-O(-)
    "product": {"name": "tfe_alkoxide", "charge": -1, "multiplicity": 1},
}

DONOR_TEMPLATES["bih"] = {
    # 1,3-dimethyl-2-phenyl-benzimidazoline (BIH); C2-H is the donor H.
    # Connectivity-faithful sketch geometry; refinement by the scan itself.
    "xyz": """17

N      0.000000    1.390000    0.000000
N      0.000000   -1.390000    0.000000
C      0.000000    0.000000    0.300000
C      1.150000    0.680000    0.900000
C      1.150000   -0.680000    0.900000
C     -1.150000    0.680000    0.900000
C     -1.150000   -0.680000    0.900000
C      2.320000    0.000000    1.550000
C     -2.320000    0.000000    1.550000
C      2.900000    0.850000    2.300000
C      2.900000   -0.850000    2.300000
C     -2.900000    0.850000    2.300000
C     -2.900000   -0.850000    2.300000
N      0.000000    2.100000    0.950000
N      0.000000   -2.100000    0.950000
H      0.000000    0.000000    1.350000
H      0.000000    0.000000   -0.650000
""",
    "n_atoms": 17,
    "donor_h_index": 15,           # C2-H
    "donor_heavy_index": 2,        # C2
    "transfer": "hat",             # HAT: [M-H] + BI radical (open shell)
    "product": {"name": "bi_radical", "charge": 0, "multiplicity": 2},
}

DONOR_TEMPLATES["phenol"] = {
    # Phenol C6H5-O-H; donor H = hydroxyl H (last atom).
    "xyz": """14

C      0.000000    1.200000    0.000000
C      1.040000    0.600000    0.000000
C      1.040000   -0.600000    0.000000
C      0.000000   -1.200000    0.000000
C     -1.040000   -0.600000    0.000000
C     -1.040000    0.600000    0.000000
H      0.000000    2.140000    0.000000
H      1.860000    1.100000    0.000000
H      1.860000   -1.100000    0.000000
H      0.000000   -2.140000    0.000000
H     -1.860000   -1.100000    0.000000
H     -1.860000    1.100000    0.000000
O      0.000000    0.000000    0.000000
H      0.000000    0.000000    0.980000
""",
    "n_atoms": 14,
    "donor_h_index": 13,
    "donor_heavy_index": 12,
    "transfer": "hat",             # HAT: [M-H] + PhO radical (open shell)
    "product": {"name": "phenoxyl_radical", "charge": 0, "multiplicity": 2},
    # PCET alternative (proton transfer, closed-shell phenolate) selectable
    # via proton_transfer_model=protic; the per-model product lives here:
    "protic_product": {"name": "phenolate", "charge": -1, "multiplicity": 1},
}

DONOR_TEMPLATES["hoac"] = {
    # Acetic acid CH3-C(=O)-O-H; donor H = hydroxyl H (index 3).
    "xyz": """8

C      0.000000    0.000000    0.000000
C      1.350000    0.000000    0.000000
O      1.900000    1.050000    0.000000
H      1.850000    1.850000    0.000000
O      1.950000   -1.050000    0.000000
C      1.900000    0.000000   -1.550000
H      3.000000    0.000000   -1.550000
H      1.500000    0.900000   -2.100000
""",
    "n_atoms": 8,
    "donor_h_index": 3,
    "donor_heavy_index": 2,
    "transfer": "protic",          # PT: [M-H] + OAc(-)
    "product": {"name": "acetate", "charge": -1, "multiplicity": 1},
}

DONOR_TEMPLATES["h3o"] = {
    # Eigen-like cluster H3O+(H2O)3 (n=3 water shells as default; the
    # coordinates sketch a tetrahedral first shell). Donor H = one H of the
    # hydronium (index 4). Never a bare H3O+ in vacuum: the explicit water
    # shell stays part of the donor and ends up in the product fragments.
    "xyz": """13

O      0.000000    0.000000    0.000000
H      0.960000    0.260000    0.000000
H     -0.260000    0.260000    0.960000
H     -0.260000    0.260000   -0.960000
H      0.000000    0.000000    -0.980000
O      1.600000    1.600000    1.600000
H      2.300000    1.000000    1.200000
H      1.000000    2.300000    1.200000
O     -1.600000   -1.600000    1.600000
H     -2.300000   -1.000000    1.200000
H     -1.000000   -2.300000    1.200000
O     -1.600000    1.600000   -1.600000
H     -2.300000    1.000000   -1.200000
""",
    "n_atoms": 13,
    "donor_charge": 1,             # H3O+ itself carries +1
    "donor_h_index": 4,
    "donor_heavy_index": 0,
    "transfer": "protic",          # PT: [M-H] + H2O(H2O)3 residual cluster
    "product": {"name": "water_cluster", "charge": 0, "multiplicity": 1},
}

DONOR_ALIASES = {
    "trifluorethanol": "tfe",
    "trifluoroethanol": "tfe",
    "tfeoh": "tfe",
    "benzimidazoline": "bih",
    "bi-dih": "bih",
    "ph": "phenol",
    "phenolat": "phenol",
    "acetic_acid": "hoac",
    "hoac": "hoac",
    "essigsaure": "hoac",
    "hydronium": "h3o",
    "h3o+": "h3o",
    "h3o(h2o)3": "h3o",
}

CANONICAL_DONOR_DIRS = {
    "tfe": "TFE_coordinator",
    "bih": "BIH_coordinator",
    "phenol": "phenol_coordinator",
    "hoac": "HOAc_coordinator",
    "h3o": "H3O+_coordinator",
}


def resolve_donor_name(raw) -> Optional[str]:
    """Map a CONTROL donor string to a canonical template key.

    Unknown values return None: the caller treats them as a custom donor
    xyz file (proton_donor_xyz must then be set).
    """
    if not isinstance(raw, str):
        return None
    key = raw.strip().lower()
    if not key or key in {"auto", "none", ""}:
        return None
    key = DONOR_ALIASES.get(key, key)
    return key if key in DONOR_TEMPLATES else None


# === XYZ helpers =============================================================

def read_xyz_robust(path: str) -> Atoms:
    """Read an xyz file, tolerating blank lines / odd separators like CO2 does."""
    atoms = read(path, format="xyz")
    return atoms


def write_xyz(atoms: Atoms, path: str) -> str:
    write(path, atoms, format="xyz")
    return path


def detect_metal_index(atoms: Atoms) -> int:
    """First transition-metal / p-block metal in the frame (CO2 convention)."""
    from ase.data import chemical_symbols
    metals = set()
    # rows 4-6 d-block + main-group metals commonly used here
    for z in range(19, 89):
        try:
            metals.add(chemical_symbols[z])
        except IndexError:
            pass
    for i, sym in enumerate(atoms.get_chemical_symbols()):
        if sym in metals:
            return i
    raise ValueError("Kein Metallatom in der Struktur gefunden.")


# === Geometry preparation ====================================================

def build_donor_atoms(donor_key: str) -> Atoms:
    """Instantiate a donor from the template table."""
    tpl = DONOR_TEMPLATES.get(donor_key)
    if tpl is None:
        raise KeyError(f"Unbekannter Donor '{donor_key}'.")
    return read_xyz_text(tpl["xyz"])


def build_donor_opt_input(donor_key: str, outdir: str, keywords: str,
                          charge: int = 0, multiplicity: int = 1,
                          PAL: int = 8, maxcore: int = 2000) -> str:
    """Write a standalone OPT input for one donor (pre-optimisation).

    The templates are sketch geometries; a relaxed scan started from them
    would inherit their strain. The donor is optimised alone, in its ground
    state, with the SAME method keywords as the later scan so the energies
    stay comparable. Returns the input path.
    """
    atoms = build_donor_atoms(donor_key)
    os.makedirs(outdir, exist_ok=True)
    inp = os.path.join(outdir, f"{donor_key}_preopt.inp")
    with open(inp, "w", newline="\n") as f:
        f.write(f"! {keywords}\n")
        f.write(f"%maxcore {maxcore}\n")
        f.write(f"%pal nprocs {PAL} end\n")
        f.write(f"* xyz {charge} {multiplicity}\n")
        for sym, pos in zip(atoms.get_chemical_symbols(), atoms.positions):
            f.write(f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
        f.write("*\n")
    return inp


def read_xyz_text(text: str) -> Atoms:
    import io
    return read(io.StringIO(text), format="xyz")


def orient_donor_toward_metal(donor: Atoms, donor_h_index: int,
                              donor_heavy_index: int,
                              metal_pos: Sequence[float],
                              distance: float = 3.0) -> Atoms:
    """Translate+rotate the donor so the donor H sits at `distance` from the
    metal along the M...H axis, with the X-H bond collinear to it.

    After this the H points AT the metal and the X-H bond stretches along the
    scan coordinate - the geometric precondition for a relaxed transfer scan.
    """
    donor = donor.copy()
    h_pos = np.asarray(donor.positions[donor_h_index], dtype=float)
    x_pos = np.asarray(donor.positions[donor_heavy_index], dtype=float)

    # unit vector along X->H (the bond to be broken/relaxed)
    bond_vec = h_pos - x_pos
    bond_len = np.linalg.norm(bond_vec)
    if bond_len < 1e-6:
        raise ValueError("Donor-Geometrie: X-H-Bindung degeneriert.")
    bond_hat = bond_vec / bond_len

    m_pos = np.asarray(metal_pos, dtype=float)
    # place H on the line metal <- H, so the H "points at" the metal
    target_h = m_pos + bond_hat * distance
    shift = target_h - h_pos
    donor.translate(shift)
    return donor


def combine_complex_donor(complex_atoms: Atoms, donor: Atoms) -> Tuple[Atoms, int]:
    """Concatenate complex + donor; returns (combined, donor_h_global_index)."""
    combined = complex_atoms.copy()
    offset = len(combined)
    for atom in donor:
        combined += atom
    # donor H is the LAST atom in every template
    return combined, offset + len(donor) - 1


def rotate_donor_about_axis(atoms: Atoms, donor_start: int, angle_deg: float,
                            axis_pos: Sequence[float],
                            axis_dir: Sequence[float]) -> Atoms:
    """Rotate donor atoms (indices >= donor_start) about axis through axis_pos."""
    out = atoms.copy()
    theta = np.deg2rad(angle_deg)
    axis = np.asarray(axis_dir, dtype=float)
    axis = axis / np.linalg.norm(axis)
    pos = np.asarray(axis_pos, dtype=float)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)
    for i in range(donor_start, len(out)):
        p = out.positions[i] - pos
        out.positions[i] = pos + R @ p
    return out


# === ORCA input builders (mirroring CO2_Coordinator6) =========================

def write_orca_sp_input(atoms: Atoms, outdir: str, tag: str,
                        keywords: str, charge: int, multiplicity: int,
                        PAL: int = 8, maxcore: int = 2000,
                        broken_sym: str = "") -> str:
    """Write an SP input for one orientation-scan point."""
    os.makedirs(outdir, exist_ok=True)
    inp = os.path.join(outdir, f"{tag}.inp")
    with open(inp, "w", newline="\n") as f:
        f.write(f"! {keywords}\n")
        f.write(f"%maxcore {maxcore}\n")
        f.write(f"%pal nprocs {PAL} end\n")
        if broken_sym:
            f.write((broken_sym if broken_sym.endswith("\n") else broken_sym + "\n"))
        f.write(f"* xyz {charge} {multiplicity}\n")
        for sym, pos in zip(atoms.get_chemical_symbols(), atoms.positions):
            f.write(f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
        f.write("*\n")
    return inp


def write_orca_scan_input(atoms: Atoms, outdir: str,
                          metal_index: int, donor_h_index: int,
                          start_distance: float, end_distance: float,
                          steps: int, keywords: str,
                          charge: int, multiplicity: int,
                          PAL: int = 8, maxcore: int = 2000,
                          broken_sym: str = "") -> str:
    """Write the relaxed M...H distance scan input (%geom Scan B)."""
    os.makedirs(outdir, exist_ok=True)
    inp = os.path.join(outdir, "scan.inp")
    with open(inp, "w", newline="\n") as f:
        f.write(f"! {keywords}\n")
        f.write(f"%maxcore {maxcore}\n")
        f.write(f"%scf\n  MaxIter 375\nend\n")
        if broken_sym:
            f.write((broken_sym if broken_sym.endswith("\n") else broken_sym + "\n"))
        f.write(f"%pal nprocs {PAL} end\n")
        f.write("%geom\n  MaxIter 400\n  Scan\n")
        f.write(f"    B  {metal_index}  {donor_h_index} = "
                f"{start_distance:.2f}, {end_distance}, {steps}\n")
        f.write("  end\nend\n")
        f.write(f"* xyz {charge} {multiplicity}\n")
        for sym, pos in zip(atoms.get_chemical_symbols(), atoms.positions):
            f.write(f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
        f.write("*\n")
    return inp


def run_orca_in(outdir: str, inp_name: str = "scan.inp",
                out_name: str = "scan.out") -> int:
    """Run ORCA inside outdir; raises if the binary is missing."""
    orca_path = shutil.which("orca")
    if orca_path is None:
        raise RuntimeError("ORCA wurde nicht gefunden! Ist es im $PATH?")
    with open(os.path.join(outdir, out_name), "w") as f:
        return subprocess.run([orca_path, inp_name], cwd=outdir,
                              stdout=f, stderr=subprocess.STDOUT).returncode


# === Scans ===================================================================

def orientation_scan(base_atoms: Atoms, donor_start_index: int,
                     donor_h_index: int, charge: int, multiplicity: int,
                     keywords: str, outdir: str = "orientation_scan",
                     PAL: int = 8, maxcore: int = 2000,
                     broken_sym: str = "",
                     angle_step_deg: int = 10, angle_range_deg: int = 180,
                     metal_index: int = 0,
                     run_calculator=None) -> Tuple[str, float, Optional[float]]:
    """Rotate ONLY the donor about the M...H axis at fixed distance and run
    SP calcs. run_calculator allows tests to inject a mock.

    metal_index locates the rotation centre (assumed first atom only when
    the frame was aligned, e.g. by CO2's align_complex).

    Returns (best_xyz_path, best_angle_deg, best_energy).
    """
    os.makedirs(outdir, exist_ok=True)
    m_pos = base_atoms.positions[metal_index]
    if run_calculator is None:
        run_calculator = _default_sp_runner
    results = []
    for ang in range(0, angle_range_deg + 1, angle_step_deg):
        rot = rotate_donor_about_axis(base_atoms, donor_start_index, ang,
                                      axis_pos=m_pos,
                                      axis_dir=base_atoms.positions[donor_h_index] - m_pos)
        tag = f"angle_{ang:03d}"
        inp = write_orca_sp_input(rot, outdir, tag, keywords, charge,
                                  multiplicity, PAL, maxcore, broken_sym)
        energy = run_calculator(outdir, tag, rot)
        results.append((ang, energy, os.path.join(outdir, tag + ".xyz")))

    finite = [(a, e) for a, e, _ in results if e is not None and np.isfinite(e)]
    if not finite:
        raise RuntimeError("Keine gueltige Energie im Orientierungsscan.")
    best_angle, best_E = min(finite, key=lambda t: t[1])
    best_xyz = os.path.join(outdir, f"angle_{best_angle:03d}.xyz")

    csv_path = os.path.join(outdir, "orientation_scan.csv")
    with open(csv_path, "w") as f:
        f.write("angle_deg,energy_Eh,relative_kcal_per_mol\n")
        ref = min(e for _, e in finite)
        for ang, energy, _ in results:
            if energy is not None and np.isfinite(energy):
                f.write(f"{ang},{energy},{(energy - ref) * 627.509:.3f}\n")
            else:
                f.write(f"{ang},,\n")
    return best_xyz, float(best_angle), float(best_E)


def _default_sp_runner(outdir, tag, atoms):
    """Real ORCA SP: write .inp, run, parse last 'Final Gibbs' or total energy."""
    write_xyz(atoms, os.path.join(outdir, tag + ".xyz"))
    run_orca_in(outdir, tag + ".inp", tag + ".out")
    out_file = os.path.join(outdir, tag + ".out")
    return parse_orca_energy(out_file)


def parse_orca_energy(out_path: str) -> Optional[float]:
    """Parse the FINAL SINGLE POINT ENERGY (Hartree) from an ORCA output."""
    if not os.path.exists(out_path):
        return None
    energy = None
    with open(out_path, "r", errors="replace") as f:
        for line in f:
            if "Final Single Point Energy" in line:
                try:
                    energy = float(line.split()[-1])
                except ValueError:
                    pass
    return energy


def relaxed_distance_scan(atoms: Atoms, metal_index: int, donor_h_index: int,
                          start_distance: float, end_distance: float,
                          steps: int, keywords: str,
                          charge: int, multiplicity: int,
                          outdir: str = "relaxed_surface_scan",
                          PAL: int = 8, maxcore: int = 2000,
                          broken_sym: str = "",
                          run_orca=None) -> str:
    """Write and (optionally) run the relaxed M...H scan; returns the .dat path."""
    os.makedirs(outdir, exist_ok=True)
    write_orca_scan_input(atoms, outdir, metal_index, donor_h_index,
                          start_distance, end_distance, steps, keywords,
                          charge, multiplicity, PAL, maxcore, broken_sym)
    if run_orca is None:
        run_orca = run_orca_in
    run_orca(outdir)
    return os.path.join(outdir, "scan.relaxscanact.dat")


def plot_scan_result(datapath: str, png_path: Optional[str] = None) -> Optional[str]:
    """Plot scan.relaxscanact.dat (distance vs E_rel) like CO2 does."""
    if not os.path.exists(datapath):
        print(f"[WARN] {datapath} nicht gefunden - kein Plot.")
        return None
    dist, energies = [], []
    with open(datapath) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    dist.append(float(parts[0]))
                    energies.append(float(parts[1]))
                except ValueError:
                    continue
    if not dist:
        return None
    ref = min(energies)
    rel = [(e - ref) * 627.509 for e in energies]
    if png_path is None:
        png_path = os.path.join(os.path.dirname(datapath),
                                "scan_relative.png")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(dist, rel, "o-", color="#1f77b4")
    ax.set_xlabel("M...H distance (Å)")
    ax.set_ylabel("ΔE (kcal/mol)")
    ax.set_title("Relaxed M...H transfer scan")
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    return png_path


# === Main flow ===============================================================

def main(control_path: str = "CONTROL.txt") -> int:
    """Standalone flow, mirroring CO2_Coordinator6.main().

    Reads a light-weight CONTROL.txt (same parser family), builds the
    <Donor>_coordinator folder for each selected donor, and drives:
    align -> place donor -> orientation scan -> relaxed scan -> plot.
    """
    from delfin.proton.chain_setup import bookkeeping, build_control

    args = _minimal_read_control_file(control_path)
    donor_raw = args.get("proton_donor", "tfe")
    donor_keys = [d.strip() for d in str(donor_raw).split(",") if d.strip()]
    xyz_in = args.get("xyz", "complex.xyz")
    charge = int(args.get("charge", 0))
    multiplicity = int(args.get("multiplicity", 1))
    PAL = int(args.get("PAL", 8))
    maxcore = int(args.get("maxcore", 2000))
    orientation_distance = float(args.get("orientation_distance", 3.0))
    rot_step_deg = int(args.get("rot_step_deg", 10))
    rot_range_deg = int(args.get("rot_range_deg", 180))
    scan_end = float(args.get("scan_end", 1.5))
    scan_steps = int(args.get("scan_steps", 15))

    from delfin.co2.CO2_Coordinator6 import build_orca_keywords, align_complex
    keywords = build_orca_keywords(args, "SP")

    for donor_key in donor_keys:
        canon = resolve_donor_name(donor_key)
        if canon is None:
            print(f"[WARN] Unbekannter Donor '{donor_key}' - uebersprungen.")
            continue
        folder = CANONICAL_DONOR_DIRS[canon]
        os.makedirs(folder, exist_ok=True)
        cwd = os.getcwd()
        os.chdir(folder)
        try:
            _run_for_donor(args, canon, xyz_in_abs=os.path.join(cwd, xyz_in),
                           charge=charge, multiplicity=multiplicity,
                           PAL=PAL, maxcore=maxcore,
                           keywords=keywords,
                           orientation_distance=orientation_distance,
                           rot_step_deg=rot_step_deg,
                           rot_range_deg=rot_range_deg,
                           scan_end=scan_end, scan_steps=scan_steps,
                           control_path=os.path.join(cwd, control_path))
        finally:
            os.chdir(cwd)
    return 0


def _run_for_donor(args, donor_key, *, xyz_in_abs, charge, multiplicity,
                   PAL, maxcore, keywords, orientation_distance,
                   rot_step_deg, rot_range_deg, scan_end, scan_steps,
                   control_path):
    """One donor's full pipeline inside its <Donor>_coordinator folder."""
    from delfin.proton.chain_setup import bookkeeping, build_control
    from delfin.co2.CO2_Coordinator6 import build_orca_keywords

    donor_tpl = DONOR_TEMPLATES[donor_key]

    # Donor pre-optimisation: the template is a sketch geometry; a scan
    # started from it would inherit its strain. A preopt result from an
    # earlier run is reused; otherwise the OPT input is (re)written and the
    # scan proceeds from the sketch WITH A WARNING (not validated chemistry).
    preopt_xyz = f"{donor_key}_preopt.xyz"
    if not os.path.exists(preopt_xyz):
        build_donor_opt_input(
            donor_key, outdir=".", keywords=keywords,
            charge=int(donor_tpl.get("donor_charge", 0)),
            multiplicity=int(donor_tpl.get("donor_multiplicity", 1)),
            PAL=PAL, maxcore=maxcore)
        print(f"[WARN] {donor_key}: no {preopt_xyz} found - {donor_key}_preopt.inp "
              "written. Run ORCA on it, rename the optimised output to "
              f"{preopt_xyz}, then rerun. Scan continues from the SKETCH "
              "geometry (not validated chemistry).")
    else:
        print(f"[INFO] {donor_key}: using pre-optimised donor geometry "
              f"{preopt_xyz}")

    donor_atoms = build_donor_atoms(donor_key)
    if os.path.exists(preopt_xyz):
        donor_atoms = read_xyz_robust(preopt_xyz)
    complex_atoms = read_xyz_robust(xyz_in_abs)
    metal_idx = detect_metal_index(complex_atoms)

    bk = bookkeeping(donor_key, complex_charge=charge,
                     complex_multiplicity=multiplicity)
    combined_charge = bk["combined_charge"]
    combined_mult = bk["combined_multiplicity"]
    print(f"[INFO] {donor_key}: combined charge={combined_charge}, "
          f"mult={combined_mult} ({bk['scheme']})")

    # place donor H toward metal
    placed = orient_donor_toward_metal(
        donor_atoms, donor_tpl["donor_h_index"], donor_tpl["donor_heavy_index"],
        metal_pos=complex_atoms.positions[metal_idx],
        distance=orientation_distance)
    combined, donor_h_global = combine_complex_donor(complex_atoms, placed)
    combined.center()
    write_xyz(combined, "complex_with_donor.xyz")

    best_xyz, best_angle, best_E = orientation_scan(
        combined, donor_start_index=len(complex_atoms),
        donor_h_index=donor_h_global, charge=combined_charge,
        multiplicity=combined_mult, keywords=keywords,
        PAL=PAL, maxcore=maxcore,
        broken_sym=str(args.get("broken_sym", "") or ""),
        angle_step_deg=rot_step_deg, angle_range_deg=rot_range_deg,
        metal_index=metal_idx)

    atoms_best = read_xyz_robust(best_xyz)
    metal_idx_best = detect_metal_index(atoms_best)
    start_distance = float(np.linalg.norm(
        atoms_best.positions[metal_idx_best]
        - atoms_best.positions[donor_h_global]))
    print(f"[INFO] Best orientation {best_angle} deg, M...H = {start_distance:.3f} A")

    scan_keywords = build_orca_keywords(args, "OPT")
    relaxed_distance_scan(atoms_best, metal_idx_best, donor_h_global,
                          start_distance=start_distance,
                          end_distance=scan_end, steps=scan_steps,
                          keywords=scan_keywords, charge=combined_charge,
                          multiplicity=combined_mult, PAL=PAL,
                          maxcore=maxcore,
                          broken_sym=str(args.get("broken_sym", "") or ""))
    plot_scan_result(os.path.join("relaxed_surface_scan", "scan.relaxscanact.dat"))

    build_control(args, donor_key, control_path=control_path,
                  outdir=os.getcwd())


def _minimal_read_control_file(path: str = "CONTROL.txt") -> Dict:
    """Light-weight key=value parser (CO2 convention, deliberately separate
    from delfin.config to stay standalone-capable)."""
    params: Dict = {}
    if not os.path.exists(path):
        print("[WARN] CONTROL.txt not found - using defaults.")
        return params
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("----"):
                continue
            if "=" in line:
                k, _, v = line.partition("=")
                params[k.strip()] = v.strip()
    return params


if __name__ == "__main__":
    import argparse
    import sys
    parser = argparse.ArgumentParser(description="Proton coordinator")
    parser.add_argument("--define", action="store_true",
                        help="Erzeuge CONTROL-Vorlage und beende.")
    cli = parser.parse_args()
    if cli.define:
        from delfin.proton.chain_setup import write_default_control
        write_default_control()
        sys.exit(0)
    sys.exit(main())
