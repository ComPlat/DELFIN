#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CO2_coordinator.py – richtet Komplex aus, platziert CO2, macht Winkel-SPs
und startet ORCA-Distanz-Scan. Parameter werden aus CONTROL.txt gelesen.
"""

import os
# --- Headless Plot Backend, bevor matplotlib importiert wird ---
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")

import re
import shutil
import subprocess
import numpy as np
from io import StringIO
from ase.io import read, write
from ase.data import covalent_radii
import matplotlib.pyplot as plt

# === Templates erzeugen (--define) ===========================================
def write_default_files(control_path="CONTROL.txt", co2_path="co2.xyz",
                        charge=None, multiplicity=None, solvent=None, metal=None,
                        overwrite=False):
    # 'metal=auto' signals runtime detection from the xyz (see main()).
    control_template = """# Input / Output
------------------------------------
xyz=input6.xyz
out=complex_aligned.xyz
co2=co2.xyz

# Charge & Multiplicity
------------------------------------
charge=[CHARGE]
multiplicity=[MULTIPLICITY]
additions=

# Solvation
------------------------------------
implicit_solvation_model=CPCM
solvent=[SOLVENT]

# Orientation Scan (single points)
------------------------------------
orientation_distance=5.0
rot_step_deg=10
rot_range_deg=180

# Method Settings
------------------------------------
functional=PBE0
disp_corr=D4
ri_jkx=RIJCOSX
aux_jk=def2/J
main_basisset=def2-SVP
metal_basisset=def2-TZVP
orientation_job=SP
scan_job=OPT

# Relaxed Distance Scan
------------------------------------
scan_end=1.6
scan_steps=25

# Alignment (0-based indices)
------------------------------------
metal=auto
metal_index=
align_bond_index=
neighbors=

# CO2 placement
------------------------------------
place_axis=z
mode=side-on
perp_axis=y
no_place_co2=false

# Resources
------------------------------------
PAL=32
maxcore=3800

# Alternative keywords (commented examples)
# orientation_job=GFN2-XTB
# scan_job=GFN2-XTB OPT
# additions=%SCF BrokenSym M,N END
"""

    # Platzhalter optional ersetzen (sonst bleiben sie wie im Template)
    repl = {
        "[CHARGE]":       str(charge) if charge is not None else "[CHARGE]",
        "[MULTIPLICITY]": str(multiplicity) if multiplicity is not None else "[MULTIPLICITY]",
        "[SOLVENT]":      solvent if solvent is not None else "[SOLVENT]",
    }
    for k, v in repl.items():
        control_template = control_template.replace(k, v)

    if metal is not None:
        metal_str = str(metal).strip()
        control_template = control_template.replace("metal=auto", f"metal={metal_str}")

    co2_xyz = """3

O      0.000000    0.000000    1.840000
C      0.000000    0.000000    3.000000
O      0.000000    0.000000    4.160000
"""

    def _write(path, text):
        import sys
        if os.path.exists(path) and not overwrite:
            print(f"[INFO] {path} existiert bereits – nichts geschrieben (nutze --force zum Überschreiben).", file=sys.stderr)
            return
        with open(path, "w", newline="\n") as f:
            f.write(text)
        print(f"[OK] geschrieben: {path}")

    _write(control_path, control_template)
    _write(co2_path, co2_xyz)


# === CONTROL.txt einlesen ===
def read_control_file(path="CONTROL.txt"):
    params = {}
    if not os.path.exists(path):
        print(f"[WARN] CONTROL.txt not found → using defaults.")
        return params

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            key, val = map(str.strip, line.split("=", 1))

            # Boolesche Werte
            if isinstance(val, str):
                low = val.lower()
                if low == "true":
                    val = True
                elif low == "false":
                    val = False
                elif val == "":
                    val = None
                else:
                    # Versuch: Zahl (int oder float)
                    try:
                        if "." in val:
                            val = float(val)
                        else:
                            val = int(val)
                    except ValueError:
                        pass  # bleibt String

            params[key] = val

    # Explizite Typanpassung
    for key in ["distance", "scan_end", "orientation_distance"]:
        if key in params and isinstance(params[key], str):
            params[key] = float(params[key])
    for key in ["scan_steps", "charge", "multiplicity", "PAL", "maxcore", "rot_step_deg", "rot_range_deg"]:
        if key in params and isinstance(params[key], str):
            params[key] = int(params[key])

    # Optional: /n durch Zeilenumbruch ersetzen
    for key in ["orca_keywords", "rot_orca_keywords", "additions"]:
        if key in params and isinstance(params[key], str):
            params[key] = params[key].replace("/n", "\n")

    return params


# === Geometrie und Rotation ===
def rot_from_vecs(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-12:
        if c > 0:
            return np.eye(3)
        axis = np.array([1., 0., 0.]) if abs(a[0]) < 0.9 else np.array([0., 1., 0.])
        v = np.cross(a, axis)
        v /= np.linalg.norm(v)
        K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + 2 * (K @ K)
    s = np.linalg.norm(v)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * ((1 - c) / (s ** 2))

def Rz(angle_deg):
    th = np.deg2rad(angle_deg)
    c, s = np.cos(th), np.sin(th)
    return np.array([[c, -s, 0.0],
                     [s,  c, 0.0],
                     [0.0, 0.0, 1.0]])

def project_to_plane(v, n):
    n = n / np.linalg.norm(n)
    return v - np.dot(v, n) * n

def principal_plane_normal(vectors):
    M = np.stack(vectors, axis=0)
    C = M.T @ M
    _, eigvecs = np.linalg.eigh(C)
    return eigvecs[:, 0] / np.linalg.norm(eigvecs[:, 0])

# === XYZ robust lesen ===
def _read_xyz_robust(path):
    try:
        return read(path)
    except UnicodeDecodeError:
        with open(path, 'r', encoding='cp1252', errors='replace') as f:
            txt = f.read()
        txt = txt.replace('–', '-').replace('—', '-').replace('−', '-')
        return read(StringIO(txt), format='xyz')

# === Metall erkennen ===
METAL_SYMBOLS = set("""
Li Be Na Mg K Ca Rb Sr Cs Ba Fr Ra
Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd
Hf Ta W Re Os Ir Pt Au Hg La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr
Al Ga In Tl Sn Pb Bi Po
""".split())

def detect_metal_index(atoms):
    c = [i for i, a in enumerate(atoms) if a.symbol in METAL_SYMBOLS]
    if not c:
        nonmetals = set(list("HBCNOFPSI") + ["Se", "Te", "I", "Br", "Cl", "F", "Ne", "Ar", "Kr", "Xe", "Rn", "Og", "He"])
        c = [i for i, a in enumerate(atoms) if a.symbol not in nonmetals]
        if not c:
            return int(np.argmax([a.number for a in atoms]))
    return c[0] if len(c) == 1 else max(c, key=lambda i: atoms[i].number)

def guess_neighbors(atoms, metal_index, scale=1.15):
    ZM = atoms[metal_index].number
    rM = covalent_radii[ZM]
    Mpos = atoms.positions[metal_index]
    neigh = []
    for i, a in enumerate(atoms):
        if i == metal_index:
            continue
        ri = covalent_radii[a.number]
        cutoff = scale * (rM + ri)
        if np.linalg.norm(atoms.positions[i] - Mpos) <= cutoff:
            neigh.append(i)
    return neigh

# === Komplex ausrichten ===
def align_complex(infile, outfile, metal_index=None, metal_symbol=None, align_bond_index=None, neighbor_indices=None):
    atoms = _read_xyz_robust(infile)
    M = metal_index if metal_index is not None else detect_metal_index(atoms)
    atoms.positions -= atoms.positions[M]

    if not neighbor_indices:
        neighbor_indices = guess_neighbors(atoms, M)
        if not neighbor_indices:
            raise ValueError("Keine Liganden erkannt.")

    ML = [atoms.positions[i] for i in neighbor_indices]
    n = principal_plane_normal(ML)
    R1 = rot_from_vecs(n, np.array([0., 0., 1.]))
    atoms.positions = atoms.positions @ R1.T
    ML_rot = [v @ R1.T for v in ML]

    pick = (np.argmax([np.linalg.norm(project_to_plane(v, [0, 0, 1])) for v in ML_rot])
            if align_bond_index is None else neighbor_indices.index(align_bond_index))
    v_pick = project_to_plane(ML_rot[pick], [0, 0, 1])
    v_pick /= np.linalg.norm(v_pick)

    K = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    phi = np.arctan2(np.dot(np.cross(v_pick, [0, 1, 0]), [0, 0, 1]), np.dot(v_pick, [0, 1, 0]))
    R2 = np.eye(3) + np.sin(phi) * K + (1 - np.cos(phi)) * (K @ K)
    atoms.positions = atoms.positions @ R2.T

    write(outfile, atoms)
    print(f"[OK] wrote {outfile}")
    return outfile

# === CO2 platzieren ===
def _axis_vector(name):
    return {"x": np.array([1, 0, 0]), "y": np.array([0, 1, 0]), "z": np.array([0, 0, 1])}[name]

def _co2_axis_center_indices(atoms):
    syms = atoms.get_chemical_symbols()
    c = [i for i, s in enumerate(syms) if s == "C"]
    o = [i for i, s in enumerate(syms) if s == "O"]
    if len(c) != 1 or len(o) != 2:
        raise ValueError("CO2 muss 1 C und 2 O enthalten.")
    axis = atoms.positions[o[0]] - atoms.positions[o[1]]
    axis /= np.linalg.norm(axis)
    center = atoms.positions[c[0]]
    return axis, center, c[0]

def place_co2_general(complex_path, co2_path, out_path, distance=5.0, place_axis='z', mode='side-on', perp_axis='y'):
    """
    Fügt CO2 an +place_axis in 'distance' Å an. Gibt zusätzlich CO2-Indizes im kombinierten System zurück.
    """
    comp = _read_xyz_robust(complex_path)
    co2 = _read_xyz_robust(co2_path)

    axis, center, c_idx_local = _co2_axis_center_indices(co2)
    co2.positions -= center  # CO2 um sein C zentrieren

    target = (_axis_vector(perp_axis) if mode == "side-on" else _axis_vector(place_axis))
    R = rot_from_vecs(axis, target)
    co2.positions = co2.positions @ R.T + distance * _axis_vector(place_axis)

    combined = comp + co2
    write(out_path, combined)
    print(f"[OK] wrote combined (complex + CO2) → {out_path}")

    # Indizes der CO2-Atome im kombinierten System
    n_comp = len(comp)
    co2_indices = list(range(n_comp, n_comp + len(co2)))
    co2_c_index_combined = co2_indices[c_idx_local]
    return out_path, co2_indices, co2_c_index_combined

# === ORCA Helpers ===
def parse_orca_energy(out_path):
    """
    Liefert Energie in Hartree. Sucht robust nach:
    - 'FINAL SINGLE POINT ENERGY' (ORCA)
    - 'TOTAL ENERGY' (xTB-Driver in ORCA)
    """
    with open(out_path, "r", errors="ignore") as f:
        lines = f.readlines()

    for line in reversed(lines):
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                return float(line.strip().split()[-1])
            except Exception:
                pass

    for line in reversed(lines):
        if "TOTAL ENERGY" in line and "MB" not in line:
            m = re.search(r"(-?\d+\.\d+)", line)
            if m:
                return float(m.group(1))

    raise RuntimeError(f"Energie in '{out_path}' nicht gefunden.")

def _format_basis_block(metal_symbol, metal_basis, keywords):
    if not metal_symbol or not metal_basis:
        return ""
    if not isinstance(metal_symbol, str) or not isinstance(metal_basis, str):
        return ""
    # Skip explicit basis override for semiempirical runs
    kw_upper = keywords.upper()
    if "GFN" in kw_upper or "XTB" in kw_upper:
        return ""
    sym = metal_symbol.strip()
    basis = metal_basis.strip()
    if not sym or not basis:
        return ""
    block = ["%basis", f"  NewGTO {sym} \"{basis}\" end", "end"]
    return "\n".join(block)


def _clean_str(value, default=""):
    if value is None:
        return default
    if isinstance(value, bool):
        return default
    text = str(value).strip()
    if not text:
        return default
    if text.startswith("[") and text.endswith("]"):
        return default
    return text


def build_orca_keywords(config, job_spec):
    """Construct ORCA keyword string from control-style settings."""
    job = _clean_str(job_spec)
    functional = _clean_str(config.get("functional", "PBE0"), "PBE0")
    tokens = []
    if functional:
        tokens.append(functional)

    functional_upper = functional.upper()
    is_gfn = "GFN" in functional_upper or "XTB" in functional_upper

    if is_gfn:
        solvent = _clean_str(config.get("solvent"))
        model = _clean_str(config.get("implicit_solvation_model", "ALPB"), "ALPB")
        if solvent:
            if model:
                tokens.append(f"{model}({solvent})")
        elif model:
            tokens.append(model)
        if job:
            tokens.append(job)
        return " ".join(tokens).strip()

    main_basis = _clean_str(config.get("main_basisset", "def2-SVP"))
    if main_basis:
        tokens.append(main_basis)
    aux_jk = _clean_str(config.get("aux_jk", "def2/J"))
    if aux_jk:
        tokens.append(aux_jk)
    disp_corr = _clean_str(config.get("disp_corr", "D4"))
    if disp_corr:
        tokens.append(disp_corr)

    solvent = _clean_str(config.get("solvent"))
    implicit_model = _clean_str(config.get("implicit_solvation_model", "CPCM"), "CPCM")
    if solvent:
        model = implicit_model if implicit_model else "CPCM"
        tokens.append(f"{model}({solvent})")

    ri_jkx = _clean_str(config.get("ri_jkx"))
    if ri_jkx:
        tokens.append(ri_jkx)
    relativity = _clean_str(config.get("relativity"))
    if relativity:
        tokens.append(relativity)
    ri_soc = _clean_str(config.get("ri_soc"))
    if ri_soc:
        tokens.append(ri_soc)

    if job:
        tokens.append(job)

    return " ".join(tokens).strip()


def write_orca_sp_input_and_run(xyz_path, outdir, orca_keywords="GFN2-XTB SP ALPB(DMF)", additions="",
                                charge=-2, multiplicity=1, PAL=8, maxcore=2000, tag="calc",
                                metal_symbol=None, metal_basis=None):
    os.makedirs(outdir, exist_ok=True)
    inp = os.path.join(outdir, f"{tag}.inp")
    out = os.path.join(outdir, f"{tag}.out")

    xyz_basename = os.path.basename(xyz_path)
    xyz_target = os.path.join(outdir, xyz_basename)

    basis_block = _format_basis_block(metal_symbol, metal_basis, orca_keywords)
    sections = [f"! {orca_keywords}"]
    if basis_block:
        sections.append(basis_block)
    sections.append(f"%maxcore {maxcore}")
    sections.append(f"%pal nprocs {PAL} end")
    if additions:
        sections.append(additions)
    sections.append("")
    sections.append(f"* xyzfile {charge} {multiplicity} {xyz_basename} *")
    input_text = "\n".join(sections) + "\n"
    with open(inp, "w") as f:
        f.write(input_text)

    # Nur kopieren, wenn Quelle und Ziel verschieden sind
    if os.path.abspath(xyz_path) != os.path.abspath(xyz_target):
        shutil.copy(xyz_path, xyz_target)

    orca_path = shutil.which("orca")
    if orca_path is None:
        raise RuntimeError("ORCA wurde nicht gefunden! Ist es im $PATH?")

    with open(out, "w") as f:
        subprocess.run([orca_path, os.path.basename(inp)], cwd=outdir,
                       stdout=f, stderr=subprocess.STDOUT, check=False)
    return out


def write_orca_input_and_run(xyz_path, metal_index, co2_c_index, start_distance, end_distance=1.7, steps=5,
                             orca_keywords="GFN2-XTB OPT ALPB(DMF)", additions="",
                             charge=-2, multiplicity=3, PAL=4, maxcore=2000,
                             metal_symbol=None, metal_basis=None):
    additions_line = additions if additions else ""
    # ORCA ist 1-basiert → ACHTUNG: falls 'metal_index'/'co2_c_index' 0-basiert sind, hier ggf. +1 setzen.
    i_orca = metal_index
    j_orca = co2_c_index

    basis_block = _format_basis_block(metal_symbol, metal_basis, orca_keywords)
    sections = [f"! {orca_keywords}"]
    if basis_block:
        sections.append(basis_block)
    sections.append(f"%maxcore {maxcore}")
    sections.append(f"%pal nprocs {PAL} end")
    if additions_line:
        sections.append(additions_line)
    sections.append("%geom")
    sections.append("  MaxIter 200")
    sections.append("  Scan")
    sections.append(f"    B  {i_orca}  {j_orca} = {start_distance:.2f}, {end_distance}, {steps}")
    sections.append("  end")
    sections.append("end")
    sections.append("")
    sections.append(f"* xyzfile {charge} {multiplicity} {os.path.basename(xyz_path)} *")
    input_text = "\n".join(sections) + "\n"
    os.makedirs("relaxed_surface_scan", exist_ok=True)
    inp_path = os.path.join("relaxed_surface_scan", "scan.inp")
    xyz_target = os.path.join("relaxed_surface_scan", os.path.basename(xyz_path))

    with open(inp_path, "w") as f:
        f.write(input_text)
    shutil.copy(xyz_path, xyz_target)

    orca_path = shutil.which("orca")
    if orca_path is None:
        raise RuntimeError("ORCA wurde nicht gefunden! Ist es im $PATH?")

    print(f"[INFO] Scan-Bindung: B {i_orca} {j_orca} (M_idx0={metal_index}, CO2_C_idx0={co2_c_index})")
    print("[INFO] Starte ORCA (dist.-scan):", orca_path)
    with open(os.path.join("relaxed_surface_scan", "scan.out"), "w") as f:
        subprocess.run([orca_path, "scan.inp"], cwd="relaxed_surface_scan",
                       stdout=f, stderr=subprocess.STDOUT)

# === Plot ===
def plot_scan_result(datapath):
    if not os.path.exists(datapath):
        print(f"[WARNING] File '{datapath}' not found – no plot generated.")
        return

    data = np.loadtxt(datapath)
    distances = data[:, 0]
    energies_kcal = data[:, 1] * 627.509  # Eh → kcal/mol

    # --- absolut ---
    plt.figure()
    plt.plot(distances, energies_kcal, marker='o')
    plt.xlabel("M–C distance [Å]")
    plt.ylabel("Energy [kcal/mol]")
    plt.title("Relaxed Surface Scan (absolute)")
    plt.grid()
    plt.tight_layout()
    plt.savefig("relaxed_surface_scan/scan_absolute.png")

    # --- relativ zum ersten Punkt (typisch 5 Å) ---
    ref_idx = 0
    ref_d = distances[ref_idx]
    ref_E = energies_kcal[ref_idx]
    rel_energies = energies_kcal - ref_E

    plt.figure()
    plt.plot(distances, rel_energies, marker='o')
    plt.xlabel("M–C distance [Å]")
    plt.ylabel("ΔE [kcal/mol] (rel. to first point)")
    plt.title(f"Relaxed Surface Scan (ΔE vs {ref_d:.2f} Å)")
    plt.grid()
    plt.tight_layout()
    plt.savefig("relaxed_surface_scan/scan_relative.png")

    print("[OK] Plots saved: scan_absolute.png, scan_relative.png")


def plot_orientation_result(csv_path, png_path):
    if not os.path.exists(csv_path):
        return
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    ang = data[:, 0]
    rel = data[:, 2]
    plt.figure()
    plt.plot(ang, rel, marker='o')
    plt.xlabel("Rotation angle about z [deg]")
    plt.ylabel("ΔE [kcal/mol] (rel. to min)")
    plt.title("CO₂ orientation scan (fixed distance)")
    plt.grid()
    plt.tight_layout()
    plt.savefig(png_path)
    plt.close()
    print(f"[OK] Orientation plot saved: {png_path}")

# === Orientation scan (0–180°, SPs) ===
def orientation_scan_at_fixed_distance(combined_xyz_path, co2_indices, charge, multiplicity,
                                       PAL, maxcore, orca_keywords, additions,
                                       angle_step_deg=10, angle_range_deg=180,
                                       metal_symbol=None, metal_basis=None):
    """
    Dreht NUR die CO2-Atome um die z-Achse (durch den Ursprung) auf ihrer Position (z=const),
    macht für jeden Winkel eine SP-Rechnung und liefert die beste Geometrie zurück.
    """
    base_atoms = _read_xyz_robust(combined_xyz_path)
    os.makedirs("orientation_scan", exist_ok=True)

    # Winkel-Liste 0..angle_range_deg inkl. Endpunkt
    angles = list(range(0, angle_range_deg + 1, angle_step_deg))
    results = []  # (angle_deg, energy_Eh, energy_kcal)

    best_energy = None
    best_angle = None
    best_xyz = None

    for ang in angles:
        atoms = base_atoms.copy()
        R = Rz(ang)
        pos = atoms.positions.copy()
        pos[co2_indices] = (pos[co2_indices] @ R.T)  # rotiere nur CO2
        atoms.positions = pos

        ang_dir = os.path.join("orientation_scan", f"ang_{ang:03d}")
        os.makedirs(ang_dir, exist_ok=True)
        xyz_path = os.path.join(ang_dir, "structure.xyz")
        write(xyz_path, atoms)

        out_path = write_orca_sp_input_and_run(
            xyz_path, ang_dir,
            orca_keywords=orca_keywords,
            additions=additions,
            charge=charge, multiplicity=multiplicity,
            PAL=PAL, maxcore=maxcore, tag="calc",
            metal_symbol=metal_symbol, metal_basis=metal_basis
        )
        try:
            E = parse_orca_energy(out_path)  # Hartree
        except Exception as e:
            print(f"[WARN] Angle {ang}°: energy parse failed: {e}")
            E = np.nan

        kcal = E * 627.509 if np.isfinite(E) else np.nan
        results.append((ang, E, kcal))

        if np.isfinite(E) and (best_energy is None or E < best_energy):
            best_energy = E
            best_angle = ang
            best_xyz = xyz_path

        if np.isfinite(E):
            print(f"[INFO] angle {ang:3d}° → E = {E:.10f} Eh")
        else:
            print(f"[INFO] angle {ang:3d}° → E = NaN")

    # Save CSV
    csv_path = os.path.join("orientation_scan", "orientation_scan.csv")
    with open(csv_path, "w") as f:
        f.write("angle_deg,energy_Eh,relative_kcal_per_mol\n")
        finite_kcal = [k for _, _, k in results if np.isfinite(k)]
        ref = min(finite_kcal) if finite_kcal else np.nan
        for ang, E, kcal in results:
            rel = (kcal - ref) if (np.isfinite(kcal) and np.isfinite(ref)) else np.nan
            f.write(f"{ang},{E if np.isfinite(E) else ''},{rel if np.isfinite(rel) else ''}\n")

    plot_orientation_result(csv_path, os.path.join("orientation_scan", "orientation_relative.png"))

    if best_xyz is None:
        raise RuntimeError("Keine gültige Energie im Orientierungsscan gefunden.")
    print(f"[OK] Best orientation: {best_angle}°")
    return best_xyz, best_angle, best_energy

# === Main ===
def main():
    args = read_control_file()

    # ---- Defaults / CONTROL-Parameter ----
    xyz_in        = args.get("xyz", "complex.xyz")
    xyz_out_align = args.get("out", "complex_aligned.xyz")
    co2_path      = args.get("co2", "co2.xyz")

    # Handle metal specification: allow manual override, otherwise auto-detect
    metal_setting = args.get("metal")
    metal_symbol = None
    if isinstance(metal_setting, str):
        metal_stripped = metal_setting.strip()
        if metal_stripped and metal_stripped.lower() not in {"auto", "[metal]"}:
            metal_symbol = metal_stripped
    elif metal_setting:
        metal_symbol = str(metal_setting).strip()

    if metal_symbol:
        pass
    else:
        try:
            atoms_input = _read_xyz_robust(xyz_in)
            metal_idx0 = detect_metal_index(atoms_input)
            metal_symbol = atoms_input[metal_idx0].symbol
            print(f"[INFO] Metall automatisch erkannt: {metal_symbol} (Index {metal_idx0})")
        except Exception as exc:
            print(f"[WARN] Metall konnte nicht automatisch aus '{xyz_in}' ermittelt werden: {exc}")
            metal_symbol = None

    args["metal"] = metal_symbol
    metal_basis = _clean_str(args.get("metal_basisset", "def2-TZVP")) or None

    def _apply_metal_placeholder(value):
        if isinstance(value, str) and metal_symbol:
            return value.replace("[METAL]", metal_symbol)
        return value

    # Orientation at fixed distance (defaults to 5.0 Å)
    orientation_distance = args.get("orientation_distance", 5.0)
    rot_step_deg = args.get("rot_step_deg", 10)     # 0,10,20,...,180
    rot_range_deg = args.get("rot_range_deg", 180)

    # Build ORCA keyword strings from Delfin-like controls
    orientation_job = args.get("orientation_job", "SP")
    scan_job = args.get("scan_job", "OPT")
    rot_orca_keywords = build_orca_keywords(args, orientation_job)
    scan_orca_keywords = build_orca_keywords(args, scan_job)

    if "rot_orca_keywords" in args:
        custom_rot = _apply_metal_placeholder(args.get("rot_orca_keywords"))
        if _clean_str(custom_rot):
            rot_orca_keywords = custom_rot
    if "orca_keywords" in args:
        custom_scan = _apply_metal_placeholder(args.get("orca_keywords"))
        if _clean_str(custom_scan):
            scan_orca_keywords = custom_scan

    if not rot_orca_keywords:
        raise ValueError("Keine gültigen ORCA-Schlüsselwörter für den Orientierungsscan gefunden. Bitte CONTROL-Einträge prüfen.")
    if not scan_orca_keywords:
        raise ValueError("Keine gültigen ORCA-Schlüsselwörter für den Distanzscan gefunden. Bitte CONTROL-Einträge prüfen.")

    additions_raw = args.get("additions", "")
    additions = _apply_metal_placeholder(additions_raw)

    if metal_symbol is None and isinstance(additions_raw, str) and "[METAL]" in additions_raw:
        print("[WARN] Platzhalter [METAL] in 'additions' nicht ersetzt. Bitte Metall manuell angeben.")

    # Charge/mult etc.
    charge = args.get("charge", -2)
    multiplicity = args.get("multiplicity", 1)
    PAL = args.get("PAL", 32)
    maxcore = args.get("maxcore", 3800)

    # Distance scan settings
    scan_end = args.get("scan_end", 1.7)
    scan_steps = args.get("scan_steps", 15)

    # --- 1) Align complex ---
    neighbors = [int(i) for i in args["neighbors"].split(",")] if args.get("neighbors") else None
    aligned = align_complex(
        xyz_in, xyz_out_align,
        args.get("metal_index"), args.get("metal"),
        align_bond_index=args.get("align_bond_index"),
        neighbor_indices=neighbors
    )

    # --- 2) Place CO2 at fixed distance on +z (default 5.0 Å) ---
    if not os.path.exists(co2_path):
        raise FileNotFoundError(f"CO2-Datei nicht gefunden: {co2_path}")

    combined_out = xyz_out_align.rsplit(".", 1)[0] + "_with_CO2.xyz"
    combined_path, co2_indices, co2_c_idx = place_co2_general(
        aligned, co2_path, combined_out,
        distance=orientation_distance,
        place_axis=args.get("place_axis", "z"),
        mode=args.get("mode", "side-on"),
        perp_axis=args.get("perp_axis", "y")
    )

    # --- 3) Orientation scan at fixed distance ---
    best_xyz_path, best_angle_deg, best_E = orientation_scan_at_fixed_distance(
        combined_xyz_path=combined_path,
        co2_indices=co2_indices,
        charge=charge, multiplicity=multiplicity,
        PAL=PAL, maxcore=maxcore,
        orca_keywords=rot_orca_keywords,
        additions=additions,
        angle_step_deg=rot_step_deg,
        angle_range_deg=rot_range_deg,
        metal_symbol=metal_symbol,
        metal_basis=metal_basis
    )

    # --- 4) Use best orientation to start relaxed distance scan ---
    atoms_best = _read_xyz_robust(best_xyz_path)
    metal_idx = detect_metal_index(atoms_best)
    metal_symbol_scan = metal_symbol or atoms_best[metal_idx].symbol
    start_distance = np.linalg.norm(atoms_best.positions[metal_idx] - atoms_best.positions[co2_c_idx])
    print(f"[INFO] Using best angle {best_angle_deg}° for distance scan. Start M–C = {start_distance:.3f} Å")

    write_orca_input_and_run(
        best_xyz_path,
        metal_idx,
        co2_c_idx,
        start_distance=start_distance,
        end_distance=scan_end,
        steps=scan_steps,
        orca_keywords=scan_orca_keywords,
        additions=additions,
        charge=charge,
        multiplicity=multiplicity,
        PAL=PAL,
        maxcore=maxcore,
        metal_symbol=metal_symbol_scan,
        metal_basis=metal_basis
    )

    # --- 5) Plot distance scan results (if present) ---
    plot_scan_result(os.path.join("relaxed_surface_scan", "scan.relaxscanact.dat"))

if __name__ == "__main__":
    import argparse, sys
    parser = argparse.ArgumentParser(description="CO2_coordinator.py")
    parser.add_argument("--define", action="store_true",
                        help="Erzeuge CONTROL.txt und co2.xyz und beende.")
    parser.add_argument("--force", action="store_true",
                        help="Vorhandene Dateien überschreiben.")
    # Optional: Platzhalter direkt befüllen
    parser.add_argument("--charge", type=int, help="ersetzt [CHARGE] im CONTROL-Template")
    parser.add_argument("--multiplicity", type=int, help="ersetzt [MULTIPLICITY] im CONTROL-Template")
    parser.add_argument("--solvent", type=str, help="ersetzt [SOLVENT] im CONTROL-Template, z.B. DMF")
    parser.add_argument("--metal", type=str, help='ersetzt [METAL] im CONTROL-Template')

    cli = parser.parse_args()

    if cli.define:
        write_default_files(charge=cli.charge,
                            multiplicity=cli.multiplicity,
                            solvent=cli.solvent,
                            metal=cli.metal,
                            overwrite=cli.force)
        sys.exit(0)

    # normaler Ablauf
    main()
