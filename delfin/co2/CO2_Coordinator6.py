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
from typing import Optional, Tuple, List, Dict, Sequence
from ase.io import read, write
from ase.data import covalent_radii
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

# === Templates erzeugen (--define) ===========================================
def write_default_files(control_path="CONTROL.txt", co2_path="co2.xyz",
                        charge=None, multiplicity=None, solvent=None, metal=None,
                        additions=None, overwrite=False):
    # 'metal=auto' signals runtime detection from the xyz (see main()).
    control_template = """# Input / Output
------------------------------------
xyz=[FILE]
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
orientation_distance=4.0
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
first_coordination_sphere_metal_basisset=no
first_coordination_sphere_scale=1.20
second_coordination_sphere_metal_basisset=no
second_coordination_sphere_scale=1.30
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
place_optimize=true
place_samples=800
place_clearance_scale=1.0
no_place_co2=false

# Resources
------------------------------------
PAL=12
maxcore=3800

# Parallelization (orientation scan only)
------------------------------------
parallel_orientation_scan=true
max_workers=4

# Alternative keywords (commented examples)
# orientation_job=GFN2-XTB
# scan_job=GFN2-XTB OPT
# additions=%SCF BrokenSym M,N END
#
# Parallelization notes:
# - Each worker uses PAL cores for ORCA
# - max_workers is auto-calculated as: total_cores / PAL
# - Example: 32 cores with PAL=32 → 1 worker (sequential)
# - Example: 32 cores with PAL=8 → 4 workers (parallel)
# - Set max_workers explicitly to override (but stay within limits!)
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

    if additions is not None:
        additions_str = str(additions).strip()
        control_template = control_template.replace("additions=", f"additions={additions_str}")

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
    for key in ["distance", "scan_end", "orientation_distance", "place_clearance_scale"]:
        if key in params and isinstance(params[key], str):
            params[key] = float(params[key])
    for key in ["scan_steps", "charge", "multiplicity", "PAL", "maxcore", "rot_step_deg",
                "rot_range_deg", "place_samples"]:
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
_SEPARATOR_CHARS = {"$", "*"}


def _is_separator_line(stripped: str) -> bool:
    return bool(stripped) and all(ch in _SEPARATOR_CHARS for ch in stripped)


def _read_xyz_robust(path):
    def _sanitize(text: str) -> str:
        lines = []
        for raw in text.splitlines():
            stripped = raw.strip()
            if _is_separator_line(stripped):
                continue
            lines.append(raw)
        return "\n".join(lines)

    def _read_from_text(text: str):
        sanitized = _sanitize(text)
        return read(StringIO(sanitized), format="xyz")

    try:
        with open(path, "r", encoding="utf-8") as f:
            txt = f.read()
        return _read_from_text(txt)
    except UnicodeDecodeError:
        with open(path, "r", encoding="cp1252", errors="replace") as f:
            txt = f.read()
        txt = txt.replace("–", "-").replace("—", "-").replace("−", "-")
        return _read_from_text(txt)
    except Exception:
        return read(path)


def _read_xyz_lines(path: str) -> Sequence[str]:
    try:
        with open(path, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except UnicodeDecodeError:
        with open(path, "r", encoding="cp1252", errors="replace") as f:
            lines = f.readlines()
    return lines


def _detect_qmmm_metadata_from_lines(lines: Sequence[str]) -> Tuple[Optional[Tuple[int, int]], Optional[str]]:
    qm_count = 0
    separator = None
    for raw in lines:
        stripped = raw.strip()
        if not stripped:
            continue
        if _is_separator_line(stripped):
            separator = stripped
            break
        qm_count += 1
    if separator and qm_count > 0:
        return (0, qm_count - 1), separator
    return None, separator


def detect_qmmm_range_from_xyz(path: str) -> Optional[Tuple[int, int]]:
    """
    Returns (start, end) indices (0-based) for QM atoms based on separator lines in XYZ.
    If no separator is present, returns None.
    """
    lines = _read_xyz_lines(path)

    return detect_qmmm_range_from_lines(lines[2:])  # skip atom count + comment


def detect_qmmm_range_from_lines(lines: Sequence[str]) -> Optional[Tuple[int, int]]:
    range_info, _ = _detect_qmmm_metadata_from_lines(lines)
    return range_info


def detect_qmmm_separator_from_xyz(path: str) -> Optional[str]:
    lines = _read_xyz_lines(path)
    return detect_qmmm_separator_from_lines(lines[2:])


def detect_qmmm_separator_from_lines(lines: Sequence[str]) -> Optional[str]:
    _, separator = _detect_qmmm_metadata_from_lines(lines)
    return separator


def _write_xyz_with_separator(atoms, path: str, qm_count: Optional[int] = None, comment: str = "",
                              separator: str = "$", insert_separator: bool = True):
    """Write XYZ file, optionally inserting separator line between QM and MM atoms when qm_count is provided."""
    n_atoms = len(atoms)
    if comment is None:
        comment = ""
    comment_line = comment.strip() if comment else "generated by CO2_Coordinator6"
    sep_line = (separator or "$").strip() or "$"
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{n_atoms}\n")
        f.write(comment_line + "\n")
        for idx, atom in enumerate(atoms):
            if insert_separator and qm_count is not None and idx == qm_count:
                f.write(sep_line + "\n")
            x, y, z = atom.position
            f.write(f"{atom.symbol:<3} {x:>14.8f} {y:>14.8f} {z:>14.8f}\n")


def _is_enabled(value) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"yes", "true", "1", "on"}


def _parse_float(value, default: float) -> float:
    if value is None or value == "":
        return default
    try:
        return float(str(value).strip())
    except (ValueError, TypeError):
        return default


def _collect_metal_indices(atoms) -> List[int]:
    indices = []
    for idx, atom in enumerate(atoms):
        if atom.symbol.capitalize() in METAL_SYMBOLS:
            indices.append(idx)
    return indices


def _first_coordination_sphere_indices(atoms, metal_indices: Sequence[int], scale: float) -> set:
    first = set()
    if not metal_indices:
        return first
    scale = max(0.0, float(scale))
    for im in metal_indices:
        Mpos = atoms.positions[im]
        rM = covalent_radii[atoms[im].number]
        for idx, atom in enumerate(atoms):
            if idx == im or idx in metal_indices:
                continue
            cutoff = scale * (rM + covalent_radii[atom.number])
            if np.linalg.norm(atoms.positions[idx] - Mpos) <= cutoff:
                first.add(idx)
    return first


def _second_coordination_sphere_indices(atoms, first_indices: Sequence[int], metal_indices: Sequence[int], scale: float) -> set:
    second = set()
    if not first_indices:
        return second
    scale = max(0.0, float(scale))
    first_set = set(first_indices)
    metal_set = set(metal_indices)
    for jf in first_indices:
        pos_f = atoms.positions[jf]
        r_f = covalent_radii[atoms[jf].number]
        for idx, atom in enumerate(atoms):
            if idx == jf or idx in first_set or idx in metal_set:
                continue
            cutoff = scale * (r_f + covalent_radii[atom.number])
            if np.linalg.norm(atoms.positions[idx] - pos_f) <= cutoff:
                second.add(idx)
    return second


def _is_semiempirical_keywords(keywords: str) -> bool:
    kw = (keywords or "").upper()
    if "QM/XTB" in kw:
        return False
    return "GFN" in kw or "XTB" in kw


def _ensure_qmmm_keyword(keywords: str, qmmm_range: Optional[Tuple[int, int]]) -> str:
    if not qmmm_range:
        return keywords
    kw = keywords or ""
    if "QM/XTB" in kw.upper():
        return kw
    return (kw + " QM/XTB").strip()


def _build_qmmm_block_lines(qmmm_range: Optional[Tuple[int, int]]) -> List[str]:
    if not qmmm_range:
        return []
    start, end = qmmm_range
    return [
        "%QMMM\n",
        # ORCA expects QMATOMS ranges to use its internal 0-based indexing.
        f"  QMATOMS {{{start}:{end}}}\n",
        "END\n",
        "END\n",
    ]


def _build_newgto_assignments(atoms, keywords: str, metal_basis: Optional[str],
                              control_args: Optional[Dict], qmmm_range: Optional[Tuple[int, int]],
                              skip_metals: bool = False) -> Dict[int, str]:
    if not metal_basis:
        return {}
    if _is_semiempirical_keywords(keywords):
        return {}

    config = control_args or {}
    basis = str(metal_basis).strip()
    if not basis:
        return {}

    qm_last = qmmm_range[1] if qmmm_range else len(atoms) - 1
    metal_indices = [idx for idx in _collect_metal_indices(atoms) if idx <= qm_last]
    if not metal_indices:
        return {}

    first_enabled = _is_enabled(config.get("first_coordination_sphere_metal_basisset"))
    second_enabled = _is_enabled(config.get("second_coordination_sphere_metal_basisset"))
    first_scale = _parse_float(config.get("first_coordination_sphere_scale"), 1.20)
    second_scale = _parse_float(config.get("second_coordination_sphere_scale"), 1.30)

    assignments: Dict[int, str] = {}
    if not skip_metals:
        assignments.update({idx: basis for idx in metal_indices})

    need_first = first_enabled or second_enabled
    first_indices = _first_coordination_sphere_indices(atoms, metal_indices, first_scale) if need_first else set()
    if first_enabled:
        for idx in first_indices:
            if idx <= qm_last:
                assignments[idx] = basis

    if second_enabled:
        second_indices = _second_coordination_sphere_indices(atoms, first_indices, metal_indices, second_scale)
        for idx in second_indices:
            if idx <= qm_last:
                assignments[idx] = basis

    return assignments


def _format_geometry_lines(atoms, keywords: str, metal_basis: Optional[str], control_args: Optional[Dict],
                           qmmm_range: Optional[Tuple[int, int]], basis_block_present: bool = False) -> List[str]:
    assignments = _build_newgto_assignments(
        atoms, keywords, metal_basis, control_args, qmmm_range, skip_metals=basis_block_present
    )
    qm_last = qmmm_range[1] if qmmm_range else None
    lines = []
    for idx, atom in enumerate(atoms):
        x, y, z = atom.position
        line = f"{atom.symbol:<3} {x:>14.8f} {y:>14.8f} {z:>14.8f}"
        if (qm_last is None or idx <= qm_last) and idx in assignments:
            line += f'   NewGTO "{assignments[idx]}" end'
        lines.append(line + "\n")
    return lines

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

def _fibonacci_sphere(n_samples):
    """
    Gleichmäßig verteilte Richtungen auf der Einheitssphäre (Fibonacci-Sampling).
    """
    n_samples = int(max(1, n_samples))
    if n_samples == 1:
        return [np.array([0.0, 0.0, 1.0])]

    offset = 2.0 / n_samples
    increment = np.pi * (3.0 - np.sqrt(5.0))
    dirs = []
    for i in range(n_samples):
        y = ((i * offset) - 1.0) + (offset / 2.0)
        r = np.sqrt(max(0.0, 1.0 - y * y))
        phi = i * increment
        x = np.cos(phi) * r
        z = np.sin(phi) * r
        dirs.append(np.array([x, y, z]))
    return dirs

def _find_max_clearance_direction(atoms, metal_index, distance, samples=800, clearance_scale=1.0):
    """
    Suche die Richtung, in der ein Punkt im Abstand 'distance' zur Metallposition
    den größten minimalen Abstand zu allen anderen Atomen besitzt.
    """
    if len(atoms) <= 1:
        return np.array([0.0, 0.0, 1.0]), np.inf

    directions = _fibonacci_sphere(samples)
    positions = atoms.positions
    other_indices = [i for i in range(len(atoms)) if i != metal_index]
    if not other_indices:
        return np.array([0.0, 0.0, 1.0]), np.inf

    other_positions = positions[other_indices]
    other_numbers = np.array([atoms[i].number for i in other_indices])
    other_radii = np.array([covalent_radii[num] for num in other_numbers])
    carbon_radius = covalent_radii[6]

    best_dir = np.array([0.0, 0.0, 1.0])
    best_clearance = -np.inf
    scale = float(max(0.0, clearance_scale))

    for direction in directions:
        direction = direction / np.linalg.norm(direction)
        candidate = direction * distance
        diffs = other_positions - candidate
        dists = np.linalg.norm(diffs, axis=1)
        if scale > 0.0:
            guard = scale * (other_radii + carbon_radius)
            clearances = dists - guard
        else:
            clearances = dists
        min_clearance = float(np.min(clearances))
        if min_clearance > best_clearance:
            best_clearance = min_clearance
            best_dir = direction

    return best_dir / np.linalg.norm(best_dir), best_clearance

def place_co2_general(complex_path, co2_path, out_path, distance=5.0, place_axis='z', mode='side-on',
                      perp_axis='y', optimize_direction=True, direction_samples=800, clearance_scale=1.0,
                      qm_count=None, qm_separator="$"):
    """
    Fügt CO2 an +place_axis in 'distance' Å an. Gibt zusätzlich CO2-Indizes im kombinierten System zurück.
    """
    comp = _read_xyz_robust(complex_path)
    co2 = _read_xyz_robust(co2_path)

    direction_samples = max(1, int(direction_samples) if direction_samples is not None else 1)
    clearance_scale = float(clearance_scale) if clearance_scale is not None else 1.0

    metal_index = None
    metal_shift = np.zeros(3)
    try:
        metal_index = detect_metal_index(comp)
    except Exception as exc:
        print(f"[WARN] Metallindex für Platzierung nicht bestimmt: {exc}")

    if metal_index is not None:
        metal_shift = comp.positions[metal_index].copy()
        comp.positions -= metal_shift

    axis_target = _axis_vector(place_axis)
    clearance = None
    place_dir = axis_target

    if optimize_direction and metal_index is not None:
        best_dir, clearance = _find_max_clearance_direction(
            comp, metal_index, distance,
            samples=direction_samples,
            clearance_scale=clearance_scale
        )
        place_dir = best_dir / np.linalg.norm(best_dir)
        print(f"[INFO] Gefundene freie Richtung: "
              f"({place_dir[0]:+.3f}, {place_dir[1]:+.3f}, {place_dir[2]:+.3f})")
        if clearance is not None and np.isfinite(clearance):
            print(f"[INFO] Optimierte CO₂-Richtung: minimale Clearance {clearance:.3f} Å")
            if clearance < 0.0:
                print("[WARN] Clearance < 0 Å – Platzierung kann zu Überlappungen führen.")
        else:
            print("[INFO] Optimierte CO₂-Richtung: keine Nachbaratome für Clearance-Berechnung.")
    elif optimize_direction and metal_index is None:
        print("[WARN] Platzierungsoptimierung nicht möglich – nutze statische place_axis.")

    axis, center, c_idx_local = _co2_axis_center_indices(co2)
    co2.positions -= center  # CO2 um sein C zentrieren

    if mode == "side-on":
        perp_target = _axis_vector(perp_axis)
        if optimize_direction and metal_index is not None:
            perp_projected = project_to_plane(perp_target, place_dir)
            if np.linalg.norm(perp_projected) < 1e-6:
                fallback = np.array([1.0, 0.0, 0.0])
                if abs(np.dot(fallback, place_dir)) > 0.9:
                    fallback = np.array([0.0, 1.0, 0.0])
                perp_projected = project_to_plane(fallback, place_dir)
            perp_target = perp_projected / np.linalg.norm(perp_projected)
        target = perp_target
    else:
        target = place_dir if (optimize_direction and metal_index is not None) else axis_target
    R = rot_from_vecs(axis, target)
    co2.positions = co2.positions @ R.T + distance * place_dir

    if metal_index is not None:
        comp.positions += metal_shift
        co2.positions += metal_shift

    combined = comp + co2
    n_comp = len(comp)
    n_co2 = len(co2)

    final_qm_count = None
    co2_indices = list(range(n_comp, n_comp + n_co2))

    if qm_count is not None:
        qm_core = int(max(0, min(int(qm_count), n_comp)))
        qm_indices = list(range(qm_core))
        mm_indices = list(range(qm_core, n_comp))
        new_order = qm_indices + co2_indices + mm_indices
        combined = combined[new_order]
        co2_indices = list(range(len(qm_indices), len(qm_indices) + n_co2))
        final_qm_count = len(qm_indices) + n_co2
    else:
        final_qm_count = None

    _write_xyz_with_separator(combined, out_path, final_qm_count, separator=qm_separator, insert_separator=False)
    print(f"[OK] wrote combined (complex + CO2) → {out_path}")

    # Indizes der CO2-Atome im kombinierten System
    co2_c_index_combined = co2_indices[c_idx_local]
    return out_path, combined, co2_indices, co2_c_index_combined, final_qm_count

# === ORCA Helpers ===
def _is_orca_calculation_complete(out_path):
    """Check if ORCA calculation is complete by looking for termination marker."""
    if not os.path.exists(out_path):
        return False
    try:
        # Check file size to avoid reading incomplete files
        if os.path.getsize(out_path) < 100:
            return False
        with open(out_path, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()
            return "ORCA TERMINATED NORMALLY" in content
    except Exception:
        return False


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
    # Inline NewGTO tags directly on atom lines to avoid duplicate blocks.
    return ""


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


def write_orca_sp_input_and_run(atoms, xyz_path, outdir, orca_keywords="GFN2-XTB SP ALPB(DMF)", additions="",
                                charge=-2, multiplicity=1, PAL=8, maxcore=2000, tag="calc",
                                metal_symbol=None, metal_basis=None, control_args=None,
                                qmmm_range: Optional[Tuple[int, int]] = None):
    os.makedirs(outdir, exist_ok=True)
    inp = os.path.join(outdir, f"{tag}.inp")
    out = os.path.join(outdir, f"{tag}.out")

    xyz_basename = os.path.basename(xyz_path)
    xyz_target = os.path.join(outdir, xyz_basename)

    keywords_line = _ensure_qmmm_keyword(orca_keywords, qmmm_range)
    basis_block = _format_basis_block(metal_symbol, metal_basis, keywords_line)

    lines = [f"! {keywords_line}\n"]
    if basis_block:
        block = basis_block if basis_block.endswith("\n") else basis_block + "\n"
        lines.append(block)
    lines.append(f"%maxcore {maxcore}\n")
    lines.append(f"%pal nprocs {PAL} end\n")
    if additions:
        add = additions if additions.endswith("\n") else additions + "\n"
        lines.append(add)

    lines.extend(_build_qmmm_block_lines(qmmm_range))
    lines.append(f"* xyz {charge} {multiplicity}\n")
    has_basis_block = bool(basis_block)
    lines.extend(_format_geometry_lines(
        atoms, keywords_line, metal_basis, control_args, qmmm_range, basis_block_present=has_basis_block
    ))
    lines.append("*\n")

    with open(inp, "w", encoding="utf-8") as f:
        f.writelines(lines)

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


def write_orca_input_and_run(atoms, xyz_path, metal_index, co2_c_index, start_distance, end_distance=1.7, steps=5,
                             orca_keywords="GFN2-XTB OPT ALPB(DMF)", additions="",
                             charge=-2, multiplicity=3, PAL=4, maxcore=2000,
                             metal_symbol=None, metal_basis=None, control_args=None,
                             qmmm_range: Optional[Tuple[int, int]] = None):
    additions_line = additions if additions else ""
    # For QM/MM workflows we keep the 0-based indexing convention throughout.
    i_orca = metal_index
    j_orca = co2_c_index

    keywords_line = _ensure_qmmm_keyword(orca_keywords, qmmm_range)
    basis_block = _format_basis_block(metal_symbol, metal_basis, keywords_line)
    sections = [f"! {keywords_line}\n"]
    if basis_block:
        block = basis_block if basis_block.endswith("\n") else basis_block + "\n"
        sections.append(block)
    sections.append(f"%maxcore {maxcore}\n")
    sections.append(f"%pal nprocs {PAL} end\n")
    if additions_line:
        add = additions_line if additions_line.endswith("\n") else additions_line + "\n"
        sections.append(add)
    sections.extend(_build_qmmm_block_lines(qmmm_range))
    sections.append("%geom\n")
    sections.append("  MaxIter 200\n")
    sections.append("  Scan\n")
    sections.append(f"    B  {i_orca}  {j_orca} = {start_distance:.2f}, {end_distance}, {steps}\n")
    sections.append("  end\n")
    sections.append("end\n")
    sections.append(f"* xyz {charge} {multiplicity}\n")
    sections.extend(_format_geometry_lines(
        atoms, keywords_line, metal_basis, control_args, qmmm_range, basis_block_present=bool(basis_block)
    ))
    sections.append("*\n")

    os.makedirs("relaxed_surface_scan", exist_ok=True)
    inp_path = os.path.join("relaxed_surface_scan", "scan.inp")
    xyz_target = os.path.join("relaxed_surface_scan", os.path.basename(xyz_path))

    with open(inp_path, "w", encoding="utf-8") as f:
        f.writelines(sections)
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
    # Use genfromtxt to handle NaN values properly
    data = np.genfromtxt(csv_path, delimiter=",", skip_header=1, filling_values=np.nan)
    ang = data[:, 0]
    rel = data[:, 2]

    # Filter out NaN values for plotting
    mask = ~np.isnan(rel)
    ang_valid = ang[mask]
    rel_valid = rel[mask]

    if len(ang_valid) == 0:
        print(f"[WARN] No valid data points to plot in {csv_path}")
        return

    plt.figure()
    plt.plot(ang_valid, rel_valid, marker='o')
    plt.xlabel("Rotation angle about z [deg]")
    plt.ylabel("ΔE [kcal/mol] (rel. to min)")
    plt.title("CO₂ orientation scan (fixed distance)")
    plt.grid()
    plt.tight_layout()
    plt.savefig(png_path)
    plt.close()
    print(f"[OK] Orientation plot saved: {png_path}")

def _calculate_single_angle(ang, base_atoms, co2_indices, charge, multiplicity, PAL, maxcore,
                           orca_keywords, additions, metal_symbol, metal_basis, control_args,
                           qmmm_range, qm_separator, recalc_mode=False):
    """
    Worker function to calculate energy for a single angle.
    Returns (angle_deg, energy_Eh, xyz_path) or (angle_deg, np.nan, None) on failure.
    """
    atoms = base_atoms.copy()
    R = Rz(ang)
    pos = atoms.positions.copy()
    pos[co2_indices] = (pos[co2_indices] @ R.T)  # rotiere nur CO2
    atoms.positions = pos

    ang_dir = os.path.join("orientation_scan", f"ang_{ang:03d}")
    os.makedirs(ang_dir, exist_ok=True)
    xyz_path = os.path.join(ang_dir, "structure.xyz")
    qm_count = (qmmm_range[1] + 1) if qmmm_range else None
    separator_line = (qm_separator or "$").strip() or "$"
    _write_xyz_with_separator(atoms, xyz_path, qm_count, separator=separator_line, insert_separator=False)

    out_path = os.path.join(ang_dir, "calc.out")

    # Check if calculation is already complete (recalc mode)
    if recalc_mode and _is_orca_calculation_complete(out_path):
        print(f"[recalc] Skipping angle {ang}° - already complete")
        try:
            E = parse_orca_energy(out_path)
            return (ang, E, xyz_path)
        except Exception as e:
            print(f"[WARN] Angle {ang}°: existing calc found but energy parse failed: {e}")
            # Fall through to re-run calculation

    out_path = write_orca_sp_input_and_run(
        atoms,
        xyz_path,
        ang_dir,
        orca_keywords=orca_keywords,
        additions=additions,
        charge=charge, multiplicity=multiplicity,
        PAL=PAL, maxcore=maxcore, tag="calc",
        metal_symbol=metal_symbol, metal_basis=metal_basis,
        control_args=control_args, qmmm_range=qmmm_range
    )
    try:
        E = parse_orca_energy(out_path)  # Hartree
        print(f"[INFO] angle {ang:3d}° → E = {E:.10f} Eh")
        return (ang, E, xyz_path)
    except Exception as e:
        print(f"[WARN] Angle {ang}°: energy parse failed: {e}")
        return (ang, np.nan, None)


# === Orientation scan (0–180°, SPs) ===
def orientation_scan_at_fixed_distance(base_atoms, combined_xyz_path, co2_indices, charge, multiplicity,
                                       PAL, maxcore, orca_keywords, additions,
                                       angle_step_deg=10, angle_range_deg=180,
                                       metal_symbol=None, metal_basis=None,
                                       control_args=None, qmmm_range=None, qm_separator="$",
                                       parallel=True, max_workers=None):
    """
    Dreht NUR die CO2-Atome um die z-Achse (durch den Ursprung) auf ihrer Position (z=const),
    macht für jeden Winkel eine SP-Rechnung und liefert die beste Geometrie zurück.

    Args:
        parallel: If True, run angle calculations in parallel (default: True)
        max_workers: Maximum number of parallel workers (default: auto-detect)
    """
    base_atoms = base_atoms.copy() if base_atoms is not None else _read_xyz_robust(combined_xyz_path)
    os.makedirs("orientation_scan", exist_ok=True)

    # Check recalc mode
    recalc_mode = os.environ.get("DELFIN_CO2_RECALC") == "1"
    if recalc_mode:
        print("[INFO] CO2 recalc mode enabled - skipping completed calculations")

    # Winkel-Liste 0..angle_range_deg inkl. Endpunkt
    angles = list(range(0, angle_range_deg + 1, angle_step_deg))
    results = []  # (angle_deg, energy_Eh, xyz_path)

    if parallel:
        # Parallel execution
        print(f"[INFO] Running orientation scan in parallel with {max_workers or 'auto'} workers")
        worker_func = partial(_calculate_single_angle,
                            base_atoms=base_atoms, co2_indices=co2_indices,
                            charge=charge, multiplicity=multiplicity,
                            PAL=PAL, maxcore=maxcore,
                            orca_keywords=orca_keywords, additions=additions,
                            metal_symbol=metal_symbol, metal_basis=metal_basis,
                            control_args=control_args, qmmm_range=qmmm_range,
                            qm_separator=qm_separator, recalc_mode=recalc_mode)

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_angle = {executor.submit(worker_func, ang): ang for ang in angles}
            for future in as_completed(future_to_angle):
                ang, E, xyz_path = future.result()
                results.append((ang, E, xyz_path))
    else:
        # Sequential execution (original behavior)
        print("[INFO] Running orientation scan sequentially")
        qm_count = (qmmm_range[1] + 1) if qmmm_range else None
        separator_line = (qm_separator or "$").strip() or "$"

        for ang in angles:
            ang, E, xyz_path = _calculate_single_angle(
                ang, base_atoms, co2_indices, charge, multiplicity, PAL, maxcore,
                orca_keywords, additions, metal_symbol, metal_basis, control_args,
                qmmm_range, qm_separator, recalc_mode
            )
            results.append((ang, E, xyz_path))

    # Sort results by angle
    results.sort(key=lambda x: x[0])

    # Find best geometry
    best_energy = None
    best_angle = None
    best_xyz = None

    for ang, E, xyz_path in results:
        if np.isfinite(E) and (best_energy is None or E < best_energy):
            best_energy = E
            best_angle = ang
            best_xyz = xyz_path

    # Save CSV
    csv_path = os.path.join("orientation_scan", "orientation_scan.csv")
    with open(csv_path, "w") as f:
        f.write("angle_deg,energy_Eh,relative_kcal_per_mol\n")
        finite_energies = [(ang, E) for ang, E, _ in results if np.isfinite(E)]
        if finite_energies:
            ref = min(E for _, E in finite_energies)
            for ang, E, _ in results:
                if np.isfinite(E):
                    kcal = E * 627.509
                    rel = (kcal - ref * 627.509)
                    f.write(f"{ang},{E},{rel}\n")
                else:
                    f.write(f"{ang},,\n")

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

    qmmm_initial = detect_qmmm_range_from_xyz(xyz_in) if os.path.exists(xyz_in) else None
    qm_atom_count = (qmmm_initial[1] + 1) if qmmm_initial else None
    qm_separator_detected = detect_qmmm_separator_from_xyz(xyz_in) if os.path.exists(xyz_in) else None
    has_dollar_separator = bool(qm_separator_detected and qm_separator_detected.strip() == "$")
    qm_separator = qm_separator_detected or "$"

    if has_dollar_separator:
        print('[WARN] "$" in Koordinaten erkannt – CPCM ist mit QM/MM nicht verwendbar.')
        current_model_raw = args.get("implicit_solvation_model")
        current_model_clean = _clean_str(current_model_raw, "CPCM")
        if current_model_clean.upper() == "CPCM":
            args["implicit_solvation_model"] = "ALPB"
            print("[INFO] implicit_solvation_model automatisch auf ALPB gesetzt.")
        elif current_model_clean:
            print(f"[INFO] Behalte implicit_solvation_model={current_model_clean}.")

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
    combined_path, combined_atoms, co2_indices, co2_c_idx, qm_atom_count = place_co2_general(
        aligned, co2_path, combined_out,
        distance=orientation_distance,
        place_axis=args.get("place_axis", "z"),
        mode=args.get("mode", "side-on"),
        perp_axis=args.get("perp_axis", "y"),
        optimize_direction=args.get("place_optimize", True),
        direction_samples=args.get("place_samples", 800),
        clearance_scale=args.get("place_clearance_scale", 1.0),
        qm_count=qm_atom_count,
        qm_separator=qm_separator
    )
    qmmm_range = (0, qm_atom_count - 1) if qm_atom_count is not None else None

    # --- 3) Orientation scan at fixed distance (optional) ---
    orientation_flag = args.get("perform_orientation_scan")
    perform_orientation_scan = True if orientation_flag is None else _is_enabled(orientation_flag)

    # Parallelization settings
    parallel_flag = args.get("parallel_orientation_scan")
    parallel_orientation = True if parallel_flag is None else _is_enabled(parallel_flag)
    max_workers = args.get("max_workers")
    if max_workers is not None and max_workers != "":
        try:
            max_workers = int(max_workers)
        except (ValueError, TypeError):
            print(f"[WARN] Invalid max_workers value '{max_workers}', using auto-detect")
            max_workers = None

    # Calculate safe max_workers to avoid exceeding total cores
    # Each worker will use PAL cores, so total_cores / PAL = safe_workers
    if parallel_orientation and max_workers is None:
        try:
            import multiprocessing
            total_cores = multiprocessing.cpu_count()
            safe_workers = max(1, total_cores // PAL)
            max_workers = safe_workers
            print(f"[INFO] Auto-detected parallelization: {max_workers} workers (each using {PAL} cores = {max_workers * PAL}/{total_cores} total cores)")
        except Exception:
            max_workers = 1
            print(f"[WARN] Could not detect CPU count, using 1 worker")
    elif parallel_orientation and max_workers is not None:
        # User specified max_workers - check if it's safe
        try:
            import multiprocessing
            total_cores = multiprocessing.cpu_count()
            required_cores = max_workers * PAL
            if required_cores > total_cores:
                safe_workers = max(1, total_cores // PAL)
                print(f"[WARN] max_workers={max_workers} × PAL={PAL} = {required_cores} cores exceeds available {total_cores} cores")
                print(f"[WARN] Reducing max_workers to {safe_workers} to stay within limits")
                max_workers = safe_workers
            else:
                print(f"[INFO] Using {max_workers} workers × {PAL} cores = {required_cores}/{total_cores} total cores")
        except Exception:
            print(f"[INFO] Using user-specified max_workers={max_workers}")

    if perform_orientation_scan:
        best_xyz_path, best_angle_deg, best_E = orientation_scan_at_fixed_distance(
            combined_atoms,
            combined_xyz_path=combined_path,
            co2_indices=co2_indices,
            charge=charge, multiplicity=multiplicity,
            PAL=PAL, maxcore=maxcore,
            orca_keywords=rot_orca_keywords,
            additions=additions,
            angle_step_deg=rot_step_deg,
            angle_range_deg=rot_range_deg,
            metal_symbol=metal_symbol,
            metal_basis=metal_basis,
            control_args=args,
            qmmm_range=qmmm_range,
            qm_separator=qm_separator,
            parallel=parallel_orientation,
            max_workers=max_workers
        )
    else:
        print("[INFO] Orientation scan disabled via CONTROL (perform_orientation_scan=false).")
        best_xyz_path = combined_path
        best_angle_deg = 0.0
        best_E = None

    # --- 4) Use best orientation to start relaxed distance scan ---
    atoms_best = _read_xyz_robust(best_xyz_path)
    metal_idx = detect_metal_index(atoms_best)
    metal_symbol_scan = metal_symbol or atoms_best[metal_idx].symbol
    start_distance = np.linalg.norm(atoms_best.positions[metal_idx] - atoms_best.positions[co2_c_idx])
    print(f"[INFO] Using best angle {best_angle_deg}° for distance scan. Start M–C = {start_distance:.3f} Å")

    write_orca_input_and_run(
        atoms_best,
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
        metal_basis=metal_basis,
        control_args=args,
        qmmm_range=qmmm_range
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
