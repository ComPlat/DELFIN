"""Results Overview: extract comparable quantities from finished calculations.

Reads only files that already exist on disk (archive/ or calc/ folders);
it never runs a program.  Every extractor returns ``None`` for a quantity
that is not present, so the table can show "not computed" instead of a
made-up value.

Data sources (verified against archive/403-429_403):
* ``DELFIN_Data.json`` — metadata, state optimizations (charge/multiplicity/
  hartree/converged), thermochemistry (Gibbs), orbitals (HOMO/LUMO/gap),
  ``computed.redox_potentials`` (E_ref, m1/m2/m3/final, gibbs per charge).
* ``CO2_coordination/relaxed_surface_scan/scan.relaxscanact.dat`` — two
  columns: distance (Angstrom), energy (hartree).  Used for the
  endergonic/exergonic comparison between the CO2-bound minimum and 4 A.
* PNGs under ``CO2_coordination/relaxed_surface_scan/*.png`` — collected
  with paths relative to the calculation folder.
"""

from __future__ import annotations

import glob
import json as _json
import os
import re
from typing import Any, Dict, List, Optional

# hartree -> kJ/mol
HARTREE_KJ_MOL = 2625.499638

DELFIN_DATA_NAME = "DELFIN_Data.json"
RSS_DIR = os.path.join("CO2_coordination", "relaxed_surface_scan")
RSS_DAT = os.path.join(RSS_DIR, "scan.relaxscanact.dat")

# The CO2-bound minimum lies around 2.0-2.1 A; the unbound reference is
# 3.5 A (4.0 A sits too close to the plain-minimum and mislabels endergonic
# cases, per user guidance 2026-10-30).
# Window around the CO2-bound minimum used for the legacy global-minimum
# column and the shape (saddle vs minimum) judgement.
RSS_MIN_WINDOW = (1.9, 2.3)
# RSS energetics: reference point (CO2 unbound side) vs CO2-bound point.
RSS_REF_DISTANCE = 3.5
RSS_BOUND_DISTANCE = 2.0
# Neighbourhood inspected to judge saddle vs true minimum around 2.0 A.
RSS_SHAPE_DISTANCES = (1.7, 1.8, 1.9, 2.1, 2.2, 2.3)

# Metric ids offered in the dashboard tab.  Values are extracted lazily.
METRICS: List[Dict[str, str]] = [
    {"id": "redox_E_ox", "label": "Oxidation potential E_ox (V vs E_ref)"},
    {"id": "redox_E_red", "label": "Reduction potential E_red (V vs E_ref)"},
    {"id": "redox_E_red_2", "label": "Reduction potential E_red2 (V vs E_ref)"},
    {"id": "redox_E_ref", "label": "Reference potential E_ref (V)"},
    {"id": "redox_method", "label": "Redox method (m1/m2/m3 flags)"},
    {"id": "mult_ground", "label": "Multiplicity ground state (charge/mult bs)"},
    {"id": "mult_oxidized", "label": "Multiplicity oxidized (charge/mult bs)"},
    {"id": "mult_reduced", "label": "Multiplicity reduced (charge/mult bs)"},
    {"id": "mult_reduced2", "label": "Multiplicity 2x reduced (charge/mult bs)"},
    {"id": "rss_delta_kj", "label": "RSS: E(2.0 A) - E(3.5 A) (kJ/mol)"},
    {"id": "rss_direction", "label": "RSS: endergonic/exergonic (2.0 vs 3.5 A)"},
    {"id": "rss_shape", "label": "RSS: minimum type (global/local/none)"},
    {"id": "rss_min_distance", "label": "RSS: distance at sweep minimum (A)"},
    {"id": "rss_min_delta_kj", "label": "RSS: E(min) - E(3.5 A) (kJ/mol)"},
    {"id": "rss_tight_min_distance", "label": "RSS: min in 1.9-2.3 A window (A)"},
    {"id": "rss_tight_delta_kj", "label": "RSS: E(tight min) - E(3.5 A) (kJ/mol)"},
    {"id": "gibbs_ground", "label": "Gibbs energy ground state (hartree)"},
    {"id": "gap_ground", "label": "HOMO-LUMO gap ground state (eV)"},
    {"id": "converged_flags", "label": "Convergence of all state optimizations"},
    {"id": "method", "label": "Method (functional/basis/solvent/charge)"},
]


def list_calculations(roots: List[str]) -> List[Dict[str, str]]:
    """Find calculation folders containing DELFIN_Data.json under *roots*.

    Returns entries with ``path`` (absolute), ``label`` and ``root``.  Both
    flat calc folders and archive ``<range>_<name>`` folders are handled.
    """
    found: Dict[str, Dict[str, str]] = {}
    for root in roots:
        root = os.path.abspath(root)
        if not os.path.isdir(root):
            continue
        for dirpath, _dirnames, filenames in os.walk(root):
            if DELFIN_DATA_NAME in filenames and os.path.join(
                dirpath, RSS_DIR
            ) != dirpath:
                rel = os.path.relpath(dirpath, root)
                label = os.path.basename(dirpath)
                found.setdefault(dirpath, {
                    "path": dirpath,
                    "label": label,
                    "rel": rel,
                    "root": root,
                })
    # Walk order is arbitrary; sort by label for a stable table.
    return sorted(found.values(), key=lambda e: (e["rel"], e["label"]))


def _load_delfin_data(calc_dir: str) -> Optional[dict]:
    path = os.path.join(calc_dir, DELFIN_DATA_NAME)
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return _json.load(fh)
    except (OSError, _json.JSONDecodeError):
        return None


def _state_opt(data: dict, *keys: str) -> Optional[dict]:
    node: Any = data
    for key in keys:
        if not isinstance(node, dict):
            return None
        node = node.get(key)
    if isinstance(node, dict):
        opt = node.get("optimization")
        if isinstance(opt, dict):
            return opt
    return None


def _format_mult(opt: Optional[dict], broken_sym: Optional[str] = None
                 ) -> Optional[str]:
    if not opt:
        return None
    charge = opt.get("charge")
    mult = opt.get("multiplicity")
    if charge is None or mult is None:
        return None
    text = f"{charge}/{mult}"
    if broken_sym:
        # Compact suffix per user request: "2/3 bs 1,1" — only the values,
        # no state prefix (the state is the column).
        m = re.search(r"brokensym\s+(-?\d+(?:\s*,\s*-?\d+)*)",
                      broken_sym, re.IGNORECASE)
        if m and m.group(1).strip():
            text += f" bs {m.group(1).strip()}"
    return text


def _redox(data: dict) -> Optional[dict]:
    computed = data.get("computed") if isinstance(data, dict) else None
    if isinstance(computed, dict):
        rp = computed.get("redox_potentials")
        if isinstance(rp, dict):
            return rp
    return None


def parse_rss_dat(calc_dir: str) -> List[tuple]:
    """Parse ``scan.relaxscanact.dat`` into [(distance_A, energy_hartree)]."""
    path = os.path.join(calc_dir, RSS_DAT)
    points: List[tuple] = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    dist = float(parts[0])
                    energy = float(parts[1])
                except ValueError:
                    continue
                points.append((dist, energy))
    except OSError:
        return []
    return points


_SMILES_METALS = {
    "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Al", "Ga", "In", "Sn", "Pb", "Bi", "Ce", "Sm", "Eu", "Gd",
}


def _smiles_from_3d(calc_dir: str):
    """Derive a kekulized SMILES from initial.xyz (metal-stripped).

    Fallback for older calcs whose DELFIN_Data.json carries no
    input.smiles (46 of 115 do).  The heavy+H graph minus metals is
      built as an RDKit molecule with the XYZ coordinates as conformer,
      then rdDetermineBonds assigns connectivity and orders; a short
      charge sweep finds the fragment charge that admits a valid
      ordering (proven on 106-138_106: charge -2 yields the fully
      kekulized complex; 106-138_109 needs -3, 106-138_108..113 need
      -2).  None if RDKit is absent or no ordering is valid for any
      charge in the sweep.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError:
        return None
    atoms = _parse_xyz(calc_dir)
    if not atoms:
        return None
    sym = [a[0] for a in atoms]
    keep = [i for i, s in enumerate(sym) if s not in _SMILES_METALS]
    if not keep:
        return None
    rw = Chem.RWMol()
    for i in keep:
        rw.AddAtom(Chem.Atom(sym[i]))
    mol = rw.GetMol()
    conf = Chem.Conformer(len(keep))
    for n, i in enumerate(keep):
        conf.SetAtomPosition(n, (float(atoms[i][1]), float(atoms[i][2]),
                                 float(atoms[i][3])))
    mol.AddConformer(conf)
    for charge in (0, 1, -1, 2, -2, 3, -3, 4, -4):
        try:
            trial = Chem.RWMol(mol)
            rdDetermineBonds.DetermineConnectivity(trial)
            rdDetermineBonds.DetermineBondOrders(trial, charge=charge)
            return Chem.MolToSmiles(trial, kekuleSmiles=True)
        except Exception:  # noqa: BLE001  (charge mismatch, valence)
            continue
    return None


def get_smiles(calc_dir: str):
    """SMILES for a calc: DELFIN_Data.json ``input.smiles`` first, then
    derived from the 3D structure (None only if neither works)."""
    path = os.path.join(calc_dir, DELFIN_DATA_NAME)
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            data = _json.load(fh)
    except (OSError, _json.JSONDecodeError):
        data = None
    value = (data or {}).get("input", {}).get("smiles")
    if value:
        return str(value)
    return _smiles_from_3d(calc_dir)


def smiles_svg_bytes(smiles: str, size: int = 260):
    """Column-2 drawing as VECTOR SVG, drawn exactly like the calc's
    SMILES.png (delfin_docx_report.py: RemoveHs, Normalizer, CoordGen
    with scaling 1.5, bondLineWidth 3, multipleBondOffset 0.15,
    stereo annotations, kekulize).  None on failure.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D, rdMolDraw2DDrawer
    except ImportError:
        return None
    try:
        from rdkit.Chem import rdCoordGen
    except ImportError:
        rdCoordGen = None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = Chem.RemoveHs(mol, updateExplicitCount=True)
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize
            mol = rdMolStandardize.Normalizer().normalize(mol)
        except Exception:
            pass
        done = False
        if rdCoordGen is not None:
            try:
                params = rdCoordGen.CoordGenParams()
                params.SetCoordgenScaling(1.5)
                params.SetMinimizerPrecision(1e-4)
                done = rdCoordGen.AddCoords(mol, params)
            except Exception:
                done = False
        if not done:
            AllChem.Compute2DCoords(
                mol, canonOrient=True, clearConfs=True, coordMap={},
                nFlipsPerSample=3, nSample=100, sampleSeed=42,
                permuteDeg4Nodes=True, bondLength=1.5)
            try:
                from rdkit.Chem import rdDepictor
                rdDepictor.StraightenDepiction(mol)
            except Exception:
                pass
        drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        opts.bondLineWidth = 3
        opts.multipleBondOffset = 0.15
        opts.addStereoAnnotation = True
        opts.addAtomIndices = False
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, kekulize=True)
        drawer.FinishDrawing()
        return drawer.GetDrawingText().encode("utf-8")
    except Exception:
        return None


def smiles_png_bytes(smiles: str, size: int = 260, fmt: str = "png"):
    """2D structure image via RDKit (PNG or vector PDF bytes); None if
    RDKit is absent or the SMILES cannot be parsed."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if fmt == "png":
        # Plain RDKit drawing, the standard the calculators browser
        # uses (tab_calculations_browser.py: Compute2DCoords + draw).
        try:
            from rdkit.Chem import rdDepictor
            rdDepictor.Compute2DCoords(mol)
        except Exception:
            pass
        drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    import numpy
    # Vector path: EXACTLY the SMILES.png pipeline of
    # delfin_docx_report.py (RemoveHs -> Normalizer -> CoordGen 1.5 /
    # 1e-4, fallback Compute2DCoords + StraightenDepiction), only the
    # RENDERER differs: the same coordinates go into the colored
    # matplotlib vector renderer (SVG-style look, vector PDF bytes).
    mol = Chem.RemoveHs(mol, updateExplicitCount=True)
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
        mol = rdMolStandardize.Normalizer().normalize(mol)
    except Exception:
        pass
    coordgen_done = False
    try:
        from rdkit.Chem import rdCoordGen
        params = rdCoordGen.CoordGenParams()
        params.SetCoordgenScaling(1.5)
        params.SetMinimizerPrecision(1e-4)
        coordgen_done = rdCoordGen.AddCoords(mol, params)
    except Exception:
        coordgen_done = False
    if not coordgen_done:
        AllChem.Compute2DCoords(
            mol, canonOrient=True, clearConfs=True, coordMap={},
            nFlipsPerSample=3, nSample=100, sampleSeed=42,
            permuteDeg4Nodes=True, bondLength=1.5)
        try:
            from rdkit.Chem import rdDepictor
            rdDepictor.StraightenDepiction(mol)
        except Exception:
            pass
    # SMILES.png kekulizes (Draw.MolToImage(..., kekulize=True)): the
    # alternating double bonds of the aromatic rings must appear.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass
    conf = mol.GetConformer()
    pos = numpy.array([[conf.GetAtomPosition(i).x,
                        conf.GetAtomPosition(i).y]
                       for i in range(mol.GetNumAtoms())])
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
              int(b.GetBondTypeAsDouble())) for b in mol.GetBonds()]
    orders = {(i, j): o for i, j, o in bonds if o >= 2}
    pairs = [(i, j) for i, j, _ in bonds]
    return _render_skeleton(sym, pos, pairs, size,
                            bond_orders=orders, fmt=fmt, colored=True)




# Covalent radii (Angstrom, Cordero 2008) for bond detection.
_XYZ_COVALENT = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "B": 0.84, "Si": 1.11, "Br": 1.20, "I": 1.39,
    "Na": 1.66, "K": 2.03, "Li": 1.28, "Mg": 1.41, "Ca": 1.76,
    "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32,
    "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Mo": 1.54,
    "Pd": 1.39, "Ag": 1.45, "Sn": 1.39, "Pt": 1.36, "Au": 1.36,
    "Hg": 1.32, "Pb": 1.46,
}
_XYZ_COVALENT_DEFAULT = 1.20
#: A bond forms when d(i,j) <= BOND_FACTOR * (r_i + r_j).
BOND_FACTOR = 1.15

# Display radii (Angstrom) per element, roughly ball-and-stick scaled.
_XYZ_RADII = {
    "H": 0.32, "C": 0.55, "N": 0.56, "O": 0.53, "F": 0.47, "P": 0.80,
    "S": 0.80, "Cl": 0.79, "B": 0.70, "Si": 0.85, "Br": 0.94, "I": 1.05,
    "Na": 0.90, "K": 1.20, "Li": 0.70, "Mg": 0.90, "Ca": 1.10,
}
# Default radii for transition metals and any element not listed.
_XYZ_RADIUS_DEFAULT = 0.95

# CPK-like colours (hex) for the common elements.
_XYZ_COLORS = {
    "H": "#f0f0f0", "C": "#404040", "N": "#3050f8", "O": "#ff2d2d",
    "F": "#90e050", "P": "#ff8000", "S": "#ffff30", "Cl": "#1ff01f",
    "B": "#ffb5b5", "Si": "#f0c8a0", "Br": "#a62929", "I": "#940094",
    "Na": "#ab5cf2", "K": "#8f40d4", "Li": "#cc80ff", "Mg": "#8aff00",
    "Ca": "#3dff00", "Ti": "#bfc2c7", "V": "#a6a6ab", "Cr": "#8a99c7",
    "Mn": "#9c7ac7", "Fe": "#e06633", "Co": "#f090a0", "Ni": "#50d050",
    "Cu": "#c88033", "Zn": "#7d80b0", "Si": "#f0c8a0", "Mo": "#54b5b5",
    "Pd": "#cdA0a0", "Ag": "#c0c0c0", "Sn": "#668f8f", "Pt": "#d0d0e0",
    "Au": "#ffd123", "Hg": "#b8b8d0", "Pb": "#575961",
}
_XYZ_COLOR_DEFAULT = "#d0a070"  # muted metal tone


def _parse_xyz(calc_dir: str):
    """Read initial.xyz -> [(symbol, x, y, z), ...] or None.

    Line 0 is a free-text title; line 1 the atom count (standard XYZ).
    ORCA's initial.xyz variant instead starts coordinates at line 1 with
    the count in the title line 0, so parse robustly: take any line whose
    first token is a symbol and next three fields are numbers.
    """
    xyz = os.path.join(calc_dir, "initial.xyz")
    try:
        with open(xyz, "r", encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
    except OSError:
        return None
    atoms = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) >= 4:
            try:
                atoms.append((parts[0].capitalize(), float(parts[1]),
                              float(parts[2]), float(parts[3])))
                continue
            except ValueError:
                pass
        if atoms:
            break  # coordinate block ended
    return atoms or None


_XYZ_MASS = {
    "H": 1.008, "He": 4.003, "Li": 6.94, "Be": 9.012, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.18,
    "Na": 22.99, "Mg": 24.305, "Al": 26.982, "Si": 28.085,
    "P": 30.974, "S": 32.06, "Cl": 35.45, "K": 39.098, "Ca": 40.078,
    "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
    "Zn": 65.38, "Br": 79.904, "I": 126.904,
}


def _rotate_inertia(sym, pts):
    """Mass-weighted principal inertia rotation: largest cross-section on z.

    The camera looks along z; the SMALLEST extent axis must become z so the
    largest cross-section faces the viewer.  eigh() sorts ascending, so
    column 0 (smallest eigenvalue) becomes z.
    """
    import numpy
    pts = numpy.asarray(pts, dtype=float)
    pts = pts - pts.mean(axis=0)  # geometry centre is enough for orientation
    weights = numpy.array([_XYZ_MASS.get(s, 12.011) for s in sym])
    inertia = (pts * weights[:, None]).T @ pts
    _eigvals, eigvecs = numpy.linalg.eigh(inertia)
    if numpy.linalg.det(eigvecs) < 0:  # keep a right-handed frame
        eigvecs[:, 2] = -eigvecs[:, 2]
    rotated = pts @ eigvecs
    rotated = rotated[:, [1, 2, 0]]
    # Centre on the projected bounding box, not the mean: the mean can
    # sit off-centre after an asymmetric rotation and then the frame
    # clips atoms on one side (user report 430-450_432).
    rotated[:, 0] -= 0.5 * (rotated[:, 0].min() + rotated[:, 0].max())
    rotated[:, 1] -= 0.5 * (rotated[:, 1].min() + rotated[:, 1].max())
    return rotated


def initial_xyz_png_bytes(calc_dir: str, size: int = 320, show_h=True,
                         bond_orders=None, fmt="png"):
    """Skeletal formula of ``initial.xyz`` (PNG bytes), or None.

    The structure is rotated onto its principal inertia axes so the largest
    cross-section faces the camera (most atoms visible), then depth-sorted
    and drawn with matplotlib (Agg) as a skeletal wireframe (bonds as black
    lines, atom symbols in gaps, H suppressed) in RDKit-2D style.  Returns None when no readable initial.xyz exists or
    numpy/matplotlib are unavailable.  It never runs a program.
    """
    atoms = _parse_xyz(calc_dir)
    if not atoms:
        return None
    try:
        import numpy
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        return None

    sym = [a[0] for a in atoms]
    pts = numpy.array([[a[1], a[2], a[3]] for a in atoms])
    try:
        rotated = _rotate_inertia(sym, pts)
    except numpy.linalg.LinAlgError:
        return None

    # Bond detection: covalent radii sum x BOND_FACTOR (Cordero 2008).
    cov = [_XYZ_COVALENT.get(s, _XYZ_COVALENT_DEFAULT) for s in sym]
    bonds = _detect_bonds(sym, cov, rotated)
    if bond_orders is None:
        # First source: the calc's own SMILES (correct even for
        # metal-coordinated heteroaromatics, calc 428).
        from delfin.analysis._smiles_orders import smiles_bond_orders
        bond_orders = smiles_bond_orders(sym, bonds, calc_dir)
    bond_orders = bond_orders if bond_orders is not None \
        else _rdkit_bond_orders(sym, bonds, rotated)

    # Skeletal convention, INDEPENDENT of the checkbox: aliphatic H
    # (C-H) are never drawn in the skeleton; protons on heteroatoms
    # (OH, NH2, COOH, SH...) always are.  The "Show hydrogen atoms"
    # checkbox steers only the ball-stick view (optimized 3D column).
    _hetero = {"O", "N", "S", "P", "Se", "B", "Si", "F", "Cl", "Br", "I"}
    hetero_h_count = {}
    for i, j in bonds:
        for h, other in ((i, j), (j, i)):
            if sym[h] == "H" and sym[other] in _hetero:
                hetero_h_count[other] = hetero_h_count.get(other, 0) + 1
    # All H atoms drop out of the geometry; their count merges into the
    # heteroatom LABEL ("NH2", "OH", "SH").  This avoids dangling H
    # bonds to nowhere and the wrong "proton on the carbonyl O" look.
    _h_suffix = {0: "", 1: "H", 2: "H2", 3: "H3"}
    label = {i: (s + _h_suffix.get(hetero_h_count.get(i, 0), "")
                 if s in _hetero else s)
             for i, s in enumerate(sym)}
    # Hydroxyl reading direction: the label sits on the side of the O
    # away from its single heavy partner.  C right of O -> text left
    # of O -> "HO"; C left of O -> "OH".
    _odeg = [0] * len(sym)
    for i, j in bonds:
        _odeg[i] += 1
        _odeg[j] += 1
    for i, s in enumerate(sym):
        if s == "O" and hetero_h_count.get(i, 0) == 1 and _odeg[i] == 2:
            _pc = next((j if k == i else k) for k, j in bonds
                       if i in (k, j) and sym[j if k == i else k] != "H")
            label[i] = ("HO" if rotated[i][0] < rotated[_pc][0]
                        else "OH")
    keep = [i for i, s in enumerate(sym) if s != "H"]
    keepset = set(keep)
    remap = {old_i: new_i for new_i, old_i in enumerate(keep)}
    sym_vis = [label[i] for i in keep]
    rotated = rotated[keep]
    bonds_vis = [(remap[i], remap[j]) for i, j in bonds
                 if i in keepset and j in keepset]
    # The upgrade indices must follow the remap or the double-bond
    # lines point at hidden atoms and silently disappear.
    bond_orders = {(remap[i], remap[j]): o
                   for (i, j), o in bond_orders.items()
                   if i in keepset and j in keepset}
    return _render_skeleton(sym_vis, rotated, bonds_vis, size,
                            bond_orders=bond_orders, fmt=fmt)


def _detect_bonds(sym, cov, pts):
    """Bond pairs (i, j) with d <= BOND_FACTOR * (cov_i + cov_j).

    A fixed cutoff joined everything within ~2.4 A and drew far too
    many bonds; per-element radii keep H-H and ring gaps unbound.
    A per-pair sanity cap removes close NON-bonded contacts that pass
    the generic factor (stacked rings, crowded metal spheres).
    """
    import numpy
    cov = numpy.asarray(cov, dtype=float)
    n = len(sym)
    _MAX = {("C", "C"): 1.80, ("C", "H"): 1.30, ("H", "H"): 1.05,
            ("N", "H"): 1.30, ("O", "H"): 1.25}
    bonds = []
    for i in range(n):
        for j in range(i + 1, n):
            d = pts[i] - pts[j]
            dist = float(numpy.sqrt(d @ d))
            if dist > BOND_FACTOR * (cov[i] + cov[j]):
                continue
            key = tuple(sorted((sym[i], sym[j])))
            if dist <= _MAX.get(key, 2.4):
                bonds.append((i, j))
    return bonds


def _render_ball_stick(sym, rotated, bonds, size, atom_scale=0.8,
                       bond_lw=1.4, fmt="png"):
    """Draw a depth-sorted ball-and-stick projection (PNG bytes).

    ``atom_scale`` scales the display radii (0.8 = classic filling),
    ``bond_lw`` the bond line width.  The dashboard 3Dmol viewer style
    corresponds to small spheres on thick sticks (molecule_viewer.py:
    sphere scale 0.28 of vdW, bond radius 0.11 A).
    """
    import numpy
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from io import BytesIO

    radii = numpy.array([_XYZ_RADII.get(s, _XYZ_RADIUS_DEFAULT) for s in sym])
    # Headroom framing: cover outermost atoms PLUS their display radii,
    # else carbonyl O (428) or balls (434) clip at the frame edge.
    disp = 0.55 * radii
    cx, cy = (0.5 * (rotated[:, 0].max() + rotated[:, 0].min()),
              0.5 * (rotated[:, 1].max() + rotated[:, 1].min()))
    rmax = float(disp.max())
    scale = size / (2.0 * (max(abs(rotated[:, 0] - cx).max() + rmax,
                               abs(rotated[:, 1] - cy).max() + rmax)))
    rotated = rotated - numpy.array([cx, cy] + [0.0]
                            * (rotated.shape[1] - 2))
    fig = plt.figure(figsize=(size / 100, size / 100), dpi=100)
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    half = 0.5 * size / scale  # data half-width in Angstrom
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_aspect("equal")

    # Bonds behind atoms, depth-sorted by their average z (far first).
    for i, j in sorted(bonds, key=lambda b: rotated[b[0]][2] + rotated[b[1]][2]):
        ax.plot([rotated[i][0], rotated[j][0]],
                [rotated[i][1], rotated[j][1]],
                color="#909090", lw=bond_lw, zorder=1, solid_capstyle="round")

    # Atoms: near atoms (large z) last, on top.
    for idx in numpy.argsort(rotated[:, 2]):
        x, y = rotated[idx][0], rotated[idx][1]
        r_px = radii[idx] * scale * atom_scale
        r_data = r_px / size * (2 * half)
        ax.add_patch(mpatches.Circle(
            (x, y), r_data,
            facecolor=_XYZ_COLORS.get(sym[idx], _XYZ_COLOR_DEFAULT),
            edgecolor="white", linewidth=0.6,
            zorder=2 + rotated[idx][2] / (2.0 * (abs(rotated[:, 2]).max()
                                                 + 1e-6))))
    from io import BytesIO
    buf = BytesIO()
    fig.savefig(buf, format=fmt, dpi=100)
    plt.close(fig)
    return buf.getvalue()


def molecular_formula(calc_dir: str):
    """Hill-ordered molecular formula of initial.xyz, e.g. C34H24N4Ni.

    C first, H second, remaining elements alphabetically; None when no
    readable initial.xyz exists.
    """
    import collections
    atoms = _parse_xyz(calc_dir)
    if not atoms:
        return None
    counts = collections.Counter(a[0] for a in atoms)
    parts = []
    for s in ("C", "H") + tuple(sorted(k for k in counts
                                      if k not in ("C", "H"))):
        n = counts.pop(s, 0)
        if n:
            parts.append(s if n == 1 else f"{s}{n}")
    return "".join(parts)


def formula_to_html(formula: str) -> str:
    """Molecular formula with true sub/superscripts for the browser.

    Trailing digits become <sub>2</sub>, a trailing charge like 2+ or
    2- becomes <sup>2+</sup>: "[Fe(bpy)3]2+" -> "[Fe(bpy)<sub>3</sub>]"
    with a <sup>2+</sup> suffix.  Formulas without indices pass through.
    """
    import re
    if not formula:
        return formula or ""
    out = re.sub(r"(?<=[A-Za-z\)])((?:\d+)+)", r"<sub>\1</sub>", formula)
    m = re.search(r"(\d*[+-])$", out)
    if m and m.group(1):
        out = out[:m.start()] + "<sup>" + m.group(1) + "</sup>"
    return out


def formula_to_latex(formula: str) -> str:
    """Molecular formula in LaTeX math: digits become _{n}, a trailing
    charge like 2+ becomes $^{2+}$."""
    import re
    if not formula:
        return formula or ""
    m = re.search(r"(\d*[+-])$", formula)
    body, charge = formula, ""
    if m and m.group(1):
        body, charge = formula[:m.start()], m.group(1)
    out = re.sub(r"(?<=[A-Za-z\)])((?:\d+)+)", r"$_{\1}$", body)
    if charge:
        out += "$^{" + charge + "}$"
    return out


def formula_to_richtext(formula: str):
    """Molecular formula as openpyxl rich text: list of (text, vertAlign).

    vertAlign 'subscript' for the index runs, None elsewhere; charges
    become 'superscript'.  Suitable for openpyxl CellRichText.
    """
    import re
    if not formula:
        return [(formula or "", None)]
    parts, i = [], 0
    for m in re.finditer(r"(?<=[A-Za-z\)])((?:\d+)+)|(\d*[+-]$)", formula):
        if m.start() > i:
            parts.append((formula[i:m.start()], None))
        if m.group(1):
            parts.append((m.group(1), "subscript"))
        else:
            parts.append((m.group(2), "superscript"))
        i = m.end()
    if i < len(formula):
        parts.append((formula[i:], None))
    return parts or [(formula, None)]


def _rdkit_bond_orders(sym, bonds, positions):
    """Bond orders via RDKit: metal-free graph, DetermineBondOrders,
    kekulize; mapped back to the original indices.

    Metals and their bonds are removed (the ligands become fragments);
    RDKit assigns orders from valences (charged fragments allowed, total
    charge 0); Kekulization produces the alternating aromatic pattern.
    Falls back to the valence heuristic on any failure.  Returns
    {(i, j): 2|3} over ORIGINAL indices; metal bonds are never upgraded.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError:
        return _infer_bond_orders(sym, bonds)
    metal = set(_XYZ_MASS) - {"H", "B", "C", "N", "O", "F", "Si", "P",
                              "S", "Cl", "Br", "I"}
    keep = [i for i, s in enumerate(sym) if s not in metal]
    if len(keep) < 2:
        return _infer_bond_orders(sym, bonds)
    kset = set(keep)
    remap = {old: new for new, old in enumerate(keep)}
    rw = Chem.RWMol()
    for i in keep:
        rw.AddAtom(Chem.Atom(sym[i]))
    for i, j in bonds:
        if i in kset and j in kset:
            rw.AddBond(remap[i], remap[j], Chem.BondType.SINGLE)
    conf = Chem.Conformer(rw.GetNumAtoms())
    for n, o in enumerate(keep):
        conf.SetAtomPosition(n, tuple(float(x) for x in positions[o]))
    conf.SetId(0)
    rw.AddConformer(conf, assignId=True)
    mol = rw.GetMol()
    # The complex may be an anion/cation (e.g. 428: Co complex -1):
    # try small charges until one fits, else fall back.
    orders_ok = False
    for charge in (1, -1, 0, 2, -2, 3, -3, 4, -4):
        try:
            probe = Chem.Mol(mol)
            rdDetermineBonds.DetermineBondOrders(
                probe, charge=charge, allowChargedFragments=True)
            Chem.SanitizeMol(probe)
            Chem.Kekulize(probe)
            mol = probe
            orders_ok = True
            break
        except Exception:
            continue
    if not orders_ok:
        return _infer_bond_orders(sym, bonds)
    upgrades = {}
    for b in range(mol.GetNumBonds()):
        bond = mol.GetBondWithIdx(b)
        name = str(bond.GetBondType())
        if name not in ("DOUBLE", "TRIPLE"):
            continue
        i = keep[bond.GetBeginAtomIdx()]
        j = keep[bond.GetEndAtomIdx()]
        upgrades[(i, j)] = 2 if name == "DOUBLE" else 3
    # Deprotonated pyrrole-type N (bonded to a metal plus exactly two
    # ring carbons, no H): the metal bond consumes its lone pair, so
    # its ring double bonds sit BETWEEN CARBONS - never C=N (user,
    # 428 pyrrole ligand).  Drop N doubles at such N; the adjacent C
    # keeps the alternation from its other edge.
    # SIX-membered rings are pyridine-type: the coordinated ring N
    # KEEPS its C=N double bond (removing it left pyridine rings with
    # 2 instead of 3 doubles, user report 106-138_108..113).  The drop
    # therefore applies only when the N sits in a FIVE-membered ring
    # (pyrrole/pyrazolate).  Ring membership is decided on the heavy-
    # atom graph by walking 3 and 4 steps back to the N itself.
    _METALS = set(_XYZ_MASS) - {"H", "B", "C", "N", "O", "F", "Si",
                                "P", "S", "Cl", "Br", "I"}
    _neigh = {}
    for a, b in bonds:
        _neigh.setdefault(a, []).append(b)
        _neigh.setdefault(b, []).append(a)

    def _in_ring_of_size(atom, size, hops):
        """True if a cycle of length `size` returns to `atom` within
        `hops` heavy-atom steps (H and metals excluded from the walk)."""
        targets = {size}
        frontier = [(atom, frozenset())]
        for _ in range(hops):
            nxt = []
            for node, path in frontier:
                for k in _neigh.get(node, ()):
                    if sym[k] in _METALS or sym[k] == "H":
                        continue
                    if k in path:
                        continue
                    if k == atom and len(path) + 1 == size:
                        return True
                    nxt.append((k, path | {node}))
            frontier = nxt
        return False

    drop = set()
    for atom, ns in _neigh.items():
        # Target the metal-coordinated deprotonated N: it must have a
        # metal neighbour, no H, and exactly two ring-carbon neighbours,
        # and sit in a five-membered ring (pyrrolate); a pyridine-type
        # N in a six-membered ring keeps its C=N double bond.
        if sym[atom] != "N" or not any(sym[k] in _METALS for k in ns):
            continue
        hs = [k for k in ns if sym[k] == "H"]
        cs = [k for k in ns if sym[k] == "C"]
        if len(hs) == 0 and len(cs) == 2 \
                and _in_ring_of_size(atom, 5, 4):
            for k in ns:
                if sym[k] != "H":
                    drop.add(tuple(sorted((atom, k))))
    upgrades = {pair: o for pair, o in upgrades.items()
                if pair not in drop}
    # DetermineBondOrders can hand a coordinating heteroatom TWO double
    # bonds (pyridine N drawn as an N+ fragment: N= and =N).  A valid
    # kekule structure is a MATCHING: every atom sits in at most one
    # double bond.  Greedy shortest-first matching: accept a candidate
    # if neither endpoint is taken yet; this keeps the true C=N/C=C
    # (shorter) bonds and repairs the alternation (user report 428).
    import numpy as _np
    def _len(pair):
        d = positions[pair[0]] - positions[pair[1]]
        return float(_np.sqrt(d @ d))
    # Metal-ligand bonds are dative: they consume a valence of the
    # ligand atom.  DetermineBondOrders sees only the metal-stripped
    # fragment, so a ligand atom that carries the metal in the FULL
    # graph got its order assigned WITHOUT that bond.  For a thiolate
    # S (metal + ring-C) that produced a quinoid "C=S" to the ring
    # (user report class C, 11-14).  But a blanket drop would also
    # strip the pyridine C=N the metal-coordinated ring N must keep
    # (class A fix).  The precise rule: drop an upgrade only when the
    # atom would become OVERVALENT in the full graph - its heavy-atom
    # degree there (metal bonds included) already reaches its normal
    # valence, so a double bond has no valence left to use.
    _MET = set(_XYZ_MASS) - {"H", "B", "C", "N", "O", "F", "Si", "P",
                             "S", "Cl", "Br", "I"}
    _VALENCE = {"C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "Si": 4, "B": 3}
    _nb_full = {}
    for a, b in bonds:
        _nb_full.setdefault(a, []).append(b)
        _nb_full.setdefault(b, []).append(a)

    def _overvalent(atom):
        val = _VALENCE.get(sym[atom])
        if val is None:
            return False
        # metal-coordinated ring N of a pyridine-type ring keeps its
        # C=N (class A): donating the lone pair leaves N+ with the
        # double bond intact.  Detect: N with a metal neighbour and
        # exactly two carbon neighbours (ring N, no substituent).
        if sym[atom] == "N":
            ns = _nb_full.get(atom, ())
            if any(sym[k] in _MET for k in ns):
                cs = [k for k in ns if sym[k] == "C"]
                hs = [k for k in ns if sym[k] == "H"]
                if len(cs) == 2 and not hs:
                    return False
        deg_heavy = sum(1 for k in _nb_full.get(atom, ())
                        if sym[k] != "H")
        if deg_heavy < val:
            return False
        return True

    upgrades = {pair: o for pair, o in upgrades.items()
                if not (_overvalent(pair[0]) or _overvalent(pair[1]))}
    taken = set()
    matched = {}
    for pair in sorted(upgrades, key=_len):
        i, j = pair
        if i in taken or j in taken:
            continue
        matched[pair] = upgrades[pair]
        taken.add(i)
        taken.add(j)
    return matched or _infer_bond_orders(sym, bonds)


def _infer_bond_orders(sym, bonds, rotated_pos=None):
    """Assign bond orders from valence deficits (heuristic, no charges).

    Each heavy atom gets its typical valence (C 4, N 3, O 2, S 2, ...);
    every single bond consumes one bond from each partner.  Deficits are
    filled greedily: two atoms each missing >=2 bonds form a double bond,
    missing >=3 a triple bond.  Metals and H never take multiple bonds.
    Returns {(i, j): 2|3} for the bonds upgraded.
    """
    import collections
    _VALENCE = {"C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "Si": 4, "B": 3}
    order = collections.defaultdict(int)   # per-atom bonds consumed
    upgrades = {}
    for i, j in bonds:
        order[i] += 1
        order[j] += 1
    # Metal bonds are dative and do not consume organic valence: a
    # pyridine N coordinated to Fe still has its valence 3 available,
    # so the C=N of the ring is recognised again (same for fluorenyl).
    _METAL = set(_XYZ_MASS) - {"H", "B", "C", "N", "O", "F", "Si",
                               "P", "S", "Cl", "Br", "I"}
    need = {}
    for i, s in enumerate(sym):
        v = _VALENCE.get(s)
        if v is None or s == "H" or s in _METAL:
            continue
        # H bonds DO consume valence (a CH2 carbon has none left);
        # only metal bonds are dative and skipped.
        used = sum(1 for i2, j2 in bonds
                   if (i2 == i or j2 == i) and sym[i2] not in _METAL
                   and sym[j2] not in _METAL)
        need[i] = max(0, v - used)
    # Aromatic/ring alternating bonds leave BOTH partners one short
    # (pyridine C=N, benzene C=C, fluorenyl), so single-deficit pairs
    # upgrade to double as well; higher deficits stay greedy.
    def bond_length(i, j):
        if rotated_pos is None:
            return 0.0
        import numpy
        di = rotated_pos[i] - rotated_pos[j]
        return float(numpy.sqrt(di @ di))

    # At most ONE double bond per atom: in a pyridine ring every C is
    # one bond short, and the greedy all-pairs pass drew double lines
    # on adjacent ring bonds (wrong alternation, user report 406).
    # One double bond per atom, preferring larger deficits and shorter
    # bonds (true C=N/C=C are shorter than their single neighbours),
    # reproduces the alternating pattern.
    # Only TWO situations deserve a drawn double bond:
    # 1. both partners sit in a common ring (aromatic alternation -
    #    pyridine C=N/C=C, benzene, fluorenyl).  Ring bonds are found
    #    with a simple cycle search over the bond graph.
    # 2. a short C=O carbonyl (bond < 1.30 A; a hydroxyl C-OH is
    #    ~1.36 A and stays single - its O is labeled OH instead).
    # Everything else (bridging single bonds, long C-N links) must NOT
    # be doubled even when both ends are one bond short - the greedy
    # all-pairs rule drew double lines everywhere (user report 406).
    METAL = set(_XYZ_MASS) - {"H", "B", "C", "N", "O", "F", "Si",
                              "P", "S", "Cl", "Br", "I"}
    # Chelate rings through the metal are coordination, not aromatic -
    # their bonds must never be doubled.
    rings = [r for r in _find_rings(sym, bonds)
             if not any(sym[a] in METAL for a in r)]
    in_same_ring = {}
    for ring in rings:
        rset = set(ring)
        for a in ring:
            for b in ring:
                if a < b:
                    in_same_ring[(a, b)] = True
    atom_upgraded = set()
    candidates = []
    for i, j in bonds:
        if i not in need or j not in need:
            continue
        if i in atom_upgraded or j in atom_upgraded:
            continue
        deficit = min(need[i], need[j])
        length = bond_length(i, j)
        carbonyl = ({sym[i], sym[j]} == {"C", "O"} and length < 1.30)
        if (i, j) not in in_same_ring and not carbonyl:
            continue
        candidates.append((deficit, -length, i, j))
    candidates.sort(reverse=True)
    for deficit, _neglen, i, j in candidates:
        if deficit >= 3:
            upgrades[(i, j)] = 3
            need[i] -= 3; need[j] -= 3
        elif deficit >= 1:
            upgrades[(i, j)] = 2
            need[i] -= 1; need[j] -= 1
        atom_upgraded.add(i)
        atom_upgraded.add(j)
    return upgrades


def _find_rings(sym, bonds):
    """Small rings (<= 8 atoms) of the bond graph, as index tuples.

    Simple DFS cycle detection over heavy-atom bonds; returns each
    shortest cycle once.  Good enough to mark aromatic systems.
    """
    import collections
    adj = collections.defaultdict(set)
    for i, j in bonds:
        if sym[i] != "H" and sym[j] != "H":
            adj[i].add(j)
            adj[j].add(i)
    rings = []
    seen = set()
    for start in adj:
        # BFS shortest cycle through start via each neighbour pair.
        for n1 in adj[start]:
            for n2 in adj[start]:
                if n1 >= n2:
                    continue
                # shortest path n1 -> n2 not through start
                prev = {start: None}
                queue = collections.deque([n1])
                dist = {n1: 0}
                parent = {}
                while queue:
                    cur = queue.popleft()
                    if cur == n2:
                        break
                    for nb in adj[cur]:
                        if nb == start:
                            continue
                        if nb not in dist:
                            dist[nb] = dist[cur] + 1
                            parent[nb] = cur
                            queue.append(nb)
                if n2 in dist and dist[n2] <= 7:
                    path = [n2]
                    while path[-1] != n1:
                        path.append(parent[path[-1]])
                    cycle = tuple(sorted(path + [start]))
                    if cycle not in seen and len(cycle) <= 8:
                        seen.add(cycle)
                        rings.append(cycle)
    return rings


def _render_skeleton(sym, rotated, bonds, size,
                     bond_orders=None, fmt="png", colored=False):
    """Wireframe of the rotated structure in RDKit-2D style (PNG bytes).

    Same geometry as the ball-and-stick column (positions kept, no
    snapping to ideal bond lengths): bonds as black lines - H-bonds thin,
    heavy-atom bonds thick - atom symbols drawn inside gaps, H symbols
    suppressed so it reads like a skeletal formula.
    """
    import numpy
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from io import BytesIO

    radii = numpy.array([_XYZ_RADII.get(s, _XYZ_RADIUS_DEFAULT) for s in sym])
    # Headroom framing: cover outermost atoms PLUS their display radii,
    # else carbonyl O (428) or balls (434) clip at the frame edge.
    disp = 0.55 * radii
    cx, cy = (0.5 * (rotated[:, 0].max() + rotated[:, 0].min()),
              0.5 * (rotated[:, 1].max() + rotated[:, 1].min()))
    rmax = float(disp.max())
    scale = size / (2.0 * (max(abs(rotated[:, 0] - cx).max() + rmax,
                               abs(rotated[:, 1] - cy).max() + rmax)))
    rotated = rotated - numpy.array([cx, cy] + [0.0]
                            * (rotated.shape[1] - 2))
    fig = plt.figure(figsize=(size / 100, size / 100), dpi=100)
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    half = 0.5 * size / scale
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_aspect("equal")

    is_h = numpy.array([s == "H" for s in sym])
    heavy_bonds = [(i, j) for i, j in bonds if not (is_h[i] or is_h[j])]
    h_bonds = [(i, j) for i, j in bonds if (is_h[i] or is_h[j])]

    # A double/triple bond REPLACES its base line: the parallel lines
    # must sit visibly beside each other.  With lw 1.8 a 0.06 A offset
    # hides entirely under the base line (bug: upgrades were computed
    # and drawn but invisible).
    _dbl = {tuple(sorted(p)) for p in bond_orders} if bond_orders else set()
    for i, j in heavy_bonds:
        if tuple(sorted((i, j))) in _dbl:
            continue
        ax.plot([rotated[i][0], rotated[j][0]],
                [rotated[i][1], rotated[j][1]],
                color="black", lw=1.8, solid_capstyle="round", zorder=1)
    for i, j in h_bonds:
        ax.plot([rotated[i][0], rotated[j][0]],
                [rotated[i][1], rotated[j][1]],
                color="black", lw=0.8, solid_capstyle="round", zorder=1)

    # Skeletal convention: C is implicit (never labelled), only
    # heteroatoms (+ metals) get a white gap and a symbol.  H labelled
    # only when show_h is on.
    _implicit = {"C"}
    # A terminal, single-bonded O carries its proton in these
    # complexes (acid OH / ester O-alkyl chain end): label it "OH"
    # instead of a bare "O" (user request, 406).
    _degree = [0] * len(sym)
    _dbl_here = {tuple(sorted(p)) for p in (bond_orders or {})}
    for i, j in bonds:
        _degree[i] += 1
        _degree[j] += 1
    zcol = rotated[:, 2] if rotated.shape[1] > 2 else numpy.zeros(len(sym))
    for idx in numpy.argsort(zcol):
        s = sym[idx]
        if s == "H":
            continue  # heteroatom protons were kept by the caller; C-H are gone
        if s == "O" and _degree[idx] <= 1:
            # A carbonyl O (double bond to its only neighbour) must
            # NOT be labeled OH - only a single-bonded hydroxyl O
            # carries the proton.
            on_double = any(
                idx in (a, b) and tuple(sorted((a, b))) in _dbl_here
                for a, b in bonds)
            # The OH/HO label requires a REAL proton on the O.  A
            # deprotonated terminal O (carboxylate O-, nitro O, phenolate
            # drawn as terminal) has no H; forcing "OH" onto it invented
            # hydrogens (user report classes E/F: "HO-N=O" on nitro,
            # "HO-CH" in the oxalyl region).
            _has_h = any(idx in (a, b) and sym[b if a == idx else a] == "H"
                         for a, b in bonds)
            if not on_double and _has_h:
                # Direction along the O-C bond, not global x: the
                # partner sits at the other end of the O's only bond.
                # The label reads away from the partner: when the C
                # sits to the right of the O the text goes left of the
                # O, so it must read "HO" (and vice versa).
                _pc = next((j if i == idx else i) for i, j in bonds
                           if idx in (i, j))
                s = "HO" if rotated[idx][0] < rotated[_pc][0] else "OH"
        if s in _implicit:
            continue
        x, y = rotated[idx][0], rotated[idx][1]
        r_data = radii[idx] * scale * 0.55 / size * (2 * half)
        ax.add_patch(mpatches.Circle(
            (x, y), r_data, facecolor="white", edgecolor="none", zorder=2))
        ax.text(x, y, s, ha="center", va="center", fontsize=11,
                color=(colored and _XYZ_COLORS.get(s, "black")
                       or "black"), zorder=3)
    # Double / triple bonds: parallel offset lines per order.
    if bond_orders:
        import math
        for (i, j), order in bond_orders.items():
            if order < 2 or order > 3:
                continue
            x1, y1 = rotated[i][0], rotated[i][1]
            x2, y2 = rotated[j][0], rotated[j][1]
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy) or 1.0
            nx, ny = -dy / length, dx / length
            gap = 0.10 * length if order == 2 else 0.22 * length
            style = dict(color="black", lw=1.2,
                         solid_capstyle="round", zorder=1)
            if order == 2:
                offs = (gap / 2, -gap / 2)
            else:
                offs = (0.0, gap, -gap)  # middle + two outer
            for off in offs:
                ax.plot([x1 + nx * off, x2 + nx * off],
                        [y1 + ny * off, y2 + ny * off], **style)
    buf = BytesIO()
    fig.savefig(buf, format=fmt, dpi=100)
    plt.close(fig)
    return buf.getvalue()


def optimize_view_direction(sym, pts3d, bonds, n_samples: int = 500):
    """Fibonacci-sphere search for the least-overlapping camera direction.

    Samples ~500 directions (Fibonacci sphere) plus the three PCA axes;
    scores each by projected atom overlap and bond crossings and breaks
    ties by the larger projected area.  Returns the view unit vector.
    """
    import numpy
    pts3d = numpy.asarray(pts3d, dtype=float)
    radii = numpy.array([_XYZ_RADII.get(s, _XYZ_RADIUS_DEFAULT)
                         for s in sym])
    # Fibonacci sphere directions.
    i = numpy.arange(n_samples, dtype=float) + 0.5
    phi = numpy.pi * (3.0 - numpy.sqrt(5.0))
    z = 1.0 - 2.0 * i / n_samples
    r = numpy.sqrt(numpy.maximum(0.0, 1.0 - z * z))
    dirs = numpy.stack([r * numpy.cos(phi * i), r * numpy.sin(phi * i), z],
                       axis=1)
    # PCA axes as candidates: the inertia-view (smallest extent on z) and
    # the two neighbours.
    pts3d_c = pts3d - pts3d.mean(axis=0)
    cov = (pts3d_c * 1.0).T @ pts3d_c  # unweighted for extent ordering
    try:
        _w, vecs = numpy.linalg.eigh(cov)
    except numpy.linalg.LinAlgError:
        vecs = numpy.eye(3)
    for k in range(3):
        dirs = numpy.vstack([dirs, vecs[:, k]])

    # Precompute pairwise projected data once per direction; scoring is
    # vectorized.  Bond-crossing penalty: count segment intersections.
    def total_score(view):
        helper = numpy.array([1.0, 0.0, 0.0])
        if abs(view[0]) > 0.9:
            helper = numpy.array([0.0, 1.0, 0.0])
        x = numpy.cross(view, helper)
        x /= numpy.linalg.norm(x)
        y = numpy.cross(view, x)
        px, py = pts3d_c @ x, pts3d_c @ y
        dx = px[:, None] - px[None, :]
        dy = py[:, None] - py[None, :]
        d = numpy.sqrt(dx * dx + dy * dy)
        iu = numpy.triu_indices(len(sym), 1)
        sums = 0.8 * (radii[:, None] + radii[None, :])
        overlap = numpy.maximum(0.0, sums[iu] - d[iu])
        s = float(overlap.sum())
        # Crossing bonds in projection: segment intersection test.
        cross = 0
        for a, b in bonds:
            for c, d2 in bonds:
                if len({a, b} & {c, d2}) == 1:
                    # 2D segment intersection (cross products).
                    p0, p1 = (px[a], py[a]), (px[b], py[b])
                    q0, q1 = (px[c], py[c]), (px[d2], py[d2])
                    d1 = (p1[0] - p0[0], p1[1] - p0[1])
                    d2v = (q1[0] - q0[0], q1[1] - q0[1])
                    denom = d1[0] * d2v[1] - d1[1] * d2v[0]
                    if abs(denom) < 1e-12:
                        continue
                    t = ((q0[0] - p0[0]) * d2v[1]
                         - (q0[1] - p0[1]) * d2v[0]) / denom
                    u = ((q0[0] - p0[0]) * d1[1]
                         - (q0[1] - p0[1]) * d1[0]) / denom
                    if 0.02 < t < 0.98 and 0.02 < u < 0.98:
                        cross += 1
        area = float(numpy.ptp(px) * numpy.ptp(py))
        return s + 2.0 * cross, area

    best = None
    for k in range(dirs.shape[0]):
        sc, area = total_score(dirs[k])
        if best is None or sc < best[0] or (sc == best[0] and area > best[1]):
            best = (sc, area, dirs[k])
    return best[2]


def optimized_3d_png_bytes(calc_dir: str, size: int = 320, show_h=False,
                           fmt="png"):
    """Structure view with the camera direction optimized for least overlap,

    drawn in the style of the dashboard 3Dmol viewer: small van-der-Waals
    spheres (scale 0.28) on thick bonds (0.11 A radius), Jmol palette.
    ~500 Fibonacci-sphere directions (+ the three PCA axes) are scored for
    projected atom overlap and bond crossings, the best faces the camera.
    Returns None on any failure.
    """
    atoms = _parse_xyz(calc_dir)
    if not atoms:
        return None
    try:
        import numpy
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        return None
    sym = [a[0] for a in atoms]
    pts = numpy.array([[a[1], a[2], a[3]] for a in atoms])
    cov = [_XYZ_COVALENT.get(s, _XYZ_COVALENT_DEFAULT) for s in sym]
    bonds = _detect_bonds(sym, cov, pts)
    try:
        view = optimize_view_direction(sym, pts, bonds)
    except Exception:
        return None
    # Rotate onto that view: z-axis along the view vector.
    z = numpy.asarray(view, dtype=float)
    z /= numpy.linalg.norm(z)
    helper = numpy.array([1.0, 0.0, 0.0])
    if abs(z[0]) > 0.9:
        helper = numpy.array([0.0, 1.0, 0.0])
    x = numpy.cross(z, helper)
    x /= numpy.linalg.norm(x)
    y = numpy.cross(z, x)
    pts_c = pts - pts.mean(axis=0)
    rotated = numpy.stack([pts_c @ x, pts_c @ y, pts_c @ z], axis=1)
    # Dashboard viewer style: small spheres on thick sticks
    # (molecule_viewer.py DEFAULT_3DMOL_ATOM_SCALE=0.28,
    # DEFAULT_3DMOL_BOND_RADIUS=0.11 vs. the classic 0.8/1.4 here).
    if not show_h:
        is_h = numpy.array([s == "H" for s in sym])
        keep = [i for i in range(len(sym)) if not is_h[i]]
        keepset = set(keep)
        remap = {old: new for new, old in enumerate(keep)}
        sym = [sym[i] for i in keep]
        rotated = rotated[keep]
        bonds = [(remap[i], remap[j]) for i, j in bonds
                 if i in keepset and j in keepset]
    return _render_ball_stick(sym, rotated, bonds, size,
                              atom_scale=0.55, bond_lw=3.2)


def rss_metrics(calc_dir: str) -> Dict[str, Any]:
    """RSS energetics: E(2.0 A) vs E(3.5 A) plus shape around the bound point.

    The unbound reference is 3.5 A, not 4.0 A: at 4.0 A the curve is already
    near a plain minimum and the sign of the difference mislabels endergonic
    coordination.  3.5 A is still on the descending branch, so
    ``delta = E(2.0) - E(3.5)`` answers "does binding pay off" cleanly.

    When the delta is endergonic the bound point can still be a true local
    minimum or sit on a shoulder/saddle.  The shape is judged from the
    neighbourhood (1.7 ... 2.3 A): if E rises on both sides of 2.0 A, it is
    a minimum; if it falls monotonically towards 2.0 A from above while the
    short side keeps falling, it reads as "no clear minimum (saddle/shoulder)".
    """
    points = parse_rss_dat(calc_dir)
    result: Dict[str, Any] = {}
    if not points:
        return result
    # Keep the first sweep only: distances decreasing from 4.0 A.
    sweep: List[tuple] = []
    prev = None
    for dist, energy in points:
        if prev is not None and dist > prev + 1e-9:
            break
        sweep.append((dist, energy))
        prev = dist
    if not sweep:
        sweep = points

    by_dist = {round(d, 4): e for d, e in sweep}
    ref = by_dist.get(round(RSS_REF_DISTANCE, 4))
    bound = by_dist.get(round(RSS_BOUND_DISTANCE, 4))
    if ref is None or bound is None:
        return result

    delta_ha = bound - ref
    result["rss_delta_kj"] = round(delta_ha * HARTREE_KJ_MOL, 2)
    result["rss_direction"] = ("exergonic" if delta_ha < 0
                          else "endergonic")
    result["rss_points"] = len(sweep)

    # Shape judgement around the bound point (only meaningful endergonic;
    # reported either way so the table stays complete).
    neighbours = {}
    for d in RSS_SHAPE_DISTANCES:
        e = by_dist.get(round(d, 4))
        if e is not None:
            neighbours[d] = e
    # Shape judgement per user specification (2026-10-30):
    # exergonic -> "global minimum"; endergonic -> E(2.1 A):
    # higher than E(2.0 A) -> "local minimum", else no minimum.
    # The old left/right neighbourhood criteria are retired.
    shape = None
    if result.get("rss_direction") == "exergonic":
        shape = "global minimum"
    elif result.get("rss_direction") == "endergonic":
        e_20 = by_dist.get(round(RSS_BOUND_DISTANCE, 4))
        e_21 = by_dist.get(round(RSS_BOUND_DISTANCE + 0.1, 4))
        if e_20 is not None and e_21 is not None:
            shape = ("local minimum" if e_21 > e_20
                     else "no minimum")
    if shape:
        result["rss_shape"] = shape
    # Legacy columns kept: global minimum of the sweep (wherever it sits).
    window = [(d, e) for d, e in sweep if d < round(RSS_REF_DISTANCE, 4) - 0.05]
    if window:
        min_d, min_e = min(window, key=lambda p: p[1])
        result["rss_min_distance"] = round(min_d, 3)
        result["rss_min_delta_kj"] = round(
            (min_e - ref) * HARTREE_KJ_MOL, 2)
    tight = [
        (d, e) for d, e in sweep if RSS_MIN_WINDOW[0] <= d <= RSS_MIN_WINDOW[1]
    ]
    if tight:
        t_d, t_e = min(tight, key=lambda p: p[1])
        result["rss_tight_min_distance"] = round(t_d, 3)
        result["rss_tight_delta_kj"] = round(
            (t_e - ref) * HARTREE_KJ_MOL, 2)
    return result


def extract(calc_dir: str, metric_ids: List[str]) -> Dict[str, Any]:
    """Extract the requested metrics from one calculation folder.

    Returns ``{metric_id: value}``; absent quantities are omitted.
    """
    data = _load_delfin_data(calc_dir)
    wanted = set(metric_ids)
    out: Dict[str, Any] = {}
    if data is None:
        # No DELFIN_Data.json: still offer RSS + PNGs which need no JSON.
        if wanted & {"rss_min_distance", "rss_delta_kj", "rss_direction",
                     "rss_shape", "rss_min_delta_kj",
                     "rss_tight_min_distance", "rss_tight_delta_kj"}:
            out.update(rss_metrics(calc_dir))
        if "rss_pngs" in wanted:
            out["rss_pngs"] = list_rss_pngs(calc_dir)
        return out

    # Extractor that scans the ORCA inputs of a calculation for
    # ``%scf BrokenSym ... end`` lines, keyed by state.  Defined here so
    # extract() and the tests share one source of truth.
    def broken_sym_map() -> Dict[str, str]:
        found: Dict[str, str] = {}
        try:
            inps = sorted(glob.glob(os.path.join(calc_dir, "*.inp")))
        except OSError:
            return found
        for inp in inps:
            stem = os.path.splitext(os.path.basename(inp))[0]
            try:
                with open(inp, "r", encoding="utf-8", errors="replace") as fh:
                    for line in fh:
                        if "brokensym" in line.lower():
                            found[stem] = line.strip()
                            break
            except OSError:
                continue
        return found

    meta = data.get("metadata") or {}
    rp = _redox(data) or {}

    def want(mid: str) -> bool:
        return mid in wanted

    if want("method"):
        out["method"] = "/".join(
            str(meta.get(k, "?")) for k in
            ("functional", "basis_set", "solvent", "charge")
        )
    if want("redox_E_ox"):
        # Two decimals only: redox potentials carry no more precision than
        # the method supports.
        v = (rp.get("final") or {}).get("E_ox")
        out["redox_E_ox"] = round(v, 2) if isinstance(v, float) else v
    if want("redox_E_red"):
        v = (rp.get("final") or {}).get("E_red")
        out["redox_E_red"] = round(v, 2) if isinstance(v, float) else v
    if want("redox_E_red_2"):
        v = (rp.get("final") or {}).get("E_red_2")
        out["redox_E_red_2"] = round(v, 2) if isinstance(v, float) else v
    if want("redox_E_ref"):
        out["redox_E_ref"] = rp.get("E_ref")
    if want("redox_method"):
        flags = rp.get("method_flags") or {}
        active = [k.lstrip("use_m") for k, v in flags.items() if v]
        out["redox_method"] = "+".join(active) if active else None
    if want("mult_ground"):
        out["mult_ground"] = _format_mult(
            _state_opt(data, "ground_state_S0"),
            broken_sym_map().get("initial"))
    if want("mult_oxidized"):
        out["mult_oxidized"] = _format_mult(
            _state_opt(data, "oxidized_state"),
            broken_sym_map().get("ox_step_1"))
    if want("mult_reduced"):
        out["mult_reduced"] = _format_mult(
            _state_opt(data, "reduced_state"),
            broken_sym_map().get("red_step_1"))
    if want("mult_reduced2"):
        red2 = (data.get("reduced_states") or {}).get("red_step_2")
        out["mult_reduced2"] = _format_mult(
            red2.get("optimization") if isinstance(red2, dict) else None,
            broken_sym_map().get("red_step_2"))
    if want("gibbs_ground"):
        thermo = (data.get("ground_state_S0") or {}).get("thermochemistry") or {}
        out["gibbs_ground"] = thermo.get("total_gibbs_free_energy_hartree")
    if want("gap_ground"):
        orbitals = (data.get("ground_state_S0") or {}).get("orbitals") or {}
        out["gap_ground"] = orbitals.get("gap_eV")
    if want("converged_flags"):
        flags = []
        states = [("S0", data.get("ground_state_S0")),
                  ("ox", data.get("oxidized_state")),
                  ("red", data.get("reduced_state"))]
        for name, payload in states:
            opt = (payload or {}).get("optimization") if payload else None
            conv = opt.get("converged") if opt else None
            flags.append((name, conv))
        for prefix, mapping in (
            ("ox", data.get("oxidized_states") or {}),
            ("red", data.get("reduced_states") or {}),
        ):
            for step, payload in sorted(mapping.items()):
                opt = (payload or {}).get("optimization") or {}
                flags.append((f"{prefix}/{step}", opt.get("converged")))
        names_true = [n for n, c in flags if c is True]
        names_false = [n for n, c in flags if c is False]
        names_unknown = [n for n, c in flags if c is not True and c is not False]
        if not names_false and not names_unknown and names_true:
            out["converged_flags"] = "true"
        else:
            parts = ["false"]
            if names_false:
                parts.append("(" + ", ".join(names_false) + " false)")
            if names_unknown:
                parts.append("(" + ", ".join(names_unknown) + " unbekannt)")
            out["converged_flags"] = " ".join(parts)
    if wanted & {"rss_min_distance", "rss_delta_kj", "rss_direction",
                 "rss_shape", "rss_min_delta_kj",
                 "rss_tight_min_distance", "rss_tight_delta_kj"}:
        out.update(rss_metrics(calc_dir))
    if want("rss_pngs"):
        out["rss_pngs"] = list_rss_pngs(calc_dir)
    return out


def list_rss_pngs(calc_dir: str) -> List[Dict[str, str]]:
    """All PNGs under the relaxed surface scan dir, relative to calc_dir."""
    pattern = os.path.join(calc_dir, RSS_DIR, "*.png")
    rels = sorted(
        os.path.relpath(p, calc_dir) for p in glob.glob(pattern)
    )
    return [{"rel": rel, "abs": os.path.join(calc_dir, rel)} for rel in rels]


def format_value(value: Any) -> str:
    """One table cell: value or a dash for "not computed"."""
    if value is None:
        return "–"
    if isinstance(value, float):
        if abs(value) >= 1e6:
            return f"{value:.4g}"
        return f"{value:.4f}".rstrip("0").rstrip(".")
    return str(value)


def sort_table(table: Dict[str, Any], column_id: str,
               descending: bool = False) -> Dict[str, Any]:
    """Return a copy of *table* sorted by one column.

    Numeric cells sort numerically (None always last), text cells
    alphabetically; the label column sorts by label.
    """
    if column_id == "label":
        key = lambda r: r["label"]  # noqa: E731
    else:
        key = lambda r: r["values"].get(column_id)  # noqa: E731
    rows = [r for r in table["rows"] if r["values"].get(column_id) is not None]
    rest = [r for r in table["rows"] if r["values"].get(column_id) is None]
    numeric = all(
        isinstance(r["values"].get(column_id), (int, float))
        for r in rows
    ) and column_id != "label"
    if numeric:
        rows.sort(key=key, reverse=descending)
    else:
        rows.sort(key=lambda r: str(key(r)), reverse=descending)
    return {**table, "rows": rows + rest}


def _json_path_get(data: Any, path: str):
    """Follow a dotted JSON path like ``computed.final.energy``; list
    indices as digits, ``[*]`` applies the following hop to every list
    element (collecting a list of values), and ``?(@.key=="value")``
    filters list elements by string equality. Returns None when any
    hop is missing."""
    def walk(cur: Any, hops: List[str]) -> Any:
        if not hops:
            return cur
        hop = hops[0]
        # Re-join hops that were split inside a ?(...) filter: the '.'
        # after '?(@' belongs to the filter (@.key), not to the path.
        # Keep consuming hops until the filter's ')]' ends.
        if '?(@' in hop and not hop.rstrip().endswith(')]'):
            while hops[1:] and not hop.endswith(')]'):
                hop += '.' + hops.pop(1)
        filt = None
        if "?(@" in hop and '=="' in hop:
            base, _, tail = hop.partition('?(@.')
            base = base.rstrip('[')
            tail = tail.rstrip(']').rstrip(')')
            conds = [c.strip() for c in tail.split('&&')]
            parsed = []
            for cond in conds:
                key, _, val = cond.partition('=="')
                val = val.split('"')[0] if val.startswith('"') else val
                parsed.append((key.strip('@. '), val.rstrip('"')))
            filt = parsed
            hop = base
        star = hop.endswith("[*]")
        if star:
            hop = hop[:-3]
        if isinstance(cur, dict):
            cur = cur.get(hop)
        elif isinstance(cur, list) and hop.isdigit():
            try:
                cur = cur[int(hop)]
            except IndexError:
                return None
        else:
            return None
        if cur is None:
            return None
        if filt is not None and isinstance(cur, list):
            cur = [el for el in cur if isinstance(el, dict) and all(
                str(el.get(key)) == val for key, val in filt)]
            if not cur:
                return None
            # The filter selected elements; a following field hop applies
            # to each selected element (or the single one, directly).
            if len(cur) == 1:
                cur = cur[0]
        if star and isinstance(cur, list):
            rest = hops[1:]
            if not rest:
                return cur
            results = []
            for el in cur:
                got = walk(el, rest)
                if got is not None:
                    results.append(got)
            return results
        return walk(cur, hops[1:])
    return walk(data, path.split("."))


def extract_custom(calc_dir: str, columns: List[Dict[str, Any]]) -> Dict[str, Any]:
    """User-defined columns, semantics from the calculations-browser
    table panel: {id, label, file, type: text|regex|jsonpath,
    pattern, occurrence: first|last|all}.

    text: the raw matched line (file must exist); regex: first capture
    group or the whole match; jsonpath: a dotted path into
    DELFIN_Data.json. occurrence selects which match is returned
    (all joins matches with ``; ``). Missing values are omitted.
    """
    out: Dict[str, Any] = {}
    for col in columns:
        cid = str(col.get("id") or "").strip()
        if not cid:
            continue
        ftype = (col.get("type") or "text").lower()
        pattern = col.get("pattern")
        occurrence = (col.get("occurrence") or "first").lower()
        value = None
        try:
            if ftype == "jsonpath":
                data = _load_delfin_data(calc_dir)
                if data is not None:
                    value = _json_path_get(data, str(pattern or ""))
            else:
                target = os.path.join(calc_dir, str(col.get("file") or ""))
                if os.path.isfile(target):
                    with open(target, "r", encoding="utf-8",
                              errors="replace") as fh:
                        content = fh.read()
                    if ftype == "regex":
                        matches = re.findall(str(pattern or ""), content)
                        if matches:
                            if occurrence == "all":
                                value = "; ".join(
                                    m if isinstance(m, str)
                                    else (m[0] if m and m[0] else "")
                                    for m in matches)
                            elif occurrence == "last":
                                m = matches[-1]
                                value = (m if isinstance(m, str)
                                         else (m[0] if m and m[0] else ""))
                            else:
                                m = matches[0]
                                value = (m if isinstance(m, str)
                                         else (m[0] if m and m[0] else ""))
                    else:  # text: every line containing the pattern
                        hits = [ln.strip() for ln in content.splitlines()
                                if str(pattern or "") in ln]
                        if hits:
                            if occurrence == "all":
                                value = "; ".join(hits)
                            elif occurrence == "last":
                                value = hits[-1]
                            else:
                                value = hits[0]
        except (OSError, re.error):
            value = None
        if value not in (None, ""):
            out[cid] = value
    return out


def build_table(calc_dirs: List[str], metric_ids: List[str],
                custom_columns: Optional[List[Dict[str, Any]]] = None
                ) -> Dict[str, Any]:
    """Rows for the dashboard table: one per calculation, columns per metric.

    Metric ids are de-duplicated, unknown ids are dropped, and the label of
    each column comes from :data:`METRICS`. ``custom_columns`` appends
    user-defined text/regex/jsonpath columns (see :func:`extract_custom`).
    """
    seen: List[str] = []
    known_ids = {m["id"] for m in METRICS}
    for mid in metric_ids:
        if mid in known_ids and mid not in seen:
            seen.append(mid)
    custom = list(custom_columns or [])
    rows: List[Dict[str, Any]] = []
    for calc_dir in calc_dirs:
        values = extract(calc_dir, seen)
        values.update(extract_custom(calc_dir, custom))
        rows.append({
            "path": calc_dir,
            "label": os.path.basename(os.path.normpath(calc_dir)),
            "values": values,
        })
    columns = [{"id": m["id"], "label": m["label"]} for m in METRICS
               if m["id"] in seen]
    columns.extend({"id": c["id"], "label": str(c.get("label") or c["id"])}
                   for c in custom if str(c.get("id") or "").strip())
    return {"columns": columns, "rows": rows}
