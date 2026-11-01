"""Conformer shape classification engine.

Continuous geometric shape descriptors for conformer ensembles
(multi-frame XYZ), for the whole molecule and for explicit atom
selections, plus discrete chemical classifications with explicit
status. Continuous shape descriptors and chemical categories are two
separate layers: shape_descriptors() never invents classes, and
classification functions return DescriptorStatus with documented
reasons whenever the input does not determine the answer.

Definitions (all exact; thresholds are caller-supplied, none validated):

  Weighted covariance S = sum_i w_i (r_i - c)(r_i - c)^T / W,
  centroid c = sum_i w_i r_i / W, W = sum_i w_i.
  Principal moments i1 <= i2 <= i3 = eigenvalues of S.

  Rg     = sqrt(t1 + t2 + t3)                       [Angstrom]
  dmax   = max pairwise atom-atom distance          [Angstrom]
  asphericity = ((t3-t2)^2 + (t3-t1)^2 + (t2-t1)^2) / (2 (t1+t2+t3)^2)
           (DOI 10.1021/ma00148a028; gyration moments t_i);
           0 = spherical, 0.5 = disc-like lamella, 1 = rod.
           NOT a planarity measure.
  eccentricity = sqrt(I3^2 - I1^2)/I3 [0,1]; 0 sphere, 1 rod
           (physical inertia moments I_i, atomic masses).
  npr1 = I1/I3 in [0,1] (1 = spherical top), npr2 = I2/I3 in [0,1]
           (RDKit NPR1/NPR2, Sauer & Schwarz, DOI 10.1021/ci025599w;
           useAtomicMasses=True). Rod: npr1 -> 0; disc: npr1 ~ 0.5.
  planarity = (t2 - t1) / (t3 - t1), [0,1]; 1 = planar lamella
           (inertia basis: I1 = W(t1+t2), I2 = W(t1+t3), I3 = W(t2+t3)
           are the sorted PHYSICAL principal moments (= RDKit PMI1/2/3)
           and the identity (t2-t1)/(t3-t1) == (I3-I2)/(I3-I1) makes
           the two forms equal, so the gyration eigenvalues t are used
           for numerical stability).
           High npr2 with low npr1 distinguishes rod from disc.
  acylindricity = (t2 - t1) / (t1 + t2 + t3); 0 = cylinder-like
  kappa2 (relative shape anisotropy, DOI 10.1021/ma00148a028) = 1.5 * sum_k
           (tk - s/3)^2 / s^2, s = t1+t2+t3, [0,1]; 0 = isotropic
           sphere-like, 1 = rod-like. NOTE: identical to Baumgaertner
           asphericity (both = 1.5*sum(t^2)/s^2 - 0.5), NOT its double.
  compacity = Rg / dmax; ~0.5 fully extended rod, smaller = compact
  foldedness = 1 - Rg / Rg_uniform-ball(dmax) with
           Rg_ball = (dmax/2) * sqrt(3/5); 0 = ball-like filling its
           own diameter, negative = sparse/extended.

  Weightings: "mass" (atomic masses), "uniform" (w=1; translation and
  rotation invariant in the exact same way), "heavy" (w=1 on non-H, 0
  on H; errors out when the selection is all-H). Results from different
  weightings are NOT comparable; each is labelled in the output.

Chemical classification works on a connectivity-only RDKit molecule
(no bond orders from XYZ): E/Z uses the dihedral of the two highest-CIP
heavy substituents across an acyclic C=C-like pair and is reported as
"geometric" with status not_assigned unless both substituents are
unambiguous; s-cis/s-trans uses the dihedral across a central
sp2-sp2 single bond. Unknown cases get status "not_assigned" + reason.

CLI:
  python -m delfin.analysis_tools.conformer_shape_classification \
      ensemble.xyz [--json OUT] [--weights mass|uniform|heavy]
      [--atoms "0-5,9"] [--ranking key:asc key:desc] [--cache DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from delfin.tools.adapters.crest_ensemble import parse_ensemble

__all__ = [
    "Ensemble", "DescriptorStatus", "MASS", "DESCRIPTOR_DOC",
    "load_ensemble", "shape_descriptors",
    "double_bond_geometry", "diene_conformation", "ez_cip_assignment",
    "rotamer_states", "custom_torsion",
    "fragment_distance", "fragment_orientation",
    "rank_conformers", "filter_conformers",
    "analyze_ensemble", "cache_key", "load_or_compute", "main",
]

MASS = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990,
    "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "Ar": 39.95, "K": 39.098, "Ca": 40.078,
}

DESCRIPTOR_DOC = {
    "rg": {"unit": "A", "definition": "sqrt(i1+i2+i3) of weighted covariance"},
    "dmax": {"unit": "A", "definition": "max pairwise atom distance in selection"},
    "asphericity": {"unit": "1", "definition": "Baumgaertner over gyration moments: 0.5*sum (ti-tj)^2/(sum t)^2, [0,1]; 0 sphere, 1 rod; not planarity"},
    "eccentricity": {"unit": "1", "definition": "sqrt(I3^2-I1^2)/I3 of physical inertia tensor, [0,1]"},
    "npr1": {"unit": "1", "definition": "I1/I3 of physical inertia tensor (atomic masses), [0,1]; 1 spherical, ~0 rod (RDKit NPR1)"},
    "npr2": {"unit": "1", "definition": "I2/I3 of physical inertia tensor (atomic masses), [0,1]; high + low npr1 = rod (RDKit NPR2)"},
    "planarity": {"unit": "1", "definition": "(I3-I2)/(I3-I1) of physical inertia tensor (RDKit PMI basis); == (t2-t1)/(t3-t1) of gyration moments (exact identity), [0,1]; 1 = planar lamella"},
    "acylindricity": {"unit": "1", "definition": "(i2-i1)/(i1+i2+i3), 0 = cylinder"},
    "kappa2": {"unit": "1", "definition": "relative shape anisotropy 1.5*sum (tk-s/3)^2/s^2, [0,1]; 0 isotropic, 1 rod (DOI 10.1021/ma00148a028; identical to asphericity)"},
    "compacity": {"unit": "1", "definition": "Rg/dmax, ~0.5 extended rod"},
    "foldedness": {"unit": "1", "definition": "1 - Rg/Rg(ball of diameter dmax)"},
}


@dataclass
class DescriptorStatus:
    """Value + unit + explicit status; never a bare None in the dark."""
    value: object | None
    unit: str | None
    status: str  # "ok" | "not_assigned" | "unknown" | "error"
    reason: str | None = None


@dataclass
class Ensemble:
    path: Path
    symbols: list[list[str]]
    positions: list[np.ndarray]
    energies: list[float | None]
    source_hash: str

    def __len__(self):
        return len(self.symbols)

    @property
    def data_id(self) -> str:
        return f"{self.path.resolve()}#{self.source_hash[:16]}"


def load_ensemble(path) -> Ensemble:
    """Read GOAT (relative kcal/mol) and CREST (hartree) via existing adapters.

    The legacy `energy_Eh` result field carries an Eh-equivalent energy;
    differences for GOAT data are converted from kcal/mol and never mistaken
    for absolute hartree values. See docs/CONFORMER_GEOMETRY_DESCRIPTORS.md.
    """
    p = Path(path)
    if p.name.endswith(".finalensemble.xyz"):
        from delfin.tools.adapters.crest_ensemble import KCAL_PER_EH
        from delfin.tools.adapters.orca_goat_opi import read_goat_ensemble
        base = p.name[:-len(".finalensemble.xyz")]
        goat = read_goat_ensemble(p.parent, base)
        if not goat["structures"]:
            raise ValueError(f"no GOAT frames found in {p}")
        symbols = [[a[0] for a in frame] for frame in goat["structures"]]
        positions = [np.asarray([a[1:] for a in frame], dtype=float)
                     for frame in goat["structures"]]
        energies = [float(e) / KCAL_PER_EH
                    for e in goat["energies_kcal_mol"]]
    else:
        frames = parse_ensemble(p)
        symbols = [list(f["symbols"]) for f in frames]
        positions = [np.asarray(f["positions"], dtype=float) for f in frames]
        energies = [f["energy_Eh"] for f in frames]
    return Ensemble(p, symbols, positions, energies,
                    hashlib.sha256(p.read_bytes()).hexdigest())


# ---------------------------------------------------------------------------
# Continuous shape descriptors
# ---------------------------------------------------------------------------

def _weights(symbols: Sequence[str], weighting: str) -> np.ndarray:
    if weighting == "mass":
        # RDKit periodic table: any valid element symbol works; an
        # INVALID symbol aborts loudly instead of silently becoming carbon.
        # NOTE: RDKit's Boost exception for unknown symbols is NOT a
        # catchable python RuntimeError in all builds, so validity is
        # checked against the symbol list BEFORE any call.
        from rdkit.Chem import GetPeriodicTable
        pt = GetPeriodicTable()
        # periodic table has 118 elements (Z=1..118); collect the symbols
        valid = set()
        for z in range(1, 119):
            try:
                valid.add(pt.GetElementSymbol(z))
            except Exception:
                pass
        out = []
        for s in symbols:
            if s not in valid:
                raise ValueError(
                    f"unknown element symbol {s!r} in selection; "
                    f"cannot weight it (checked against the RDKit "
                    f"periodic table)")
            out.append(pt.GetAtomicWeight(pt.GetAtomicNumber(s)))
        return np.array(out)
    if weighting == "uniform":
        return np.ones(len(symbols))
    if weighting == "heavy":
        w = np.array([0.0 if s == "H" else 1.0 for s in symbols])
        if w.sum() == 0:
            raise ValueError("heavy weighting on an all-hydrogen selection")
        return w
    raise ValueError(f"unknown weighting {weighting!r} (mass|uniform|heavy)")


def shape_descriptors(coords, symbols, weighting: str = "mass"
                      ) -> dict[str, DescriptorStatus]:
    """Continuous shape descriptors of one geometry under one weighting."""
    coords = np.asarray(coords, dtype=float)
    try:
        w = _weights(symbols, weighting)
    except ValueError as exc:
        return {k: DescriptorStatus(None, None, "error", str(exc))
                for k in DESCRIPTOR_DOC}
    if coords.ndim != 2 or coords.shape[0] != len(w) or coords.shape[0] == 0:
        return {k: DescriptorStatus(None, None, "error",
                                    f"coords {coords.shape} vs {len(w)} atoms")
                for k in DESCRIPTOR_DOC}
    W = w.sum()
    c = (w[:, None] * coords).sum(axis=0) / W
    d = coords - c
    S = (w[:, None] * d).T @ d / W
    # Physical inertia tensor (atomic masses): I = W tr(S) E - W S.
    # RDKit NPR1/NPR2 (Sauer/Schwarz 2003, useAtomicMasses=True) and the
    # molecular eccentricity are DEFINED on these principal moments;
    # the asphericity of Baumgaertner (RDKit eq. 11) is defined on the
    # moments of the gyration matrix (= eigenvalues of S).
    It = W * (np.trace(S) * np.eye(3) - S)
    evI = np.sort(np.linalg.eigvalsh(It))
    I1, I2, I3 = (float(x) for x in evI)
    ev = np.sort(np.linalg.eigvalsh(S))
    t1, t2, t3 = (float(x) for x in ev)

    diff = coords[:, None, :] - coords[None, :, :]
    dmax = float(np.sqrt((diff ** 2).sum(-1)).max())
    tr = t1 + t2 + t3
    rg = math.sqrt(max(tr, 0.0))
    denom = 2.0 * (t1 + t2 + t3) ** 2
    asph = (((t3 - t2) ** 2 + (t3 - t1) ** 2 + (t2 - t1) ** 2) / denom
            if denom > 0 else None)
    ecc = (math.sqrt(max(0.0, (I3 ** 2 - I1 ** 2))) / I3 if I3 > 0 else None)
    npr1 = I1 / I3 if I3 > 0 else None
    npr2 = I2 / I3 if I3 > 0 else None
    plan = (t2 - t1) / (t3 - t1) if t3 > t1 else None
    acyl = (t2 - t1) / tr if tr > 0 else None
    if tr > 0:
        kappa2 = 1.5 * sum((ik - tr / 3.0) ** 2 for ik in ev) / tr ** 2
    else:
        kappa2 = None
    comp = rg / dmax if dmax > 0 else None
    folded = (1.0 - rg / ((dmax / 2.0) * math.sqrt(3.0 / 5.0))
              if dmax > 0 else None)

    def st(x, unit="1"):
        return DescriptorStatus(None if x is None else round(float(x), 6),
                                unit, "ok")

    return {
        "rg": st(rg, "A"), "dmax": st(dmax, "A"),
        "asphericity": st(asph), "eccentricity": st(ecc),
        "npr1": st(npr1), "npr2": st(npr2),
        "planarity": st(plan), "acylindricity": st(acyl),
        "kappa2": st(kappa2), "compacity": st(comp),
        "foldedness": st(folded),
    }


# ---------------------------------------------------------------------------
# Connectivity / chemical classification helpers
# ---------------------------------------------------------------------------

def _conn_mol(symbols: Sequence[str], coords: np.ndarray):
    """Connectivity-only RDKit mol from one geometry (no bond orders)."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError("rdkit required for chemical classification") from exc
    block = f"{len(symbols)}\n\n"
    block += "\n".join(f"{s} {x:.6f} {y:.6f} {z:.6f}"
                       for s, (x, y, z) in zip(symbols, coords)) + "\n"
    try:
        from rdkit.Chem import rdDetermineBonds
    except ImportError:
        rdDetermineBonds = None
    mol = Chem.MolFromXYZBlock(block)
    if mol is None:
        raise ValueError("RDKit could not parse the geometry")
    if rdDetermineBonds is not None:
        rdDetermineBonds.DetermineConnectivity(mol)  # connectivity, no bond orders
    Chem.FastFindRings(mol)  # ring perception fallback without bond orders
    try:
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_SETCONJUGATION)
    except Exception:
        pass  # conjugation flags are best-effort; connectivity stands
    return mol


def _dihedral(p0, p1, p2, p3) -> float:
    """Dihedral in degrees, or math.nan when it is geometrically undefined
    (coincident adjacent atoms make the rotation axis degenerate). Callers
    that present the value to a user must translate NaN into an explicit
    error/not_assigned status -- a NaN must never be shown as a number."""
    b0 = np.asarray(p0) - np.asarray(p1)
    b1 = np.asarray(p2) - np.asarray(p1)
    b2 = np.asarray(p3) - np.asarray(p2)
    norm_b1 = np.linalg.norm(b1)
    if norm_b1 == 0.0:
        return math.nan
    b1n = b1 / norm_b1
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))


def _sp2_like(conn_mol, coords, idx) -> bool:
    """sp2-like carbon by geometry: 2 heavy neighbours, angle 110-130 deg."""
    at = conn_mol.GetAtomWithIdx(idx)
    if at.GetSymbol() != "C":
        return False
    heavy = [n.GetIdx() for n in at.GetNeighbors() if n.GetSymbol() != "H"]
    if len(heavy) != 2:
        return False
    b1 = coords[heavy[0]] - coords[idx]
    b2 = coords[heavy[1]] - coords[idx]
    cosang = np.dot(b1, b2) / (np.linalg.norm(b1) * np.linalg.norm(b2))
    ang = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
    return 110.0 <= ang <= 130.0


def double_bond_geometry(conn_mol, coords) -> list[dict[str, object]]:
    """Geometric E/Z-like assignment across acyclic C=C-like bonds.

    For every acyclic bond between two sp2-like carbons, take the first
    heavy neighbour on each side (hydrogen neighbours are SKIPPED: the
    reference atom is the first non-hydrogen substituent, so a terminal
    =CH2 carries no steric reference and the bond yields no dihedral)
    and compute the dihedral across the bond. |dih| < 45 deg ->
    "Z-like", > 135 deg -> "E-like", else status "not_assigned".
    This is GEOMETRIC ONLY: without bond orders and CIP priorities no
    true stereochemical configuration is proven.
    """
    from rdkit import Chem
    Chem.FastFindRings(conn_mol)
    ri = conn_mol.GetRingInfo()
    out = []
    for bond in conn_mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if ri.IsBondInRingOfSize(bond.GetIdx(), 6):
            continue
        if not (_sp2_like(conn_mol, coords, a)
                and _sp2_like(conn_mol, coords, b)):
            continue
        na = [n.GetIdx() for n in conn_mol.GetAtomWithIdx(a).GetNeighbors()
              if n.GetIdx() != b and n.GetSymbol() != "H"]
        nb = [n.GetIdx() for n in conn_mol.GetAtomWithIdx(b).GetNeighbors()
              if n.GetIdx() != a and n.GetSymbol() != "H"]
        if not na or not nb:
            continue
        d = _dihedral(coords[na[0]], coords[a], coords[b], coords[nb[0]])
        ang = abs(d)
        if ang < 45.0:
            stereo, status, reason = "Z-like", "assigned", None
        elif ang > 135.0:
            stereo, status, reason = "E-like", "assigned", None
        else:
            stereo, status, reason = None, "not_assigned", (
                f"dihedral {d:.1f} deg in ambiguous range 45-135")
        out.append({"bond": [a, b], "dihedral_deg": round(d, 1),
                    "stereo": stereo, "status": status, "reason": reason,
                    "kind": "geometric (no CIP, no bond orders)"})
    return out


def diene_conformation(conn_mol, coords) -> list[dict[str, object]]:
    """s-cis / s-trans state across central sp2-sp2 single bonds.

    Dihedral of the outer double-bond partners across the bond:
    |dih| < 60 -> s-cis, > 120 -> s-trans, else not_assigned.
    """
    from rdkit import Chem
    Chem.FastFindRings(conn_mol)
    ri = conn_mol.GetRingInfo()
    out = []
    for bond in conn_mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if ri.IsBondInRingOfSize(bond.GetIdx(), 6):
            continue
        if not (_sp2_like(conn_mol, coords, a)
                and _sp2_like(conn_mol, coords, b)):
            continue
        na = [n.GetIdx() for n in conn_mol.GetAtomWithIdx(a).GetNeighbors()
              if n.GetIdx() != b and n.GetSymbol() != "H"]
        nb = [n.GetIdx() for n in conn_mol.GetAtomWithIdx(b).GetNeighbors()
              if n.GetIdx() != a and n.GetSymbol() != "H"]
        if not na or not nb:
            continue
        d = _dihedral(coords[na[0]], coords[a], coords[b], coords[nb[0]])
        ang = abs(d)
        if ang < 60.0:
            state, status = "s-cis", "assigned"
        elif ang > 120.0:
            state, status = "s-trans", "assigned"
        else:
            state, status = None, "not_assigned"
        out.append({"bond": [a, b], "dihedral_deg": round(d, 1),
                    "state": state, "status": status})
    return out


def ez_cip_assignment(symbols, coords, charge=0) -> list[dict[str, object]]:
    """Chemically defined E/Z for stereogenic C=C from bond orders + CIP.

    Uses rdDetermineBonds.DetermineBonds(charge) for BOND ORDERS and
    Chem.AssignStereochemistryFrom3D for CIP priorities. This is a
    RECONSTRUCTION from coordinates (no independent reference structure);
    it is chemically defined (bond orders + CIP) but not externally
    validated. Aromatic-ring 1.5-bonds are not stereogenic and are
    excluded by RDKit itself (STEREOE/Z only on real alkenes).

    Status per double bond:
      "assigned"   -- CIP stereo flags were set for this bond
      "not_assigned"-- RDKit left the bond stereovague (e.g. ring
                       fusion, symmetric substitution)
      "unknown"    -- bond-order perception failed entirely
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError as exc:
        raise ImportError("rdkit required for CIP E/Z assignment") from exc
    block = f"{len(symbols)}\n\n"
    block += "\n".join(f"{s} {x:.6f} {y:.6f} {z:.6f}"
                       for s, (x, y, z) in zip(symbols, coords)) + "\n"
    mol = Chem.MolFromXYZBlock(block)
    if mol is None:
        return [{"bond": None, "stereo": None, "status": "unknown",
                 "reason": "XYZ block could not be parsed by RDKit",
                  "kind": "coordinate reconstruction, not independently validated"}]
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=charge)
    except Exception as exc:  # RDKit raises Boost exceptions, not python
        return [{"bond": None, "stereo": None, "status": "unknown",
                 "reason": f"bond-order perception failed: {exc}",
                  "kind": "coordinate reconstruction, not independently validated"}]
    try:
        Chem.AssignStereochemistryFrom3D(mol)
    except Exception as exc:
        return [{"bond": None, "stereo": None, "status": "unknown",
                 "reason": f"stereochemistry assignment failed: {exc}",
                 "kind": "coordinate reconstruction, not independently validated"}]
    out = []
    for b in mol.GetBonds():
        if b.GetBondTypeAsDouble() <= 1.5:
            continue
        if b.GetBeginAtom().GetSymbol() != "C" or \
                b.GetEndAtom().GetSymbol() != "C":
            continue
        stereo = str(b.GetStereo())
        a, bb = sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        if stereo in ("STEREOE", "STEREOZ"):
            out.append({"bond": [a, bb],
                        "stereo": "E" if stereo == "STEREOE" else "Z",
                        "status": "assigned", "reason": None,
                        "kind": "coordinate reconstruction, not independently validated"})
        else:
            out.append({"bond": [a, bb], "stereo": None,
                        "status": "not_assigned",
                        "reason": f"RDKit stereo state {stereo}",
                        "kind": "coordinate reconstruction, not independently validated"})
    if not out:
        out.append({"bond": None, "stereo": None, "status": "not_assigned",
                    "reason": "no acyclic C=C in the perceived structure",
                    "kind": "coordinate reconstruction, not independently validated"})
    return out


def rotamer_states(conn_mol, coords, reference=None) -> list[dict[str, object]]:
    """Torsions about rotatable non-ring heavy-atom single bonds.

    Optional reference geometry adds the circular difference (in
    [-180, 180]) per torsion, so conformer-vs-reference rotamer shifts
    are readable directly.
    """
    out = []
    for bond in conn_mol.GetBonds():
        if bond.GetIsInRing():
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if conn_mol.GetAtomWithIdx(a).GetSymbol() == "H" or \
           conn_mol.GetAtomWithIdx(b).GetSymbol() == "H":
            continue
        na = [n.GetIdx() for n in conn_mol.GetAtomWithIdx(a).GetNeighbors()
              if n.GetIdx() != b and n.GetSymbol() != "H"]
        nb = [n.GetIdx() for n in conn_mol.GetAtomWithIdx(b).GetNeighbors()
              if n.GetIdx() != a and n.GetSymbol() != "H"]
        if not na or not nb:
            continue
        d = _dihedral(coords[na[0]], coords[a], coords[b], coords[nb[0]])
        row = {"bond": [a, b], "dihedral_deg": round(d, 1),
               "status": "assigned"}
        if reference is not None:
            rd = _dihedral(reference[na[0]], reference[a], reference[b],
                           reference[nb[0]])
            diff = (d - rd + 180.0) % 360.0 - 180.0
            row["delta_vs_reference_deg"] = round(diff, 1)
        out.append(row)
    return out


def custom_torsion(coords, i0, i1, i2, i3) -> DescriptorStatus:
    """User-specified dihedral of four atom indices."""
    n = len(coords)
    if not all(0 <= i < n for i in (i0, i1, i2, i3)) or len({i0, i1, i2, i3}) != 4:
        return DescriptorStatus(None, "deg", "error",
                                "torsion needs 4 distinct valid atom indices")
    value = _dihedral(coords[i0], coords[i1], coords[i2], coords[i3])
    if math.isnan(value):
        return DescriptorStatus(None, "deg", "error",
                                "torsion is undefined: coincident atoms "
                                "i1 and i2 give a degenerate rotation axis")
    return DescriptorStatus(round(value, 1), "deg", "ok")


def fragment_distance(coords, sel_a, sel_b) -> DescriptorStatus:
    """Minimum inter-fragment heavy-atom distance in Angstrom."""
    a = np.asarray(coords)[list(sel_a)]
    b = np.asarray(coords)[list(sel_b)]
    if a.size == 0 or b.size == 0:
        return DescriptorStatus(None, "A", "error", "empty fragment selection")
    diff = a[:, None, :] - b[None, :, :]
    return DescriptorStatus(round(float(np.sqrt((diff ** 2).sum(-1)).min()), 4),
                            "A", "ok")


def fragment_orientation(coords, symbols, sel_a, sel_b) -> DescriptorStatus:
    """Angle (deg) between fragment principal axes (i3 eigenvectors).

    Axis of a fragment = eigenvector of its weighted covariance with the
    LARGEST eigenvalue (longest principal extent). Returns the angle
    between the two axes in [0, 90] deg; near 0/180 = parallel (stacked
    or co-aligned), near 90 = perpendicular.
    """
    try:
        w = _weights([symbols[i] for i in sel_a], "mass")
    except ValueError as exc:
        return DescriptorStatus(None, "deg", "error", str(exc))
    ca = np.asarray(coords)[list(sel_a)]
    cb = np.asarray(coords)[list(sel_b)]
    if len(sel_a) < 3 or len(sel_b) < 3:
        return DescriptorStatus(None, "deg", "not_assigned",
                                "axis needs >= 3 atoms per fragment")

    def axis(c, w):
        W = w.sum()
        d = c - (w[:, None] * c).sum(axis=0) / W
        S = (w[:, None] * d).T @ d / W
        vals, vecs = np.linalg.eigh(S)
        return vecs[:, int(np.argmax(vals))]

    u, v = axis(ca, w), axis(cb, w)
    ang = math.degrees(math.acos(min(1.0, abs(float(np.dot(u, v))))))
    return DescriptorStatus(round(ang, 1), "deg", "ok")


# ---------------------------------------------------------------------------
# Selection parsing, ranking, filtering
# ---------------------------------------------------------------------------

def parse_selection(spec: str, n_atoms: int) -> list[int]:
    """"0-5,9,12-14" -> sorted unique index list (0-based)."""
    out = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-")
            idx = range(int(lo), int(hi) + 1)
        else:
            idx = [int(part)]
        for i in idx:
            if not 0 <= i < n_atoms:
                raise ValueError(f"atom index {i} out of range 0..{n_atoms - 1}")
            out.add(i)
    if not out:
        raise ValueError("empty atom selection")
    return sorted(out)


def rank_conformers(values, reverse: bool = False) -> list[int]:
    """Deterministic ranking: ties broken by lower conformer index.

    None values sort last regardless of direction.
    """
    have = [(v, i) for i, v in enumerate(values) if v is not None]
    have.sort(key=lambda t: t[1])          # ties keep ascending index order
    have.sort(key=lambda t: t[0], reverse=reverse)  # stable on value only
    none_idx = [i for i, v in enumerate(values) if v is None]
    return [i for _, i in have] + none_idx


def filter_conformers(rows: Sequence[dict[str, object]],
                      predicate) -> list[int]:
    """Indices of rows matching predicate; a conformer may match many."""
    return [i for i, row in enumerate(rows) if predicate(row)]


# ---------------------------------------------------------------------------
# Whole-ensemble analysis, JSON, cache
# ---------------------------------------------------------------------------

def cache_key(path: Path, weighting: str, selections: Sequence[str],
              definition_version: str = "6") -> str:
    """Cache identity: source hash + weighting + selections + version.

    Version 6 invalidates GOAT energy fields parsed with CREST hartree rules;
    Version 2 invalidates every cache entry written before the inertia-
    tensor correction (NPR1/NPR2/eccentricity/asphericity definitions);
    version 3 invalidates entries written before the kappa2 (relative
    shape anisotropy) correction; version 5 invalidates entries written
    before the CIP-based E/Z assignment (ez_cip) and the ring-bond
    heuristic fix.
    """
    h = hashlib.sha256()
    h.update(hashlib.sha256(Path(path).read_bytes()).digest())
    h.update(weighting.encode())
    h.update("|".join(sorted(set(selections))).encode())
    h.update(definition_version.encode())
    return h.hexdigest()


def load_or_compute(path, weighting="mass", selections=None, cache_dir=None):
    """Load cached JSON result if present, else compute and cache it."""
    sels = selections or ["all"]
    key = cache_key(path, weighting, sels)
    cache_path = None
    if cache_dir:
        cache_path = Path(cache_dir) / f"shape-classification-{key[:24]}.json"
        if cache_path.exists():
            return json.loads(cache_path.read_text()), cache_path
    result = analyze_ensemble(path, weighting=weighting, selections=sels)
    result["cache_key"] = key
    if cache_path is not None:
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(result, indent=1))
        except OSError:
            # Read-only cache dir (HPC home): results stay valid, just not
            # persisted. Never fail the analysis because of the cache.
            cache_path = None
    return result, cache_path


def analyze_ensemble(path, weighting="mass", selections=None) -> dict:
    """Full analysis: descriptors for all frames and selections + classes.

    The input coordinates are never modified; each analysis works on the
    loaded copies. Frame indices in the output are the original ones.
    """
    ens = load_ensemble(path)
    sels = selections or ["all"]
    sel_map: dict[str, list[int]] = {"all": list(range(len(ens.symbols[0])))}
    for spec in sels:
        if spec != "all":
            sel_map[spec] = parse_selection(spec, len(ens.symbols[0]))

    rows = []
    for k in range(len(ens)):
        coords = ens.positions[k]
        symbols = ens.symbols[k]
        entry = {
            "conformer_id": k,
            "frame_index": k,
            "energy_Eh": ens.energies[k],
            "descriptors": {},
            "chemical_classes": {},
        }
        for name, sel in sel_map.items():
            sc = coords[sel] if name != "all" else coords
            sy = [symbols[i] for i in sel] if name != "all" else symbols
            entry["descriptors"][name] = {
                key: asdict(v) for key, v in
                shape_descriptors(sc, sy, weighting).items()}
            if name == "all":
                mol = _conn_mol(symbols, coords)
                entry["chemical_classes"]["double_bonds"] = \
                    double_bond_geometry(mol, coords)
                entry["chemical_classes"]["diene_conformations"] = \
                    diene_conformation(mol, coords)
                entry["chemical_classes"]["ez_cip"] = \
                    ez_cip_assignment(symbols, coords)
        rows.append(entry)

    return {
        "input_file": str(Path(path).resolve()),
        "data_id": ens.data_id,
        "source_sha256": ens.source_hash,
        "weighting": weighting,
        "selections": {name: sel for name, sel in sel_map.items()},
        "descriptor_definitions": DESCRIPTOR_DOC,
        "caveats": [
            "XYZ has no bond orders: E/Z-like and s-cis/s-trans are geometric; "
            "CIP E/Z is coordinate-reconstructed, not independently validated.",
            "Descriptors from different weightings are not comparable.",
            "No thresholds are validated; classification into bins is the "
            "caller's responsibility.",
        ],
        "frames": rows,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Conformer shape classification (see module docstring)")
    ap.add_argument("ensemble", help="multi-frame XYZ file")
    ap.add_argument("--json", dest="json_out", help="write JSON result here")
    ap.add_argument("--weights", default="mass",
                    choices=["mass", "uniform", "heavy"])
    ap.add_argument("--atoms", action="append", default=[],
                    help="atom selection spec like 0-5,9; repeatable")
    ap.add_argument("--cache", help="cache directory for results")
    ap.add_argument("--ranking", action="append", default=[],
                    metavar="KEY[:asc|:desc]",
                    help="rank frames by descriptor KEY in selection 'all'")
    args = ap.parse_args(argv)

    selections = ["all"] + list(args.atoms)
    if args.cache:
        result, _ = load_or_compute(args.ensemble, args.weights,
                                    selections, args.cache)
    else:
        result = analyze_ensemble(args.ensemble, weighting=args.weights,
                                  selections=selections)
    for spec in args.ranking:
        key, _, order = spec.partition(":")
        vals = [f["descriptors"]["all"].get(key, {}).get("value")
                for f in result["frames"]]
        idx = rank_conformers(vals, reverse=(order == "desc"))
        result.setdefault("rankings", {})[spec] = idx
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(result, indent=1))
    else:
        print(json.dumps(result["frames"][0]["descriptors"]["all"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
