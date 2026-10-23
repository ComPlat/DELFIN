"""ERDBEBEN final pass — re-seat a planar-COLLAPSED ligand fragment from a clean ISOLATED embed.

The decisive finding (2026-07-23): the metal-context whole-complex ETKDG collapses rigid ligand cages
(AQIBAE bis-phosphite: one cage flattens to a plane, sp3 tetra-volume ~0), while the SAME ligand embedded
in ISOLATION from the metal builds 3D in 20/20 ETKDG seeds.  A spring-relax (``_bond_decollapse``) only
spreads the superimposed atoms (sev 5.8->4.6) but cannot pop the flattened cage out of its planar local
minimum -- only a FRESH embed reconstructs the tetrahedral geometry.

This runs LAST (on the final frame list, after every builder/UFF/bandage), so it catches a collapsed
fragment no matter which builder produced it.  Per frame:
  1. perceive non-metal fragments (geometric bonds);
  2. for a fragment with a planar-collapsed sp3 centre, embed that fragment ALONE (clean, 3D);
  3. Kabsch-align the clean embed onto the frame's metal-anchored DONOR atoms (they are correctly placed),
     resolving the residual rotational DOF (2-donor chelate) by minimising clash with the rest of the
     complex;
  4. splice the clean fragment coordinates (heavy + H) back in.
Per-frame ROLLBACK: keep only if the collapse strictly drops AND no metal-donor bond breaks AND
inter-atom clashes do not increase.  Bit-exact when no fragment is collapsed (early return) or the flag is
off.  Requires frame atom order == mol atom order (DELFIN convention; guarded by an atom-count match, else
the frame is skipped).  RDKit-only, no CCDC -- the fragment topology comes from ``mol`` (SMILES)."""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional

import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    _HAVE_RDKIT = True
except Exception:                                                 # pragma: no cover
    Chem = None
    _HAVE_RDKIT = False

_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Mg", "Ca", "Sr", "Ba", "Al", "Ga", "In", "Tl", "Sn", "Pb", "Bi", "Sb",
    "Li", "Na", "K", "Rb", "Cs", "Be",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Th", "U",
}
_COV = {"H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Si": 1.11, "P": 1.07,
        "S": 1.05, "Cl": 1.02, "As": 1.19, "Br": 1.20, "I": 1.39, "Se": 1.20, "Te": 1.38}

_VOL_MIN = 0.40          # sp3 tetra-volume floor (collapsed below this)
_MD_FACTOR = 1.35        # metal-donor bond detection factor (x covalent sum)
_MD_TOL = 0.25           # allowed M-D lengthening after re-seat (A) before rollback
_CLASH_MIN = 2.2         # heavy-heavy non-bonded clash floor (A)


def _cov(s: str) -> float:
    return _COV.get(s, 0.95)


def _parse(xyz: str):
    """Parse an xyz block to (syms, P).  ROBUST: skips non-atom lines (a leading count line, or an
    'AQIBAE frame0 pucker' comment line) so it works whether or not the frame carries an xyz header --
    the pipeline passes full frames, the self-test passes atom-blocks."""
    syms: List[str] = []
    pts: List[List[float]] = []
    for ln in xyz.strip().split("\n"):
        p = ln.split()
        if len(p) >= 4:
            try:
                x, y, z = float(p[1]), float(p[2]), float(p[3])
            except (ValueError, IndexError):
                continue                                          # comment / non-coordinate line
            syms.append(p[0])
            pts.append([x, y, z])
    return syms, np.array(pts, dtype=float)


def _write(syms: List[str], P: np.ndarray) -> str:
    return "\n".join("%-4s %12.6f %12.6f %12.6f" % (syms[i], P[i, 0], P[i, 1], P[i, 2])
                     for i in range(len(syms))) + "\n"


def _adjacency(syms: List[str], P: np.ndarray):
    n = len(syms)
    adj: List[List[int]] = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(P[i] - P[j]))
            if 0.3 < d < 1.3 * (_cov(syms[i]) + _cov(syms[j])):
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _nonmetal_fragments(syms: List[str], adj) -> List[List[int]]:
    n = len(syms)
    seen: set = set()
    frags: List[List[int]] = []
    for s in range(n):
        if s in seen or syms[s] in _METALS:
            continue
        comp: List[int] = []
        stack = [s]
        while stack:
            u = stack.pop()
            if u in seen or syms[u] in _METALS:
                continue
            seen.add(u)
            comp.append(u)
            for v in adj[u]:
                if v not in seen and syms[v] not in _METALS:
                    stack.append(v)
        frags.append(sorted(comp))
    return frags


def _mol_fragments(mol, n_atoms: int) -> List[List[int]]:
    """Non-metal connected components in the TOPOLOGY (mol) -- gives the COMPLETE ligand fragment (e.g. a
    bis-phosphite = both cages + linker as one) regardless of the collapsed geometry, which a geometric
    fragmentation would split into disconnected pieces (the '0 donors' failure).  mol index == frame index."""
    seen: set = set()
    frags: List[List[int]] = []
    for s in range(n_atoms):
        a = mol.GetAtomWithIdx(s)
        if s in seen or a.GetSymbol() in _METALS:
            continue
        comp: List[int] = []
        stack = [s]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            au = mol.GetAtomWithIdx(u)
            if au.GetSymbol() in _METALS:
                continue
            seen.add(u)
            comp.append(u)
            for nb in au.GetNeighbors():
                if nb.GetIdx() not in seen and nb.GetSymbol() not in _METALS:
                    stack.append(nb.GetIdx())
        frags.append(sorted(comp))
    return frags


def _sp3_nbrs(mol, ai: int, fs: set, syms: List[str], P: np.ndarray) -> List[int]:
    """Bonded neighbours of ai that are ALSO within bonding distance (GEOMETRIC).  Matches the eye's
    find_planar_collapse (which uses geometric adjacency, not mol topology): a centre with a STRETCHED bond
    has <4 geometric neighbours -> NOT a collapse.  Aligning the reseat trigger with the eye stops the
    false-fire on e.g. EKAKIK's telluroether-crown centre (a stretched-bond sp3 the eye considers fine)."""
    out = []
    for nbj in mol.GetAtomWithIdx(ai).GetNeighbors():
        j = nbj.GetIdx()
        if j in fs and 0.3 < float(np.linalg.norm(P[j] - P[ai])) < 1.3 * (_cov(syms[ai]) + _cov(syms[j])):
            out.append(j)
    return out


def _fragment_collapsed(mol, frag: List[int], syms: List[str], P: np.ndarray) -> bool:
    """True if the fragment (frame indices == mol indices) has a planar-collapsed sp3 centre."""
    fs = set(frag)
    for ai in frag:
        a = mol.GetAtomWithIdx(ai)
        if a.GetSymbol() in _METALS or a.GetIsAromatic():
            continue
        if a.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        nb = _sp3_nbrs(mol, ai, fs, syms, P)
        if len(nb) < 4 or sum(1 for j in nb if syms[j] != "H") < 2:
            continue
        near = sorted(nb, key=lambda j: float(np.sum((P[j] - P[ai]) ** 2)))[:4]
        v = [P[j] for j in near]
        vol = abs(float(np.dot(np.cross(v[1] - v[0], v[2] - v[0]), v[3] - v[0]))) / 6.0
        if vol < _VOL_MIN:
            return True
    return False


def _submol_collapsed(m, conf) -> bool:
    for a in m.GetAtoms():
        if a.GetIsAromatic() or a.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        nb = list(a.GetNeighbors())
        if len(nb) < 4 or sum(1 for x in nb if x.GetSymbol() != "H") < 2:
            continue
        pa = conf.GetAtomPosition(a.GetIdx())
        near = sorted(nb, key=lambda x: (conf.GetAtomPosition(x.GetIdx()).x - pa.x) ** 2
                      + (conf.GetAtomPosition(x.GetIdx()).y - pa.y) ** 2
                      + (conf.GetAtomPosition(x.GetIdx()).z - pa.z) ** 2)[:4]
        v = [np.array([conf.GetAtomPosition(x.GetIdx()).x, conf.GetAtomPosition(x.GetIdx()).y,
                       conf.GetAtomPosition(x.GetIdx()).z]) for x in near]
        vol = abs(float(np.dot(np.cross(v[1] - v[0], v[2] - v[0]), v[3] - v[0]))) / 6.0
        if vol < _VOL_MIN:
            return True
    return False


def _collapsed_free_set(mol, frag: List[int], syms: List[str], orig_P: np.ndarray) -> set:
    """Frag atoms in the COLLAPSE region -- planar-collapsed sp3 centres + their bonded substituents.
    GRAFT re-embed (2026-07-25, user chose the by-construction root fix): ONLY these embed freely; every
    other fragment atom is constrained (coordMap) to its ORIGINAL position, so the CCDC-anchored distal
    geometry (biaryl twist, pyramid, ring fold, poly) is PRESERVED while the collapsed cage pops."""
    fs = set(frag)
    free: set = set()
    for ai in frag:
        a = mol.GetAtomWithIdx(ai)
        if a.GetSymbol() in _METALS or a.GetIsAromatic() or a.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        nb = _sp3_nbrs(mol, ai, fs, syms, orig_P)                 # geometric-aligned (see _sp3_nbrs)
        if len(nb) < 4 or sum(1 for j in nb if syms[j] != "H") < 2:
            continue
        near = sorted(nb, key=lambda j: float(np.sum((orig_P[j] - orig_P[ai]) ** 2)))[:4]
        v = [orig_P[j] for j in near]
        vol = abs(float(np.dot(np.cross(v[1] - v[0], v[2] - v[0]), v[3] - v[0]))) / 6.0
        if vol < _VOL_MIN:
            free.add(ai)
            free.update(nb)                       # centre + substituents must move to become tetrahedral
    return free


def _embed_fragment_clean(mol, frag: List[int], donors: List[int], target_donor_P: np.ndarray,
                          n_seeds: int = 20, orig_P: Optional[np.ndarray] = None,
                          syms: Optional[List[str]] = None) -> Optional[Dict[int, np.ndarray]]:
    """Extract the fragment (with its explicit H) as a sub-mol, embed it ISOLATED, return {frame_idx: xyz}
    for the NON-COLLAPSED embed whose donor-donor geometry best matches the frame's (so the Kabsch
    alignment onto the metal-anchored donors keeps the M-D bite -- the metal enforces a bite the isolated
    ligand does not naturally hold, so a random fold breaks M-D; we pick the fold that matches).  None if
    all embeds fail/collapse."""
    amap: Dict[int, int] = {}
    rw = Chem.RWMol()
    for i in frag:
        a = mol.GetAtomWithIdx(i)
        na = Chem.Atom(a.GetAtomicNum())
        na.SetFormalCharge(a.GetFormalCharge())
        na.SetNoImplicit(True)
        amap[i] = rw.AddAtom(na)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in amap and j in amap:
            rw.AddBond(amap[i], amap[j], b.GetBondType())
    sub0 = rw.GetMol()
    sub = None
    for _neut in (False, True):                                   # try as-is, then neutralise dative [P+]/[N+]
        m = Chem.RWMol(sub0)
        if _neut:
            for a in m.GetAtoms():
                a.SetFormalCharge(0)
        mm = m.GetMol()
        try:
            Chem.SanitizeMol(mm)
            sub = mm
            break
        except Exception:
            try:
                Chem.SanitizeMol(mm, sanitizeOps=(Chem.SanitizeFlags.SANITIZE_ALL
                                                  & ~Chem.SanitizeFlags.SANITIZE_KEKULIZE
                                                  & ~Chem.SanitizeFlags.SANITIZE_PROPERTIES))
                sub = mm
                break
            except Exception:
                continue
    if sub is None:
        return None
    inv = {v: k for k, v in amap.items()}
    # GRAFT (user root fix 2026-07-25): constrain the NON-collapsed fragment atoms to their ORIGINAL
    # positions so the embed PRESERVES the CCDC-anchored distal geometry (biaryl/pyramid/fold/poly) and
    # only pops the collapse.  None -> full fresh embed (old behaviour); every failure below falls back to it.
    coord_map = None
    if orig_P is not None and syms is not None:
        try:
            from rdkit.Geometry import Point3D
            _free = _collapsed_free_set(mol, frag, syms, orig_P)
            if 0 < len(_free) < len(frag):
                coord_map = {amap[fi]: Point3D(float(orig_P[fi][0]), float(orig_P[fi][1]), float(orig_P[fi][2]))
                             for fi in frag if fi not in _free and fi in amap}
                if len(coord_map) < 3:            # too few anchors -> not a meaningful constraint
                    coord_map = None
        except Exception:
            coord_map = None
    tgt_dd = None
    if len(donors) >= 2:
        tgt_dd = np.linalg.norm(target_donor_P[:, None, :] - target_donor_P[None, :, :], axis=-1)
    best_map = None
    best_score = 1e18
    for seed in range(1, n_seeds + 1):
        m2 = Chem.Mol(sub)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        try:
            _ok = -1
            if coord_map:
                try:
                    _ok = AllChem.EmbedMolecule(m2, params, coordMap=coord_map)   # GRAFT: preserve distal geom
                except TypeError:
                    _ok = AllChem.EmbedMolecule(m2, randomSeed=seed, coordMap=coord_map)
                if _ok != 0:                                                       # graft failed -> fresh embed
                    m2 = Chem.Mol(sub)
                    _ok = AllChem.EmbedMolecule(m2, params)
            else:
                _ok = AllChem.EmbedMolecule(m2, params)
            if _ok != 0:
                continue
            conf = m2.GetConformer()
        except Exception:
            continue
        if _submol_collapsed(m2, conf):
            continue
        cmap = {inv[si]: np.array([conf.GetAtomPosition(si).x, conf.GetAtomPosition(si).y,
                                   conf.GetAtomPosition(si).z], float)
                for si in range(sub.GetNumAtoms()) if si in inv}
        if tgt_dd is not None:                                    # score by donor-donor bite match
            dp = np.array([cmap[d] for d in donors])
            edd = np.linalg.norm(dp[:, None, :] - dp[None, :, :], axis=-1)
            score = float(np.sqrt(np.mean((edd - tgt_dd) ** 2)))
        else:
            score = 0.0
        if score < best_score:
            best_score = score
            best_map = cmap
            if best_score < 0.25:                                # good enough -> stop early
                break
    return best_map


def _kabsch(A: np.ndarray, B: np.ndarray):
    """Rotation R + translation t so that R@A_i + t best matches B_i (least squares)."""
    ca, cb = A.mean(0), B.mean(0)
    H = (A - ca).T @ (B - cb)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cb - R @ ca


def _metal_donors(mol, frag, syms, P, metal_idxs):
    """Fragment atoms bonded to a metal.  Uses the TOPOLOGY (mol) first -- robust to a collapse that pushed
    a donor away from the metal (the geometric M-D distance then misses it, the '0 donors' failure mode) --
    with a geometric fallback."""
    metal_set = set(metal_idxs)
    out = [i for i in frag if mol.GetAtomWithIdx(i).GetAtomicNum() != 1
           and any(nb.GetIdx() in metal_set for nb in mol.GetAtomWithIdx(i).GetNeighbors())]
    if out:
        return out
    for i in frag:                                                # geometric fallback
        if syms[i] == "H":
            continue
        for m in metal_idxs:
            if float(np.linalg.norm(P[i] - P[m])) < _MD_FACTOR * (_cov(syms[i]) + _cov(syms[m])):
                out.append(i)
                break
    return out


def _clash_count(syms, P, frag_set, exclude_pairs=None) -> int:
    """Heavy-heavy contacts < _CLASH_MIN between the fragment and the REST (inter-fragment clash)."""
    other = [k for k in range(len(syms)) if k not in frag_set and syms[k] != "H" and syms[k] not in _METALS]
    if not other:
        return 0
    Q = P[other]
    n = 0
    for i in frag_set:
        if syms[i] == "H":
            continue
        d = np.linalg.norm(Q - P[i], axis=1)
        n += int(np.sum(d < _CLASH_MIN))
    return n


def _pyramidalisation(center, neighbors) -> float:
    """360 - sum of the three X-centre-X angles (deg).  ~0 = planar sp2, ~31.5 = ideal tetrahedral."""
    c = np.asarray(center, float)
    ns = [np.asarray(x, float) for x in neighbors]
    ssum = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            v1 = ns[i] - c
            v2 = ns[j] - c
            dn = float(np.linalg.norm(v1) * np.linalg.norm(v2))
            if dn < 1e-8:
                return 0.0
            cosang = max(-1.0, min(1.0, float(np.dot(v1, v2)) / dn))
            ssum += math.degrees(math.acos(cosang))
    return 360.0 - ssum


def _worsens_planarity(mol, frag, syms, P_orig, P_cand,
                       planar_thresh: float = 15.0, worsen_margin: float = 15.0) -> bool:
    """ROLLBACK GUARD (2026-07-24): True if the candidate PYRAMIDALISES a fragment centre that the TOPOLOGY
    says is planar (sp2 / aromatic) and that WAS planar in the original.  The isolated re-embed lacks the
    metal context that keeps a metal-induced-planar donor flat (amide/imine N, sp2 C, cyclometalated aryl
    carbanion) and can build it pyramidal -- a hybridisation regression the collapse/M-D/clash rollback
    misses (FECJIJ/JIMTIM/MAKKOC, erdbeben10k).  GEOMETRIC criterion (NOT RDKit hybridisation, which the
    isolated embed shares and mis-perceives for metal-bonded/carbanion centres -> JIMTIM): guard ANY
    3-coordinate centre that was PLANAR in the original.  4-coordinate sp3 collapse -- the reseat's actual
    job -- is excluded by len(nb)==3, so this never blocks the collapse fix.  Over-blocking (skipping a
    legit pyramidalisation) only reduces coverage, never regresses = never-worse-safe."""
    try:
        for a in frag:
            at = mol.GetAtomWithIdx(a)
            nb = [n.GetIdx() for n in at.GetNeighbors()]
            if len(nb) != 3:                      # only 3-coordinate (sp2-like) centres; 4-coord sp3 = collapse job
                continue
            dev_orig = _pyramidalisation(P_orig[a], [P_orig[j] for j in nb])
            if dev_orig > planar_thresh:          # already non-planar (sp3 amine) in the original -> not our concern
                continue
            dev_cand = _pyramidalisation(P_cand[a], [P_cand[j] for j in nb])
            if dev_cand > dev_orig + worsen_margin:   # a PLANAR 3-coord centre got pyramidalised -> reject
                return True
        return False
    except Exception:
        return False


def _coordination_regressed(mol, syms, P_cand, metal_idxs) -> bool:
    """ROLLBACK GUARD (2026-07-25, COMPLETENESS floor): True if the candidate has a metal contact with an
    atom that is NOT a metal-neighbour in the mol TOPOLOGY -> a SPURIOUS new M-L bond that raises CN and
    alters the coordination ISOMER (RANFOE: the re-embed drifted a thiourea backbone N to 2.36 A from Cu ->
    CN3 trigonal-planar became CN4 tetrahedron; the crystal is CN3, so the CCDC isomer is lost).  The M-D
    rollback preserves the ORIGINAL donors but does not stop a NEW atom entering.  Checked against the SMILES
    graph (the coordination authority) with a GENEROUS metal-aware cutoff (metal covalent radius floored at
    1.3 A -- _COV has no metals -> defaults 0.95, too tight; the eye's donor cutoff is looser)."""
    try:
        for m in metal_idxs:
            mol_nbrs = set(n.GetIdx() for n in mol.GetAtomWithIdx(m).GetNeighbors())
            r_m = max(_cov(syms[m]), 1.3)
            Pm = P_cand[m]
            for i in range(len(syms)):
                if i == m or syms[i] == "H" or i in metal_idxs or i in mol_nbrs:
                    continue
                if float(np.linalg.norm(P_cand[i] - Pm)) < 1.4 * (_cov(syms[i]) + r_m):
                    return True
        return False
    except Exception:
        return False


def _local_defect_score(mol, syms, P, frag, frag_set, metal_idxs) -> float:
    """Self-contained per-fragment defect severity (NO CCDC, NO eye): collapsed sp3 centres + impossible
    heavy-heavy short contacts (d/covsum<0.65) + crushed metal contacts (d/covsum<0.6) involving the
    fragment.  ROLLBACK NET-IMPROVEMENT check (2026-07-25): the candidate's score must be STRICTLY LOWER
    than the original's, else the reseat is not a net per-frame win and is rejected (ILEDUB: the graft popped
    a collapse but nudged the holistic quality up on a no_valid system it cannot rescue).  Ends the
    quality-regression tail by construction -- the reseat only ever LOWERS a frame's defects."""
    score = 0.0
    for ai in frag:
        a = mol.GetAtomWithIdx(ai)
        if a.GetSymbol() in _METALS or a.GetIsAromatic() or a.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        nb = _sp3_nbrs(mol, ai, frag_set, syms, P)                # geometric-aligned (see _sp3_nbrs)
        if len(nb) < 4 or sum(1 for j in nb if syms[j] != "H") < 2:
            continue
        near = sorted(nb, key=lambda j: float(np.sum((P[j] - P[ai]) ** 2)))[:4]
        v = [P[j] for j in near]
        if abs(float(np.dot(np.cross(v[1] - v[0], v[2] - v[0]), v[3] - v[0]))) / 6.0 < _VOL_MIN:
            score += 2.0
    for i in frag:
        if syms[i] == "H":
            continue
        Pi = P[i]
        for j in range(len(syms)):
            if j == i or j in frag_set or syms[j] == "H":
                continue
            cs = _cov(syms[i]) + _cov(syms[j])
            if cs <= 0:
                continue
            r = float(np.linalg.norm(Pi - P[j])) / cs
            if syms[j] in _METALS or syms[i] in _METALS:
                if r < 0.6:
                    score += 1.5
            elif r < 0.65:
                score += 1.0
    return score


def _max_angle_dev(P, ai: int, near: List[int], ideal: float = 109.47) -> float:
    """Max deviation (deg) of the neighbour-CENTRE-neighbour angles at ai from the ideal (109.47 = sp3)."""
    c = P[ai]
    m = 0.0
    for i in range(len(near)):
        for j in range(i + 1, len(near)):
            v1 = P[near[i]] - c
            v2 = P[near[j]] - c
            dn = float(np.linalg.norm(v1) * np.linalg.norm(v2))
            if dn < 1e-8:
                continue
            ang = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2)) / dn))))
            m = max(m, abs(ang - ideal))
    return m


def _worsens_sp3(mol, frag, syms, P_orig, P_cand, margin: float = 10.0) -> bool:
    """ROLLBACK GUARD (2026-07-26, ANGLE-based + RELATIVE to MATCH the eye's methyl_broken / sp3 hyb-angle,
    which measures the H-C-H / X-C-X angle deviation, NOT tetra-volume).  True if a 4-coordinate sp3 centre
    is SIGNIFICANTLY MORE angle-distorted in the candidate than in the original (max angle-dev from 109.47
    rises by > margin).  Catches HERVOQ (C11 methyl 5->26 deg, C3 backbone 49->64 deg) the tetra-volume
    version missed.  RELATIVE (not an absolute floor) so a collapsed centre the reseat is FIXING (large dev
    in orig -> small in cand) is never blocked -- only a centre the reseat makes WORSE."""
    fs = set(frag)
    try:
        for a in frag:
            at = mol.GetAtomWithIdx(a)
            if at.GetSymbol() in _METALS or at.GetIsAromatic() or at.GetHybridization() != Chem.HybridizationType.SP3:
                continue
            nb = _sp3_nbrs(mol, a, fs, syms, P_orig)
            if len(nb) < 4:
                continue
            near = sorted(nb, key=lambda j: float(np.sum((P_orig[j] - P_orig[a]) ** 2)))[:4]
            if _max_angle_dev(P_cand, a, near) > _max_angle_dev(P_orig, a, near) + margin:
                return True                       # the reseat made this sp3 centre's angles WORSE -> reject
        return False
    except Exception:
        return False


def _worsens_bonds(mol, frag, syms, P_orig, P_cand, margin: float = 0.15) -> bool:
    """ROLLBACK GUARD (2026-07-26, MATCHES the eye's org_bond / M-L length axes): True if a bond within the
    fragment is SIGNIFICANTLY MORE off its covalent-sum length in the candidate than in the original (|d -
    covsum| rises by > margin A).  Relative -> a bond the reseat un-compresses (collapse fix) is never
    blocked; only a bond the reseat stretches/compresses is rejected."""
    fs = set(frag)
    try:
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if i not in fs or j not in fs or syms[i] == "H" or syms[j] == "H":
                continue
            cs = _cov(syms[i]) + _cov(syms[j])
            off_o = abs(float(np.linalg.norm(P_orig[i] - P_orig[j])) - cs)
            off_c = abs(float(np.linalg.norm(P_cand[i] - P_cand[j])) - cs)
            if off_c > off_o + margin:
                return True
        return False
    except Exception:
        return False


def _reseat_frame(mol, xyz: str) -> Optional[str]:
    """Re-seat every collapsed fragment in one frame from a clean isolated embed; None if nothing changed."""
    _tr = os.environ.get("DELFIN_TRACE_SEATING", "0") == "1"

    def _t(msg):
        if _tr:
            try:
                import sys as _s
                print("[SEATING] RESEAT " + msg, file=_s.stderr, flush=True)
            except Exception:
                pass

    # capture any leading xyz header (count + comment) so a re-seated frame keeps the exact input format
    _hdr: List[str] = []
    for _ln in xyz.split("\n"):
        _pp = _ln.split()
        _isatom = False
        if len(_pp) >= 4:
            try:
                float(_pp[1]); float(_pp[2]); float(_pp[3]); _isatom = True
            except (ValueError, IndexError):
                pass
        if _isatom:
            break
        _hdr.append(_ln)

    syms, P = _parse(xyz)
    if len(syms) != mol.GetNumAtoms():
        _t("skip: natoms frame=%d mol=%d (order guard)" % (len(syms), mol.GetNumAtoms()))
        return None                                               # atom-order assumption broken -> skip
    if any(syms[i] != mol.GetAtomWithIdx(i).GetSymbol() for i in range(len(syms))):
        _t("skip: element-sequence mismatch (frame order != mol order)")
        return None                                               # order not confirmed -> skip (safe)
    metal_idxs = [i for i in range(len(syms)) if syms[i] in _METALS]
    if not metal_idxs:
        _t("skip: no metal in frame (natoms=%d, first syms=%s)" % (len(syms), syms[:6]))
        return None
    changed = False
    P = P.copy()
    _frags = _mol_fragments(mol, len(syms))
    _t("frame: %d frags, %d collapsed, molSP3=%d" % (
        len(_frags), sum(1 for f in _frags if _fragment_collapsed(mol, f, syms, P)),
        sum(1 for a in mol.GetAtoms() if a.GetHybridization() == Chem.HybridizationType.SP3)))
    for frag in _frags:
        if not _fragment_collapsed(mol, frag, syms, P):
            continue
        donors = _metal_donors(mol, frag, syms, P, metal_idxs)
        if len(donors) < 1:
            _t("collapsed frag (%d atoms) but 0 donors" % len(frag))
            continue
        clean = _embed_fragment_clean(mol, frag, donors, np.array([P[d] for d in donors], float),
                                      orig_P=P, syms=syms)
        if clean is None or any(d not in clean for d in donors):
            _t("collapsed frag (%d atoms, %d donors): isolated embed FAILED" % (len(frag), len(donors)))
            continue
        _t("collapsed frag (%d atoms, %d donors): embed OK, aligning" % (len(frag), len(donors)))
        A = np.array([clean[d] for d in donors], float)
        B = np.array([P[d] for d in donors], float)
        frag_set = set(frag)
        base_clash = _clash_count(syms, P, frag_set)
        base_defect = _local_defect_score(mol, syms, P, frag, frag_set, metal_idxs)
        # candidate placements: Kabsch on donors; for a 2-donor chelate resolve the residual rotation
        # about the donor-donor axis by minimising clash with the rest of the complex.
        best_P = None
        best_key = None
        R0, t0 = _kabsch(A, B) if len(donors) >= 2 else (np.eye(3), B[0] - clean[donors[0]])
        base_place = {i: (R0 @ clean[i] + t0) for i in frag}
        if len(donors) == 2:
            axis = B[1] - B[0]
            an = float(np.linalg.norm(axis))
            u = axis / an if an > 1e-8 else np.array([0.0, 0.0, 1.0])
            piv = 0.5 * (B[0] + B[1])
            angles = [2.0 * math.pi * k / 18.0 for k in range(18)]
        else:
            u = None
            angles = [0.0]
        for th in angles:
            place = {}
            if u is not None and th != 0.0:
                c, s = math.cos(th), math.sin(th)
                for i in frag:
                    v = base_place[i] - piv
                    place[i] = (v * c + np.cross(u, v) * s + np.outer(v @ u, u) * (1.0 - c)) + piv
            else:
                place = dict(base_place)
            cand = P.copy()
            for i in frag:
                cand[i] = place[i]
            # rollback guards: M-D preserved, collapse resolved, clash not worse
            md_ok = all(float(np.linalg.norm(cand[d] - P[m])) < _MD_FACTOR * (_cov(syms[d]) + _cov(syms[m])) + _MD_TOL
                        for d in donors for m in metal_idxs
                        if float(np.linalg.norm(P[d] - P[m])) < _MD_FACTOR * (_cov(syms[d]) + _cov(syms[m])))
            if not md_ok:
                continue
            if _fragment_collapsed(mol, frag, syms, cand):
                continue
            if _worsens_planarity(mol, frag, syms, P, cand):
                continue                                          # reject: pyramidalises an sp2/planar centre
            if _worsens_sp3(mol, frag, syms, P, cand):
                continue                                          # reject: distorts a clean sp3 (methyl/backbone)
            if _worsens_bonds(mol, frag, syms, P, cand):
                continue                                          # reject: stretches/compresses a bond (org_bond/M-L)
            if _coordination_regressed(mol, syms, cand, metal_idxs):
                continue                                          # reject: spurious new M-L bond (CN/donor change)
            # DONOR-PRESERVATION (2026-07-26, self-contained proxy for the CRYSTAL-anchored poly_cshm axis
            # which cannot be rebuilt in DELFIN without CCDC): the coordination polyhedron is defined by the
            # donor positions -- if the reseat holds every donor within 0.15 A of its original spot, the
            # polyhedron (and its CShM to ANY reference incl the crystal) is preserved.  Rejects the residual
            # bite-mismatch donor drift that shifted EKAKIK's poly_cshm 0.07.  The reseat aligns donors to
            # the original anyway (Kabsch), so a clean bite-matched reseat passes; only a drift is rejected.
            if any(float(np.linalg.norm(cand[d] - P[d])) > 0.15 for d in donors):
                continue                                          # reject: a donor drifted -> poly shift
            if _local_defect_score(mol, syms, cand, frag, frag_set, metal_idxs) >= base_defect:
                continue                                          # reject: not a strict net per-frame improvement
            clash = _clash_count(syms, cand, frag_set)
            key = clash
            if best_key is None or key < best_key:
                best_key, best_P = key, cand
        if best_P is not None and best_key is not None and best_key <= base_clash:
            P = best_P
            changed = True
    if not changed:
        return None
    body = _write(syms, P)
    return ("\n".join(_hdr) + "\n" + body) if _hdr else body


def _frame_topology_valid(mol, syms, P) -> bool:
    """A frame is 'good enough to protect' if it has NO collapsed sp3 fragment AND every heavy-heavy mol
    bond is realised (not broken/compressed).  ROOT never-worse gate (2026-07-25): if ANY frame of a
    system is valid, the system already realises its CCDC geometry (isomer/biaryl/pyramid/coordination) in
    that frame, so the reseat must NEVER touch the system -- the fresh embed of a sibling collapsed frame
    disturbs that axis (RANFOE/SAQCOC/SUKKAL had CLEAN frames + a collapsed sibling and lost their CCDC
    axis when the sibling was reseated).  Only FULLY-collapsed systems (AQIBAE: all 32 frames collapsed,
    no valid frame) are rescued -- there is no realised CCDC axis to lose."""
    if P is None or len(syms) != mol.GetNumAtoms():
        return False
    try:
        for frag in _mol_fragments(mol, len(syms)):
            if _fragment_collapsed(mol, frag, syms, P):
                return False
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if i >= len(syms) or j >= len(syms) or syms[i] == "H" or syms[j] == "H":
                continue
            d = float(np.linalg.norm(P[i] - P[j]))
            cs = _cov(syms[i]) + _cov(syms[j])
            if d < 0.55 * cs or d > 1.4 * cs:                     # compressed or broken bond
                return False
        return True
    except Exception:
        return False


def correct_results(mol, results):
    """Re-seat collapsed fragments across all frames.  Per-frame: keep the re-seated frame only when it
    changed (a collapse was resolved under the rollback guards); otherwise the original is kept."""
    if not _HAVE_RDKIT or not results or mol is None:
        return results
    # DELFIN frames carry explicit H; make the topology mol explicit-H too so frame index == mol index
    # (verified per frame by the element-sequence guard in _reseat_frame -> a mismatch safely skips).
    # Conditional: if mol already has explicit H, AddHs would be a no-op or (worse) mis-add -> use as-is.
    try:
        has_h = any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
        molH = mol if has_h else Chem.AddHs(mol)
    except Exception:
        molH = mol
    # ROOT never-worse gate (2026-07-25): only rescue FULLY-collapsed systems (no valid frame).  If ANY
    # frame already realises the CCDC geometry (topology-valid), NEVER touch the system -- the fresh embed
    # of a collapsed sibling frame would disturb the axis that frame realises (RANFOE CN3->CN4, SAQCOC
    # biaryl, SUKKAL pyramid all had clean frames + a collapsed sibling).  AQIBAE = all frames collapsed
    # -> rescued.  Keeps the no_valid rescues (+valid), removes the CCDC-axis losses by construction.
    try:
        for item in results:
            _xyz = item[0] if isinstance(item, (tuple, list)) else item
            _s, _P = _parse(_xyz)
            if _frame_topology_valid(molH, _s, _P):
                return results                                    # system has a good frame -> skip entirely
    except Exception:
        pass
    out = []
    n_fixed = 0
    for item in results:
        # DELFIN pipeline frames are (xyz, label) tuples; the self-test passes bare xyz strings -> both.
        is_tuple = isinstance(item, (tuple, list))
        xyz = item[0] if is_tuple else item
        try:
            new = _reseat_frame(molH, xyz)
        except Exception:
            new = None
        if new is not None:
            out.append(((new,) + tuple(item[1:])) if is_tuple else new)
            n_fixed += 1
        else:
            out.append(item)
    if n_fixed and os.environ.get("DELFIN_TRACE_SEATING", "0") == "1":
        try:
            import sys as _sys
            print("[SEATING] ISOLATED_RESEAT re-seated %d/%d frames" % (n_fixed, len(results)),
                  file=_sys.stderr, flush=True)
        except Exception:
            pass
    return out


# ------------------------------------------------------------------ #
# self-test (full visibility, bypasses the loop/--debug plumbing):
#   python _isolated_reseat.py <archive.xyz> "<smiles>"
# ------------------------------------------------------------------ #
if __name__ == "__main__":
    import sys
    os.environ["DELFIN_TRACE_SEATING"] = "1"
    archive, smi = sys.argv[1], sys.argv[2]
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        try:
            Chem.SanitizeMol(mol, sanitizeOps=(Chem.SanitizeFlags.SANITIZE_ALL
                                               & ~Chem.SanitizeFlags.SANITIZE_KEKULIZE))
        except Exception:
            pass
    molH = Chem.AddHs(mol)
    print("mol heavy=%d  molH=%d" % (mol.GetNumAtoms(), molH.GetNumAtoms()))
    # parse archive -> atom-block-only frames (drop count + comment lines)
    lines = open(archive).read().splitlines()
    frames = []
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if s.isdigit():
            n = int(s)
            atoms = lines[i + 2:i + 2 + n]
            frames.append("\n".join(atoms))
            i += 2 + n
        else:
            i += 1
    print("frames parsed: %d" % len(frames))
    if frames:
        syms0, _P0 = _parse(frames[0])
        print("frame0 atoms=%d  molH atoms=%d" % (len(syms0), molH.GetNumAtoms()))
        if len(syms0) == molH.GetNumAtoms():
            mism = [i for i in range(len(syms0)) if syms0[i] != molH.GetAtomWithIdx(i).GetSymbol()]
            print("element-order MISMATCHES: %d  first: %s" % (len(mism), mism[:8]))
            print("  frame0 first 12 syms:", syms0[:12])
            print("  molH   first 12 syms:", [molH.GetAtomWithIdx(k).GetSymbol() for k in range(12)])
        else:
            print("ATOM-COUNT MISMATCH (frame != molH) -> heavy-map needed")
    out = correct_results(mol, frames)
    n_changed = sum(1 for a, b in zip(frames, out) if a != b)
    print("=> correct_results changed %d/%d frames" % (n_changed, len(frames)))
    if len(sys.argv) > 3:                                         # 3rd arg = output path -> write multiframe xyz
        with open(sys.argv[3], "w") as fh:
            for k, blk in enumerate(out):
                lns = [x for x in blk.split("\n") if x.strip()]
                fh.write("%d\nreseat frame%d\n%s\n" % (len(lns), k, "\n".join(lns)))
        print("wrote %s" % sys.argv[3])
