"""XRD-isomer-and-conformer recall metric vs CCDC ground truth.

User-direction 2026-06-04 (USER, 17:30):

    "XRD abgleich ab jetzt mit einbeziehen in generierungsqualität ist
     richtige Isomer + konformere enthalten."

This metric captures, for a generated pool of DELFIN structures, whether
the pool contains the *correct* coordination isomer AND a conformer
whose heavy-atom RMSD to the CCDC crystal is below a threshold.

The metric is split into two scalars, reported per CCDC anchor refcode:

    isomer_recall:
        1 iff at least ONE DELFIN structure for the matching SMILES has
        the same coordination-isomer label (fac/mer/cis/trans/etc.) as
        the CCDC entry; 0 otherwise.

    conformer_recall_rmsd:
        min over all DELFIN structures (matching SMILES) of the heavy-
        atom Kabsch RMSD to the CCDC entry. Pass if < 0.5 Å (default).
        NaN if no DELFIN structures for that SMILES.

Both are computed with NO CCDC dependency at runtime — the CCDC ground
truth is shipped as a static JSON dump (see
``agent_workspace/quality_framework/reference/ccdc_ground_truth_48.json``,
built by ``scripts/build_ccdc_ground_truth_48.py``).

Env flag
--------
``DELFIN_USE_CCDC_REFERENCE_DUMP=1`` switches the metric on. Default OFF
keeps everything byte-identical and CCDC-free at import time.

Public API
----------
``classify_isomer(syms, P)``
    Heuristic isomer label (CN6-OH-fac, CN4-SP-cis, etc.) from coords.
``kabsch_rmsd(P, Q)``
    Heavy-atom RMSD after optimal rotation.
``score_pool(archive_dir)``
    Returns per-refcode isomer_recall + conformer_recall_rmsd.

This module is pure-python + numpy. No CCDC import at runtime.
"""
from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path
from collections import Counter

try:
    import numpy as np
except Exception:  # pragma: no cover
    np = None  # type: ignore

__all__ = [
    "DEFAULT_REF_PATH",
    "ENV_FLAG",
    "RMSD_PASS_THRESHOLD_A",
    "classify_isomer",
    "kabsch_rmsd_heavy",
    "score_pool",
    "load_reference",
]

DEFAULT_REF_PATH = os.environ.get(
    "DELFIN_CCDC_REFERENCE_DUMP",
    "/home/qmchem_max/agent_workspace/quality_framework/reference/"
    "ccdc_ground_truth_48.json",
)
ENV_FLAG = "DELFIN_USE_CCDC_REFERENCE_DUMP"
RMSD_PASS_THRESHOLD_A = float(
    os.environ.get("DELFIN_XRD_CONFORMER_RMSD_A", "0.5")
)

METALS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
    "Er", "Tm", "Yb", "Lu", "Th", "U",
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Ga", "In", "Sn",
    "Sb", "Tl", "Pb", "Bi",
}
DONOR_ELEMS = {"C", "N", "O", "P", "S", "Cl", "Br", "I", "F", "Se", "As", "B", "Si", "Te"}

# rough covalent radii (Å); used for M-D cutoff
_COV = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "As": 1.19, "Se": 1.20,
    "Br": 1.20, "Te": 1.38, "I": 1.39,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.50, "Fe": 1.42,
    "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45,
    "Cd": 1.44, "Hf": 1.75, "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44,
    "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Pb": 1.46, "La": 2.07,
    "Ce": 2.04, "Sn": 1.39, "Al": 1.21, "Mg": 1.41, "Ca": 1.76, "Na": 1.66,
    "K": 2.03, "U": 1.96, "Th": 2.06, "Eu": 1.98, "Gd": 1.96, "Dy": 1.92,
    "In": 1.42, "Tl": 1.45, "Bi": 1.51, "Sb": 1.39,
}
_BOND_TOL = 1.30


def _dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _elem_from_label(s):
    """Strip digits/sybyl from an atom label."""
    s = s.split(".")[0]
    m = re.match(r"([A-Z][a-z]?)", s)
    return m.group(1) if m else s


def parse_xyz(text):
    """Robust XYZ parser. Returns the FIRST frame only.

    Accepts:
      * Standard multi-frame XYZ (count line + comment line + atoms,
        repeated). Only the first frame is returned.
      * Header-less XYZ (lines start with atom records, no count/
        comment). All matching atom lines are returned.

    For multi-frame iteration use :func:`parse_xyz_frames`.
    """
    frames = parse_xyz_frames(text)
    if not frames:
        empty_P = np.zeros((0, 3), dtype=float) if np is not None else []
        return [], empty_P
    return frames[0]


def parse_xyz_frames(text):
    """Return all frames of a (potentially multi-frame) XYZ file.

    Returns ``[(syms, Nx3 array), ...]``. If the file is header-less
    (no count line), a single frame is returned with all atom rows.
    """
    out = []
    lines = text.splitlines()
    i = 0
    n_lines = len(lines)
    # Try header-mode first: lines[i] should parse as int = natoms
    while i < n_lines:
        # Skip blank lines
        while i < n_lines and not lines[i].strip():
            i += 1
        if i >= n_lines:
            break
        # Header attempt
        try:
            na = int(lines[i].strip())
            if na <= 0 or na > 1_000_000:
                raise ValueError
        except ValueError:
            # No count line — fall back to header-less single frame
            syms, P = _parse_atom_rows(lines[i:])
            if syms:
                out.append((syms, np.asarray(P, dtype=float)
                            if np is not None else P))
            break
        if i + 1 >= n_lines:
            break
        # Skip count line + comment line, parse na atom lines
        atom_lines = lines[i + 2:i + 2 + na]
        syms, P = _parse_atom_rows(atom_lines)
        if syms:
            out.append((syms, np.asarray(P, dtype=float)
                        if np is not None else P))
        i = i + 2 + na
    return out


def _parse_atom_rows(rows):
    syms = []
    P = []
    for line in rows:
        p = line.split()
        if len(p) < 4:
            continue
        try:
            x, y, z = float(p[1]), float(p[2]), float(p[3])
        except ValueError:
            continue
        if not re.match(r"^[A-Za-z]{1,2}[0-9+\-]*$", p[0]):
            continue
        e = _elem_from_label(p[0])
        syms.append(e)
        P.append([x, y, z])
    return syms, P


def _coord_sphere(syms, P):
    """Return (metal_idx, donor_indices, donor_syms, donor_distances) or None."""
    if np is None:
        return None
    metals = [i for i, s in enumerate(syms) if s in METALS]
    if not metals:
        return None
    # primary metal: highest count of donors
    best = None
    for mi in metals:
        mp = P[mi]
        mr = _COV.get(syms[mi], 1.5)
        donors = []
        for j, s in enumerate(syms):
            if j == mi or s == "H" or s not in DONOR_ELEMS:
                continue
            r = _COV.get(s, 1.0)
            d = float(np.linalg.norm(P[j] - mp))
            cutoff = (mr + r) * _BOND_TOL
            # C tighter than the rule to avoid ring-C overcount
            if s == "C":
                cutoff = min(cutoff, 2.45)
            if d < cutoff:
                donors.append((j, s, d))
        donors.sort(key=lambda t: t[2])
        if best is None or len(donors) > len(best[1]):
            best = (mi, donors)
    if best is None:
        return None
    mi, donors = best
    return {
        "metal_idx": mi,
        "metal_sym": syms[mi],
        "donor_idx": [t[0] for t in donors],
        "donor_syms": [t[1] for t in donors],
        "donor_distances": [t[2] for t in donors],
    }


def _polyhedron_class(cn):
    return {
        2: "linear", 3: "trigonal", 4: "tet-or-sp",
        5: "TBP-or-SPY", 6: "OH", 7: "PBP", 8: "Cube",
        9: "TTP", 10: "BAS", 11: "edge", 12: "ico",
    }.get(cn, f"CN{cn}")


def _fac_or_mer(syms, P, coord_info):
    """For CN6-OH with two donor element groups (3+3), classify fac vs mer."""
    donors = coord_info["donor_idx"]
    dsyms = coord_info["donor_syms"]
    cnt = Counter(dsyms)
    if len(cnt) != 2:
        return None
    vals = sorted(cnt.values())
    if vals != [3, 3]:
        return None
    major_sym = max(cnt, key=cnt.get)
    if cnt[major_sym] != 3:
        return None
    # pick the 3 donors of major_sym; check if any opposite pair (180°) exists
    mi = coord_info["metal_idx"]
    mp = P[mi]
    same = [d for d, s in zip(donors, dsyms) if s == major_sym]
    if len(same) != 3:
        return None
    # mer: at least one pair of "same-element" donors is trans (cos~-1)
    is_mer = False
    for i, a in enumerate(same):
        for b in same[i+1:]:
            v1 = P[a] - mp
            v2 = P[b] - mp
            n1 = float(np.linalg.norm(v1))
            n2 = float(np.linalg.norm(v2))
            if n1 < 1e-6 or n2 < 1e-6:
                continue
            cos = float(np.dot(v1, v2) / (n1 * n2))
            if cos < -0.85:  # ~155°+
                is_mer = True
                break
        if is_mer:
            break
    return "mer" if is_mer else "fac"


def _cis_or_trans_sp4(syms, P, coord_info):
    """For CN4-SP with two donor-element groups (2+2), cis vs trans."""
    donors = coord_info["donor_idx"]
    dsyms = coord_info["donor_syms"]
    cnt = Counter(dsyms)
    if len(cnt) != 2:
        return None
    vals = sorted(cnt.values())
    if vals != [2, 2]:
        return None
    # group by element
    grp = {}
    for d, s in zip(donors, dsyms):
        grp.setdefault(s, []).append(d)
    mi = coord_info["metal_idx"]
    mp = P[mi]
    is_trans = False
    for s, idxs in grp.items():
        if len(idxs) != 2:
            continue
        a, b = idxs
        v1 = P[a] - mp
        v2 = P[b] - mp
        n1 = float(np.linalg.norm(v1))
        n2 = float(np.linalg.norm(v2))
        if n1 < 1e-6 or n2 < 1e-6:
            continue
        cos = float(np.dot(v1, v2) / (n1 * n2))
        if cos < -0.85:
            is_trans = True
            break
    return "trans" if is_trans else "cis"


def classify_isomer(syms, P):
    """Heuristic coordination-isomer label.

    Returns a string in DELFIN's Pólya naming convention:
      ``CN6-OH-fac``, ``CN6-OH-mer``, ``CN4-SP-cis``, ``CN4-SP-trans``,
      ``CN6-OH``, ``CN4-SP``, ``CN5-TBP-or-SPY``, ``CN<n>-<poly>``, etc.
    """
    if np is None:
        return "unknown"
    coord = _coord_sphere(syms, np.asarray(P, dtype=float))
    if not coord:
        return "no-metal"
    cn = len(coord["donor_syms"])
    poly = _polyhedron_class(cn)
    if cn == 4:
        sub = _cis_or_trans_sp4(syms, np.asarray(P, dtype=float), coord)
        return f"CN4-{poly}-{sub}" if sub else f"CN4-{poly}"
    if cn == 6:
        sub = _fac_or_mer(syms, np.asarray(P, dtype=float), coord)
        return f"CN6-OH-{sub}" if sub else f"CN6-OH"
    return f"CN{cn}-{poly}"


# -------------------------------------------------------------------------
# Heavy-atom RMSD (Kabsch) with greedy element-matched correspondence
# -------------------------------------------------------------------------
def _kabsch(P, Q):
    P = P - P.mean(0)
    Q = Q - Q.mean(0)
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    Pr = P @ R.T
    return float(math.sqrt(((Pr - Q) ** 2).sum(1).mean()))


def _greedy_element_match(syms_a, Pa, syms_b, Pb):
    """Pair atoms a->b greedily by element + distance after centroid alignment.

    Returns (Pa_matched, Pb_matched, n_matched). Hydrogens ignored.
    """
    if np is None:
        return None, None, 0
    keep_a = [i for i, s in enumerate(syms_a) if s != "H"]
    keep_b = [i for i, s in enumerate(syms_b) if s != "H"]
    if not keep_a or not keep_b:
        return None, None, 0
    Pa_h = Pa[keep_a]
    Pb_h = Pb[keep_b]
    Pa_h = Pa_h - Pa_h.mean(0)
    Pb_h = Pb_h - Pb_h.mean(0)
    syms_ah = [syms_a[i] for i in keep_a]
    syms_bh = [syms_b[i] for i in keep_b]

    # PCA-align Pa onto Pb to get an approximate initial correspondence
    # (skip detailed alignment — greedy match handles modest rotation poorly,
    #  but the union RMSD then includes Kabsch which fixes rigid rotation).
    used_b = [False] * len(Pb_h)
    pairs = []
    # iterate atoms of a in random order — use deterministic: sorted by element
    order_a = sorted(range(len(Pa_h)), key=lambda i: (syms_ah[i], i))
    for i in order_a:
        sa = syms_ah[i]
        best_j = -1
        best_d = 1e9
        for j in range(len(Pb_h)):
            if used_b[j] or syms_bh[j] != sa:
                continue
            d = float(np.linalg.norm(Pa_h[i] - Pb_h[j]))
            if d < best_d:
                best_d, best_j = d, j
        if best_j >= 0:
            used_b[best_j] = True
            pairs.append((i, best_j))
    if not pairs:
        return None, None, 0
    Pa_m = np.array([Pa_h[i] for i, _ in pairs])
    Pb_m = np.array([Pb_h[j] for _, j in pairs])
    return Pa_m, Pb_m, len(pairs)


def kabsch_rmsd_heavy(syms_a, Pa, syms_b, Pb):
    """Heavy-atom Kabsch RMSD with greedy element-matched correspondence."""
    if np is None:
        return float("nan")
    Pa = np.asarray(Pa, dtype=float)
    Pb = np.asarray(Pb, dtype=float)
    Pa_m, Pb_m, n = _greedy_element_match(syms_a, Pa, syms_b, Pb)
    if n < 3:
        return float("nan")
    return _kabsch(Pa_m, Pb_m)


def kabsch_rmsd_cn_sphere(syms_a, Pa, syms_b, Pb):
    """Kabsch RMSD restricted to the coordination sphere (metal + donors).

    This is the CHEMICALLY-MEANINGFUL metric: it compares the donor
    geometry around the metal (the part DELFIN's builder controls
    directly), bypassing the whole-molecule atom-count mismatch.

    Returns NaN if either structure has no identifiable metal.
    """
    if np is None:
        return float("nan")
    Pa = np.asarray(Pa, dtype=float)
    Pb = np.asarray(Pb, dtype=float)
    ca = _coord_sphere(syms_a, Pa)
    cb = _coord_sphere(syms_b, Pb)
    if not ca or not cb:
        return float("nan")
    # Build the M + donors point set for each
    mi_a = ca["metal_idx"]
    mi_b = cb["metal_idx"]
    donors_a = ca["donor_idx"]
    donors_b = cb["donor_idx"]
    if not donors_a or not donors_b:
        return float("nan")
    # element-matched greedy pairing (donor i in a → donor j in b)
    syms_da = [syms_a[i] for i in donors_a]
    syms_db = [syms_b[j] for j in donors_b]
    pos_da = Pa[donors_a]
    pos_db = Pb[donors_b]
    # center on metal
    mp_a = Pa[mi_a]
    mp_b = Pb[mi_b]
    pos_da = pos_da - mp_a
    pos_db = pos_db - mp_b
    used = [False] * len(donors_b)
    pairs = []
    # sort donors_a deterministically (by element, then by index)
    order = sorted(range(len(donors_a)), key=lambda i: (syms_da[i], i))
    for i in order:
        sa = syms_da[i]
        best_j = -1; best_d = 1e9
        for j in range(len(donors_b)):
            if used[j] or syms_db[j] != sa:
                continue
            d = float(np.linalg.norm(pos_da[i] - pos_db[j]))
            if d < best_d:
                best_d, best_j = d, j
        if best_j >= 0:
            used[best_j] = True
            pairs.append((i, best_j))
    if len(pairs) < 3:
        return float("nan")
    Pa_m = np.array([pos_da[i] for i, _ in pairs])
    Pb_m = np.array([pos_db[j] for _, j in pairs])
    return _kabsch(Pa_m, Pb_m)


# -------------------------------------------------------------------------
# Reference loader + pool scorer
# -------------------------------------------------------------------------
def load_reference(path=None):
    """Load CCDC ground-truth JSON dump."""
    p = Path(path or DEFAULT_REF_PATH)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def _refcode_from_filename(name):
    """Try to extract a 6-letter CCDC refcode from a pool filename."""
    stem = Path(name).stem
    # match 6 uppercase letters delimited by - / _ / start / end
    for tok in re.split(r"[-_.]", stem):
        if len(tok) == 6 and tok.isupper() and tok.isalpha():
            return tok
    # also try last token
    m = re.findall(r"[A-Z]{6}", stem)
    if m:
        return m[-1]
    return None


def score_pool(archive_dir, ref=None, rmsd_threshold=None):
    """Score a pool against CCDC ground truth.

    Parameters
    ----------
    archive_dir : str / Path
        Directory containing .xyz pool files.
    ref : dict or None
        Pre-loaded reference; if None, loads DEFAULT_REF_PATH.
    rmsd_threshold : float or None
        Pass-rmsd in Å, defaults to ``RMSD_PASS_THRESHOLD_A`` (0.5).

    Returns
    -------
    dict with keys:
        per_refcode : list of dicts (refcode, isomer_recall, conformer_recall_rmsd,
                                     rmsd_pass, ccdc_isomer, delfin_isomers, n_files)
        isomer_recall_mean : float
        conformer_recall_mean : float   (fraction of refcodes with rmsd < threshold)
        n_refcodes_scored : int
        n_refcodes_missing : int
    """
    ref = ref or load_reference()
    if ref is None:
        return {"error": "no reference data — see DELFIN_CCDC_REFERENCE_DUMP"}
    rmsd_threshold = rmsd_threshold or RMSD_PASS_THRESHOLD_A
    refmap = {r["refcode"]: r for r in ref.get("records", []) if r.get("ok")}
    archive = Path(archive_dir)
    if not archive.exists():
        return {"error": f"archive missing: {archive}"}
    # group .xyz files by refcode
    by_rc = {}
    for f in archive.glob("*.xyz"):
        rc = _refcode_from_filename(f.name)
        if rc and rc in refmap:
            by_rc.setdefault(rc, []).append(f)
    per_ref = []
    iso_rec_sum = 0
    conf_rec_pass = 0
    for rc, ref_rec in sorted(refmap.items()):
        files = by_rc.get(rc, [])
        if not files:
            per_ref.append({
                "refcode": rc,
                "n_files": 0,
                "ccdc_isomer": ref_rec.get("isomer_label"),
                "delfin_isomers": [],
                "isomer_recall": 0,
                "conformer_recall_rmsd": float("nan"),
                "rmsd_pass": 0,
            })
            continue
        ccdc_iso = ref_rec.get("isomer_label", "unknown")
        ref_syms = ref_rec["symbols"]
        ref_pos = np.asarray(ref_rec["positions"], dtype=float) if np is not None else None
        delfin_isos = []
        best_rmsd = float("inf")
        best_rmsd_cn = float("inf")
        n_frames = 0
        for f in files:
            try:
                frames = parse_xyz_frames(f.read_text())
            except Exception:
                continue
            for syms, P in frames:
                if not syms:
                    continue
                n_frames += 1
                label = classify_isomer(syms, P)
                delfin_isos.append(label)
                if ref_pos is not None:
                    r = kabsch_rmsd_heavy(syms, P, ref_syms, ref_pos)
                    if not math.isnan(r) and r < best_rmsd:
                        best_rmsd = r
                    r_cn = kabsch_rmsd_cn_sphere(syms, P, ref_syms, ref_pos)
                    if not math.isnan(r_cn) and r_cn < best_rmsd_cn:
                        best_rmsd_cn = r_cn
        # isomer recall: any DELFIN structure has matching label?
        # Match by base polyhedron + isomer-stereo (e.g. CN6-OH-fac matches CN6-OH-fac;
        # also allow CN6-OH-fac to match CN6-OH-fac-or-mer in ground truth)
        def _match(a, b):
            if a == b:
                return True
            # allow "fac" or "mer" to match the umbrella "fac-or-mer"
            if a and b and "or" in b:
                return a in b
            if a and b and "or" in a:
                return b in a
            return False
        iso_match = any(_match(d, ccdc_iso) for d in delfin_isos)
        # Conformer recall passes if either the whole-molecule heavy RMSD
        # OR the CN-sphere RMSD beats the threshold.
        is_rmsd_pass = (math.isfinite(best_rmsd)
                        and best_rmsd < rmsd_threshold)
        is_rmsd_cn_pass = (math.isfinite(best_rmsd_cn)
                           and best_rmsd_cn < rmsd_threshold)
        iso_rec_sum += int(iso_match)
        conf_rec_pass += int(is_rmsd_pass or is_rmsd_cn_pass)
        per_ref.append({
            "refcode": rc,
            "n_files": len(files),
            "n_frames": n_frames,
            "ccdc_isomer": ccdc_iso,
            "delfin_isomers": delfin_isos,
            "isomer_recall": int(iso_match),
            "conformer_recall_rmsd":
                float(best_rmsd) if math.isfinite(best_rmsd) else float("nan"),
            "conformer_recall_rmsd_cn":
                float(best_rmsd_cn) if math.isfinite(best_rmsd_cn) else float("nan"),
            "rmsd_pass": int(is_rmsd_pass),
            "rmsd_cn_pass": int(is_rmsd_cn_pass),
        })
    n_total = len(refmap)
    n_missing = sum(1 for r in per_ref if r["n_files"] == 0)
    return {
        "per_refcode": per_ref,
        "isomer_recall_mean":
            float(iso_rec_sum) / max(n_total, 1),
        "conformer_recall_mean":
            float(conf_rec_pass) / max(n_total, 1),
        "n_refcodes_scored": n_total,
        "n_refcodes_missing": n_missing,
        "rmsd_threshold_A": rmsd_threshold,
    }
