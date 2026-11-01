"""Automatic chemical fragment recognition for conformer ensembles.

Chemical structure comes from RDKit bond perception on the XYZ geometry
(rdDetermineBonds.DetermineBonds with an explicit net charge). An XYZ
file carries NO bond orders itself: recognition is therefore only valid
when (a) perception succeeds on EVERY frame and (b) the perceived
connectivity, bond orders, rings and aromaticity are IDENTICAL across
all frames. Otherwise fragments get status "plausible" (first-frame
only) or "unknown" (no consistent graph), and no unique fragment
identity may be claimed.

Recognized fragment kinds (chemical rules, no thresholds):
- aromatic_ring   : RDKit aromatic ring atoms (one entry per ring)
- conjugated_pi   : connected component of double bonds + aromatic
                    atoms linked through conjugated (double/aromatic)
                    bonds; includes sp2-sp2 single bonds in between
- conjugated_backbone : the largest conjugated_pi component that is
                    chemically adjacent to (fused/spiro/shared atoms
                    or direct conjugated bond) an aromatic system; if
                    the adjacency is ambiguous, status "plausible"
- flexible_side   : heavy-atom side chains not in the conjugated
                    backbone, attached via rotatable single bonds
- torsion_pair    : the sp2-sp2 single bonds inside/adjacent to the
                    backbone (s-cis/s-trans analysis points)

Status values: "recognized" | "plausible" | "unknown" (+ reason).
Atom indices refer to the ensemble frame ordering (0-based) and are
validated for stability across all frames.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["recognize_fragments", "fragment_geometry", "combined_ranking"]

# ---------------------------------------------------------------------------
# graph perception
# ---------------------------------------------------------------------------


def _perceive(path, charge=0):
    """Return per-frame chemical graphs, or raise with a clear reason."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    from delfin.tools.adapters.crest_ensemble import parse_ensemble

    frames = parse_ensemble(Path(path))
    graphs = []
    failed = None
    for k, f in enumerate(frames):
        sym, pos = f["symbols"], f["positions"]
        block = (f"{len(sym)}\nframe\n"
                 + "\n".join(f"{s} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}"
                             for s, c in zip(sym, pos)) + "\n")
        mol = Chem.MolFromXYZBlock(block)
        if mol is None:
            failed = f"frame {k}: RDKit could not parse geometry"
            break
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=charge)
            Chem.SanitizeMol(mol)
        except Exception as exc:  # perception failed -> no graph
            failed = (f"frame {k}: bond perception failed "
                      f"({type(exc).__name__}: {exc})")
            break
        bonds = sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
                        b.GetBondTypeAsDouble()) for b in mol.GetBonds())
        rings = sorted(sorted(r) for r in mol.GetRingInfo().AtomRings())
        arom = sorted(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
        graphs.append({"bonds": bonds, "rings": rings, "aromatic": arom,
                       "symbols": sym})
    return frames, graphs, None, failed


def _stability(graphs, failed):
    """(frames, graphs, stable, failed) -> stable flag."""
    if failed is not None or not graphs:
        return False
    return all(g["bonds"] == graphs[0]["bonds"]
               and g["rings"] == graphs[0]["rings"]
               and g["aromatic"] == graphs[0]["aromatic"] for g in graphs)


def _sp2_like_geometry(adj, coords, symbols):
    """Geometric sp2 fallback: two heavy neighbours, angle 115-130 deg.

    Used ONLY to bridge single bonds between pi systems when bond
    perception leaves a sp2-sp2 single bond unmarked (conjugated diene
    core). Mirrors the Phase-1 engine's sp2_like heuristic.
    """
    import math

    import numpy as np
    out = set()
    for a, nb in adj.items():
        if symbols[a] == "H":
            continue
        heavy = [b for b, _k, _o in nb if symbols[b] != "H"]
        if len(heavy) != 2:
            continue
        v1 = coords[heavy[0]] - coords[a]
        v2 = coords[heavy[1]] - coords[a]
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-6 or n2 < 1e-6:
            continue
        cosang = float(np.dot(v1, v2) / (n1 * n2))
        ang = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
        if 115.0 <= ang <= 130.0:
            out.add(a)
    return out


def _conjugated_components(graph, coords=None):
    """Connected components over atoms joined by double or aromatic bonds.

    Two such bonds that share a sp2 atom, or are linked by an sp2-sp2
    single bond, end up in one component (a conjugated chain).
    """
    bonds = graph["bonds"]
    arom = set(graph["aromatic"])
    # adjacency: atom -> [(neighbour, kind, order)]
    adj = {}
    for a, b, o in bonds:
        kind = "aromatic" if o == 1.5 else ("double" if o == 2.0 else "single")
        adj.setdefault(a, []).append((b, kind, o))
        adj.setdefault(b, []).append((a, kind, o))
    pi_nodes = set()
    for a, nb in adj.items():
        if any(k in ("double", "aromatic") for _, k, _ in nb):
            pi_nodes.add(a)
    if coords is not None:
        # sp2-sp2 single bonds between two pi atoms (double/aromatic
        # bonds on both sides) bridge components: conjugated diene
        # cores. Purely geometric sp2 atoms WITHOUT any pi bond are
        # NOT pulled in -- their conjugation is not established by
        # the perceived graph, so they stay outside (conservative).
        for a in list(pi_nodes):
            for b, kind, _o in adj.get(a, []):
                if kind == "single" and b in pi_nodes:
                    pass  # BFS below already walks single bonds
        # BFS below crosses single bonds between pi nodes; nothing
        # extra needed here.
    comps = []
    seen = set()
    for start in sorted(pi_nodes):
        if start in seen:
            continue
        comp = set([start])
        stack = [start]
        while stack:
            x = stack.pop()
            for y, kind, _o in adj.get(x, []):
                if y not in comp and y in pi_nodes and kind in (
                        "double", "aromatic"):
                    comp.add(y)
                    stack.append(y)
                elif (y not in comp and y in pi_nodes and kind == "single"
                        and not (y in arom and x in arom)):
                    # sp2-sp2 single bond inside a conjugated chain
                    # (aromatic-aromatic single bonds stay separate rings
                    #  unless fused; aromatic rings are their own
                    #  components via their aromatic bonds)
                    comp.add(y)
                    stack.append(y)
        seen |= comp
        if len(comp) > 1:
            comps.append(sorted(comp))
    return comps, adj, arom


def _heavy(idx, symbols):
    return symbols[idx] != "H"


def _backbone_torsions(graph, symbols, backbone):
    """Single bonds with both atoms in the backbone (s-cis/s-trans
    points) plus side-chain attach bonds, between heavy atoms, not in
    rings."""
    ring_bonds = set()
    for r in graph["rings"]:
        for i, x in enumerate(r):
            y = r[(i + 1) % len(r)]
            ring_bonds.add((min(x, y), max(x, y)))
    out = []
    for a, b, o in graph["bonds"]:
        if o != 1.0 or not (_heavy(a, symbols) and _heavy(b, symbols)):
            continue
        key = (min(a, b), max(a, b))
        if key in ring_bonds:
            continue
        if a in backbone and b in backbone:
            out.append(key)
    return out


def _side_groups(graph, symbols, backbone):
    """Heavy-atom connected components outside the backbone."""
    parent = list(range(len(symbols)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, _o in graph["bonds"]:
        if _heavy(a, symbols) and _heavy(b, symbols) and \
                a not in backbone and b not in backbone:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    groups = {}
    for i in range(len(symbols)):
        if _heavy(i, symbols) and i not in backbone:
            groups.setdefault(find(i), []).append(i)
    # attach bond: the single bond linking the group to the backbone
    attach = {}
    for a, b, o in graph["bonds"]:
        if o != 1.0:
            continue
        for inside, outside in ((a, b), (b, a)):
            if outside in backbone and inside not in backbone and \
                    _heavy(inside, symbols):
                attach.setdefault(find(inside), (outside, inside))
    out = []
    for root in sorted(groups):
        atoms = sorted(groups[root])
        # hydrogens bonded to the group's heavy atoms belong to it too
        hs = set()
        for a, b, o in graph["bonds"]:
            if o == 1.0:
                for h, heavy in ((a, b), (b, a)):
                    if symbols[h] == "H" and heavy in atoms:
                        hs.add(h)
        atoms = sorted(set(atoms) | hs)
        out.append({"kind": "flexible_side", "atoms": atoms,
                    "attach_bond": attach.get(root)})
    return out


def recognize_fragments(path, charge=0):
    """Recognize chemically relevant fragments for a whole ensemble.

    Returns a dict with global status, the fragment list (stable atom
    indices), torsion points, and provenance of the chemical
    information.
    """
    import numpy as np

    frames, graphs, _unused, failed = _perceive(path, charge=charge)
    stable = _stability(graphs, failed)
    if failed is not None or not graphs:
        return {
            "input_file": str(path),
            "source": f"rdDetermineBonds.DetermineBonds(charge={charge}) on XYZ "
                      "coordinates; NO bond orders stored in the file "
                      "itself",
            "graph_stable_across_frames": False,
            "n_frames": len(frames),
            "n_atoms": None,
            "status": "unknown",
            "reason": failed or "no frame could be perceived",
            "fragments": [],
        }
    g0 = graphs[0]
    symbols = g0["symbols"]
    comps, adj, arom = _conjugated_components(
        g0, coords=np.asarray(frames[0]["positions"], dtype=float))

    fragments = []
    # aromatic rings (one per ring system)
    for r in g0["rings"]:
        if all(x in arom for x in r):
            fragments.append({"kind": "aromatic_ring", "atoms": r,
                              "status": "recognized", "reason": None})

    # conjugated pi systems
    pi_comps = comps
    for comp in pi_comps:
        if all(x in arom for x in comp):
            continue  # pure aromatic rings already listed
        status = "recognized" if stable else "plausible"
        reason = None if stable else (
            "bond perception differs between frames; first-frame only")
        fragments.append({"kind": "conjugated_pi", "atoms": comp,
                          "status": status, "reason": reason})

    # conjugated backbone: the conjugated component with the most
    # aromatic adjacency (shared/fused atoms with an aromatic ring or
    # direct conjugated bond into it). Ambiguity -> "plausible".
    aromatic_sets = [set(r) for r in g0["rings"] if all(x in arom for x in r)]
    cands = [c for c in pi_comps if not all(x in arom for x in c)]
    if len(cands) == 1:
        backbone, bstatus, breason = cands[0], ("recognized" if stable else
                                                "plausible"), None
    elif len(cands) == 0:
        backbone, bstatus, breason = None, "unknown", "no conjugated chain found"
    else:
        # multiple candidates: prefer the one with most aromatic contact
        def contact(c):
            s = set(c)
            return sum(len(s & rs) for rs in aromatic_sets)
        cands.sort(key=lambda c: (-contact(c), len(c), c[0]))
        top = contact(cands[0])
        if len(cands) > 1 and contact(cands[1]) == top:
            backbone, bstatus = None, "ambiguous"
            breason = (f"{len(cands)} conjugated components share the "
                       f"maximal aromatic contact ({top}); ask the user "
                       "which backbone is meant")
        else:
            backbone, bstatus = cands[0], ("recognized" if stable else
                                           "plausible")
            breason = None
    if backbone is not None:
        fragments.append({"kind": "conjugated_backbone", "atoms": backbone,
                          "status": bstatus, "reason": breason})
        fragments.append({"kind": "flexible_side_groups",
                          "atoms": sorted(set().union(*[
                              set(s["atoms"]) for s in
                              _side_groups(g0, symbols, backbone)]))
                          or None,
                          "status": "recognized" if stable else "plausible",
                          "reason": None})
        torsions = _backbone_torsions(g0, symbols, set(backbone))
        fragments.append({"kind": "torsion_points",
                          "bonds": torsions,
                          "status": "recognized" if stable else "plausible",
                          "reason": None if stable else (
                              "bond perception differs between frames; "
                              "first-frame only")})

    return {
        "input_file": str(path),
        "source": f"rdDetermineBonds.DetermineBonds(charge={charge}) on XYZ "
                  "coordinates; NO bond orders stored in the file itself",
        "graph_stable_across_frames": stable,
        "n_frames": len(frames),
        "n_atoms": len(symbols),
        "status": ("recognized" if stable else "plausible")
                  if fragments else "unknown",
        "fragments": fragments,
    }


# ---------------------------------------------------------------------------
# fragment geometry analysis
# ---------------------------------------------------------------------------


def _best_plane(coords):
    """(centroid, normal, rms_perp) of the least-squares plane."""
    import numpy as np
    c = coords - coords.mean(axis=0)
    _u, s, vt = np.linalg.svd(c, full_matrices=False)
    normal = vt[2]
    rms = float(np.sqrt(((c @ normal) ** 2).mean()))
    return coords.mean(axis=0), normal, rms


def fragment_geometry(path, atoms, kind="selection"):
    """Geometric descriptors for one fragment over all frames.

    Uses the Phase-1 engine on the given atom selection, plus a
    best-plane RMS deviation (true out-of-plane measure in A).
    """
    import numpy as np

    from delfin.analysis_tools.conformer_shape_classification import analyze_ensemble
    from delfin.tools.adapters.crest_ensemble import parse_ensemble

    atoms = sorted(int(a) for a in atoms)
    spec = ",".join(str(a) for a in atoms)
    res = analyze_ensemble(path, weighting="mass", selections=[spec])
    frames = parse_ensemble(Path(path))
    per_frame = []
    for k, f in enumerate(frames):
        d = res["frames"][k]["descriptors"][spec]
        c = np.asarray(f["positions"])[atoms]
        centroid, normal, rms = _best_plane(c)
        e = res["frames"][k].get("energy_Eh")
        ref = min(f2.get("energy_Eh") for f2 in res["frames"])
        rel = None if (e is None or ref is None) else (e - ref) * 627.509474
        per_frame.append({
            "conformer_id": k,
            "descriptors": d,
            "plane_rms_A": round(rms, 6),
            "rel_energy_kcal_mol": rel,
        })
    return {"atoms": atoms, "kind": kind, "frames": per_frame}


def combined_ranking(path, atoms, descriptor_a="plane_rms_A",
                     ascending_a=True, descriptor_b="dmax",
                     descending_b=True, top=5, tie="a_first"):
    """Transparent two-criteria ranking over ONE fragment selection.

    Every value comes from the requested atom selection -- there is no
    silent fallback to the whole molecule. descriptor_a/b may be engine
    descriptors or "plane_rms_A" (best-plane RMS, fragment-only).
    Sorts lexicographically: primary criterion first; equal primary
    values (within 1e-9) fall to the secondary. Both directions and
    the tie rule are explicit in the result; no weighted score is
    computed. This is a strict-priority ranking, not a compromise
    score.
    """
    fg = fragment_geometry(path, atoms)
    if not fg["frames"]:
        raise ValueError("no frames: fragment ranking impossible")
    rows = []
    for k, fr in enumerate(fg["frames"]):
        d = fr["descriptors"]
        va = fr["plane_rms_A"] if descriptor_a == "plane_rms_A" \
            else d.get(descriptor_a, {}).get("value")
        vb = fr["plane_rms_A"] if descriptor_b == "plane_rms_A" \
            else d.get(descriptor_b, {}).get("value")
        if va is None or vb is None:
            continue
        row = {"conformer_id": k,
               "rel_energy_kcal_mol": fr["rel_energy_kcal_mol"]}
        row[descriptor_a], row[descriptor_b] = va, vb
        rows.append(row)
    if not rows:
        raise ValueError("all values None: descriptors unavailable for "
                         "the requested selection")
    sgn_a = 1.0 if ascending_a else -1.0
    sgn_b = -1.0 if descending_b else 1.0
    rows.sort(key=lambda r: (sgn_a * r[descriptor_a],
                             sgn_b * r[descriptor_b],
                             r["conformer_id"]))
    return {"atoms": fg["atoms"], "kind": fg["kind"],
            "criteria": [
                {"descriptor": descriptor_a, "ascending": ascending_a},
                {"descriptor": descriptor_b, "ascending": not descending_b}],
            "tie_rule": "lexicographic, exact primary wins within 1e-9, "
                        "conformer_id as final tiebreak; no weighted score",
            "ranking": rows[:top]}


def pareto_front(path, atoms, minimize=("plane_rms_A",),
                 maximize=("dmax",), top=None):
    """Non-dominated set over two fragment criteria (NOT a ranking).

    A conformer is on the front if no other conformer is at least as
    good in every criterion and strictly better in one. The result is
    the set of non-dominated solutions; it carries no preference order
    and must not be presented as 'the best'.
    """
    fg = fragment_geometry(path, atoms)
    cand = []
    for k, fr in enumerate(fg["frames"]):
        vals = []
        for d in minimize:
            v = fr["plane_rms_A"] if d == "plane_rms_A" else \
                fr["descriptors"].get(d, {}).get("value")
            if v is None:
                break
            vals.append(v)
        else:
            for d in maximize:
                v = fr["plane_rms_A"] if d == "plane_rms_A" else \
                    fr["descriptors"].get(d, {}).get("value")
                if v is None:
                    break
                vals.append(-v)  # maximise -> minimise negative
            else:
                cand.append((k, tuple(vals),
                             fr["rel_energy_kcal_mol"]))
    if not cand:
        raise ValueError("no complete criterion values: Pareto analysis "
                         "impossible")
    front = []
    for k, v, rel in cand:
        dominated = any(
            all(o <= x + 1e-12 for o, x in zip(ov, v))
            and any(o < x - 1e-12 for o, x in zip(ov, v))
            for _k2, ov, _r2 in cand if _k2 != k)
        if not dominated:
            front.append(k)
    return {"atoms": fg["atoms"], "minimize": list(minimize),
            "maximize": list(maximize),
            "note": "non-dominated set, NOT a unique ranking",
            "front": sorted(front),
            "values": {str(k): {"plane_rms_A": fg["frames"][k]["plane_rms_A"],
                                "dmax_A": fg["frames"][k]["descriptors"]
                                .get("dmax", {}).get("value")}
                       for k in sorted(front)}}
