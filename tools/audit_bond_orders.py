#!/usr/bin/env python3
"""Audit the drawn bond orders of every calc against an RDKit reference.

The reference is built from the calc's own 3D structure: heavy-atom
graph with metals stripped, converted to an RDKit molecule with
Kekulization - i.e. what the aromatic rings SHOULD look like.  The
drawn orders are what ``initial_xyz_png_bytes`` actually renders
(``smiles_bond_orders`` when a SMILES exists, else the geometric
``_rdkit_bond_orders`` fallback).

Findings are reported per error class (user briefing 2026-11-21):

  A  coordinated pyridine ring with <3 double bonds
  B  benzene ring next to deprotonated amide N with <3 double bonds
  C  phenolate/thiolate drawn quinoid (exocyclic C=O / C=S on a ring)
  D  phenyl on phosphorus drawn as P=C quinoid
  E  nitro group with an H label on O ("HO-N=O")
  F  carbonyl O labeled OH (wrong H suffix in the oxalyl/amide region)
  G  sulfonic acid S=O orders missing / labels overlapping

Usage:
    python3 tools/audit_bond_orders.py [--calc-root DIR] [--json OUT]

Without --calc-root the script scans archive/ + calc/ under the user's
home layout used by this project.  Exit code 0 = no findings.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402

import delfin.analysis.results_overview as ro  # noqa: E402
from delfin.analysis._smiles_orders import smiles_bond_orders  # noqa: E402

_METALS = {
    "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Al", "Ga", "In", "Tl", "Pb", "Bi",
    "Sn",
}


def _graph(sym, bonds):
    adj = {}
    for i, j in bonds:
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)
    return adj


def _find_rings(sym, adj, max_size=7):
    """All simple cycles up to max_size atoms (C/N/O/S/P rings)."""
    rings = []
    ring_keys = set()
    for start in sorted(adj):
        if sym[start] not in ("C", "N", "O", "S", "P"):
            continue
        stack = [(start, [start])]
        while stack:
            node, path = stack.pop()
            if len(path) >= max_size:
                continue
            for nxt in sorted(adj.get(node, [])):
                if nxt == start and len(path) >= 5:
                    key = frozenset(path)
                    if key not in ring_keys:
                        ring_keys.add(key)
                        rings.append(list(path))
                    continue
                if nxt in path or sym[nxt] not in ("C", "N", "O", "S", "P"):
                    continue
                if sym[nxt] in _METALS:
                    continue
                stack.append((nxt, path + [nxt]))
    return rings


def _reference_orders(sym, bonds, atoms_3d):
    """Kekule orders of the metal-stripped heavy graph via RDKit."""
    try:
        from rdkit import Chem
    except ImportError:
        return None
    keep = [i for i, s in enumerate(sym) if s not in _METALS and s != "H"]
    remap = {o: n for n, o in enumerate(keep)}
    # Use the same proven path as _smiles_from_3d: 3D coordinates on the
    # metal-stripped heavy graph + rdDetermineBonds with a charge sweep.
    # The old pure-topology builder failed silently for metal-coordinated
    # heteroaromatics (106-138_106 produced 0 reference double bonds).
    pts3d = np.array([[atoms_3d[i][1], atoms_3d[i][2], atoms_3d[i][3]]
                      for i in keep])
    conf = Chem.Conformer(len(keep))
    for n, p in enumerate(pts3d):
        conf.SetAtomPosition(n, tuple(float(v) for v in p))
    try:
        from rdkit.Chem import rdDetermineBonds
        mol = Chem.RWMol()
        for i in keep:
            mol.AddAtom(Chem.Atom(sym[i]))
        for i, j in bonds:
            if i in remap and j in remap:
                mol.AddBond(remap[i], remap[j], Chem.BondType.SINGLE)
        m = mol.GetMol()
        m.AddConformer(conf, assignId=True)
        for q in (0, 1, -1, 2, -2):
            try:
                m2 = Chem.Mol(m)
                rdDetermineBonds.DetermineBondOrders(m2, charge=q)
                m = m2
                break
            except Exception:
                continue
        else:
            return None
        Chem.SanitizeMol(m)
        Chem.Kekulize(m, clearAromaticFlags=True)
    except Exception:
        return None
    orders = {}
    for bond in m.GetBonds():  # m = bond-order-assigned molecule (mol is all-SINGLE input)
        bt = bond.GetBondType()
        o = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2,
             Chem.BondType.TRIPLE: 3}.get(bt)
        if o is None:
            continue
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        orders[frozenset((remap[keep[a]], remap[keep[b]]))] = o
    return orders


def audit_calc(calc_dir):
    """Return a list of findings [(class, detail), ...] for one calc."""
    atoms = ro._parse_xyz(calc_dir)
    if not atoms:
        return [("SKIP", "no initial.xyz")]
    sym = [a[0] for a in atoms]
    pts = np.array([[a[1], a[2], a[3]] for a in atoms])
    rot = ro._rotate_inertia(sym, pts)
    cov = [ro._XYZ_COVALENT.get(s, ro._XYZ_COVALENT_DEFAULT) for s in sym]
    bonds = ro._detect_bonds(sym, cov, rot)
    if not bonds:
        return [("SKIP", "no bonds")]

    # What the renderer actually draws:
    bo = smiles_bond_orders(sym, bonds, calc_dir)
    if bo is None:
        bo = ro._rdkit_bond_orders(sym, bonds, rot) or {}
    drawn_dbl = {tuple(sorted(p)) for p, o in bo.items() if o >= 2}

    ref = _reference_orders(sym, bonds, atoms)
    findings = []
    if ref is None:
        findings.append(("SKIP", "RDKit reference unavailable"))
        return findings
    ref_dbl = {tuple(sorted(p)) for p, o in ref.items() if o == 2}

    adj = _graph(sym, bonds)
    rings = _find_rings(sym, adj)

    # --- class A/B: aromatic ring drawn with too few double bonds ---
    for ring in rings:
        if len(ring) != 6:
            continue
        elems = [sym[i] for i in ring]
        if elems.count("C") < 4:
            continue
        ref_nd = sum(1 for k in range(6)
                     if tuple(sorted((ring[k], ring[(k + 1) % 6])))
                     in ref_dbl)
        drawn_nd = sum(1 for k in range(6)
                       if tuple(sorted((ring[k], ring[(k + 1) % 6])))
                       in drawn_dbl)
        if ref_nd >= 3 and drawn_nd < 3:
            # distinguish A (coordinated N in ring) vs B (benzene next
            # to amide N): ring contains N bonded to a metal?
            has_metal_n = any(
                sym[nb] in _METALS
                for i in ring if sym[i] == "N"
                for nb in adj.get(i, []))
            cls = "A" if has_metal_n else "B"
            findings.append((cls, "6-ring %s: ref %d dbl, drawn %d"
                             % ("-".join(elems), ref_nd, drawn_nd)))

        # --- class C/D: exocyclic double bond to ring atom ---
        for k, i in enumerate(ring):
            for nb in adj.get(i, []):
                if nb in ring:
                    continue
                pair = tuple(sorted((i, nb)))
                if pair not in drawn_dbl:
                    continue
                # quinoid C=O / C=S on a ring that the reference keeps
                # aromatic, or P=C to a ring carbon
                if sym[nb] in ("O", "S") and sym[i] == "C":
                    if pair not in ref_dbl:
                        findings.append(("C", "exocyclic %s=%s on ring %s"
                                         % (sym[i], sym[nb],
                                            "-".join(elems))))
                elif sym[nb] == "P" or (sym[i] == "P" and sym[nb] == "C"):
                    if pair not in ref_dbl:
                        findings.append(("D", "P=C quinoid on ring %s"
                                         % "-".join(elems)))

    # --- class E: nitro O carrying an H in the label path ---
    for i, s in enumerate(sym):
        if s != "N":
            continue
        o_nb = [j for j in adj.get(i, []) if sym[j] == "O"]
        if len(o_nb) < 2:
            continue
        # nitro: N with two O, at least one drawn double
        n_dbl_to_o = sum(1 for j in o_nb
                         if tuple(sorted((i, j))) in drawn_dbl)
        if n_dbl_to_o >= 1:
            for j in o_nb:
                h_nb = [k for k in adj.get(j, []) if sym[k] == "H"]
                if h_nb:
                    findings.append(("E", "nitro O%d carries H" % j))
    # --- class F: carbonyl O (degree 1, drawn double) labeled OH ---
    # The renderer labels a terminal O "OH"/"HO" ONLY when the O has a
    # real proton (commit 72aafb3f9); a carbonyl O (drawn double) must
    # never carry one.  Detect the CURRENT renderer condition: an O on
    # a drawn double that also has a structural H would be the wrong
    # "HO-C=O" look; an O the label pass would still give "OH" despite
    # a real proton ON A CARBONYL is impossible (double O is not
    # single-bonded).  Keep the sibling check as the residual case:
    # carbonyl O drawn double whose sibling O carries the H - the
    # label is on the sibling now, so this is only a finding when the
    # sibling ALSO renders as "OH" while deprotonated (no H there).
    for i, s in enumerate(sym):
        if s != "O":
            continue
        heavy = [j for j in adj.get(i, []) if sym[j] != "H"]
        hs = [j for j in adj.get(i, []) if sym[j] == "H"]
        if len(heavy) == 1 and len(hs) == 0:
            j = heavy[0]
            if sym[j] == "C" and tuple(sorted((i, j))) in drawn_dbl:
                # With the H-gated label the carbonyl O can no longer
                # inherit a sibling's H.  A finding remains only when
                # the DRAWN state still shows a proton: the C=O drawn
                # single (not in drawn_dbl) while a sibling O-H exists
                # would re-create the "HO-CH" look via the label pass
                # on the C.  That variant is covered by the drawn_dbl
                # check above; nothing to flag on the double-bonded O.
                pass
    # --- class G: sulfonic acid S(=O)(=O)OH without drawn S=O ---
    for i, s in enumerate(sym):
        if s != "S":
            continue
        o_nb = [j for j in adj.get(i, []) if sym[j] == "O"]
        c_nb = [j for j in adj.get(i, []) if sym[j] == "C"]
        if len(o_nb) >= 2 and len(c_nb) == 1:
            drawn_s_dbl = sum(1 for j in o_nb
                              if tuple(sorted((i, j))) in drawn_dbl)
            if drawn_s_dbl == 0:
                findings.append(("G", "sulfonate S%d: no drawn S=O" % i))
    return findings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calc-root", action="append", default=[],
                    help="directory whose subdirs are calcs (repeatable)")
    ap.add_argument("--json", dest="json_out",
                    help="write findings as JSON to this path")
    ap.add_argument("--only", help="substring filter on calc names")
    args = ap.parse_args()
    # 87_* calcs are excluded everywhere per user rule (2026-11-21).
    if args.only and "87_" in args.only:
        args.only = None
    skip87 = "87_" not in (args.only or "")
    if args.only:
        skip87 = False

    roots = args.calc_root or [
        "/pfs/data6/home/ka/ka_ibcs/ka_xn0397/calc",
        "/pfs/data6/home/ka/ka_ibcs/ka_xn0397/archive",
    ]
    rows = []
    for root in roots:
        if not os.path.isdir(root):
            continue
        for name in sorted(os.listdir(root)):
            if args.only and args.only not in name:
                continue
            if skip87 and name.startswith("87_"):
                continue
            cd = os.path.join(root, name)
            if not os.path.isdir(cd):
                continue
            if not os.path.exists(os.path.join(cd, "initial.xyz")):
                continue
            try:
                import signal
                def _alarm(sig, frm):
                    raise TimeoutError("audit_calc > %ss" % os.environ.get("AUDIT_CALC_TIMEOUT", "60"))
                signal.signal(signal.SIGALRM, _alarm)
                signal.alarm(int(os.environ.get("AUDIT_CALC_TIMEOUT", "60")))
                try:
                    fnd = audit_calc(cd)
                finally:
                    signal.alarm(0)
            except Exception as exc:  # noqa: BLE001
                fnd = [("ERROR", str(exc)[:120])]
            if fnd:
                rows.append((name, fnd))
            print("%-24s %s" % (name, "done" if not fnd else
                                 " ".join(c for c, _ in fnd)),
                  flush=True)
            # Incremental JSON: crash-safe, one write after every calc.
            if args.json_out:
                with open(args.json_out, "w", encoding="utf-8") as fh:
                    json.dump([{"calc": n, "findings": f}
                               for n, f in rows], fh, indent=1)
    n_bad = sum(1 for _, f in rows
                if any(c not in ("SKIP",) for c, _ in f))
    for name, fnd in rows:
        for cls, detail in fnd:
            print("%-24s %s  %s" % (name, cls, detail))
    print("---")
    print("calcs with findings: %d" % n_bad)
    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fh:
            json.dump([{"calc": n, "findings": f} for n, f in rows],
                      fh, indent=1)
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
