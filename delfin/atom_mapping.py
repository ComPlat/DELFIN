#!/usr/bin/env python3
"""Universal atom-order mapping between two XYZ structures.

Reorders the atoms of a *target* structure so that they match the atom order of
a *reference* structure.  The two structures may be

  * two conformers / two exports of the *same* molecule (pure renumbering), or
  * the educt and the product of a reaction (a few bonds formed / broken),

including transition-metal, lanthanide and actinide complexes.  The method is
element-agnostic: covalent radii come from the full periodic table and metals
are recognised by element category, so *any* molecule or complex is handled
without editing a curated list.

Method
------
1.  Connectivity graphs from covalent radii, with bonds to *metals* EXCLUDED.
    The coordination sphere of a metal changes during a reaction and every
    distance cutoff would give an unstable bond count there, so metal bonds are
    simply left out.  Metals are unique enough to map on element alone.
2.  Automatic detection of the bond edits (formed / broken bonds) that make the
    two connectivity graphs isomorphic.
3.  Joint Weisfeiler-Lehman colour refinement -> topologically unique atoms.
4.  Kabsch superposition on the unique ("anchor") atoms.
5.  Hungarian assignment *inside* each degenerate colour class (methyl H, CH2 H,
    equivalent aryl rings, ...) -> the geometrically sensible choice.
6.  Verification that the result really is a graph isomorphism, with a VF2
    fallback (restricted by the WL colours) if it is not.

Requires: numpy, scipy, networkx, rdkit (all already DELFIN dependencies).

CLI:
    python -m delfin.atom_mapping reference.xyz target.xyz -o reordered.xyz \
        [--align] [--map mapping.txt] [--max-edits N] [-v]

Library:
    from delfin.atom_mapping import map_atoms
    result = map_atoms(sym_ref, xyz_ref, sym_tgt, xyz_tgt)
    order = result['order']          # target index for each reference atom
"""

from __future__ import annotations

import argparse
import itertools
from collections import Counter, defaultdict

import numpy as np

# ---------------------------------------------------------------------------
# Covalent radii -- periodic-table-complete.
#
# Primary source: RDKit's built-in periodic table (covers every element).
# DELFIN's own tuned metal/donor radii (delfin.manta.polyhedra.COV) are layered
# on top so we reuse DELFIN data where it exists.  A finite last-resort fallback
# guarantees an unknown / dummy label can never raise.
# ---------------------------------------------------------------------------
_FALLBACK_RADIUS = 1.5  # A -- only reached for genuinely unknown labels

try:  # DELFIN's own metal/donor radii (optional, used as an override)
    from delfin.manta.polyhedra import COV as _DELFIN_COV
except Exception:  # pragma: no cover - polyhedra is always importable in-tree
    _DELFIN_COV = {}

_PT = None


def _periodic_table():
    global _PT
    if _PT is None:
        from rdkit import Chem
        _PT = Chem.GetPeriodicTable()
    return _PT


_COV_CACHE: dict = {}


def cov_radius(symbol: str) -> float:
    """Single-bond covalent radius (A) for an element symbol, universal.

    DELFIN's tuned value wins where defined; otherwise RDKit's periodic-table
    value; otherwise a finite fallback.  Cached; case-insensitive on input.
    """
    sym = symbol.strip().capitalize()
    if sym in _COV_CACHE:
        return _COV_CACHE[sym]
    r = _DELFIN_COV.get(sym)
    if r is None:
        try:
            pt = _periodic_table()
            z = pt.GetAtomicNumber(sym)
            r = float(pt.GetRcovalent(z)) if z > 0 else None
        except Exception:
            r = None
    if not r or r <= 0.0:
        r = _FALLBACK_RADIUS
    _COV_CACHE[sym] = float(r)
    return float(r)


# ---------------------------------------------------------------------------
# Metal classification -- by element *category*, not a hand-maintained list.
#
# Everything that is NOT a non-metal / noble gas / covalent-bonding metalloid is
# treated as a metal, whose bonds are omitted from the connectivity graph.  This
# automatically covers transition metals, lanthanides, actinides, alkali /
# alkaline-earth and post-transition metals.
# ---------------------------------------------------------------------------
_NON_METALS = frozenset({
    'H', 'He',
    'B', 'C', 'N', 'O', 'F', 'Ne',
    'Si', 'P', 'S', 'Cl', 'Ar',
    'Ge', 'As', 'Se', 'Br', 'Kr',
    'Sb', 'Te', 'I', 'Xe',
    'Po', 'At', 'Rn',
})


def is_metal(symbol: str) -> bool:
    """True if ``symbol`` should be treated as a metal (its bonds are excluded
    from the connectivity graph).  Boron and the classic metalloids that form
    ordinary covalent ligand bonds (Si, Ge, As, Sb, Te, Po) are NOT metals."""
    return symbol.strip().capitalize() not in _NON_METALS


# ---------------------------------------------------------------------------
# XYZ I/O (small, self-contained -- keeps the CLI usable standalone)
# ---------------------------------------------------------------------------
def read_xyz(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    comment = lines[1] if len(lines) > 1 else ''
    syms, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        syms.append(p[0].capitalize())
        xyz.append([float(v) for v in p[1:4]])
    return syms, np.asarray(xyz, dtype=float), comment


def build_xyz(syms, xyz, comment=''):
    xyz = np.asarray(xyz, dtype=float)
    out = ['%d' % len(syms), str(comment)]
    for s, r in zip(syms, xyz):
        out.append('%-2s %15.8f %15.8f %15.8f' % (s, r[0], r[1], r[2]))
    return '\n'.join(out) + '\n'


def write_xyz(path, syms, xyz, comment=''):
    with open(path, 'w') as fh:
        fh.write(build_xyz(syms, xyz, comment))


# ---------------------------------------------------------------------------
# Connectivity
# ---------------------------------------------------------------------------
def connectivity(syms, xyz, tol=1.15):
    """Covalent connectivity graph.  Bonds involving a metal are omitted,
    because a metal's coordination sphere changes during the reaction and would
    destroy the graph isomorphism."""
    import networkx as nx

    n = len(syms)
    G = nx.Graph()
    radii = [cov_radius(s) for s in syms]
    metal = [is_metal(s) for s in syms]
    for i, s in enumerate(syms):
        G.add_node(i, el=s)
    xyz = np.asarray(xyz, dtype=float)
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    for i in range(n):
        if metal[i]:
            continue
        for j in range(i + 1, n):
            if metal[j]:
                continue
            if d[i, j] < tol * (radii[i] + radii[j]):
                G.add_edge(i, j)
    return G


def metal_contacts(syms, xyz, tol=1.35):
    """Atoms in the first coordination sphere of any metal -> reaction centre."""
    out = set()
    xyz = np.asarray(xyz, dtype=float)
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    radii = [cov_radius(s) for s in syms]
    for i, s in enumerate(syms):
        if not is_metal(s):
            continue
        for j in range(len(syms)):
            if j != i and d[i, j] < tol * (radii[i] + radii[j]):
                out.add(j)
    return out


# ---------------------------------------------------------------------------
# Bond-edit detection
# ---------------------------------------------------------------------------
def _env_signature(G, syms):
    return {i: (syms[i], tuple(sorted(syms[j] for j in G[i]))) for i in G.nodes}


def _candidate_atoms(G, syms, xyz, other_sig_counter):
    """Chemically 'suspicious' atoms: in a metal coordination sphere, or with a
    local environment whose count differs between the two structures."""
    sig = _env_signature(G, syms)
    mine = Counter(sig.values())
    bad = {s for s in set(mine) | set(other_sig_counter)
           if mine.get(s, 0) != other_sig_counter.get(s, 0)}
    cand = {i for i, s in sig.items() if s in bad}
    return cand | metal_contacts(syms, xyz)


def _elem_match(a, b):
    return a['el'] == b['el']


def find_bond_edits(G1, s1, x1, G2, s2, x2, max_edits=2):
    """Smallest set of bonds to add to G1 and/or G2 so the two graphs become
    isomorphic.  Returns (edits1, edits2) as lists of (i, j) edges."""
    import networkx as nx

    if nx.is_isomorphic(G1, G2, node_match=_elem_match):
        return [], []

    c1 = _candidate_atoms(G1, s1, x1, Counter(_env_signature(G2, s2).values()))
    c2 = _candidate_atoms(G2, s2, x2, Counter(_env_signature(G1, s1).values()))

    def pairs(G, cand, xyz, rmax=3.6):
        out = []
        xyz = np.asarray(xyz, dtype=float)
        d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
        for i, j in itertools.combinations(sorted(cand), 2):
            if not G.has_edge(i, j) and d[i, j] < rmax:
                out.append((i, j))
        return out

    p1, p2 = pairs(G1, c1, x1), pairs(G2, c2, x2)

    for k in range(1, max_edits + 1):
        for k1 in range(k + 1):
            k2 = k - k1
            for e1 in itertools.combinations(p1, k1):
                H1 = G1.copy()
                H1.add_edges_from(e1)
                for e2 in itertools.combinations(p2, k2):
                    H2 = G2.copy()
                    H2.add_edges_from(e2)
                    if H1.number_of_edges() != H2.number_of_edges():
                        continue
                    if nx.is_isomorphic(H1, H2, node_match=_elem_match):
                        return list(e1), list(e2)
    raise RuntimeError('no isomorphism found within %d bond edit(s)' % max_edits)


# ---------------------------------------------------------------------------
# Weisfeiler-Lehman colour refinement (joint alphabet over both graphs)
# ---------------------------------------------------------------------------
def wl_colours(G1, s1, G2, s2, rounds=30):
    c1 = {i: s1[i] for i in G1}
    c2 = {i: s2[i] for i in G2}
    for _ in range(rounds):
        n1 = {v: (c1[v], tuple(sorted(c1[u] for u in G1[v]))) for v in G1}
        n2 = {v: (c2[v], tuple(sorted(c2[u] for u in G2[v]))) for v in G2}
        alph = {c: k for k, c in
                enumerate(sorted(set(n1.values()) | set(n2.values()), key=str))}
        n1 = {v: alph[c] for v, c in n1.items()}
        n2 = {v: alph[c] for v, c in n2.items()}
        stable = (len(set(n1.values())) == len(set(c1.values())) and
                  len(set(n2.values())) == len(set(c2.values())))
        c1, c2 = n1, n2
        if stable:
            break
    return c1, c2


# ---------------------------------------------------------------------------
# Kabsch superposition
# ---------------------------------------------------------------------------
def kabsch_rt(P, Q):
    """Rotation + translation that puts P onto Q (both N x 3, same order)."""
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, pc, qc


def apply_rt(X, R, pc, qc):
    return (np.asarray(X, dtype=float) - pc) @ R.T + qc


def is_isomorphism(mapping, G1, s1, G2, s2):
    if len(set(mapping.values())) != len(mapping):
        return False
    if G1.number_of_edges() != G2.number_of_edges():
        return False
    if any(s1[i] != s2[mapping[i]] for i in mapping):
        return False
    return all(G2.has_edge(mapping[u], mapping[v]) for u, v in G1.edges())


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def map_atoms(sym_ref, xyz_ref, sym_tgt, xyz_tgt, max_edits=2, verbose=False):
    """Map the atoms of the *target* onto the *reference* atom order.

    Parameters
    ----------
    sym_ref, sym_tgt : list[str]      element symbols
    xyz_ref, xyz_tgt : (N, 3) arrays  Cartesian coordinates
    max_edits : int                   max number of bonds formed/broken to try
    verbose : bool                    print progress

    Returns
    -------
    dict with keys:
      order            list[int] -- target index for each reference atom i
      bond_edits_ref   list[(i,j,el_i,el_j)] bonds added to the reference graph
      bond_edits_tgt   list[(i,j,el_i,el_j)] bonds added to the target graph
      n_bond_edits     int
      rmsd             float -- RMSD after superposition on the full mapping
      method           'xyzmap-verified' | 'xyzmap-vf2' | 'failed'
      verified         bool
      reordered_target_xyz  str -- XYZ block of the target in reference order
    """
    sym_ref = [s.capitalize() for s in sym_ref]
    sym_tgt = [s.capitalize() for s in sym_tgt]
    x1 = np.asarray(xyz_ref, dtype=float)
    x2 = np.asarray(xyz_tgt, dtype=float)

    if Counter(sym_ref) != Counter(sym_tgt):
        raise ValueError(
            'different molecular formulae: %s vs %s'
            % (dict(Counter(sym_ref)), dict(Counter(sym_tgt)))
        )
    n = len(sym_ref)

    G1 = connectivity(sym_ref, x1)
    G2 = connectivity(sym_tgt, x2)

    e1, e2 = find_bond_edits(G1, sym_ref, x1, G2, sym_tgt, x2, max_edits)
    H1, H2 = G1.copy(), G2.copy()
    H1.add_edges_from(e1)
    H2.add_edges_from(e2)

    def _fmt(es, syms):
        return [(i, j, syms[i], syms[j]) for i, j in es]

    if verbose:
        print('bond edits: +%s in reference, +%s in target'
              % (_fmt(e1, sym_ref), _fmt(e2, sym_tgt)))

    c1, c2 = wl_colours(H1, sym_ref, H2, sym_tgt)
    cls1, cls2 = defaultdict(list), defaultdict(list)
    for i in range(n):
        cls1[c1[i]].append(i)
        cls2[c2[i]].append(i)
    if sorted(map(len, cls1.values())) != sorted(map(len, cls2.values())):
        raise RuntimeError('WL colour classes are incompatible between structures')
    if verbose:
        sings = sum(1 for v in cls1.values() if len(v) == 1)
        print('WL classes: %d, singletons: %d/%d' % (len(cls1), sings, n))

    # anchor on topologically unique atoms
    anchor = {v[0]: cls2[c][0] for c, v in cls1.items() if len(v) == 1}
    if len(anchor) < 3:
        raise RuntimeError('too few unique anchor atoms for superposition')
    ke = sorted(anchor)
    R, pc, qc = kabsch_rt(x2[[anchor[i] for i in ke]], x1[ke])
    x2a = apply_rt(x2, R, pc, qc)

    from scipy.optimize import linear_sum_assignment

    best, best_rmsd = None, np.inf
    for _ in range(10):
        m = {}
        for c, ve in cls1.items():
            vp = cls2[c]
            if len(ve) == 1:
                m[ve[0]] = vp[0]
                continue
            cost = np.linalg.norm(x1[ve][:, None, :] - x2a[vp][None, :, :], axis=-1)
            row, col = linear_sum_assignment(cost)
            for a, b in zip(row, col):
                m[ve[a]] = vp[b]
        if not is_isomorphism(m, H1, sym_ref, H2, sym_tgt):
            break
        order = [m[i] for i in range(n)]
        R, pc, qc = kabsch_rt(x2[order], x1)
        x2a = apply_rt(x2, R, pc, qc)
        rmsd = float(np.sqrt(((x2a[order] - x1) ** 2).sum(1).mean()))
        if rmsd < best_rmsd - 1e-8:
            best, best_rmsd = order, rmsd
        else:
            break

    method = 'xyzmap-verified' if best is not None else 'failed'

    # fallback / safety net: VF2 restricted by WL colours (both element + colour)
    if best is None:
        import networkx as nx

        if verbose:
            print('per-class assignment inconsistent -> VF2 fallback')
        # bake the WL colour into the node attributes so VF2 respects it
        for i in range(n):
            H1.nodes[i]['wl'] = c1[i]
            H2.nodes[i]['wl'] = c2[i]
        nm = lambda a, b: a['el'] == b['el'] and a['wl'] == b['wl']
        gm = nx.algorithms.isomorphism.GraphMatcher(H1, H2, node_match=nm)
        for k, mm in enumerate(gm.isomorphisms_iter()):
            order = [mm[i] for i in range(n)]
            R, pc, qc = kabsch_rt(x2[order], x1)
            rmsd = float(np.sqrt(((apply_rt(x2, R, pc, qc)[order] - x1) ** 2).sum(1).mean()))
            if rmsd < best_rmsd:
                best, best_rmsd = order, rmsd
                method = 'xyzmap-vf2'
            if k > 50000:
                break

    if best is None:
        raise RuntimeError('could not construct a valid atom mapping')

    verified = is_isomorphism({i: best[i] for i in range(n)}, H1, sym_ref, H2, sym_tgt)
    if not verified:
        raise RuntimeError('constructed mapping failed isomorphism verification')

    R, pc, qc = kabsch_rt(x2[best], x1)
    disp = np.linalg.norm(apply_rt(x2, R, pc, qc)[best] - x1, axis=1)
    if verbose:
        moved = sum(1 for i in range(n) if best[i] != i)
        print('method=%s  RMSD=%.3f A (max %.3f A on atom %d)  lines moved %d/%d'
              % (method, best_rmsd, disp.max(), int(disp.argmax()) + 1, moved, n))

    reordered_syms = [sym_tgt[j] for j in best]
    reordered_xyz = x2[best]
    return {
        'order': list(best),
        'bond_edits_ref': _fmt(e1, sym_ref),
        'bond_edits_tgt': _fmt(e2, sym_tgt),
        'n_bond_edits': len(e1) + len(e2),
        'rmsd': best_rmsd,
        'method': method,
        'verified': verified,
        'displacement': disp.tolist(),
        'reordered_target_xyz': build_xyz(reordered_syms, reordered_xyz,
                                          comment='Reordered to reference numbering'),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Reorder a target XYZ onto a reference XYZ atom order '
                    '(renumbering or reaction educt<->product, metal-aware).')
    ap.add_argument('reference')
    ap.add_argument('target')
    ap.add_argument('-o', '--out', default='reordered.xyz')
    ap.add_argument('--align', action='store_true',
                    help='also rigid-body superimpose the output onto the reference')
    ap.add_argument('--max-edits', type=int, default=2)
    ap.add_argument('--map', default=None, help='write a mapping table to this file')
    ap.add_argument('-v', '--verbose', action='store_true')
    a = ap.parse_args(argv)

    s1, x1, _ = read_xyz(a.reference)
    s2, x2, cm = read_xyz(a.target)
    res = map_atoms(s1, x1, s2, x2, max_edits=a.max_edits, verbose=a.verbose)
    order = res['order']

    xnew = x2[order]
    if a.align:
        R, pc, qc = kabsch_rt(xnew, x1)
        xnew = apply_rt(xnew, R, pc, qc)
    write_xyz(a.out, [s2[j] for j in order], xnew, cm)

    print('method=%s  bond_edits=%d  rmsd=%.3f A  -> %s'
          % (res['method'], res['n_bond_edits'], res['rmsd'], a.out))
    for i, j, ei, ej in res['bond_edits_ref']:
        print('  formed/broken bond (reference): %s%d-%s%d' % (ei, i + 1, ej, j + 1))
    for i, j, ei, ej in res['bond_edits_tgt']:
        print('  formed/broken bond (target):    %s%d-%s%d' % (ei, i + 1, ej, j + 1))

    if a.map:
        R, pc, qc = kabsch_rt(x2[order], x1)
        disp = np.linalg.norm(apply_rt(x2, R, pc, qc)[order] - x1, axis=1)
        with open(a.map, 'w') as fh:
            fh.write('# new_index  element  old_index  displacement/A\n')
            for i, j in enumerate(order):
                fh.write('%5d  %-2s  %5d  %7.3f\n' % (i + 1, s2[j], j + 1, disp[i]))


if __name__ == '__main__':
    main()
