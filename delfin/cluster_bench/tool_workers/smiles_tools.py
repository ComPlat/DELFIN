"""rxembed, Molassembler and MetalloGen: one build per tool, run in the tool's own interpreter.

Shared by the dashboard / CONTROL build (``delfin.common.external_builders``, which imports this
module in the worker it starts) and by the cluster batch worker (``smiles_worker.py``, which
imports it as a sibling file).  Standard library only at import time, so it loads in each tool's
environment whatever else that environment has (Molassembler's has no RDKit), and it imports
nothing else from DELFIN.

Each ``run_<tool>`` returns ``(frames, extra)``: ``frames`` a list of
``{'label', 'symbols', 'coords' (Angstrom), 'energy' (or None)}`` in the tool's order, ``extra`` a
dict with what the run reports about itself.  A SMILES the tool cannot take raises
:class:`NotExpressible`; nothing ever falls back to another builder.

rxembed (github.com/aligfellow/rxembed, MIT; pinned commit 9c0ba5a6, version 0.2.0)
    ``rx.embed(rx.isomers(smiles))`` exactly as its README's "Metal Complexes" section: every
    stereoisomer and every coordination arrangement on rxembed's default polyhedron of the
    coordination number, its default conformer count and seed.  Python >= 3.12, RDKit >= 2026.9.1.

Molassembler (Sobez & Reiher, JCIM 2020, 60, 3884; scine-molassembler 3.0.1, BSD-3)
    ``io.experimental.from_smiles_multiple``; the component that holds a metal (else the largest)
    is taken.  Its SMILES reader does not parse dative bonds (``->``), so such a SMILES is not
    expressible; write the complex charge-separated.  Every stereoisomer is enumerated by assigning
    each unassigned atom and bond stereopermutator in turn (the molecule is copied after each
    assignment because ranking can change); each leaf is canonicalised and kept only when no
    earlier isomer compares equal (meso forms once), at most ``max_isomers`` (1000) distinct ones
    and 50000 leaves -- the same enumeration as the paper benchmark's runner.  Each isomer gets
    ``dg.generate_ensemble(isomer, num_confs, seed + i, Configuration())``; coordinates come in Bohr
    and are converted to Angstrom.  Molassembler picks one polyhedron per centre from its graph
    model (Pt(II)/Pd(II) with four donors: tetrahedral).

MetalloGen (Lee, Park, Park & Kim, JCIM 2025, 65, 11878; github.com/kyunghoonlee777/MetalloGen,
BSD-3; pinned commit fa26e2ad)
    Its default calculator is xTB through Gaussian (``xtb-gaussian``); DELFIN uses only its ORCA
    calculator (``-c orca``: ORCA's ``! XTB2``, GFN2-xTB through ORCA's ``otool_xtb``), which needs
    no Gaussian.  MetalloGen builds the one arrangement its m-SMILES writes (donor -> numbered
    site of a named geometry); the m-SMILES are written here from DELFIN's split of the SMILES,
    one per distinct arrangement of the donor sites on each geometry (proper rotations of the
    polyhedron, automorphisms of each ligand and swaps of identical ligands identify
    arrangements, mirror images are kept; donors of one ligand at most four bonds apart never sit
    more than 135 degrees apart; at most ``max_arrangements`` per geometry, ``truncated`` says
    when that cut).  Each m-SMILES is built as MetalloGen's command line builds it (``run.main``:
    its embedding scales, ``TMCGenerator.sample_conformer``, Butina clustering at 0.5 A, final
    relaxation with the calculator, ``num_conformer`` cap), with ``random`` and ``numpy.random``
    seeded to ``seed + index`` first; energies are the calculator's (Hartree, GFN2-xTB).  Paths
    with blanks or shell characters are refused: MetalloGen hands them to the shell unquoted.
"""

from __future__ import annotations

import copy
import itertools
import math
import os
import shutil
import sys
import tempfile

BOHR_TO_ANGSTROM = 0.529177210903

#: Defaults, the same for the dashboard, CONTROL and the cluster batch.
RXEMBED_DEFAULTS = {"n": None, "seed": None}
MOLASSEMBLER_DEFAULTS = {"num_confs": 10, "max_isomers": 1000, "seed": 42}
METALLOGEN_DEFAULTS = {"geometries": "common", "num_conformer": 10, "final_relax": 1,
                       "max_arrangements": 48, "calculator": "orca", "seed": 42}

#: MetalloGen geometries (its globalvars names, typos included) per number of donor sites.
#: ``common``: the polyhedra the other constructors build; ``all``: every one MetalloGen knows.
METALLOGEN_GEOMS = {
    "common": {2: ["2_linear"], 3: ["3_trigonal_planar"],
               4: ["4_square_planar", "4_tetrahedral"],
               5: ["5_trigonal_bipyramidal", "5_square_pyramidal"],
               6: ["6_octahedral"], 7: ["7_pentagonal_bipyramidal"],
               8: ["8_squre_antiprismatic"]},
    "all": None,  # MetalloGen's own CN_known_geometries_dict, read in its environment
}

#: Donors of one ligand this many bonds apart or closer must not sit trans (> 135 degrees).
CHELATE_TRANS_MAX_BONDS = 4

_ST_METALS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am",
}


class NotExpressible(ValueError):
    """A SMILES the tool cannot take, with the reason."""


def st_read_xyz_frames(path):
    """``[(comment, symbols, coords)]`` of a multi-frame xyz file."""
    out = []
    with open(path) as fh:
        lines = fh.read().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        syms, xyz = [], []
        for row in lines[i + 2:i + 2 + n]:
            p = row.split()
            syms.append(p[0])
            xyz.append([float(p[1]), float(p[2]), float(p[3])])
        out.append((lines[i + 1] if i + 1 < len(lines) else "", syms, xyz))
        i += 2 + n
    return out


# ---------------------------------------------------------------------------
# rxembed
# ---------------------------------------------------------------------------

def run_rxembed(smiles, work, options=None):
    import rxembed as rx

    opts = dict(RXEMBED_DEFAULTS, **(options or {}))
    kw = {k: opts[k] for k in ("n", "seed") if opts.get(k) is not None}
    try:
        species = rx.isomers(smiles)
    except (ValueError, RuntimeError) as exc:
        raise NotExpressible("rxembed cannot read the SMILES: %s" % str(exc)[:300])
    res = rx.embed(species, **kw)
    d = tempfile.mkdtemp(prefix="rxembed_", dir=work)
    # EnsembleSet.dump: one multi-frame xyz per species, named "sp_<index>_<species name>".
    paths = res.dump(os.path.join(d, "sp.xyz"))
    frames = []
    for p in (paths if isinstance(paths, (list, tuple)) else [paths]):
        if not os.path.exists(p) or not os.path.getsize(p):
            continue
        label = os.path.basename(p)[:-4]
        label = label[3:] if label.startswith("sp_") else label
        for j, (_c, syms, xyz) in enumerate(st_read_xyz_frames(p)):
            frames.append({"label": "%s-conf%d" % (label, j), "symbols": syms, "coords": xyz,
                           "energy": None})
    errors = [str(e)[:200] for e in (getattr(res, "errors", ()) or ())]
    return frames, {"n_species": len(species), "n_species_failed": len(errors),
                    "species_errors": errors[:20]}


# ---------------------------------------------------------------------------
# Molassembler
# ---------------------------------------------------------------------------

def _masm_element(m, i):
    return str(m.graph.element_type(i)).split(".")[-1]


def _masm_unassigned(m):
    for p in m.stereopermutators.atom_stereopermutators():
        if p.assigned is None and p.num_assignments > 1:
            return p.placement, p.num_assignments
    for p in m.stereopermutators.bond_stereopermutators():
        if p.assigned is None and p.num_assignments > 1:
            return p.placement, p.num_assignments
    return None


#: Raw assignment leaves visited before giving up (duplicates included), as in the benchmark runner.
MOLASSEMBLER_MAX_LEAVES = 50000


def _masm_enumerate(m, out, cap, seen=None, leaves=None):
    """Every distinct stereoisomer, depth-first (identical to the paper benchmark's mas_one.py).

    Assigning stereopermutators one after another reaches the same isomer more than once (meso:
    FC(Cl)C(F)Cl gave 4 instead of 3), so each leaf is canonicalised on a copy and kept only if no
    earlier isomer compares equal; the cap counts distinct isomers.
    """
    if seen is None:
        seen, leaves = [], [0]
    if len(out) >= cap:
        return
    u = _masm_unassigned(m)
    if u is None:
        leaves[0] += 1
        if leaves[0] > MOLASSEMBLER_MAX_LEAVES:
            raise NotExpressible("more than %d assignment leaves" % MOLASSEMBLER_MAX_LEAVES)
        key = copy.copy(m)
        key.canonicalize()
        if not any(key == s for s in seen):
            seen.append(key)
            out.append(m)
        return
    for k in range(u[1]):
        c = copy.copy(m)
        c.assign_stereopermutator(u[0], k)
        _masm_enumerate(c, out, cap, seen, leaves)
        if len(out) >= cap:
            return


def run_molassembler(smiles, work, options=None):
    if "->" in smiles or "<-" in smiles:
        raise NotExpressible("Molassembler's SMILES reader does not parse dative bonds (->); "
                             "write the complex charge-separated, e.g. [NH3+][Pt-2]([NH3+])(Cl)Cl")
    import scine_molassembler as masm

    opts = dict(MOLASSEMBLER_DEFAULTS, **(options or {}))
    # The enumeration recurses once per stereopermutator.
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    try:
        mols = masm.io.experimental.from_smiles_multiple(smiles)
    except RuntimeError as exc:
        raise NotExpressible("Molassembler cannot read the SMILES: %s"
                             % (str(exc).strip()[:200] or "parsing failure"))
    with_metal = [m for m in mols
                  if any(_masm_element(m, i) in _ST_METALS for i in range(m.graph.V))]
    mol = max(with_metal or mols, key=lambda x: x.graph.V)
    cap = int(opts["max_isomers"])
    isomers = []
    _masm_enumerate(mol, isomers, cap + 1)
    if len(isomers) > cap:
        raise NotExpressible("more than %d stereoisomers" % cap)
    cfg = masm.dg.Configuration()
    frames, n_err = [], 0
    for i, iso in enumerate(isomers):
        els = [_masm_element(iso, a) for a in range(iso.graph.V)]
        ens = masm.dg.generate_ensemble(iso, int(opts["num_confs"]), int(opts["seed"]) + i, cfg)
        for j, pos in enumerate(ens):
            if isinstance(pos, masm.dg.Error):
                n_err += 1
                continue
            frames.append({"label": "iso%d-conf%d" % (i, j), "symbols": els,
                           "coords": [[x * BOHR_TO_ANGSTROM, y * BOHR_TO_ANGSTROM,
                                       z * BOHR_TO_ANGSTROM] for x, y, z in pos],
                           "energy": None})
    return frames, {"n_isomers": len(isomers), "n_dg_errors": n_err,
                    "n_components": len(mols)}


# ---------------------------------------------------------------------------
# MetalloGen: DELFIN's split -> m-SMILES, one per arrangement
# ---------------------------------------------------------------------------

def _st_unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return [c / n for c in v]


def _st_dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def _st_orthonormal(a, b):
    """Orthonormal frame (rows) from two non-collinear unit vectors."""
    e1 = _st_unit(a)
    e2 = _st_unit([b[i] - _st_dot(b, e1) * e1[i] for i in range(3)])
    e3 = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
          e1[0] * e2[1] - e1[1] * e2[0]]
    return [e1, e2, e3]


def site_rotations(vectors, tol=1e-3):
    """Every permutation of the sites a proper rotation of the polyhedron makes.

    ``P[i] = j``: the rotation takes site i onto site j.  Mirror images are not identified.
    """
    v = [_st_unit(x) for x in vectors]
    n = len(v)
    if n == 1:
        return [(0,)]
    j0 = next((j for j in range(1, n) if abs(abs(_st_dot(v[0], v[j])) - 1.0) > 1e-6), None)
    if j0 is None:                     # every site on one line (linear): identity and the flip
        return sorted({tuple(range(n)), tuple(reversed(range(n)))})
    cos0 = _st_dot(v[0], v[j0])
    f1 = _st_orthonormal(v[0], v[j0])
    perms = set()
    for a in range(n):
        for b in range(n):
            if a == b or abs(_st_dot(v[a], v[b]) - cos0) > tol:
                continue
            f2 = _st_orthonormal(v[a], v[b])
            # R = F2^T F1 takes the frame of (v0, vj0) onto that of (va, vb).
            R = [[sum(f2[k][r] * f1[k][c] for k in range(3)) for c in range(3)] for r in range(3)]
            perm = []
            for i in range(n):
                w = [sum(R[r][c] * v[i][c] for c in range(3)) for r in range(3)]
                hit = [j for j in range(n) if sum((w[k] - v[j][k]) ** 2 for k in range(3)) < 1e-4]
                if len(hit) != 1:
                    break
                perm.append(hit[0])
            if len(perm) == n and len(set(perm)) == n:
                perms.add(tuple(perm))
    return sorted(perms)


def _mg_group_perms(m, groups):
    """The permutations of a ligand's donor sites that an automorphism of the ligand makes.

    From ``GetSubstructMatches(mol, mol, uniquify=False, useChirality=False)``: each match maps
    atom i onto match[i]; it is a symmetry of the bound ligand only when it maps every donor site
    onto a donor site.  ``p[g]`` = the site group g goes to.  The identity is always among them.
    """
    index = {tuple(g): k for k, g in enumerate(groups)}
    perms = {tuple(range(len(groups)))}
    for match in m.GetSubstructMatches(m, uniquify=False, useChirality=False, maxMatches=100000):
        p = []
        for g in groups:
            k = index.get(tuple(sorted(match[a] for a in g)))
            if k is None:
                break
            p.append(k)
        if len(p) == len(groups):
            perms.add(tuple(p))
    return sorted(perms)


def metallogen_ligands(spec):
    """The donor sites of the split: ``[(ligand index, [donor atoms], group index), ...]`` -- a
    hapto group (donors of one ligand bonded to each other) is one site -- and per-ligand RDKit
    data, with the donor-site permutations its automorphisms make (``gperms``)."""
    from rdkit import Chem

    sites, ligs = [], []
    for li, lig in enumerate(spec["ligands"]):
        m = Chem.MolFromSmiles(lig["smiles"])
        if m is None:
            raise NotExpressible("ligand %s cannot be read by the RDKit of MetalloGen" % lig["smiles"])
        cs = set(lig["coordList"])
        seen, groups = set(), []
        for d in sorted(cs):
            if d in seen:
                continue
            comp, stack = [], [d]
            seen.add(d)
            while stack:
                a = stack.pop()
                comp.append(a)
                for nb in m.GetAtomWithIdx(a).GetNeighbors():
                    if nb.GetIdx() in cs and nb.GetIdx() not in seen:
                        seen.add(nb.GetIdx())
                        stack.append(nb.GetIdx())
            groups.append(sorted(comp))
        ligs.append({"mol": m, "smiles": Chem.MolToSmiles(m), "groups": groups,
                     "gperms": _mg_group_perms(m, groups), "dist": Chem.GetDistanceMatrix(m)})
        for gi, g in enumerate(groups):
            sites.append((li, g, gi))
    return sites, ligs


def _mg_canonical_key(assign, sites, ligs, perm, gchoice):
    """One donor-site -> polyhedron-site assignment, as seen after one rotation and one
    automorphism per ligand.

    ``assign[k]`` = polyhedron site of donor site k; ``gchoice[li]`` = the donor-site permutation
    of ligand li applied.  A polyhedron site is named by its ligand, the (permuted) donor site of
    that ligand and the ligand instance; identical instances are interchangeable, so an instance
    is numbered by where its first donor appears.
    """
    n = len(assign)
    at_site = [None] * n
    for k, s in enumerate(assign):
        at_site[perm[s]] = k
    inst, key = {}, []
    for s in range(n):
        li, _g, gi = sites[at_site[s]]
        name = ligs[li]["smiles"]
        if li not in inst:
            inst[li] = sum(1 for x in inst if ligs[x]["smiles"] == name)
        key.append((name, gchoice[li][gi], inst[li]))
    return tuple(key)


def metallogen_arrangements(spec, vectors, cap=None):
    """Distinct assignments of the donor sites to the sites of one polyhedron.

    ``([assign, ...], sites, ligs, truncated)`` with ``assign[k]`` = 0-based polyhedron site of
    donor site k (order of :func:`metallogen_ligands`).  Two assignments are one arrangement when
    a proper rotation of the polyhedron, an automorphism of each ligand and a swap of identical
    ligands turn one into the other (mirror images stay apart).  Donors of one ligand at most
    ``CHELATE_TRANS_MAX_BONDS`` bonds apart never sit more than 135 degrees apart.
    ``truncated`` is True when *cap* stopped the enumeration before a further distinct one.
    """
    sites, ligs = metallogen_ligands(spec)
    n = len(sites)
    if n != len(vectors):
        raise ValueError("%d donor sites for a %d-site polyhedron" % (n, len(vectors)))
    v = [_st_unit(x) for x in vectors]
    close = []
    for a in range(n):
        for b in range(a + 1, n):
            if sites[a][0] != sites[b][0]:
                continue
            d = ligs[sites[a][0]]["dist"]
            if min(d[i][j] for i in sites[a][1] for j in sites[b][1]) <= CHELATE_TRANS_MAX_BONDS:
                close.append((a, b))
    rots = site_rotations(v)
    choices = list(itertools.product(*[lig["gperms"] for lig in ligs]))
    limit = math.cos(math.radians(135.0))
    seen, out, truncated = set(), [], False
    for assign in itertools.permutations(range(n)):
        if any(_st_dot(v[assign[a]], v[assign[b]]) < limit for a, b in close):
            continue
        key = min(_mg_canonical_key(assign, sites, ligs, p, c) for p in rots for c in choices)
        if key in seen:
            continue
        if cap and len(out) >= cap:
            truncated = True
            break
        seen.add(key)
        out.append(list(assign))
    return out, sites, ligs, truncated


def metallogen_msmiles(spec, geometry, assign, sites, ligs):
    """The m-SMILES of one arrangement: ``[M+ox]|ligand|...|geometry``; sites are 1-based."""
    from rdkit import Chem

    per_lig = {}
    for k, (li, g, _gi) in enumerate(sites):
        for a in g:
            per_lig.setdefault(li, {})[a] = assign[k] + 1
    ox = int(spec.get("metal_ox") or 0)
    parts = ["[%s+%d]" % (spec["metal"], ox) if ox > 0 else "[%s]" % spec["metal"]]
    for li, lig in enumerate(ligs):
        m = Chem.Mol(lig["mol"])
        for a, site in per_lig.get(li, {}).items():
            m.GetAtomWithIdx(a).SetAtomMapNum(site)
        parts.append(Chem.MolToSmiles(m, canonical=False))
    parts.append(geometry)
    return "|".join(parts)


def metallogen_jobs(spec, options=None, known=None, vectors=None, report=None):
    """``[(geometry, arrangement index, m-SMILES), ...]`` for every geometry that fits.

    *known* / *vectors*: MetalloGen's ``CN_known_geometries_dict`` and
    ``known_geometries_vector_dict`` (read from it in its environment).  *report*, a dict, gets
    ``truncated`` (whether ``max_arrangements`` cut any geometry short) and
    ``truncated_geometries``.
    """
    opts = dict(METALLOGEN_DEFAULTS, **(options or {}))
    if opts["geometries"] not in METALLOGEN_GEOMS:
        raise ValueError("unknown MetalloGen geometry set %r (known: %s)"
                         % (opts["geometries"], ", ".join(sorted(METALLOGEN_GEOMS))))
    sites, _ligs = metallogen_ligands(spec)
    n = len(sites)
    table = METALLOGEN_GEOMS[opts["geometries"]]
    geoms = list((known or {}).get(n, [])) if table is None else list(table.get(n, []))
    geoms = [g for g in geoms if g in (vectors or {}) and len(vectors[g]) == n]
    if not geoms:
        raise NotExpressible("MetalloGen has no geometry for %d donor sites (a hapto ligand is "
                             "one site)" % n)
    jobs, cut = [], []
    for geom in geoms:
        arrs, sites, ligs, truncated = metallogen_arrangements(
            spec, vectors[geom], cap=int(opts["max_arrangements"]))
        if truncated:
            cut.append(geom)
        for k, assign in enumerate(arrs):
            jobs.append((geom, k, metallogen_msmiles(spec, geom, assign, sites, ligs)))
    if report is not None:
        report["truncated"] = bool(cut)
        report["truncated_geometries"] = cut
    if not jobs:
        raise NotExpressible("no arrangement of the donors fits %s (a chelate cannot span "
                             "trans positions)" % ", ".join(geoms))
    return jobs


def _metallogen_one(msmiles, work, opts, orca):
    """One m-SMILES through MetalloGen as its command line builds it (run.main)."""
    import numpy as np
    from MetalloGen import om
    from MetalloGen.Calculator import orca as mg_orca
    from MetalloGen.run import TMCGenerator, clean_working_directory
    from MetalloGen.utils import clustering

    os.makedirs(work, exist_ok=True)
    complex_ = om.get_om_from_modified_smiles(msmiles)
    calc = mg_orca.Orca(command=orca)
    calc.change_working_directory(work)
    nc = int(opts["num_conformer"])
    if nc < 10:
        scales = [0.70, 0.80, 0.90, 1.00, 1.10]
    else:
        scales = [0.7 + ((0.8 + 0.1 * (nc - 10)) / (nc - 2)) * k for k in range(int(nc / 2))]
    mols = []
    for scale in scales:
        mols.extend(TMCGenerator(calc, scale, True).sample_conformer(complex_))
    if not mols:
        return []
    symbols = [a.get_element() for a in complex_.get_atom_list()]
    clusters, _D = clustering.cluster_conformers_butina(
        [m.get_coordinate_list() for m in mols], symbols, cutoff=0.5, exclude_H=True)
    out = []
    for i, mol in enumerate(mols[c[0]] for c in clusters):
        if mol.energy == 1e6:
            continue
        if int(opts["final_relax"]):
            calc.optimize_geometry(mol, file_name="final_relax_%d" % i, initial_hessian=None,
                                   save_directory=None)
            clean_working_directory(work)
        syms, xyz = [], []
        for atom in mol.atom_list:
            p = atom.get_content().split()
            syms.append(p[0])
            xyz.append([float(p[1]), float(p[2]), float(p[3])])
        e = mol.energy
        out.append((syms, xyz, None if e is None or not np.isfinite(e) else float(e)))
        if len(out) >= nc:
            break
    return out


#: Characters MetalloGen's os.system() calls (``mkdir``, ``mv``, ``orca x.com > x.log``) would
#: read as shell syntax; a path with one of them is refused instead of run.
_SHELL_UNSAFE = set(" \t\n\r;&|$`'\"<>()*?![]{}#~\\")


def _mg_check_path(what, path):
    bad = sorted(set(str(path)) & _SHELL_UNSAFE)
    if bad:
        raise NotExpressible("MetalloGen passes the %s to the shell unquoted, and %r contains %s; "
                             "use a path without blanks or shell characters"
                             % (what, str(path), " ".join(repr(c) for c in bad)))


def _mg_append_partial(path, frames):
    """Append frames to a multi-frame xyz as soon as an arrangement is done (kept on a time-out)."""
    with open(path, "a") as fh:
        for f in frames:
            fh.write("%d\n%s E=%s\n" % (len(f["symbols"]), f["label"],
                                        "na" if f["energy"] is None else "%.6f" % f["energy"]))
            for s, (x, y, z) in zip(f["symbols"], f["coords"]):
                fh.write("%-2s %14.6f %14.6f %14.6f\n" % (s, x, y, z))


def run_metallogen(spec, work, options=None):
    opts = dict(METALLOGEN_DEFAULTS, **(options or {}))
    if opts.get("calculator", "orca") != "orca":
        raise NotExpressible("DELFIN runs MetalloGen only with its ORCA calculator (GFN2-xTB "
                             "through ORCA); its default xtb_gaussian route needs Gaussian")
    orca = opts.get("orca") or shutil.which("orca")
    if not orca or not os.access(orca, os.X_OK):
        raise NotExpressible("MetalloGen's Gaussian-free route needs ORCA (orca and its "
                             "otool_xtb); no ORCA executable was found")
    orca = os.path.abspath(orca)
    _mg_check_path("ORCA path", orca)
    _mg_check_path("working directory", os.path.abspath(work))
    os.environ["PATH"] = os.path.dirname(orca) + os.pathsep + os.environ.get("PATH", "")
    import random

    import numpy as np
    from MetalloGen import globalvars as gv

    report = {}
    jobs = metallogen_jobs(spec, opts, gv.CN_known_geometries_dict,
                           gv.known_geometries_vector_dict, report)
    frames, per = [], {}
    for idx, (geom, k, msmiles) in enumerate(jobs):
        tag = "%s-arr%d" % (geom, k)
        # MetalloGen draws its embeddings from random and numpy.random; one seed per m-SMILES
        # makes the embeddings repeat (the xTB/ORCA relaxations can still differ in the last
        # digits, so a frame count can change when a near-duplicate falls either side of the
        # 0.5 A clustering cut-off).
        random.seed(int(opts["seed"]) + idx)
        np.random.seed(int(opts["seed"]) + idx)
        try:
            got = _metallogen_one(msmiles, os.path.join(work, tag), opts, orca)
        except (Exception, SystemExit) as exc:   # MetalloGen also exit()s
            per[tag] = "%s: %s" % (type(exc).__name__, str(exc)[:160])
            continue
        per[tag] = "%d conformers (%s)" % (len(got), msmiles)
        done = [{"label": "%s-conf%d" % (tag, j), "symbols": syms, "coords": xyz, "energy": e}
                for j, (syms, xyz, e) in enumerate(got)]
        frames.extend(done)
        if opts.get("partial_xyz") and done:
            _mg_append_partial(opts["partial_xyz"], done)
    return frames, {"per_geometry": per, "m_smiles": [j[2] for j in jobs],
                    "truncated": report.get("truncated", False),
                    "truncated_geometries": report.get("truncated_geometries", [])}


#: Which runner takes the SMILES as written and which takes DELFIN's split (a spec).
RUNNERS = {"rxembed": run_rxembed, "molassembler": run_molassembler,
           "metallogen": run_metallogen}
TAKES_SPEC = {"metallogen"}
