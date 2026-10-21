"""Ligand conformers of the FF-free constructor: cached MMFF-free relaxation, degenerate symmetrisation, clash count, torsion and joint declash frames, guarded refinement and sphere flex.

Moved verbatim from delfin/manta/assemble_complex.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem

from delfin.manta.assemble_donor_plane import (
    _ffree_flag,
)
from delfin.manta.assemble_ligand_embed import (
    SEED,
)


def _ligand_3d_from_mol(frag_mol):
    """Embed a ligand fragment mol (heavy-atom indices preserved under AddHs)."""
    m = Chem.AddHs(frag_mol)
    if AllChem.EmbedMolecule(m, randomSeed=SEED) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(m)
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    return syms, m.GetConformer().GetPositions(), m


# Process-local conformer memo (deterministic embed -> cacheable).  Keyed by the
# fragment's canonical SMILES + k, so the ENSEMBLE builder re-uses one ligand embed
# across all its variant builds instead of re-running ETKDG per variant (the embed of
# a large flexible η-substituent dominates the per-build cost).  Transparent: the same
# input always returns the same (deterministic) result, so behaviour is unchanged —
# this is a speedup, not a semantic change.  Bounded to avoid unbounded growth.
_CONF_CACHE = {}


_CONF_CACHE_MAX = 256


def _relax_confs_ffree(m, cids):
    """Minimise each conformer with U_total instead of MMFF.

    Routed through ``variational_refine`` rather than a second copy of the minimiser:
    ONE minimiser, one place.  The XYZ round-trip is cheap against L-BFGS.
    ``enable_global_pg=False`` because Tier D is the point group of the WHOLE molecule
    and this is only a cut-out fragment; Tier B/C (Morgan equivalence + graph orbits)
    stay on and are exactly what keeps a symmetric ligand's conformer symmetric.
    A conformer the refiner declines is left untouched -- never-worse per conformer.
    """
    try:
        from delfin.manta._variational_refiner import variational_refine
    except Exception:
        return
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    for cid in cids:
        try:
            conf = m.GetConformer(cid)
            pos = conf.GetPositions()
            xyz = f"{len(syms)}\nconf\n" + "\n".join(
                f"{s:4s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}"
                for s, p in zip(syms, pos))
            new_xyz, rep = variational_refine(xyz, m, class_label="no_metal",
                                              enable_global_pg=False)
            if rep.get("fallback_used", True):
                continue
            n_set = 0
            for ln in new_xyz.splitlines():
                parts = ln.split()
                if len(parts) == 4 and n_set < len(syms):
                    conf.SetAtomPosition(n_set, (float(parts[1]), float(parts[2]),
                                                 float(parts[3])))
                    n_set += 1
        except Exception:
            continue


def _ligand_confs_from_mol(frag_mol, k=10):
    """UNIVERSAL multi-conformer generation for a ligand (deterministic): K diverse
    ETKDG conformers (fixed seed, single-thread) + MMFF.  Returns (syms, [coords],
    mol).  Used to pick the clash-minimal conformer per ligand at placement — a
    fundamental Layer-3 mechanism applied to every ligand, not a per-case patch.

    Result is memoised by canonical SMILES + k (deterministic embed): repeated calls
    for the same ligand (e.g. the ensemble builder's variant loop) skip the costly
    re-embed.  The cached coords/mol are NOT mutated by any caller."""
    # ⚠️ ENTRY TRACE.  The first mesomery smoke test delivered SIX built systems
    # and ZERO trace lines -- the default sat behind the MMFF block and was never
    # reached.  That has three causes which cannot be separated without measuring:
    # the function is not called at all, it serves from the cache, or it returns
    # earlier.  My static proof "assemble_from_config calls it at
    # :4381" was TRUE and still not sufficient -- for the third time today the
    # confusion of "the site exists" with "this run arrives there".
    # This one line answers it, and it costs nothing when the path is empty.
    _entp = os.environ.get("DELFIN_MESO_TRACE", "")
    if _entp and _entp != "0":
        try:
            with open(_entp, "a") as _fh:
                _fh.write("[ENTER] _ligand_confs_from_mol natoms=%d k=%d\n"
                          % (frag_mol.GetNumAtoms(), int(k)))
        except Exception:
            pass
    key = None
    try:
        key = (Chem.MolToSmiles(frag_mol), int(k))
    except Exception:
        key = None
    if key is not None and key in _CONF_CACHE:
        return _CONF_CACHE[key]
    m = Chem.AddHs(frag_mol)
    try:
        cids = list(AllChem.EmbedMultipleConfs(m, numConfs=k, randomSeed=SEED,
                                               numThreads=1))
    except Exception:
        # Default (byte-identical): re-raise the original unguarded behaviour.
        if os.environ.get("DELFIN_FFFREE_KEKULIZE_SPLIT", "0") != "1":
            raise
        # Kekulize-robust retry (same root as the decompose split): a cleaved
        # aromatic-N⁺ ligand (triazolide / pyridinium / scorpionate cap) carries
        # ARTEFACT positive charges on ring N from the metal-dative SMILES encoding,
        # which leave the ring unkekulizable -> EmbedMultipleConfs raises and the
        # WHOLE complex falls to legacy.  Neutralise those aromatic-N⁺ formal charges
        # (geometry-only: connectivity + donor indices unchanged) so a full sanitize
        # kekulizes and the fragment embeds.  Verified: scorpionate cap embeds.
        cids = []
        try:
            fm = Chem.RWMol(frag_mol)
            for a in fm.GetAtoms():
                if (a.GetIsAromatic() and a.GetSymbol() == "N"
                        and a.GetFormalCharge() > 0):
                    a.SetFormalCharge(0)
            m = fm.GetMol()
            Chem.SanitizeMol(m)
            m = Chem.AddHs(m)
            cids = list(AllChem.EmbedMultipleConfs(m, numConfs=k, randomSeed=SEED,
                                                   numThreads=1))
        except Exception:
            cids = []
    if not cids:
        try:
            _emb_ok = AllChem.EmbedMolecule(m, randomSeed=SEED) == 0
        except Exception:
            if os.environ.get("DELFIN_FFFREE_KEKULIZE_SPLIT", "0") != "1":
                raise                  # byte-identical: original unguarded behaviour
            _emb_ok = False            # kekulize-robust: clean None instead of crash
        if not _emb_ok:
            if key is not None:
                _CONF_CACHE[key] = None
            return None
        cids = [0]
    # FF-FREE CONFORMER RELAX (DELFIN_FFREE_CONF_RELAX, default OFF -> byte-identical).
    # MMFF is the last real force field on the conformer axis, and it has NO metal
    # parameters -- it relaxes the conformers of a COORDINATED ligand with a model that
    # does not know the metal exists, the same defect class as the
    # "UFFTYPER: Unrecognized atom type: Pd+2" this build prints.  The fragment is cut
    # metal-free here, so the functional's existing "no_metal" preset fits exactly:
    # k_topology = 0, k_A = 0, and bond / signature-angle / torsion / clash / symmetry
    # carry the geometry.  ETKDG above is NOT touched: it is distance geometry, not a
    # force field, and it stays the generator.
    if _ffree_flag("CONF_RELAX"):
        _relax_confs_ffree(m, cids)
    else:
        try:
            AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=1)
        except Exception:
            pass
    _symmetrize_degenerate(m, cids)
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    out = (syms, [np.array(m.GetConformer(c).GetPositions(), float) for c in cids], m)
    if key is not None and len(_CONF_CACHE) < _CONF_CACHE_MAX:
        _CONF_CACHE[key] = out
    return out


def _symmetrize_degenerate(m, cids):
    """Set degenerate bond pairs to ONE length -- carboxylate, nitro, amidinate.

    ===== THE KEKULE DRAWING IS NOT A GEOMETRY ==================================
    A SMILES draws a carboxylate as ``C(=O)[O-]`` -- one double and one single
    bond.  In reality both C-O are equally long.  Measured 18.08. against the
    HIT distribution (only systems with ``org_bond_realized == false``):

        nitro,       worst bond N-O:  24.79 % against 1.94 %  = 12.75x
        carboxylate, worst bond C-O:  24.82 % against 7.00 %  =  3.55x

    And the direction is unambiguous: for nitro, **29 of 30 are TOO LONG**, mean
    deviation 0.233 Angstrom, of which **19 on the drawn SINGLE bond**.
    That is the Kekule signature, undisguised.  It also doubles
    ``pyramidal_sp2`` (factor 2.41) and thereby hits the three-times-measured
    pi root from the build side.
    Reach: 780 of 10000 systems hard-degenerate (7.8 %).

    ⚠ WHY HERE AND NOT IN A CORRECTOR.  This function runs in the METAL-FREE
    cut ligand fragment, BEFORE placement -- it is part of the construction of
    the ligand scaffold, not a post-correction on the finished complex.
    The obvious place ``refine._precompute_arom_targets`` would have reach
    near ZERO: its never-worse guard excludes everything that touches the
    coordination sphere -- and 281 of 301 carboxylates are on the metal.

    ⚠ WHY THE AROMATICS AXIS DOES NOT ALREADY DO IT.  ``AROM_SEAT`` demands RDKit
    aromaticity or a geometric 5/6-ring at all three seating sites.  A
    carboxylate is neither; the acac chelate ring would be a six-ring, but carries
    the metal and drops out.  Reach on these groups: exactly 0.

    THE RULE is element-agnostic, so that it is not pinned to patterns:
    for every heavy atom X its TERMINAL heavy neighbours are grouped by
    element; if a group has at least two members AND
    different bond orders, it is degenerate and all its bonds
    get the MEAN length.  That hits carboxylate, nitro, nitrate, sulfonate,
    phosphonate and amidinate, without any of them appearing in the code.  If the
    orders are already equal, nothing is touched -- then the drawer has already
    expressed the symmetry.

    DELFIN_FFFREE_MESOMERY_SEAT (default 0 -> byte-identical).
    """
    if os.environ.get("DELFIN_FFFREE_MESOMERY_SEAT", "0") != "1":
        return
    # ⚠️ THE TRACE MUST STAND BEFORE THE try.  In the first draft it sat inside -- an
    # exception during the group search would thereby have landed in `except: pass`
    # WITHOUT writing a line, and could not have been distinguished from "no group
    # found".  A silent exception path is exactly the gap at which the fire census
    # invented findings on 14.08.
    _mtp = os.environ.get("DELFIN_MESO_TRACE", "")

    def _mtrace(txt):
        if not _mtp or _mtp == "0":
            return
        try:
            with open(_mtp, "a") as _fh:
                _fh.write("[MESO] %s\n" % txt)
        except Exception:
            pass

    try:
        groups = []
        for a in m.GetAtoms():
            if a.GetSymbol() == "H":
                continue
            by_el = {}
            for b in a.GetBonds():
                nb = b.GetOtherAtom(a)
                if nb.GetSymbol() == "H":
                    continue
                # terminal: no heavy neighbour other than X
                if sum(1 for x in nb.GetNeighbors() if x.GetSymbol() != "H") != 1:
                    continue
                by_el.setdefault(nb.GetSymbol(), []).append(
                    (nb.GetIdx(), float(b.GetBondTypeAsDouble())))
            for _el, mem in by_el.items():
                if len(mem) < 2:
                    continue
                if len({round(o, 2) for _i, o in mem}) < 2:
                    continue          # already drawn symmetrically -> nothing to do
                groups.append((a.GetIdx(), [i for i, _o in mem]))
        # THE SAME TRACE DISCIPLINE AS IN THE ASSERTION PROTOCOL: a byte-identical
        # build here too has several indistinguishable causes -- the function does not
        # run, it finds no group, it aborts, or it finds one and the
        # bonds are already equally long.
        if not groups:
            _mtrace("natoms=%d groups=0 moved=0.0" % m.GetNumAtoms())
            return
        _worst = 0.0
        for c in cids:
            conf = m.GetConformer(c)
            for cen, terms in groups:
                pc = np.array(conf.GetAtomPosition(cen), float)
                vecs, lens = [], []
                for t in terms:
                    v = np.array(conf.GetAtomPosition(t), float) - pc
                    n = float(np.linalg.norm(v))
                    if n < 1e-6:
                        vecs, lens = [], []
                        break
                    vecs.append(v / n)
                    lens.append(n)
                if not lens:
                    continue
                tgt = float(sum(lens) / len(lens))
                for t, v, l0 in zip(terms, vecs, lens):
                    _worst = max(_worst, abs(tgt - l0))
                    p = pc + v * tgt
                    conf.SetAtomPosition(t, (float(p[0]), float(p[1]), float(p[2])))
        _mtrace("natoms=%d groups=%d conf=%d worst_shift=%.4f"
                % (m.GetNumAtoms(), len(groups), len(cids), _worst))
    except Exception as _mx:
        # a default that does not take effect must cost nothing -- but it must not
        # stay silent either, otherwise the abort reads like "nothing found".
        _mtrace("ABBRUCH %s" % type(_mx).__name__)


def _clash_count(Q, existing, syms_Q, syms_ex):
    """# heavy/H pairs between block Q and existing atoms closer than 0.7*(vdW sum)."""
    if len(existing) == 0:
        return 0
    from delfin.manta.refine import _vdw
    c = 0
    for a in range(len(Q)):
        for b in range(len(existing)):
            d = float(np.linalg.norm(Q[a] - existing[b]))
            if d < 0.70 * (_vdw(syms_Q[a]) + _vdw(syms_ex[b])):
                c += 1
    return c


def _ligand_block_bonds(lmol, offset, donor_local):
    """True connectivity of one ligand block for the #308 torsion relaxer.

    Returns ``(offset, [(li, lj), ...], donor_local)`` where the local (i,j) bond
    pairs come directly from the ligand mol (atom order preserved through assembly,
    AddHs included), so the relaxer never has to GUESS bonds from distance on a
    crowded complex (where two ligands at a fortuitous bonding distance would be
    mis-read as covalently bonded).  Returns ``None`` on any failure (relaxer then
    falls back to geometric perception)."""
    try:
        lb = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in lmol.GetBonds()]
        return (int(offset), lb, int(donor_local))
    except Exception:
        return None


def _torsion_relax_frame(out_syms, P, fixed, block_specs):
    """Apply the env-gated #308 torsion-space clash relax to one assembled frame,
    threading the true per-ligand connectivity (``block_specs`` = list of
    ``(offset, lmol, donor_local)``).  No-op when the flag is unset; never raises."""
    try:
        from delfin.manta import torsion_relax as _TR
        bp = None
        if block_specs:
            blocks = [bb for bb in (_ligand_block_bonds(m, off, dl)
                                    for (off, m, dl) in block_specs) if bb is not None]
            if blocks:
                bp = _TR.bonds_from_blocks(0, blocks)
        return np.asarray(_TR.relax_if_enabled(out_syms, P, fixed, bond_pairs=bp),
                          dtype=float)
    except Exception:
        return P


def _joint_declash_frame(out_syms, P, fixed, block_specs, geom=None):
    """Apply the env-gated JOINT global INTER-LIGAND heavy-heavy declash to one
    assembled frame (``DELFIN_FFFREE_JOINT_DECLASH``), threading the true
    per-ligand connectivity (``block_specs`` = list of ``(offset, lmol,
    donor_local)``).  Runs AFTER #308 torsion-relax and BEFORE the self-gate so a
    declashed class-B build passes ``_build_is_clean``.  No-op when the flag is
    unset; never raises."""
    try:
        from delfin.manta import joint_declash as _JD
        bp = None
        if block_specs:
            blocks = [bb for bb in (_ligand_block_bonds(m, off, dl)
                                    for (off, m, dl) in block_specs) if bb is not None]
            if blocks:
                bp = _JD._TR.bonds_from_blocks(0, blocks)
        return np.asarray(_JD.declash_if_enabled(out_syms, P, fixed, geom=geom, bond_pairs=bp),
                          dtype=float)
    except Exception:
        return P


def _refine_guarded(out_syms, P, fixed):
    """``refine()`` with the seating's assertion before and after it.

    ⚠️ WHY THIS FUNCTION EXISTS -- a mistake of mine, measured on 18.08. and
    recorded here so that nobody repeats it.  I had first written the guard
    INLINE at ONE call site (``assemble_heteroleptic_from_mols``) and
    then checked it on the 19 measured OC-6 cases: 19 of 19 byte-identical.  That
    looked like "the assertion holds".  The positive control refuted it: with
    tolerance 0.0001 Angstrom -- where EVERY relaxation must trigger -- the build
    stayed identical as well.  So the block never ran.  The FF-free chelate builder goes
    through ``assemble_from_config`` -> ``_finish_config_frame``, a DIFFERENT function
    with its OWN refine call site.
    I had checked the line for reachability and not the FUNCTION -- the same
    construction as ``ISOLATED_SEAT``, which I documented in others on the same day.
    ⇒ The guard belongs at ALL four refine call sites, hence in ONE function.

    Default OFF -> byte-identical: without the switch this is exactly the old
    ``try: P = refine(...) except: pass`` block.
    """
    _assert_on = os.environ.get("DELFIN_FFFREE_ASSERT_ENFORCE", "0") == "1"
    _assertion, _P_before = None, None
    if _assert_on:
        try:
            from delfin.manta import _frame_assertions as _FA
            # The BUILDER's set, not my reconstruction of it: refine() receives
            # `fixed` as a promise, so exactly that is the contract it has to keep.
            _assertion = _FA.derive((list(out_syms), P), frozen=fixed)
            _P_before = P.copy()
        except Exception:
            _assertion = _P_before = None
    try:
        from delfin.manta.refine import refine as _refine
        P = _refine(out_syms, P, fixed)
    except Exception:
        pass
    if _assertion is not None and _P_before is not None:
        try:
            from delfin.manta import _frame_assertions as _FA
            _v = _FA.violations(_assertion, (list(out_syms), P))
            # ⚠️ A TRACE, BECAUSE A BYTE COMPARISON DOES NOT DECIDE HERE.
            # The first smoke test showed "identical" -- and that has three causes
            # indistinguishable to the naked eye: (a) the block does not run, (b)
            # derive() returns None, (c) refine() moves nothing, in which case the
            # rollback is a no-op.  Exactly this confusion let the fire census
            # invent findings on 14.08.  The trace separates them:
            #   derived=1 rules out (b), moved=... rules out (c), broke=1 is the hit.
            # DELFIN_ASSERT_TRACE=<path>, otherwise silent and free.
            _tp = os.environ.get("DELFIN_ASSERT_TRACE", "")
            if _tp and _tp != "0":
                try:
                    _mv = float(np.max(np.linalg.norm(P - _P_before, axis=1)))
                except Exception:
                    _mv = -1.0
                try:
                    with open(_tp, "a") as _fh:
                        _fh.write("[ASSERT] n=%d derived=1 moved=%.4f broke=%d %s\n"
                                  % (len(out_syms), _mv,
                                     1 if (_v and _v.get("any_broken")) else 0,
                                     "" if not _v else
                                     "md=%d planar=%d frozen=%d trans=%d" % (
                                         _v.get("md_broken", 0), _v.get("planar_broken", 0),
                                         _v.get("frozen_moved", 0),
                                         _v.get("trans_lost_metals", 0))))
                except Exception:
                    pass
            if _v is not None and _v.get("any_broken"):
                P = _P_before              # rollback: the claim weighs more
                # ===== DIAGNOSIS ONLY, NEVER IN PRODUCTION =======================
                # Question this answers: the full run over the 19
                # octahedral cases broke the assertion SIX TIMES (109 calls, all
                # six trans) -- and still delivered 19 of 19 byte-identical
                # archive files.  Two explanations, both plausible:
                #   (a) the affected frames are rejected by the self-gate anyway, they
                #       stand in NO arm in the archive;
                #   (b) a pass AFTER refine restores the twisted state.
                # Both produce byte-identical archives, so a byte comparison cannot
                # separate them -- the same trap as in the first smoke test.
                # DELFIN_ASSERT_MARK_ROLLBACK=1 additionally shifts the rolled-back frame
                # by 100 Angstrom.  Such a frame is geometrically
                # impossible and falls through every gate.  If the archive is THEN still
                # byte-identical, the frame was never in it -> case (a).
                # If it changes, the rollback does very well reach the archive ->
                # case (b), and the cause lies behind refine.
                if os.environ.get("DELFIN_ASSERT_MARK_ROLLBACK", "0") == "1":
                    P = P + np.array([100.0, 0.0, 0.0], float)
        except Exception:
            pass
    return P


def _sphere_flex_frame(out_syms, P, fixed, block_specs):
    """Apply the env-gated soft coordination-sphere clash relax to one assembled
    frame (``DELFIN_FFFREE_SPHERE_FLEX``).  Donors are soft-restrained (not frozen)
    so the sphere can BREATHE a few hundredths of an Angstrom to open the residual
    inter-ligand heavy-heavy contacts that a frozen-donor refine (and pure M-D-axis
    rotation) cannot.  Threads the true connectivity; runs AFTER joint-declash and
    BEFORE the self-gate.  No-op when the flag is unset; never raises."""
    try:
        from delfin.manta import sphere_flex as _SF
        from delfin.manta import torsion_relax as _TR
        bp = None
        if block_specs:
            blocks = [bb for bb in (_ligand_block_bonds(m, off, dl)
                                    for (off, m, dl) in block_specs) if bb is not None]
            if blocks:
                bp = _TR.bonds_from_blocks(0, blocks)
        return np.asarray(_SF.flex_if_enabled(out_syms, P, fixed, bond_pairs=bp),
                          dtype=float)
    except Exception:
        return P
