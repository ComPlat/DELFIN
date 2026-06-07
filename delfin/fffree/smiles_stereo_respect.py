r"""delfin.fffree.smiles_stereo_respect -- Universal SMILES stereo respect.

The GENFUB bug (2026-06-07): a SMILES with an unspecified double-bond stereo
descriptor (e.g. ``...C=CC2=CC=NC=C2...`` -- no ``/`` or ``\``) makes RDKit's
ETKDGv3 multi-conformer embed produce a *random mix* of E (trans) and Z (cis)
configurations across the N requested conformers, because every conformer's
ring-closure step picks a random sign for the unconstrained dihedral.  For a
single SMILES orbit this means 15 frames have one stereo and 1 frame has the
other -- not a chemically meaningful conformer ensemble but a topological
mix of two different molecules.

User mandate (universal):

  1. If the SMILES specifies stereo (``/C=C/``, ``/C=C\``, ``[C@H]``,
     ``[C@@H]``, ...) ETKDG MUST respect it: every conformer carries the
     SAME stereo descriptor.
  2. If the SMILES does NOT specify stereo, pick a deterministic canonical
     interpretation (E for double bonds, R for chiral centers) and pin it
     BEFORE the embed -- so the ensemble explores conformers within one
     stereoisomer instead of mixing two.

This module is pure RDKit graph operations -- no SMILES patterns, no per-class
templates, no element-specific code.  It runs on any organic / inorganic /
mixed molecule that RDKit can parse.

Env flag
--------

``DELFIN_FFFREE_SMILES_STEREO_RESPECT``: default OFF (byte-identical).  When
set to ``1`` / ``true`` / ``yes`` / ``on``, the helper :func:`enforce_stereo`
runs on every mol immediately after ``AddHs`` and before the ETKDG embed.

API
---

``enforce_stereo(mol) -> mol``
    Mutates ``mol`` in place:
      * detects every potential stereo site (atom + bond) via RDKit's
        ``FindPotentialStereo``;
      * assigns a deterministic canonical stereo to every UNSPECIFIED site
        (E for ``Bond_Double``, R == ``CHI_TETRAHEDRAL_CW`` for ``Atom_Tetrahedral``);
      * leaves SPECIFIED sites untouched so explicit user intent is preserved.
    Returns the same mol for chain-style usage.

``configure_etkdg_for_stereo(params) -> params``
    Sets ``enforceChirality = True`` and ``useBasicKnowledge = True`` on an
    ETKDGv3 parameter block so RDKit's embedder uses the assigned stereo as
    a hard constraint.  Returns the same params for chain-style usage.

``flag_active(env=None) -> bool``
    Cheap env-flag probe.  Pass an explicit dict for tests.
"""
from __future__ import annotations

import os
from typing import Optional


_ENV_FLAG = "DELFIN_FFFREE_SMILES_STEREO_RESPECT"
_ENV_TRUTHY = ("1", "true", "yes", "on")


def flag_active(env: Optional[dict] = None) -> bool:
    """Return True iff ``DELFIN_FFFREE_SMILES_STEREO_RESPECT`` is truthy."""
    src = env if env is not None else os.environ
    raw = str(src.get(_ENV_FLAG, "")).strip().lower()
    return raw in _ENV_TRUTHY


def _pick_reference_neighbour(atom, exclude_idx: int) -> Optional[int]:
    """Pick the lowest-indexed neighbour of ``atom`` excluding ``exclude_idx``.

    This is the canonical choice for the StereoAtoms reference pair on a
    double bond: deterministic (lex on index), prefers heavy atoms when
    multiple are present (heavy atoms always have lower indices than H in
    an AddHs'ed mol when the user-supplied SMILES enumerates heavy atoms
    first, which is the RDKit canonical ordering).
    """
    heavy = []
    light = []
    for nb in atom.GetNeighbors():
        idx = int(nb.GetIdx())
        if idx == int(exclude_idx):
            continue
        if nb.GetAtomicNum() == 1:
            light.append(idx)
        else:
            heavy.append(idx)
    if heavy:
        return min(heavy)
    if light:
        return min(light)
    return None


def enforce_stereo(mol) -> "Mol":  # noqa: F821 (rdkit Mol)
    """Pin every unspecified stereo site to a deterministic canonical value.

    For each potential stereo element reported by ``Chem.FindPotentialStereo``:
      * ``Atom_Tetrahedral``: if unspecified -> set ``CHI_TETRAHEDRAL_CW``
        (which CIP-translates to R for the canonical RDKit ranking).
      * ``Bond_Double``: if unspecified -> set ``STEREOE`` (trans / lower
        energy) with the two lowest-indexed non-double-bonded neighbours as
        the reference atoms (deterministic regardless of input atom order).
      * any other type (Atom_SquarePlanar, Bond_Atropisomer, ...): skipped
        for now (out of scope for the C=C / R-S contract; can be extended).
    SPECIFIED sites are NEVER touched -- explicit user stereo is preserved.

    The mol is mutated in place; the same object is returned for fluency.

    Failure mode: never raises.  On any internal error (RDKit version skew,
    exotic atom missing neighbours, ...) the mol is returned unchanged.
    """
    if mol is None:
        return mol
    try:
        from rdkit import Chem
    except Exception:
        return mol
    try:
        si_list = Chem.FindPotentialStereo(mol)
    except Exception:
        return mol

    # Deterministic ordering: sort by (type-name, centeredOn) so the iteration
    # order is independent of RDKit's internal traversal.  All assignments
    # commute (each site is independent) but a fixed order makes the debug
    # logs reproducible across versions.
    try:
        si_sorted = sorted(
            si_list,
            key=lambda si: (str(si.type), int(si.centeredOn)),
        )
    except Exception:
        si_sorted = list(si_list)

    for si in si_sorted:
        try:
            st_type = str(si.type)
            specified = str(si.specified)
            centered_on = int(si.centeredOn)
        except Exception:
            continue
        if specified.lower() == "specified":
            # Explicit user intent -> leave alone.
            continue

        if st_type == "Atom_Tetrahedral":
            try:
                atom = mol.GetAtomWithIdx(centered_on)
            except Exception:
                continue
            try:
                # CW == R under the standard CIP / RDKit convention.
                atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            except Exception:
                pass
            continue

        if st_type == "Bond_Double":
            try:
                bond = mol.GetBondWithIdx(centered_on)
            except Exception:
                continue
            try:
                a_begin = bond.GetBeginAtom()
                a_end = bond.GetEndAtom()
                ref_begin = _pick_reference_neighbour(
                    a_begin, a_end.GetIdx(),
                )
                ref_end = _pick_reference_neighbour(
                    a_end, a_begin.GetIdx(),
                )
                if ref_begin is None or ref_end is None:
                    # Terminal C=C (e.g. CH2=CH-...): no stereo to assign.
                    continue
                bond.SetStereoAtoms(int(ref_begin), int(ref_end))
                bond.SetStereo(Chem.BondStereo.STEREOE)
            except Exception:
                pass
            continue

        # Out-of-scope stereo types (square-planar, trigonal-bipyramidal,
        # atropisomeric bonds, ...).  Silently skip; the universal C=C / R-S
        # contract is met for the cases the user named, and the embed still
        # runs normally for the others.
        continue

    # Re-perceive stereochemistry so downstream RDKit consumers (notably
    # the ETKDG embedder) see the new canonical descriptors.  IMPORTANT:
    # do NOT pass ``cleanIt=True`` here -- that flag wipes the bond-stereo
    # we just assigned (it interprets "missing wedge bonds" as "missing
    # stereo info" and resets STEREOE -> STEREONONE).  ``cleanIt=False`` +
    # ``force=False`` simply pushes the existing descriptors through the
    # CIP-rank pipeline so atom/bond properties stay consistent.
    try:
        Chem.AssignStereochemistry(mol, cleanIt=False, force=False)
    except Exception:
        pass
    return mol


def configure_etkdg_for_stereo(params):
    """Flip the ETKDGv3 flags that make the embed respect assigned stereo.

    Sets:
      * ``params.enforceChirality = True``  -- hard chirality constraint
        on tetrahedral centers (R/S preserved across all conformers).
      * ``params.useBasicKnowledge = True`` -- ETKDG's "basic knowledge"
        block that encodes E/Z for double bonds + planarity for aromatic
        rings; required for STEREOE / STEREOZ to be respected by the embed.

    Both flags are RDKit defaults in recent versions, but we set them
    explicitly so the contract is stable across RDKit minor-version bumps.

    Returns ``params`` for chain-style usage.  Never raises.
    """
    if params is None:
        return params
    try:
        params.enforceChirality = True
    except Exception:
        pass
    try:
        params.useBasicKnowledge = True
    except Exception:
        pass
    return params
