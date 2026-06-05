"""delfin.fffree.embed_fallback -- Universal Embed + FF-free polish fallback.

Mission F2 (User 2026-06-05): when fffree's constructive path can't build a
SMILES (unknown class, exotic geometry, organic, ...), instead of falling
through to the legacy UFF pipeline we run a deterministic RDKit ETKDGv3
embed and -- if the molecule contains a recognised metal centre -- finish
with the SAME GRIP polish the constructive path uses.

The result is the same ``[(xyz_string, label), ...]`` shape produced by
``_fffree_isomers`` so the converter dispatcher is byte-identical for the
*caller* regardless of which path generated the coordinates.

Determinism contract
--------------------
* ``ETKDGv3`` with fixed ``randomSeed = 42``
* ``useRandomCoords = False``
* No threading: ``numThreads = 1``
* Two consecutive calls with the same SMILES + same env return
  byte-identical XYZ blocks (verified by the F2 test suite).

Selectable fallback mode
------------------------
The dispatcher in :mod:`delfin.smiles_converter` reads
``DELFIN_FFFREE_FALLBACK_MODE`` (see :func:`resolve_fallback_mode` for the
accepted values) and routes accordingly.  The polish step here is itself
gated by the value: ``grip`` polishes coordinates, ``uff`` skips polish
(legacy semantics), ``none`` is never reached (the dispatcher returns ``[]``
before invoking us), ``both`` returns both variants.

This module never raises -- failure modes (RDKit unavailable, embed fails,
polish fails) all degrade silently to a sensible result so the converter
caller is not forced to handle exotic errors.
"""
from __future__ import annotations

import os
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

# Deterministic ETKDG seed (F2 contract -- never change without bumping the
# byte-identity proof).  Each conformer index uses a +k offset of this so the
# user can extend the ensemble without disturbing earlier conformers.
ETKDG_SEED = 42

# Valid env values for DELFIN_FFFREE_FALLBACK_MODE.
_VALID_MODES = ("grip", "uff", "none", "both")


def resolve_fallback_mode(env: Optional[dict] = None) -> str:
    """Resolve the active fallback mode from the environment.

    Parameters
    ----------
    env : mapping, optional
        Pass an explicit dict for tests; defaults to ``os.environ``.

    Returns
    -------
    str
        One of ``"grip"``, ``"uff"``, ``"none"``, ``"both"``.  Unknown
        values silently coerce to ``"uff"`` (the default-OFF byte-identical
        legacy behaviour); unset = ``"uff"`` as well.
    """
    src = env if env is not None else os.environ
    raw = str(src.get("DELFIN_FFFREE_FALLBACK_MODE", "")).strip().lower()
    if raw in _VALID_MODES:
        return raw
    return "uff"


def _xyz_block(syms: Sequence[str], P: np.ndarray) -> str:
    """Canonical headerless XYZ atom block (matches converter_backend._xyz)."""
    return "\n".join(
        f"{s:4s} {float(x):12.6f} {float(y):12.6f} {float(z):12.6f}"
        for s, (x, y, z) in zip(syms, P)
    )


def _detect_metal_donors(mol) -> Tuple[Optional[int], List[int]]:
    """Return (metal_atom_idx, donor_atom_idxs) for the polish call.

    Falls back to ``(None, [])`` when no recognised metal is present (which
    means the polish step is skipped and we just emit the raw embed).
    """
    try:
        from delfin._bond_decollapse import _is_metal
    except Exception:
        return None, []
    metal_idx: Optional[int] = None
    for atom in mol.GetAtoms():
        if _is_metal(atom.GetSymbol()):
            metal_idx = atom.GetIdx()
            break
    if metal_idx is None:
        return None, []
    donors = [
        nbr.GetIdx()
        for nbr in mol.GetAtomWithIdx(metal_idx).GetNeighbors()
        if not _is_metal(nbr.GetSymbol())
    ]
    return metal_idx, donors


def _maybe_grip_polish(
    coords: np.ndarray,
    mol,
    *,
    geom_hint: str = "",
) -> Tuple[np.ndarray, str]:
    """Optionally run grip_polish on ``coords``.

    The polish is silently skipped when:
      * no metal atom is present (organic SMILES -- nothing to polish around)
      * grip_polish import / call fails
      * the polish would mutate frozen / donor atoms in a non-finite way

    Returns
    -------
    (np.ndarray, str)
        The (possibly polished) coordinates and a short status tag for the
        label suffix (``"grip"`` if polish accepted, ``"raw"`` otherwise).
    """
    metal_idx, donors = _detect_metal_donors(mol)
    if metal_idx is None or not donors:
        return coords, "raw"
    try:
        from delfin.fffree.grip_polish import grip_polish
    except Exception:
        return coords, "raw"
    try:
        P_polished = grip_polish(
            coords,
            mol,
            metal=int(metal_idx),
            donors=list(donors),
            geom=geom_hint,
        )
        P_polished = np.asarray(P_polished, dtype=float)
        if P_polished.shape == coords.shape and np.all(np.isfinite(P_polished)):
            return P_polished, "grip"
    except Exception:
        pass
    return coords, "raw"


def embed_isomers(
    smiles: str,
    *,
    max_isomers: int = 16,
    polish: str = "grip",
) -> Optional[List[Tuple[str, str]]]:
    """Generate ETKDGv3 embed isomers + (optionally) GRIP polish.

    Parameters
    ----------
    smiles : str
        Input SMILES.
    max_isomers : int
        Upper bound on conformers; clamped to ``[1, 16]`` so very large
        molecules don't explode the wall-clock.  The polish step is run
        per accepted embed.
    polish : str
        ``"grip"`` (default): FF-free GRIP polish; ``"raw"``: skip polish
        and emit the bare ETKDG coordinates; ``"both"``: emit both variants
        per embed (label suffix distinguishes them).  Any other value falls
        back to ``"raw"`` semantics.

    Returns
    -------
    list of (xyz, label) tuples, or ``None`` if RDKit cannot parse the
    SMILES or the embed produced no conformers.  Never raises.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = Chem.AddHs(mol)
    except Exception:
        return None

    n_embed = max(1, min(int(max_isomers), 16))

    # Deterministic ETKDGv3 parameter block.  ``useRandomCoords = False`` is
    # essential for the byte-identity contract; ``numThreads = 1`` prevents
    # non-deterministic thread scheduling on multi-core boxes.
    params = AllChem.ETKDGv3()
    params.randomSeed = ETKDG_SEED
    params.useRandomCoords = False
    params.numThreads = 1
    try:
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_embed, params=params)
    except Exception:
        return None
    cids = list(cids)
    if not cids:
        return None

    syms = [a.GetSymbol() for a in mol.GetAtoms()]

    polish_mode = polish if polish in ("grip", "raw", "both") else "raw"
    out: List[Tuple[str, str]] = []
    for cid in cids:
        try:
            conf = mol.GetConformer(cid)
            coords = np.asarray(conf.GetPositions(), dtype=float)
        except Exception:
            continue
        if coords.size == 0 or not np.all(np.isfinite(coords)):
            continue
        if polish_mode == "raw":
            out.append((_xyz_block(syms, coords), f"embed-conf{cid}-raw"))
            continue
        # grip and both: try polish; "both" additionally emits the raw variant.
        if polish_mode == "both":
            out.append((_xyz_block(syms, coords), f"embed-conf{cid}-raw"))
        P_out, status = _maybe_grip_polish(coords, mol)
        suffix = "grip" if status == "grip" else "raw"
        out.append((_xyz_block(syms, P_out), f"embed-conf{cid}-{suffix}"))

    return out or None


def grip_embed_fallback(
    smiles: str,
    *,
    max_isomers: int = 16,
) -> Optional[List[Tuple[str, str]]]:
    """Convenience: ETKDG + GRIP polish (the production F2 default)."""
    return embed_isomers(smiles, max_isomers=max_isomers, polish="grip")


def raw_embed_fallback(
    smiles: str,
    *,
    max_isomers: int = 16,
) -> Optional[List[Tuple[str, str]]]:
    """Convenience: ETKDG only -- no polish (forensic A/B counterpart)."""
    return embed_isomers(smiles, max_isomers=max_isomers, polish="raw")


def both_embed_fallback(
    smiles: str,
    *,
    max_isomers: int = 16,
) -> Optional[List[Tuple[str, str]]]:
    """Convenience: emit both raw and grip-polished variants per conformer."""
    return embed_isomers(smiles, max_isomers=max_isomers, polish="both")
