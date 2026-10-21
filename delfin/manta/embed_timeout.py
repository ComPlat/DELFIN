"""Timeout-guarded RDKit embedding (single and multi-conformer) of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import threading
from typing import List, Optional

from delfin.common.logging import get_logger
from delfin.manta.converter_flags import (
    _DEFAULT_EMBED_SEED,
    _EMBED_TIMEOUT,
    _MULTIEMBED_TIMEOUT,
    _MULTIEMBED_TIMEOUT_OVERRIDE,
    _deterministic_mode,
)
from delfin.manta.ml_tables import (
    AllChem,
)

logger = get_logger("delfin.smiles_converter")


def _embed_with_timeout(mol, params=None, timeout: Optional[float] = None):
    """Run AllChem.EmbedMolecule with a timeout guard.

    Returns the embed result (0 on success, -1 on failure/timeout).
    The *mol* is modified in-place on success (same as EmbedMolecule).

    When *params* is ``None`` a fixed-seed :class:`EmbedParameters` is
    substituted so the embed stays deterministic; a bare
    ``AllChem.EmbedMolecule(mol)`` would otherwise use RDKit's wall-clock
    default seed (``randomSeed = -1``).
    """
    if timeout is None:
        timeout = _EMBED_TIMEOUT

    # Determinism guard: never let a caller fall through to RDKit's
    # wall-clock-seeded default embed.
    if params is None:
        params = AllChem.EmbedParameters()
        params.randomSeed = _DEFAULT_EMBED_SEED

    # Fast path: skip timeout machinery for very small molecules
    # BUT always use timeout for highly connected molecules (many ring bonds)
    # which can cause ETKDG to hang even at low atom counts (e.g. borane cages).
    n_atoms = mol.GetNumAtoms()
    n_bonds = mol.GetNumBonds()
    highly_connected = n_bonds > 2 * n_atoms  # cage/cluster topology
    if n_atoms < 40 and not highly_connected:
        if params is not None:
            return AllChem.EmbedMolecule(mol, params)
        return AllChem.EmbedMolecule(mol)

    # Scale timeout with mol size: large hexadentate / fused-ring ligands
    # (Cu-salen-biphep, Fe-terpy-NMe2, Pt-phos-terpy) need much more than
    # the default 10 s to finish a single ETKDG embedding pass.  Without
    # this scaling, every attempt times out, the background thread is
    # orphaned (still running ETKDG in C++), the pipeline retries with
    # a new seed, and orphan threads accumulate until the whole
    # subprocess hits the caller's wall-clock timeout.
    if n_atoms > 50:
        timeout = max(timeout, 10.0 + 0.5 * (n_atoms - 50))

    result = [-1]
    exc_holder = [None]

    def _do_embed():
        try:
            if params is not None:
                result[0] = AllChem.EmbedMolecule(mol, params)
            else:
                result[0] = AllChem.EmbedMolecule(mol)
        except Exception as e:
            exc_holder[0] = e

    t = threading.Thread(target=_do_embed, daemon=True)
    t.start()
    # Determinism: under the master switch, wait for the embed to finish
    # (no wall-clock cutoff) so the result never depends on timing/CPU load.
    t.join(timeout=None if _deterministic_mode() else timeout)

    if t.is_alive():
        logger.debug(
            "EmbedMolecule timed out after %.1fs for mol with %d atoms",
            timeout, mol.GetNumAtoms(),
        )
        # Thread cannot be killed; it will finish eventually in the background.
        # Return failure so the caller moves on to the next strategy.
        return -1

    if exc_holder[0] is not None:
        raise exc_holder[0]

    return result[0]


def _embed_multiple_confs_with_timeout(
    mol,
    num_confs: int,
    params,
    timeout: Optional[float] = None,
) -> List[int]:
    """Run AllChem.EmbedMultipleConfs with a hard wall-clock timeout.

    The embed runs in a daemon thread; if it hasn't finished within *timeout*
    seconds we return whatever conformers have already been added to *mol*.
    RDKit continues executing in the background (threads are not killable),
    but the caller is unblocked immediately.
    """
    if timeout is None:
        # Honour the multi-sigma V2 thread-local override (default OFF).
        # Operator-set ``DELFIN_MULTIEMBED_TIMEOUT`` still applies via the
        # module-level constant when no override is active.
        try:
            _ovr = getattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value", None)
        except Exception:
            _ovr = None
        timeout = float(_ovr) if _ovr is not None else _MULTIEMBED_TIMEOUT

    # Scale timeout + num_confs with mol size.  For huge fused-ring
    # ligands (>80 atoms) the default 25 s default cannot finish even
    # a single conformer; shrinking num_confs while extending timeout
    # lets one attempt finish and avoids the orphan-thread pile-up.
    n_atoms = mol.GetNumAtoms()
    if n_atoms > 80:
        timeout = max(timeout, 25.0 + 1.0 * (n_atoms - 80))
        num_confs = min(num_confs, 3)
    elif n_atoms > 50:
        timeout = max(timeout, 25.0 + 0.5 * (n_atoms - 50))
        num_confs = min(num_confs, 5)

    before = {int(c.GetId()) for c in mol.GetConformers()}
    result: List[int] = []
    exc_holder: List[Optional[BaseException]] = [None]

    def _do_embed():
        try:
            ids = list(
                AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
            )
            result.extend(ids)
        except BaseException as e:
            exc_holder[0] = e

    t = threading.Thread(target=_do_embed, daemon=True)
    t.start()
    # Determinism: under the master switch, wait for the embed to finish
    # (no wall-clock cutoff) so the result never depends on timing/CPU load.
    t.join(timeout=None if _deterministic_mode() else timeout)

    if t.is_alive():
        logger.debug(
            "EmbedMultipleConfs timed out after %.1fs (n_atoms=%d, n_bonds=%d); "
            "returning partial conformer set.",
            timeout, mol.GetNumAtoms(), mol.GetNumBonds(),
        )
        after = {int(c.GetId()) for c in mol.GetConformers()}
        return sorted(after - before)

    if exc_holder[0] is not None:
        raise exc_holder[0]

    return list(result)


def _make_random_embed_params(seed: int = 42):
    """Create minimal embed params that skip ETKDG torsion/knowledge checks.

    These are much faster for complex metal ring systems where the ETKDG
    distance bounds matrix computation can hang.
    """
    params = AllChem.EmbedParameters()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.useBasicKnowledge = False
    params.useExpTorsionAnglePrefs = False
    params.enforceChirality = False
    return params
