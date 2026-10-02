"""The embed wall-clock join follows the determinism master switch.

The mission brief's finding: Ir(ppy)2(acac) TP-6-conf2, Cd-histidine CN=7 and
D-AQIWAZ build differently on the cluster vs. the GitHub slow-tests runner
(67fc8f9e marked all three non-strict). The mechanism this file pins:

  _embed_with_timeout and _embed_multiple_confs_with_timeout run the RDKit
  embed in a daemon thread and join it with a WALL-CLOCK timeout
  (delfin/manta/embed_timeout.py:84-91 and :185-190). On timeout the caller
  rotates seeds / takes a partial conformer set, so the emitted geometry —
  and through it the topology-gate verdict — depends on how fast the machine
  was at that moment. Under DELFIN_DETERMINISTIC=1 the join is unbounded
  (converter_flags.py:132-142), removing exactly that load dependence.

These tests assert that the switch actually governs both joins, by forcing a
50 ms timeout on a quaterphenyl (24 heavy C, 42 atoms with H: above the
n_atoms>=40 threading threshold at embed_timeout.py:73-75, below the >50-atom
timeout-scaling floor at :81-83 and :127-133 — the raw probe timeout only
survives in that 40..50 window). Legacy mode must return the timed-out
result (-1 / partial set); deterministic mode must wait and return the full
embed. If the switch ever stops governing a join, a loaded machine will
silently drop heavy conformers again — the environment-dependent build class
this branch was sent out to diagnose.

Molecule choice is asserted, not assumed: 40..50 atoms, checked before the
probe runs.
"""
from __future__ import annotations

import os

import pytest


def _clear_delfin_env():
    saved = {k: v for k, v in os.environ.items() if k.startswith("DELFIN_")}
    for k in list(saved):
        del os.environ[k]
    return saved


def _restore(saved):
    for k, v in saved.items():
        os.environ[k] = v
    for k in [k for k in os.environ if k.startswith("DELFIN_")
              and k not in saved]:
        del os.environ[k]


# Quaterphenyl: 24 heavy C, 42 atoms with H — the 40..50-atom window where
# the join runs with the RAW caller timeout (no scaling floor) but the
# threading path is taken.
_QUATERPHENYL = "c1ccccc1-c1ccccc1-c1ccccc1-c1ccccc1"


def _probe_mol():
    from rdkit import Chem

    mol = Chem.AddHs(Chem.MolFromSmiles(_QUATERPHENYL))
    assert mol is not None
    assert 40 <= mol.GetNumAtoms() <= 50, (
        "probe needs 40..50 atoms: >=40 for the threading path, <=50 to "
        "stay under the timeout-scaling floor that would override the "
        "forced tiny timeout")
    return mol


def test_embed_join_ignores_the_wall_clock_under_deterministic_mode():
    """_embed_with_timeout under DELFIN_DETERMINISTIC=1 waits for the embed
    even when the caller passed a tiny timeout; legacy mode honours it."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from delfin.manta.embed_timeout import _embed_with_timeout

    mol = _probe_mol()
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True

    saved = _clear_delfin_env()
    try:
        # ORDER IS FIXED: legacy first. The deterministic call warms RDKit's
        # process-wide ETKDG tables — measured on this RDKit (2025.09.6),
        # a warm quaterphenyl embed takes <50 ms, so a legacy probe run
        # after a deterministic one would finish "too fast" and measure
        # nothing. Cold legacy must time out; warm deterministic must wait.
        os.environ.pop("DELFIN_DETERMINISTIC", None)
        legacy = _embed_with_timeout(Chem.Mol(mol), params, timeout=0.05)
        assert legacy == -1, (
            "legacy embed finished within 50 ms — the probe molecule is "
            "too easy on this machine and the wall-clock join was not "
            "exercised; the legacy assertion measures nothing here"
        )

        os.environ["DELFIN_DETERMINISTIC"] = "1"
        det = _embed_with_timeout(Chem.Mol(mol), params, timeout=0.05)
        assert det == 0, (
            "embed FAILED under DELFIN_DETERMINISTIC=1 with a 50 ms caller "
            "timeout — the master switch no longer makes the join "
            "unbounded (embed_timeout.py:84-91); a loaded machine would "
            "drop heavy embeds again and the environment-dependent build "
            "class returns"
        )
    finally:
        _restore(saved)


def test_multiembed_join_yields_the_full_set_under_deterministic_mode():
    """_embed_multiple_confs_with_timeout: legacy + 50 ms timeout returns a
    partial conformer set; DELFIN_DETERMINISTIC=1 returns all 5."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from delfin.manta.embed_timeout import _embed_multiple_confs_with_timeout

    mol = _probe_mol()

    def _params():
        p = AllChem.ETKDGv3()
        p.randomSeed = 42
        p.useRandomCoords = True
        p.clearConfs = False
        return p

    saved = _clear_delfin_env()
    try:
        os.environ.pop("DELFIN_DETERMINISTIC", None)
        legacy = _embed_multiple_confs_with_timeout(
            Chem.Mol(mol), 5, _params(), timeout=0.05)
        assert len(legacy) < 5, (
            "legacy multiembed placed all 5 conformers within 50 ms — "
            "probe too easy on this machine, measured nothing"
        )

        os.environ["DELFIN_DETERMINISTIC"] = "1"
        det = _embed_multiple_confs_with_timeout(
            Chem.Mol(mol), 5, _params(), timeout=0.05)
        assert len(det) == 5, (
            "only %d/5 conformers under DELFIN_DETERMINISTIC=1 — the "
            "master switch no longer governs "
            "_embed_multiple_confs_with_timeout's join "
            "(embed_timeout.py:185-190)" % len(det)
        )
    finally:
        _restore(saved)
