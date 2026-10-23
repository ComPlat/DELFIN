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

The embed itself is slowed to a fixed 0.3 s (``_slow_embed``): a real
quaterphenyl embed finished inside 50 ms on the GitHub runner, so a probe
that relies on the machine being slow measures the machine, not the join.
The slowed call still runs RDKit's real embed afterwards.
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


_EMBED_DELAY_S = 0.3


@pytest.fixture(autouse=True)
def _slow_embed(monkeypatch):
    """Every embed takes at least 0.3 s, six times the 50 ms join timeout,
    whatever the machine: the legacy join must give up, the deterministic
    one must wait."""
    import time

    from delfin.manta import embed_timeout

    real = embed_timeout.AllChem

    class _SlowAllChem:
        def __getattr__(self, name):
            return getattr(real, name)

        @staticmethod
        def EmbedMolecule(*args, **kwargs):
            time.sleep(_EMBED_DELAY_S)
            return real.EmbedMolecule(*args, **kwargs)

        @staticmethod
        def EmbedMultipleConfs(*args, **kwargs):
            time.sleep(_EMBED_DELAY_S)
            return real.EmbedMultipleConfs(*args, **kwargs)

    monkeypatch.setattr(embed_timeout, "AllChem", _SlowAllChem())


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
        # The embed is held at 0.3 s by _slow_embed, so legacy must time
        # out and deterministic must wait, in either order and on any
        # machine.
        os.environ.pop("DELFIN_DETERMINISTIC", None)
        legacy = _embed_with_timeout(Chem.Mol(mol), params, timeout=0.05)
        assert legacy == -1, (
            "legacy embed returned a result although the embed took "
            "0.3 s against a 50 ms join — the wall-clock join no longer "
            "gives up in legacy mode"
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
            "legacy multiembed returned all 5 conformers although the "
            "embed took 0.3 s against a 50 ms join — the join no longer "
            "gives up in legacy mode"
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
