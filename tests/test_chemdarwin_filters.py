"""Tests for the ChemDarwin mutation validity guards.

Covers the three always-on, seed-relative guards added to the mutation pipeline
in ``delfin.dashboard.tab_chemdarwin`` — fragmentation, radicals/open valences
and runaway growth — plus the complex-reassembly fragmentation guard.
"""

from __future__ import annotations

import pytest

pytest.importorskip("ipywidgets")
Chem = pytest.importorskip("rdkit.Chem")

from delfin.dashboard import tab_chemdarwin as cd


def _benign(smi):
    """Aromatic CH -> C-Cl substitution; used as the per-ligand mutation func."""
    return cd.apply_custom_reaction_iter(smi, "[cH:1]>>[c:1]Cl", iterations=1, keep_rings=True)


# --------------------------------------------------------------------------
# Organic path (apply_custom_reaction_iter)
# --------------------------------------------------------------------------

def test_benign_substitution_passes():
    res = cd.apply_custom_reaction_iter(
        "c1ccccc1", "[cH:1]>>[c:1]Cl", iterations=1, keep_rings=True
    )
    assert res
    assert all("." not in s for s, _ in res)


def test_radical_products_filtered():
    # Removing an H from a methyl leaves an open-valence carbon radical.
    with pytest.raises(ValueError, match="radicals"):
        cd.apply_custom_reaction_iter(
            "CCC", "[CH3:1]>>[CH2:1]", iterations=1, keep_rings=False
        )


def test_growth_limit_filters_oversized():
    # seed has 1 heavy atom -> cap = round(1 * 3.0) = 3; product has 11 -> rejected.
    with pytest.raises(ValueError, match="oversized"):
        cd.apply_custom_reaction_iter(
            "C", "[CH4:1]>>[C:1]CCCCCCCCCC", iterations=1, keep_rings=False
        )


def test_growth_limit_can_be_disabled():
    res = cd.apply_custom_reaction_iter(
        "C", "[CH4:1]>>[C:1]CCCCCCCCCC",
        iterations=1, keep_rings=False, max_growth_factor=0,
    )
    assert res
    assert any(m.GetNumHeavyAtoms() > 3 for _, m in res)


# --------------------------------------------------------------------------
# Complex path (mutate_complex_ligands)
# --------------------------------------------------------------------------

_COMPLEX_SEED = "[Fe+2](<-n1ccccc1)<-n1ccccc1"


def test_complex_path_produces_only_connected_complexes():
    res = cd.mutate_complex_ligands(_COMPLEX_SEED, _benign)
    assert res
    assert all("." not in s for _, s in res)


def test_complex_fragmentation_guard_drops_disconnected(monkeypatch):
    """If reassembly yields a disconnected complex, the guard must drop it."""
    orig = cd.reassemble_complex_from_mols

    def fragmenting(*args, **kwargs):
        combo = orig(*args, **kwargs)
        if combo is None:
            return None
        rw = Chem.RWMol(combo)
        rw.InsertMol(Chem.MolFromSmiles("c1ccccc1"))  # add a disconnected fragment
        return rw.GetMol()

    monkeypatch.setattr(cd, "reassemble_complex_from_mols", fragmenting)
    res = cd.mutate_complex_ligands(_COMPLEX_SEED, _benign)
    assert res == []
