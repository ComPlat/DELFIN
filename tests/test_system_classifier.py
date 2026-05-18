"""Tests for delfin._system_classifier (Welle-5m-Y).

Universal 8-axis classifier — pure-graph features, no SMILES regex.
Coverage: one test per axis + contrast SMILES descriptors.
"""

from __future__ import annotations

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem

from delfin._system_classifier import classify_complex_system


def _mol(smiles: str):
    """Parse SMILES tolerant of metal complex notation."""
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    assert m is not None, f"Failed to parse: {smiles}"
    return m


# ---------------------------------------------------------------------------
# Contrast SMILES from the brief (with `Ph` shorthand expanded — RDKit
# does not parse the `Ph` atom shortcut, this is a transcription quirk
# of the brief, not a SMILES feature).
# ---------------------------------------------------------------------------
SMILES_YIRQIC = (
    "[O+]#[C][Re-5]1([Br])([C]#[O+])"
    "([P+](c2ccccc2)(c2ccccc2)c2ccccc2)"
    "[C+]2N(c3ccccc3)C=CN2C2=CC=CC=[N+]21"
)
SMILES_YIVROM_LIKE = (
    # Tris-bidentate Fe(en)3 - style: three N,N chelates, all sigma.
    "[Fe-3]123([N]CC[N]1)([N]CC[N]2)[N]CC[N]3"
)
SMILES_SIGMA_HOMO = "[Cl][Ni]([Cl])([N])([N])"
SMILES_FERROCENE = "[CH-]1C=CC=C1.[CH-]1C=CC=C1.[Fe+2]"
SMILES_HEXACYANIDO = "[N]#C[Fe-4](C#[N])(C#[N])(C#[N])(C#[N])C#[N]"


# ---------------------------------------------------------------------------
# Axis 1 — class
# ---------------------------------------------------------------------------
def test_class_no_metal():
    d = classify_complex_system(_mol("CCO"))
    assert d["class"] == "no_metal"
    assert d["cell_key"] == "no_metal"


def test_class_sigma_vs_hapto():
    d_sigma = classify_complex_system(_mol(SMILES_SIGMA_HOMO))
    d_hapto = classify_complex_system(_mol(SMILES_FERROCENE))
    assert d_sigma["class"] == "sigma"
    # Ferrocene as two separate fragments + bare Fe+2 does not have
    # graph-level M-C edges, so it lands as sigma with CN 0.  Use a
    # bonded ferrocene representation to test hapto axis specifically.
    # Bonded Cp-like ring: 4 M-C edges → hapto group (η4)
    d_bonded_cp = classify_complex_system(_mol("C12=C3C4=C1[Fe]234"))
    assert d_bonded_cp["class"] == "hapto"
    assert d_bonded_cp["hapticity_max"] >= 4


# ---------------------------------------------------------------------------
# Axis 2 — CN
# ---------------------------------------------------------------------------
def test_cn_hexacoord_homo():
    d = classify_complex_system(_mol(SMILES_HEXACYANIDO))
    assert d["CN"] == 6
    assert d["donor_heterogeneity"] == "homo"


def test_cn_sigma_cn4():
    d = classify_complex_system(_mol(SMILES_SIGMA_HOMO))
    assert d["CN"] == 4


# ---------------------------------------------------------------------------
# Axis 3 — donor_heterogeneity (homo / bi / tri / maximal-asym)
# ---------------------------------------------------------------------------
def test_heterogeneity_axis():
    # homo (CN-Fe only): one element class
    assert classify_complex_system(_mol(SMILES_HEXACYANIDO))[
        "donor_heterogeneity"
    ] == "homo"
    # bi: Cl + N donors at Ni
    assert classify_complex_system(_mol(SMILES_SIGMA_HOMO))[
        "donor_heterogeneity"
    ] == "bi"
    # YIRQIC: C(CO) + Br + P + C(carbene) + N → C/Br/P/N = 4 classes → maximal-asym
    yir = classify_complex_system(_mol(SMILES_YIRQIC))
    assert yir["donor_heterogeneity"] == "maximal-asym"


# ---------------------------------------------------------------------------
# Axis 4 — chelate_pattern
# ---------------------------------------------------------------------------
def test_chelate_all_mono():
    d = classify_complex_system(_mol(SMILES_SIGMA_HOMO))
    assert d["chelate_pattern"] == "all-mono"


def test_chelate_tris_bid():
    d = classify_complex_system(_mol(SMILES_YIVROM_LIKE))
    # Three N,O bidentate ligands → tris-bid
    assert d["chelate_pattern"] in ("tris-bid", "pincer")
    # Three distinct ligand fragments × 2 donors each → no fragment
    # contributes >=3 donors, so pincer should not fire.
    assert d["chelate_pattern"] == "tris-bid"


# ---------------------------------------------------------------------------
# Axis 5 — metal_block
# ---------------------------------------------------------------------------
def test_metal_block_axis():
    # Ni (Z=28) → 3d
    assert classify_complex_system(_mol(SMILES_SIGMA_HOMO))[
        "metal_block"
    ] == "3d"
    # Fe (Z=26) → 3d
    assert classify_complex_system(_mol(SMILES_HEXACYANIDO))[
        "metal_block"
    ] == "3d"
    # Re (Z=75) → 5d
    assert classify_complex_system(_mol(SMILES_YIRQIC))[
        "metal_block"
    ] == "5d"
    # Sn (Z=50) → main-group
    d_sn = classify_complex_system(_mol("[Cl][Sn]([Cl])([Cl])[Cl]"))
    assert d_sn["metal_block"] == "main-group"
    # Eu (Z=63) → f-block
    d_eu = classify_complex_system(_mol("[Cl][Eu]([Cl])[Cl]"))
    assert d_eu["metal_block"] == "f-block"


# ---------------------------------------------------------------------------
# Axis 6 — q_magnitude
# ---------------------------------------------------------------------------
def test_q_magnitude_axis():
    assert classify_complex_system(_mol(SMILES_YIRQIC))["q_magnitude"] == 5
    assert classify_complex_system(_mol(SMILES_HEXACYANIDO))[
        "q_magnitude"
    ] == 4
    assert classify_complex_system(_mol(SMILES_SIGMA_HOMO))[
        "q_magnitude"
    ] == 0


# ---------------------------------------------------------------------------
# Axis 7 — hapticity_max
# ---------------------------------------------------------------------------
def test_hapticity_axis():
    assert classify_complex_system(_mol(SMILES_SIGMA_HOMO))[
        "hapticity_max"
    ] == 0
    d_cp = classify_complex_system(_mol("C12=C3C4=C1[Fe]234"))
    assert d_cp["hapticity_max"] >= 4


# ---------------------------------------------------------------------------
# Axis 8 — bulky_flag (tBu graph signature: sp3-C with 3 sp3-C nbrs)
# ---------------------------------------------------------------------------
def test_bulky_flag_axis():
    # No bulky substituents
    assert (
        classify_complex_system(_mol(SMILES_SIGMA_HOMO))["bulky_flag"]
        is False
    )
    # P-tBu3 ligand on Pd: sp3 C with 3 sp3-C nbrs near P donor
    d_tbu = classify_complex_system(
        _mol("[Cl][Pd]([Cl])[P](C(C)(C)C)(C(C)(C)C)C(C)(C)C")
    )
    assert d_tbu["bulky_flag"] is True


# ---------------------------------------------------------------------------
# cell_key stability + format
# ---------------------------------------------------------------------------
def test_cell_key_format():
    d = classify_complex_system(_mol(SMILES_SIGMA_HOMO))
    parts = d["cell_key"].split("_")
    # class _ CNn _ het _ chelate _ block _ qN _ etaN _ bulkyflag
    assert parts[0] == "sigma"
    assert parts[1] == "CN4"
    assert parts[5] == "q0"
    assert parts[6] == "eta0"
    assert parts[7] == "nobulky"


def test_cell_key_yirqic_smoke():
    d = classify_complex_system(_mol(SMILES_YIRQIC))
    assert d["class"] == "sigma"
    assert d["metal_block"] == "5d"
    assert d["q_magnitude"] == 5
    assert d["hapticity_max"] == 0
    # cell key should be stable, deterministic
    k1 = d["cell_key"]
    k2 = classify_complex_system(_mol(SMILES_YIRQIC))["cell_key"]
    assert k1 == k2


# ---------------------------------------------------------------------------
# Null safety
# ---------------------------------------------------------------------------
def test_none_mol_returns_no_metal():
    d = classify_complex_system(None)
    assert d["class"] == "no_metal"
    assert d["cell_key"] == "no_metal"
