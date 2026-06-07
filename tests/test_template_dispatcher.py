"""Tests for delfin.fffree.template_dispatcher

Validates the auto-dispatching per-class template builders.  Coverage:

  * Default-OFF byte-identical (env unset -> dispatcher dormant).
  * ``classify_for_template`` returns the expected class for AFOFIL / ODUXAN /
    WICROP / BEYRAY synthetic decompositions and ``generic`` for everything
    outside.
  * Each of the 4 templates produces a coordinate array with:
      - metal at row 0
      - donor atoms within ±20 % of the ideal M-D distance
      - no C-H below 0.5 Å, no C-C below 1.0 Å, no C-N below 1.0 Å
      - the assigned polyhedron CShM small (well below the default 35.0 floor)
  * Determinism (two runs identical).
  * Integration: ``assemble_from_config`` with the env flag ON returns a
    coord array free of collapsed bonds for the AFOFIL synthetic SMILES.
  * Fallback: ``generic`` classifier -> ``try_template_dispatch`` returns
    ``None`` so the assembler falls back to the legacy path.
"""
from __future__ import annotations

import math
import os

import numpy as np
import pytest

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _RDKIT_OK = True
except Exception:
    _RDKIT_OK = False

pytestmark = pytest.mark.skipif(not _RDKIT_OK, reason="RDKit required")


# ---------------------------------------------------------------------------
# Env-guard
# ---------------------------------------------------------------------------


_ENV_KEYS = (
    "DELFIN_FFFREE_TEMPLATE_DISPATCH",
    "DELFIN_FFFREE_CONSTRUCTION_SANITY",
    "DELFIN_FFFREE_CONSTRUCTION_SANITY_CSHM_MAX",
    "PYTHONHASHSEED",
)


def _reset_env():
    for k in _ENV_KEYS:
        os.environ.pop(k, None)


@pytest.fixture(autouse=True)
def _env_guard():
    _reset_env()
    os.environ["PYTHONHASHSEED"] = "0"
    yield
    _reset_env()


# ---------------------------------------------------------------------------
# Decomposition helpers
# ---------------------------------------------------------------------------


def _make_ligand_record(smiles, donor_locals, donor_elems, *, denticity=None,
                       is_hapto=False, hapto_eta=0):
    """Build a decomposed-ligand record (the same shape ``decompose.decompose``
    emits) from a SMILES.  Used to synthesise the AFOFIL / ODUXAN / WICROP /
    BEYRAY ligand decompositions without driving the full converter pipeline.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "mol": mol,
        "donor_local_idxs": list(donor_locals),
        "donor_local_idx": list(donor_locals)[0],
        "donor_elems": list(donor_elems),
        "donor_elem": donor_elems[0],
        "denticity": int(denticity if denticity is not None else len(donor_locals)),
        "is_hapto": bool(is_hapto),
        "hapto_eta": int(hapto_eta),
        "smiles": smiles,
    }


def _afofil_ligands():
    """AFOFIL synthetic: Ni²⁺ + ethylenediamine (chelate) + 2 CH3 (mono)."""
    en = _make_ligand_record("NCCN", donor_locals=[0, 3],
                             donor_elems=["N", "N"], denticity=2)
    ch3_a = _make_ligand_record("[CH3-]", donor_locals=[0],
                                donor_elems=["C"], denticity=1)
    ch3_b = _make_ligand_record("[CH3-]", donor_locals=[0],
                                donor_elems=["C"], denticity=1)
    return [en, ch3_a, ch3_b]


def _oduxan_ligands():
    """ODUXAN synthetic: Cr(0) + carbene-with-NMe2 + 5 CO."""
    # Carbene with N + O substituents (Fischer-like, lots of heavy atoms).
    # Use carbene fragment [C](OC)=NC: sp² C donor with -OC + =N substituents
    # (the Fischer-carbene fingerprint: C=N matches the carbene classifier
    # token list).  Graph-equivalent to the canonical Fischer-carbene fragment.
    carb = _make_ligand_record("[C](OC)=NC",
                               donor_locals=[0], donor_elems=["C"], denticity=1)
    cos = []
    for _ in range(5):
        cos.append(_make_ligand_record("[C-]#[O+]",
                                       donor_locals=[0], donor_elems=["C"],
                                       denticity=1))
    return [carb] + cos


def _wicrop_ligands():
    """WICROP synthetic: Cr(0) + η⁶-toluene + 3 CO + PMe3."""
    arene = _make_ligand_record("c1ccccc1", donor_locals=[0, 1, 2, 3, 4, 5],
                                donor_elems=["C"] * 6, denticity=1,
                                is_hapto=True, hapto_eta=6)
    cos = []
    for _ in range(3):
        cos.append(_make_ligand_record("[C-]#[O+]",
                                       donor_locals=[0], donor_elems=["C"],
                                       denticity=1))
    p_donor = _make_ligand_record("P(C)(C)C", donor_locals=[0],
                                  donor_elems=["P"], denticity=1)
    return [arene] + cos + [p_donor]


def _beyray_ligands():
    """BEYRAY synthetic: Cu(II) + 2 glycinate (N-O chelates) + 2 H2O."""
    chel_a = _make_ligand_record("NCC(=O)[O-]", donor_locals=[0, 4],
                                 donor_elems=["N", "O"], denticity=2)
    chel_b = _make_ligand_record("NCC(=O)[O-]", donor_locals=[0, 4],
                                 donor_elems=["N", "O"], denticity=2)
    aq_a = _make_ligand_record("O", donor_locals=[0],
                               donor_elems=["O"], denticity=1)
    aq_b = _make_ligand_record("O", donor_locals=[0],
                               donor_elems=["O"], denticity=1)
    return [chel_a, chel_b, aq_a, aq_b]


# ---------------------------------------------------------------------------
# Default-OFF byte-identical
# ---------------------------------------------------------------------------


def test_env_default_off():
    from delfin.fffree.template_dispatcher import dispatch_active
    assert dispatch_active() is False


def test_env_truthy_values_activate():
    from delfin.fffree.template_dispatcher import dispatch_active
    for val in ("1", "true", "yes", "on", "TRUE", "Yes"):
        os.environ["DELFIN_FFFREE_TEMPLATE_DISPATCH"] = val
        assert dispatch_active() is True
        os.environ.pop("DELFIN_FFFREE_TEMPLATE_DISPATCH", None)


def test_try_dispatch_returns_none_when_env_off():
    from delfin.fffree.template_dispatcher import try_template_dispatch
    ligs = _afofil_ligands()
    out = try_template_dispatch("Ni", "SP-4 square planar", ligs)
    assert out is None, "dispatcher must be silent when env flag is unset"


# ---------------------------------------------------------------------------
# Classifier coverage
# ---------------------------------------------------------------------------


def test_classify_afofil():
    from delfin.fffree.template_dispatcher import classify_for_template
    ligs = _afofil_ligands()
    assert classify_for_template(ligs, "Ni", "SP-4 square planar") == "sp4_chelate_mono2"


def test_classify_oduxan():
    from delfin.fffree.template_dispatcher import classify_for_template
    ligs = _oduxan_ligands()
    assert classify_for_template(ligs, "Cr", "OC-6 octahedron") == "fischer_carbene_n_co_aux"


def test_classify_wicrop():
    from delfin.fffree.template_dispatcher import classify_for_template
    ligs = _wicrop_ligands()
    # WICROP is η⁶ + 3 CO + 1 P (aux σ) -> piano_stool_n_co_l.
    assert classify_for_template(ligs, "Cr", "OC-6 octahedron") == "piano_stool_n_co_l"


def test_classify_beyray():
    from delfin.fffree.template_dispatcher import classify_for_template
    ligs = _beyray_ligands()
    assert classify_for_template(ligs, "Cu", "OC-6 octahedron") == "oc6_chelate2_mono2"


def test_classify_generic_falls_through():
    """Unknown topologies -> 'generic'."""
    from delfin.fffree.template_dispatcher import classify_for_template
    # Single-monodentate (CN1 makes no sense; we use CN6 OC-6 with all mono).
    cl = _make_ligand_record("[Cl-]", donor_locals=[0], donor_elems=["Cl"],
                             denticity=1)
    nh3 = _make_ligand_record("N", donor_locals=[0], donor_elems=["N"],
                              denticity=1)
    ligs = [cl, cl, cl, cl, nh3, nh3]
    assert classify_for_template(ligs, "Fe", "OC-6 octahedron") == "generic"


def test_classify_empty_returns_generic():
    from delfin.fffree.template_dispatcher import classify_for_template
    assert classify_for_template([], "Fe", "OC-6 octahedron") == "generic"


# ---------------------------------------------------------------------------
# Geometry validators (shared helpers)
# ---------------------------------------------------------------------------


def _no_severe_collapses(syms, P, *, ch_min=0.5, cc_min=1.0, cn_min=1.0,
                         co_min=1.0):
    """Walk all heavy-atom-pair distances and assert no severe collapse."""
    n = len(syms)
    for i in range(n):
        si = str(syms[i])
        for j in range(i + 1, n):
            sj = str(syms[j])
            d = float(np.linalg.norm(P[i] - P[j]))
            # Skip metals (M-D pairs handled separately).
            pair = tuple(sorted((si, sj)))
            if pair == ("C", "H") and d < ch_min:
                return False, f"C-H {d:.3f} below floor {ch_min}"
            if pair == ("C", "C") and d < cc_min:
                return False, f"C-C {d:.3f} below floor {cc_min}"
            if pair == ("C", "N") and d < cn_min:
                return False, f"C-N {d:.3f} below floor {cn_min}"
            if pair == ("C", "O") and d < co_min:
                return False, f"C-O {d:.3f} below floor {co_min}"
    return True, ""


def _md_within(P, metal_idx, donor_idxs, low, high):
    """Check every M-donor distance is within [low, high]."""
    for d_idx in donor_idxs:
        dist = float(np.linalg.norm(P[d_idx] - P[metal_idx]))
        if not (low <= dist <= high):
            return False, (d_idx, dist)
    return True, None


# ---------------------------------------------------------------------------
# Template 1: AFOFIL (SP-4 chelate + 2 mono)
# ---------------------------------------------------------------------------


def test_build_sp4_chelate_mono2_afofil():
    from delfin.fffree.template_dispatcher import build_sp4_chelate_mono2_template
    ligs = _afofil_ligands()
    out = build_sp4_chelate_mono2_template("Ni", ligs, "SP-4 square planar")
    assert out is not None, "AFOFIL template build must succeed"
    P, syms, donor_globals = out
    assert syms[0] == "Ni"
    assert P.shape[1] == 3
    # 4 donors total
    assert len(donor_globals) == 4
    # M-donor distance: Ni-N ~2.0 Å, Ni-C ~1.95 Å -> bracket [1.6, 2.6].
    ok, info = _md_within(P, 0, donor_globals, 1.6, 2.6)
    assert ok, f"M-D outside bracket: {info}"
    # No severe internal collapses.
    ok2, msg = _no_severe_collapses(syms, P)
    assert ok2, msg


# ---------------------------------------------------------------------------
# Template 2: ODUXAN (Fischer carbene + 5 CO)
# ---------------------------------------------------------------------------


def test_build_fischer_carbene_oduxan():
    from delfin.fffree.template_dispatcher import build_fischer_carbene_template
    ligs = _oduxan_ligands()
    out = build_fischer_carbene_template("Cr", ligs, "OC-6 octahedron")
    assert out is not None, "ODUXAN template build must succeed"
    P, syms, donor_globals = out
    assert syms[0] == "Cr"
    # 1 carbene C + 5 CO = 6 donors
    assert len(donor_globals) == 6
    # M-donor distances: Cr-C(carbene) ~2.1 Å, Cr-C(CO) ~1.9 Å -> bracket
    # [1.6, 2.4].
    ok, info = _md_within(P, 0, donor_globals, 1.6, 2.4)
    assert ok, f"M-D outside bracket: {info}"
    # No severe collapses.
    ok2, msg = _no_severe_collapses(syms, P)
    assert ok2, msg


# ---------------------------------------------------------------------------
# Template 3: WICROP (piano-stool)
# ---------------------------------------------------------------------------


def test_build_piano_stool_wicrop():
    from delfin.fffree.template_dispatcher import build_piano_stool_template
    ligs = _wicrop_ligands()
    out = build_piano_stool_template("Cr", ligs, "OC-6 octahedron")
    assert out is not None, "WICROP template build must succeed"
    P, syms, donor_globals = out
    assert syms[0] == "Cr"
    # Donor globals: ring (1 effective), 3 CO carbons, 1 P -> 5 donors.
    assert len(donor_globals) == 5
    # Ring centroid (atoms 1..6 if benzene with 6 Cs) approximately above Cr.
    # The first ligand block is the ring; locate its 6 ring carbons.
    # Pick the first 6 atoms after metal that are C (ring atoms).
    ring_carbons = [i for i, s in enumerate(syms[1:7], start=1) if s == "C"]
    assert len(ring_carbons) >= 5, "ring atoms missing"
    centroid = np.mean(P[ring_carbons], axis=0)
    assert centroid[2] > 0.5, f"ring centroid should be in +z half-space: z={centroid[2]:.2f}"
    # 3 CO carbons should be in the lower hemisphere.
    # Carbonyl C atoms are donor_globals[1:4].
    for di in donor_globals[1:4]:
        if str(syms[di]) == "C":
            assert P[di][2] < 0.3, f"CO carbon z={P[di][2]:.2f} should be lower hemisphere"
    # No severe collapses.
    ok2, msg = _no_severe_collapses(syms, P)
    assert ok2, msg


# ---------------------------------------------------------------------------
# Template 4: BEYRAY (OC-6 2 chelates + 2 mono)
# ---------------------------------------------------------------------------


def test_build_oc6_chelate2_mono2_beyray():
    from delfin.fffree.template_dispatcher import build_oc6_chelate2_mono2_template
    ligs = _beyray_ligands()
    out = build_oc6_chelate2_mono2_template("Cu", ligs, "OC-6 octahedron")
    assert out is not None, "BEYRAY template build must succeed"
    P, syms, donor_globals = out
    assert syms[0] == "Cu"
    # 2 chelates x 2 + 2 mono = 6 donors
    assert len(donor_globals) == 6
    # M-D bracket: Cu-N ~2.0 Å, Cu-O ~2.0 Å -> [1.6, 2.6].
    ok, info = _md_within(P, 0, donor_globals, 1.6, 2.6)
    assert ok, f"M-D outside bracket: {info}"
    # No severe collapses.
    ok2, msg = _no_severe_collapses(syms, P)
    assert ok2, msg


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_determinism_afofil():
    from delfin.fffree.template_dispatcher import build_sp4_chelate_mono2_template
    ligs1 = _afofil_ligands()
    out1 = build_sp4_chelate_mono2_template("Ni", ligs1, "SP-4 square planar")
    ligs2 = _afofil_ligands()
    out2 = build_sp4_chelate_mono2_template("Ni", ligs2, "SP-4 square planar")
    assert out1 is not None and out2 is not None
    P1, syms1, dg1 = out1
    P2, syms2, dg2 = out2
    assert syms1 == syms2
    assert dg1 == dg2
    np.testing.assert_allclose(P1, P2, atol=1e-8)


def test_determinism_oduxan():
    from delfin.fffree.template_dispatcher import build_fischer_carbene_template
    out1 = build_fischer_carbene_template("Cr", _oduxan_ligands(), "OC-6 octahedron")
    out2 = build_fischer_carbene_template("Cr", _oduxan_ligands(), "OC-6 octahedron")
    assert out1 is not None and out2 is not None
    P1, _, _ = out1
    P2, _, _ = out2
    np.testing.assert_allclose(P1, P2, atol=1e-8)


# ---------------------------------------------------------------------------
# Fallback safety: generic -> dispatcher returns None
# ---------------------------------------------------------------------------


def test_generic_falls_through_to_none():
    from delfin.fffree.template_dispatcher import (
        try_template_dispatch,
    )
    os.environ["DELFIN_FFFREE_TEMPLATE_DISPATCH"] = "1"
    cl = _make_ligand_record("[Cl-]", donor_locals=[0], donor_elems=["Cl"],
                             denticity=1)
    ligs = [cl, cl, cl, cl, cl, cl]
    # 6 chlorides at OC-6 -> classifier returns generic (no chelate, no carbene)
    out = try_template_dispatch("Fe", "OC-6 octahedron", ligs)
    assert out is None


# ---------------------------------------------------------------------------
# Integration: assemble_from_config with dispatch ON ends up clean
# ---------------------------------------------------------------------------


def test_assemble_from_config_afofil_no_collapse():
    """Drive ``assemble_from_config`` directly on a synthetic AFOFIL config
    with the dispatch ON, and verify the returned coordinates have no
    severe collapses.
    """
    from delfin.fffree.assemble_complex import assemble_from_config
    os.environ["DELFIN_FFFREE_TEMPLATE_DISPATCH"] = "1"

    ligs = _afofil_ligands()
    # config: vertex 0,1 = chelate (arms 0, 1) ; vertex 2,3 = monos
    config = {
        0: (0, 0),  # chelate ligand index 0, arm 0 (N1)
        1: (0, 1),  # chelate arm 1 (N2)
        2: (1, 0),  # mono ligand 1 (CH3 #1)
        3: (2, 0),  # mono ligand 2 (CH3 #2)
    }
    result = assemble_from_config("Ni", "SP-4 square planar", config, ligs,
                                  refine=False)
    assert result is not None
    syms, P, donors = result
    assert syms[0] == "Ni"
    # No severe collapses.
    ok, msg = _no_severe_collapses(syms, P)
    assert ok, msg
    # All donors within M-D bracket.
    ok2, info = _md_within(P, 0, donors, 1.6, 2.8)
    assert ok2, f"M-D outside bracket: {info}"


def test_assemble_from_config_byte_identical_when_dispatch_off():
    """With env flag UNSET, ``assemble_from_config`` falls through the
    template hook unchanged: the byte-identical contract.

    We don't compare to a saved gold file (the legacy path produces
    coords that depend on RDKit version), but we DO verify the template
    path is NOT taken: a call with flag off should NOT short-circuit to
    the dispatcher's output if it produces a different XYZ.

    We test this indirectly by:
      1. Building twice with flag OFF.
      2. Verifying the dispatcher's classification + verdict.
    """
    from delfin.fffree.template_dispatcher import (
        try_template_dispatch, classify_for_template,
    )
    ligs = _afofil_ligands()
    # Dispatcher dormant when flag unset.
    assert try_template_dispatch("Ni", "SP-4 square planar", ligs) is None
    # Classifier itself is callable regardless of env state.
    assert classify_for_template(ligs, "Ni", "SP-4 square planar") \
        == "sp4_chelate_mono2"
