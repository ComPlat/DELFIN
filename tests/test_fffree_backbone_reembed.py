"""Tests for the core-preserving backbone re-embed conformer source
(delfin.manta.backbone_reembed, env DELFIN_FFFREE_BACKBONE_REEMBED).

Validates the four hard contracts of the feature:
  1. byte-identical when the flag is OFF (default) -> the recall path is untouched,
  2. the #82/#100 core invariant: metal + ALL donors stay native (<=0.05 A) in every
     re-embedded frame (the whole point -- whole-complex DG destroyed the core),
  3. the flag ON adds backbone-DIVERSE frames (heavy backbone genuinely re-folds),
  4. determinism: same input twice ON -> byte-identical output.

All graph-based / universal (no refcode- or SMILES-specific assertions).
"""
import os
import numpy as np
import pytest

from delfin.manta._bond_decollapse import _is_metal as bd_is_metal

# CN_EXTEND is required for the small CN2 test ligands to reach the FF-free path; the
# re-embed itself is independent of it.  RIGID_HAPTO/SIGMA stay default (off here).
_BASE_ENV = {"DELFIN_FFFREE_CN_EXTEND": "1"}

# A monodentate Werner complex (CN6, flexible-armed donors) and a chelate (CN6) so
# both grafting paths (single-donor pin + multidentate bite-constrained embed) run.
_MONO = "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]"             # CoCl3(NH3)3, monodentate
_FLEX = "CCO[Co](OCC)(OCC)([NH3])([NH3])[NH3]"                # flexible ethoxy arms
_CHELATE = "[NH2]CC[NH2][Co]1([NH2]CC[NH2]1)([Cl])[Cl]"      # bis(en) chelate (CN6)


def _isomers(smiles):
    from delfin.manta import converter_backend as CB
    return CB._fffree_isomers(smiles)


def _parse(xyz):
    syms = [ln.split()[0] for ln in xyz.splitlines()]
    P = np.array([[float(x) for x in ln.split()[1:4]] for ln in xyz.splitlines()])
    return syms, P


def _base_of(label):
    """The parent (native) base label of a re-embed frame: strip the trailing
    '-reembedN' suffix.  Each isomer/config has its OWN base frame + donor vertices,
    so a re-embed frame must be compared to ITS parent, not to an arbitrary frame."""
    i = label.rfind("-reembed")
    return label[:i] if i >= 0 else label


def _donors(syms, P):
    import delfin.manta._bond_decollapse as bd
    from delfin.manta import metal_sphere_builder as MSB
    m = next(i for i, s in enumerate(syms) if bd._is_metal(s))
    don = [j for j in range(len(syms))
           if j != m and syms[j] != "H"
           and float(np.linalg.norm(P[j] - P[m])) < 1.02 * MSB.md_distance(syms[m], syms[j])]
    return m, don


@pytest.fixture(autouse=True)
def _clean_env():
    saved = dict(os.environ)
    os.environ.pop("DELFIN_FFFREE_BACKBONE_REEMBED", None)
    for k, v in _BASE_ENV.items():
        os.environ[k] = v
    yield
    os.environ.clear()
    os.environ.update(saved)


def _sig(rr):
    if rr is None:
        return None
    return tuple((lab, xyz) for xyz, lab in rr)


def test_byte_identical_when_flag_off():
    """Flag unset == flag '0': no extra frames, output identical (path untouched)."""
    os.environ.pop("DELFIN_FFFREE_BACKBONE_REEMBED", None)
    a = _sig(_isomers(_MONO))
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "0"
    b = _sig(_isomers(_MONO))
    assert a == b
    # and no frame carries a -reembed label
    assert a is not None
    assert not any("reembed" in lab for lab, _ in a)


def test_flag_on_adds_reembed_frames():
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "0"
    off = _isomers(_FLEX)
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "1"
    on = _isomers(_FLEX)
    assert off is not None and on is not None
    n_re = sum("reembed" in lab for _, lab in on)
    assert n_re >= 1, "flag ON must add at least one backbone re-embed frame"
    assert len(on) > len(off)


def test_core_invariant_metal_and_donors_native():
    """The #82/#100 guard: in EVERY re-embed frame the metal and all donor atoms stay
    within 0.05 A of their PARENT base frame's positions (per isomer/config)."""
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "1"
    rr = _isomers(_FLEX)
    assert rr is not None
    bases = {lab: xyz for xyz, lab in rr if "reembed" not in lab}
    assert bases
    n_re = 0
    for xyz, lab in rr:
        if "reembed" not in lab:
            continue
        n_re += 1
        parent = _base_of(lab)
        assert parent in bases, f"orphan re-embed frame {lab}"
        bsyms, bP = _parse(bases[parent])
        bm, donors = _donors(bsyms, bP)
        assert donors, "expected at least one detected donor"
        syms, P = _parse(xyz)
        assert len(syms) == len(bsyms)
        m = next(i for i, s in enumerate(syms) if bd_is_metal(s))
        assert float(np.linalg.norm(P[m] - bP[bm])) <= 0.05
        for dgi in donors:
            assert float(np.linalg.norm(P[dgi] - bP[dgi])) <= 0.05, (
                f"donor {dgi} drifted in {lab}")
    assert n_re >= 1


def test_backbone_actually_refolds():
    """A re-embed frame must move heavy BACKBONE atoms (non-core) vs its parent --
    otherwise it is not a new conformer source, just a duplicate."""
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "1"
    rr = _isomers(_FLEX)
    bases = {lab: _parse(xyz) for xyz, lab in rr if "reembed" not in lab}
    moved = 0.0
    for xyz, lab in rr:
        if "reembed" not in lab:
            continue
        bsyms, bP = bases[_base_of(lab)]
        bm, donors = _donors(bsyms, bP)
        donors = set(donors)
        _, P = _parse(xyz)
        heavy = [i for i, s in enumerate(bsyms)
                 if s != "H" and i != bm and i not in donors]
        if heavy:
            mv = float(np.sqrt(((P[heavy] - bP[heavy]) ** 2).sum(axis=1).mean()))
            moved = max(moved, mv)
    assert moved > 0.10, "re-embed must produce a genuinely distinct backbone fold"


def test_deterministic():
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "1"
    a = _sig(_isomers(_FLEX))
    b = _sig(_isomers(_FLEX))
    assert a == b


def test_never_nonfinite():
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "1"
    for smi in (_MONO, _FLEX, _CHELATE):
        rr = _isomers(smi)
        assert rr is not None
        for xyz, _ in rr:
            _, P = _parse(xyz)
            assert np.all(np.isfinite(P))


def test_chelate_byte_identical_off_and_core_preserved_on():
    """Chelate path: OFF byte-identical; ON keeps the chelate donors native (the
    bite-constrained embed must not drift the coordinating N atoms)."""
    import delfin.manta._bond_decollapse as bd
    from delfin.manta import metal_sphere_builder as MSB
    os.environ.pop("DELFIN_FFFREE_BACKBONE_REEMBED", None)
    a = _sig(_isomers(_CHELATE))
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "0"
    b = _sig(_isomers(_CHELATE))
    assert a == b
    os.environ["DELFIN_FFFREE_BACKBONE_REEMBED"] = "1"
    rr = _isomers(_CHELATE)
    assert rr is not None and any("reembed" in lab for _, lab in rr)
    bases = {lab: _parse(xyz) for xyz, lab in rr if "reembed" not in lab}
    for xyz, lab in rr:
        if "reembed" not in lab:
            continue
        bsyms, bP = bases[_base_of(lab)]
        bm, donors = _donors(bsyms, bP)
        syms, P = _parse(xyz)
        m = next(i for i, s in enumerate(syms) if bd_is_metal(s))
        assert float(np.linalg.norm(P[m] - bP[bm])) <= 0.05
        for dgi in donors:
            assert float(np.linalg.norm(P[dgi] - bP[dgi])) <= 0.05
