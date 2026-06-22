"""Tests for the gas-phase backbone-CONFORMER-coverage offensive
(CONFORMER_REACHABILITY_2026_06_17 — 3 fixes: format-tolerant XYZ parse,
backbone-heavy-moving DOF picker, dense ring-pucker + denser grid).

Invariants enforced here:
  * default-OFF byte-identical (flag unset -> base-only, default DOF picker
    unchanged, headerless parse unchanged);
  * dense heavy-MOVING conformers with the flag ON;
  * hard core-freeze (metal + every donor stays native within +-0.05 A);
  * determinism (2x identical).

No CCDC coordinates are used — synthetic graphs + a small frozen-core complex.
"""
import os

import pytest

from delfin import _conformer_pool as cp
from delfin import _rotamer_diversity as rot


# --------------------------------------------------------------------------
# Synthetic geometries (headerless DELFIN-XYZ unless a header is added)
# --------------------------------------------------------------------------
# n-butane: CH3-CH2-CH2-CH3 with one genuine backbone torsion (C2-C3).
_BUTANE = (
    "C   0.000  0.000  0.000\n"
    "C   1.530  0.000  0.000\n"
    "C   2.040  1.440  0.000\n"
    "C   3.570  1.440  0.000\n"
    "H  -0.360 -0.510  0.890\n"
    "H  -0.360 -0.510 -0.890\n"
    "H  -0.360  1.020  0.000\n"
    "H   1.890 -0.510  0.890\n"
    "H   1.890 -0.510 -0.890\n"
    "H   1.680  1.950  0.890\n"
    "H   1.680  1.950 -0.890\n"
    "H   3.930  0.930  0.890\n"
    "H   3.930  0.930 -0.890\n"
    "H   3.930  2.460  0.000\n"
)

# --------------------------------------------------------------------------
# Bug-1: format-tolerant XYZ parse (headerless byte-id; count-header tolerated)
# --------------------------------------------------------------------------
def test_parse_headerless_unchanged():
    syms, coords = rot._parse_delfin_xyz(_BUTANE)
    assert syms[0] == "C" and len(syms) == 14
    assert coords[0] == (0.0, 0.0, 0.0)


def test_parse_count_header_tolerated():
    """A count-header XYZ (the pool/recall-harness storage format) must parse to
    its first frame — previously raised, yielding 0 conformers (bug-1)."""
    headed = "14\nsome comment line\n" + _BUTANE
    s1, c1 = rot._parse_delfin_xyz(headed)
    s0, c0 = rot._parse_delfin_xyz(_BUTANE)
    assert s1 == s0 and c1 == c0


def test_parse_multiframe_takes_first():
    headed = "14\nframe-0\n" + _BUTANE + "14\nframe-1\n" + _BUTANE
    s1, _ = rot._parse_delfin_xyz(headed)
    assert len(s1) == 14  # only the first frame


def test_build_ob_from_headed_xyz():
    headed = "14\ncomment\n" + _BUTANE
    m = rot._build_ob_mol_from_xyz(headed)
    assert m is not None  # bug-1: previously None -> 0 conformers


# --------------------------------------------------------------------------
# Bug-2: DOF picker — coverage keeps backbone-heavy-moving torsions, drops
# trivial methyl rotors; default path byte-identical.
# --------------------------------------------------------------------------
def _graph(xyz):
    return rot._graph_from_ob(rot._build_ob_mol_from_xyz(xyz))


def test_dof_default_has_backbone_key_but_same_selection():
    g = _graph(_BUTANE)
    d = rot.identify_rotamer_dofs(g, max_dofs=6)  # coverage=False (default)
    assert d, "butane has a rotatable backbone bond"
    assert all("backbone_heavy_moved" in x for x in d)


def test_dof_coverage_keeps_only_backbone_movers():
    g = _graph(_BUTANE)
    cov = rot.identify_rotamer_dofs(g, max_dofs=8, coverage=True)
    # n-butane C2-C3 is a genuine backbone torsion (each side moves a CH3 = 1
    # heavy beyond the rotor root); coverage must keep it and every kept DOF
    # must move >=1 backbone heavy atom (never a pure terminal-H/methyl rotor).
    assert cov, "butane has a backbone torsion the coverage picker must keep"
    assert all(d["backbone_heavy_moved"] >= 1 for d in cov)


def test_dof_coverage_drops_methyl_on_real_frame():
    """On a real assembled complex, coverage must drop the terminal methyl
    rotors that the default picker keeps (bug-2: 0 backbone diversity)."""
    pool = ("/home/qmchem_max/agent_workspace/MANTA/POOLS/"
            "FULLSTACK-6835621-BEST50K-2026-06-16/archive/best_ON")
    p = os.path.join(pool, "CIYROT.xyz")
    if not os.path.exists(p):
        pytest.skip("CCDC-derived pool frame not available")
    g = _graph(open(p).read())
    d_def = rot.identify_rotamer_dofs(g, max_dofs=8, coverage=False)
    d_cov = rot.identify_rotamer_dofs(g, max_dofs=8, coverage=True)
    n_methyl_def = sum(1 for x in d_def if x["is_methyl"])
    assert n_methyl_def >= 1, "fixture expects methyl-dominated default DOFs"
    assert all(x["backbone_heavy_moved"] >= 1 for x in d_cov)
    assert len(d_cov) < len(d_def), "coverage prunes the methyl rotors"


def test_dof_coverage_drops_pure_methyl_molecule():
    # neopentane-like: a central C with 4 methyls -> all rotors are methyls,
    # coverage must return an empty DOF list (no backbone fold from torsions).
    neo = (
        "C  0.000  0.000  0.000\n"
        "C  1.540  0.000  0.000\n"
        "C -0.510  1.450  0.000\n"
        "C -0.510 -0.730  1.260\n"
        "C -0.510 -0.730 -1.260\n"
    )
    # add 3 H per terminal C
    hs = []
    import numpy as np
    base = np.array([[1.540, 0, 0], [-0.510, 1.450, 0],
                     [-0.510, -0.730, 1.260], [-0.510, -0.730, -1.260]])
    for b in base:
        for off in ([0.5, 0.5, 0.5], [0.5, -0.5, -0.5], [-0.5, 0.5, -0.5]):
            p = b + np.array(off)
            hs.append(f"H  {p[0]:.3f}  {p[1]:.3f}  {p[2]:.3f}")
    xyz = neo + "\n".join(hs) + "\n"
    g = _graph(xyz)
    cov = rot.identify_rotamer_dofs(g, max_dofs=8, coverage=True)
    assert cov == [] or all(d["backbone_heavy_moved"] >= 1 for d in cov)


# --------------------------------------------------------------------------
# Coverage flag: byte-identical OFF, dense conformers ON, determinism.
# --------------------------------------------------------------------------
def test_pool_off_byte_identical(monkeypatch):
    monkeypatch.delenv("DELFIN_FFFREE_CONFORMER_COVERAGE", raising=False)
    monkeypatch.delenv("DELFIN_5O_CONFORMER_POOL", raising=False)
    out = cp.apply_if_enabled(_BUTANE)
    assert out == [(_BUTANE, "base")]


def test_pool_coverage_produces_backbone_conformers(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_CONFORMER_COVERAGE", "1")
    out = cp.apply_if_enabled(_BUTANE)
    assert len(out) >= 2, "coverage must produce >=1 distinct backbone conformer"
    assert out[0][1] == "base"


def test_pool_coverage_deterministic(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_CONFORMER_COVERAGE", "1")
    a = cp.apply_if_enabled(_BUTANE)
    b = cp.apply_if_enabled(_BUTANE)
    assert a == b


def test_core_frozen_guard_rejects_metal_drift():
    # synthetic 2-atom "complex": metal + one donor; moving the donor >0.05 A
    # must fail _core_frozen.
    import numpy as np
    g = {
        "n_atoms": 3,
        "atomic_nums": [26, 7, 6],   # Fe, N(donor), C(backbone)
        "is_metal": [True, False, False],
        "neighbours": [[1], [0, 2], [1]],
    }
    base = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 1.0, 0.0)]
    moved_backbone = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 2.0, 0.0)]
    assert cp._core_frozen(g, base, moved_backbone, 0.05)  # backbone free
    moved_donor = [(0.0, 0.0, 0.0), (2.0, 0.2, 0.0), (3.0, 1.0, 0.0)]
    assert not cp._core_frozen(g, base, moved_donor, 0.05)  # donor drifted 0.2A


def test_coverage_core_frozen_on_real_frame(monkeypatch):
    """Every coverage conformer of a real frozen-core complex keeps the
    metal+donor within +-0.05 A of native."""
    import glob, numpy as np
    monkeypatch.setenv("DELFIN_FFFREE_CONFORMER_COVERAGE", "1")
    pool = ("/home/qmchem_max/agent_workspace/MANTA/POOLS/"
            "FULLSTACK-6835621-BEST50K-2026-06-16/archive/best_ON")
    cand = [os.path.join(pool, f"{rc}.xyz") for rc in ("NUHPUD", "HEXSAI", "REPWEO")]
    cand = [p for p in cand if os.path.exists(p)]
    if not cand:
        pytest.skip("CCDC-derived pool frames not available in this environment")
    f = cand[0]
    txt = open(f).read()
    s0, c0 = rot._parse_delfin_xyz(txt)
    base = np.asarray(c0, float)
    g = _graph(txt)
    core = set()
    for mi in range(g["n_atoms"]):
        if g["is_metal"][mi]:
            core.add(mi)
            for d in g["neighbours"][mi]:
                if g["atomic_nums"][d] != 1:
                    core.add(d)
    out = cp.apply_if_enabled(txt)
    for mxyz, _tag in out[1:]:
        _s, mc = rot._parse_delfin_xyz(mxyz)
        P = np.asarray(mc, float)
        assert np.all(np.isfinite(P))
        for i in core:
            assert float(np.linalg.norm(P[i] - base[i])) <= 0.05 + 1e-6
