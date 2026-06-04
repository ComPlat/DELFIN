"""POLYA-COVERAGE-FIX-v1 — universal chelate-isomer counter.

Forensik on the b00f9a0 voll-pool isocoverage (2026-06-04, n=10647) showed
that 1746 SMILES (16.4%) sit at 0% coverage with theoretical >= 2 isomers.
The dominant single contributor is "chelate counting only works for the
octahedron", which silently dropped TPR-6 trigonal-prism tris-chelates and
returned ``NotImplementedError`` for every other geometry.

This test module pins:

* env-OFF byte-identical to HEAD bcf56f8 (only octahedron supported,
  every other geometry raises ``NotImplementedError`` from
  ``count_chelate_isomers``; ``enumerate_chelate_configs("trigonal_prism", ...)``
  still raises ``KeyError`` because TPR-6 is not in ``_GEOM_KEY_TO_SHAPE``);
* env-ON (``DELFIN_FFFREE_POLYA_COVERAGE_FIX_v1=1``) unlocks the universal
  path with textbook chelate-isomer counts for SP-4, T-4, TBP-5, SPY-5,
  TPR-6, PBP-7, SAP-8, TTP-9, plus asymmetric-bidentate (e.g. glycinate)
  and mixed-monodentate-label support.

The two paths share no module-level state, so toggling the env-flag mid-
process is safe (it is read at every ``count_chelate_isomers`` call).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Ensure deterministic dict iteration / hash ordering even when pytest is
# invoked outside the wrapper script.
os.environ.setdefault("PYTHONHASHSEED", "0")

# Allow running the test file standalone without ``pip install -e .``.
_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from delfin.fffree import polya_isomer_count as PIC  # noqa: E402


_ENV_FLAG = "DELFIN_FFFREE_POLYA_COVERAGE_FIX_v1"


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """Each test starts with the env-flag cleared so default-OFF behaviour
    is verifiable; tests opt in by setting the flag themselves."""
    monkeypatch.delenv(_ENV_FLAG, raising=False)
    yield


# --------------------------------------------------------------------- #
# Group 1 — DEFAULT-OFF BYTE-IDENTICAL PARITY (HEAD bcf56f8)            #
# --------------------------------------------------------------------- #

def test_env_off_octahedron_unchanged():
    """The octahedron path is the only one HEAD supports; env-OFF results
    must match the textbook + previously-emitted values bit-exact."""
    assert PIC.count_chelate_isomers("octahedron", 1) == 1   # M(NN)X4
    assert PIC.count_chelate_isomers("octahedron", 2) == 3   # M(NN)2X2 (cis-Δ/Λ + trans)
    assert PIC.count_chelate_isomers("octahedron", 3) == 2   # M(NN)3 (Δ/Λ)


@pytest.mark.parametrize("geom", [
    "square_planar", "tetrahedron",
    "trigonal_bipyramid", "square_pyramid",
    "trigonal_prism", "pentagonal_bipyramid",
    "square_antiprism", "tricapped_trigonal_prism",
])
def test_env_off_nonoctahedron_raises(geom):
    """HEAD bcf56f8: ``count_chelate_isomers`` raises NotImplementedError
    for every non-octahedron geometry — the v1 fix must preserve this
    when its env-flag is unset."""
    with pytest.raises(NotImplementedError):
        PIC.count_chelate_isomers(geom, 1)


def test_env_off_trigonal_prism_enumerate_raises_keyerror():
    """The historical bug we are unlocking: TPR-6 silently fails because
    ``_GEOM_KEY_TO_SHAPE`` had no entry, ``_chelate_cis_edges`` falls back
    to ``_ANTIPODE_FULL`` and KeyErrors out.  Env-OFF must keep raising —
    that's the exact behaviour the FF-free builder (converter_backend
    line ~813) wraps in ``except Exception: return None``."""
    specs = [{"type": "B", "denticity": 2}] * 3
    with pytest.raises(KeyError):
        PIC.enumerate_chelate_configs("trigonal_prism", specs)


def test_env_off_other_enumerate_paths_still_work():
    """Non-TPR geometries with proper Pólya groups + ref_vectors should
    keep enumerating exactly as today (no behaviour change from v1)."""
    # OC-6 tris-chelate
    configs = PIC.enumerate_chelate_configs(
        "octahedron", [{"type": "B", "denticity": 2}] * 3,
    )
    assert len(configs) == 2
    # SP-4 mono-chelate + 2 mono X
    configs = PIC.enumerate_chelate_configs(
        "square_planar",
        [{"type": "B", "denticity": 2},
         {"type": "X", "denticity": 1},
         {"type": "X", "denticity": 1}],
    )
    assert len(configs) == 1


# --------------------------------------------------------------------- #
# Group 2 — ENV-ON TEXTBOOK CHELATE COUNTS                              #
# --------------------------------------------------------------------- #

def test_env_on_octahedron_path_preserved(monkeypatch):
    """Even with the flag on, the octahedron path retains its original
    edge-combinatorial algorithm (the universal helper is consulted only
    for non-octahedron geometries).  This guards against accidental
    divergence between the two implementations."""
    monkeypatch.setenv(_ENV_FLAG, "1")
    assert PIC.count_chelate_isomers("octahedron", 1) == 1
    assert PIC.count_chelate_isomers("octahedron", 2) == 3
    assert PIC.count_chelate_isomers("octahedron", 3) == 2


@pytest.mark.parametrize("geom,n_chelate,expected", [
    # Square planar M(NN)X2 — only cis placement possible.
    ("square_planar", 1, 1),
    # Tetrahedron — bite spans every edge equivalently.
    ("tetrahedron", 1, 1),                         # M(NN)X2
    ("tetrahedron", 2, 1),                         # M(NN)2
    # Trigonal bipyramid — the bite occupies an ax-eq edge (eq-eq is
    # disallowed at 120°), so M(NN)X3 collapses to a single orbit and
    # M(NN)2X must place the second chelate on the complementary ax-eq.
    ("trigonal_bipyramid", 1, 1),
    ("trigonal_bipyramid", 2, 2),
    # Square pyramid — apical-basal vs basal-basal: 2 distinct orbits.
    ("square_pyramid", 1, 2),
    # Trigonal prism tris-chelate — the case that silently dropped at
    # HEAD (KeyError into except Exception in converter_backend); the
    # eclipsed prism has Δ/Λ helical isomers exactly as OC-6 does.
    ("trigonal_prism", 3, 2),
    # Pentagonal bipyramid — eq-eq adjacent + ax-eq: 2 mono-chelate orbits.
    ("pentagonal_bipyramid", 1, 2),
    # Square antiprism — symmetry D4 makes 1 chelate give 3 orbits.
    ("square_antiprism", 1, 3),
])
def test_env_on_textbook_counts(monkeypatch, geom, n_chelate, expected):
    monkeypatch.setenv(_ENV_FLAG, "1")
    got = PIC.count_chelate_isomers(geom, n_chelate)
    assert got == expected, f"{geom} n={n_chelate}: got {got}, expected {expected}"


def test_env_on_asymmetric_bidentate_oc6(monkeypatch):
    """Glycinate-like M(N^O)3 with distinguishable N vs O arms.  Counting
    enantiomers as distinct: fac-Δ, fac-Λ, mer-Δ, mer-Λ = 4."""
    monkeypatch.setenv(_ENV_FLAG, "1")
    n = PIC.count_chelate_isomers_universal(
        "octahedron", 3, n_monodentate_fillers=0, asymmetric=True,
    )
    assert n == 4


def test_env_on_mixed_monodentate_labels_sp4(monkeypatch):
    """M(NN)(X)(Y) on square planar: chelate spans one cis-edge, X/Y on
    the other two vertices.  Up to D4, this gives a single orbit (X and
    Y are forced into the two remaining adjacent positions)."""
    monkeypatch.setenv(_ENV_FLAG, "1")
    n = PIC.count_chelate_isomers_universal(
        "square_planar", 1, n_monodentate_fillers=2,
        monodentate_labels=("X", "Y"),
    )
    assert n == 1


def test_env_on_mixed_monodentate_labels_oc6(monkeypatch):
    """M(NN)2(X)(Y) on octahedron: X-trans-Y vs X-cis-Y(Δ) vs X-cis-Y(Λ)
    = 3 orbits — same count as homo M(NN)2X2 because X<>Y disambiguation
    doesn't add orbits when the two chelates already pin Δ/Λ symmetry."""
    monkeypatch.setenv(_ENV_FLAG, "1")
    n = PIC.count_chelate_isomers_universal(
        "octahedron", 2, n_monodentate_fillers=2,
        monodentate_labels=("X", "Y"),
    )
    assert n == 3


# --------------------------------------------------------------------- #
# Group 3 — INVARIANTS                                                  #
# --------------------------------------------------------------------- #

def test_env_on_counts_match_enumerate_chelate_configs(monkeypatch):
    """The universal counter must agree with ``enumerate_chelate_configs``
    (single source of truth contract): every chelate count is the length
    of the corresponding enumerated config list."""
    monkeypatch.setenv(_ENV_FLAG, "1")
    cases = [
        ("square_planar", 1, 2, False),
        ("tetrahedron", 2, 0, False),
        ("trigonal_bipyramid", 2, 1, False),
        ("square_pyramid", 1, 3, False),
        ("trigonal_prism", 3, 0, False),
        ("octahedron", 3, 0, True),    # asym
    ]
    for geom, n_chel, n_mono, asym in cases:
        count = PIC.count_chelate_isomers_universal(
            geom, n_chel, n_mono, asymmetric=asym,
        )
        specs = []
        for _ in range(n_chel):
            specs.append({"type": "B", "denticity": 2, "asym": asym})
        for k in range(n_mono):
            specs.append({"type": "X", "denticity": 1})
        configs = PIC.enumerate_chelate_configs(geom, specs)
        assert count == len(configs), (
            f"{geom} n_chel={n_chel} n_mono={n_mono} asym={asym}: "
            f"count={count}, enumerate={len(configs)}"
        )


def test_env_on_deterministic_across_calls(monkeypatch):
    """The universal counter must give identical answers across repeated
    calls (no RNG, no hash-iteration order leak)."""
    monkeypatch.setenv(_ENV_FLAG, "1")
    for geom, n_chel in [
        ("square_planar", 1),
        ("trigonal_prism", 3),
        ("square_antiprism", 4),
        ("octahedron", 3),
    ]:
        first = PIC.count_chelate_isomers(geom, n_chel)
        for _ in range(4):
            assert PIC.count_chelate_isomers(geom, n_chel) == first


def test_env_on_input_validation(monkeypatch):
    """Universal counter must reject impossible combinations early."""
    monkeypatch.setenv(_ENV_FLAG, "1")
    with pytest.raises(ValueError):
        # 2 chelates + 3 mono = 7 != CN6 octahedron
        PIC.count_chelate_isomers_universal("octahedron", 2, 3)
    with pytest.raises(ValueError):
        PIC.count_chelate_isomers_universal("square_planar", -1, 6)
    with pytest.raises(ValueError):
        PIC.count_chelate_isomers_universal(
            "octahedron", 2, 2, monodentate_labels=("X",),  # length 1 != 2
        )


def test_list_supported_geometries_covers_core_polyhedra(monkeypatch):
    """The introspection helper must advertise every core polyhedron
    referenced by the FF-free builder, both env-OFF and env-ON (the list
    reads ``_GROUPS`` directly + the lazy f-block discovery, neither of
    which depends on the v1 flag)."""
    supported = set(PIC.list_supported_chelate_geometries())
    must_have = {
        "octahedron", "square_planar", "tetrahedron",
        "trigonal_bipyramid", "square_pyramid",
        "trigonal_prism", "pentagonal_bipyramid",
        "square_antiprism", "tricapped_trigonal_prism",
    }
    missing = must_have - supported
    assert not missing, f"missing geometries: {missing}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
