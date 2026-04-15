"""Dedup + energy-window filter tests for GUPPY (M2)."""
from typing import List

from delfin import guppy_sampling
from delfin.guppy_sampling import RunResult


def _mk_result(energy: float, coords: List[str], run_idx: int = 1,
               label: str = "isomer", source: str = "isomer") -> RunResult:
    # RunResult = (energy, natoms, coords, run_idx, label, source)
    return (energy, len(coords), coords, run_idx, label, source)


def _tetrahedron(scale: float = 1.0) -> List[str]:
    return [
        f"C {scale * 0.000:.6f} {scale * 0.000:.6f} {scale * 0.000:.6f}",
        f"C {scale * 1.000:.6f} {scale * 1.000:.6f} {scale * 0.000:.6f}",
        f"C {scale * 1.000:.6f} {scale * 0.000:.6f} {scale * 1.000:.6f}",
        f"C {scale * 0.000:.6f} {scale * 1.000:.6f} {scale * 1.000:.6f}",
    ]


def _tetrahedron_translated(dx: float = 10.0) -> List[str]:
    """Translation must not register as a new structure (rotation invariance)."""
    return [
        f"C {dx + 0.0:.6f} 0.0 0.0",
        f"C {dx + 1.0:.6f} 1.0 0.0",
        f"C {dx + 1.0:.6f} 0.0 1.0",
        f"C {dx + 0.0:.6f} 1.0 1.0",
    ]


def test_empty_input_returns_empty() -> None:
    assert guppy_sampling._filter_xtb_candidates(
        [], rmsd_cutoff=0.3, energy_window_eh=0.04
    ) == []


def test_energy_window_drops_high_energy() -> None:
    results = [
        _mk_result(-10.0, _tetrahedron(), run_idx=1),
        _mk_result(-9.95, _tetrahedron(scale=1.3), run_idx=2),  # ~31 kcal/mol above min
    ]
    kept = guppy_sampling._filter_xtb_candidates(
        results, rmsd_cutoff=0.0, energy_window_eh=25.0 / 627.509474
    )
    assert len(kept) == 1
    assert kept[0][3] == 1


def test_translation_treated_as_duplicate() -> None:
    """A translated copy is identical by sorted-pair-distance fingerprint."""
    results = [
        _mk_result(-10.00, _tetrahedron(), run_idx=1),
        _mk_result(-9.99, _tetrahedron_translated(dx=25.0), run_idx=2),
    ]
    kept = guppy_sampling._filter_xtb_candidates(
        results, rmsd_cutoff=0.3, energy_window_eh=0.0
    )
    assert len(kept) == 1
    assert kept[0][3] == 1  # lowest-energy wins


def test_distinct_structures_both_kept() -> None:
    results = [
        _mk_result(-10.00, _tetrahedron(scale=1.0), run_idx=1),
        _mk_result(-9.95, _tetrahedron(scale=1.5), run_idx=2),  # different bond lengths
    ]
    kept = guppy_sampling._filter_xtb_candidates(
        results, rmsd_cutoff=0.1, energy_window_eh=0.0
    )
    assert len(kept) == 2


def test_lowest_energy_wins_within_cluster() -> None:
    results = [
        _mk_result(-10.05, _tetrahedron(), run_idx=3),
        _mk_result(-10.00, _tetrahedron_translated(dx=5.0), run_idx=1),
        _mk_result(-9.90, _tetrahedron_translated(dx=10.0), run_idx=2),
    ]
    # Pre-sort as run_sampling does before calling the filter.
    results.sort(key=lambda r: (r[0], r[3]))
    kept = guppy_sampling._filter_xtb_candidates(
        results, rmsd_cutoff=0.3, energy_window_eh=0.0
    )
    assert len(kept) == 1
    assert kept[0][3] == 3  # lowest energy


def test_sorted_pair_distances_ignores_hydrogens() -> None:
    coords_with_h = [
        "C 0.0 0.0 0.0",
        "C 1.0 0.0 0.0",
        "H 2.0 0.0 0.0",
        "H -1.0 0.0 0.0",
    ]
    heavy = guppy_sampling._extract_heavy_atom_coords(coords_with_h)
    assert len(heavy) == 2
    assert all(sym == "C" for sym, *_ in heavy)
    fp = guppy_sampling._sorted_pair_distances(heavy)
    assert len(fp) == 1
    assert abs(fp[0] - 1.0) < 1e-9


def test_zero_cutoffs_disable_filtering() -> None:
    results = [
        _mk_result(-10.0, _tetrahedron(), run_idx=1),
        _mk_result(-9.9, _tetrahedron_translated(), run_idx=2),
    ]
    kept = guppy_sampling._filter_xtb_candidates(
        results, rmsd_cutoff=0.0, energy_window_eh=0.0,
        constitution_dedup=False,
    )
    assert len(kept) == 2


# ---------------------------------------------------------------------------
# M2b: constitution-fingerprint dedup (conformers of same constitution collapse)
# ---------------------------------------------------------------------------

def test_constitution_dedup_collapses_same_fingerprint(monkeypatch) -> None:
    """Two candidates with identical constitution FP → only lowest-energy kept."""
    # We stub the fingerprint helper so the test doesn't need RDKit / metal mol.
    shared_fp = ("metal_Fe", ("N", "N"), ("cis", "trans"))
    fps = iter([shared_fp, shared_fp, ("different_fp",)])

    def fake_fp(mol_template, coords):
        return next(fps)

    monkeypatch.setattr(
        guppy_sampling, "_constitution_fingerprint_from_coords", fake_fp
    )

    results = [
        _mk_result(-10.00, _tetrahedron(), run_idx=1, source="isomer"),
        _mk_result(-9.95, _tetrahedron(scale=1.2), run_idx=2, source="isomer"),
        _mk_result(-9.90, _tetrahedron(scale=1.4), run_idx=3, source="isomer"),
    ]
    kept = guppy_sampling._filter_xtb_candidates(
        results, rmsd_cutoff=0.3, energy_window_eh=0.0,
        mol_template=object(),  # sentinel, fake_fp ignores it
        constitution_dedup=True,
    )
    # Two distinct fingerprints -> two kept (runs 1 and 3); run 2 collapsed.
    assert [r[3] for r in kept] == [1, 3]


def test_constitution_dedup_falls_back_to_rmsd_when_fp_none(monkeypatch) -> None:
    """When fingerprint unavailable, geometric RMSD is used as fallback."""
    monkeypatch.setattr(
        guppy_sampling, "_constitution_fingerprint_from_coords",
        lambda mt, c: None,
    )

    results = [
        _mk_result(-10.00, _tetrahedron(), run_idx=1),
        _mk_result(-9.99, _tetrahedron_translated(dx=25.0), run_idx=2),
    ]
    kept = guppy_sampling._filter_xtb_candidates(
        results, rmsd_cutoff=0.3, energy_window_eh=0.0,
        mol_template=None, constitution_dedup=True,
    )
    assert len(kept) == 1
    assert kept[0][3] == 1


def test_constitution_dedup_logs_unique_count(monkeypatch, caplog) -> None:
    import logging
    fps = iter([("A",), ("A",), ("B",)])
    monkeypatch.setattr(
        guppy_sampling, "_constitution_fingerprint_from_coords",
        lambda mt, c: next(fps),
    )
    results = [
        _mk_result(-10.0, _tetrahedron(), run_idx=1),
        _mk_result(-9.9, _tetrahedron(), run_idx=2),
        _mk_result(-9.8, _tetrahedron(), run_idx=3),
    ]
    caplog.set_level(logging.INFO)
    guppy_sampling._filter_xtb_candidates(
        results, rmsd_cutoff=0.3, energy_window_eh=0.0,
        mol_template=object(), constitution_dedup=True,
    )
    messages = "\n".join(rec.getMessage() for rec in caplog.records)
    assert "unique constitutions: 2" in messages
    assert "same-constitution" in messages
