"""Tests for delfin/analysis_tools/conformer_clustering.py.

Covers:
- multi-XYZ parser (energies, malformed input, truncation)
- EDM-eigenvalue features (rotation/translation invariance, H exclusion)
- PCA + KMeans + silhouette scan (EnAn reference parity, degenerate cases)
- representative selection (lowest energy, tie -> smaller index)
- output writers (CSV/JSON/XYZ/PNG)
- CLI subcommand and CONTROL-driven pipeline hook
- direct parity against the archived Ensemble Analyzer reference package
  (optional: skipped when ``ensemble_analyzer`` is not importable)
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from delfin.analysis_tools.conformer_clustering import (  # noqa: E402
    ClusteringResult,
    ConformerRecord,
    ConformerClusteringError,
    FamilyInfo,
    cluster_ensemble,
    compute_edm_eigenvalues,
    compute_features,
    compute_quality_checks,
    parse_energy_from_comment,
    read_finalensemble_xyz,
    run_clustering,
    silhouette_scan,
    write_cluster_assignments_csv,
    write_cluster_summary_json,
    write_clustered_representatives_xyz,
)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

ELEM4 = ("C", "C", "H", "H")


def _write_ensemble(path: Path, frames) -> Path:
    """Write frames [(symbols, coords, energy), ...] as multi-XYZ."""
    lines = []
    for symbols, coords, energy in frames:
        lines.append(str(len(symbols)))
        lines.append("" if energy is None else f"{energy:.10f}")
        for sym, (x, y, z) in zip(symbols, coords):
            lines.append(f"{sym} {x:.8f} {y:.8f} {z:.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _two_family_frames(
    n_per_family: int = 10,
    separation: float = 5.0,
    seed: int = 0,
    symbols: tuple = ELEM4,
):
    """Synthetic ensemble: two well-separated Gaussian atom clouds."""
    rng = np.random.default_rng(seed)
    frames = []
    for family in range(2):
        center = family * separation
        for _ in range(n_per_family):
            coords = center + rng.normal(0.0, 0.05, (len(symbols), 3))
            energy = -100.0 + 0.5 * family + float(rng.uniform(0.0, 0.1))
            frames.append((list(symbols), coords, energy))
    return frames


def _records(frames):
    return [
        ConformerRecord(index=i, elements=list(sym), coords=np.asarray(c), comment="", energy_hartree=e)
        for i, (sym, c, e) in enumerate(frames)
    ]


# ---------------------------------------------------------------------------
# parser
# ---------------------------------------------------------------------------


class TestParser:
    def test_reads_frames_energies_and_comments(self, tmp_path):
        frames = [
            (("C", "H", "H"), np.array([[0.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0]]), -12.3456),
            (("C", "H", "H"), np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]), None),
        ]
        path = _write_ensemble(tmp_path / "ens.xyz", frames)
        path.write_text(
            path.read_text().replace(
                "-12.3456000000",
                "-12.3456 some conformer comment",
            )
        )
        confs = read_finalensemble_xyz(path)
        assert len(confs) == 2
        assert confs[0].index == 0
        assert confs[0].energy_hartree == pytest.approx(-12.3456)
        assert confs[1].energy_hartree is None
        assert "comment" in confs[0].comment
        assert confs[0].elements == ["C", "H", "H"]
        assert confs[0].coords.shape == (3, 3)

    def test_truncated_file_raises(self, tmp_path):
        path = tmp_path / "broken.xyz"
        path.write_text("5\ncomment\nC 0 0 0\n", encoding="utf-8")
        with pytest.raises(ConformerClusteringError, match="[Tt]runcated"):
            read_finalensemble_xyz(path)

    def test_malformed_atom_line_raises(self, tmp_path):
        path = tmp_path / "broken2.xyz"
        path.write_text("2\nc\nC 0 0 0\nN 1 2\n", encoding="utf-8")
        with pytest.raises(ConformerClusteringError, match="atom line"):
            read_finalensemble_xyz(path)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(ConformerClusteringError, match="not found"):
            read_finalensemble_xyz(tmp_path / "nope.xyz")

    def test_empty_file_raises(self, tmp_path):
        path = tmp_path / "empty.xyz"
        path.write_text("", encoding="utf-8")
        with pytest.raises(ConformerClusteringError):
            read_finalensemble_xyz(path)

    def test_goat_style_energy_line(self, tmp_path):
        path = tmp_path / "goat.xyz"
        path.write_text(
            "1\n-154.872431 Eh\nC 0.0 0.0 0.0\n",
            encoding="utf-8",
        )
        (confs,) = read_finalensemble_xyz(path)
        assert confs.energy_hartree == pytest.approx(-154.872431)


# ---------------------------------------------------------------------------
# EDM eigenvalue features
# ---------------------------------------------------------------------------


class TestEdmFeatures:
    def _geom(self):
        return np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

    def test_translation_invariance(self):
        g = self._geom()
        shifted = g + np.array([3.0, -2.0, 5.0])
        a = compute_edm_eigenvalues([("C", "C", "C")], [g], include_h=True)
        b = compute_edm_eigenvalues([("C", "C", "C")], [shifted], include_h=True)
        assert np.allclose(a, b, atol=1e-10)

    def test_rotation_invariance(self):
        g = self._geom()
        angle = np.deg2rad(37.0)
        rot = np.array(
            [[np.cos(angle), -np.sin(angle), 0.0],
             [np.sin(angle), np.cos(angle), 0.0],
             [0.0, 0.0, 1.0]]
        )
        a = compute_edm_eigenvalues([("C", "C", "C")], [g], include_h=True)
        b = compute_edm_eigenvalues([("C", "C", "C")], [g @ rot.T], include_h=True)
        assert np.allclose(a, b, atol=1e-10)

    def test_eigenvalues_match_numpy_reference(self):
        from scipy.spatial import distance_matrix

        g = self._geom()
        expected = np.linalg.eig(distance_matrix(g, g))[0].real
        got = compute_edm_eigenvalues([("C", "C", "C")], [g], include_h=True)[0]
        assert np.allclose(np.sort(got), np.sort(expected), atol=1e-10)

    def test_hydrogen_exclusion(self):
        g = np.array([[0.0, 0, 0], [1.0, 0, 0], [5.0, 0, 0], [5.0, 1.0, 0]])
        symbols = ["C", "H", "C", "H"]
        with_h = compute_edm_eigenvalues(symbols, [g], include_h=True)
        without_h = compute_edm_eigenvalues(symbols, [g], include_h=False)
        assert without_h.shape[1] == 2
        assert with_h.shape[1] == 4

    def test_feature_matrix_shape(self):
        frames = _two_family_frames(n_per_family=5)
        feats = compute_features(_records(frames), include_h=True)
        assert feats.shape == (10, 4)


# ---------------------------------------------------------------------------
# PCA / KMeans / silhouette scan
# ---------------------------------------------------------------------------


class TestClustering:
    def test_two_families_recovered(self):
        frames = _two_family_frames(n_per_family=10)
        result = cluster_ensemble(_records(frames), k="auto")
        assert result.k_chosen >= 2
        assert len(result.families) == result.k_chosen
        members = sorted(m for fam in result.families for m in fam.members)
        assert members == list(range(20))
        # both families present with disjoint members
        fam_sets = [frozenset(fam.members) for fam in result.families]
        assert len(set(fam_sets)) == 2

    def test_fixed_k(self):
        frames = _two_family_frames(n_per_family=10)
        result = cluster_ensemble(_records(frames), k=2)
        assert result.k_chosen == 2
        assert result.k_values == [2]

    def test_fixed_k_out_of_range(self):
        frames = _two_family_frames(n_per_family=3)
        with pytest.raises(ConformerClusteringError, match="out of range"):
            cluster_ensemble(_records(frames), k=6)

    def test_n_equals_1(self):
        frames = _two_family_frames(n_per_family=1)
        result = cluster_ensemble(_records([frames[0]]), k="auto")
        assert result.k_chosen == 1
        assert len(result.families) == 1
        assert any("one conformer" in w.lower() for w in result.warnings)

    def test_identical_conformers_single_family(self):
        frame = _two_family_frames(n_per_family=1)[0]
        frames = [frame] * 6
        result = cluster_ensemble(_records(frames), k="auto")
        assert len(result.families) == 1
        assert any("degenerate" in w.lower() for w in result.warnings)

    def test_unknown_method_rejected(self):
        frames = _two_family_frames(n_per_family=2)
        with pytest.raises(ConformerClusteringError, match="method"):
            cluster_ensemble(_records(frames), k=2, method="tfd")

    def test_determinism_same_seed(self):
        frames = _two_family_frames(n_per_family=8, seed=3)
        recs = _records(frames)
        r1 = cluster_ensemble(recs, k="auto", seed=42)
        r2 = cluster_ensemble(recs, k="auto", seed=42)
        assert r1.k_chosen == r2.k_chosen
        assert [fam.members for fam in r1.families] == [fam.members for fam in r2.families]


# ---------------------------------------------------------------------------
# representative selection
# ---------------------------------------------------------------------------


class TestRepresentatives:
    def test_lowest_energy_is_representative(self):
        frames = _two_family_frames(n_per_family=6)
        result = cluster_ensemble(_records(frames), k=2)
        for fam in result.families:
            member_energies = [(result.energies_hartree[m], m) for m in fam.members]
            best_energy, best_idx = min(member_energies)
            assert fam.representative_index == best_idx

    def test_energy_tie_breaks_to_smaller_index(self):
        # Two families at identical geometry; every member of family A has
        # energy -10.0, family B -5.0.  Within A all energies are equal, so
        # the representative must be the smallest original index.
        coords = np.array([[0.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])
        frames = []
        for family, energy in ((0, -10.0), (1, -5.0)):
            for _ in range(4):
                frames.append((["C"] * 4, coords + family * 10.0, energy))
        result = cluster_ensemble(_records(frames), k=2)
        for fam in result.families:
            idx = fam.members
            assert fam.representative_index == min(idx)


# ---------------------------------------------------------------------------
# silhouette scan boundaries (reference: min_k = max(int(N*0.1), 2),
# max_k = int(N*0.8), degenerate -> min_k with score 0.0)
# ---------------------------------------------------------------------------


class TestSilhouetteScan:
    def test_k_bounds(self):
        rng = np.random.default_rng(1)
        feats = rng.normal(0, 1, (30, 3))
        k_chosen, k_values, scores, _ = silhouette_scan(feats, seed=42)
        # min_k = max(int(30*0.1), 2) = 3 ; max_k = int(30*0.8) = 24
        assert k_values[0] == 3
        assert k_values[-1] == 24
        assert k_chosen in k_values
        assert len(scores) == len(k_values)

    def test_small_ensemble_uses_min_k(self):
        # N=4: min_k = 2, max_k = int(3.2) = 3 -> normal scan over [2, 3]
        rng = np.random.default_rng(2)
        feats = rng.normal(0, 1, (4, 2))
        k_chosen, k_values, _, _ = silhouette_scan(feats, seed=42)
        assert k_values[0] == 2
        assert k_chosen in k_values

    def test_degenerate_scan_returns_min_k_zero_score(self):
        # N=3: min_k = 2, max_k = int(2.4) = 2 < ... no: 2 >= 2 -> normal scan.
        # N=2: min_k = 2, max_k = int(1.6) = 1 < 2 -> degenerate path.
        feats = np.array([[0.0, 0.0], [1.0, 1.0]])
        k_chosen, k_values, scores, _ = silhouette_scan(feats, seed=42)
        assert k_chosen == 2
        assert k_values == [2]
        assert scores == [0.0]


# ---------------------------------------------------------------------------
# output writers
# ---------------------------------------------------------------------------


class TestOutputs:
    def _run(self, tmp_path):
        frames = _two_family_frames(n_per_family=8)
        path = _write_ensemble(tmp_path / "ens.xyz", frames)
        return run_clustering(path, out_dir=tmp_path / "out", k="auto", write_plots=False)

    def test_all_files_written(self, tmp_path):
        summary = self._run(tmp_path)
        written = {Path(p).name for p in summary["files"].values()}
        assert {
            "clustered_representatives.xyz",
            "cluster_assignments.csv",
            "cluster_summary.json",
            "silhouette_scan.csv",
            "pca_coordinates.csv",
        } <= written

    def test_assignments_csv_content(self, tmp_path):
        summary = self._run(tmp_path)
        with open(summary["files"]["cluster_assignments.csv"], newline="") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 16
        assert {int(r["cluster"]) for r in rows} == {1, 2}
        assert sum(int(r["is_representative"]) for r in rows) == 2
        # relative energies: the global minimum must be 0.0
        assert min(float(r["relative_energy_kcal_mol"]) for r in rows) == pytest.approx(0.0, abs=1e-9)

    def test_representatives_xyz_structure(self, tmp_path):
        summary = self._run(tmp_path)
        text = Path(summary["files"]["clustered_representatives.xyz"]).read_text()
        # One frame per representative: each frame is "<natoms>\n<comment>\n"
        # followed by exactly <natoms> atom lines.  Count frames via the
        # digit-only headers and the total atom lines via those headers.
        lines = text.splitlines()
        headers = [
            i for i, line in enumerate(lines)
            if line.strip().isdigit()
        ]
        assert len(headers) == 2  # one frame per family
        n_atoms_from_headers = sum(int(lines[i].strip()) for i in headers)
        # Atom lines are the non-header, non-comment lines; the comment of
        # every frame is the line right after a header.
        comment_idx = {i + 1 for i in headers}
        atom_line_count = sum(
            1 for i, line in enumerate(lines)
            if i not in set(headers) and i not in comment_idx and line.strip()
        )
        assert atom_line_count == n_atoms_from_headers
        assert atom_line_count == 8  # 2 reps x 4 atoms

    def test_summary_json_fields(self, tmp_path):
        summary = self._run(tmp_path)
        data = json.loads(Path(summary["files"]["cluster_summary.json"]).read_text())
        for key in (
            "method", "n_conformers", "n_families", "k_scan",
            "pca", "kmeans_params", "reference", "random_seed",
        ):
            assert key in data
        assert data["method"] == "EnAn-style EDM-eigenvalue PCA + KMeans"
        assert data["kmeans_params"]["n_init"] == "auto"
        assert data["random_seed"] == 42

    def test_k_fixed_summary_has_no_scan(self, tmp_path):
        path = _write_ensemble(tmp_path / "ens2.xyz", _two_family_frames(n_per_family=8))
        summary = run_clustering(path, out_dir=tmp_path / "out2", k=3, write_plots=False)
        data = json.loads(Path(summary["files"]["cluster_summary.json"]).read_text())
        assert data["k_scan"]["mode"] == "fixed_or_degenerate"
        assert data["k_scan"]["k_chosen"] == 3


# ---------------------------------------------------------------------------
# CLI + pipeline hook
# ---------------------------------------------------------------------------


class TestCliAndHook:
    def test_cli_subcommand(self, tmp_path, capsys):
        path = _write_ensemble(tmp_path / "ens.xyz", _two_family_frames(n_per_family=8))
        from delfin.analysis_tools.conformer_clustering import run_cli

        rc = run_cli([str(path), "--out", str(tmp_path / "cli_out"), "-k", "2"])
        assert rc == 0
        out = capsys.readouterr().out
        assert "Selected families:    2 (k=2)" in out
        assert "All 2 representatives written successfully." in out
        assert (tmp_path / "cli_out" / "cluster_summary.json").exists()

    def test_cli_invalid_k_errors_cleanly(self, tmp_path, capsys):
        path = _write_ensemble(tmp_path / "ens.xyz", _two_family_frames(n_per_family=4))
        from delfin.analysis_tools.conformer_clustering import run_cli

        rc = run_cli([str(path), "--out", str(tmp_path / "cli_out"), "-k", "1"])
        assert rc != 0
        assert "out of range" in capsys.readouterr().err

    def test_maybe_run_from_config_disabled(self, tmp_path):
        from delfin.analysis_tools.conformer_clustering import maybe_run_from_config

        assert maybe_run_from_config({"conformer_clustering": "no"}, tmp_path) is None

    def test_maybe_run_from_config_clusters_finalensemble(self, tmp_path, caplog):
        import logging

        from delfin.analysis_tools.conformer_clustering import maybe_run_from_config

        goat_dir = tmp_path / "def2-SVP_GOAT"
        goat_dir.mkdir()
        _write_ensemble(goat_dir / "molecule.finalensemble.xyz", _two_family_frames(n_per_family=8))
        config = {
            "conformer_clustering": "yes",
            "conformer_clustering_method": "enan",
            "conformer_clusters": "auto",
            "conformer_cluster_exclude_h": "no",
        }
        with caplog.at_level(logging.INFO):
            summary = maybe_run_from_config(config, tmp_path, logger_=logging.getLogger(__name__))
        assert summary is not None
        assert summary["n_families"] >= 2
        assert (goat_dir / "conformer_clustering" / "cluster_summary.json").exists()

    def test_maybe_run_from_config_extended_scans_k2(self, tmp_path, caplog):
        import json
        import logging

        from delfin.analysis_tools.conformer_clustering import maybe_run_from_config

        goat_dir = tmp_path / "def2-SVP_GOAT"
        goat_dir.mkdir()
        _write_ensemble(goat_dir / "molecule.finalensemble.xyz", _two_family_frames(n_per_family=18))
        config = {
            "conformer_clustering": "yes",
            "conformer_clustering_method": "enan",
            "conformer_clusters": "auto",
            "conformer_cluster_exclude_h": "no",
        }
        # reference (extended off): scan starts at max(int(0.1*n),2)=3 for n=36
        with caplog.at_level(logging.INFO):
            ref = maybe_run_from_config(config, tmp_path, logger_=logging.getLogger(__name__))
        ref_js = json.loads((goat_dir / "conformer_clustering" / "cluster_summary.json").read_text())
        assert ref_js["k_scan"]["k_values"][0] == 3

        # extended on: k=2 is scanned as well
        config["conformer_cluster_extended"] = "yes"
        with caplog.at_level(logging.INFO):
            ext = maybe_run_from_config(config, tmp_path, logger_=logging.getLogger(__name__))
        ext_js = json.loads((goat_dir / "conformer_clustering" / "cluster_summary.json").read_text())
        assert ext_js["k_scan"]["extended"] is True
        assert ext_js["k_scan"]["k_values"][0] == 2

    def test_maybe_run_from_config_no_ensemble_is_silent(self, tmp_path, caplog):
        import logging

        from delfin.analysis_tools.conformer_clustering import maybe_run_from_config

        config = {"conformer_clustering": "yes"}
        with caplog.at_level(logging.INFO):
            summary = maybe_run_from_config(config, tmp_path, logger_=logging.getLogger(__name__))
        assert summary is None

    def test_hook_after_goat_call_sites(self):
        # The two GOAT call sites in the pipeline must be followed by the
        # clustering hook (task requirement: run directly on the produced
        # finalensemble.xyz).
        from pathlib import Path as _Path

        src = (_Path(__file__).resolve().parent.parent / "delfin" / "workflows" / "pipeline.py").read_text()
        assert src.count("maybe_run_from_config(") >= 2


# ---------------------------------------------------------------------------
# EnAn reference parity (optional import)
# ---------------------------------------------------------------------------


class TestEnAnParity:
    """Direct comparison against the archived Ensemble Analyzer package.

    Skipped when the reference is not importable; the archived source is
    pinned in ``reference/enan/ensemble_analyzer-main`` (Zenodo DOI
    10.5281/zenodo.18255912, GitHub andre-cloud/ensemble_analyzer).
    """

    @pytest.fixture(autouse=True)
    def _require_reference(self):
        """Skip parity tests unless the archived reference submodules import.

        Unlike ``pytest.importorskip("ensemble_analyzer")`` at class level —
        which succeeds for the bare namespace package even when its source has
        been deleted and would here wipe the *whole module's* collection — this
        guards the deep submodule import. The failure mode we guard against is
        real: the reference source tree was emptied once (only ``__pycache__``
        left), and a bare top-level import cannot detect that.
        """
        pytest.importorskip("ensemble_analyzer._clustering.cluster_manager")
        pytest.importorskip("ensemble_analyzer._clustering.cluster_config")

    def test_eigenvalue_features_match_reference(self):
        import inspect

        from ensemble_analyzer._clustering.cluster_manager import ClusteringManager

        rng = np.random.default_rng(7)
        geoms = [rng.normal(0, 1, (5, 3)) for _ in range(4)]
        symbols = ["C", "N", "O", "H", "H"]
        # call the reference routine directly
        manager = ClusteringManager.__new__(ClusteringManager)  # no logger needed
        ref = manager._calculate_distance_matrix_eigenvalues(
            geometries=np.array(geoms), atoms=np.array(symbols), include_H=True
        )
        ours = compute_edm_eigenvalues(symbols, geoms, include_h=True)
        assert np.allclose(ref, ours, atol=1e-10)

        ref_no_h = manager._calculate_distance_matrix_eigenvalues(
            geometries=np.array(geoms), atoms=np.array(symbols), include_H=False
        )
        ours_no_h = compute_edm_eigenvalues(symbols, geoms, include_h=False)
        assert np.allclose(ref_no_h, ours_no_h, atol=1e-10)

    def test_silhouette_scan_matches_reference(self):
        from unittest.mock import MagicMock

        from ensemble_analyzer._clustering.cluster_manager import ClusteringManager

        rng = np.random.default_rng(11)
        feats = np.vstack([rng.normal(c, 0.1, (6, 2)) for c in (0.0, 3.0, 6.0)])
        manager = ClusteringManager(
            logger=MagicMock(),
            config=type("C", (), {"random_state": 42, "include_H": True})(),
        )
        ref_k, ref_range, ref_scores = manager.find_optimal_clusters(feats)
        our_k, our_values, our_scores, _ = silhouette_scan(feats, seed=42)
        assert our_k == ref_k
        assert our_values == list(ref_range)
        for ours, theirs in zip(our_scores, ref_scores):
            assert ours == pytest.approx(theirs, abs=1e-9)

    def test_full_pipeline_matches_reference(self):
        from unittest.mock import MagicMock

        from ensemble_analyzer._clustering.cluster_config import ClusteringConfig
        from ensemble_analyzer._clustering.cluster_manager import ClusteringManager
        from ensemble_analyzer.conformer.conformer import Conformer
        from ensemble_analyzer.conformer.energy_data import EnergyRecord, EnergyStore

        rng = np.random.default_rng(13)
        n_per_family = 12
        ref_confs = []
        our_frames = []
        for family, center in enumerate((0.0, 4.0)):
            for i in range(n_per_family):
                geom = center + rng.normal(0, 0.05, (4, 3))
                energy = -100.0 - family * 2.0 - i * 0.01
                ref_conf = Conformer(number=i, geom=geom, atoms=("C", "C", "H", "H"), raw=True)
                store = EnergyStore()
                store.add(1, EnergyRecord(E=energy))
                ref_conf.energies = store
                ref_confs.append(ref_conf)
                our_frames.append((("C", "C", "H", "H"), geom, energy))

        manager = ClusteringManager(logger=MagicMock(), config=ClusteringConfig(include_H=True))
        ref_result = manager._execute_pca_pipeline(ref_confs, n_clusters=None, include_H=True)
        labels_ref = np.asarray(ref_result.clusters)

        ours = cluster_ensemble(
            [
                ConformerRecord(
                    index=i,
                    elements=list(our_frames[i][0]),
                    coords=np.asarray(our_frames[i][1]),
                    comment="",
                    energy_hartree=our_frames[i][2],
                )
                for i in range(len(our_frames))
            ],
            k="auto",
        )
        labels_ours = np.zeros(len(our_frames), dtype=int)
        for fam in ours.families:
            for m in fam.members:
                labels_ours[m] = fam.family_id
        # same k and an identical partition up to label permutation.
        # NB: sort by the sorted element list — sorted() on frozensets uses
        # the subset relation and gives an unstable order for disjoint sets.
        ref_sets = sorted(
            (sorted(np.where(labels_ref == c)[0].tolist()) for c in set(labels_ref))
        )
        our_sets = sorted(
            (sorted(np.where(labels_ours == c)[0].tolist()) for c in set(labels_ours))
        )
        assert ref_sets == our_sets


# ---------------------------------------------------------------------------
# New: energy-parser robustness, extended k-scan, representative roundtrip,
# quality checks (incl. stereochemistry), seed stability, missing energies.
# ---------------------------------------------------------------------------


def _chiral_geometry():
    """A genuinely chiral 5-atom geometry (mirror NOT superposable by rotation)."""
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.2, 0.0],
            [-0.3, -0.4, 1.0],
            [0.2, -0.3, -0.9],
        ],
        dtype=float,
    )


class TestEnergyParserExtended:
    def test_energy_field_wins_over_ids(self):
        # Cluster/conformer ids must never be parsed as the energy.
        assert parse_energy_from_comment(
            "Family 1 | conformer 0 | E=-58.2220784074 Eh"
        ) == pytest.approx(-58.2220784074)
        assert parse_energy_from_comment(
            "E=-58.2220784074 Eh family=1 conformer=0"
        ) == pytest.approx(-58.2220784074)

    def test_goat_raw_lead_number(self):
        # GOAT writes the energy as the first number.
        assert parse_energy_from_comment(
            "-62.9834912942 converged=true"
        ) == pytest.approx(-62.9834912942)

    def test_id_only_line_is_none(self):
        assert parse_energy_from_comment("Family 3 conformer 7 no energy") is None
        assert parse_energy_from_comment("family=2 conformer=5") is None
        assert parse_energy_from_comment("") is None

    def test_energy_field_variant(self):
        assert parse_energy_from_comment("energy: -12.3456") == pytest.approx(-12.3456)
        assert parse_energy_from_comment("Energy = -1.0") == pytest.approx(-1.0)


class TestRepresentativeRoundtrip:
    def test_exactly_one_rep_per_family(self, tmp_path):
        path = _write_ensemble(tmp_path / "ens.xyz", _two_family_frames(n_per_family=6))
        summary = run_clustering(path, out_dir=tmp_path / "out", write_plots=False)
        reps = [fam.representative_index for fam in summary["result"].families]
        # one representative per family, all distinct
        assert len(reps) == summary["n_families"]
        assert len(set(reps)) == len(reps)
        # every representative is a true member of its own family
        for fam in summary["result"].families:
            assert fam.representative_index in fam.members

    def test_roundtrip_representatives_xyz(self, tmp_path):
        path = _write_ensemble(tmp_path / "ens.xyz", _two_family_frames(n_per_family=6))
        out_dir = tmp_path / "out"
        summary = run_clustering(path, out_dir=out_dir, write_plots=False)
        rep_path = Path(summary["files"]["clustered_representatives.xyz"])
        text = rep_path.read_text()
        reps = read_finalensemble_xyz(rep_path)
        # one frame per family
        assert len(reps) == summary["n_families"]
        # every frame documents its ORIGINAL conformer id and family id
        for rec in reps:
            assert rec.conformer_id is not None
            assert rec.family_id is not None
        # energies + original ids round-trip against the in-memory results
        expected = {
            fam.representative_index: (
                summary["result"].energies_hartree[fam.representative_index],
                fam.family_id,
            )
            for fam in summary["result"].families
        }
        for rec in reps:
            orig_energy, fam_id = expected[rec.conformer_id]
            assert rec.family_id == fam_id
            if orig_energy is not None:
                assert rec.energy_hartree == pytest.approx(orig_energy, abs=1e-9)
            else:
                assert rec.energy_hartree is None
        # no representative is silently dropped
        assert sorted(rec.conformer_id for rec in reps) == sorted(expected)

    def test_original_file_not_modified(self, tmp_path):
        path = _write_ensemble(tmp_path / "ens.xyz", _two_family_frames(n_per_family=6))
        original = path.read_text()
        run_clustering(path, out_dir=tmp_path / "out", write_plots=False)
        assert path.read_text() == original


class TestExtendedKScan:
    def _feats_36(self):
        rng = np.random.default_rng(3)
        return np.vstack([rng.normal(c, 0.1, (12, 2)) for c in (0.0, 3.0, 6.0)])

    def test_extended_includes_k2(self):
        feats = self._feats_36()
        _, ref_vals, _, _ = silhouette_scan(feats, seed=42)
        _, ext_vals, _, ext_w = silhouette_scan(feats, seed=42, extended=True)
        assert ref_vals[0] == 3  # reference starts at k=3 for n=36
        assert ext_vals[0] == 2  # extended also tests k=2
        assert any("extended" in w.lower() for w in ext_w)

    def test_reference_mode_unchanged(self):
        feats = self._feats_36()
        k_ref, vals_ref, scores_ref, _ = silhouette_scan(feats, seed=42)
        # default (extended=False) is exactly the reference range
        assert vals_ref[0] == max(int(36 * 0.1), 2) == 3
        assert all(s is None or isinstance(s, float) for s in scores_ref)
        assert k_ref in vals_ref


class TestQualityChecks:
    def _records_from_frames(self, frames):
        return _records(frames)

    def test_mirror_twin_flagged(self):
        # A conformer that is the mirror image of the family representative
        # passes the EDM-eigenvalue (mirror-insensitive) descriptor, so it
        # can land in the same family; compute_quality_checks must flag it.
        tet = _chiral_geometry()
        mirror = _mirror_of(tet)
        conformers = _records(
            [
                (["C", "C", "C", "C", "C"], tet, -100.0),
                (["C", "C", "C", "C", "C"], tet, -100.05),
                (["C", "C", "C", "C", "C"], mirror, -100.10),
            ]
        )
        # Realistic single-family result: EDM eigenvalues cannot tell the
        # mirror image apart, so all three share one family with member 0
        # (lowest energy) as representative.
        result = ClusteringResult(
            method="enan",
            input_file="",
            n_conformers=3,
            include_h=True,
            feature_dimension=5,
            eigenvalues=np.zeros((3, 5)),
            pca_scores=np.zeros((3, 2)),
            explained_variance_ratio=np.zeros(2),
            k_values=[1],
            silhouette_scores=[None],
            k_chosen=1,
            kmeans_params={"n_init": "auto", "random_state": 42},
            seed=42,
            labels=np.zeros(3, dtype=int),
            families=[
                FamilyInfo(
                    family_id=1,
                    representative_index=0,
                    members=[0, 1, 2],
                    rep_energy_hartree=-100.0,
                    rep_rel_energy_kcal_mol=0.0,
                )
            ],
            energies_hartree=[-100.0, -100.05, -100.10],
        )
        checks = compute_quality_checks(result, conformers)
        # the mirror twin member (index 2) must be flagged in its family
        flagged = [
            i
            for fam in checks["intra_family"]
            for i in fam["mirror_twin_member_indices"]
        ]
        assert 2 in flagged
        assert any("mirror" in w.lower() for w in checks["warnings"])

    def test_inhomogeneous_family_warning(self):
        # A member with a genuinely DIFFERENT shape (not a pure translation,
        # which Kabsch alignment makes invisible) far from the representative
        # -> homogeneity warning.
        base = _chiral_geometry()
        # structurally distinct: elongated chain geometry, not a shifted copy
        far = np.array(
            [
                [0.0, 0.0, 0.0],
                [6.0, 0.0, 0.0],
                [12.0, 0.0, 0.0],
                [18.0, 1.0, 0.0],
                [24.0, 0.0, 1.0],
            ],
            dtype=float,
        )
        conformers = _records(
            [
                (["C", "C", "C", "C", "C"], base, -100.0),
                (["C", "C", "C", "C", "C"], far, -99.0),
            ]
        )
        result = ClusteringResult(
            method="enan",
            input_file="",
            n_conformers=2,
            include_h=True,
            feature_dimension=5,
            eigenvalues=np.zeros((2, 5)),
            pca_scores=np.zeros((2, 2)),
            explained_variance_ratio=np.zeros(2),
            k_values=[1],
            silhouette_scores=[None],
            k_chosen=1,
            kmeans_params={"n_init": "auto", "random_state": 42},
            seed=42,
            labels=np.zeros(2, dtype=int),
            families=[
                FamilyInfo(
                    family_id=1,
                    representative_index=0,
                    members=[0, 1],
                    rep_energy_hartree=-100.0,
                    rep_rel_energy_kcal_mol=0.0,
                )
            ],
            energies_hartree=[-100.0, -99.0],
        )
        checks = compute_quality_checks(result, conformers)
        max_rmsd = checks["intra_family"][0]["rmsd_to_rep_A_max"]
        assert max_rmsd is not None and max_rmsd > 2.0  # far member detected
        assert any("exceeds 2.0" in w for w in checks["warnings"])


class TestStabilityAndMissing:
    def test_same_seed_stable(self):
        rng = np.random.default_rng(5)
        frames = []
        for fam, center in enumerate((0.0, 4.0)):
            for i in range(12):
                coords = center + rng.normal(0, 0.05, (4, 3))
                frames.append((list("CCHH"), coords, -100.0 + fam))
        a = cluster_ensemble(_records(frames), k="auto", seed=42)
        b = cluster_ensemble(_records(frames), k="auto", seed=42)
        # identical partition up to family numbering
        set_a = {frozenset(f.members) for f in a.families}
        set_b = {frozenset(f.members) for f in b.families}
        assert set_a == set_b

    def test_missing_energy_handled(self, tmp_path):
        frames_wo = [
            (list("CCHH"), np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float), None),
            (list("CCHH"), np.array([[0.1, 0, 0], [1.1, 0, 0], [0.1, 1, 0], [0.1, 0, 1]], float), None),
        ]
        path = _write_ensemble(tmp_path / "ens.xyz", frames_wo)
        # frames with None energy -> write nothing in the energy field
        text = path.read_text()
        assert "None" not in text
        summary = run_clustering(path, out_dir=tmp_path / "out", write_plots=False)
        assert summary["n_families"] >= 1
        # representatives either carry an energy or are documented as unknown
        for fam in summary["result"].families:
            assert fam.representative_index in fam.members

    def test_translation_vs_distinct_geometries(self):
        # Pure translation of the SAME geometry yields identical EDM
        # eigenvalues (that is why shifted copies do not validate distinct
        # structure classes).  A genuinely different shape must differ.
        base = np.array([[0.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0]])
        shifted = base + np.array([3.0, -1.0, 2.0])
        dist = np.array([[0.0, 0, 0], [2.0, 0, 0], [0, 0, 3.0]])  # different shape
        same = compute_features(_records([
            (["C", "C", "C"], base, None),
            (["C", "C", "C"], shifted, None),
        ]))
        assert np.allclose(same[0], same[1], atol=1e-10)  # translation-invariant
        diff = compute_features(_records([
            (["C", "C", "C"], base, None),
            (["C", "C", "C"], dist, None),
        ]))
        assert not np.allclose(diff[0], diff[1], atol=1e-6)  # shape change detectable


def _mirror_of(coords):
    return np.asarray(coords, dtype=float) * np.array([-1.0, 1.0, 1.0])
