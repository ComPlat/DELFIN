"""Unit tests for :mod:`delfin.fffree.realism_ranking`.

Covers:

* default OFF: env-flag unset → no behaviour change (byte-identical for any
  pipeline integration; this test asserts the gate function and weight
  invariants only — the module is pure / does not auto-run on import).
* weights sum to 1 (default and after env overrides)
* determinism: two runs over the same input → identical ranking
* signal extraction handles missing / NaN robustly
* rank_xyz_group: best-frame appears at rank 0
* hard gates push failures to the back
* empty group → empty output
* single-frame group → rank 0 with finite score
* all-NaN signals → preserves original order
* archive batch ranking: JSONL sidecar produced, ordered by rank
* archive batch ranking: rewrite-xyz preserves original via sidecar
* env weight override picks up new values
* polya / burnside normalisation: True / 1.0 produces lower (better) score
* ties broken deterministically by frame index
* labels round-trip into the ranking JSONL
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pytest

from delfin.fffree import realism_ranking as RR


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
@pytest.fixture
def env_clean(monkeypatch):
    """Reset all realism-ranking env-flags."""
    for k in list(os.environ):
        if k.startswith("DELFIN_FFFREE_REALISM"):
            monkeypatch.delenv(k, raising=False)
    yield


def _write_multixyz(path: Path, frames):
    """frames = list of (label, atoms)  where atoms is list[str]."""
    with path.open("w") as f:
        for label, atoms in frames:
            f.write(f"{len(atoms)}\n")
            f.write(f"commit=test smi={path.stem} label={label}\n")
            for a in atoms:
                f.write(a + "\n")


# ---------------------------------------------------------------------------
# Gate / config invariants
# ---------------------------------------------------------------------------
def test_default_off(env_clean):
    """Master env-flag default OFF — emission pipeline must not auto-sort."""
    assert RR.realism_sort_active() is False


def test_master_env_on(env_clean, monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_REALISM_SORT", "1")
    assert RR.realism_sort_active() is True


def test_weights_sum_to_one(env_clean):
    w = RR.get_weights()
    assert math.isclose(sum(w.values()), 1.0, abs_tol=1e-9)
    assert all(v >= 0.0 for v in w.values())


def test_weights_env_override(env_clean, monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_REALISM_WEIGHT_MOGUL", "0.8")
    w = RR.get_weights()
    assert math.isclose(sum(w.values()), 1.0, abs_tol=1e-9)
    # Mogul should dominate after the override+normalisation step
    # (raw 0.8 vs other defaults summing to 0.7 → 0.8/1.5 ≈ 0.53).
    assert w["mogul"] > 0.5
    assert w["mogul"] == max(w.values())


def test_weights_all_zero_falls_back(env_clean, monkeypatch):
    for k in RR._DEFAULT_WEIGHTS:
        monkeypatch.setenv(f"DELFIN_FFFREE_REALISM_WEIGHT_{k.upper()}", "0")
    w = RR.get_weights()
    # Fallback to defaults (sum 1.0).
    assert math.isclose(sum(w.values()), 1.0, abs_tol=1e-9)
    assert w == RR._DEFAULT_WEIGHTS


# ---------------------------------------------------------------------------
# Signal extraction
# ---------------------------------------------------------------------------
def test_extract_signals_complete():
    md = {
        "mogul_anomaly_count": 3,
        "cshm": 2.0,
        "inter_ligand_clash_count": 0,
        "hh_clash_count": 0,
        "grip_final_loss": 1.5,
        "polya_complete": True,
        "burnside_coverage": 0.9,
        "topology_ok": True,
        "build_time_clash_ok": True,
    }
    sig = RR.extract_signals(md)
    assert sig["mogul"] == 3.0
    assert sig["cshm"] == 2.0
    assert sig["polya"] == 0.0  # True → inverted to 0 (best)
    # burnside_coverage 0.9 → 0.1 after inversion (close to best)
    assert math.isclose(sig["burnside"], 0.1, abs_tol=1e-9)
    assert sig["topology_ok"] is True


def test_extract_signals_missing_returns_nan():
    sig = RR.extract_signals({})
    for k in ("mogul", "cshm", "inter_clash", "hh_clash", "grip_loss",
              "polya", "burnside"):
        assert math.isnan(sig[k]), f"{k} should be NaN when absent"
    # Hard gates default True (pass).
    assert sig["topology_ok"] is True
    assert sig["build_time_clash_ok"] is True


def test_extract_signals_handles_garbage():
    md = {
        "mogul_anomaly_count": "not-a-number",
        "cshm": float("inf"),
        "polya_complete": "complete",
        "burnside_coverage": 1.5,  # out-of-range → clamped
    }
    sig = RR.extract_signals(md)
    assert math.isnan(sig["mogul"])
    assert math.isnan(sig["cshm"])
    assert sig["polya"] == 0.0  # "complete" string → True → 0
    # 1.5 clamped to 1.0 then inverted → 0.0
    assert sig["burnside"] == 0.0


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------
def test_rank_empty():
    assert RR.rank_xyz_group([]) == []


def test_rank_single_frame_returns_zero():
    out = RR.rank_xyz_group([{"mogul_anomaly_count": 7}])
    assert len(out) == 1
    assert out[0][0] == 0
    assert math.isfinite(out[0][1])


def test_rank_best_frame_at_rank_zero():
    """A frame with strictly-best metrics on every signal must rank-0."""
    frames = [
        {  # frame 0: middling
            "mogul_anomaly_count": 5,
            "cshm": 2.0,
            "inter_ligand_clash_count": 1,
            "hh_clash_count": 1,
            "grip_final_loss": 1.0,
        },
        {  # frame 1: best on every axis
            "mogul_anomaly_count": 0,
            "cshm": 0.1,
            "inter_ligand_clash_count": 0,
            "hh_clash_count": 0,
            "grip_final_loss": 0.0,
        },
        {  # frame 2: worst on every axis
            "mogul_anomaly_count": 10,
            "cshm": 5.0,
            "inter_ligand_clash_count": 4,
            "hh_clash_count": 3,
            "grip_final_loss": 4.0,
        },
    ]
    ranking = RR.rank_xyz_group(frames)
    assert ranking[0][0] == 1  # best frame → rank 0
    assert ranking[-1][0] == 2  # worst frame → last
    # Scores monotonically increase.
    scores = [s for _, s in ranking]
    assert scores == sorted(scores)


def test_hard_gate_failure_pushed_to_back():
    frames = [
        {  # broken topology — should rank last despite great soft signals
            "mogul_anomaly_count": 0,
            "cshm": 0.0,
            "topology_ok": False,
        },
        {  # legitimate good frame
            "mogul_anomaly_count": 2,
            "cshm": 1.0,
        },
        {  # mediocre but passing
            "mogul_anomaly_count": 5,
            "cshm": 3.0,
        },
    ]
    ranking = RR.rank_xyz_group(frames)
    # Hard-gate failure (frame 0) must NOT rank 0.
    assert ranking[0][0] != 0
    # Frame 0 should be last because of the 1e3 penalty.
    assert ranking[-1][0] == 0


def test_all_nan_preserves_insertion_order():
    """When no signals are present every frame gets neutral rank → ties
    broken by original index (stable sort)."""
    frames = [{}, {}, {}, {}]
    ranking = RR.rank_xyz_group(frames)
    assert [orig for orig, _ in ranking] == [0, 1, 2, 3]


def test_ties_broken_deterministically():
    """Identical signals → original index order preserved."""
    frames = [
        {"mogul_anomaly_count": 3, "cshm": 1.0},
        {"mogul_anomaly_count": 3, "cshm": 1.0},
        {"mogul_anomaly_count": 3, "cshm": 1.0},
    ]
    ranking = RR.rank_xyz_group(frames)
    assert [orig for orig, _ in ranking] == [0, 1, 2]


def test_determinism_across_runs():
    frames = [
        {"mogul_anomaly_count": 5, "cshm": 2.0, "hh_clash_count": 1},
        {"mogul_anomaly_count": 0, "cshm": 0.5, "hh_clash_count": 0},
        {"mogul_anomaly_count": 3, "cshm": 1.5, "hh_clash_count": 2},
        {"mogul_anomaly_count": 8, "cshm": 3.0, "hh_clash_count": 4},
    ]
    r1 = RR.rank_xyz_group(frames)
    r2 = RR.rank_xyz_group(frames)
    assert r1 == r2


def test_compute_realism_isolated_mode():
    """Isolated-mode score is finite and monotone in raw signal."""
    md_good = {"mogul_anomaly_count": 0, "cshm": 0.0}
    md_bad = {"mogul_anomaly_count": 50, "cshm": 50.0}
    s_good = RR.compute_realism_score(md_good)
    s_bad = RR.compute_realism_score(md_bad)
    assert math.isfinite(s_good)
    assert math.isfinite(s_bad)
    assert s_good < s_bad


# ---------------------------------------------------------------------------
# Archive batch processing
# ---------------------------------------------------------------------------
def test_rank_archive_sidecar(tmp_path):
    archive = tmp_path / "archive"
    archive.mkdir()
    # Two synthetic XYZ files, three frames each.
    _write_multixyz(
        archive / "smi-A.xyz",
        [
            ("isomer-1", ["C 0 0 0", "H 1 0 0"]),
            ("isomer-2", ["C 0 0 0", "H 1 0 0"]),
        ],
    )
    _write_multixyz(
        archive / "smi-B.xyz",
        [
            ("conf-1", ["N 0 0 0", "H 1 0 0", "H 0 1 0"]),
        ],
    )

    # Detector JSONL with mogul anomaly counts.
    mogul = archive.parent / "archive_mogul_v3.jsonl"
    with mogul.open("w") as f:
        f.write(json.dumps({
            "file": "smi-A.xyz",
            "header": "commit=test smi=smi-A label=isomer-1",
            "n_anomalies": 5,
        }) + "\n")
        f.write(json.dumps({
            "file": "smi-A.xyz",
            "header": "commit=test smi=smi-A label=isomer-2",
            "n_anomalies": 0,
        }) + "\n")
        f.write(json.dumps({
            "file": "smi-B.xyz",
            "header": "commit=test smi=smi-B label=conf-1",
            "n_anomalies": 2,
        }) + "\n")

    out = tmp_path / "ranking.jsonl"
    stats = RR.rank_archive(archive, out, metric_jsonls=[mogul])
    assert stats["n_files"] == 2
    assert stats["n_frames_total"] == 3
    assert out.exists()
    records = [json.loads(l) for l in out.read_text().splitlines() if l.strip()]
    # 2 SMILES → 2 records
    assert len(records) == 2
    by_smi = {r["smiles_id"]: r for r in records}
    # smi-A: isomer-2 (0 anomalies) should rank 0, isomer-1 last.
    assert by_smi["smi-A"]["ranking"][0]["label"] == "isomer-2"
    assert by_smi["smi-A"]["ranking"][-1]["label"] == "isomer-1"


def test_rank_archive_rewrite_xyz_preserves_original(tmp_path):
    archive = tmp_path / "archive"
    archive.mkdir()
    _write_multixyz(
        archive / "test.xyz",
        [
            ("worse", ["C 0 0 0"]),
            ("better", ["C 0 0 0"]),
        ],
    )
    mogul = archive.parent / "archive_mogul_v3.jsonl"
    with mogul.open("w") as f:
        f.write(json.dumps({
            "file": "test.xyz",
            "header": "commit=t smi=test label=worse",
            "n_anomalies": 10,
        }) + "\n")
        f.write(json.dumps({
            "file": "test.xyz",
            "header": "commit=t smi=test label=better",
            "n_anomalies": 0,
        }) + "\n")
    out = tmp_path / "ranking.jsonl"
    RR.rank_archive(archive, out, metric_jsonls=[mogul], rewrite_xyz=True)

    # XYZ rewritten — first frame now the "better" one.
    text = (archive / "test.xyz").read_text()
    first_comment_line = text.splitlines()[1]
    assert "label=better" in first_comment_line

    # Original order preserved.
    orig = json.loads((archive / "test.orig_order.json").read_text())
    assert orig["orig_labels"] == ["worse", "better"]


def test_rank_archive_empty_directory(tmp_path):
    archive = tmp_path / "empty"
    archive.mkdir()
    out = tmp_path / "ranking.jsonl"
    stats = RR.rank_archive(archive, out, metric_jsonls=[])
    assert stats["n_files"] == 0
    assert stats["n_frames_total"] == 0


def test_polya_complete_lowers_score():
    """Setting polya_complete=True must give a STRICTLY lower score than
    polya_complete=False for an otherwise identical frame."""
    f_good = {"polya_complete": True}
    f_bad = {"polya_complete": False}
    ranking = RR.rank_xyz_group([f_bad, f_good])
    # Good frame must rank 0.
    assert ranking[0][0] == 1


def test_burnside_coverage_higher_is_better():
    f_low = {"burnside_coverage": 0.1}
    f_high = {"burnside_coverage": 0.95}
    ranking = RR.rank_xyz_group([f_low, f_high])
    assert ranking[0][0] == 1  # 0.95 coverage → rank 0


# ---------------------------------------------------------------------------
# New env-flags: TUNED weights + SOFT gates (Mission 6 of 2026-06-04 opt)
# ---------------------------------------------------------------------------
def test_tuned_weights_env_flag_default_off(env_clean):
    """``DELFIN_FFFREE_REALISM_TUNED_WEIGHTS`` unset → default weights."""
    assert RR.tuned_weights_active() is False
    w = RR.get_weights()
    # Default profile: mogul=0.30, cshm=0.25, …
    assert abs(w["mogul"] - 0.30) < 1e-9
    assert abs(w["cshm"] - 0.25) < 1e-9


def test_tuned_weights_env_flag_swaps_profile(env_clean, monkeypatch):
    """Setting the flag swaps to the tuned 4-signal-focused profile."""
    monkeypatch.setenv(RR.TUNED_WEIGHTS_ENV, "1")
    assert RR.tuned_weights_active() is True
    w = RR.get_weights()
    # Tuned profile heavily weights hh_clash + inter_clash, drops polya
    # / burnside / grip_loss.
    assert w["hh_clash"] > 0.25
    assert w["inter_clash"] > 0.25
    assert w["polya"] < 1e-6
    assert w["burnside"] < 1e-6


def test_soft_gates_env_flag_default_off(env_clean):
    assert RR.soft_gates_active() is False


def test_soft_gates_env_flag_changes_ranking(env_clean, monkeypatch):
    """With soft-gate mode, a topology-failed frame can still rank ahead
    of a gate-pass frame iff its other soft signals dominate.  Without
    soft-gate mode the 1e3 hard penalty forces it to the back."""
    frames = [
        # Frame A: bad soft signals, gate passes
        {"cshm": 10.0, "mogul_anomaly_count": 10, "topology_ok": True,
         "build_time_clash_ok": True},
        # Frame B: good soft signals, gate FAILS
        {"cshm": 0.1, "mogul_anomaly_count": 0, "topology_ok": False,
         "build_time_clash_ok": True},
    ]
    # Hard gates ON (default): A ranks 0, B is pushed to back.
    rank_hard = RR.rank_xyz_group(frames)
    assert rank_hard[0][0] == 0  # gate-pass frame first
    # Soft gates ON: B's strong soft signals overpower the small gate
    # penalty; B should now rank 0.
    monkeypatch.setenv(RR.SOFT_GATES_ENV, "1")
    rank_soft = RR.rank_xyz_group(frames)
    assert rank_soft[0][0] == 1  # B (better soft) now leads


def test_tuned_and_soft_combination_is_stable(env_clean, monkeypatch):
    """Both flags ON produce a deterministic, finite ranking."""
    monkeypatch.setenv(RR.TUNED_WEIGHTS_ENV, "1")
    monkeypatch.setenv(RR.SOFT_GATES_ENV, "1")
    frames = [
        {"cshm": 1.0, "hh_clash_count": 0, "topology_ok": True,
         "build_time_clash_ok": True},
        {"cshm": 2.0, "hh_clash_count": 5, "topology_ok": True,
         "build_time_clash_ok": True},
    ]
    rank = RR.rank_xyz_group(frames)
    assert len(rank) == 2
    assert all(math.isfinite(s) for _, s in rank)
    # First frame must beat the second under tuned weights (h-clash
    # gets very high weight in the tuned profile).
    assert rank[0][0] == 0
