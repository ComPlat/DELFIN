"""Adversarial tests for HX1 P1/P2 (reviewer nacht-s10).

Each test drives the real public path of delfin/agent/harness_loop.py.
A red test here is a FINDING for the builder; green ones are committed
as regression guards for the reviewer branch.
"""
import math

import pytest

from delfin.agent.harness_loop import (
    Candidate,
    HarnessCandidateError,
    contrast_supported,
    pick_weakest_family,
    validate_candidate,
)


def _candidate(target, old="a", new="b", **kw):
    return Candidate(
        kind=kw.pop("kind", "invariant"),
        target=target,
        old_text=old,
        new_text=new,
        observation_task_id=kw.pop("observation_task_id", "t1"),
        observation=kw.pop("observation", "agent picked the wrong tool"),
        **kw,
    )


# --- family selection ------------------------------------------------------

def test_nan_score_is_refused_not_silently_skipped():
    # A NaN score must not silently fall out of the headroom list and
    # let the picker choose a different family without recording why:
    # NaN means "no measurement", and acting on it is fabrication.
    scores = {"generalist": 60.0, "classic_opt_freq": float("nan")}
    with pytest.raises(ValueError):
        pick_weakest_family(scores)


def test_all_nan_scores_refused():
    with pytest.raises(ValueError):
        pick_weakest_family({"a": float("nan"), "b": float("nan")})


# --- candidate validation: the deny-list is not enough ----------------------

def test_traversal_target_refused():
    # A text-diff candidate must not be able to reach OUTSIDE the repo.
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_candidate("../../etc/hosts.txt"))


def test_absolute_target_refused():
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_candidate("/etc/hosts.txt"))


def test_guidance_file_outside_owned_dirs_is_refused():
    # The loop may edit skill/guidance text only; an .md file in an
    # arbitrary repo location is NOT in its scope (allow-list, not a
    # deny-list).
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_candidate("docs/random_note.md"))


def test_noop_candidate_refused():
    # An empty diff passes "text-diff" shape but has no content to
    # probe; probing it burns budget for a guaranteed no contrast.
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_candidate(
            "delfin/tools/guides/tool_guidance.md", old="same", new="same"))


# --- contrast: real cases ---------------------------------------------------

def test_contrast_gain_exactly_at_band_not_supported():
    # Strictly-above semantics: a gain equal to the band is inside the
    # noise band and must not be reported as supported.
    res = contrast_supported(
        ["t1", "t2"], [80.0, 90.0], [88.0, 98.0], spread=0.0)
    assert res["band"] == pytest.approx(8.0)
    assert res["mean_gain"] == pytest.approx(8.0)
    assert res["supported"] is False


def test_contrast_nan_gain_not_supported():
    # NaN anywhere in the pairing must not produce supported=True; a
    # NaN mean compared with > is False, but the verdict must also not
    # look like a clean "no contrast" — it must be refused.
    with pytest.raises(ValueError):
        contrast_supported(
            ["t1", "t2"], [80.0, float("nan")], [90.0, 90.0], spread=0.0)


def test_contrast_mismatched_lengths_refused():
    with pytest.raises(ValueError):
        contrast_supported(["t1"], [80.0, 90.0], [90.0], spread=0.0)


def test_contrast_single_task_real_band_respected():
    # One task with spread 20 (bench_watch's measured run-to-run noise):
    # band = max(8, 1.5*20) = 30, so a +20 gain must NOT be supported.
    res = contrast_supported(["t1"], [60.0], [80.0], spread=20.0)
    assert res["band"] == pytest.approx(30.0)
    assert res["supported"] is False
