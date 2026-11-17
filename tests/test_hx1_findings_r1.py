"""HX1 fixes for reviewer nacht-s10's findings F1-F5 (round 1).

F1/F2: validate_candidate uses an ALLOW-LIST, not a deny-list:
repo-relative guidance paths under an enumerated allowed root only.
Traversal (..) and absolute paths are refused; a .md outside the
allowed roots is refused.
F3: a no-op candidate (old_text == new_text or empty new_text) is
refused — probing a zero-diff burns budget for a guaranteed
no-contrast.
F4: pick_weakest_family refuses NaN/None scores with ValueError
(silently skipping them manufactured a false LookupError claim).
F5: contrast_supported refuses non-finite reference/probe values with
ValueError (a NaN gain must never read as a clean "no contrast").
"""
from __future__ import annotations

import math

import pytest

from delfin.agent.harness_loop import (
    Candidate,
    HarnessCandidateError,
    RoundConfig,
    contrast_supported,
    pick_weakest_family,
    run_round,
    validate_candidate,
)


def _cand(**overrides):
    fields = dict(
        kind="invariant",
        target="delfin/agent/pack/skills/tool-guidance/SKILL.md",
        old_text="Old guidance text.",
        new_text="New guidance text.",
        observation_task_id="t1",
        observation="something went wrong",
    )
    fields.update(overrides)
    return Candidate(**fields)


# F1: traversal / absolute paths
def test_traversal_target_refused() -> None:
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_cand(target="../../etc/hosts.txt"))


def test_absolute_target_refused() -> None:
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_cand(target="/etc/hosts.txt"))


# F2: allow-list roots
def test_guidance_file_outside_owned_dirs_is_refused() -> None:
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_cand(target="docs/random_note.md"))


def test_allowed_skill_txt_still_passes() -> None:
    # No refusal — the ordinary, correct candidate must still work.
    validate_candidate(_cand())


# F3: no-op
def test_noop_candidate_refused() -> None:
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_cand(new_text="Old guidance text."))


def test_empty_new_text_refused() -> None:
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_cand(new_text=""))


# F4: NaN scores
def test_nan_score_is_refused_not_silently_skipped() -> None:
    with pytest.raises(ValueError):
        pick_weakest_family({"workflow": float("nan"), "esd": 93.0})


def test_all_nan_scores_refused() -> None:
    with pytest.raises(ValueError):
        pick_weakest_family({"a": float("nan"), "b": float("nan")})


def test_none_score_is_refused() -> None:
    with pytest.raises(ValueError):
        pick_weakest_family({"workflow": None, "esd": 93.0})


# F5: NaN in contrast samples
def test_contrast_nan_gain_not_supported() -> None:
    with pytest.raises(ValueError):
        contrast_supported(["t1", "t2"], [50.0, float("nan")],
                           [60.0, 70.0], spread=10.0)


def test_contrast_nan_probe_refused() -> None:
    with pytest.raises(ValueError):
        contrast_supported(["t1"], [50.0], [math.inf], spread=10.0)


# F6 (reviewer round 2): an empty candidate set is refused, not a
# silent noise_band=None "measured round".
def test_empty_candidate_set_refused() -> None:
    base = RoundConfig(
        family_id="workflow", candidate_edits={}, paired_tasks=["t1"],
        reference=[40.0], probe_quality={}, spread=20.0)
    with pytest.raises(ValueError):
        run_round(base, dry_run=True)
