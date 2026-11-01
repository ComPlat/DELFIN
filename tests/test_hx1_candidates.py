"""HX1 Phase 2 red test: candidate edits are TEXT diffs to allowed
targets only, and each cites an observation.

Contract under test: delfin/agent.harness_loop.Candidate / validate_candidate
- a candidate is a text diff (old_text/new_text) plus a named observation
  (task id + what went wrong) and a kind (invariant|domain);
- a candidate whose target path ends in a code extension (.py, .sh,
  .ipynb) or lands under a protected path (delfin/agent/pack/,
  delfin/tools/, any */sandbox*, api_client, engine.py, …) is REFUSED
  (HarnessCandidateError) before probing — the brief's first red test;
- a candidate without an observation (empty task id or empty what-went-
  wrong) is refused — "the report cites an observation for every
  decision";
- an unknown kind is refused (one invariant, one domain per brief).
"""
from __future__ import annotations

import pytest

from delfin.agent.harness_loop import (
    Candidate,
    HarnessCandidateError,
    validate_candidate,
)


def _ok_candidate(**overrides):
    fields = dict(
        kind="invariant",
        target="delfin/agent/pack/skills/tool-guidance/SKILL.md",
        old_text="Before writing an input file, ask the user.",
        new_text="Before writing an input file, read CONTROL.txt first.",
        observation_task_id="hx1_ref_003",
        observation="Agent wrote an ORCA input without reading CONTROL.txt "
                    "and guessed method=basis from memory; task failed "
                    "grading on the method key.",
    )
    fields.update(overrides)
    return Candidate(**fields)


def test_valid_candidate_passes() -> None:
    validate_candidate(_ok_candidate())


def test_code_target_is_refused() -> None:
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_ok_candidate(target="delfin/agent/tools.py"))


def test_protected_path_is_refused() -> None:
    with pytest.raises(HarnessCandidateError):
        validate_candidate(
            _ok_candidate(target="delfin/agent/api_client.py"))
    with pytest.raises(HarnessCandidateError):
        validate_candidate(
            _ok_candidate(target="delfin/agent/pack/knowledge/SKILL.md"))


def test_missing_observation_is_refused() -> None:
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_ok_candidate(observation_task_id=""))
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_ok_candidate(observation="  "))


def test_unknown_kind_is_refused() -> None:
    with pytest.raises(HarnessCandidateError):
        validate_candidate(_ok_candidate(kind="security"))


def test_domain_kind_is_valid() -> None:
    validate_candidate(_ok_candidate(kind="domain"))
