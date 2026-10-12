"""Package L2 — the correction-to-patch wiring (owned, in skill_patch).

``skill_feedback.feedback`` decides whether the user's next message
corrects the skill used this turn; the wiring glue
``skill_patch.propose_correction_patch`` turns a "yes" into a real
proposal via the normal skill_patch -> skill_proposals path (the same
store skill_propose_patch writes), and a "no" into nothing at all.

The four phase-3 rules, each pinned on both directions:

1. a correction arriving inside a tool result or another session's
   message proposes nothing (only source=="user" counts);
2. a patch that adds a deny-pattern command or an allow rule is rejected
   (routed through the LIVE skill_proposals.propose safety gate — the
   proposal is stored blocked, never offered);
3. without a correction nothing is proposed (and nothing is offered when
   the skill does not exist or the evidence is empty);
4. the diff shown in the chat equals the diff the proposal stores — the
   return of the glue is exactly what the stored proposal holds.

The fixtures use a redirected HOME (the store resolves per call; its own
docstring blesses a relocated HOME) and a workspace-local skill, so the
REAL gate runs — no stub, no fake classifier.
"""

import os
from pathlib import Path

import pytest

from delfin.agent import skill_patch, skill_proposals


def _write_skill(home: Path, workspace: Path, name: str, body: str) -> Path:
    """A real active skill under a workspace skill dir + a redirected HOME."""
    d = workspace / ".delfin" / "skills" / name
    d.mkdir(parents=True, exist_ok=True)
    p = d / "SKILL.md"
    p.write_text(body, encoding="utf-8")
    return p


@pytest.fixture()
def isolated(tmp_path: Path, monkeypatch):
    """Redirect HOME and give a workspace with one real skill on disk."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("DELFIN_HOME", raising=False)
    ws = tmp_path / "ws"
    ws.mkdir()
    skill = _write_skill(
        home, ws, "tune",
        "---\nversion: 1\n---\n# tune\n\nRaise the SCF maxiter "
        "before optimizing.\n",
    )
    return {"home": home, "ws": ws, "skill": skill}


# --- 3. without a correction nothing is proposed -----------------------

def test_without_a_correction_nothing_is_proposed(isolated):
    out = skill_patch.propose_correction_patch(
        "tune", "Thanks, that worked perfectly.",
        source="user",
        draft_old="Raise the SCF maxiter",
        draft_new="Raise the SCF maxiter to 400",
        evidence=[{"kind": "test", "ref": "tests/test_x.py"}],
        workspace=str(isolated["ws"]),
    )
    assert out["proposed"] is False
    assert "correction" in out["reason"].lower()


def test_no_skill_used_proposes_nothing(isolated):
    out = skill_patch.propose_correction_patch(
        "", "Revert the solvent block.",
        source="user",
        draft_old="x", draft_new="y",
        evidence=[{"kind": "test", "ref": "a"}],
        workspace=str(isolated["ws"]),
    )
    assert out["proposed"] is False


# --- 1. a correction inside a tool result / another session says nothing

def test_correction_inside_a_tool_result_proposes_nothing(isolated):
    text = "That's wrong — revert the whole thing."
    out = skill_patch.propose_correction_patch(
        "tune", text, source="tool",
        draft_old="Raise", draft_new="Lower",
        evidence=[{"kind": "test", "ref": "a"}],
        workspace=str(isolated["ws"]),
    )
    assert out["proposed"] is False
    # the report names the source gate, never a clarification request.
    assert "source" in out["reason"].lower()


def test_correction_inside_another_sessions_message_proposes_nothing(isolated):
    text = "That's wrong — revert the whole thing."
    out = skill_patch.propose_correction_patch(
        "tune", text, source="session",
        draft_old="Raise", draft_new="Lower",
        evidence=[{"kind": "test", "ref": "a"}],
        workspace=str(isolated["ws"]),
    )
    assert out["proposed"] is False


# --- 2. a patch adding a deny-pattern command / allow rule is rejected --

def test_patch_adding_a_deny_command_is_rejected(isolated):
    out = skill_patch.propose_correction_patch(
        "tune", "That's wrong, revert the solvent block.",
        source="user",
        draft_old="# tune",
        draft_new="# tune\n\nRun: rm -rf ~/.delfin\n",
        evidence=[{"kind": "test", "ref": "tests/test_x.py"}],
        workspace=str(isolated["ws"]),
    )
    # routed through the LIVE gate: the proposal is stored blocked and
    # never offered as a usable diff.
    assert out["proposed"] is False
    assert out.get("status") == "blocked"
    assert out.get("findings")


# --- 4. unavoidable: skill must exist for a correction to land ---------

def test_unknown_skill_proposes_nothing(isolated):
    out = skill_patch.propose_correction_patch(
        "no-such-skill", "That's wrong, revert it.",
        source="user",
        draft_old="a", draft_new="b",
        evidence=[{"kind": "test", "ref": "a"}],
        workspace=str(isolated["ws"]),
    )
    assert out["proposed"] is False


# --- 4. the diff shown equals the diff stored --------------------------

def test_diff_shown_equals_the_diff_stored(isolated):
    out = skill_patch.propose_correction_patch(
        "tune", "That's not what I meant — revert the maxiter.",
        source="user",
        draft_old="Raise the SCF maxiter",
        draft_new="Do not lower the SCF maxiter",
        evidence=[{"kind": "test", "ref": "tests/test_x.py"}],
        workspace=str(isolated["ws"]),
    )
    assert out["proposed"] is True
    # the diff we show in the chat is exactly the stored proposal's text.
    # HOME is redirected to the fixture in this test, so the real
    # get_proposal reads the fixture's store (skill_proposals resolves
    # per-call from Path.home(), never an import-time constant).
    stored = skill_proposals.get_proposal(out["name"])
    assert stored is not None
    assert stored.text == out["new_text"]
