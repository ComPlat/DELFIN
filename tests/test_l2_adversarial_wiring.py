"""Adversarial tests for the L2 wiring glue at the fail-closed seam.

Reviewer nacht-s17 against builder nacht-s13, phase-3 glue landed in the
fail-closed commit 9f861f5b.

skill_feedback.feedback was made fail-closed by F2 (default source
"unknown"); this suite attacks the GLUE boundary the same way. The glue
``skill_patch.propose_correction_patch`` takes its own ``source`` default,
so the question is whether omitting the source at the glue level quietly
re-opens the F2 hole (treating an unknown message as the user's own).

The security rule (waves 14/15, memory-poisoning surface): learned
content comes ONLY from the user's own next message on THIS session. A
caller that does not state a source must therefore propose nothing, just
like at the pure-function level — the default must not silently claim the
user's identity.

The fixtures redirect HOME and give a workspace skill exactly like the
builder's wiring tests, so the real gate and the real proposal store run.
"""
from pathlib import Path

import pytest

from delfin.agent import skill_patch

CORRECTION = "That's wrong — revert the whole thing."


def _write_skill(home: Path, workspace: Path, name: str, body: str) -> Path:
    d = workspace / ".delfin" / "skills" / name
    d.mkdir(parents=True, exist_ok=True)
    p = d / "SKILL.md"
    p.write_text(body, encoding="utf-8")
    return p


@pytest.fixture()
def isolated(tmp_path: Path, monkeypatch):
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


def test_omitted_source_at_the_glue_proposes_nothing(isolated):
    # F2 at the glue boundary: a correction-looking message with NO
    # source must propose nothing — the glue must not default to "user".
    out = skill_patch.propose_correction_patch(
        "tune", CORRECTION,
        draft_old="Raise the SCF maxiter",
        draft_new="Do not lower the SCF maxiter",
        evidence=[{"kind": "test", "ref": "tests/test_x.py"}],
        workspace=str(isolated["ws"]),
    )
    assert out["proposed"] is False
    # it must fail on the SOURCE gate, not on a missing-skill or
    # not-a-correction excuse: the reported reason names the source.
    assert "source" in out["reason"].lower()
