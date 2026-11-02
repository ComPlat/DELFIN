"""Package SW1 phase 3: swarm triage + merge.

Red test first. Control: tests/test_sw1_triage_merge.py, 9 cases, all 9
red on the previous commit (ModuleNotFoundError: delfin.agent.swarm.
triage).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.swarm.triage import (
    TriageError,
    classify,
    triage_answer,
    triage_list,
)
from delfin.agent.swarm.merge import merge_plan


@pytest.fixture()
def swarm(tmp_path, monkeypatch):
    """A swarm dir with a member set, isolated state."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "swarm.yaml").write_text(
        "goal: G\npackages:\n"
        "  - name: core\n"
        "    owned: [\"delfin/pkg/**\"]\n"
        "    task: t\n",
        encoding="utf-8")
    (tmp_path / "delfin" / "pkg").mkdir(parents=True)
    return tmp_path


def test_partner_message_from_same_swarm_is_auto_answerable(swarm):
    verdict = classify({"member": "core-builder", "from": "core-reviewer",
                        "text": "review done, sent APPROVED with marks",
                        "kind": "session_message",
                        "roster": ["core-builder", "core-reviewer"]})
    assert verdict["action"] == "auto_answer"
    assert verdict["reason"].startswith("plain partner message")


def test_push_request_escalates_to_the_user(swarm):
    verdict = classify({"member": "core-builder", "from": "git",
                        "text": "push to origin/main refused -- allow?",
                        "kind": "approval"})
    assert verdict["action"] == "escalate"


def test_permission_request_escalates(swarm):
    verdict = classify({"member": "core-builder", "from": "gate",
                        "text": "rm -rf .swarm/ui-builder not on the "
                                "allow list", "kind": "approval"})
    assert verdict["action"] == "escalate"


def test_own_test_output_is_auto_answerable(swarm):
    verdict = classify({"member": "core-builder", "from": "core-builder",
                        "text": "gate: 12 passed", "kind": "test_output",
                        "roster": ["core-builder", "core-reviewer"]})
    assert verdict["action"] == "auto_answer"


def test_triage_list_is_read_only_and_sorted(swarm):
    queue = [
        {"id": "b", "member": "core-builder", "from": "git",
         "text": "push refused", "kind": "approval"},
        {"id": "a", "member": "core-reviewer", "from": "core-builder",
         "text": "hello", "kind": "session_message",
         "roster": ["core-builder", "core-reviewer"]},
    ]
    listed = triage_list(queue)
    # Escalations first, then auto-answerables; stable within class.
    assert [i["id"] for i in listed] == ["b", "a"]
    assert listed[0]["action"] == "escalate"


def test_triage_answer_logs_and_answers(swarm, capsys, tmp_path,
                                         monkeypatch):
    monkeypatch.setenv("HOME", str(swarm))
    item = {"id": "a", "member": "core-builder", "from": "core-reviewer",
            "text": "APPROVED, marks 5/5/5", "kind": "session_message",
            "roster": ["core-builder", "core-reviewer"]}
    verdict = classify(item)
    assert verdict["action"] == "auto_answer"
    result = triage_answer(item, verdict, answer="ack")
    assert result["answered"] is True
    logged = Path(result["logged"]).read_text(encoding="utf-8")
    assert "auto_answer" in logged and "core-builder" in logged
    # Nothing is written inside the workspace (T4/D14).
    assert not (swarm / ".swarm" / "triage_log.jsonl").exists()


def test_triage_answer_refuses_escalated_items(swarm):
    item = {"id": "b", "member": "core-builder", "from": "git",
            "text": "push refused", "kind": "approval",
            "roster": ["core-builder", "core-reviewer"]}
    verdict = classify(item)
    with pytest.raises(TriageError, match="user"):
        triage_answer(item, verdict, answer="yes push")


def test_merge_plan_prints_the_command_and_runs_nothing(swarm):
    plan = merge_plan("core", base="main")
    assert plan["merge_command"].startswith("git merge")
    assert "core" in plan["merge_command"]
    assert plan["executed"] is False


def test_merge_plan_requires_a_package_name(swarm):
    with pytest.raises(TriageError, match="package"):
        merge_plan("", base="main")
