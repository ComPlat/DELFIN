"""Verification as mechanics, not prose-parsing.

Four gates hardened together:
1. structured verdict (report_verdict tool) — beats malformed prose, so a
   rejecting critic can no longer auto-continue on a formatting slip;
2. test-evidence ledger — a PASS claim without a recorded test execution is
   downgraded to UNTESTED, and recorded failures veto a prose PASS;
3. test-tamper gate — editing a test file that was RED earlier in the turn
   records a security event and demands justification;
4. honest auto-verify — a slow suite falls back to a scoped run instead of
   silently turning verification off.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent.engine import AgentEngine


# --- FIX 1a: tolerant prose fallback ---------------------------------------


def test_extract_status_tolerates_trailing_text():
    assert AgentEngine.extract_status_field(
        "**status:** reject — see findings") == "reject"
    assert AgentEngine.extract_status_field(
        "**status:** approve_with_risks (minor nits)") == "approve_with_risks"
    # Exact old format still extracts.
    assert AgentEngine.extract_status_field("**status:** approve") == "approve"


def test_extract_status_accepts_word_variants():
    assert AgentEngine.extract_status_field("**Status:** Rejected") == "reject"
    assert AgentEngine.extract_status_field("**status:** approved") == "approve"
    assert AgentEngine.extract_status_field("no verdict here at all") == ""


def test_rejecting_critic_with_trailing_text_no_longer_auto_continues():
    """The original regex was end-anchored: '**status:** reject — …'
    extracted '' and the gate returned continue — it failed OPEN."""
    out = textwrap.dedent("""\
        ## REVIEW

        ### Critical findings
        1. [engine.py:10] — data loss — **Fix:** revert

        **status:** reject — see findings above
    """)
    action, gate_type, _msg = AgentEngine.evaluate_role_gate(
        "critic_agent", out)
    assert action == "pause"
    assert gate_type == "risk"


# --- FIX 1b: structured verdict beats prose ---------------------------------


def test_structured_reject_beats_missing_prose_status():
    action, gate_type, _ = AgentEngine.evaluate_role_gate(
        "critic_agent", "free-form review with no status line",
        {"status": "reject"})
    assert action == "pause"
    assert gate_type == "risk"


def test_structured_verdict_wins_over_prose():
    # The tool call is the mechanical channel — it overrides the prose line.
    action, _, _ = AgentEngine.evaluate_role_gate(
        "critic_agent", "**status:** reject", {"status": "approve"})
    assert action == "continue"


def test_reviewer_structured_reject_pauses():
    action, gate_type, _ = AgentEngine.evaluate_role_gate(
        "reviewer_agent", "looks fine overall", {"status": "reject"})
    assert action == "pause"
    assert gate_type == "review"


def test_invalid_structured_verdict_falls_back_to_prose():
    action, _, _ = AgentEngine.evaluate_role_gate(
        "critic_agent", "**status:** reject", {"status": "banana"})
    assert action == "pause"


# --- FIX 1c: report_verdict tool registration + executor ---------------------


def test_report_verdict_tool_is_advertised():
    schemas = {t["function"]["name"]: t["function"] for t in A._DOC_TOOLS_OPENAI}
    assert "report_verdict" in schemas
    fn = schemas["report_verdict"]
    assert fn["parameters"]["required"] == ["status"]
    assert set(fn["parameters"]["properties"]["status"]["enum"]) == {
        "approve", "approve_with_risks", "reject"}


def test_report_verdict_executor_echoes_structured_json():
    result = A._doc_executor.execute("report_verdict", {
        "status": "reject",
        "criteria": [{"name": "c1", "state": "fail"},
                     {"name": "c2", "state": "nonsense"}],
        "evidence": "pytest exit 1",
    })
    data = json.loads(result)
    assert data["status"] == "reject"
    assert data["recorded"] is True
    assert data["criteria"][0] == {"name": "c1", "state": "FAIL"}
    # Lenient on sloppy criteria: unknown states become UNTESTED, never lost.
    assert data["criteria"][1]["state"] == "UNTESTED"


def test_report_verdict_executor_rejects_invalid_status():
    result = A._doc_executor.execute("report_verdict", {"status": "maybe"})
    assert result.lstrip().startswith('{"error"')


def test_report_verdict_always_allowed_for_locked_down_roles():
    assert A._tool_denied_for_role("dashboard_agent", "report_verdict") is False


# --- FIX 2: test-evidence ledger ---------------------------------------------


def _gate_engine(test_out="", verdict=None, evidence=None):
    eng = MagicMock()
    eng.role_outputs = {"test_agent": test_out} if test_out else {}
    eng.role_verdicts = {"test_agent": verdict} if verdict else {}
    eng.role_test_evidence = (
        {"test_agent": evidence} if evidence is not None else {})
    return eng


def test_pass_without_evidence_downgraded_to_untested():
    from delfin.dashboard.tab_agent import _check_acceptance_gate
    out = "## TEST REPORT\n1. crit — PASS\n\n**status:** approve\n"
    verdict = _check_acceptance_gate(_gate_engine(test_out=out))
    assert "UNTESTED" in verdict
    assert "no test execution" in verdict


def test_evidence_with_failures_beats_prose_pass():
    from delfin.dashboard.tab_agent import _check_acceptance_gate
    out = "## TEST REPORT\n1. crit — PASS\n\n**status:** approve\n"
    evidence = [{"tool": "run_tests", "command": "", "exit_code": 1,
                 "status": "failed", "passed": 3, "failed": 2, "ts": 0.0}]
    verdict = _check_acceptance_gate(_gate_engine(test_out=out,
                                                  evidence=evidence))
    assert verdict.startswith("❌")
    assert "test evidence" in verdict


def test_pass_with_green_evidence_stays_green():
    from delfin.dashboard.tab_agent import _check_acceptance_gate
    out = "## TEST REPORT\n1. crit — PASS\n\n**status:** approve\n"
    evidence = [{"tool": "run_tests", "command": "", "exit_code": 0,
                 "status": "ok", "passed": 5, "failed": 0, "ts": 0.0}]
    verdict = _check_acceptance_gate(_gate_engine(test_out=out,
                                                  evidence=evidence))
    assert verdict.startswith("✅")


def test_skipped_auto_verify_is_not_execution_evidence():
    from delfin.dashboard.tab_agent import _check_acceptance_gate
    out = "**status:** approve\n1. crit — PASS\n"
    evidence = [{"tool": "auto_verify", "command": "", "exit_code": None,
                 "status": "skipped", "passed": 0, "failed": 0, "ts": 0.0}]
    verdict = _check_acceptance_gate(_gate_engine(test_out=out,
                                                  evidence=evidence))
    assert "UNTESTED" in verdict


def test_structured_verdict_beats_prose_token_counting():
    from delfin.dashboard.tab_agent import _check_acceptance_gate
    # Prose full of PASS tokens, but the STRUCTURED verdict is reject.
    out = "PASS PASS PASS **status:** approve\n"
    verdict_dict = {"status": "reject",
                    "criteria": [{"name": "c1", "state": "FAIL"}]}
    verdict = _check_acceptance_gate(
        _gate_engine(test_out=out, verdict=verdict_dict))
    assert verdict.startswith("❌")
    assert "1 FAIL" in verdict


def test_observe_run_tests_records_evidence_and_red_files():
    evidence: list = []
    red: set = set()
    result = json.dumps({
        "status": "failed", "exit_code": 1,
        "summary": {"passed": 2, "failed": 1, "errors": 0},
        "failures": [{"node_id": "tests/test_x.py::test_y",
                      "message": "boom"}],
    })
    note = A._observe_test_evidence(evidence, red, "run_tests",
                                    {"target": "tests/"}, result)
    assert note == ""
    assert evidence[0]["failed"] == 1 and evidence[0]["exit_code"] == 1
    assert "tests/test_x.py" in red


def test_observe_bash_pytest_records_evidence():
    evidence: list = []
    red: set = set()
    result = json.dumps({
        "exit_code": 1, "stdout":
            "FAILED tests/test_mod.py::test_val - AssertionError\n"
            "1 failed, 3 passed in 0.2s\n",
        "stderr": "", "command": "python -m pytest -q",
    })
    A._observe_test_evidence(evidence, red, "bash",
                             {"command": "python -m pytest -q"}, result)
    assert evidence[0]["passed"] == 3 and evidence[0]["failed"] == 1
    assert "tests/test_mod.py" in red


def test_observe_non_test_bash_records_nothing():
    evidence: list = []
    red: set = set()
    A._observe_test_evidence(evidence, red, "bash",
                             {"command": "ls -la"}, '{"exit_code": 0}')
    assert evidence == [] and red == set()


def test_echo_pytest_cannot_clear_red_state():
    """`echo pytest` (exit 0) must neither count as a test run nor clear
    the tamper gate's red set — the runner must be the COMMAND."""
    evidence: list = []
    red = {"tests/test_mod.py"}
    A._observe_test_evidence(
        evidence, red, "bash", {"command": "echo pytest"},
        json.dumps({"exit_code": 0, "stdout": "pytest", "stderr": ""}))
    assert evidence == []
    assert red == {"tests/test_mod.py"}


def test_pytest_version_exit0_does_not_clear_red_state():
    evidence: list = []
    red = {"tests/test_mod.py"}
    A._observe_test_evidence(
        evidence, red, "bash", {"command": "pytest --version"},
        json.dumps({"exit_code": 0, "stdout": "pytest 8.0.0", "stderr": ""}))
    # Recorded (it IS a pytest invocation) but ran nothing → red survives.
    assert red == {"tests/test_mod.py"}


def test_pytest_version_is_not_execution_evidence_for_the_gate():
    from delfin.dashboard.tab_agent import _check_acceptance_gate
    out = "**status:** approve\n1. crit — PASS\n"
    evidence = [{"tool": "bash", "command": "pytest --version",
                 "exit_code": 0, "status": "ok",
                 "passed": 0, "failed": 0, "ts": 0.0}]
    verdict = _check_acceptance_gate(_gate_engine(test_out=out,
                                                  evidence=evidence))
    assert "UNTESTED" in verdict


def test_engine_copies_verdict_and_evidence_per_role(tmp_path):
    """The engine snapshots the client's per-turn mechanics per role so the
    acceptance gate can consult the TEST role's turn at cycle end."""
    agent_dir = tmp_path / "pack"
    (agent_dir / "shared").mkdir(parents=True)
    (agent_dir / "agents").mkdir()
    for name in ("delfin_context.md", "work_cycle_rules.md",
                 "goal_decomposition_rules.md", "universal_input_template.md",
                 "minimal_final_verdict.md"):
        (agent_dir / "shared" / name).write_text("# x")
    for name in ("session_manager.md", "builder_agent.md", "test_agent.md"):
        (agent_dir / "agents" / name).write_text("# role")
    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "quick.md").write_text("# quick mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: quick
            file: modes/quick.md
            route:
              - session_manager
              - builder_agent
              - test_agent
    """))

    from delfin.agent.api_client import StreamEvent

    def fake_stream(system, messages, max_tokens=4096, session_id="",
                    thinking_budget=0):
        yield StreamEvent(type="text_delta", text="done")
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=fake_stream)
    client._last_structured_verdict = {"status": "reject"}
    client._test_evidence = [{"tool": "run_tests", "failed": 1,
                              "status": "failed", "ts": 0.0}]

    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick",
                             pack_dir=tmp_path)
    engine.stream_response("run the tests")

    role = engine.current_role
    assert engine._last_structured_verdict == {"status": "reject"}
    assert engine.role_verdicts[role] == {"status": "reject"}
    assert engine.role_test_evidence[role][0]["failed"] == 1


# --- FIX 3: test-tamper gate -------------------------------------------------


def _seed_red(evidence, red):
    result = json.dumps({
        "status": "failed", "exit_code": 1,
        "summary": {"passed": 0, "failed": 1, "errors": 0},
        "failures": [{"node_id": "tests/test_mod.py::test_val"}],
    })
    A._observe_test_evidence(evidence, red, "run_tests", {"target": ""}, result)


def test_tamper_gate_records_event_and_demands_justification(monkeypatch):
    events: list = []
    monkeypatch.setattr(
        A, "_record_security_event",
        lambda kind, tool, detail, blocked=True: events.append(
            (kind, tool, detail, blocked)))
    evidence: list = []
    red: set = set()
    _seed_red(evidence, red)

    note = A._observe_test_evidence(
        evidence, red, "edit_file",
        {"path": "tests/test_mod.py", "old": "== 2", "new": "== 3"},
        '{"status": "ok"}')
    assert "TEST-TAMPER GATE" in note
    assert "why the test (not the code) was wrong" in note.lower() \
        or "why the test (not the" in note
    assert events and events[0][0] == "test_tamper"
    assert events[0][3] is False        # recorded, not blocked


def test_tamper_gate_covers_any_test_file_while_red(monkeypatch):
    monkeypatch.setattr(A, "_record_security_event",
                        lambda *a, **k: None)
    evidence: list = []
    red: set = set()
    _seed_red(evidence, red)
    # A DIFFERENT test file than the failing one still needs justification.
    note = A._observe_test_evidence(
        evidence, red, "write_file",
        {"path": "tests/test_other.py", "content": "def test(): pass"},
        '{"status": "ok"}')
    assert "TEST-TAMPER GATE" in note


def test_tamper_gate_ignores_source_edits(monkeypatch):
    monkeypatch.setattr(A, "_record_security_event",
                        lambda *a, **k: None)
    evidence: list = []
    red: set = set()
    _seed_red(evidence, red)
    note = A._observe_test_evidence(
        evidence, red, "edit_file", {"path": "pkg/mod.py"},
        '{"status": "ok"}')
    assert note == ""


def test_tamper_gate_catches_apply_patch_diffs(monkeypatch):
    monkeypatch.setattr(A, "_record_security_event",
                        lambda *a, **k: None)
    evidence: list = []
    red: set = set()
    _seed_red(evidence, red)
    diff = ("--- a/tests/test_mod.py\n+++ b/tests/test_mod.py\n"
            "@@ -1 +1 @@\n-assert val() == 2\n+assert val() == 3\n")
    note = A._observe_test_evidence(
        evidence, red, "apply_patch", {"diff": diff}, '{"status": "ok"}')
    assert "TEST-TAMPER GATE" in note


def test_green_full_run_clears_red_state(monkeypatch):
    monkeypatch.setattr(A, "_record_security_event",
                        lambda *a, **k: None)
    evidence: list = []
    red: set = set()
    _seed_red(evidence, red)
    assert red
    green = json.dumps({"exit_code": 0,
                        "stdout": "4 passed in 0.1s", "stderr": ""})
    A._observe_test_evidence(evidence, red, "bash",
                             {"command": "python -m pytest -q"}, green)
    assert red == set()
    # After green, a test edit no longer needs justification.
    note = A._observe_test_evidence(
        evidence, red, "edit_file", {"path": "tests/test_mod.py"},
        '{"status": "ok"}')
    assert note == ""


def test_green_scoped_run_only_clears_its_own_file(monkeypatch):
    monkeypatch.setattr(A, "_record_security_event",
                        lambda *a, **k: None)
    evidence: list = []
    red = {"tests/test_a.py", "tests/test_b.py"}
    green = json.dumps({"exit_code": 0, "stdout": "1 passed", "stderr": ""})
    A._observe_test_evidence(
        evidence, red, "bash",
        {"command": "python -m pytest tests/test_a.py -q"}, green)
    assert red == {"tests/test_b.py"}


# --- FIX 4: honest auto-verify (see also tests/test_auto_verify.py) ---------


def test_paths_from_diff_extracts_edited_files():
    diff = ("--- a/pkg/mod.py\n+++ b/pkg/mod.py\n@@ -1 +1 @@\n-x\n+y\n"
            "--- /dev/null\n+++ b/pkg/new.py\n@@ -0,0 +1 @@\n+z\n")
    assert A._paths_from_diff(diff) == ["pkg/mod.py", "pkg/new.py"]


def test_apply_patch_and_skip_notices_are_wired():
    """Source-level wiring assertions (same style as test_loop_is_wired):
    apply_patch feeds the verify gate, and skipped/exhausted verification is
    surfaced visibly instead of silently reading as clean."""
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    assert 'fn_name == "apply_patch"' in src
    assert "_paths_from_diff(" in src
    assert "Auto-verify skipped:" in src
    assert "NOT machine-verified" in src
    assert "Auto-verify gave up after" in src
    assert "_observe_test_evidence(" in src


def test_verdict_capture_is_wired_into_the_tool_loop():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    assert "_last_structured_verdict = _sv" in src or \
        "self._last_structured_verdict = _sv" in src
    assert "self._test_evidence" in src
    assert "self._red_test_files" in src


def test_gate_role_prompts_mention_report_verdict():
    root = Path(__file__).resolve().parent.parent
    for name in ("critic_agent", "test_agent", "reviewer_agent"):
        body = (root / "delfin" / "agent" / "pack" / "agents"
                / f"{name}.md").read_text(encoding="utf-8")
        assert "report_verdict" in body, f"{name}.md lacks report_verdict note"
