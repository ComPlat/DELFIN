"""A verification claim in the answer is checked against this turn's runs.

"tests pass", "5 passed", "alle Tests grün" -- and, in a turn that edited
code, "verified" / "funktioniert jetzt" -- need an executed, fresh,
relevant run that was not bought by editing the failing test
(verification_claims.judge over turn_verification's account). An
unbacked claim gets one visible "[verify] Caveat" line, no correction
turn. Driven through AgentEngine.stream_response with a scripted backend
whose account is prepared by hand.

Offline before wiring: 1,063 benchmark answers, 3 flagged, 3 true by
reading; field reports, 249 aligned turns with 102 claim turns, 5
flagged, 2 true (the three false ones: a pass restated from an earlier
turn, a coverage the replay could not see without the files, a phase
summary).
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import channels
from delfin.agent import turn_verification as tv
from delfin.agent import verification_claims as vcl


@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    (lite_dir / "modes").mkdir(parents=True)
    (lite_dir / "modes" / "solo.md").write_text("# solo mode")
    (lite_dir / "manifest.yaml").write_text(
        "pack_name: DELFIN_AGENT_LITE\nversion: 1\nmodes:\n  - id: solo\n"
        "    file: modes/solo.md\n    route:\n      - session_manager\n")
    return tmp_path


def _ok(cmd, code=0, out=""):
    return json.dumps({"exit_code": code, "stdout": out, "stderr": "",
                       "command": cmd})


def _engine(agent_tree, reply, account):
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine
    fake = MagicMock()
    fake._observed_files_session = set()
    fake._turn_verification = account

    def _stream(*a, **k):
        yield StreamEvent(type="text_delta", text=reply)
        yield StreamEvent(type="message_delta", output_tokens=5,
                          cost_usd=0.0)

    fake.stream_message = MagicMock(side_effect=_stream)
    with patch("delfin.agent.engine.create_client", return_value=fake):
        eng = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick",
                          pack_dir=agent_tree)
    return eng, fake


def _account(ws, *steps):
    a = tv.TurnVerification(str(ws))
    for tool, args, res in steps:
        a.note(tool, args, res)
    return a


_REFUSED = '{"error": "bash: \'go test ./...\' is not on the auto-allow list"}'


def test_a_green_claim_after_a_refused_run_gets_the_caveat(agent_tree):
    acc = _account(agent_tree,
                   ("edit_file", {"path": "netcfg/port.go"}, "Edited netcfg/port.go"),
                   ("bash", {"command": "go test ./...", "cwd": "netcfg"}, _REFUSED))
    channels.reset()
    eng, fake = _engine(agent_tree,
                        "Der Fix ist fertig; alle Tests grün.", acc)
    out = eng.stream_response("fix the port check")
    assert "[verify] Caveat: these verification claims are not backed" in out
    assert "no test run in this turn passed" in out
    assert fake.stream_message.call_count == 1          # no correction turn
    assert channels.snapshot().get("test_claim_flagged") == 1
    channels.reset()


def test_a_green_claim_after_a_covering_pass_is_left_alone(agent_tree):
    acc = _account(agent_tree,
                   ("edit_file", {"path": "netcfg/port.go"}, "Edited netcfg/port.go"),
                   ("bash", {"command": "go test ./...", "cwd": "netcfg"},
                    _ok("go test ./...", 0, "ok netcfg")))
    eng, _ = _engine(agent_tree, "Alle Tests grün.", acc)
    out = eng.stream_response("fix the port check")
    assert "[verify]" not in out


def test_a_pass_after_editing_the_failing_test_backs_nothing(agent_tree):
    (agent_tree / "tests").mkdir()
    acc = _account(agent_tree,
                   ("bash", {"command": "pytest -q"}, _ok("pytest -q", 1, "1 failed")),
                   ("edit_file", {"path": "tests/test_t.py", "old_string": "== 2",
                                  "new_string": "== 3"}, "Edited tests/test_t.py"),
                   ("bash", {"command": "pytest -q"}, _ok("pytest -q", 0, "4 passed")))
    eng, _ = _engine(agent_tree, "All tests pass now (4 passed).", acc)
    out = eng.stream_response("make the tests green")
    assert "an edit to the failing test itself" in out


def test_an_honest_answer_gets_no_caveat(agent_tree):
    acc = _account(agent_tree,
                   ("edit_file", {"path": "netcfg/port.go"}, "Edited netcfg/port.go"),
                   ("bash", {"command": "go test ./...", "cwd": "netcfg"}, _REFUSED))
    eng, _ = _engine(agent_tree,
                     "Geändert, aber ungetestet: go test wurde nicht freigegeben.",
                     acc)
    assert "[verify]" not in eng.stream_response("fix the port check")


def test_a_client_without_an_account_changes_nothing(agent_tree):
    eng, _ = _engine(agent_tree, "All tests pass.", None)
    assert eng.stream_response("status?") == "All tests pass."


def test_a_summary_turn_without_tests_or_edits_is_not_judged(agent_tree):
    eng, _ = _engine(agent_tree, "Tests are green (10/10), nothing touches main.",
                     _account(agent_tree))
    assert "[verify]" not in eng.stream_response("where are we?")


def test_the_import_of_the_edited_module_counts_as_touching_it(tmp_path):
    (tmp_path / "pkg").mkdir()
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "test_api.py").write_text(
        "from pkg.engine_core import run\n")
    acc = _account(tmp_path,
                   ("edit_file", {"path": "pkg/engine_core.py"}, "Edited"),
                   ("bash", {"command": "pytest tests/test_api.py -q"},
                    _ok("pytest tests/test_api.py -q", 0, "3 passed")))
    assert acc.unverified_edit() is None
    assert acc.unbacked_claims("3 passed, all green.") == []


def test_a_finished_background_job_is_its_run(tmp_path):
    acc = _account(tmp_path,
                   ("edit_file", {"path": "m.py"}, "Edited m.py"),
                   ("bash_background", {"command": "pytest -q"},
                    json.dumps({"status": "started", "job_id": "j1",
                                "command": "pytest -q"})),
                   ("bash_status", {"job_id": "j1"},
                    json.dumps({"job_id": "j1", "running": False,
                                "exit_code": 0, "command": "pytest -q",
                                "cwd": "."})))
    assert acc.unverified_edit() is None
    assert acc.unbacked_claims("All tests pass.") == []


def test_a_nested_correction_turn_continues_the_account(tmp_path):
    a = _account(tmp_path, ("edit_file", {"path": "m.py"}, "Edited m.py"))
    b = _account(tmp_path, ("bash", {"command": "pytest -q"},
                            _ok("pytest -q", 0, "1 passed")))
    a.extend(b)
    assert a.unverified_edit() is None


def test_the_caveat_names_each_claim_and_its_reason():
    v = vcl.Verdict(vcl.Claim("test", "**Alle 4 Tests grün.**", 4), "stale")
    text = vcl.caveat([v])
    assert "'Alle 4 Tests grün.'" in text
    assert "came before a later edit" in text
    assert vcl.caveat([]) == ""


def test_fixing_a_test_the_turn_wrote_itself_is_not_gaming(tmp_path):
    acc = _account(tmp_path,
                   ("write_file", {"path": "tests/test_words.py"},
                    "File created: tests/test_words.py"),
                   ("bash", {"command": "pytest -q"}, _ok("pytest -q", 1, "1 failed")),
                   ("edit_file", {"path": "tests/test_words.py", "old_string": "a b c …",
                                  "new_string": "a b c…"}, "Edited tests/test_words.py"),
                   ("bash", {"command": "pytest -q"}, _ok("pytest -q", 0, "5 passed")))
    assert tv.gamed_runs(acc.events) == set()
    assert acc.unbacked_claims("5 passed, all green.") == []


def test_a_cut_argument_takes_the_path_from_the_result(tmp_path):
    acc = _account(tmp_path, ("write_file", '{"content": "x', "File created: pkg/new_mod.py"))
    assert acc.events[0]["path"] == "pkg/new_mod.py"
    assert acc.events[0]["created"] is True


@pytest.mark.parametrize("text,kinds", [
    ("Targeted: 32/32 pass.", ["test"]),
    ("Alle vier Fundstellen verifiziert: subagents.py:173, a.py:9.", []),
    ("Zeilenangaben geprüft und verifiziert.", []),
])
def test_claim_forms(text, kinds):
    assert [c.kind for c in vcl.find_claims(text)] == kinds
