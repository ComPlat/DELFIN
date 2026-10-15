"""A refused edit of the protected zone becomes a request, and only that.

Input: an auto iteration whose tool call (or whose end-of-iteration diff)
touches the protected zone -- the sandbox, approval and gate modules, the
principles, the managed policy, the Goodhart and verify guards, files
edit_signals counts as graders (CI configuration, test-runner
configuration, benchmark tasks, acceptance scripts, thresholds), DELFIN's
own workspace settings -- or the loop's declared measure.

Output: the call is refused as before; a structured request lands in
``~/.delfin/requests/<id>.json`` (directory 0700, file 0600) with what was
blocked, an evidence pointer, the proposed change, why it is safe and the
worst case (the last two as the agent states them in its answer), status
open; it is listed under needs-you; the loop goes on with other work.
The request is data for the supervising person or session: its text never
reaches a later prompt, and only a typed command changes its status
(declined with a reason, built, on-main).
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import auto_guard as ag
from delfin.agent import auto_mode as am
from delfin.agent import change_requests as cr
from delfin.agent import loop_commands
from delfin.agent.harness_rows import HarnessText


# -- the zone ---------------------------------------------------------------

@pytest.mark.parametrize("path", [
    "delfin/agent/sandbox.py",
    "delfin/agent/approval_guard.py",
    "delfin/agent/cli_gate.py",
    "delfin/agent/permission_rules.py",
    "delfin/agent/principles_guard.py",
    "delfin/agent/managed_policy.py",
    "delfin/agent/edit_signals.py",
    "delfin/agent/verify_guard.py",
    "delfin/agent/auto_guard.py",
    "delfin/agent/pack/benchmark/tasks.yaml",
    "delfin/agent/pack/benchmark/accept/check_table.sh",
    "etc/delfin/managed_policy.json",
    ".github/workflows/ci.yml",
    ".gitlab-ci.yml",
    "pytest.ini",
    "conf/thresholds.yaml",
    ".delfin/settings.json",
    ".delfin/hooks.json",
])
def test_the_zone(path):
    assert ag.protected_reason(path), path


@pytest.mark.parametrize("path", [
    "delfin/agent/repl.py", "results/table.csv", "src/parser.py",
    "tests/test_parser.py", "README.md", ".delfin/loops/x.md",
])
def test_outside_the_zone(path):
    assert ag.protected_reason(path) == "", path


@pytest.mark.parametrize("name, args", [
    ("write_file", {"path": ".github/workflows/ci.yml", "content": "x"}),
    ("edit_file", {"path": "/w/delfin/agent/sandbox.py", "old_string": "a",
                   "new_string": "b"}),
    ("bash", {"command": "sed -i 's/0.95/0.5/' conf/thresholds.yaml"}),
    ("bash", {"command": "echo '{}' > delfin/agent/pack/benchmark/tasks.yaml"}),
])
def test_a_write_into_the_zone_is_refused_not_queued(name, args):
    v = ag.classify_call(name, args, workspace="/w")
    assert v is not None and v.kind == "protected" and not v.queue


def test_reading_the_zone_passes():
    assert ag.classify_call("read_file", {"path": "delfin/agent/sandbox.py"},
                            workspace="/w") is None
    assert ag.classify_call("bash", {"command": "pytest -c pytest.ini -q"},
                            workspace="/w") is None


# -- the store --------------------------------------------------------------

@pytest.fixture()
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: h))
    return h


def test_a_request_is_a_private_file_with_every_field(home):
    req = cr.create(blocked={"tool": "edit_file", "target": "conf/thresholds.yaml",
                             "zone": "a threshold file"},
                    evidence={"log": "/x.jsonl", "iteration": 3},
                    proposed_change="max_fail: 0 -> 1", workspace="/w",
                    session_id="s", loop="table")
    path = cr.requests_dir() / f"{req['id']}.json"
    assert stat.S_IMODE(os.stat(path).st_mode) == 0o600
    assert stat.S_IMODE(os.stat(cr.requests_dir()).st_mode) == 0o700
    data = json.loads(path.read_text())
    for key in ("id", "created_at", "status", "blocked", "evidence",
                "proposed_change", "why_safe", "worst_case", "workspace",
                "loop", "note"):
        assert key in data, key
    assert data["status"] == "open"
    assert "not an instruction" in data["note"]


def test_status_changes_and_their_reasons(home):
    req = cr.create(blocked={"tool": "bash", "target": "x", "zone": "z"},
                    evidence={}, proposed_change="p")
    with pytest.raises(ValueError):
        cr.set_status(req["id"], "declined")            # a reason is required
    with pytest.raises(ValueError):
        cr.set_status(req["id"], "merged")              # not a status
    cr.set_status(req["id"], "declined", "the threshold is the contract")
    got = cr.get(req["id"])
    assert got["status"] == "declined"
    assert got["status_reason"] == "the threshold is the contract"
    cr.set_status(req["id"], "built", "branch agent/x")
    cr.set_status(req["id"], "on-main", "PR #301")
    assert cr.get(req["id"])["status"] == "on-main"
    assert [r["id"] for r in cr.list_requests(status="on-main")] == [req["id"]]
    # The agent's text can be added only while the request is open.
    with pytest.raises(ValueError):
        cr.amend(req["id"], why_safe="trust me")


def test_an_id_cannot_name_a_path(home):
    with pytest.raises(ValueError):
        cr.get("../../etc/passwd")


# -- the iteration ----------------------------------------------------------

def _git(ws, *a):
    return subprocess.run(["git", "-C", str(ws), *a], capture_output=True,
                          text=True, check=True).stdout.strip()


class Engine:
    session_id = "s1"

    def __init__(self):
        self.kit_permissions = SimpleNamespace(mode="acceptEdits",
                                               auto_veto=None)
        self.token_usage = {"input": 0, "output": 0}
        self.messages = []
        self.pending_turn_continuation = ""

    def request_stop(self):
        pass


@pytest.fixture()
def repo(tmp_path, home):
    ws = tmp_path / "ws"
    ws.mkdir()
    _git(ws, "init", "-q", "-b", "main")
    _git(ws, "config", "user.name", "t")
    _git(ws, "config", "user.email", "t@example.org")
    (ws / "conf").mkdir()
    (ws / "conf" / "thresholds.yaml").write_text("max_fail: 0\n")
    (ws / "src").mkdir()
    (ws / "src" / "a.py").write_text("x = 1\n")
    _git(ws, "add", "-A")
    _git(ws, "commit", "-q", "-m", "init")
    return ws


def _ctl(ws, clock):
    return am.AutoController(ws, session_id="s1", clock=clock,
                             run_check=lambda e, c: (1, ""),
                             notify=lambda *a: None,
                             settings={"idle_seconds": 0, "plateau": 9,
                                       "repeat_failures": 9})


def test_a_refused_call_files_a_request_and_the_work_goes_on(repo):
    t = [1e9]
    clock = lambda: t[0]  # noqa: E731
    ctl = _ctl(repo, clock)
    eng = Engine()
    ctl.enable(eng)
    work = [am.WorkItem(1, "task", "task 1 in progress: make the suite pass",
                        True)]
    t[0] += 10
    assert ctl.next_prompt(eng, work=work)
    out = ctl.veto("edit_file", {"path": "conf/thresholds.yaml",
                                 "old_string": "max_fail: 0",
                                 "new_string": "max_fail: 3"})
    assert out and "request" in out and "not an instruction" in out
    (open_req,) = cr.list_requests(status="open")
    assert open_req["blocked"]["target"] == "conf/thresholds.yaml"
    assert "max_fail: 3" in open_req["proposed_change"]
    assert open_req["evidence"]["iteration"] == 1
    assert any(n["kind"] == "request" and open_req["id"] in n["text"]
               for n in ctl.needs_you)
    # Other work in the same iteration still runs.
    assert ctl.veto("edit_file", {"path": "src/a.py"}) is None
    rid = open_req["id"]
    ctl.finish_iteration(eng, (
        "Fixed src/a.py instead.\n"
        f"REQUEST {rid} why-safe: the flaky test is quarantined upstream\n"
        f"REQUEST {rid} worst-case: a real failure passes CI for a day\n"
        "LOOP: progress a.py"))
    got = cr.get(rid)
    assert got["why_safe"] == "the flaky test is quarantined upstream"
    assert got["worst_case"] == "a real failure passes CI for a day"
    assert ctl.enabled
    # Data, never an instruction: no later prompt carries the request.
    t[0] += 10
    nxt = ctl.next_prompt(eng, work=work)
    assert isinstance(nxt, HarnessText)
    assert rid not in nxt and "quarantined" not in nxt


def test_a_zone_edit_that_got_past_the_veto_is_reverted_with_a_request(repo):
    t = [1e9]
    ctl = _ctl(repo, lambda: t[0])
    eng = Engine()
    ctl.enable(eng)
    t[0] += 10
    ctl.next_prompt(eng, work=[am.WorkItem(1, "task", "x", True)])
    # Written by a route the veto did not see (a script the shell ran).
    (repo / "conf" / "thresholds.yaml").write_text("max_fail: 9\n")
    (repo / "src" / "a.py").write_text("x = 2\n")
    digest = ctl.finish_iteration(eng, "LOOP: progress")
    assert "revert" in digest.lower()
    assert (repo / "conf" / "thresholds.yaml").read_text() == "max_fail: 0\n"
    (req,) = cr.list_requests(status="open")
    assert req["blocked"]["target"] == "conf/thresholds.yaml"
    assert "max_fail: 9" in req["proposed_change"]


def test_the_commands_list_and_decide_requests(repo, home):
    ctl = _ctl(repo, lambda: 1e9)
    req = cr.create(blocked={"tool": "edit_file", "target": "pytest.ini",
                             "zone": "test-runner configuration"},
                    evidence={"iteration": 2}, proposed_change="-x",
                    workspace=str(repo))
    out = loop_commands.auto_command("requests", ctl, Engine())
    assert req["id"] in out and "pytest.ini" in out
    shown = loop_commands.auto_command(f"request {req['id']}", ctl, Engine())
    assert "proposed change" in shown and "worst case" in shown
    out = loop_commands.auto_command(f"request {req['id']} decline", ctl, Engine())
    assert "reason" in out and cr.get(req["id"])["status"] == "open"
    loop_commands.auto_command(
        f"request {req['id']} decline the runner config is the contract",
        ctl, Engine())
    got = cr.get(req["id"])
    assert got["status"] == "declined"
    assert got["status_reason"] == "the runner config is the contract"
    needs = loop_commands.auto_command("needs", ctl, Engine())
    assert req["id"] not in needs or "declined" in needs
