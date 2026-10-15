"""run_tests runs the test verb of the target's stack and says when nothing ran.

Diagnosis 2026-10-15: run_tests was pytest-only ("run_tests ist
pytest-basiert, passt nicht zu Go" in 4 traces), and it answered
``status: ok`` for 0 passed / 5 skipped. Now:

* stack_tests.plan picks cargo / go / npm / node / pytest per target from
  the target and the marker files above it;
* a non-pytest group runs as one bash call in its project, through the
  bash gate (unasked in the cage, asked or refused outside it);
* the summary counts passed, failed, errors, skipped and xfailed from the
  runner's own output, and 0 passed with nothing failed is
  ``nothing_ran``, for pytest too;
* several targets in different stacks run as several groups, summed.
"""

from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import stack_tests as st
from delfin.agent import test_runner as TR
from delfin.agent import turn_verification as tv
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    d = tmp_path / "ws"
    (d / "rs" / "src").mkdir(parents=True)
    (d / "rs" / "Cargo.toml").write_text('[package]\nname = "rs"\n')
    (d / "rs" / "src" / "lib.rs").write_text("")
    (d / "rs" / "tests").mkdir()
    (d / "rs" / "tests" / "conv.rs").write_text("")
    (d / "gomod" / "netcfg").mkdir(parents=True)
    (d / "gomod" / "go.mod").write_text("module gomod\n")
    (d / "gomod" / "netcfg" / "port.go").write_text("package netcfg\n")
    (d / "ts" / "test").mkdir(parents=True)
    (d / "ts" / "package.json").write_text(json.dumps(
        {"scripts": {"test": "node --test test/"}}))
    (d / "js" / "test").mkdir(parents=True)
    (d / "js" / "package.json").write_text(json.dumps({"name": "js"}))
    (d / "js" / "test" / "a.test.js").write_text("")
    (d / "py" / "tests").mkdir(parents=True)
    (d / "py" / "pyproject.toml").write_text("[project]\nname='p'\n")
    (d / "py" / "tests" / "test_a.py").write_text("def test_a():\n    pass\n")
    return d


@pytest.mark.parametrize("target,framework,project,command", [
    ("rs", "cargo", "rs", "cargo test --offline"),
    ("rs/src/lib.rs", "cargo", "rs", "cargo test --offline"),
    ("rs/tests/conv.rs", "cargo", "rs", "cargo test --offline --test conv"),
    ("gomod", "go", "gomod", "go test -v ./..."),
    ("gomod/netcfg/port.go", "go", "gomod", "go test -v ./netcfg/..."),
    ("ts", "npm", "ts", "npm test"),
    ("js/test/a.test.js", "node", "js", "node --test test/a.test.js"),
    ("py/tests/test_a.py", "pytest", "", ""),
])
def test_the_plan_names_the_stack_of_the_target(ws, target, framework,
                                                project, command):
    groups = st.plan(ws, target)
    assert [g.framework for g in groups] == [framework]
    if project:
        assert groups[0].project.resolve() == (ws / project).resolve()
    assert groups[0].command == command


def test_no_target_and_no_marker_is_pytest_as_before(tmp_path):
    assert [g.framework for g in st.plan(tmp_path, "")] == ["pytest"]


def test_several_targets_in_two_stacks_are_two_groups(ws):
    groups = st.plan(ws, "rs gomod")
    assert sorted(g.framework for g in groups) == ["cargo", "go"]


@pytest.mark.parametrize("framework,code,output,want", [
    ("cargo", 0, "test result: ok. 3 passed; 0 failed; 1 ignored; 0 measured",
     {"passed": 3, "failed": 0, "skipped": 1, "errors": 0}),
    ("cargo", 101, "test result: FAILED. 2 passed; 1 failed; 0 ignored;",
     {"passed": 2, "failed": 1}),
    ("cargo", 101, "error[E0425]: cannot find value `x`", {"errors": 1}),
    ("go", 1, "--- PASS: TestA (0.00s)\n--- FAIL: TestB (0.00s)\nFAIL",
     {"passed": 1, "failed": 1}),
    ("npm", 0, "ℹ tests 4\nℹ pass 3\nℹ fail 0\nℹ skipped 1\nℹ todo 0",
     {"passed": 3, "skipped": 1}),
    ("npm", 1, " Tests  1 failed | 2 passed (3)", {"passed": 2, "failed": 1}),
    ("npm", 1, "Tests:       1 failed, 4 passed, 5 total",
     {"passed": 4, "failed": 1}),
    ("npm", 0, "  6 passing (10ms)\n  1 pending", {"passed": 6, "skipped": 1}),
])
def test_each_runner_output_is_counted(framework, code, output, want):
    got = st.summarise(framework, code, output)
    for k, v in want.items():
        assert got[k] == v, (k, got)


def test_zero_passed_is_nothing_ran_not_ok():
    assert st.status_of({"passed": 0, "failed": 0, "errors": 0,
                         "skipped": 5}, 0) == "nothing_ran"
    assert st.status_of({"passed": 3, "failed": 0, "errors": 0}, 0) == "ok"
    assert st.status_of({"passed": 3, "failed": 1, "errors": 0}, 1) == "failed"


def test_pytest_with_every_test_skipped_is_nothing_ran():
    repo = Path(tempfile.mkdtemp(prefix="rt_skip_"))
    try:
        (repo / "test_s.py").write_text(
            "import pytest\n\n"
            "@pytest.mark.skip(reason='no gpu')\ndef test_a():\n    pass\n\n"
            "@pytest.mark.skip(reason='no gpu')\ndef test_b():\n    pass\n")
        r = TR.run_tests(repo, target="test_s.py")
        assert r["status"] == "nothing_ran", r
        assert r["summary"]["skipped"] == 2 and r["summary"]["passed"] == 0
        assert "no test passed" in r["note"]
    finally:
        shutil.rmtree(repo, ignore_errors=True)


def test_an_expected_failure_is_counted_beside_the_passes():
    repo = Path(tempfile.mkdtemp(prefix="rt_xf_"))
    try:
        (repo / "test_x.py").write_text(
            "import pytest\n\ndef test_ok():\n    pass\n\n"
            "@pytest.mark.xfail(reason='known')\ndef test_known():\n"
            "    assert False\n")
        r = TR.run_tests(repo, target="test_x.py")
        assert r["status"] == "ok", r
        assert r["summary"]["passed"] == 1
        assert r["summary"]["xfailed"] == 1
    finally:
        shutil.rmtree(repo, ignore_errors=True)


def _bash_double(monkeypatch, answer):
    """The bash call run_tests makes, answered by ``answer(args)``; every
    other tool runs for real."""
    seen: list[dict] = []
    real = A._DocToolExecutor.execute

    def fake(self, name, arguments, permissions=None):
        if name == "bash":
            seen.append(dict(arguments))
            return json.dumps(answer(arguments))
        return real(self, name, arguments, permissions)

    monkeypatch.setattr(A._DocToolExecutor, "execute", fake)
    return seen


def _run(ws, target, **kw):
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    return json.loads(_doc_executor.execute(
        "run_tests", {"target": target, **kw}, perms))


def test_a_cargo_target_runs_cargo_in_its_crate(ws, monkeypatch):
    seen = _bash_double(monkeypatch, lambda a: {
        "exit_code": 0, "elapsed_s": 0.3, "command": a["command"],
        "stdout": "test result: ok. 4 passed; 0 failed; 0 ignored;",
        "stderr": ""})
    r = _run(ws, "rs")
    assert seen[0]["command"] == "cargo test --offline"
    assert Path(seen[0]["cwd"]).resolve() == (ws / "rs").resolve()
    assert r["framework"] == "cargo" and r["status"] == "ok"
    assert r["summary"]["passed"] == 4 and r["ran_in"] == "rs"


def test_two_stacks_are_summed_and_the_worst_status_wins(ws, monkeypatch):
    def answer(a):
        if a["command"].startswith("cargo"):
            return {"exit_code": 0, "command": a["command"], "stderr": "",
                    "stdout": "test result: ok. 2 passed; 0 failed; 0 ignored;"}
        return {"exit_code": 1, "command": a["command"], "stderr": "",
                "stdout": "--- PASS: TestA (0.00s)\n--- FAIL: TestB (0.00s)"}
    _bash_double(monkeypatch, answer)
    r = _run(ws, "rs gomod")
    assert r["framework"] == "multi" and r["status"] == "failed"
    assert r["summary"]["passed"] == 3 and r["summary"]["failed"] == 1
    assert len(r["runs"]) == 2


def _cage(monkeypatch, kind: str) -> None:
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: (kind, "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: kind == "bwrap")
    monkeypatch.setattr(A, "_landlock_functional", lambda: kind == "bwrap")


def test_outside_the_cage_a_stack_run_is_refused_not_reported_green(
        ws, monkeypatch):
    _cage(monkeypatch, "none")
    r = _run(ws, "rs")
    assert r["status"] == "refused", r
    assert "UNTESTED" in r["error"]
    assert "summary" not in r


def test_the_account_scopes_a_stack_run_to_its_project():
    acc = tv.TurnVerification("")
    acc.note("edit_file", {"path": "rs/src/lib.rs"}, "Edited rs/src/lib.rs")
    acc.note("run_tests", {"target": "gomod"}, json.dumps({
        "status": "ok", "framework": "go", "ran_in": "gomod",
        "exit_code": 0, "summary": {"passed": 2, "failed": 0, "errors": 0}}))
    assert acc.unverified_edit() == "rs/src/lib.rs"
    acc.note("run_tests", {"target": "rs"}, json.dumps({
        "status": "ok", "framework": "cargo", "ran_in": "rs",
        "exit_code": 0, "summary": {"passed": 4, "failed": 0, "errors": 0}}))
    assert acc.unverified_edit() is None


def test_a_nothing_ran_result_is_not_a_check_of_the_edit():
    acc = tv.TurnVerification("")
    acc.note("edit_file", {"path": "py/m.py"}, "Edited py/m.py")
    acc.note("run_tests", {"target": "py"}, json.dumps({
        "status": "nothing_ran", "framework": "pytest", "exit_code": 0,
        "summary": {"passed": 0, "failed": 0, "errors": 0, "skipped": 5}}))
    assert acc.unverified_edit() == "py/m.py"
