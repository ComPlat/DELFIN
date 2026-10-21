"""The verification signals use the measure-edit detection when it exists.

Interfaces with the experiment-discipline modules (edit_signals,
paired_compare; on main since #269), both used through adapters that
degrade to this bundle's own checks when the module cannot be imported.
Both paths are tested: the real modules, and an import that fails.

* edit_signals.measure_edits(session, after_red_only=True): a test file it
  names makes every later pass gamed (turn_verification.account), on top
  of the account's own rule;
* edit_signals.turn_caveat: when its caveat already names the changed
  check, the test-claim caveat does not repeat it;
* paired_compare.provenance: the benchmark run header uses its stamp as
  the base, and fills the same fields itself without it.

counter_test (on this base, pytest-only, session-level) and the account's
gamed_runs (per turn, any stack) state one rule -- a test file edited
between its red and its green run is no control -- and the pin below
makes them agree on a shared case.
"""

from __future__ import annotations

import json
import sys
import types

from delfin.agent import benchmark as bm
from delfin.agent import counter_test as ct
from delfin.agent import turn_verification as tv


def _ok(cmd, code=0, out=""):
    return json.dumps({"exit_code": code, "stdout": out, "stderr": "",
                       "command": cmd})


def _account(*steps, session=""):
    a = tv.TurnVerification("", session=session)
    for tool, args, res in steps:
        a.note(tool, args, res)
    return a


# An additive edit to the red test file: the account's own rule lets it
# pass (a regression test added beside the red one).
_STEPS = (
    ("bash", {"command": "pytest -q"}, _ok("pytest -q", 1, "1 failed")),
    ("edit_file", {"path": "tests/test_t.py", "old_string": "def test_a():",
                   "new_string": "def test_a():\n    x = 1"},
     "Edited tests/test_t.py"),
    ("bash", {"command": "pytest -q"}, _ok("pytest -q", 0, "3 passed")),
)


def _fake_edit_signals(monkeypatch, paths):
    mod = types.ModuleType("delfin.agent.edit_signals")
    seen = {}

    def measure_edits(session, *, after_red_only=False):
        seen["call"] = (session, after_red_only)
        return [types.SimpleNamespace(path=p, kind="assertion widened")
                for p in paths]

    mod.measure_edits = measure_edits
    monkeypatch.setitem(sys.modules, "delfin.agent.edit_signals", mod)
    import delfin.agent as pkg
    monkeypatch.setattr(pkg, "edit_signals", mod, raising=False)
    return seen


def _without(monkeypatch, name):
    """The import of delfin.agent.<name> fails, as on a tree without it."""
    import delfin.agent as pkg
    monkeypatch.setitem(sys.modules, f"delfin.agent.{name}", None)
    monkeypatch.delattr(pkg, name, raising=False)


def test_without_the_detection_the_own_rule_decides(monkeypatch):
    _without(monkeypatch, "edit_signals")
    assert tv.external_measure_paths("s1") is None
    acc = _account(*_STEPS, session="s1")
    assert acc.unbacked_claims("All tests pass (3 passed).") == []


def test_the_real_detection_names_a_widened_assertion(tmp_path):
    from delfin.agent import edit_signals as es
    session = "verify-iface-real"
    es.reset_turn(session)
    path = str(tmp_path / "tests" / "test_t.py")
    es.note_command(session, "pytest -q", "1 failed, 2 passed")
    found = es.note_change(session, path,
                           "def test_a():\n    assert f() == 2\n",
                           "def test_a():\n    assert f() in (2, 3)\n")
    assert es.measure_edits(session, after_red_only=True), found
    assert tv.external_measure_paths(session) == {"test_t.py"}
    acc = _account(*_STEPS, session=session)
    assert [v.reason for v in acc.unbacked_claims(
        "All tests pass (3 passed).")] == ["gamed"]
    es.reset_turn(session)


def test_a_file_the_detection_names_makes_the_later_pass_gamed(monkeypatch):
    seen = _fake_edit_signals(monkeypatch, ["/ws/tests/test_t.py"])
    acc = _account(*_STEPS, session="s1")
    verdicts = acc.unbacked_claims("All tests pass (3 passed).")
    assert [v.reason for v in verdicts] == ["gamed"]
    assert seen["call"] == ("s1", True)


def test_a_detection_that_names_another_file_changes_nothing(monkeypatch):
    _fake_edit_signals(monkeypatch, ["/ws/tests/test_other.py"])
    acc = _account(*_STEPS, session="s1")
    assert acc.unbacked_claims("All tests pass (3 passed).") == []


def test_counter_test_and_the_account_agree_on_an_edited_red_test():
    # A red run, the red test file rewritten (not additively), a green run.
    steps = (
        ("bash", {"command": "pytest -q tests/test_t.py"},
         _ok("pytest -q tests/test_t.py", 1,
             "FAILED tests/test_t.py::test_a\n1 failed")),
        ("edit_file", {"path": "tests/test_t.py", "old_string": "== 2",
                       "new_string": "== 3"}, "Edited tests/test_t.py"),
        ("bash", {"command": "pytest -q tests/test_t.py"},
         _ok("pytest -q tests/test_t.py", 0, "1 passed")),
    )
    acc = _account(*steps)
    assert tv.gamed_runs(acc.events) == {2.0}
    runs = ct.normalize_test_runs([
        {"command": "pytest -q tests/test_t.py",
         "output": "FAILED tests/test_t.py::test_a\n1 failed",
         "status": "failed", "failed": 1, "ts": 1.0},
        {"command": "pytest -q tests/test_t.py", "output": "1 passed",
         "status": "ok", "failed": 0, "ts": 3.0}])
    # The fix turn changed the product too; assess judges only then.
    verdict = ct.assess([{"path": "pkg/m.py", "ts": 2.0},
                         {"path": "tests/test_t.py", "ts": 2.0}], runs)
    assert verdict and verdict.startswith("no control run")


def test_the_run_header_has_the_generic_stamp_fields_without_the_module(
        monkeypatch):
    _without(monkeypatch, "paired_compare")
    head = bm.run_provenance(model="m")
    assert head.get("stamp_schema") is None
    for key in ("taken_at", "host", "platform", "machine", "python", "code",
                "instrument", "settings", "process_start"):
        assert key in head, key
    assert {"head", "dirty_hash"} <= set(head["git"])


def test_the_run_header_is_the_real_paired_compare_stamp():
    from delfin.agent import paired_compare as pc
    head = bm.run_provenance(model="m")
    assert head["stamp_schema"] == pc.SCHEMA
    assert head["code"] and head["instrument"]["sha256"]


def test_the_run_header_has_the_generic_stamp_fields():
    head = bm.run_provenance(model="m")
    for key in ("taken_at", "host", "platform", "machine", "python", "code",
                "instrument", "settings", "process_start"):
        assert key in head, key
    assert {"head", "dirty_hash"} <= set(head["git"])


def test_the_run_header_uses_paired_compare_when_present(monkeypatch):
    mod = types.ModuleType("delfin.agent.paired_compare")
    calls = {}

    def provenance(**kw):
        calls.update(kw)
        return {"schema": 7, "platform": "p", "git": {"toplevel": "/t",
                                                      "dirty_hash": "d"},
                "code": {"sha256": "c"}, "env": {"DELFIN_X": "1"}}

    mod.provenance = provenance
    monkeypatch.setitem(sys.modules, "delfin.agent.paired_compare", mod)
    import delfin.agent as pkg
    monkeypatch.setattr(pkg, "paired_compare", mod, raising=False)
    head = bm.run_provenance(model="m")
    assert head["platform"] == "p" and head["env"] == {"DELFIN_X": "1"}
    assert head["stamp_schema"] == 7 and head["schema"] == 1
    assert head["git"]["dirty_hash"] == "d" and "head" in head["git"]
    assert calls["prefixes"] == ("DELFIN_",)
    assert any(p.endswith("benchmark.py") for p in calls["code"])
