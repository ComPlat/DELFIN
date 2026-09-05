"""A size/mtime diff of a shared tree cannot name a culprit.

``run_suite`` fingerprints the checkout before the tasks and diffs it
after, then prints "this run changed N file(s) OUTSIDE the fixture
workspaces". The diff is true about the tree and the sentence is a claim
about the run, and the two are not the same thing whenever anything else
writes into the checkout while the suite is going.

Observed 2026-08-13: an office run reported 569 changed files, every one
of the shape ``MagicMock/mock._permissions.workspace/<id>/.delfin/...``.
The benchmark never writes such a path. They were transient leaks of a
pytest suite running in the same checkout, created after the fingerprint
was taken and gone before the diff. 569 lines of noise is also how a real
warning gets skipped.

The run does have a record of its own writes -- the audit log, the same
source ``_denials_during`` already reads -- so the report can name what
it actually touched and say plainly that the rest is unattributable.

``None`` for that record means the log could not be read. It must not
collapse to "wrote nothing", which would relabel every changed file as
somebody else's.
"""

from __future__ import annotations

import pytest

from delfin.agent import benchmark_runner as br


# ---------------------------------------------------------------------------
# Splitting the diff by what the trace can vouch for
# ---------------------------------------------------------------------------

def test_a_file_this_run_wrote_is_named_as_this_run_s():
    mine, other = br.attributed_changes(
        ["delfin/agent/engine.py"], {"delfin/agent/engine.py"})
    assert mine == ["delfin/agent/engine.py"]
    assert other == []


def test_a_file_nothing_in_the_trace_wrote_is_kept_apart():
    mine, other = br.attributed_changes(
        ["MagicMock/mock._permissions.workspace/x/.delfin/bash_jobs.json.lock"],
        set())
    assert mine == []
    assert other == [
        "MagicMock/mock._permissions.workspace/x/.delfin/bash_jobs.json.lock"]


def test_nothing_is_counted_twice():
    changed = ["a.py", "b.py", "c.py"]
    mine, other = br.attributed_changes(changed, {"b.py"})
    assert sorted(mine + other) == sorted(changed)
    assert set(mine) & set(other) == set()


def test_a_write_that_changed_nothing_is_not_invented_into_the_report():
    """The agent rewrote a file with identical bytes, or wrote inside a
    fixture workspace the guard then restored. Neither is a change to the
    checkout and neither belongs in a report about one."""
    mine, other = br.attributed_changes(["a.py"], {"a.py", "b.py", "c.py"})
    assert mine == ["a.py"]
    assert other == []


def test_an_unreadable_trace_is_not_an_empty_one():
    """``None`` means nobody looked. Treating it as "wrote nothing" would
    move every changed file into the other column and accuse a process
    that may not exist."""
    mine, other = br.attributed_changes(["a.py"], None)
    assert mine is None
    assert other is None


# ---------------------------------------------------------------------------
# What the run prints
# ---------------------------------------------------------------------------

def test_a_clean_checkout_says_nothing():
    assert br.format_checkout_change_report([], set()) == ""
    assert br.format_checkout_change_report([], None) == ""


def test_the_run_s_own_writes_are_reported_as_its_own():
    text = br.format_checkout_change_report(["delfin/agent/engine.py"],
                                            {"delfin/agent/engine.py"})
    assert "this run wrote" in text
    assert "delfin/agent/engine.py" in text


def test_what_it_cannot_attribute_is_reported_as_exactly_that():
    text = br.format_checkout_change_report(["stray.json"], set())
    assert "stray.json" in text
    # The sentence that was missing: the run does not claim these.
    assert "nothing in this run" in text.lower()
    assert "another process" in text.lower()


def test_the_two_groups_are_not_merged_into_one_number():
    text = br.format_checkout_change_report(
        ["mine.py", "theirs.json"], {"mine.py"})
    assert "1 file(s) this run wrote" in text
    assert "1 file(s) changed" in text


def test_without_a_trace_the_report_names_its_own_limit():
    """The old message, plus the sentence it was missing. It must NOT
    claim the run made the changes when it could not check."""
    text = br.format_checkout_change_report(["a.py", "b.py"], None)
    assert "2 file(s)" in text
    assert "could not be read" in text or "cannot tell" in text
    assert "this run wrote" not in text


def test_a_long_list_is_capped_and_says_how_many_it_dropped():
    """569 lines is not a warning anybody reads twice."""
    changed = [f"f{i}.json" for i in range(60)]
    text = br.format_checkout_change_report(changed, set())
    assert "… and 40 more" in text
    assert text.count("\n") < 40


def test_the_run_s_own_writes_are_never_silently_dropped_by_the_cap():
    """The cap must not hide the group that actually matters behind 20
    lines of the group that does not."""
    changed = [f"noise{i}.json" for i in range(40)] + ["real.py"]
    text = br.format_checkout_change_report(changed, {"real.py"})
    assert "real.py" in text


# ---------------------------------------------------------------------------
# Where the record comes from
# ---------------------------------------------------------------------------

def test_the_written_paths_come_back_relative_to_the_checkout(tmp_path,
                                                              monkeypatch):
    """The fingerprint is keyed by relative path; the audit log records
    absolute ones. A comparison across the two would match nothing and
    quietly report every change as unattributable."""
    root = tmp_path / "checkout"
    (root / "delfin").mkdir(parents=True)

    def _fake_report(**kw):
        return {"files_written": [{"path": str(root / "delfin" / "x.py")}],
                "window": {"records": 1}}

    monkeypatch.setattr("delfin.agent.audit_log.build_changes_report",
                        _fake_report)
    got = br.paths_this_run_wrote(root, since_ts="2026-01-01T00:00:00")
    assert got == {"delfin/x.py"}


def test_a_write_outside_the_checkout_is_not_mapped_into_it(tmp_path,
                                                            monkeypatch):
    root = tmp_path / "checkout"
    root.mkdir()

    def _fake_report(**kw):
        return {"files_written": [{"path": "/etc/passwd"},
                                  {"path": str(root / "in.py")}],
                "window": {"records": 2}}

    monkeypatch.setattr("delfin.agent.audit_log.build_changes_report",
                        _fake_report)
    assert br.paths_this_run_wrote(root, since_ts="x") == {"in.py"}


def test_an_empty_audit_window_is_unreadable_not_empty(tmp_path, monkeypatch):
    """The skeleton that comes back on failure has zero records and an
    empty list — the same shape as a run that wrote nothing."""
    monkeypatch.setattr(
        "delfin.agent.audit_log.build_changes_report",
        lambda **kw: {"files_written": [], "window": {"records": 0}})
    assert br.paths_this_run_wrote(tmp_path, since_ts="x") is None


def test_a_window_that_saw_records_can_say_the_run_wrote_nothing(
        tmp_path, monkeypatch):
    monkeypatch.setattr(
        "delfin.agent.audit_log.build_changes_report",
        lambda **kw: {"files_written": [], "window": {"records": 12}})
    assert br.paths_this_run_wrote(tmp_path, since_ts="x") == set()


def test_a_broken_audit_log_does_not_take_the_run_down(tmp_path, monkeypatch):
    def _boom(**kw):
        raise RuntimeError("log unreadable")

    monkeypatch.setattr("delfin.agent.audit_log.build_changes_report", _boom)
    assert br.paths_this_run_wrote(tmp_path, since_ts="x") is None


# ---------------------------------------------------------------------------
# The suite uses it
# ---------------------------------------------------------------------------

def test_the_suite_prints_the_attributed_report(tmp_path, monkeypatch, capsys):
    """The mechanism has to be on the path the run actually takes."""
    seen: dict = {}

    def _fake_format(changed, written):
        seen["changed"], seen["written"] = list(changed), written
        return "REPORT-MARKER"

    monkeypatch.setattr(br, "format_checkout_change_report", _fake_format)
    monkeypatch.setattr(br, "changed_outside_workspaces",
                        lambda root, before: ["stray.json"])
    monkeypatch.setattr(br, "checkout_fingerprint", lambda root: {"a": (1, 1)})
    monkeypatch.setattr(br, "paths_this_run_wrote",
                        lambda root, since_ts: {"other.py"})

    monkeypatch.chdir(tmp_path)
    br.run_suite([], model="test-model")
    out = capsys.readouterr().out
    assert "REPORT-MARKER" in out
    assert seen["changed"] == ["stray.json"]
    assert seen["written"] == {"other.py"}


def test_a_suite_that_changed_nothing_prints_nothing(tmp_path, monkeypatch,
                                                     capsys):
    monkeypatch.setattr(br, "changed_outside_workspaces",
                        lambda root, before: [])
    monkeypatch.setattr(br, "checkout_fingerprint", lambda root: {"a": (1, 1)})
    monkeypatch.chdir(tmp_path)
    br.run_suite([], model="test-model")
    assert "OUTSIDE" not in capsys.readouterr().out


@pytest.mark.parametrize("name", ["attributed_changes",
                                  "format_checkout_change_report",
                                  "paths_this_run_wrote"])
def test_the_pieces_are_importable_by_name(name):
    assert callable(getattr(br, name))
