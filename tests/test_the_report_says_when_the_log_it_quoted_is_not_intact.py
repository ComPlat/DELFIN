"""The changes report says when the log it quoted is not intact.

Input: an audit log written by ``audit_log.append`` under ``tmp_path`` (the
conftest points the log's default path into the test's own directory), then
altered the way a process outside DELFIN would alter it -- a line edited, a
tail cut off, the file emptied.

Output: ``audit_log.build_changes_report``'s new ``integrity`` element, the
text ``format_changes_report`` renders, and the answer the real entry point
``list_changes_made`` gives through ``_DocToolExecutor.execute``.

Semantics:

* ``integrity`` reports what ``verify`` said about the SAME bytes the records
  came from: ``{"ok", "files": [{"name", "ok", "problems", "notes"}]}``.
* an intact log says nothing extra; a log that is not intact puts
  ``NOT INTACT: <file>`` on the FIRST line, above the changes, because a list
  read under a false premise is worse than no list.
* a log with a live writer is NOT "not intact": an append writes the line and
  then its head file, so the log can be a few bytes longer than the head for
  an instant. Measured on this host before the fix: 16.6 % of reads of a log
  with a live writer reported a ``head`` problem, against 0.1 % after (all of
  them a torn last line, a separate race). ``_just_grew`` accepts
  that gap only when the bytes past the marker continue the chain.
* a rewrite that keeps the byte count but changes a line, and a log longer
  than the head whose extra bytes do NOT continue the chain, stay findings.

The baseline this is measured against: at ``origin/main`` (``14a08f5f9``)
the report had no ``integrity`` key at all and rendered the tampered records
as an ordinary list. Run against that commit, 10 of these 11 fail; the one
that passes is the never-worse guard (growth that does not continue the
chain is still a finding), which is what it is for -- it must pass on both
sides of the change.
"""
from __future__ import annotations

import inspect
import json
from pathlib import Path

from delfin.agent import audit_log
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _log(tmp_path: Path) -> Path:
    return tmp_path / "state" / "audit.log"


def _write(path: Path, n: int = 4, sid: str = "s1") -> None:
    for i in range(n):
        audit_log.append(
            audit_log.make_record(tool="write_file", decision="ok",
                                  path=f"file{i}.py", session_id=sid),
            log_path=path)


def _edit_line(path: Path, index: int, old: str, new: str) -> None:
    """Change one line in place, leaving the line count and byte edges."""
    lines = path.read_text(encoding="utf-8").splitlines()
    assert old in lines[index], lines[index]
    lines[index] = lines[index].replace(old, new, 1)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# An intact log: the report is what it was, and nothing extra is claimed
# ---------------------------------------------------------------------------

def test_an_intact_log_reports_itself_intact_and_warns_nothing(tmp_path):
    """The control for every warning below: same log, nobody touched it."""
    path = _log(tmp_path)
    _write(path)

    report = audit_log.build_changes_report("s1", log_path=path)
    text = audit_log.format_changes_report(report)

    assert report["integrity"]["ok"] is True, report["integrity"]
    assert [f["name"] for f in report["integrity"]["files"]] == ["audit.log"]
    assert "NOT INTACT" not in text
    assert "could not be checked" not in text
    # Still an answer, not a warning that swallowed the answer.
    assert [f["path"] for f in report["files_written"]] == [
        f"file{i}.py" for i in range(4)]
    assert "Files written:" in text


# ---------------------------------------------------------------------------
# A tampered log: the report says so, first
# ---------------------------------------------------------------------------

def test_an_edited_record_is_not_quoted_as_if_it_were_the_record(tmp_path):
    """A line rewritten in place changes the answer; the report must say so."""
    path = _log(tmp_path)
    _write(path)
    _edit_line(path, 1, "file1.py", "secrets.py")

    report = audit_log.build_changes_report("s1", log_path=path)
    text = audit_log.format_changes_report(report)

    assert report["integrity"]["ok"] is False, report["integrity"]
    kinds = {p["kind"] for p in report["integrity"]["files"][-1]["problems"]}
    assert "chain" in kinds, kinds
    assert text.splitlines()[0].startswith("**NOT INTACT: audit.log**"), text
    # The tampered record is still shown -- marked, not silently dropped.
    # Hiding it would answer "nothing changed" about a file that was touched.
    assert "Files written:" in text


def test_a_cut_tail_is_named_before_the_changes_it_hid(tmp_path):
    """The last two writes were removed; the report must not read as complete."""
    path = _log(tmp_path)
    _write(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(lines[:2]) + "\n", encoding="utf-8")

    report = audit_log.build_changes_report("s1", log_path=path)
    text = audit_log.format_changes_report(report)

    assert report["integrity"]["ok"] is False
    assert any(p["kind"] == "head"
               for p in report["integrity"]["files"][-1]["problems"]), report
    first = text.splitlines()[0]
    assert first.startswith("**NOT INTACT: audit.log**"), text
    assert text.index(first) < text.index("### Changes made")


def test_an_emptied_log_does_not_answer_no_recorded_changes(tmp_path):
    """Empty is the shape a wiped log has, and also the shape of a quiet week.

    The report cannot tell them apart from the records alone -- which is
    exactly why the verdict has to travel with it.
    """
    path = _log(tmp_path)
    _write(path)
    path.write_text("", encoding="utf-8")

    report = audit_log.build_changes_report("s1", log_path=path)
    text = audit_log.format_changes_report(report)

    assert report["integrity"]["ok"] is False
    assert text.startswith("**NOT INTACT: audit.log**"), text
    assert "No recorded changes" in text


# ---------------------------------------------------------------------------
# The real entry point -- through the executor, the way a model reaches it
# ---------------------------------------------------------------------------

def _perms(ws: Path, sid: str = "") -> KitToolPermissions:
    p = KitToolPermissions(workspace=ws, mode="default")
    p.task_session_id = sid
    return p


def test_list_changes_made_carries_the_warning_through_the_executor(tmp_path,
                                                                   monkeypatch):
    """The tool, not the function: this is what an agent actually reads."""
    log = tmp_path / "audit.log"
    monkeypatch.setattr(audit_log, "_default_log_path", lambda: log)
    ws = tmp_path / "ws"
    ws.mkdir()
    target = ws / "note.txt"
    target.write_text("one\n", encoding="utf-8")

    ex = _DocToolExecutor()
    perms = _perms(ws, "s-tool")
    ex.execute("write_file", {"path": str(target), "content": "two\n"}, perms)
    intact = ex.execute("list_changes_made", {}, perms)
    assert "note.txt" in intact and "NOT INTACT" not in intact, intact

    _edit_line(log, 0, "note.txt", "other.txt")

    tampered = ex.execute("list_changes_made", {}, perms)
    assert tampered.splitlines()[0].startswith("**NOT INTACT: audit.log**"), (
        tampered)
    assert "other.txt" in tampered, "the warning replaced the answer"


def test_every_file_the_report_quoted_carries_its_own_verdict(tmp_path):
    """A report quoting two files must not pass one file's verdict for both.

    The active log alone answers "no recorded changes" on the first day of a
    new ISO week, so the report falls back to the rotated file and quotes
    BOTH. An intact active log beside a tampered rotated one is the case that
    a single ``ok`` cannot express: the answer is built from a file whose
    chain does not hold, and the report has to say which file failed.
    """
    log = tmp_path / "audit.log"
    rotated = tmp_path / "audit-2026-W20.log"
    _write(rotated, 3, sid="s-old")
    _edit_line(rotated, 1, "file1.py", "secrets.py")
    _write(log, 2, sid="s1")

    # Asking for the old session: the active log holds no record of it, so
    # the report falls back to the rotated file and quotes both.
    report = audit_log.build_changes_report("s-old", log_path=log)
    integrity = report["integrity"]
    by_name = {f["name"]: f for f in integrity["files"]}

    assert set(by_name) == {"audit-2026-W20.log", "audit.log"}, by_name
    assert by_name["audit.log"]["ok"] is True, by_name
    assert by_name["audit-2026-W20.log"]["ok"] is False, by_name
    assert integrity["ok"] is False, "one bad quoted file made the report intact"
    text = audit_log.format_changes_report(report)
    assert text.splitlines()[0].startswith("**NOT INTACT"), text


# ---------------------------------------------------------------------------
# Never-worse: what must NOT be reported as tampering
# ---------------------------------------------------------------------------

def test_a_log_longer_than_its_head_because_a_write_was_in_flight(tmp_path):
    """The head is written AFTER the line; the gap between them is not a cut.

    Constructed rather than raced: the state is exactly the one an append
    leaves for an instant, and it is deterministic here.
    """
    path = _log(tmp_path)
    _write(path)
    head = audit_log._read_head(path)
    assert head is not None

    # One more chained record, with the head file left where it was -- the
    # bytes past the marker are that record, and its prev is the marker.
    audit_log.append(audit_log.make_record(tool="bash", decision="ok",
                                           command="echo late", session_id="s1"),
                     log_path=path)
    audit_log._write_head(path, head["last"], head["size"])

    report = audit_log.verify(path)
    assert report["ok"], report["problems"]
    assert not report["problems"]
    assert any("in flight" in n for n in report["notes"]), report["notes"]


def test_bytes_appended_without_the_chain_are_still_a_finding(tmp_path):
    """The other half: growth that does NOT continue the chain is tampering."""
    path = _log(tmp_path)
    _write(path)
    head = audit_log._read_head(path)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps({"tool": "bash", "decision": "ok",
                             "prev": "0" * 64}) + "\n")
    audit_log._write_head(path, head["last"], head["size"])

    report = audit_log.verify(path)
    assert not report["ok"], report
    assert any(p["kind"] in ("head", "chain") for p in report["problems"]), report


def test_a_symlinked_log_is_refused_rather_than_followed(tmp_path):
    """``_scan`` reads with ``_read_private``: the private-file check applies.

    Before this the report followed the link with ``read_text`` -- a log
    pointed at another file was read as if it were the audit log.
    """
    real = tmp_path / "elsewhere.log"
    _write(real)
    link = tmp_path / "link.log"
    link.symlink_to(real)

    report = audit_log.build_changes_report("s1", log_path=link)

    assert report["integrity"]["ok"] is False
    assert report["files_written"] == [], report["files_written"]


# ---------------------------------------------------------------------------
# The deletion: the two names that nothing reached are gone
# ---------------------------------------------------------------------------

def test_the_two_readers_of_the_log_cannot_disagree(tmp_path):
    """``read_last_n`` and the report's reader are one parser, not two.

    Measured at the baseline: a log with a line holding ``[1, 2, 3]`` and one
    holding ``42`` made ``read_last_n`` return ``['dict', 'list', 'int']``
    where the report's reader returned ``['dict']`` -- the same file, two
    answers. Whatever they return now, it has to be the same answer.
    """
    path = _log(tmp_path)
    _write(path, 2)
    with path.open("a", encoding="utf-8") as fh:
        fh.write("[1, 2, 3]\n")          # valid JSON, not a record
        fh.write("42\n")                 # valid JSON, not a record

    tail = audit_log.read_last_n(10, log_path=path)
    # _scan is the one reader: records and the verification of the
    # same bytes. The _parse_log wrapper this line used to call was
    # left behind by the consolidation -- its docstring claimed two
    # call sites read a rotated log through it, and both of them
    # call _scan directly, so nothing in delfin/ named it and only
    # this test kept it alive.
    whole, _ = audit_log._scan(path)

    assert tail == whole, (tail, whole)
    assert all(isinstance(r, dict) for r in tail), tail
    # n is a count, not a slice bound that wraps: n=0 is no records.
    assert audit_log.read_last_n(0, log_path=path) == []


def test_the_unreached_diff_helpers_are_gone():
    """Nothing called them, so their behaviour is not covered by a test.

    ``test_a_value_that_is_not_shipped_has_no_value_without_its_table`` is the
    precedent: a name that must not come back is pinned by its absence.
    """
    assert not hasattr(audit_log, "_diff_line_counts")

    params = inspect.signature(audit_log.make_record).parameters
    assert "diff_lines_added" not in params
    assert "diff_lines_removed" not in params
