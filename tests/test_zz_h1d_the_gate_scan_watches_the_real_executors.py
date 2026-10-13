"""The gate scan runs on the real tool executors and watches for change.

``gate_coverage`` reads ``api_client.py`` and asks, for every ``_execute_*``
method that touches the filesystem, whether a read or write gate runs before
the access. On its first real run it pointed at ``grep_file``, which did read
files outside the workspace through a symlink (fixed separately).

The scan is not finished: on today's file it flags eight executors, and for
most of them the flag is the scan's own imprecision (a sink inside a shared
helper is charged to the wrong executor; a read-only git call on a gated
folder). Those eight are listed below as OPEN -- not as approved exceptions.
What this file makes of the scan today is an alarm in two directions:

* an executor that is not on the list and gets flagged fails the test -- a new
  tool, or a change to an old one, that reads or writes before a gate;
* an executor on the list that is no longer flagged fails the test too, so
  the list shrinks when the scan or the code improves and cannot rot.

And it pins that the scan is not blind: for seven executors that are clean
today, taking their gate calls out of the source makes the scan flag them.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.agent import gate_coverage as gc

SOURCE = Path(gc.__file__).with_name("api_client.py")

#: Flagged on the real file today. OPEN means nobody has shown the flag to be
#: harmless; it is a finding to work through, kept visible here.
OPEN = {
    "_execute_grep_file": "real: read a symlink's target outside the searched "
                          "folder; fixed in the tool, the scan does not read "
                          "the containment check as a gate yet",
    "_execute_pinned": "open: sinks of shared helpers are charged here",
    "_execute_list_files": "open: lists names under a gated folder",
    "_execute_run_tests": "open: unverified",
    "_execute_apply_patch": "open: reads pre/post images under "
                            "_run_permission_gate",
    "_execute_exit_worktree": "open: git on a gated folder",
    "_execute_worktree_merge": "open: git on a gated folder",
    "_execute_make_plot": "open: unverified",
}

#: Clean today, and proven to be watched: without their gate calls the scan
#: flags them.
WATCHED = (
    "_execute_read_file", "_execute_edit_sheet", "_execute_fill_pdf_form",
    "_execute_fill_docx_template", "_execute_notebook_edit",
    "_execute_bash", "_execute_bash_background",
)

_GATE_CALL = re.compile(
    r"self\.(_check_read_access|_gate_write_path|_resolve_in_workspace|"
    r"_run_permission_gate|_gate_[a-z_]+|_check_read_[a-z_]+)\(")


@pytest.fixture(scope="module")
def source() -> str:
    return SOURCE.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def report(source):
    return gc.analyze(source, str(SOURCE))


def test_the_scan_sees_the_executors_and_their_file_access(report):
    assert len(report.executors) >= 60
    assert sum(1 for e in report.executors if e.sinks) >= 15


def test_no_executor_outside_the_open_list_is_flagged(report):
    flagged = {e.name for e in report.ungated()}
    new = sorted(flagged - set(OPEN))
    assert not new, (
        f"{new}: file access with no read or write gate before it. Gate it, "
        "or show why the scan is wrong and say so in OPEN.")


def test_the_open_list_holds_nothing_that_is_clean_now(report):
    flagged = {e.name for e in report.ungated()}
    stale = sorted(set(OPEN) - flagged)
    assert not stale, f"{stale}: no longer flagged, take them off OPEN"


def test_the_open_list_names_only_executors_that_exist(report):
    assert not report.stale_exemptions(OPEN)


@pytest.mark.parametrize("name", WATCHED)
def test_taking_the_gate_out_of_a_clean_executor_gets_it_flagged(
        name, source, report):
    finding = next(e for e in report.executors if e.name == name)
    assert finding.sinks and not finding.ungated_sinks, (
        f"{name} is expected to be clean on the real file")
    lines = source.split("\n")
    for line_no in sorted({g.line for g in finding.gates}):
        lines[line_no - 1] = _GATE_CALL.sub("self._NOGATE(", lines[line_no - 1])
    flagged = {e.name for e in gc.analyze("\n".join(lines), "x").ungated()}
    # Others may be flagged with it: several executors share one gate that
    # sits in a common helper, and taking it out uncovers all of them.
    assert name in flagged
