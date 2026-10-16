"""The gate scan on the real ``api_client.py``.

Input: the source of ``delfin/agent/api_client.py``.  Semantics: every
``_execute_*`` method whose file access has no gate before it fails the
test, unless

* the access is listed in ``EXPLAINED`` by executor, function and call,
  with the reason it needs no gate (a path the harness chooses, not the
  model), or
* the executor is listed in ``OPEN``: flagged, not shown harmless or
  shown harmless only in a way the scan cannot read, one line each.

Both lists are checked in the other direction too: an entry that no longer
matches anything fails, so neither can rot.  ``WATCHED`` pins that the scan
is not blind for the executors that are clean: with their gate calls taken
out of the source, each is flagged.

Cost: one parse of the file, about 0.4 s; each WATCHED case parses again.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.agent import gate_coverage as gc

SOURCE = Path(gc.__file__).with_name("api_client.py")

_SIBLINGS = ("existence of the same relative path in the repository's own "
             ".delfin/worktrees, after the gate passed and the file was "
             "missing; reached from the not-found message")
_GITIGNORE = "reads <workspace>/.gitignore; root is perms.workspace"
_GREP_WALK = ("files of the walk of search_path, which _check_read_access "
              "gated; a file symlink whose target leaves the searched folder "
              "is skipped and one to a secret name is skipped, before "
              "_looks_binary or open; os.walk does not follow directory "
              "links. The walk itself is a flow the scan does not follow")

_OWNED_SOURCE = ("`git show HEAD:./<name>` of the path the secret-deny "
                 "check is deciding, to compare its HEAD blob with the "
                 "session's journalled pre-image; the bytes are hashed, "
                 "never returned, and a False only keeps the refusal")

_REVIEWED = ("`git ls-files --error-unmatch` and `git status --porcelain` "
             "of the path the secret-deny check is deciding (the same "
             "call site as _session_owned_source), through safe_git; only "
             "git's answer is used, and a False keeps the refusal. Seen "
             "by the scan since it reads safe_git calls; the earlier "
             "`import subprocess as _sp` spelling was not")

#: (executor, function holding the access, sink kind) -> why it needs no
#: gate.  Each was driven through the real executor in a tmp workspace.
EXPLAINED = {
    ("_execute_pinned", "_ensure_loaded", "read_text"):
        "the documentation index at get_default_index_path()",
    ("_execute_read_file", "_on_another_sessions_branch", "iterdir"):
        _SIBLINGS,
    ("_execute_edit_file", "_on_another_sessions_branch", "iterdir"):
        _SIBLINGS,
    ("_execute_multi_edit", "_on_another_sessions_branch", "iterdir"):
        _SIBLINGS,
    ("_execute_notebook_read", "_on_another_sessions_branch", "iterdir"):
        _SIBLINGS,
    ("_execute_notebook_edit", "_on_another_sessions_branch", "iterdir"):
        _SIBLINGS,
    ("_execute_grep_file", "_session_owned_source", "safe_git.run"):
        _OWNED_SOURCE,
    ("_execute_grep_file", "_is_reviewed_project_file", "safe_git.run"):
        _REVIEWED,
    ("_execute_apply_patch", "_session_owned_source", "safe_git.run"):
        _OWNED_SOURCE,
    ("_execute_apply_patch", "_is_reviewed_project_file", "safe_git.run"):
        _REVIEWED,
    ("_execute_code_nav", "_session_owned_source", "safe_git.run"):
        _OWNED_SOURCE,
    ("_execute_code_nav", "_is_reviewed_project_file", "safe_git.run"):
        _REVIEWED,
    ("_execute_grep_file", "_gitignore_skip_dirs", "read_text"): _GITIGNORE,
    ("_execute_list_files", "_gitignore_skip_dirs", "read_text"): _GITIGNORE,
    ("_execute_list_files", "_iter_scan_files", "os.walk"):
        "walks perms.workspace or the folder _resolve_in_workspace returned "
        "(reached through `base = resolved`, which the scan does not "
        "follow); an outside path and a symlink to outside are refused",
    ("_execute_run_tests", "_git_metadata_roots", "safe_git.run"):
        "git rev-parse in perms.workspace",
    ("_execute_grep_file", "_looks_binary", "open"): _GREP_WALK,
    ("_execute_grep_file", "_looks_binary", "read"): _GREP_WALK,
    ("_execute_grep_file", "_execute_grep_file", "open"): _GREP_WALK,
    ("_execute_apply_patch", "_note_edit_signals_raw", "read_bytes"):
        "the post-image of a file the write path has just changed and is "
        "journalling (_capture_raw_change, after _gate_write_path ran on "
        "it), at most 400 kB, read to classify the change for "
        "edit_signals; the bytes are never returned",
    ("_execute_apply_patch", "_execute_apply_patch", "read_bytes"):
        "pre-images under perms.workspace, read after _run_permission_gate "
        "ran _gate_write_path on every file of the diff; both sit under "
        "`if not check_only:`, a correlation the scan does not see. A "
        "`../` path and a symlink to outside are refused in both modes",
}

#: Executor -> why it is still flagged.
OPEN: dict[str, str] = {
    "_execute_task_update":
        "_check_completion reads the session's own tool trace "
        "(tool_trace.read(task_session_id), a harness path with a "
        "sanitised id); not driven through the executor here yet",
}

#: Executors allowed to be flagged or not while a pending change lands.
#: Empty: an entry here hides a flag, so a known one belongs in OPEN.
EXPECTED_SOON: set[str] = set()

#: Clean, with file access, and flagged once their gate calls are removed.
WATCHED = (
    "_execute_read_file", "_execute_edit_sheet", "_execute_fill_pdf_form",
    "_execute_fill_docx_template", "_execute_write_file",
    "_execute_edit_file", "_execute_multi_edit", "_execute_notebook_read",
    "_execute_notebook_edit", "_execute_bash", "_execute_bash_background",
    "_execute_exit_worktree",
    "_execute_worktree_merge", "_execute_make_plot",
)

_GATE_CALL = re.compile(r"self\.(%s)\(" % "|".join(gc.GATES))


def _unexplained(report) -> dict[str, list]:
    out: dict[str, list] = {}
    for ex in report.executors:
        left = [r for r in ex.ungated_sinks
                if (ex.name, r.sink.owner, r.sink.kind) not in EXPLAINED]
        if left:
            out[ex.name] = left
    return out


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
    new = {name: [f"{r.sink.owner}:{r.sink.line} {r.sink.kind}"
                  for r in rows]
           for name, rows in _unexplained(report).items()
           if name not in OPEN and name not in EXPECTED_SOON}
    assert not new, (
        f"{new}: file access with no read or write gate before it. Gate it, "
        "or show why it needs none and say so in EXPLAINED or OPEN.")


def test_the_open_list_holds_nothing_that_is_clean_now(report):
    stale = sorted(set(OPEN) - set(_unexplained(report)))
    assert not stale, f"{stale}: no longer flagged, take them off OPEN"
    assert not report.stale_exemptions(OPEN)


def test_every_explained_access_still_exists_ungated(report):
    seen = {(ex.name, r.sink.owner, r.sink.kind)
            for ex in report.executors for r in ex.ungated_sinks}
    assert not set(EXPLAINED) - seen


@pytest.mark.parametrize("name", WATCHED)
def test_taking_the_gate_out_of_a_clean_executor_gets_it_flagged(
        name, source, report):
    finding = next(e for e in report.executors if e.name == name)
    assert finding.sinks and name not in _unexplained(report), (
        f"{name} is expected to be clean on the real file")
    lines = source.split("\n")
    # The gate lines it reaches in helpers, and every gate in its own body
    # (a nested callback's gate is the callback's, not in finding.gates).
    own = range(finding.line, finding.end_line + 1)
    for line_no in {g.line for g in finding.gates} | set(own):
        lines[line_no - 1] = _GATE_CALL.sub("self._NOGATE(", lines[line_no - 1])
    # Others may be flagged with it where they share a gate in a helper.
    assert name in _unexplained(gc.analyze("\n".join(lines), "x"))
