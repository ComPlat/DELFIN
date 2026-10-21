"""The gate scan on fixture source: which accesses it sees, and when a gate
counts as running before them.

``gate_coverage.analyze`` parses source text and never imports it, so every
case here is a few lines of class text.  Each case is either a sink that
must be reported (``_execute_x`` in the flagged set) or its gated mirror
that must not be.
"""
from __future__ import annotations

import pytest

from delfin.agent import gate_coverage as gc


def _src(*body: str, helpers: tuple[str, ...] = ()) -> str:
    lines = ["import os, shutil, subprocess, tempfile",
             "from pathlib import Path",
             "class T:"]
    lines += ["    " + h for h in helpers]
    lines.append("    def _execute_x(self, arguments, perms=None):")
    lines += ["        " + b for b in body]
    return "\n".join(lines) + "\n"


def _flagged(src: str) -> set[str]:
    return {e.name for e in gc.analyze(src).ungated()}


GATE_P = "self._check_read_access(perms, p)"


# ---------------------------------------------------------------------------
# What is a sink
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("access", [
    "open(p).read()",
    "os.open(p, os.O_RDONLY)",
    "p.read_text()",
    "p.write_text('x')",
    "p.open('w')",
    "p.unlink()",
    "p.mkdir(parents=True)",
    "p.rmdir()",
    "p.touch()",
    "p.rename(arguments['q'])",
    "list(p.iterdir())",
    "list(p.rglob('*'))",
    "os.listdir(p)",
    "os.replace(p, arguments['q'])",
    "shutil.copy(p, arguments['q'])",
    "shutil.move(p, arguments['q'])",
    "shutil.rmtree(p)",
    "subprocess.run(['cat', arguments['path']])",
    "subprocess.run(['cat', str(p)])",
    "subprocess.run(arguments['cmd'])",
    "subprocess.run(['ls'], cwd=p)",
    "subprocess.run(['ls'], executable=p)",
    "_text_files.write_text_file(p, 'x')",
    "_text_files.read_text_file(p)",
    "safe_git.run(['ls-files'], p)",
    "safe_git.run(['show', str(p)], '.')",
    "safe_git.check_output(['status'], cwd=p)",
    "safe_git.out(p, 'status')",
])
def test_an_ungated_access_is_reported(access):
    assert "_execute_x" in _flagged(_src("p = arguments['path']", access))


def test_an_argv_of_constants_and_refs_is_no_file_access():
    assert _flagged(_src("ref = arguments.get('ref') or 'HEAD'",
                         "subprocess.run(['git', 'log', ref])")) == set()


def test_safe_git_with_a_constant_directory_and_refs_is_no_file_access():
    assert _flagged(_src("ref = arguments.get('ref') or 'HEAD'",
                         "safe_git.run(['log', ref], '.')")) == set()


def test_a_function_that_is_not_an_executor_is_not_reported():
    report = gc.analyze(_src("return 1", helpers=(
        "def _not_an_executor(self, a):",
        "    return open(a['path'])")))
    assert report.names() == ["_execute_x"]
    assert report.ungated() == []


# ---------------------------------------------------------------------------
# Order and path
# ---------------------------------------------------------------------------

def test_a_gate_before_the_access_on_the_same_path_covers_it():
    assert _flagged(_src("p = arguments['path']", GATE_P,
                         "with open(p) as f:", "    return f.read()")) == set()


def test_a_gate_after_the_access_does_not_cover_it():
    report = gc.analyze(_src("p = arguments['path']",
                             "text = open(p).read()", GATE_P))
    (ex,) = report.ungated()
    assert all(r.gate_line is None for r in ex.ungated_sinks)


@pytest.mark.parametrize("access", [
    "open(b).read()", "b.write_text('x')", "b.unlink()",
    "os.replace(p, b)", "shutil.copy(p, b)",
])
def test_a_gate_on_another_path_does_not_cover_it(access):
    assert "_execute_x" in _flagged(_src(
        "p = arguments['path']", "b = arguments['other']", GATE_P, access))


def test_str_and_path_wrappers_name_the_same_path():
    assert _flagged(_src(
        "wt = Path(arguments['path'])",
        "err = self._worktree_path_refusal('path', wt, perms)",
        "subprocess.check_output(['git', '-C', str(wt), 'status'])")) == set()


def test_the_path_a_gate_returns_counts_as_checked():
    assert _flagged(_src(
        "resolved, err = self._resolve_in_workspace(arguments['path'], perms)",
        "if err:", "    return err",
        "resolved.write_text('x')")) == set()


def test_a_handle_opened_from_a_path_method_is_that_path():
    assert _flagged(_src(
        "p = arguments['path']",
        "if perms is not None:", "    " + GATE_P,
        "with p.open() as fh:", "    return fh.read()")) == set()


def test_a_gate_that_is_not_a_path_gate_covers_nothing():
    assert "_execute_x" in _flagged(_src(
        "p = arguments['path']", "self._gate_mcp_tool('x', arguments)",
        "return open(p)"))


def test_a_gate_without_a_path_position_covers_every_path():
    assert _flagged(_src(
        "self._run_permission_gate('x', arguments, perms)",
        "Path(arguments['path']).write_text('x')")) == set()


# ---------------------------------------------------------------------------
# Control flow
# ---------------------------------------------------------------------------

def test_a_gate_in_one_branch_does_not_cover_the_other():
    for gate_branch, sink_branch in ((GATE_P, "return open(p)"),
                                     ("return open(p)", GATE_P)):
        assert "_execute_x" in _flagged(_src(
            "p = arguments['path']", "if perms is not None:",
            "    " + gate_branch, "else:", "    " + sink_branch))


@pytest.mark.parametrize("head", [
    ("if arguments.get('strict'):", "    " + GATE_P),
    ("if perms is not None and arguments.get('strict'):", "    " + GATE_P),
    ("for _ in arguments.get('items', []):", "    " + GATE_P),
    ("try:", "    " + GATE_P, "except Exception:", "    pass"),
])
def test_a_gate_under_a_condition_does_not_cover_what_follows(head):
    assert "_execute_x" in _flagged(_src(
        "p = arguments['path']", *head, "return open(p)"))


def test_a_gate_under_perms_is_not_none_covers_what_follows():
    assert _flagged(_src(
        "p = arguments['path']", "if perms is not None:", "    " + GATE_P,
        "if p.is_dir():", "    return list(p.iterdir())",
        "return open(p)")) == set()


def test_a_gate_and_access_in_one_branch_cover():
    assert _flagged(_src(
        "p = arguments['path']", "if arguments.get('x'):",
        "    " + GATE_P, "    return open(p)")) == set()


def test_a_gate_in_a_nested_function_does_not_cover_the_outer_access():
    assert "_execute_x" in _flagged(_src(
        "p = arguments['path']",
        "def _later():", "    " + GATE_P,
        "return open(p)"))


def test_a_callback_is_judged_by_its_own_gate():
    assert _flagged(_src(
        "def _reader(q):",
        "    self._check_read_access(perms, q)",
        "    with q.open() as fh:", "        return fh.read()",
        "return arguments.get('render')(_reader)")) == set()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def test_an_access_in_a_helper_is_charged_to_the_executor():
    assert "_execute_x" in _flagged(_src(
        "return self._read_all(arguments['path'])",
        helpers=("def _read_all(self, p):", "    return open(p)")))


def test_a_gate_on_the_argument_covers_the_helper_access():
    assert _flagged(_src(
        "p = arguments['path']", GATE_P, "return self._read_all(p)",
        helpers=("def _read_all(self, q):", "    return open(q).read()"))) \
        == set()


def test_another_executor_is_not_charged_to_its_caller():
    flagged = _flagged(_src(
        "return self._execute_y(arguments)",
        helpers=("def _execute_y(self, arguments):",
                 "    return open(arguments['path']).read()")))
    assert flagged == {"_execute_y"}


def test_the_reads_inside_a_gate_are_not_the_executors():
    assert _flagged(_src(
        "t = arguments['path']",
        "err = self._gate_write_path(t, perms, 'x', arguments)",
        "t.write_text('x')",
        helpers=("def _gate_write_path(self, path, perms, name, args):",
                 "    return self._preview(path)",
                 "def _preview(self, q):",
                 "    return q.read_text()"))) == set()


# ---------------------------------------------------------------------------
# Tempfile staging
# ---------------------------------------------------------------------------

def test_a_mkstemp_staging_file_is_no_user_access():
    assert _flagged(_src(
        "_fd, _name = tempfile.mkstemp(suffix='.txt')",
        "with os.fdopen(_fd, 'w') as _fh:", "    _fh.write('staged')",
        "_p = Path(_name)", "_p.read_bytes()", "_p.unlink()")) == set()


def test_a_named_temporary_file_is_no_user_access():
    assert _flagged(_src(
        "with tempfile.NamedTemporaryFile(mode='w') as fh:",
        "    fh.write('x')")) == set()


def test_fdopen_of_a_descriptor_not_from_tempfile_is_an_access():
    assert "_execute_x" in _flagged(_src(
        "fd = arguments['fd']",
        "with os.fdopen(fd, 'w') as fh:", "    fh.write('x')"))


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------

def test_a_name_on_a_list_that_is_no_executor_is_stale():
    report = gc.analyze(_src("return 1"))
    assert report.stale_exemptions(["_execute_ghost", "_execute_x"]) == [
        "_execute_ghost"]


def test_the_scan_does_not_import_the_module_it_reads():
    own = open(gc.__file__, encoding="utf-8").read()
    assert "import api_client" not in own and "from .api_client" not in own
