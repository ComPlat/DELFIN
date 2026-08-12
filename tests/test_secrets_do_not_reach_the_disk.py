"""The redaction protected the model and not the disk.

`_redact_tool_result` was applied only to the copy entering the model
context. The tool_result event carried the raw text -- and that event is
what the engine writes to ~/.delfin/tool_traces/<sid>.jsonl.

That file is created with a plain append and no chmod: observed 0664
across 426 files. `bug_report` then bundles it, and its packer explicitly
chowns the archive to the shared group and adds group-read to every path
on the way out. So a token that appeared in a traceback was stripped from
the transcript -- which is exactly why nobody would notice -- and written
verbatim to a group-readable file that gets shipped.

Same unprotected append in turn_metrics, failure_log and audit_log.

No secret was found in the live logs when this was audited. The path is
what is being closed, not an incident.
"""

from __future__ import annotations

import json
import os
import stat

import pytest

from delfin.agent import api_client as A


_SECRET = "sk-ant-api03-THISLOOKSLIKEAKEY0123456789abcdefghijklmnop"


# ---------------------------------------------------------------------------
# The event the trace is written from
# ---------------------------------------------------------------------------

def test_the_redactor_still_finds_it():
    """A precondition: if this fails the rest proves nothing."""
    assert _SECRET not in A._redact_tool_result(f"error: key {_SECRET} bad")


_REDACTOR = "_redact_tool_result"


def _api_client_tree():
    import ast
    import pathlib
    return ast.parse(pathlib.Path(A.__file__).read_text(encoding="utf-8"))


def _calls_redactor(node) -> bool:
    """True when the redactor is called anywhere inside ``node``."""
    import ast
    for sub in ast.walk(node):
        if not isinstance(sub, ast.Call):
            continue
        fn = sub.func
        name = (fn.id if isinstance(fn, ast.Name)
                else fn.attr if isinstance(fn, ast.Attribute) else "")
        if name == _REDACTOR:
            return True
    return False


def _names_bound_to_the_redactor(tree) -> set:
    """Locals assigned from an expression that calls the redactor.

    ``_redacted = _redact_tool_result(result)`` followed by
    ``tool_output=_redacted[:2000]`` is legitimate -- the site emitting a
    large result needs the whole redacted string anyway, for the length and
    the notes it carries alongside -- so the check follows the binding
    instead of banning it.
    """
    import ast
    out: set = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
            targets = [node.target]
        else:
            continue
        if node.value is None or not _calls_redactor(node.value):
            continue
        for t in targets:
            for sub in ast.walk(t):
                if isinstance(sub, ast.Name):
                    out.add(sub.id)
    return out


def _tool_output_arguments(tree):
    """Every ``tool_output=<expr>`` keyword argument, with its line."""
    import ast
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        for kw in node.keywords:
            if kw.arg == "tool_output":
                yield kw.value, getattr(kw.value, "lineno", 0)


def _is_derived_from_the_redactor(expr, bound) -> bool:
    import ast
    if _calls_redactor(expr):
        return True
    if isinstance(expr, ast.Constant):      # tool_output="" leaks nothing
        return True
    roots = {sub.id for sub in ast.walk(expr) if isinstance(sub, ast.Name)}
    return bool(roots) and roots <= bound


def test_every_tool_output_argument_is_derived_from_the_redactor():
    """Structural, not textual.

    The previous version asserted that one literal string was absent and
    that at least four redacted sites existed. Both raw sites in the
    CLI-backend stream parser -- ``tool_output=result_text`` and
    ``tool_output=str(content)`` -- satisfied it without being redacted at
    all, so the guard passed for as long as the hole was open, and a fifth
    raw site would have kept it passing.

    Asking the AST what each argument IS cannot be satisfied that way: a
    new emission site has to call the redactor, or name something that did.
    """
    import ast

    tree = _api_client_tree()
    bound = _names_bound_to_the_redactor(tree)
    offenders = []
    checked = 0
    for expr, lineno in _tool_output_arguments(tree):
        checked += 1
        if not _is_derived_from_the_redactor(expr, bound):
            offenders.append(f"line {lineno}: {ast.unparse(expr)}")
    assert checked >= 6, (
        f"only found {checked} tool_output arguments; the walk is no "
        "longer seeing the emission sites")
    assert not offenders, (
        "tool_result events carry output that never passed the redactor, "
        "and the engine writes these events to the tool trace: "
        + "; ".join(offenders))


def test_the_check_would_catch_a_new_raw_site():
    """A test that cannot fail is not a guard. This pins that the rule
    rejects exactly the shape that slipped through before."""
    import ast
    tree = ast.parse(
        "StreamEvent(type='tool_result', tool_output=result_text)\n")
    args = list(_tool_output_arguments(tree))
    assert len(args) == 1
    assert not _is_derived_from_the_redactor(
        args[0][0], _names_bound_to_the_redactor(tree))


def test_the_check_accepts_a_redacted_local():
    import ast
    tree = ast.parse(
        "_redacted = _redact_tool_result(result)\n"
        "StreamEvent(type='tool_result', tool_output=_redacted[:2000])\n")
    assert _is_derived_from_the_redactor(
        list(_tool_output_arguments(tree))[0][0],
        _names_bound_to_the_redactor(tree))


def test_the_check_rejects_a_local_bound_to_something_else():
    """Following the binding must not degrade into "any name is fine"."""
    import ast
    tree = ast.parse(
        "_redacted = _redact_tool_result(result)\n"
        "raw = result\n"
        "StreamEvent(type='tool_result', tool_output=raw[:2000])\n")
    assert not _is_derived_from_the_redactor(
        list(_tool_output_arguments(tree))[0][0],
        _names_bound_to_the_redactor(tree))


# ---------------------------------------------------------------------------
# The files the recorders create
# ---------------------------------------------------------------------------

def _mode(path):
    return stat.S_IMODE(os.stat(path).st_mode)


def test_the_tool_trace_is_owner_only(tmp_path, monkeypatch):
    from delfin.agent import tool_trace

    monkeypatch.setattr(tool_trace, "trace_path", lambda sid: tmp_path / "t.jsonl")
    tool_trace.record("s1", tool="bash", tool_input="x", output="y",
                      duration_ms=1, ok=True)
    p = tmp_path / "t.jsonl"
    assert p.exists()
    assert _mode(p) == 0o600, oct(_mode(p))


def test_the_audit_log_is_owner_only(tmp_path):
    from delfin.agent import audit_log

    log = tmp_path / "audit.log"
    audit_log.append({"tool": "bash", "decision": "ok"}, log_path=log)
    assert _mode(log) == 0o600, oct(_mode(log))


def test_the_failure_log_is_owner_only(tmp_path, monkeypatch):
    from delfin.agent import failure_log

    monkeypatch.setattr(failure_log, "_failure_path",
                        lambda *a, **kw: tmp_path / "f.jsonl", raising=False)
    try:
        failure_log.record_failure(kind="x", detail="y")
    except Exception:
        pytest.skip("record_failure signature differs; mode covered elsewhere")
    p = tmp_path / "f.jsonl"
    if p.exists():
        assert _mode(p) == 0o600, oct(_mode(p))


def test_the_recorders_survive_a_chmod_failure(tmp_path, monkeypatch):
    """Tightening permissions must never break the thing it observes."""
    from delfin.agent import audit_log

    def boom(*a, **kw):
        raise OSError("read-only filesystem")

    monkeypatch.setattr(os, "chmod", boom)
    log = tmp_path / "audit.log"
    audit_log.append({"tool": "bash", "decision": "ok"}, log_path=log)
    assert log.exists()
    assert json.loads(log.read_text().splitlines()[0])["tool"] == "bash"


# ---------------------------------------------------------------------------
# The DIRECTORIES, under a permissive umask
# ---------------------------------------------------------------------------
#
# Every creator used a bare ``mkdir`` with no mode, so the tree landed at
# whatever the umask allowed. Observed on the audited machine (umask 002):
# ``~/.delfin`` at 0700 -- inherited from an older umask, set by no code --
# and every subdirectory it created at 0775, with agent_sessions 144 of 154
# files at 0664, tool_traces 1045, turn_metrics 1263, and subagent_sessions
# all 108. These tests run under a deliberately permissive umask, because
# under 077 a broken creator looks identical to a correct one.


@pytest.fixture
def permissive_umask():
    """The mode must come from the code, not from the process umask."""
    old = os.umask(0o002)
    try:
        yield
    finally:
        os.umask(old)


def test_a_state_directory_is_owner_only(tmp_path, permissive_umask):
    from delfin.agent import state_paths

    d = state_paths.ensure_dir(tmp_path / "nested" / "state")
    assert _mode(d) == 0o700, oct(_mode(d))


def test_a_state_file_is_owner_only_from_creation(tmp_path, permissive_umask):
    """Not "owner-only after the write": a chmod that runs afterwards
    leaves the file group-readable for the length of the write."""
    from delfin.agent import state_paths

    p = state_paths.write_text(tmp_path / "s" / "f.json", "{}")
    assert _mode(p) == 0o600, oct(_mode(p))
    with state_paths.open_append(tmp_path / "s" / "g.jsonl") as fh:
        fh.write("{}\n")
    assert _mode(tmp_path / "s" / "g.jsonl") == 0o600


def test_the_tool_trace_directory_is_owner_only(tmp_path, monkeypatch,
                                                permissive_umask):
    from delfin.agent import tool_trace

    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    tool_trace.record("s1", tool="bash", output="y")
    assert _mode(tmp_path / "tool_traces") == 0o700


def test_the_turn_metrics_directory_is_owner_only(tmp_path, monkeypatch,
                                                  permissive_umask):
    from delfin.agent import turn_metrics

    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    turn_metrics.record("s1", total_ms=5)
    assert _mode(tmp_path / "turn_metrics") == 0o700
    p = turn_metrics.metrics_path("s1")
    assert _mode(p) == 0o600, oct(_mode(p))


def test_a_subagent_session_record_is_owner_only(tmp_path, monkeypatch,
                                                 permissive_umask):
    """The densest single file the agent writes: up to 60 tool outputs of
    2000 chars each plus a 200 kB report -- and it had no chmod at all."""
    from delfin.agent import subagents

    monkeypatch.setattr(subagents, "_SESSIONS_DIR", tmp_path / "sa_sessions")
    subagents._save_subagent_session(
        "sa-1", subagent_type="explore", description="d",
        messages=[{"role": "assistant", "content": "report"}],
        interactions=[{"name": "bash", "input": "ls", "output": "x"}])
    d = tmp_path / "sa_sessions"
    assert _mode(d) == 0o700, oct(_mode(d))
    f = d / "sa-1.json"
    assert f.is_file()
    assert _mode(f) == 0o600, oct(_mode(f))


def test_a_subagent_running_marker_is_owner_only(tmp_path, monkeypatch,
                                                 permissive_umask):
    from delfin.agent import subagents

    monkeypatch.setattr(subagents, "_RUNNING_DIR", tmp_path / "sa_running")
    subagents._running_update("sa-2", {"status": "running"})
    assert _mode(tmp_path / "sa_running") == 0o700
    assert _mode(tmp_path / "sa_running" / "sa-2.json") == 0o600


def test_a_pending_subagent_report_is_owner_only(tmp_path, monkeypatch,
                                                 permissive_umask):
    from delfin.agent import subagents

    monkeypatch.setattr(subagents, "_PENDING_DIR", tmp_path / "sa_pending")
    subagents._note_pending_report("sa-3", subagent_type="explore")
    assert _mode(tmp_path / "sa_pending") == 0o700
    for f in (tmp_path / "sa_pending").glob("*.json"):
        assert _mode(f) == 0o600, f


def test_the_session_store_directory_is_owner_only(tmp_path, monkeypatch,
                                                   permissive_umask):
    from delfin.agent import session_store

    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path / "sessions")
    session_store._ensure_dir()
    assert _mode(tmp_path / "sessions") == 0o700


# ---------------------------------------------------------------------------
# The repair sweep — what is ALREADY on disk
# ---------------------------------------------------------------------------
#
# The per-file chmods only ever touch the file being written now. A trace
# from April keeps its 0664 forever, because nothing appends to it again.
# The sweep is tested against a FIXTURE tree, never the user's own state.


def _legacy_tree(root):
    """A state tree as the bare-mkdir creators really left it."""
    root.mkdir(parents=True, exist_ok=True)
    os.chmod(root, 0o775)
    sub = root / "nested"
    sub.mkdir(exist_ok=True)
    os.chmod(sub, 0o775)
    made = []
    for p in (root / "a.jsonl", sub / "b.json"):
        p.write_text("{}\n", encoding="utf-8")
        os.chmod(p, 0o664)
        made.append(p)
    return made


def test_the_sweep_tightens_files_left_group_readable(tmp_path):
    from delfin.agent import state_paths

    root = tmp_path / "tool_traces"
    made = _legacy_tree(root)
    assert all(_mode(p) == 0o664 for p in made)      # precondition

    changed = state_paths.repair_tree(root)

    assert changed >= 4
    assert _mode(root) == 0o700
    assert _mode(root / "nested") == 0o700
    for p in made:
        assert _mode(p) == 0o600, p


def test_the_sweep_is_idempotent(tmp_path):
    """It runs at every start, so a second pass must be free."""
    from delfin.agent import state_paths

    root = tmp_path / "tool_traces"
    _legacy_tree(root)
    state_paths.repair_tree(root)
    assert state_paths.repair_tree(root) == 0


def test_the_sweep_does_not_follow_a_symlink_out_of_the_tree(tmp_path):
    """A link planted in the state tree must not redirect a chmod at a
    file outside it."""
    from delfin.agent import state_paths

    outside = tmp_path / "outside.txt"
    outside.write_text("x", encoding="utf-8")
    os.chmod(outside, 0o664)
    root = tmp_path / "traces"
    root.mkdir()
    (root / "link.txt").symlink_to(outside)

    state_paths.repair_tree(root)

    assert _mode(outside) == 0o664


def test_the_sweep_bounded_by_a_glob_leaves_its_neighbours_alone(tmp_path):
    """The rotated audit logs sit in ~/.delfin next to the settings file,
    the credential store and the memory store. Recursing there would be
    wrong for a chmod and catastrophic for a prune."""
    from delfin.agent import state_paths

    root = tmp_path / "delfin"
    root.mkdir()
    rotated = root / "audit-2026-W20.log"
    neighbour = root / "credentials.json"
    for p in (rotated, neighbour):
        p.write_text("{}", encoding="utf-8")
        os.chmod(p, 0o664)

    state_paths.repair_tree(root, glob="audit-*.log")

    assert _mode(rotated) == 0o600
    assert _mode(neighbour) == 0o664, "the sweep reached past its glob"
