"""A file written by a shell command is in the change journal.

A command the bash gate admits without a question because the cage holds
it (``_caged_command_may_run_unasked``) may write workspace files through
an interpreter: ``python3 -c "open('x','w').write(...)"``, a script, a
compiler. Those writes do not pass through write_file, so before this
they were not journalled: ``undo_changes`` could not take them back, the
completion check (``_journal_changes``) did not count them as written,
and the self-written detection (``_own_written_file``) did not know the
session had written the script it was refusing.

Contract (``delfin.agent.shell_changes``): every foreground ``bash`` and
every ``bash_background`` command is bracketed by a (path, size,
mtime_ns, inode) snapshot of the writable roots. Each file created,
modified or deleted under them becomes one journal entry -- created
files are removed by undo, modified and deleted ones restored from the
pre-image held for the snapshot. A pre-image over the per-file cap is
recorded as a change without undo, and the entry says so. A root over
the file-count cap is not tracked, and the bash result says so once.

The caged decision is pinned and the argv replaced by a plain shell, so
the journal path runs on every host; one test drives the real cage
where it works.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pytest

from conftest import requires_a_working_sandbox
from delfin.agent import api_client as A
from delfin.agent import change_journal as CJ
from delfin.agent.api_client import KitToolPermissions, _doc_executor

# `rm -rf /tmp/...` as hex: a script carrying it trips the payload scan,
# and a command text carrying it would be refused before it ran.
_TRIPS_HEX = "import os\nos.system('rm -rf /tmp/never-there-cj')\n".encode(
).hex()


# The path is assembled at run time, so no reading of the command text
# finds the target: only the snapshot of the roots can.
_PY_WRITE = ("python3 -c \"import os; open(os.path.join('.', 'notes.txt'), "
             "'w').write('{text}' + chr(10))\"")


@pytest.fixture
def box(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", lambda: home)
    ws = tmp_path / "ws"
    ws.mkdir()
    # The gate admits the command because the cage holds it ...
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: True)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    # ... and the executor runs it in a plain shell, so the journal path
    # is exercised on hosts where bubblewrap cannot run.
    seen: list[bool] = []

    def _argv(cmd, run_cwd, perms, mode=None, extra_write=(),
              caged_unasked=False):
        seen.append(caged_unasked)
        return ["/bin/bash", "-c", cmd]

    monkeypatch.setattr(A, "_bash_isolation_argv", _argv)
    perms = KitToolPermissions(mode="acceptEdits", workspace=str(ws))
    perms.task_session_id = f"cj-{tmp_path.name}"[:40]
    return ws, perms, seen


def _bash(perms, cmd: str, tool: str = "bash") -> dict:
    return json.loads(_doc_executor.execute(
        tool, {"command": cmd, "description": "d"}, perms))


def _write(perms, ws: Path, rel: str, content: str) -> None:
    out = str(_doc_executor.execute(
        "write_file", {"path": str(ws / rel), "content": content}, perms))
    assert '"error"' not in out, out
    assert (ws / rel).read_text() == content


def _records(perms, tool=("bash", "bash_background")) -> list[dict]:
    return [r for r in CJ.list_changes(perms.task_session_id)
            if r.get("tool") in tool]


def _undo(perms) -> dict:
    return json.loads(_doc_executor.execute(
        "undo_changes", {"scope": "session"}, perms))


def test_a_caged_python_write_is_journalled_and_undone(box):
    ws, perms, seen = box
    (ws / "notes.txt").write_text("original\n")
    out = _bash(perms, _PY_WRITE.format(text="replaced"))
    assert "error" not in out, out
    assert seen and seen[-1] is True, "the command was not admitted caged"
    assert (ws / "notes.txt").read_text() == "replaced\n"
    recs = _records(perms)
    assert [Path(r["path"]).name for r in recs] == ["notes.txt"], recs

    res = _undo(perms)
    assert str(ws / "notes.txt") in res.get("reverted", []), res
    assert (ws / "notes.txt").read_text() == "original\n"


def test_a_file_a_script_created_is_removed_by_undo(box):
    ws, perms, _ = box
    _write(perms, ws, "gen.py",
           "open('results.csv','w').write('a,b\\n1,2\\n')\n")
    out = _bash(perms, "python3 gen.py")
    assert "error" not in out, out
    assert (ws / "results.csv").is_file()
    recs = _records(perms)
    assert [(Path(r["path"]).name, r["created"]) for r in recs] == [
        ("results.csv", True)], recs

    _undo(perms)
    assert not (ws / "results.csv").exists()


def test_a_file_a_script_deleted_is_restored_by_undo(box):
    ws, perms, _ = box
    (ws / "old.dat").write_bytes(b"\x00\x01keep me\xff")
    _write(perms, ws, "drop.py", "import os\nos.remove('old.dat')\n")
    out = _bash(perms, "python3 drop.py")
    assert "error" not in out, out
    assert not (ws / "old.dat").exists()
    recs = _records(perms)
    assert [(Path(r["path"]).name, bool(r.get("deleted"))) for r in recs] \
        == [("old.dat", True)], recs

    res = _undo(perms)
    assert str(ws / "old.dat") in res.get("reverted", []), res
    assert (ws / "old.dat").read_bytes() == b"\x00\x01keep me\xff"


def test_self_written_recognises_a_script_written_through_the_shell(box):
    ws, perms, _ = box
    _write(perms, ws, "gen.py",
           f"open('x.py','w').write(bytes.fromhex('{_TRIPS_HEX}')"
           ".decode())\n")
    out = _bash(perms, "python3 gen.py")
    assert "error" not in out, out
    assert (ws / "x.py").is_file()
    refused = _doc_executor.execute(
        "bash", {"command": "python3 x.py", "description": "d"}, perms)
    assert "trip the safety scan" in refused, refused
    assert "a file this session wrote" in refused, refused


def test_the_completion_check_sees_the_written_file(box):
    ws, perms, _ = box
    _write(perms, ws, "gen.py",
           "open('results.csv','w').write('a,b\\n1,2\\n')\n")
    before = A._journal_changes(perms.task_session_id)
    assert not [c for c in before if c["path"].endswith("results.csv")]
    control = A.check_completion_claim("Write results.csv", changes=before)
    assert control["verdict"] != "verified", control

    _bash(perms, "python3 gen.py")
    after = A._journal_changes(perms.task_session_id)
    assert [c for c in after if c["path"].endswith("results.csv")
            and c["created"]], after
    verdict = A.check_completion_claim("Write results.csv", changes=after)
    assert verdict["verdict"] == "verified", verdict


def test_an_unchanged_run_records_nothing(box):
    ws, perms, _ = box
    (ws / "data.txt").write_text("a b c\n")
    out = _bash(perms, 'python3 -c "print(open(\'data.txt\').read())"')
    assert "a b c" in out.get("stdout", ""), out
    assert _records(perms) == []


def test_git_internals_and_bytecode_are_not_journalled(box):
    ws, perms, _ = box
    (ws / ".git").mkdir()
    (ws / "mod.py").write_text("X = 1\n")
    _write(perms, ws, "w.py",
           "import os, mod\n"
           "open(os.path.join('.git', 'probe'), 'w').write('x')\n"
           "open('kept.txt', 'w').write('y')\n")
    _bash(perms, "python3 w.py")
    assert (ws / "__pycache__").is_dir()
    names = [Path(r["path"]).name for r in _records(perms)]
    assert names == ["kept.txt"], names


def test_an_oversize_pre_image_is_recorded_without_undo(box, monkeypatch):
    from delfin.agent import shell_changes as SC
    monkeypatch.setattr(SC, "MAX_PRE_IMAGE_BYTES", 64)
    ws, perms, _ = box
    (ws / "big.txt").write_text("x" * 200)
    _write(perms, ws, "grow.py", "open('big.txt','a').write('more')\n")
    _bash(perms, "python3 grow.py")
    recs = _records(perms)
    assert [Path(r["path"]).name for r in recs] == ["big.txt"], recs
    assert recs[0].get("truncated") and "too large" in recs[0].get(
        "reason", ""), recs
    res = _undo(perms)
    skipped = [s for s in res.get("skipped", [])
               if s.get("path") == str(ws / "big.txt")]
    assert skipped and "too large" in skipped[0]["reason"], res
    assert (ws / "big.txt").read_text().endswith("more")


def test_a_huge_tree_is_capped_and_said_once(box, monkeypatch):
    from delfin.agent import shell_changes as SC
    monkeypatch.setattr(SC, "MAX_FILES_PER_ROOT", 50)
    ws, perms, _ = box
    (ws / "many").mkdir()
    for i in range(60):
        (ws / "many" / f"f{i}.txt").write_text(str(i))
    _write(perms, ws, "gen.py", "open('out.txt','w').write('z')\n")
    first = _bash(perms, "python3 gen.py")
    assert "too many files to track" in json.dumps(first), first
    second = _bash(perms, "python3 gen.py")
    # The reason once; afterwards a short reminder, so the model is not
    # left reading the second result as journalled.
    assert "too many files to track" not in json.dumps(second), second
    assert "still not journalled" in second.get("untracked_changes", "")
    assert _records(perms) == []


def _counting_slow_scandir(monkeypatch, SC, delay: float) -> list[str]:
    """A slow filesystem stand-in: every directory listing costs
    ``delay`` s. Returns the list of directories listed."""
    listed: list[str] = []
    real = os.scandir

    def _slow(path="."):
        listed.append(str(path))
        time.sleep(delay)
        return real(path)

    monkeypatch.setattr(SC.os, "scandir", _slow)
    return listed


def test_a_slow_root_is_not_walked_again_on_the_next_command(box,
                                                            monkeypatch):
    from delfin.agent import shell_changes as SC
    monkeypatch.setattr(SC, "MAX_WALK_SECONDS", 0.2)
    ws, perms, _ = box
    for d in range(12):
        (ws / f"d{d}").mkdir()
        (ws / f"d{d}" / "x.txt").write_text("x")
    _write(perms, ws, "gen.py", "open('out.txt','w').write('z')\n")
    listed = _counting_slow_scandir(monkeypatch, SC, 0.05)

    t0 = time.monotonic()
    first = _bash(perms, "python3 gen.py")
    first_s = time.monotonic() - t0
    assert "too many files to track" in json.dumps(first), first
    # The walk stops at its budget instead of finishing: 13 listings
    # at 0.05 s would be 0.65 s per walk.
    assert len(listed) < 13, listed
    assert first_s < 5, first_s

    listed.clear()
    second = _bash(perms, "python3 gen.py")
    assert listed == [], "the remembered root was walked again"
    assert "still not journalled" in second.get("untracked_changes", ""), \
        second
    assert _records(perms) == []


def test_a_remembered_root_is_walked_again_after_the_recheck_time(
        box, monkeypatch):
    from delfin.agent import shell_changes as SC
    monkeypatch.setattr(SC, "MAX_FILES_PER_ROOT", 5)
    ws, perms, _ = box
    for i in range(8):
        (ws / f"f{i}.txt").write_text(str(i))
    _write(perms, ws, "gen.py", "open('out.txt','w').write('z')\n")
    _bash(perms, "python3 gen.py")
    assert _records(perms) == []
    # Below the bound again, and the memory has expired.
    for i in range(8):
        (ws / f"f{i}.txt").unlink()
    monkeypatch.setattr(SC, "RECHECK_SECONDS", 0.0)
    (ws / "out.txt").unlink()
    _bash(perms, "python3 gen.py")
    assert [Path(r["path"]).name for r in _records(perms)] == ["out.txt"]


def test_calc_and_read_only_roots_are_never_walked(box, monkeypatch):
    from delfin.agent import shell_changes as SC
    ws, perms, _ = box
    calc = ws / "calc"
    calc.mkdir()
    (calc / "orca.out").write_text("big output\n")
    ro = ws.parent / "ro"
    ro.mkdir()
    perms.confirm_write_dirs = (calc.resolve(),)
    perms.read_only_workspace_dirs = (ro.resolve(),)
    perms.extra_workspace_dirs = (ro.resolve(),)
    listed = _counting_slow_scandir(monkeypatch, SC, 0.0)
    _write(perms, ws, "gen.py",
           "open('calc/new.out','w').write('n')\n"
           "open('kept.txt','w').write('k')\n")
    _bash(perms, "python3 gen.py")
    assert not [d for d in listed
                if str(calc.resolve()) in d or str(ro.resolve()) in d], listed
    assert [Path(r["path"]).name for r in _records(perms)] == ["kept.txt"]


def test_a_background_job_is_journalled_when_it_ends(box):
    ws, perms, seen = box
    _write(perms, ws, "gen.py",
           "open('bg.out','w').write('done')\n")
    out = _bash(perms, "python3 gen.py", tool="bash_background")
    assert out.get("status") == "started", out
    deadline = time.monotonic() + 30
    while time.monotonic() < deadline:
        if _records(perms):
            break
        time.sleep(0.1)
    recs = _records(perms)
    assert [(r["tool"], Path(r["path"]).name, r["created"]) for r in recs] \
        == [("bash_background", "bg.out", True)], recs
    _undo(perms)
    assert not (ws / "bg.out").exists()


@requires_a_working_sandbox
def test_through_the_real_cage(tmp_path, monkeypatch):
    """No argv pinned: the command runs under bubblewrap where it works."""
    if not A._landlock_functional():
        pytest.skip("no Landlock ABI on this kernel: the caged rule does "
                    "not admit here")
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", lambda: home)
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "notes.txt").write_text("original\n")
    perms = KitToolPermissions(mode="acceptEdits", workspace=str(ws))
    perms.task_session_id = f"cjreal-{os.getpid()}"
    out = _bash(perms, _PY_WRITE.format(text="caged"))
    assert "error" not in out, out
    assert (ws / "notes.txt").read_text() == "caged\n"
    assert [Path(r["path"]).name for r in _records(perms)] == ["notes.txt"]
    _undo(perms)
    assert (ws / "notes.txt").read_text() == "original\n"
