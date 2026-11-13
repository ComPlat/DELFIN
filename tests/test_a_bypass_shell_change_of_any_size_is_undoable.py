"""In bypass, every change a foreground shell command makes in the
workspace is undoable, whatever the file's size and however it was made.

The shell change journal holds pre-images up to 1 MiB per file and 64 MiB
in all; the trash rule (``unrestorable``) reads only the targets a command
names. A deletion inside ``python -c`` of a 3 MiB file passed both and was
gone. Contract (``overlay_workspace``): the command runs over an overlay
of the roots; afterwards its changes are applied to the real roots, every
original it replaced or deleted moved (renamed) into
``.delfin-trash/<UTC date>/cmd-<id>/``, and each change is journalled so
``undo_changes`` renames it back.

Where the overlay cannot be mounted (no unprivileged user namespaces, no
bubblewrap, a filesystem without user xattrs) the tests assert the
fallback instead: the trash rule refuses what nothing could restore, and
a command is never run without its cage.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import change_journal as CJ
from delfin.agent import overlay_workspace as OW
from delfin.agent import unrestorable as U
from delfin.agent.api_client import KitToolPermissions, _doc_executor

MIB = 1024 * 1024


def _blob(n: int, seed: int) -> bytes:
    return bytes((i * 31 + seed) % 251 for i in range(256)) * (n // 256)


BIG = _blob(3 * MIB, 1)          # over the journal's 1 MiB per file
OTHER = _blob(3 * MIB, 2)


@pytest.fixture
def box(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    monkeypatch.setattr(Path, "home", lambda: home)
    ws = tmp_path / "ws"
    (ws / "tree" / "sub").mkdir(parents=True)
    (ws / "big.dat").write_bytes(BIG)
    (ws / "other.dat").write_bytes(OTHER)
    (ws / "tree" / "a.txt").write_text("a\n")
    (ws / "tree" / "sub" / "b.bin").write_bytes(_blob(2 * MIB, 3))
    perms = KitToolPermissions(mode="bypassPermissions", workspace=str(ws))
    perms.task_session_id = f"ov-{os.getpid()}-{tmp_path.name}"[:40]
    return ws.resolve(), perms


def _bash(perms, cmd: str) -> dict:
    return json.loads(_doc_executor.execute(
        "bash", {"command": cmd, "description": "d"}, perms))


def _read(perms, path: Path) -> None:
    """The read the shell's overwrite guard asks for before a named
    target is replaced; the trash rule is what the test is about."""
    perms.read_tracker[str(path.resolve())] = path.stat().st_mtime


def _undo(perms) -> dict:
    return json.loads(_doc_executor.execute(
        "undo_changes", {"scope": "session"}, perms))


def _overlaid(perms) -> bool:
    return bool(A._overlay_roots(perms)[0])


def _tree(ws: Path) -> dict:
    out = {}
    for p in sorted(ws.rglob("*")):
        if ".delfin-trash" in p.parts:
            continue
        out[str(p.relative_to(ws))] = (p.read_bytes() if p.is_file()
                                       else "dir")
    return out


def test_python_deletion_overwrite_and_rm_r_are_undone(box):
    ws, perms = box
    before = _tree(ws)
    if not _overlaid(perms):
        # The fallback: a named target nothing can restore is refused.
        _read(perms, ws / "big.dat")
        out = _bash(perms, "rm big.dat")
        assert "refused" in out.get("error", ""), out
        assert (ws / "big.dat").read_bytes() == BIG
        return
    for cmd in ("python3 -c \"import os; os.remove('big' + '.dat')\"",
                "python3 -c \"open('other' + '.dat', 'wb').write(b'new')\"",
                "rm -r tree"):
        out = _bash(perms, cmd)
        assert out.get("exit_code") == 0, out
    assert not (ws / "big.dat").exists()
    assert (ws / "other.dat").read_bytes() == b"new"
    assert not (ws / "tree").exists()
    recs = [r for r in CJ.list_changes(perms.task_session_id)
            if r.get("tool") == "bash"]
    assert {(Path(r["path"]).name, bool(r.get("deleted")), r["entry"])
            for r in recs} == {("big.dat", True, "file"),
                               ("other.dat", False, "file"),
                               ("tree", True, "dir")}, recs
    assert all(r.get("overlay") and not r.get("truncated") for r in recs)

    res = _undo(perms)
    assert not res.get("conflicts"), res
    assert _tree(ws) == before


def test_the_original_is_moved_not_copied(box):
    ws, perms = box
    if not _overlaid(perms):
        return
    ino = (ws / "big.dat").stat().st_ino
    out = _bash(perms, "python3 -c \"import os; os.remove('big' + '.dat')\"")
    assert out.get("exit_code") == 0, out
    rec = [r for r in CJ.list_changes(perms.task_session_id)
           if r.get("tool") == "bash"][-1]
    assert os.stat(rec["trash"]).st_ino == ino
    assert "/.delfin-trash/" in rec["trash"] and "/cmd-" in rec["trash"]


def test_an_undo_does_not_overwrite_the_persons_newer_file(box):
    ws, perms = box
    if not _overlaid(perms):
        return
    out = _bash(perms, "python3 -c \"open('other' + '.dat', 'wb').write(b'cmd')\"")
    assert out.get("exit_code") == 0, out
    (ws / "other.dat").write_bytes(b"person")
    res = _undo(perms)
    assert res.get("conflicts"), res
    assert (ws / "other.dat").read_bytes() == b"person"


def test_the_trash_is_read_only_and_the_layers_unseen(box):
    ws, perms = box
    if not _overlaid(perms):
        return
    _bash(perms, "python3 -c \"import os; os.remove('big' + '.dat')\"")
    out = _bash(perms, "rm -rf .delfin-trash/*/cmd-*; echo rc=$?")
    assert "rc=0" not in out.get("stdout", ""), out
    out = _bash(perms, "find . > /dev/null; echo find=$?; "
                       "ls -A .delfin-trash/.layers | wc -l")
    assert out.get("stdout", "").split() == ["find=0", "0"], out
    res = _undo(perms)
    assert (ws / "big.dat").read_bytes() == BIG, res


def test_a_git_ref_change_is_undone_and_objects_are_not_journalled(box):
    ws, perms = box
    if not _overlaid(perms):
        return
    env = {"PATH": "/usr/bin:/bin", "HOME": str(ws),
           "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.org",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.org"}
    for args in (["init", "-q", "-b", "main"], ["add", "tree"],
                 ["commit", "-q", "-m", "one"]):
        subprocess.run(["git", "-C", str(ws), *args], check=True,
                       capture_output=True, env=env)
    head = (ws / ".git" / "refs" / "heads" / "main").read_text()
    out = _bash(perms, "git -c user.name=t -c user.email=t@example.org "
                       "commit -q --allow-empty -m two && git log --oneline "
                       "| wc -l")
    assert out.get("stdout", "").strip() == "2", out
    names = [os.path.relpath(r["path"], ws)
             for r in CJ.list_changes(perms.task_session_id)
             if r.get("tool") == "bash"]
    assert ".git/refs/heads/main" in names, names
    assert not [n for n in names if n.startswith(".git/objects/")], names
    _undo(perms)
    assert (ws / ".git" / "refs" / "heads" / "main").read_text() == head


def _over(ws: Path, script: str):
    """Run ``script`` in ws over a real overlay; (plan, upper) before the
    apply, or None where no overlay mounts on this host."""
    roots, left = OW.eligible_roots([str(ws)])
    if not roots:
        assert left[str(ws)]
        return None
    plan = OW.prepare(roots)
    argv = ["bwrap", "--ro-bind", "/", "/", "--dev", "/dev", "--proc",
            "/proc", "--bind", str(ws), str(ws), "--chdir", str(ws),
            "/bin/sh", "-c", script]
    proc = subprocess.run(OW.wrap(argv, plan), capture_output=True,
                          text=True)
    assert proc.returncode == 0, proc.stderr
    assert not OW.setup_failed(plan)
    return plan, Path(plan.layers[0].upper)


def test_a_directory_removed_and_made_again_holds_only_the_new(tmp_path,
                                                               monkeypatch):
    """The opaque marker is set by this kernel in a user-namespace mount and
    read by the apply: nothing of the old directory survives in it."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ws = tmp_path / "ws"
    (ws / "d" / "sub").mkdir(parents=True)
    (ws / "d" / "old.txt").write_text("old\n")
    (ws / "d" / "sub" / "s.txt").write_text("s\n")
    ws = Path(os.path.realpath(ws))
    run = _over(ws, "rm -r d && mkdir d && echo x > d/new")
    if run is None:
        return
    plan, upper = run
    assert os.getxattr(upper / "d", "user.overlay.opaque") == b"y"
    out = OW.apply(plan)
    assert sorted(p.name for p in (ws / "d").iterdir()) == ["new"]
    trashed = {Path(c.trash).relative_to(plan.layers[0].trash).as_posix()
               for c in out.changes if c.trash}
    assert trashed == {"d/old.txt", "d/sub"}, out.changes
    CJ.record_overlay_changes("opq", [
        {"tool": "bash", "path": c.path, "kind": c.kind, "entry": c.entry,
         "root": c.root, "trash": c.trash, "post_hash": c.post_hash,
         "post_sig": c.post_sig, "post_tree": c.post_tree}
        for c in out.changes if c.journal])
    res = CJ.revert("opq", scope="session", workspace=ws)
    assert not res["conflicts"], res
    assert sorted(p.name for p in (ws / "d").iterdir()) == ["old.txt", "sub"]
    assert (ws / "d" / "sub" / "s.txt").read_text() == "s\n"


def test_a_directory_renamed_is_undone_as_one(tmp_path, monkeypatch):
    """No redirect_dir in a user namespace: mv of a directory copies it into
    the upper layer and whites the old name out; both are applied, and the
    undo puts the old name back."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ws = tmp_path / "ws"
    (ws / "a" / "x").mkdir(parents=True)
    (ws / "a" / "x" / "f").write_bytes(b"f" * 5000)
    (ws / "a" / "g").write_text("g\n")
    ws = Path(os.path.realpath(ws))
    run = _over(ws, "mv a b")
    if run is None:
        return
    plan, _upper = run
    out = OW.apply(plan)
    assert not (ws / "a").exists()
    assert (ws / "b" / "x" / "f").read_bytes() == b"f" * 5000
    assert {(Path(c.path).name, c.kind) for c in out.changes} == {
        ("a", "deleted"), ("b", "created")}
    CJ.record_overlay_changes("mv", [
        {"tool": "bash", "path": c.path, "kind": c.kind, "entry": c.entry,
         "root": c.root, "trash": c.trash, "post_hash": c.post_hash,
         "post_sig": c.post_sig, "post_tree": c.post_tree}
        for c in out.changes])
    res = CJ.revert("mv", scope="session", workspace=ws)
    assert not res["conflicts"], res
    assert (ws / "a" / "x" / "f").read_bytes() == b"f" * 5000
    assert not (ws / "b").exists()


def test_the_trash_over_its_budget_loses_its_oldest_commands(tmp_path,
                                                             monkeypatch):
    ws = tmp_path / "ws"
    ws.mkdir()
    ws = Path(os.path.realpath(ws))
    undo = ws / ".delfin-trash" / ".undo"
    for day, cmd, size in (("2026-11-01", "cmd-090000-1-1", 3000),
                           ("2026-11-02", "cmd-090000-1-2", 3000),
                           ("2026-11-03", "cmd-090000-1-3", 3000)):
        (undo / day / cmd).mkdir(parents=True)
        (undo / day / cmd / "f").write_bytes(b"x" * size)
    events = []
    monkeypatch.setattr(OW, "_budget_bytes", lambda root, used: used - 1)
    monkeypatch.setattr(OW, "_record_event",
                        lambda kind, detail: events.append((kind, detail)))
    OW._USAGE.pop(str(ws), None)
    notes = OW.keep_within_budget(str(ws), keep=str(undo / "2026-11-03"))
    left = sorted(p.parent.name for p in undo.glob("*/cmd-*"))
    assert left == ["2026-11-02", "2026-11-03"]
    assert events and events[0][0] == "trash_evicted"
    assert notes and "removed" in notes[0]


def test_a_change_made_while_the_command_ran_is_kept_beside(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "f.txt").write_text("old\n")
    ws = Path(os.path.realpath(ws))
    roots, left = OW.eligible_roots([str(ws)])
    if not roots:
        assert left[str(ws)]
        return
    plan = OW.prepare(roots)
    argv = ["bwrap", "--ro-bind", "/", "/", "--dev", "/dev", "--proc",
            "/proc", "--bind", str(ws), str(ws), "--chdir", str(ws),
            "/bin/sh", "-c", "echo cmd > f.txt"]
    proc = subprocess.run(OW.wrap(argv, plan), capture_output=True,
                          text=True)
    assert proc.returncode == 0, proc.stderr
    (ws / "f.txt").write_text("person\n")     # after the start
    out = OW.apply(plan)
    assert (ws / "f.txt").read_text() == "person\n"
    side = ws / f"f.txt.delfin-conflict-{plan.cmd_id}"
    assert side.read_text() == "cmd\n"
    assert any("changed by someone else" in n for n in out.notes)


def test_a_mount_failure_runs_the_old_way_with_the_trash_rule(box, monkeypatch):
    ws, perms = box
    if not _overlaid(perms):
        return
    fail = ["/bin/sh", "-c",
            "echo 'delfin-overlay: simulated' >&2; exit 97", "--"]
    monkeypatch.setattr(OW, "wrap", lambda argv, plan: fail + list(argv))
    monkeypatch.setattr(OW, "_RUNTIME_FAILURE", "")
    _read(perms, ws / "big.dat")
    out = _bash(perms, "rm big.dat")
    assert "refused" in out.get("error", ""), out
    assert (ws / "big.dat").read_bytes() == BIG
    assert OW._RUNTIME_FAILURE == "simulated"
    # From now on the gate itself refuses, and nothing is overlaid.
    assert not _overlaid(perms)


def test_a_command_imitating_the_failure_does_not_switch_it_off(box):
    ws, perms = box
    if not _overlaid(perms):
        return
    out = _bash(perms, "echo run >> runs.log; "
                       "echo 'delfin-overlay: fake' >&2; exit 97")
    assert out.get("exit_code") == 97, out
    assert (ws / "runs.log").read_text() == "run\n"      # ran once
    assert _overlaid(perms)


def test_the_gate_lets_a_covered_target_pass(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "big.dat").write_bytes(b"x" * (2 * MIB))
    perms = KitToolPermissions(mode="bypassPermissions", workspace=str(ws))
    assert U.refusal("rm big.dat", ws, perms)
    assert U.refusal("rm big.dat", ws, perms, covered=lambda p: True) is None
    (ws / ".delfin-trash").mkdir(exist_ok=True)
    (ws / ".delfin-trash" / "old").write_bytes(b"x" * (2 * MIB))
    assert U.refusal("rm .delfin-trash/old", ws, perms,
                     covered=lambda p: True)


def test_overlay_records_do_not_count_toward_the_pre_image_cap(tmp_path,
                                                               monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.setattr(CJ, "MAX_ENTRIES_PER_SESSION", 3)
    sid = "cap"
    for i in range(3):
        CJ.record_deletion(sid, tool="bash", path=str(tmp_path / f"p{i}"),
                           pre_bytes=b"x")
    CJ.record_overlay_changes(sid, [
        {"tool": "bash", "path": str(tmp_path / f"o{i}"), "kind": "created",
         "entry": "file"} for i in range(10)])
    recs = CJ.list_changes(sid)
    assert not [r for r in recs if r.get("dropped")], recs
    assert len(recs) == 13
