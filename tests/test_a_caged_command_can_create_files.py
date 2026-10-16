"""A command admitted unasked by the cage creates, removes and renames files
in its workspace roots, and the read-only overlay still holds.

Input: a command ``_caged_command_may_run_unasked`` admits, driven through
``_DocToolExecutor().execute("bash", ...)`` in ``default`` mode on a
workspace that is a git repository and holds protected files at the top
level (``.github``, ``.gitattributes``) and deeper (``.delfin/settings.json``,
``delfin/agent/api_client.py``).

Output measured on the host after the command:

* a new file and directory at the top level of the workspace, inside a
  directory on the way to a protected file (``.delfin/``,
  ``delfin/agent/``), a removed and a renamed top-level file -- all done;
* every ``_caged_readonly_paths`` entry byte-identical, nothing planted
  inside a protected directory (``.git/hooks/x`` among them), and no
  directory on the way to a protected path renamed away (renaming it and
  recreating the path would leave a fresh, writable copy of that path).

The same for a granted ``~/.delfin`` holding the memory stores.

The protected file names are assembled at run time inside a script the
gate never reads: a literal protected path is refused by the gate before
the cage, which is not what is measured here.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest

from conftest import requires_a_working_sandbox
from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def landlock_host():
    if not A._landlock_functional():
        pytest.skip("no Landlock ABI on this kernel: the caged rule does "
                    "not admit here")


@pytest.fixture
def repo(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    ws = tmp_path / "ws"
    ws.mkdir()
    subprocess.run(["git", "init", "-q", str(ws)], check=True,
                   capture_output=True)
    (ws / ".git" / "hooks").mkdir(exist_ok=True)
    (ws / ".delfin").mkdir()
    (ws / ".delfin" / "settings.json").write_text("{}")
    (ws / ".github").mkdir()
    (ws / ".github" / "ci.yml").write_text("on: push\n")
    (ws / ".gitattributes").write_text("* text=auto\n")
    (ws / "delfin" / "agent").mkdir(parents=True)
    (ws / "delfin" / "agent" / "api_client.py").write_text("# agent\n")
    (ws / "data.txt").write_text("a b c\n")
    (ws / "old.txt").write_text("old\n")
    return ws


def _perms(ws: Path, **kw) -> KitToolPermissions:
    return KitToolPermissions(mode="default", workspace=str(ws), **kw)


def _run(perms: KitToolPermissions, script: Path, *args: str) -> dict:
    """Run ``script`` through a shell script: ``bash x.sh`` is admitted by
    the caged rule, while ``python3 x.py`` is on the auto-allow list and
    runs in the ordinary cage without the read-only overlay."""
    sh = script.with_suffix(".sh")
    sh.write_text(" ".join(["python3", script.name, *args]) + "\n")
    seen: list[bool] = []
    real = A._bash_isolation_argv

    def _spy(cmd, run_cwd, perms, mode=None, extra_write=(),
             caged_unasked=False):
        seen.append(caged_unasked)
        return real(cmd, run_cwd, perms, mode, extra_write, caged_unasked)

    A._bash_isolation_argv = _spy
    try:
        out = json.loads(_doc_executor.execute(
            "bash", {"command": f"bash {sh.name}"}, perms))
    finally:
        A._bash_isolation_argv = real
    assert "not on the auto-allow list" not in out.get("error", ""), out
    assert seen == [True], (seen, out)
    return out


def _snapshot(paths) -> dict:
    """Every file and directory beneath ``paths``, with file digests."""
    snap: dict[str, str] = {}
    for p in paths:
        p = Path(p)
        items = [p] + (sorted(p.rglob("*")) if p.is_dir() else [])
        for q in items:
            snap[str(q)] = ("dir" if q.is_dir() else
                            hashlib.sha256(q.read_bytes()).hexdigest())
    return snap


_ATTACK = r'''
import json, os, sys
targets, roots = json.loads(open(sys.argv[1]).read())
done = []
def attempt(name, fn):
    try:
        fn()
        done.append(name)
    except OSError:
        pass
def ancestors(p):
    out = []
    root = max((r for r in roots if p.startswith(r + "/")), key=len)
    d = os.path.dirname(p)
    while d != root and d.startswith(root + "/"):
        out.append(d)
        d = os.path.dirname(d)
    return out
for t in targets:
    if os.path.isdir(t):
        attempt("plant " + t, lambda: open(os.path.join(t, "planted"), "w").write("x"))
        attempt("mkdir " + t, lambda: os.mkdir(os.path.join(t, "planted.d")))
        for base, dirs, files in os.walk(t):
            for f in files:
                fp = os.path.join(base, f)
                attempt("append " + fp, lambda: open(fp, "a").write("x"))
    else:
        attempt("append " + t, lambda: open(t, "a").write("x"))
        attempt("unlink " + t, lambda: os.unlink(t))
    attempt("rename " + t, lambda: os.rename(t, t + ".moved"))
    for a in ancestors(t):
        attempt("rename " + a, lambda: os.rename(a, a + ".moved"))
print("DONE " + json.dumps(done))
'''


def _attack(perms: KitToolPermissions, ws: Path, targets, roots) -> list:
    (ws / "targets.json").write_text(json.dumps(
        [list(targets), [str(Path(r).resolve()) for r in roots]]))
    (ws / "attack.py").write_text(_ATTACK)
    out = _run(perms, ws / "attack.py", "targets.json")
    line = [ln for ln in out.get("stdout", "").splitlines()
            if ln.startswith("DONE ")]
    assert line, out
    return json.loads(line[0][5:])


@requires_a_working_sandbox
def test_a_caged_command_creates_files_in_a_git_workspace(repo,
                                                          landlock_host):
    (repo / "work.py").write_text(
        "import os\n"
        "open('new.txt', 'w').write('x')\n"
        "os.mkdir('newdir')\n"
        "open('newdir/f.txt', 'w').write('y')\n"
        "open('.delfin/notes.txt', 'w').write('n')\n"
        "open('delfin/agent/helper.py', 'w').write('h')\n"
        "os.remove('data.txt')\n"
        "os.rename('old.txt', 'renamed.txt')\n"
        "print('WORK-OK')\n")
    out = _run(_perms(repo), repo / "work.py")
    assert "WORK-OK" in out.get("stdout", ""), out
    assert (repo / "new.txt").read_text() == "x"
    assert (repo / "newdir" / "f.txt").read_text() == "y"
    assert (repo / ".delfin" / "notes.txt").read_text() == "n"
    assert (repo / "delfin" / "agent" / "helper.py").read_text() == "h"
    assert not (repo / "data.txt").exists()
    assert (repo / "renamed.txt").read_text() == "old\n"


@requires_a_working_sandbox
def test_every_read_only_path_still_refuses_the_caged_command(repo,
                                                              landlock_host):
    perms = _perms(repo)
    targets = A._caged_readonly_paths(perms)
    real = {str(Path(p).resolve()) for p in targets}
    for rel in (".git", ".github", ".gitattributes",
                ".delfin/settings.json", "delfin/agent/api_client.py"):
        assert str((repo / rel).resolve()) in real, (rel, targets)
    before = _snapshot(targets)
    done = _attack(perms, repo, targets, perms.all_workspace_roots())
    assert done == []
    assert _snapshot(targets) == before
    assert not list(repo.rglob("*.moved"))
    assert not (repo / ".git" / "hooks" / "planted").exists()


@requires_a_working_sandbox
def test_the_home_memory_stores_hold_and_their_parents_stay_usable(
        repo, landlock_host):
    home = Path(os.environ["HOME"])
    base = home / ".delfin"
    (base / "memory").mkdir(parents=True)
    (base / "memory" / "MEMORY.md").write_text("- a\n")
    (base / "projects" / "p1" / "memory").mkdir(parents=True)
    (base / "projects" / "p1" / "memory" / "n.md").write_text("n\n")
    perms = _perms(repo, extra_workspace_dirs=(base.resolve(),))
    targets = A._caged_readonly_paths(perms)
    real = {str(Path(p).resolve()) for p in targets}
    assert str((base / "memory").resolve()) in real
    assert str((base / "projects" / "p1" / "memory").resolve()) in real
    stores = [p for p in targets if p.startswith(str(base.resolve()))]
    before = _snapshot(stores)
    done = _attack(perms, repo, stores, perms.all_workspace_roots())
    assert done == []
    assert _snapshot(stores) == before
    # The directory on the way to a project store takes new files.
    (repo / "notes.py").write_text(
        f"open({str(base / 'projects' / 'p1' / 'notes.txt')!r}, 'w')"
        ".write('k')\nprint('NOTES-OK')\n")
    out = _run(perms, repo / "notes.py")
    assert "NOTES-OK" in out.get("stdout", ""), out
    assert (base / "projects" / "p1" / "notes.txt").read_text() == "k"


def test_the_pins_are_the_directories_on_the_way(tmp_path):
    ws = tmp_path / "ws"
    (ws / "a" / "b").mkdir(parents=True)
    (ws / "a" / "b" / "f").write_text("x")
    (ws / "top").write_text("x")
    (ws / ".git").mkdir()
    real = tmp_path / "elsewhere"
    real.mkdir()
    (real / "g").write_text("x")
    (ws / "link").symlink_to(real, target_is_directory=True)
    pins, unpinned = A._mount_pins(
        [str(ws / "a" / "b" / "f"), str(ws / "top"), str(ws / ".git"),
         str(tmp_path / "outside"), str(ws / "link" / "g")],
        [str(ws)])
    assert pins == [str(ws / "a"), str(ws / "a" / "b")]
    assert unpinned == []          # outside and the link's target: no root
    # A chain through something that is not a plain directory is not
    # pinned; its caller keeps the Landlock carve-out for it.
    pins, unpinned = A._mount_pins([str(ws / "gone" / "f")], [str(ws)])
    assert pins == [] and unpinned == [str(ws / "gone" / "f")]


def test_an_unpinnable_path_keeps_the_landlock_carve_out(repo, monkeypatch):
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    monkeypatch.setattr(A, "_BWRAP_RUNTIME_FAILURE", None)
    gone = str(repo / "vanished" / "settings.json")
    real_ro = A._caged_readonly_paths
    monkeypatch.setattr(A, "_caged_readonly_paths",
                        lambda perms: real_ro(perms) + [gone])
    argv = A._bash_isolation_argv("true", repo, _perms(repo),
                                  caged_unasked=True)
    ro = [argv[i + 1] for i, a in enumerate(argv) if a == "--read-only"]
    assert ro == [gone]
    pinned = [argv[i + 1] for i, a in enumerate(argv)
              if a in ("--bind-try", "--bind")]
    for d in (".delfin", "delfin", "delfin/agent"):
        assert str((repo / d).resolve()) in pinned
    # Shallow first, the git metadata binds included: bwrap binds from the
    # host view and recursively, so a later ancestor bind would hide a pin.
    depth = [p.rstrip("/").count("/") for p in pinned]
    assert depth == sorted(depth)


@requires_a_working_sandbox
def test_a_refused_write_says_it_is_the_cage(repo, landlock_host):
    (repo / "touch.py").write_text(
        "p = '.githu' + 'b/ci.yml'\nopen(p, 'a').write('x')\n")
    out = _run(_perms(repo), repo / "touch.py")
    assert out.get("exit_code") != 0
    assert "Read-only file system" in out.get("stderr", "")
    note = out.get("cage_note", "")
    assert note.startswith("cage: ") and ".github" in note, out
    # A command that fails for its own reasons gets no such line.
    (repo / "fail.py").write_text("raise SystemExit(3)\n")
    out = _run(_perms(repo), repo / "fail.py")
    assert out.get("exit_code") == 3 and "cage_note" not in out
