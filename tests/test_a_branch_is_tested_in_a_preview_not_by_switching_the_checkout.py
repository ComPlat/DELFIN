"""A branch is tested in a second dashboard from its own worktree.

Input: a shared checkout on ``main`` with uncommitted work, a branch to
test. Output: a worktree at ``<parent>/delfin-preview/<slug>`` (detached at
the head when the branch is held elsewhere), a server started from it with
``PYTHONPATH=<worktree>`` and a scrubbed environment, its token read back
from the launcher's token file, and a registry record; closing it ends the
process group and removes the worktree when it is clean, and asks when it
is not. The shared checkout's branch, head, files and stash are unchanged
throughout.

Universal: the repository is made here and carries its own ``delfin``
package whose ``cli_voila`` is a stand-in server (a socket and one child
process in the same group), so no Voila, no network and no fixed port.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

from delfin.dashboard import branch_preview as BP

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}

_FAKE_LAUNCHER = textwrap.dedent('''
    import argparse, json, os, pathlib, socket, subprocess, sys
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int)
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--strict-port", action="store_true")
    a = ap.parse_args()
    out = pathlib.Path(os.environ["FAKE_OUT"])
    out.mkdir(parents=True, exist_ok=True)
    (out / "env.json").write_text(json.dumps(dict(os.environ)))
    tokdir = pathlib.Path(os.environ["XDG_RUNTIME_DIR"])
    tok = tokdir / f"voila-token-{os.getpid()}-abc"
    fd = os.open(str(tok), os.O_WRONLY | os.O_CREAT, 0o600)
    os.write(fd, b"secret-token-123"); os.close(fd)
    print(f"Read it with:  cat {tok}", flush=True)
    # an old launcher overwrites the note that says where the dashboard runs
    note = pathlib.Path.home() / ".delfin" / "dashboard_here.json"
    note.parent.mkdir(parents=True, exist_ok=True)
    note.write_text(json.dumps({"pid": os.getpid(), "from": "preview"}))
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(600)"])
    (out / "child.pid").write_text(str(child.pid))
    s = socket.socket(); s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    s.bind(("127.0.0.1", a.port)); s.listen(5)
    while True:
        c, _ = s.accept()
        req = c.recv(4096).decode("latin-1")
        if "token=secret-token-123" in req.split("\\r\\n", 1)[0]:
            c.sendall(b"HTTP/1.1 302 Found\\r\\nLocation: /voila/render/"
                      b"delfin/notebooks/delfin_dashboard.ipynb\\r\\n"
                      b"Content-Length: 0\\r\\n\\r\\n")
        else:
            c.sendall(b"HTTP/1.1 403 Forbidden\\r\\nContent-Length: 0\\r\\n\\r\\n")
        c.close()
''')


def _git(cwd, *args, check=True) -> str:
    out = subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=str(cwd), env=_ENV, check=check, capture_output=True, text=True)
    return out.stdout.strip()


@pytest.fixture
def setup(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    runtime = tmp_path / "run"
    runtime.mkdir(mode=0o700)
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    root = tmp_path / "projects"
    root.mkdir()
    co = root / "DELFIN"
    _git(tmp_path, "clone", "-q", str(origin), str(co))
    _git(co, "checkout", "-q", "-b", "main")
    (co / "delfin").mkdir()
    (co / "delfin" / "__init__.py").write_text("")
    (co / "delfin" / "cli_voila.py").write_text(_FAKE_LAUNCHER)
    (co / ".gitignore").write_text("delfin_voila_runtime/\n__pycache__/\n")
    _git(co, "add", "-A")
    _git(co, "commit", "-q", "-m", "init")
    _git(co, "push", "-q", "-u", "origin", "main")
    _git(co, "branch", "feat")
    _git(co, "checkout", "-q", "feat")
    (co / "feature.txt").write_text("feature\n")
    _git(co, "add", "feature.txt")
    _git(co, "commit", "-q", "-m", "feature")
    _git(co, "checkout", "-q", "main")
    (co / "work-in-progress.txt").write_text("not committed\n")
    note = home / ".delfin" / "dashboard_here.json"
    note.parent.mkdir()
    note.write_text('{"pid": 1, "from": "the running dashboard"}')
    env = {k: v for k, v in os.environ.items()}
    env.update({
        "HOME": str(home), "XDG_RUNTIME_DIR": str(runtime),
        "FAKE_OUT": str(tmp_path / "out"),
        "JUPYTER_TOKEN": "the-main-dashboards-token",
        "OPENAI_API_KEY": "sk-should-not-pass",
        "MY_SERVICE_TOKEN": "nor-this",
        "VSCODE_IPC_HOOK_CLI": "/tmp/vscode.sock",
        "BROWSER": "/usr/bin/browser.sh",
        "DELFIN_LIFELINE_PID": "12345",
        "PYTHONPATH": str(co),                 # the shared checkout: must go
        "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_SYSTEM": os.devnull,
    })
    started = []
    yield co, env, tmp_path, started
    for p in started:
        try:
            os.killpg(p, 9)
        except OSError:
            pass


def _state(co):
    return (_git(co, "branch", "--show-current"), _git(co, "rev-parse", "HEAD"),
            _git(co, "stash", "list"), (co / "work-in-progress.txt").read_text())


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    try:
        st = Path(f"/proc/{pid}/stat").read_text().split()
        return st[2] != "Z"
    except OSError:
        return True


def test_preview_runs_from_its_own_worktree_and_closes_cleanly(setup):
    co, env, tmp, started = setup
    before = _state(co)
    ok, msg, p = BP.start(co, "feat", python=sys.executable, base_env=env,
                          wait_s=30)
    assert ok, msg
    started.append(p.pid)
    wt = Path(p.worktree)
    assert wt == co.parent / "delfin-preview" / "feat"
    assert wt.resolve() in Path(p.delfin_file).resolve().parents
    assert p.attached and _git(wt, "branch", "--show-current") == "feat"
    assert p.token == "secret-token-123"
    # the address is the one the server redirects to, not a guess: the
    # launcher serves the packaged notebook in place when it lies under
    # the root (measured in the end-to-end run: the guessed
    # delfin_voila_runtime/ path answered 404)
    assert p.url_path == "/voila/render/delfin/notebooks/delfin_dashboard.ipynb"
    assert _state(co) == before                     # nothing moved here
    seen = json.loads((tmp / "out" / "env.json").read_text())
    assert seen["PYTHONPATH"] == str(wt.resolve())
    for gone in ("JUPYTER_TOKEN", "OPENAI_API_KEY", "MY_SERVICE_TOKEN",
                 "VSCODE_IPC_HOOK_CLI", "BROWSER", "DELFIN_LIFELINE_PID"):
        assert gone not in seen, gone
    assert seen["DELFIN_NO_BROWSER"] == "1"
    assert seen[BP.PREVIEW_ENV_BRANCH] == "feat"
    assert seen["HOME"] == env["HOME"]
    # the running dashboard's note is back in place
    note = Path(env["HOME"]) / ".delfin" / "dashboard_here.json"
    assert "the running dashboard" in note.read_text()
    listed = BP.list_previews(co)
    assert [x.slug for x in listed] == ["feat"]
    rec = BP._record_path("feat")
    assert oct(rec.stat().st_mode & 0o777) == "0o600"
    child = int((tmp / "out" / "child.pid").read_text())
    assert _alive(p.pid) and _alive(child)

    ok, msg, ask = BP.stop("feat")
    assert ok and not ask, msg
    for _ in range(50):
        if not _alive(p.pid) and not _alive(child):
            break
        time.sleep(0.1)
    assert not _alive(p.pid) and not _alive(child)
    assert not wt.exists()
    assert BP.list_previews(co) == []
    assert _state(co) == before
    assert _git(co, "rev-parse", "--verify", "feat")      # branch kept


def test_a_held_branch_is_previewed_detached(setup):
    co, env, _tmp, started = setup
    ok, msg, p = BP.start(co, "main", python=sys.executable, base_env=env,
                          wait_s=30)
    assert ok, msg
    started.append(p.pid)
    assert not p.attached
    assert _git(p.worktree, "branch", "--show-current") == ""
    assert _git(p.worktree, "rev-parse", "HEAD") == _git(co, "rev-parse", "main")
    BP.stop(p.slug)


def test_a_preview_with_changes_asks_before_its_worktree_goes(setup):
    co, env, _tmp, started = setup
    ok, msg, p = BP.start(co, "feat", python=sys.executable, base_env=env,
                          wait_s=30)
    assert ok, msg
    started.append(p.pid)
    (Path(p.worktree) / "edited-in-preview.txt").write_text("keep me?\n")
    ok, msg, ask = BP.stop(p.slug)
    assert ok and ask and "edited-in-preview.txt" in msg
    assert Path(p.worktree).exists()
    ok, msg = BP.finish_close(p.slug, discard=False)
    assert ok and Path(p.worktree).exists()
    # asked again, and this time discarded: logged with the head commit
    BP._save_pending(p)
    ok, msg = BP.finish_close(p.slug, discard=True)
    assert ok, msg
    assert not Path(p.worktree).exists()
    from delfin.dashboard import branch_cleanup as BC
    rows = [json.loads(x) for x in BC._log_path().read_text().splitlines()]
    assert rows[-1]["kind"] == "worktree" and rows[-1]["sha"]
    assert "discarded" in rows[-1]["reason"]


def test_the_import_check_refuses_code_from_elsewhere(setup, tmp_path):
    co, env, _tmp, _started = setup
    target, _ = BP.resolve_target(co, "feat")
    ok, _msg, wt, _att = BP.ensure_worktree(co, target)
    assert ok
    good = BP.preview_env(env, wt, target)
    assert BP.verify_import(sys.executable, good, wt)[0]
    bad = dict(good, PYTHONPATH=str(co))            # the shared checkout
    ok, file, msg = BP.verify_import(sys.executable, bad, wt)
    assert not ok and "Refused" in msg and str(co) in file


def test_no_stash_or_checkout_in_the_shared_checkout(setup):
    """The control for the rule: the module never names stash, and every
    checkout it runs is in the preview worktree (by reading the source)."""
    src = Path(BP.__file__).read_text()
    assert '"stash"' not in src and "'stash'" not in src
    for line in src.splitlines():
        if '"checkout"' in line:
            assert "BL._git(path," in line, line
