"""The agent's own source edits and scratch directories are not refused as
secrets or as a destructive delete; the credential stores still are.

Measured on the shared bug archive, 2026-09-16 .. 10-14 (refusals in
tool_trace, 33 sessions):

  secret-deny     15 calls / 7 sessions; DELFIN's source module
                  delfin/agent/credentials.py among them: `git diff`,
                  `git checkout`, `git reset` of it refused while the
                  session was editing it, then 4 more calls reworking
                  the staging around the file.
  rm deny-pattern 11 calls / 7 sessions, mostly `rm -rf` of a scratch
                  directory the session had made itself (`.tmpdir`,
                  `bench_trial/`), e.g. `mkdir -p .tmpdir && TMPDIR=...
                  git push ...; rm -rf .tmpdir`.

The name glob ``credentials*`` already let a committed, UNMODIFIED source
module through (_reviewed_source_named_like_a_secret). A module the
session itself is editing is modified, so the exception stopped exactly
when the agent needed to diff, stage or restore its own work.

  tracked .py, every change since HEAD this session's    readable
  tracked .py with a change nobody journalled            refused
  the session's edit on top of someone else's            refused
  ~/.delfin/credentials.json, ~/.ssh, .env, *.key ...     refused
  rm -rf of a directory whose files this session made     allowed
  rm -rf of a directory the same command creates          allowed
  rm -rf of anything else (tracked, foreign, outside,
    a glob, a symlink)                                    refused
"""

from __future__ import annotations

import os
import subprocess

import pytest

from delfin.agent import change_journal as cj
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _git(ws, *args):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    subprocess.run(["git", "-C", str(ws), "-c", "user.name=t", "-c",
                    "user.email=t@example.org", *args], check=True,
                   capture_output=True, env=env)


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(cj, "_undo_root", lambda: tmp_path / "undo")
    root = tmp_path / "ws"
    (root / "pkg").mkdir(parents=True)
    (root / "pkg" / "credentials.py").write_text(
        "def load():\n    return None\n", encoding="utf-8")
    (root / "pkg" / "other.py").write_text("x = 1\n", encoding="utf-8")
    _git(root, "init", "-q")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", "c")
    perms = KitToolPermissions(workspace=root, mode="bypassPermissions")
    perms.task_session_id = "sess-own"
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms

    def gate(cmd):
        return eng._run_permission_gate("bash", {"command": cmd}, perms)

    def read(path):
        # The read tools' own check (read_file, grep_file), where the
        # secret-deny for a read is decided.
        return eng._check_read_access(perms, root / path, path)
    return root, perms, gate, read


def _agent_edits(root, rel, new):
    """An edit as the write tools journal it: the write, then the record
    (record_change hashes the post-image from the file itself)."""
    p = root / rel
    old = p.read_text(encoding="utf-8")
    p.write_text(new, encoding="utf-8")
    cj.record_change("sess-own", tool="edit_file", path=str(p.resolve()),
                     old_text=old, new_text=new)


def _refused_as_secret(out):
    return bool(out) and "secret" in str(out)


# -- a source module named like a secret ------------------------------------

def test_the_sessions_own_edit_of_credentials_py_can_be_diffed(ws):
    root, perms, gate, read = ws
    _agent_edits(root, "pkg/credentials.py",
                 "def load():\n    return 'from the store'\n")
    for cmd in ("git diff HEAD pkg/credentials.py",
                "git checkout HEAD -- pkg/credentials.py",
                "git reset HEAD pkg/credentials.py"):
        assert not _refused_as_secret(gate(cmd)), cmd
    assert not _refused_as_secret(read("pkg/credentials.py"))


def test_an_edit_nobody_journalled_stays_refused(ws):
    root, perms, gate, read = ws
    (root / "pkg" / "credentials.py").write_text(
        "API_KEY = 'typed by the user'\n", encoding="utf-8")
    assert _refused_as_secret(gate("git diff HEAD pkg/credentials.py"))
    assert _refused_as_secret(read("pkg/credentials.py"))


def test_the_sessions_edit_on_top_of_a_foreign_one_stays_refused(ws):
    root, perms, gate, read = ws
    (root / "pkg" / "credentials.py").write_text(
        "API_KEY = 'typed by the user'\n", encoding="utf-8")
    _agent_edits(root, "pkg/credentials.py",
                 "API_KEY = 'typed by the user'\n# note\n")
    assert _refused_as_secret(gate("git diff HEAD pkg/credentials.py"))


@pytest.mark.parametrize("cmd", [
    "cat ~/.delfin/credentials.json",
    "cat ~/.ssh/id_rsa",
    "cat .env",
    "cat credentials.json",
    "head secrets.yaml",
    "cat ~/.git-credentials",
])
def test_the_credential_stores_stay_refused(ws, cmd):
    root, perms, gate, read = ws
    (root / ".env").write_text("K=V\n", encoding="utf-8")
    (root / "credentials.json").write_text("{}\n", encoding="utf-8")
    (root / "secrets.yaml").write_text("k: v\n", encoding="utf-8")
    assert _refused_as_secret(gate(cmd)), cmd


# -- rm -rf of the session's own scratch -------------------------------------

def _agent_creates(root, rel, text="x\n"):
    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    cj.record_change("sess-own", tool="bash", path=str(p.resolve()),
                     old_text=None, new_text=text)


def _refused_as_delete(out):
    return bool(out) and "deny-pattern" in str(out)


def test_rm_rf_of_the_sessions_scratch_dir_runs(ws):
    root, perms, gate, read = ws
    _agent_creates(root, "bench_trial/run.sh")
    _agent_creates(root, "bench_trial/out/log.txt")
    assert not _refused_as_delete(gate("rm -rf bench_trial"))
    assert not _refused_as_delete(gate("rm -rf bench_trial/ && ls"))


def test_rm_rf_of_a_dir_the_command_itself_creates_runs(ws):
    root, perms, gate, read = ws
    cmd = ("mkdir -p .tmpdir && TMPDIR=$PWD/.tmpdir make check 2>&1; "
           "rc=$?; rm -rf .tmpdir; exit $rc")
    assert not _refused_as_delete(gate(cmd))


def test_an_empty_dir_can_go(ws):
    root, perms, gate, read = ws
    (root / "empty").mkdir()
    assert not _refused_as_delete(gate("rm -rf empty"))


@pytest.mark.parametrize("setup,cmd", [
    ("foreign", "rm -rf data"),           # a file nobody journalled
    ("tracked", "rm -rf pkg"),            # tracked sources
    ("mixed", "rm -rf bench_trial data"), # one own, one foreign
    ("own", "rm -rf bench_*"),            # a glob
    ("own", "rm -rf /tmp/elsewhere"),     # outside the workspace
    ("own", "rm -rf ~"),
    ("link", "rm -rf linkdir/"),          # a symlink
    ("moved", "mkdir -p d && mv pkg d/ && rm -rf d"),
    ("own", "rm -rf bench_trial && sudo ls"),  # another deny still fires
])
def test_everything_else_stays_refused(ws, setup, cmd):
    root, perms, gate, read = ws
    if setup in ("own", "mixed"):
        _agent_creates(root, "bench_trial/run.sh")
    if setup in ("foreign", "mixed"):
        (root / "data").mkdir()
        (root / "data" / "results.csv").write_text("1\n", encoding="utf-8")
    if setup == "link":
        target = root.parent / "outside"
        target.mkdir()
        (root / "linkdir").symlink_to(target, target_is_directory=True)
    out = gate(cmd)
    assert out and ("deny-pattern" in out or "rejected" in out), (cmd, out)
