"""``delfin-agent ops status``: one read-only overview of DELFIN's machinery.

Input under test: a temporary repository with a main checkout, a
worktree detached at the branch the main checkout holds, a worktree
whose directory is gone, an agent worktree no process works in, two
``index.lock`` files (one with a git process of the user in that
checkout, one without), a stand-in process table (``/proc`` layout) and
a temporary HOME with session locks and a kept-session record.

Asserted: the rows and the ``problem -> safe next step`` lines; the
``--json`` schema; that the permission gate auto-allows the command and
its flags while refusing what rides along; and that nothing in the tree
changes (an mtime/size snapshot of every file before and after).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from delfin.agent import ops_status

HOST = "node-a"
NOW = 2_000_000_000.0
UID = 4242
OTHER_UID = 4343


def _git(cwd: Path, *args: str) -> str:
    env = dict(os.environ)
    env.update(GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@example.invalid",
               GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@example.invalid",
               GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1")
    out = subprocess.run(["git", *args], cwd=str(cwd), env=env,
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()


def _proc(root: Path, pid: int, *, argv: list[str], cwd: Path,
          uid: int = UID, ppid: int = 1, environ: dict | None = None) -> None:
    d = root / str(pid)
    d.mkdir(parents=True)
    name = Path(argv[0]).name[:15]
    (d / "status").write_text(
        f"Name:\t{name}\nPPid:\t{ppid}\n"
        f"Uid:\t{uid}\t{uid}\t{uid}\t{uid}\n")
    (d / "cmdline").write_bytes(b"\0".join(a.encode() for a in argv) + b"\0")
    (d / "environ").write_bytes(b"\0".join(
        f"{k}={v}".encode() for k, v in (environ or {}).items()) + b"\0")
    os.symlink(str(cwd), str(d / "cwd"))


@pytest.fixture
def world(tmp_path):
    origin = tmp_path / "origin.git"
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(tmp_path, "init", "-q", "--bare", str(origin))
    _git(repo, "init", "-q", "-b", "main")
    (repo / ".gitignore").write_text(".delfin/\n.agentx/\n")
    (repo / "a.txt").write_text("a\n")
    _git(repo, "add", ".gitignore", "a.txt")
    _git(repo, "commit", "-q", "-m", "one")
    _git(repo, "remote", "add", "origin", str(origin))
    _git(repo, "push", "-q", "-u", "origin", "main")
    (repo / "a.txt").write_text("a\nb\n")
    _git(repo, "commit", "-q", "-am", "two")          # ahead 1

    dash = tmp_path / "dash"
    _git(repo, "worktree", "add", "-q", "-b", "dash-tmp", str(dash))
    _git(dash, "checkout", "-q", "--detach", "main")  # held by repo

    gone = tmp_path / "gone"
    _git(repo, "worktree", "add", "-q", "-b", "gone", str(gone))
    shutil.rmtree(gone)

    agent_wt = repo / ".agentx" / "worktrees" / "agent-1"
    _git(repo, "worktree", "add", "-q", "-b", "agent-1", str(agent_wt))
    session_wt = repo / ".delfin" / "worktrees" / "s1"
    _git(repo, "worktree", "add", "-q", "-b", "s1", str(session_wt))

    common = repo / ".git"
    (common / "index.lock").write_text("")
    (common / "worktrees" / "dash" / "index.lock").write_text("")

    proc = tmp_path / "proc"
    (proc / "self").mkdir(parents=True)
    (proc / "self" / "status").write_text("Name:\tpython\nNSpid:\t77\n")
    _proc(proc, 4100, argv=["/opt/env/bin/python", "/opt/env/bin/voila",
                            "--port=8866", "--no-browser", "delfin.ipynb"],
          cwd=dash, environ={"PYTHONPATH": str(dash),
                             "DELFIN_SESSION_GRACE_SECONDS": "90",
                             "SECRET_TOKEN": "do-not-print"})
    _proc(proc, 4101, argv=["/opt/env/bin/python", "-m", "ipykernel_launcher",
                            "-f", "/run/user/x/jupyter/kernel-abc123.json"],
          cwd=dash, ppid=4100)
    _proc(proc, 4200, argv=["git", "status"], cwd=repo)
    _proc(proc, 4300, argv=["/opt/env/bin/voila", "--port=9999"],
          cwd=repo, uid=OTHER_UID)

    home = tmp_path / "home"
    sessions = home / ".delfin" / "agent_sessions"
    sessions.mkdir(parents=True)
    (sessions / "s-live.lock").write_text(json.dumps(
        {"pid": 4101, "ts": NOW - 10, "host": HOST}))
    (sessions / "s-dead.lock").write_text(json.dumps(
        {"pid": 999999, "ts": NOW - 10, "host": HOST}))
    (sessions / "s-old.lock").write_text(json.dumps(
        {"pid": 4100, "ts": NOW - 7200, "host": HOST}))
    (sessions / "s-away.lock").write_text(json.dumps(
        {"pid": 5, "ts": NOW - 10, "host": "node-b"}))
    kept = home / ".delfin" / "kept_sessions"
    kept.mkdir(parents=True)
    (kept / "work.json").write_text(json.dumps(
        {"session_name": "work", "kernel_id": "abc123", "pid": 4101,
         "host": HOST, "updated_at": NOW - 5}))
    return {"tmp": tmp_path, "repo": repo.resolve(), "dash": dash.resolve(),
            "gone": gone.resolve(), "agent_wt": agent_wt.resolve(),
            "session_wt": session_wt.resolve(), "proc": proc, "home": home}


def _collect(w):
    return ops_status.collect(repo=w["dash"], proc_root=w["proc"], uid=UID,
                              home=w["home"], host=HOST, now=NOW)


def _by_path(status, path):
    for row in status["checkouts"]:
        if row["path"] == str(path):
            return row
    raise AssertionError(f"{path} not in {[r['path'] for r in status['checkouts']]}")


def test_checkouts_kinds_branch_and_upstream(world):
    st = _collect(world)
    assert st["repository"] == str(world["repo"])
    main = _by_path(st, world["repo"])
    assert main["kind"] == "main" and main["branch"] == "main"
    assert main["upstream"] == "origin/main"
    assert (main["ahead"], main["behind"]) == (1, 0)
    assert main["dirty"] is False

    dash = _by_path(st, world["dash"])
    assert dash["branch"] == "" and dash["detached"] is True
    assert dash["detached_from"] == "main"
    assert dash["held_by"] == str(world["repo"])
    assert dash["same_commit_as_branch"] is True
    assert dash["dashboards"] == [4100]

    assert _by_path(st, world["agent_wt"])["kind"] == "agent"
    assert _by_path(st, world["session_wt"])["kind"] == "session"
    gone = _by_path(st, world["gone"])
    assert gone["stale"] is True and gone["exists"] is False


def test_the_held_branch_line_names_the_safe_next_step(world):
    text = ops_status.format_text(_collect(world))
    repo, dash = world["repo"], world["dash"]
    line = next(ln for ln in text.splitlines() if "is held by" in ln)
    assert f"main is held by {repo}" in line
    assert "this dashboard (pid 4100) runs detached at the same commit" in line
    assert "→" in line
    assert f"git -C {repo} checkout --detach" in line
    assert f"git -C {dash} switch main" in line
    assert "detached @ " in text and f"(from main, held by {repo})" in text


def test_index_locks_with_and_without_a_git_process(world):
    st = _collect(world)
    locks = {Path(lk["path"]).parent.name: lk
             for lk in st["leftovers"]["index_locks"]}
    assert locks[".git"]["git_pids"] == [4200]
    assert locks[".git"]["stale"] is False
    assert locks["dash"]["git_pids"] == [] and locks["dash"]["stale"] is True
    text = ops_status.format_text(st)
    assert any("index.lock" in ln and "4200" in ln and "wait" in ln
               for ln in text.splitlines())
    stale_line = next(ln for ln in text.splitlines()
                      if "stale" in ln and "worktrees/dash/index.lock" in ln)
    assert "rm " in stale_line and "→" in stale_line


def test_stale_and_orphaned_worktrees_are_flagged(world):
    st = _collect(world)
    left = st["leftovers"]
    assert str(world["gone"]) in left["stale_worktrees"]
    assert str(world["agent_wt"]) in left["orphaned_agent_worktrees"]
    assert str(world["session_wt"]) not in left["orphaned_agent_worktrees"]
    text = ops_status.format_text(st)
    assert any(str(world["gone"]) in ln and "git -C" in ln and "worktree prune" in ln
               for ln in text.splitlines())


def test_dashboards_are_this_users_only_with_port_checkout_and_kernels(world):
    st = _collect(world)
    assert [d["pid"] for d in st["dashboards"]] == [4100]
    d = st["dashboards"][0]
    assert d["port"] == 8866
    assert d["cwd"] == str(world["dash"])
    assert d["checkout"] == str(world["dash"])
    assert d["pythonpath"] == [str(world["dash"])]
    assert d["grace_s"] == 90.0
    assert d["kernels"] == [{"kernel_id": "abc123", "pid": 4101,
                             "kept_as": "work", "turn_running": False}]
    dumped = json.dumps(st)
    assert "do-not-print" not in dumped
    assert "9999" not in json.dumps(st["dashboards"])


def test_session_locks_live_stale_and_elsewhere(world):
    st = _collect(world)
    locks = {lk["session_id"]: lk for lk in st["sessions"]["locks"]}
    assert locks["s-live"]["state"] == "live"
    assert locks["s-live"]["kernel_id"] == "abc123"
    assert locks["s-dead"]["state"] == "stale"
    assert "not running" in locks["s-dead"]["why"]
    assert locks["s-old"]["state"] == "stale"
    assert "older" in locks["s-old"]["why"]
    assert locks["s-away"]["state"] == "elsewhere"
    assert sorted(st["leftovers"]["stale_session_locks"]) == sorted(
        [locks["s-dead"]["path"], locks["s-old"]["path"]])
    assert st["sessions"]["kept"][0]["session_name"] == "work"
    assert st["sessions"]["kept"][0]["alive"] is True


def test_json_schema(world):
    st = _collect(world)
    assert st["schema"] == ops_status.SCHEMA
    assert set(st) >= {"schema", "generated_at", "host", "repository",
                       "process_view", "checkouts", "dashboards", "sessions",
                       "leftovers", "findings"}
    for f in st["findings"]:
        assert set(f) == {"kind", "problem", "next_step"}
    json.dumps(st)


def _snapshot(root: Path) -> dict:
    out = {}
    for dirpath, dirnames, filenames in os.walk(root):
        for name in dirnames + filenames:
            p = Path(dirpath) / name
            try:
                s = p.lstat()
            except OSError:
                continue
            out[str(p)] = (s.st_mtime_ns, s.st_size)
    return out


def test_nothing_is_written(world, monkeypatch, capsys):
    before = _snapshot(world["tmp"])
    st = _collect(world)
    ops_status.format_text(st)
    monkeypatch.setenv("HOME", str(world["home"]))
    from delfin.agent import cli
    rc = cli.main(["ops", "status", "--json", "--repo", str(world["dash"])])
    assert rc == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed["repository"] == str(world["repo"])
    rc = cli.main(["ops", "status", "--repo", str(world["repo"])])
    assert rc == 0
    assert "→" in capsys.readouterr().out
    assert _snapshot(world["tmp"]) == before


# ---------------------------------------------------------------------------
# The gate: allowed without a question in every mode that runs bash, and
# nothing may ride along on the allowance.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", ["default", "diff_approval", "acceptEdits",
                                  "bypassPermissions"])
@pytest.mark.parametrize("cmd", [
    "delfin-agent ops status",
    "delfin-agent ops status --json",
    "delfin-agent ops status --repo . --json",
    "delfin-agent ops status --help",
])
def test_the_gate_allows_it(tmp_path, mode, cmd):
    from delfin.agent import boundaries
    ws = tmp_path / "ws"
    ws.mkdir()
    said = boundaries.explain("bash", {"command": cmd}, cwd=ws,
                              permission_mode=mode, asker="nobody")
    assert said["verdict"] == "allow", said
    assert said["reason_code"] == "allowed", said


def test_the_allowance_does_not_widen_the_read_roots(tmp_path):
    """--repo naming a path outside the workspace still meets the read
    gate: the auto-allow entry decides only that nobody is asked."""
    from delfin.agent import boundaries
    ws = tmp_path / "ws"
    ws.mkdir()
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    said = boundaries.explain(
        "bash", {"command": f"delfin-agent ops status --repo {outside}"},
        cwd=ws, permission_mode="default", asker="nobody")
    assert said["verdict"] == "refuse", said
    assert said["reason_code"] == "outside_ws", said


@pytest.mark.parametrize("cmd", [
    "delfin-agent ops status --json --write",
    "delfin-agent ops statusx",
    "delfin-agent ops status --repo $(touch x)",
    "delfin-agent ops status --repo `touch x`",
    "delfin-agent ops status > out.txt",
    "delfin-agent ops status --json | tee out.txt",
    "delfin-agent ops stop",
])
def test_nothing_rides_along_on_the_pattern(cmd):
    from delfin.agent.api_client import _DEFAULT_BASH_AUTO_ALLOW
    ops = [p for p in _DEFAULT_BASH_AUTO_ALLOW if "ops" in p]
    assert len(ops) == 1
    assert re.search(ops[0], cmd) is None, cmd


def test_a_compound_is_judged_by_its_other_half(tmp_path):
    from delfin.agent import boundaries
    ws = tmp_path / "ws"
    ws.mkdir()
    said = boundaries.explain("bash", {"command": "delfin-agent ops status; rm -f x"},
                              cwd=ws, permission_mode="default", asker="nobody")
    assert said["reason_code"] != "allowed", said


def test_the_parser_has_no_writing_flag():
    from delfin.agent import cli
    parser = cli.build_parser()
    sub = next(a for a in parser._actions
               if a.__class__.__name__ == "_SubParsersAction")
    ops = sub.choices["ops"]
    ops_sub = next(a for a in ops._actions
                   if a.__class__.__name__ == "_SubParsersAction")
    status = ops_sub.choices["status"]
    flags = sorted(o for a in status._actions for o in a.option_strings)
    assert flags == sorted(["--json", "--repo", "-h", "--help"])
