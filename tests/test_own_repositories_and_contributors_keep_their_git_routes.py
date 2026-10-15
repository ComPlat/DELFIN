"""The shared-checkout rule touches DELFIN's checkout and nothing else.

People run the agent on DELFIN and on their own repositories. Driven
through the real executor and gate (``_DocToolExecutor.execute``), with
repositories made here, a bare ``origin`` each, and a stand-in ``gh`` on
PATH that records its arguments:

(a) in a repository that is not the shared DELFIN checkout (nor a worktree
    of it) ``git switch``, ``git switch -c`` and ``git checkout -b`` run
    exactly as before: no refusal, no question;
(b) ``git push origin main`` there (git role maintainer for that remote,
    the user's message asked for the push) is asked once by default and
    runs; with a per-repository allow rule (``permissions.repos``) it runs
    without asking -- and that allowance reaches neither another
    repository nor the DELFIN checkout. Without the user's request a push
    is asked (or refused unattended) whatever the rules say: that is the
    push gate's one-request grant, unchanged here;
(c) in the DELFIN shared checkout a branch move is refused with the route
    message, and the route itself -- ``enter_worktree(base_ref=
    "origin/main")``, edit, commit, ``git push -u origin
    HEAD:<user>/<topic>``, ``gh pr create --head <user>/<topic>`` -- runs
    in default and acceptEdits, asking only where pushing and opening a
    pull request ask today, and the shared checkout stays on main.

The origins and the stand-in gh are added to the session's roots, as a
GitHub remote is reached over the network: a file remote outside the
sandbox's roots is not visible to a sandboxed command.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from delfin import user_settings as US
from delfin.agent import api_client as A
from delfin.agent import managed_policy as MP
from delfin.agent import permission_rules as PR
from delfin.agent import shared_checkout as SC
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


def _git(cwd, *args) -> str:
    env = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
           "GIT_CONFIG_SYSTEM": os.devnull}
    return subprocess.run(["git", *args], cwd=str(cwd), env=env, check=True,
                          capture_output=True, text=True).stdout.strip()


class Asker:
    def __init__(self):
        self.asked: list = []

    def __call__(self, name, args, preview):
        self.asked.append(str((args or {}).get("command") or name))
        return True


def _repo(tmp, name) -> tuple[Path, Path]:
    origin = tmp / f"{name}.git"
    _git(tmp, "init", "-q", "--bare", "-b", "main", str(origin))
    ws = tmp / name
    _git(tmp, "clone", "-q", str(origin), str(ws))
    _git(ws, "checkout", "-q", "-b", "main")
    _git(ws, "config", "user.email", "u@example.invalid")
    _git(ws, "config", "user.name", "u")
    (ws / "f.txt").write_text("x\n")
    _git(ws, "add", "f.txt")
    _git(ws, "commit", "-qm", "init")
    _git(ws, "push", "-q", "-u", "origin", "main")
    return ws, origin


@pytest.fixture
def world(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    settings = home / US.SETTINGS_FILE_NAME
    monkeypatch.setattr(
        US, "get_settings_path",
        lambda base_path=None: (Path(base_path).expanduser() if base_path
                                else settings))
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(tmp_path / "no-policy.json"),))
    monkeypatch.setattr(A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", os.devnull)
    MP.clear_cache()
    PR._CACHE.clear()
    own, own_origin = _repo(tmp_path, "own")
    other, other_origin = _repo(tmp_path, "other")
    delfin, delfin_origin = _repo(tmp_path, "delfin-shared")
    monkeypatch.setenv(SC.ENV, str(delfin))
    fake = tmp_path / "bin"
    fake.mkdir()
    (fake / "gh").write_text(
        f"#!{sys.executable}\nimport json, os, sys\n"
        "with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"
        " 'calls.jsonl'), 'a') as fh:\n"
        "    fh.write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "if sys.argv[1:3] == ['pr', 'create']:\n"
        "    print('https://github.com/o/r/pull/7')\n")
    (fake / "gh").chmod(0o755)
    monkeypatch.setenv("PATH", f"{fake}{os.pathsep}{os.environ['PATH']}")
    w = {"tmp": tmp_path, "settings": settings, "own": own,
         "own_origin": own_origin, "other": other, "other_origin": other_origin,
         "delfin": delfin, "delfin_origin": delfin_origin, "fake": fake}
    yield w
    MP.clear_cache()
    PR._CACHE.clear()


def _settings(w, **block):
    w["settings"].write_text(json.dumps(block))
    PR._CACHE.clear()


def _perms(w, ws, mode="default"):
    p = KitToolPermissions(workspace=ws, mode=mode, confirm_callback=Asker())
    for extra in ("own_origin", "other_origin", "delfin_origin", "fake"):
        p.add_extra_dir(w[extra])
    return p


def _run(perms, cmd, cwd=None):
    asker = perms.confirm_callback
    n = len(asker.asked)
    args = {"command": cmd}
    if cwd:
        args["cwd"] = str(cwd)
    out = str(_DocToolExecutor().execute("bash", args, perms))
    return out, asker.asked[n:]


def _ran(out) -> bool:
    return '"exit_code": 0' in out


# -- (a) ----------------------------------------------------------------------

@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
def test_a_branch_moves_in_an_own_repository_run_as_before(world, mode):
    _settings(world)
    ws = world["own"]
    perms = _perms(world, ws, mode)
    # The questions are the ones the base commit (bad5a550) asks, measured
    # with this same test there: `git switch -c <new>` is asked once, the
    # approved shape `git switch -c pr/<name> origin/main` and the other
    # three run unasked. This change adds none and removes none.
    for cmd, branch, asks in (("git switch -c feat-a", "feat-a", True),
                              ("git checkout -b feat-b", "feat-b", False),
                              ("git switch main", "main", False),
                              ("git checkout feat-a", "feat-a", False),
                              ("git switch -c pr/x origin/main", "pr/x", False)):
        out, asked = _run(perms, cmd)
        assert _ran(out) and asked == ([cmd] if asks else []), (cmd, out, asked)
        assert "shared checkout" not in out
        assert _git(ws, "branch", "--show-current") == branch


# -- (b) ----------------------------------------------------------------------

def _maintainer_of(w, *origins, allow_in=None):
    roles = {str(Path(o).resolve()): "maintainer" for o in origins}
    block = {"agent": {"git_roles": roles}}
    if allow_in is not None:
        block["permissions"] = {"repos": {str(Path(allow_in).resolve()): {
            "allow": ["bash(git push origin main)"]}}}
    _settings(w, **block)


def _asked_for_a_push(perms):
    """The user's message asked for the push: the one-request grant the
    push gate wants before any push (unchanged by this change)."""
    A._grant_push_from(perms, "push it to main", new_request=True)
    return perms


def _new_commit(ws, name):
    (ws / name).write_text(name)
    _git(ws, "add", name)
    _git(ws, "commit", "-qm", name)
    return _git(ws, "rev-parse", "HEAD")


def test_b_a_push_to_main_in_an_own_repository_is_asked_once(world):
    _maintainer_of(world, world["own_origin"], world["other_origin"])
    sha = _new_commit(world["own"], "g.txt")
    out, asked = _run(_asked_for_a_push(_perms(world, world["own"])),
                      "git push origin main")
    assert asked == ["git push origin main"], (out, asked)
    assert _ran(out), out
    assert _git(world["own_origin"], "rev-parse", "main") == sha


def test_b_an_allow_rule_for_the_repository_runs_it_unasked_there_only(world):
    # A permissions.repos key starting with "/" names the checkout's
    # location (permission_rules.repo_entry_matches), so the rule is
    # written for the own checkout's path.
    _maintainer_of(world, world["own_origin"], world["other_origin"],
                   allow_in=world["own"])
    sha = _new_commit(world["own"], "g.txt")
    out, asked = _run(_asked_for_a_push(_perms(world, world["own"])),
                      "git push origin main")
    assert asked == [] and _ran(out), (out, asked)
    assert _git(world["own_origin"], "rev-parse", "main") == sha
    # another repository: asked, as without the rule
    other_sha = _new_commit(world["other"], "h.txt")
    out, asked = _run(_asked_for_a_push(_perms(world, world["other"])),
                      "git push origin main")
    assert asked == ["git push origin main"], (out, asked)
    assert _git(world["other_origin"], "rev-parse", "main") == other_sha
    # the DELFIN checkout (contributor there): not pushed, rule or not
    before = _git(world["delfin_origin"], "rev-parse", "main")
    _new_commit(world["delfin"], "i.txt")
    out, asked = _run(_asked_for_a_push(_perms(world, world["delfin"])),
                      "git push origin main")
    assert not _ran(out) or asked, (out, asked)
    assert _git(world["delfin_origin"], "rev-parse", "main") == before


# -- (c) ----------------------------------------------------------------------

def test_c_a_branch_move_in_the_shared_checkout_names_the_route(world):
    _settings(world)
    out, asked = _run(_perms(world, world["delfin"]),
                      "git switch -c pr/x origin/main")
    assert "shared checkout" in out and not _ran(out)
    # the message names the route the next test drives
    for part in ('enter_worktree(base_ref=\\"origin/main\\")',
                 "git push -u origin HEAD:<user>/<topic>",
                 "gh pr create --head <user>/<topic>"):
        assert part in out, (part, out)
    assert asked == []
    assert _git(world["delfin"], "branch", "--show-current") == "main"


@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
def test_c_the_pull_request_route_runs_from_a_worktree(world, mode):
    _settings(world)
    perms = _perms(world, world["delfin"], mode)
    ex = _DocToolExecutor()
    info = json.loads(ex.execute("enter_worktree", {"base_ref": "origin/main"},
                                 perms))
    wt = info["path"]
    assert "error" not in info, info
    out = ex.execute("write_file", {"path": f"{wt}/new.txt", "content": "n\n"},
                     perms)
    assert "File created" in str(out), out
    out, asked = _run(perms, "git add new.txt && git commit -m 'Add new'", wt)
    assert _ran(out) and asked == [], (out, asked)
    out, asked = _run(perms, "git push -u origin HEAD:contrib/add-new", wt)
    assert _ran(out) and asked == ["git push -u origin HEAD:contrib/add-new"], \
        (out, asked)
    out, asked = _run(perms,
                      "gh pr create --base main --head contrib/add-new --fill", wt)
    assert _ran(out) and "pull/7" in out, out
    assert asked == ["gh pr create --base main --head contrib/add-new --fill"]
    calls = [json.loads(x) for x in
             (world["fake"] / "calls.jsonl").read_text().splitlines()]
    assert ["pr", "create", "--base", "main", "--head", "contrib/add-new",
            "--fill"] in calls
    assert "contrib/add-new" in _git(world["delfin_origin"], "branch", "--list")
    assert _git(world["delfin"], "branch", "--show-current") == "main"
    assert _git(world["delfin"], "status", "--porcelain") == ""
