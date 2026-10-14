"""The user's allow rules are scoped to a repository.

Input: ``permissions.repos`` in the user's settings file, keyed like
``agent.git_roles`` (a normalised remote URL, else the top-level path) or
by a path glob; the global ``permissions.allow`` list and
``permissions.global_allow``; tool calls through the real executor.

Pinned: a per-repository allow answers in that repository and in no other
(remote key and path key); a global allow still answers everywhere, and
the start notice names it with the way to scope it; ``global_allow:
false`` turns the global list off while the per-repository lists stay;
deny and ask stay global; a per-repository deny or ask applies in its
repository only; a project file cannot add a ``repos`` block; the managed
``only_allow`` bounds a per-repository allow; ``/permissions`` and
``capabilities`` name the repository and the scope of every allow.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import pytest

from conftest import child_env
from delfin import user_settings as US
from delfin.agent import api_client as A
from delfin.agent import managed_policy as MP
from delfin.agent import permission_rules as PR
from delfin.agent import security_events
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


class Asker:
    def __init__(self, answer: bool):
        self.answer = answer
        self.asked: list = []

    def __call__(self, name, args, preview):
        self.asked.append(name)
        return self.answer


def _git(tmp_path, where, *args):
    subprocess.run(["git", *args], cwd=str(where), check=True,
                   capture_output=True, env=child_env(tmp_path))


def _repo(tmp_path, name, remote=None):
    ws = tmp_path / name
    ws.mkdir()
    _git(tmp_path, ws, "init", "-q")
    if remote:
        _git(tmp_path, ws, "remote", "add", "origin", remote)
    run = ws / "run.sh"
    run.write_text("#!/bin/sh\necho RUN_MARKER\n")
    run.chmod(0o755)
    return ws


@pytest.fixture
def env(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    settings = home / US.SETTINGS_FILE_NAME
    monkeypatch.setattr(
        US, "get_settings_path",
        lambda base_path=None: (Path(base_path).expanduser() if base_path
                                else settings))
    etc = tmp_path / "etc" / "delfin"
    etc.mkdir(parents=True)
    etc.chmod(0o755)
    etc.parent.chmod(0o755)
    policy = etc / "managed_policy.json"
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(policy),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE",
                        frozenset({0, os.getuid()}))
    monkeypatch.setattr(MP, "_ISOLATION_PROBE", lambda backend="": True)
    monkeypatch.setattr(A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)
    a = _repo(tmp_path, "a", "git@GitHub.com:me/a.git")
    b = _repo(tmp_path, "b", "https://github.com/me/b")
    local = _repo(tmp_path, "local")
    MP.clear_cache()
    PR._CACHE.clear()
    security_events.clear()
    yield {"settings": settings, "policy": policy, "a": a, "b": b,
           "local": local, "tmp": tmp_path}
    MP.clear_cache()
    PR._CACHE.clear()


def _user(env, **block):
    env["settings"].write_text(json.dumps({"permissions": block}))
    PR._CACHE.clear()


def _asked(ws) -> tuple:
    """(ran, asked) for ``./run.sh`` in ``ws`` with a 'no' answer."""
    perms = KitToolPermissions(workspace=ws, mode="default")
    perms.confirm_callback = Asker(False)
    out = _DocToolExecutor().execute("bash", {"command": "./run.sh"}, perms)
    return "RUN_MARKER" in out, bool(perms.confirm_callback.asked)


def test_control_without_a_rule_the_command_is_asked(env):
    _user(env)
    assert _asked(env["a"]) == (False, True)


def test_a_repo_allow_by_remote_key_applies_there_and_nowhere_else(env):
    _user(env, repos={"github.com/me/a": {"allow": ["bash(./run.sh)"]}})
    assert _asked(env["a"]) == (True, False)
    assert _asked(env["b"]) == (False, True)
    assert _asked(env["local"]) == (False, True)
    rs = PR.load(env["a"])
    assert rs.repo_key == "github.com/me/a"
    assert [(r.text, r.scope) for r in rs.allow] == [
        ("bash(./run.sh)", "github.com/me/a")]
    assert PR.load(env["b"]).allow == ()


def test_a_repo_key_may_be_written_as_a_url(env):
    _user(env, repos={"https://github.com/me/b.git": {
        "allow": ["bash(./run.sh)"]}})
    assert _asked(env["b"]) == (True, False)
    assert _asked(env["a"]) == (False, True)


def test_a_repo_without_a_remote_is_keyed_by_its_top_level_path(env):
    top = str(env["local"].resolve())
    _user(env, repos={top: {"allow": ["bash(./run.sh)"]}})
    assert PR.load(env["local"]).repo_key == top
    assert _asked(env["local"]) == (True, False)
    assert _asked(env["a"]) == (False, True)
    # A subdirectory of the checkout is the same repository.
    sub = env["local"] / "sub"
    sub.mkdir()
    assert PR.load(sub).repo_key == top


def test_a_path_glob_key_matches_the_checkout_location(env):
    _user(env, repos={str(env["tmp"]) + "/[ab]": {}, str(env["tmp"]) + "/a*":
                      {"allow": ["bash(./run.sh)"]}})
    assert _asked(env["a"]) == (True, False)
    assert _asked(env["b"]) == (False, True)


def test_a_global_allow_still_applies_and_the_start_says_so(env):
    _user(env, allow=["bash(./run.sh)"])
    assert _asked(env["a"]) == (True, False)
    assert _asked(env["b"]) == (True, False)
    notice = PR.scope_notice(PR.load(env["a"]))
    assert notice and "bash(./run.sh)" in notice
    assert "permissions.repos" in notice and "global_allow" in notice
    from delfin.agent import cli as agent_cli
    engine = type("E", (), {"workspace": env["a"],
                            "client": type("C", (), {"model": "m"})(),
                            "provider": "openai"})()
    args = argparse.Namespace(backend="", provider="", model="", mode="",
                              cwd=str(env["a"]), effort="",
                              permission_mode="", settings_defaults=False)
    notes = agent_cli._bounding_notices(args, engine)
    assert any("bash(./run.sh)" in n and "every repository" in n
               for n in notes), notes


def test_an_explicit_global_allow_true_is_kept_without_the_notice(env):
    _user(env, allow=["bash(./run.sh)"], global_allow=True)
    assert _asked(env["b"]) == (True, False)
    assert PR.scope_notice(PR.load(env["b"])) is None


def test_global_allow_false_turns_the_global_list_off(env):
    _user(env, allow=["bash(./run.sh)"], global_allow=False,
          repos={"github.com/me/a": {"allow": ["bash(./run.sh)"]}})
    assert _asked(env["a"]) == (True, False)
    assert _asked(env["b"]) == (False, True)
    rs = PR.load(env["b"])
    assert rs.allow == ()
    assert [r.text for r in rs.allow_off] == ["bash(./run.sh)"]
    assert PR.scope_notice(rs) is None
    assert "global_allow is false" in PR.describe(rs)


def test_a_global_allow_that_is_not_a_bool_is_an_error_and_turns_it_off(env):
    _user(env, allow=["bash(./run.sh)"], global_allow="no")
    rs = PR.load(env["a"])
    assert rs.allow == ()
    assert any("global_allow" in line for line in PR.error_lines(rs))


def test_deny_and_ask_stay_global(env):
    _user(env, deny=["bash(./run.sh)"],
          repos={"github.com/me/a": {"allow": ["bash(./run.sh)"]}})
    for ws in (env["a"], env["b"], env["local"]):
        assert _asked(ws)[0] is False
    _user(env, ask=["bash(./run.sh)"],
          repos={"github.com/me/a": {"allow": ["bash(./run.sh)"]}})
    assert _asked(env["a"]) == (False, True)


def test_a_repo_deny_or_ask_applies_in_that_repo_only(env):
    _user(env, allow=["bash(./run.sh)"], global_allow=True,
          repos={"github.com/me/a": {"deny": ["bash(./run.sh)"]},
                 "github.com/me/b": {"ask": ["bash(./run.sh)"]}})
    assert _asked(env["a"]) == (False, False)
    assert _asked(env["b"]) == (False, True)
    assert _asked(env["local"]) == (True, False)


def test_a_project_file_cannot_add_a_repos_block(env):
    _user(env)
    p = env["a"] / ".delfin" / "settings.json"
    p.parent.mkdir()
    p.write_text(json.dumps({"permissions": {
        "repos": {"github.com/me/a": {"allow": ["bash(./run.sh)"]}},
        "global_allow": True}}))
    PR._CACHE.clear()
    rs = PR.load(env["a"])
    assert rs.allow == ()
    assert any("cannot grant" in line for line in PR.error_lines(rs))
    assert _asked(env["a"]) == (False, True)


def test_a_broken_repos_block_is_reported_and_grants_nothing(env):
    _user(env, repos={"github.com/me/a": ["bash(./run.sh)"],
                      "not a key": {"allow": ["bash(./run.sh)"]},
                      "github.com/me/b": {"allow": ["bash(a ; b)"],
                                          "frob": []}})
    rs = PR.load(env["a"])
    assert rs.allow == ()
    lines = "\n".join(PR.error_lines(rs))
    assert "github.com/me/a" in lines and "not a key" in lines
    assert "frob" in lines and "bash(a ; b)" in lines


def test_managed_only_allow_bounds_a_repo_allow(env):
    env["policy"].write_text(json.dumps(
        {"permissions": {"only_allow": ["bash(ls:*)"]}}))
    env["policy"].chmod(0o644)
    MP.clear_cache()
    _user(env, repos={"github.com/me/a": {"allow": ["bash(./run.sh)"]}})
    assert _asked(env["a"]) == (False, True)
    env["policy"].write_text(json.dumps(
        {"permissions": {"only_allow": ["bash(./run.sh)"]}}))
    MP.clear_cache()
    PR._CACHE.clear()
    assert _asked(env["a"]) == (True, False)


def test_permissions_and_capabilities_name_the_repository_and_scope(env):
    _user(env, allow=["bash(make:*)"],
          repos={"github.com/me/a": {"allow": ["bash(./run.sh)"]},
                 "github.com/me/b": {"allow": ["bash(pytest:*)"]}})
    text = PR.describe(PR.load(env["a"]))
    assert "repository github.com/me/a" in text
    assert "bash(./run.sh)  [user, repository github.com/me/a]" in text
    assert "bash(make:*)  [user, every repository]" in text
    assert "bash(pytest:*)" in text and "not this repository" in text
    from delfin.agent import boundaries
    recs = boundaries._rules_summary(env["a"])["rules"]
    scoped = [r for r in recs if r["rule"] == "bash(./run.sh)"]
    assert scoped and scoped[0]["scope"] == "github.com/me/a"
    glob = [r for r in recs if r["rule"] == "bash(make:*)"]
    assert glob and glob[0]["scope"] == "every repository"


def test_the_repo_key_is_the_git_role_key(env):
    from delfin.agent import git_role
    for ws in (env["a"], env["b"], env["local"]):
        assert PR.load(ws).repo_key == git_role.repo_key(ws)
