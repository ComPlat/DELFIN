"""The managed policy's rules are permission rules: one parser, one matcher,
evaluated first.

Input: a managed policy file injected through
``managed_policy._PATHS_OVERRIDE`` (trusted through
``_ADMIN_UIDS_OVERRIDE``), the user's rules in a settings file under a
temporary HOME, and tool calls through ``_DocToolExecutor.execute``,
``_mcp_call_refusal``, ``cli_gate.decide`` and the bash gate's last tier
(``_run_permission_gate`` with a caged admission).

Pinned:

* ``managed_policy`` has no matcher of its own; its rules are
  ``permission_rules.Rule`` objects with ``managed=True``, loaded ahead of
  the user's and the project's, and match as widely (chains, ``sh -c``,
  wrappers, non-adjacent words).
* every deny outranks any ask; within a list the managed rule answers.
* ``only_allow`` bounds the user's allow rules (bash, read, tool), the
  auto-allow list and the caged admission.
* a managed deny or ask holds at the caged tier, which the entry points
  do not pass through when the gate is driven directly, and is reported
  to a program as ``permission_denied`` with ``managed_policy_deny`` /
  ``managed_policy_ask``.
* ``capabilities`` and ``/permissions`` list managed and user rules
  together, each with its origin and source file.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from delfin import user_settings as US
from delfin.agent import api_client as A
from delfin.agent import managed_policy as MP
from delfin.agent import permission_rules as PR
from delfin.agent import security_events
from delfin.agent import turn_outcome as TO
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


class Asker:
    def __init__(self, answer: bool):
        self.answer = answer
        self.asked: list = []

    def __call__(self, name, args, preview):
        self.asked.append((name, dict(args or {}), preview))
        return self.answer


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
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "notes.txt").write_text("NOTES_BODY\n")
    run = ws / "run.sh"
    run.write_text("#!/bin/sh\necho RUN_MARKER\n")
    run.chmod(0o755)
    far = tmp_path / "far"
    far.mkdir()
    (far / "data.txt").write_text("FAR_MARKER\n")
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
    # The question an off-list command gets is what these tests observe;
    # the caged tier is driven by its own tests below.
    monkeypatch.setattr(A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)
    MP.clear_cache()
    PR._CACHE.clear()
    security_events.clear()
    yield {"ws": ws, "policy": policy, "settings": settings, "far": far,
           "home": home, "tmp": tmp_path}
    MP.clear_cache()
    PR._CACHE.clear()


def _managed(env, **perms):
    env["policy"].write_text(json.dumps({"permissions": perms}))
    env["policy"].chmod(0o644)
    MP.clear_cache()
    PR._CACHE.clear()


def _user(env, **lists):
    env["settings"].write_text(json.dumps({"permissions": lists}))
    PR._CACHE.clear()


def _perms(env, mode="default", answer=None):
    perms = KitToolPermissions(workspace=env["ws"], mode=mode)
    if answer is not None:
        perms.confirm_callback = Asker(answer)
    return perms


def _run(perms, name, args):
    return _DocToolExecutor().execute(name, args, perms)


def _err(result: str) -> str:
    try:
        data = json.loads(result)
    except Exception:
        return ""
    return str(data.get("error") or "") if isinstance(data, dict) else ""


# ---------------------------------------------------------------------------
# One parser, one matcher
# ---------------------------------------------------------------------------

def test_the_managed_policy_has_no_matcher_of_its_own():
    for name in ("friction", "denies_path", "allow_may_grant_bash",
                 "path_matches", "_hits", "_bash_hits", "_call_paths",
                 "_words", "_glob_regex"):
        assert not hasattr(MP, name), name
    assert MP.parse_rule is PR.parse_rule
    assert MP.Rule is PR.Rule


def test_managed_rules_load_ahead_of_the_users(env):
    _managed(env, deny=["bash(sudo)"], ask=["web_search"],
             only_allow=["bash(ls:*)"])
    _user(env, deny=["bash(git push)"], ask=["bash(make)"],
          allow=["bash(pytest:*)"])
    rs = PR.load(env["ws"])
    assert [(r.text, r.managed) for r in rs.deny] == [
        ("bash(sudo)", True), ("bash(git push)", False)]
    assert [(r.text, r.managed) for r in rs.ask] == [
        ("web_search", True), ("bash(make)", False)]
    assert rs.deny[0].source == str(env["policy"])
    assert [[r.text for r in g] for g in rs.only_allow] == [["bash(ls:*)"]]


@pytest.mark.parametrize("cmd", [
    "git push", "true && git push origin x", "sh -c 'git push'",
    "/usr/bin/git push", "timeout 5 git --no-pager push",
    "cd sub && git push", "env GIT_TRACE=1 git push", "x=$(git push)",
])
def test_a_managed_deny_matches_as_widely_as_a_user_deny(env, cmd):
    _managed(env, deny=["bash(git push)"])
    managed = PR.friction(PR.load(env["ws"]), "bash", {"command": cmd},
                          env["ws"])
    assert managed is not None and managed[1].managed, cmd
    MP.clear_cache()
    env["policy"].unlink()
    _user(env, deny=["bash(git push)"])
    user = PR.friction(PR.load(env["ws"]), "bash", {"command": cmd},
                       env["ws"])
    assert user is not None and not user[1].managed, cmd
    _managed(env, deny=["bash(git push)"])
    err = _err(_run(_perms(env, "bypassPermissions"), "bash",
                    {"command": cmd}))
    assert "managed policy" in err, (cmd, err)


def test_a_managed_ask_is_put_before_a_users_allow_answers(env):
    _managed(env, ask=["bash(./run.sh)"])
    _user(env, allow=["bash(./run.sh:*)"])
    perms = _perms(env, answer=False)
    out = _run(perms, "bash", {"command": "./run.sh"})
    assert "RUN_MARKER" not in out
    assert perms.confirm_callback.asked
    assert "MANAGED POLICY" in perms.confirm_callback.asked[0][2]


def test_a_user_deny_outranks_a_managed_ask(env):
    _managed(env, ask=["bash(./run.sh)"])
    _user(env, deny=["bash(./run.sh)"])
    perms = _perms(env, answer=True)
    err = _err(_run(perms, "bash", {"command": "./run.sh"}))
    assert "permission rule" in err, err
    assert perms.confirm_callback.asked == []


def test_a_managed_deny_answers_before_a_user_deny(env):
    _managed(env, deny=["bash(./run.sh)"])
    _user(env, deny=["bash(./run.sh)"])
    err = _err(_run(_perms(env, "bypassPermissions"), "bash",
                    {"command": "./run.sh"}))
    assert "managed policy" in err and "administrator" in err, err


def test_the_mcp_route_decides_managed_and_user_rules(env):
    _managed(env, deny=["mcp__srv__drop"])
    _user(env, deny=["mcp__srv__wipe"])
    perms = _perms(env, "bypassPermissions")
    managed = json.loads(A._mcp_call_refusal("mcp__srv__drop", {}, perms))
    user = json.loads(A._mcp_call_refusal("mcp__srv__wipe", {}, perms))
    assert "managed policy" in managed["error"]
    assert "permission rule" in user["error"]


def test_the_cli_backend_hook_decides_managed_rules(env, monkeypatch):
    from delfin.agent import cli_gate, kit_settings
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        env["home"] / ".delfin" / "settings.json")
    _managed(env, deny=["bash(git push)"], ask=["read(notes.txt)"])
    token = cli_gate.new_token()
    cli_gate.write_state(token, workspace=str(env["ws"]),
                         mode="bypassPermissions", can_ask=False)
    ev = {"tool_name": "Bash",
          "tool_input": {"command": "timeout 5 git --no-pager push"},
          "cwd": str(env["ws"])}
    assert "managed policy" in (cli_gate.decide(token, ev) or "")
    ev = {"tool_name": "Read",
          "tool_input": {"file_path": str(env["ws"] / "notes.txt")},
          "cwd": str(env["ws"])}
    out = cli_gate.decide(token, ev) or ""
    assert "managed policy" in out and "nobody can be asked" in out, out


# ---------------------------------------------------------------------------
# only_allow bounds every grant
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd,inside", [
    ("ls -la", True), ("ls && pwd", True), ("ls | sort", False),
    ("ls; rm notes.txt", False), ("env ls", False), ("ls $(id)", False),
    ("ls > out.txt", False), ("sh -c ls", False), ("pwd", True),
])
def test_only_allow_admits_each_part_of_a_chain_and_nothing_else(env, cmd,
                                                                 inside):
    _managed(env, only_allow=["bash(ls:*)", "bash(pwd)"])
    assert PR.bash_within_only_allow(PR.load(env["ws"]), cmd) is inside


def test_without_only_allow_nothing_is_bounded(env):
    _managed(env, deny=["bash(sudo)"])
    assert PR.bash_within_only_allow(PR.load(env["ws"]), "make all")


def test_a_user_bash_allow_outside_only_allow_grants_nothing(env):
    _managed(env, only_allow=["bash(ls:*)"])
    _user(env, allow=["bash(./run.sh:*)"])
    perms = _perms(env, answer=False)
    out = _run(perms, "bash", {"command": "./run.sh"})
    assert "RUN_MARKER" not in out
    assert perms.confirm_callback.asked, "the allow rule answered"
    out = _run(_perms(env), "bash", {"command": "./run.sh"})
    assert "RUN_MARKER" not in out


def test_a_user_bash_allow_inside_only_allow_still_answers(env):
    _managed(env, only_allow=["bash(./run.sh:*)"])
    _user(env, allow=["bash(./run.sh:*)"])
    perms = _perms(env, answer=False)
    out = _run(perms, "bash", {"command": "./run.sh"})
    assert "RUN_MARKER" in out, out
    assert not perms.confirm_callback.asked


def test_a_user_read_allow_outside_only_allow_grants_nothing(env):
    _user(env, allow=[f"read({env['far']}/**)"])
    target = {"path": str(env["far"] / "data.txt")}
    assert "FAR_MARKER" in _run(_perms(env, answer=False), "read_file", target)
    _managed(env, only_allow=["bash(ls:*)"])
    perms = _perms(env, answer=False)
    assert "FAR_MARKER" not in _run(perms, "read_file", target)
    assert perms.confirm_callback.asked
    _managed(env, only_allow=[f"read({env['far']}/**)"])
    perms = _perms(env, answer=False)
    assert "FAR_MARKER" in _run(perms, "read_file", target)
    assert not perms.confirm_callback.asked


def test_a_user_tool_allow_outside_only_allow_answers_nothing(env):
    _user(env, allow=["mcp__srv__send"])
    rs_user = PR.load(env["ws"])
    assert PR.allows(rs_user, "tool", "mcp__srv__send", {}, env["ws"])
    _managed(env, only_allow=["bash(ls:*)"])
    rs = PR.load(env["ws"])
    assert PR.allows(rs, "tool", "mcp__srv__send", {}, env["ws"]) is None
    _managed(env, only_allow=["mcp__srv"])
    rs = PR.load(env["ws"])
    assert PR.allows(rs, "tool", "mcp__srv__send", {}, env["ws"]) is not None


def test_only_allow_bounds_the_auto_allow_list(env):
    _managed(env, only_allow=["bash(ls:*)"])
    perms = _perms(env)
    assert perms.matches_bash_auto_allow("ls -la")
    assert not perms.matches_bash_auto_allow("cat notes.txt")


# ---------------------------------------------------------------------------
# The caged tier
# ---------------------------------------------------------------------------

#: The real predicate, taken before the fixture pins it to False.
_REAL_CAGED = A._caged_command_may_run_unasked


def _cage(monkeypatch):
    """A host where the cage holds a command (as in the caged tests)."""
    monkeypatch.setattr(A, "_caged_command_may_run_unasked", _REAL_CAGED)
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: True)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)


def _gate(env, cmd, attended):
    asked: list = []

    def _cb(tool, args, preview):
        asked.append(preview)
        return False

    perms = KitToolPermissions(mode="default", workspace=str(env["ws"]),
                               confirm_callback=_cb if attended else None)
    admitted: list = []
    verdict = A._doc_executor._run_permission_gate(
        "bash", {"command": cmd}, perms, caged_admission=admitted)
    return verdict, asked, admitted


@pytest.mark.parametrize("attended", [True, False])
def test_control_the_cage_admits_the_command_without_a_policy(
        env, monkeypatch, attended):
    _cage(monkeypatch)
    verdict, asked, admitted = _gate(env, "./run.sh", attended)
    assert verdict is None and asked == [] and admitted == [True]


@pytest.mark.parametrize("attended", [True, False])
def test_a_managed_deny_refuses_at_the_caged_tier(env, monkeypatch, attended):
    _cage(monkeypatch)
    _managed(env, deny=["bash(./run.sh)"])
    verdict, asked, admitted = _gate(env, "./run.sh", attended)
    assert verdict and "managed policy" in verdict, verdict
    assert asked == [] and admitted == []
    kinds = [e.kind for e in security_events.recent(10)]
    assert "managed_policy_deny" in kinds


@pytest.mark.parametrize("attended", [True, False])
def test_a_managed_ask_keeps_the_command_off_the_cage(env, monkeypatch,
                                                      attended):
    _cage(monkeypatch)
    _managed(env, ask=["bash(./run.sh)"])
    verdict, asked, admitted = _gate(env, "./run.sh", attended)
    assert admitted == []
    assert verdict, "the command ran unasked"
    assert bool(asked) is attended


def test_only_allow_keeps_a_command_off_the_cage(env, monkeypatch):
    _cage(monkeypatch)
    _managed(env, only_allow=["bash(ls:*)"])
    verdict, asked, admitted = _gate(env, "./run.sh", True)
    assert admitted == [] and asked, (verdict, asked)
    assert not _REAL_CAGED("./run.sh", _perms(env), "default")


# ---------------------------------------------------------------------------
# Reported as data
# ---------------------------------------------------------------------------

def test_a_managed_refusal_is_a_permission_denied_event(env):
    _managed(env, deny=["bash(touch)"], ask=["bash(./run.sh)"])
    events: list = []
    with TO.TurnObserver(emit=events.append) as obs:
        _run(_perms(env, "bypassPermissions"), "bash",
             {"command": "touch x"})
        obs.saw_tool_use("bash")
        _run(_perms(env, "bypassPermissions"), "bash",
             {"command": "./run.sh"})
        summary = obs.summary(error="", engine_stop_reason="")
    denied = [e for e in events if e["type"] == "permission_denied"]
    codes = [e["reason_code"] for e in denied]
    assert codes == ["managed_policy_deny", "managed_policy_ask"], denied
    assert [e["needs_approval"] for e in denied] == [False, True]
    assert summary["stop_reason"] == "needs_approval"


def test_the_managed_reason_codes_are_known_and_labelled():
    labels = security_events.known_kinds()
    for code in ("managed_policy_deny", "managed_policy_ask"):
        assert code in TO.REASON_CODES and code in labels
    assert "managed_policy_ask" in TO.APPROVABLE_REASON_CODES
    assert "managed_policy_deny" not in TO.APPROVABLE_REASON_CODES


# ---------------------------------------------------------------------------
# Shown together, with sources
# ---------------------------------------------------------------------------

def test_capabilities_list_managed_and_user_rules_with_sources(env):
    from delfin.agent import boundaries
    _managed(env, deny=["bash(sudo)"], only_allow=["bash(ls:*)"])
    _user(env, deny=["bash(git push)"], allow=["bash(pytest:*)"])
    caps = boundaries.capabilities(cwd=env["ws"])
    rules = caps["permission_rules"]["rules"]
    got = {(r["list"], r["rule"], r["origin"], r["source"]) for r in rules}
    assert ("deny", "bash(sudo)", "managed", str(env["policy"])) in got
    assert ("only_allow", "bash(ls:*)", "managed", str(env["policy"])) in got
    assert ("deny", "bash(git push)", "user", str(env["settings"])) in got
    assert ("allow", "bash(pytest:*)", "user", str(env["settings"])) in got
    assert [r["origin"] for r in rules if r["list"] == "deny"] == [
        "managed", "user"]
    text = boundaries.format_capabilities(caps)
    assert f"bash(sudo)  [managed: {env['policy']}]" in text, text
    assert f"bash(git push)  [user: {env['settings']}]" in text, text
    assert text.count("bash(sudo)") == 1, text


def test_permissions_lists_managed_and_user_rules_with_sources(env):
    from delfin.agent import repl_commands as rc
    _managed(env, deny=["bash(sudo)"], only_allow=["bash(ls:*)"])
    _user(env, deny=["bash(git push)"])
    for line in ("/permissions", "/permissions rules"):
        out = rc.dispatch(line, rc.ReplContext(workspace=env["ws"])).output
        assert f"bash(sudo)  [managed {env['policy']}]" in out, out
        assert "bash(git push)  [user]" in out, out
        assert "only_allow" in out and "bash(ls:*)" in out, out
        assert str(env["settings"]) in out and str(env["policy"]) in out
