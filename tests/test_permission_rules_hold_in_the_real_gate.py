"""Structured permission rules (``permissions.deny/ask/allow``) through the
real executor.

Input: rules in the user's settings file and in a project's
``.delfin/settings.json``; tool calls through ``_DocToolExecutor.execute``
(and ``_mcp_call_refusal`` for MCP names, ``cli_gate.decide`` for the CLI
backend). Pinned: each rule kind, deny > ask > allow, deny and ask in every
mode, allow skipping only the question for the exact call and never a hard
refusal, a project allow ignored, a chained or redirected command never
covered, the settings file unwritable through every write route, the
startup notice for a broken rule, and the read-only ``/permissions`` view.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from delfin import user_settings as US
from delfin.agent import api_client as A
from delfin.agent import permission_rules as PR
from delfin.agent import security_events
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

MODES = ("plan", "default", "acceptEdits", "bypassPermissions")


class Asker:
    """A confirm callback that records every question and answers one way."""

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
    (ws / "src").mkdir(parents=True)
    (ws / "src" / "a.py").write_text("x = 1\n")
    (ws / "notes.txt").write_text("NOTES_BODY\n")
    run = ws / "run.sh"
    run.write_text("#!/bin/sh\necho RUN_MARKER\n")
    run.chmod(0o755)
    far = tmp_path / "far"
    far.mkdir()
    (far / "data.txt").write_text("FAR_MARKER\n")
    PR._CACHE.clear()
    security_events.clear()
    return {"home": home, "settings": settings, "ws": ws, "far": far,
            "tmp": tmp_path}


def _rules(env, **lists):
    env["settings"].write_text(json.dumps({"permissions": lists}))
    PR._CACHE.clear()


def _project_rules(env, **lists):
    p = env["ws"] / ".delfin" / "settings.json"
    p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps({"permissions": lists}))
    PR._CACHE.clear()


def _perms(env, mode="default", answer=None, **over):
    perms = KitToolPermissions(workspace=env["ws"], mode=mode, **over)
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
    if isinstance(data, dict):
        return str(data.get("error") or "")
    return ""


# ---------------------------------------------------------------------------
# Parsing and loading
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text,kind", [
    ("bash(git push)", "bash"), ("bash(pytest:*)", "bash"),
    ("read(/data/**)", "read"), ("write(src/**)", "write"),
    ("web_fetch(domain:example.org)", "web_fetch"),
    ("mcp__srv__tool", "tool"), ("mcp__srv__*", "tool"), ("mcp__srv", "tool"),
    ("web_search", "tool"),
])
def test_every_rule_kind_parses(text, kind):
    rule = PR.parse_rule(text)
    assert isinstance(rule, PR.Rule), rule
    assert rule.kind == kind


@pytest.mark.parametrize("text", [
    "", "bash()", "bash(a && b)", "bash(echo > f)", "frob(x)",
    "web_fetch(example.org)", "web_fetch(domain:exa mple)", "mcp__",
    "two words", 7,
])
def test_a_broken_rule_is_no_rule(text):
    assert isinstance(PR.parse_rule(text), str)


def test_invalid_rules_are_reported_and_grant_nothing(env):
    _rules(env, deny=["frob(x)", "bash(git push)"], allow=["bash(a ; b)"],
           ask="not a list")
    rs = PR.load(env["ws"])
    assert [r.text for r in rs.deny] == ["bash(git push)"]
    assert rs.allow == () and rs.ask == ()
    lines = PR.error_lines(rs)
    assert any(line.startswith("ERROR") and "frob(x)" in line for line in lines)
    assert any("bash(a ; b)" in line for line in lines)
    assert any("'ask' is not a list" in line for line in lines)


def test_an_unreadable_settings_file_grants_nothing(env):
    env["settings"].write_text("{not json")
    PR._CACHE.clear()
    rs = PR.load(env["ws"])
    assert rs.is_empty()
    assert any("not valid JSON" in line for line in PR.error_lines(rs))


def test_a_project_allow_is_ignored_and_said(env):
    _project_rules(env, allow=["bash(./run.sh:*)"], deny=["web_search"])
    rs = PR.load(env["ws"])
    assert rs.allow == ()
    assert [r.text for r in rs.deny] == ["web_search"]
    assert any("cannot grant" in line for line in PR.error_lines(rs))
    perms = _perms(env, answer=False)
    _run(perms, "bash", {"command": "./run.sh"})
    assert perms.confirm_callback.asked, "a project allow skipped the question"


# ---------------------------------------------------------------------------
# deny: every kind, every mode
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", [
    "git push origin topic", "echo hi && git push", "sh -c 'git push'",
    "timeout 5 git --no-pager push", "x=$(git push)",
])
def test_a_bash_deny_holds_in_every_mode_and_every_spelling(env, mode, cmd):
    _rules(env, deny=["bash(git push)"])
    perms = _perms(env, mode=mode, answer=True)
    err = _err(_run(perms, "bash", {"command": cmd}))
    assert "bash(git push)" in err, err
    assert not perms.confirm_callback.asked


@pytest.mark.parametrize("mode", MODES)
def test_a_read_deny_holds_for_the_tool_and_the_shell(env, mode):
    _rules(env, deny=["read(notes.txt)"])
    perms = _perms(env, mode=mode, answer=True)
    r = _run(perms, "read_file", {"path": "notes.txt"})
    assert "read(notes.txt)" in _err(r) and "NOTES_BODY" not in r, r
    r = _run(perms, "bash", {"command": "cat notes.txt"})
    assert "read(notes.txt)" in _err(r) and "NOTES_BODY" not in r, r


@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
def test_a_write_deny_holds_for_the_tool_and_the_shell(env, mode):
    _rules(env, deny=["write(src/**)"])
    perms = _perms(env, mode=mode, answer=True)
    _run(perms, "write_file", {"path": "src/new.py", "content": "y\n"})
    _run(perms, "bash", {"command": "echo y > src/new2.py"})
    _run(perms, "edit_file", {"path": "src/a.py", "old_string": "x = 1",
                              "new_string": "x = 2"})
    assert not (env["ws"] / "src" / "new.py").exists()
    assert not (env["ws"] / "src" / "new2.py").exists()
    assert (env["ws"] / "src" / "a.py").read_text() == "x = 1\n"
    # A write rule does not stop reading.
    assert "x = 1" in _run(perms, "read_file", {"path": "src/a.py"})


def test_a_read_deny_reaches_a_path_the_tool_resolves_itself(env):
    """grep over the workspace resolves the file; the deny reaches it."""
    _rules(env, deny=["read(notes.txt)"])
    perms = _perms(env, mode="bypassPermissions")
    ex = _DocToolExecutor()
    assert ex._check_read_access(perms, env["ws"] / "notes.txt") is not None
    assert ex._check_read_access(perms, env["ws"] / "src" / "a.py") is None


@pytest.mark.parametrize("mode", MODES)
def test_a_web_fetch_domain_deny_covers_subdomains(env, mode):
    _rules(env, deny=["web_fetch(domain:example.org)"])
    perms = _perms(env, mode=mode, answer=True)
    err = _err(_run(perms, "web_fetch", {"url": "https://docs.example.org/x"}))
    assert "web_fetch(domain:example.org)" in err, err


@pytest.mark.parametrize("mode", MODES)
def test_a_bare_tool_deny_holds_in_every_mode(env, mode):
    _rules(env, deny=["web_search"])
    perms = _perms(env, mode=mode, answer=True)
    err = _err(_run(perms, "web_search", {"query": "x"}))
    assert "web_search" in err and "permission rule" in err, err


@pytest.mark.parametrize("rule", ["mcp__srv__do", "mcp__srv__*", "mcp__srv"])
def test_an_mcp_deny_holds_on_the_mcp_route(env, rule):
    _rules(env, deny=[rule])
    perms = _perms(env, mode="bypassPermissions", answer=True)
    err = _err(A._mcp_call_refusal("mcp__srv__do", {}, perms) or "{}")
    assert rule in err, err
    assert A._mcp_call_refusal("mcp__other__do", {}, perms) is None


def test_a_deny_is_recorded_as_a_labelled_event(env):
    _rules(env, deny=["web_search"])
    _run(_perms(env, mode="bypassPermissions"), "web_search", {"query": "x"})
    kinds = {e.kind for e in security_events.recent(200)}
    assert "permission_rule_deny" in kinds
    assert "permission_rule_deny" in security_events.known_kinds()


# ---------------------------------------------------------------------------
# ask: adds a question in every mode; nobody to ask is a refusal
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
def test_an_ask_rule_asks_in_every_mode(env, mode):
    _rules(env, ask=["read(notes.txt)"])
    perms = _perms(env, mode=mode, answer=True)
    r = _run(perms, "read_file", {"path": "notes.txt"})
    assert "NOTES_BODY" in r, r
    assert len(perms.confirm_callback.asked) == 1
    assert "read(notes.txt)" in perms.confirm_callback.asked[0][2]


@pytest.mark.parametrize("mode", MODES)
def test_an_ask_declined_or_unanswerable_refuses(env, mode):
    _rules(env, ask=["read(notes.txt)"])
    perms = _perms(env, mode=mode, answer=False)
    r = _run(perms, "read_file", {"path": "notes.txt"})
    assert _err(r) and "NOTES_BODY" not in r, r
    r = _run(_perms(env, mode=mode), "read_file", {"path": "notes.txt"})
    assert "nobody can be asked" in _err(r), r


def test_an_approved_ask_is_not_asked_twice(env):
    """default mode would ask about ./run.sh anyway: one question, not two."""
    _rules(env, ask=["bash(./run.sh)"])
    perms = _perms(env, answer=True)
    r = _run(perms, "bash", {"command": "./run.sh"})
    assert "RUN_MARKER" in r, r
    assert len(perms.confirm_callback.asked) == 1


def test_an_approved_ask_does_not_lift_a_hard_refusal(env):
    _rules(env, ask=["bash(rm:*)"])
    perms = _perms(env, answer=True)
    (env["ws"] / "build").mkdir()
    err = _err(_run(perms, "bash", {"command": "rm -rf build"}))
    assert "deny-pattern" in err, err
    assert (env["ws"] / "build").exists()


# ---------------------------------------------------------------------------
# precedence
# ---------------------------------------------------------------------------

def test_deny_wins_over_ask_and_allow(env):
    _rules(env, deny=["bash(./run.sh)"], ask=["bash(./run.sh)"],
           allow=["bash(./run.sh:*)"])
    perms = _perms(env, answer=True)
    assert "bash(./run.sh)" in _err(_run(perms, "bash", {"command": "./run.sh"}))
    assert not perms.confirm_callback.asked


def test_ask_wins_over_allow(env):
    _rules(env, ask=["bash(./run.sh)"], allow=["bash(./run.sh:*)"])
    perms = _perms(env, mode="bypassPermissions", answer=False)
    assert _err(_run(perms, "bash", {"command": "./run.sh"}))
    assert len(perms.confirm_callback.asked) == 1


# ---------------------------------------------------------------------------
# allow: skips the question for exactly that call
# ---------------------------------------------------------------------------

def test_a_bash_allow_skips_the_question(env):
    _rules(env, allow=["bash(./run.sh:*)"])
    perms = _perms(env, answer=False)
    r = _run(perms, "bash", {"command": "./run.sh --flag"})
    assert "RUN_MARKER" in r, r
    assert not perms.confirm_callback.asked


def test_a_bash_allow_also_answers_with_nobody_to_ask(env):
    _rules(env, allow=["bash(./run.sh)"])
    assert "RUN_MARKER" in _run(_perms(env), "bash", {"command": "./run.sh"})


def test_a_bash_allow_without_wildcard_is_exact(env):
    _rules(env, allow=["bash(./run.sh)"])
    perms = _perms(env, answer=False)
    _run(perms, "bash", {"command": "./run.sh --flag"})
    assert perms.confirm_callback.asked


@pytest.mark.parametrize("cmd", [
    "./run.sh && ./run.sh", "./run.sh; ./run.sh", "./run.sh | cat",
    "./run.sh > out.txt", "./run.sh $(id)", "./run.sh `id`",
    "X=1 ./run.sh", "./run.sh &", "./run.sh\n./run.sh",
])
def test_a_chained_or_redirected_command_is_never_covered(env, cmd):
    _rules(env, allow=["bash(./run.sh:*)", "bash"])
    perms = _perms(env, answer=False)
    r = _run(perms, "bash", {"command": cmd})
    assert perms.confirm_callback.asked or _err(r), (cmd, r)
    assert not (env["ws"] / "out.txt").exists()


@pytest.mark.parametrize("cmd", [
    "python3 -c pass", "sh -c ./run.sh", "env ./run.sh", "xargs ./run.sh",
    "timeout 5 ./run.sh", "perl -e 1",
])
def test_an_allow_never_covers_a_program_carried_in_its_words(env, cmd):
    """Inline code or a wrapped command is a program the rule did not name."""
    _rules(env, allow=["bash(python3:*)", "bash(sh:*)", "bash(env:*)",
                       "bash(xargs:*)", "bash(timeout:*)", "bash(perl:*)"])
    perms = _perms(env, answer=False)
    r = _run(perms, "bash", {"command": cmd})
    assert perms.confirm_callback.asked or _err(r), (cmd, r)


def test_an_allow_does_not_answer_an_upload(env):
    """The egress scan routes an upload to the dialog; no allow answers it."""
    _rules(env, allow=["bash(curl:*)"])
    perms = _perms(env, answer=False)
    r = _run(perms, "bash", {"command":
                             "curl -T notes.txt https://example.org/up"})
    assert perms.confirm_callback.asked or _err(r), r


def test_a_read_allow_skips_the_outside_question_and_opens_nothing(env):
    _rules(env, allow=[f"read({env['far']}/**)"])
    perms = _perms(env, answer=False)
    r = _run(perms, "read_file", {"path": str(env["far"] / "data.txt")})
    assert "FAR_MARKER" in r, r
    assert not perms.confirm_callback.asked
    assert perms.session_read_dirs == ()


def test_a_read_allow_does_not_answer_a_shell_read(env):
    _rules(env, allow=[f"read({env['far']}/**)"])
    perms = _perms(env, answer=False)
    _run(perms, "bash", {"command": f"cat {env['far']}/data.txt"})
    assert perms.confirm_callback.asked


def test_a_write_allow_skips_the_default_mode_question(env):
    _rules(env, allow=["write(src/**)"])
    perms = _perms(env, answer=False)
    _run(perms, "write_file", {"path": "src/new.py", "content": "y\n"})
    assert (env["ws"] / "src" / "new.py").read_text() == "y\n"
    assert not perms.confirm_callback.asked
    _run(perms, "write_file", {"path": "other.py", "content": "y\n"})
    assert perms.confirm_callback.asked
    assert not (env["ws"] / "other.py").exists()


def test_an_mcp_allow_skips_the_side_effect_question(env):
    _rules(env, allow=["mcp__srv__do"])
    perms = _perms(env, answer=False)
    assert A._mcp_call_refusal("mcp__srv__do", {}, perms) is None
    assert not perms.confirm_callback.asked
    assert A._mcp_call_refusal("mcp__srv__other", {}, perms) is not None
    assert perms.confirm_callback.asked


# ---------------------------------------------------------------------------
# allow cannot lift a hard refusal
# ---------------------------------------------------------------------------

def _git_repo(env):
    import subprocess
    from conftest import child_env
    e = child_env(env["tmp"])
    e["GIT_CONFIG_NOSYSTEM"] = "1"
    for args in (["init", "-q"], ["remote", "add", "origin",
                                  "git@github.com:me/repo.git"]):
        subprocess.run(["git", *args], cwd=env["ws"], env=e, check=True,
                       capture_output=True)


HARD = [
    # (allow rules, mode, tool, args, what must not have happened)
    ("secret", ["read(**)"], "default", "read_file", {"path": ".env"}),
    ("outside write", None, "acceptEdits", "write_file",
     {"path": "{far}/new.txt", "content": "x"}),
    ("locked scope", None, "default", "read_file", {"path": "{far}/data.txt"}),
    ("protected config", ["write(**)"], "default", "write_file",
     {"path": ".delfin/settings.json", "content": "{}"}),
    ("self-modification", ["write(**)"], "default", "write_file",
     {"path": "delfin/agent/api_client.py", "content": "x"}),
    ("deny pattern", ["bash(rm:*)"], "default", "bash",
     {"command": "rm -rf src"}),
    ("git config program", ["bash(git config:*)"], "default", "bash",
     {"command": "git config core.hooksPath hooks"}),
    ("unrequested push", ["bash(git push:*)"], "default", "bash",
     {"command": "git push origin topic"}),
    ("refused earlier", ["bash(./run.sh:*)"], "default", "bash",
     {"command": "./run.sh"}),
]


@pytest.mark.parametrize("case", HARD, ids=[c[0] for c in HARD])
@pytest.mark.parametrize("attended", [False, True])
def test_an_allow_cannot_lift_a_hard_refusal(env, case, attended):
    label, allow, mode, tool, args = case
    far = env["far"]
    allow = allow or [f"read({far}/**)", f"write({far}/**)"]
    _rules(env, allow=allow)
    (env["ws"] / ".env").write_text("SECRET_MARKER=1\n")
    (env["ws"] / "delfin" / "agent").mkdir(parents=True)
    (env["ws"] / "delfin" / "agent" / "api_client.py").write_text("orig\n")
    if label in ("git config program", "unrequested push"):
        _git_repo(env)
    args = {k: (v.replace("{far}", str(far)) if isinstance(v, str) else v)
            for k, v in args.items()}
    perms = _perms(env, mode=mode, answer=False if attended else None,
                   lock_workspace=(label == "locked scope"))
    if label == "refused earlier":
        perms.record_denied_action("bash", args["command"])
    r = _run(perms, tool, args)
    asked = perms.confirm_callback.asked if attended else []
    # Refused, or put to the person (who said no): never let through unseen.
    assert _err(r) or asked, (label, r)
    assert "SECRET_MARKER" not in r and "FAR_MARKER" not in r and "RUN_MARKER" not in r
    assert not (far / "new.txt").exists()
    assert (env["ws"] / "src").is_dir()
    assert not (env["ws"] / ".delfin" / "settings.json").exists()
    assert (env["ws"] / "delfin" / "agent" / "api_client.py").read_text() == "orig\n"


@pytest.mark.parametrize("label,tool,args", [
    ("protected config", "write_file",
     {"path": ".delfin/settings.local.json", "content": "{}"}),
    ("self-modification", "write_file",
     {"path": "delfin/agent/api_client.py", "content": "x"}),
])
def test_an_allow_leaves_the_protected_question_in_place(env, label, tool, args):
    _rules(env, allow=["write(**)"])
    (env["ws"] / "delfin" / "agent").mkdir(parents=True)
    (env["ws"] / "delfin" / "agent" / "api_client.py").write_text("orig\n")
    perms = _perms(env, mode="acceptEdits", answer=False)
    _run(perms, tool, args)
    assert perms.confirm_callback.asked, label


# ---------------------------------------------------------------------------
# The settings file stays the user's
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
def test_no_write_route_reaches_the_settings_file(env, mode):
    target = env["settings"]
    target.write_text(json.dumps({"permissions": {"allow": [
        "write(~/**)", "bash(tee:*)", "bash(cp:*)", "bash(sed:*)",
        "bash(python3:*)", "edit_file", "write_file", "apply_patch"]}}))
    PR._CACHE.clear()
    before = target.read_text()
    planted = json.dumps({"permissions": {"allow": ["bash"]}})
    (env["ws"] / "planted.json").write_text(planted)
    perms = _perms(env, mode=mode, answer=True)
    name = US.SETTINGS_FILE_NAME
    attempts = [
        ("write_file", {"path": str(target), "content": planted}),
        ("edit_file", {"path": str(target), "old_string": "write",
                       "new_string": "bash"}),
        ("multi_edit", {"path": str(target), "edits": [
            {"old_string": "write", "new_string": "bash"}]}),
        ("apply_patch", {"diff": (
            f"--- {target}\n+++ {target}\n@@ -1 +1 @@\n-{before}\n+{planted}\n")}),
        ("bash", {"command": f"echo '{planted}' > ~/{name}"}),
        ("bash", {"command": f"tee ~/{name} < planted.json"}),
        ("bash", {"command": f"cp planted.json {target}"}),
        ("bash", {"command": f"sed -i s/write/bash/ {target}"}),
        ("bash", {"command": f"cd ~ && cp {env['ws']}/planted.json {name}"}),
        ("bash", {"command": (
            "python3 -c \"import pathlib; pathlib.Path.home().joinpath("
            f"'{name}').write_text('{{}}')\"")}),
    ]
    granted = json.loads(before)["permissions"]
    for tool, args in attempts:
        _run(perms, tool, args)
        # The gate's own settings read may merge shipped defaults into the
        # file (user_settings.load_settings); the rules must be untouched.
        assert json.loads(target.read_text())["permissions"] == granted, (
            tool, args)


def test_the_settings_file_name_is_a_protected_path(env):
    perms = _perms(env)
    assert perms.matches_path_protected(US.SETTINGS_FILE_NAME)
    assert perms.matches_path_protected(f"sub/{US.SETTINGS_FILE_NAME}")


# ---------------------------------------------------------------------------
# The CLI backend's hook asks the same rules
# ---------------------------------------------------------------------------

def test_the_cli_backend_hook_applies_deny_and_ask(env, monkeypatch):
    from delfin.agent import cli_gate, kit_settings
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        env["home"] / ".delfin" / "settings.json")
    _rules(env, deny=["bash(git push)"], ask=["read(notes.txt)"])
    token = cli_gate.new_token()
    cli_gate.write_state(token, workspace=str(env["ws"]),
                         mode="bypassPermissions", can_ask=False)
    ev = {"tool_name": "Bash", "tool_input": {"command": "git push"},
          "cwd": str(env["ws"])}
    assert "bash(git push)" in (cli_gate.decide(token, ev) or "")
    ev = {"tool_name": "Read",
          "tool_input": {"file_path": str(env["ws"] / "notes.txt")},
          "cwd": str(env["ws"])}
    assert "nobody can be asked" in (cli_gate.decide(token, ev) or "")


# ---------------------------------------------------------------------------
# What the user sees
# ---------------------------------------------------------------------------

def test_permissions_shows_the_rules_and_the_file(env):
    from delfin.agent import repl_commands as rc
    _rules(env, deny=["bash(git push)", "frob(x)"], allow=["bash(pytest:*)"])
    _project_rules(env, ask=["web_search"], allow=["bash"])
    out = rc.dispatch("/permissions", rc.ReplContext(workspace=env["ws"])).output
    for text in ("bash(git push)", "bash(pytest:*)", "web_search",
                 str(env["settings"]), "frob(x)", "cannot grant"):
        assert text in out, (text, out)
    assert "approval is" in out
    # Read-only: the view writes nothing.
    assert "frob(x)" in env["settings"].read_text()


def test_permissions_rules_is_the_view_alone(env):
    from delfin.agent import repl_commands as rc
    _rules(env, deny=["bash(git push)"])
    out = rc.dispatch("/permissions rules",
                      rc.ReplContext(workspace=env["ws"])).output
    assert "bash(git push)" in out and "could not switch" not in out, out


def test_a_broken_rule_is_named_at_startup(env):
    import argparse
    from delfin.agent import cli as agent_cli
    _rules(env, deny=["frob(x)"], allow=["bash(a;b)"])
    engine = type("E", (), {"client": type("C", (), {"model": "m"})(),
                            "provider": "kit", "backend": "api",
                            "workspace": env["ws"]})()
    args = argparse.Namespace(
        backend="", provider="", model="", mode="", cwd=str(env["ws"]),
        effort="", permission_mode="", settings_defaults=False,
        allowed_tools="", disallowed_tools="")
    notes = "\n".join(agent_cli._bounding_notices(args, engine))
    assert "frob(x)" in notes and "ERROR" in notes, notes
    assert "bash(a;b)" in notes, notes


def test_no_rules_no_startup_line(env):
    import argparse
    from delfin.agent import cli as agent_cli
    engine = type("E", (), {"client": type("C", (), {"model": "m"})(),
                            "provider": "kit", "backend": "api"})()
    args = argparse.Namespace(
        backend="", provider="", model="", mode="", cwd=str(env["ws"]),
        effort="", permission_mode="", settings_defaults=False,
        allowed_tools="", disallowed_tools="")
    assert agent_cli._bounding_notices(args, engine) == []


def test_no_rules_changes_nothing(env):
    """Without a permissions block the gate decides exactly as before."""
    perms = _perms(env, answer=False)
    _run(perms, "bash", {"command": "./run.sh"})
    assert len(perms.confirm_callback.asked) == 1
    assert os.path.exists(env["ws"] / "notes.txt")
