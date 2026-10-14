"""The administrator's managed policy through the real gate.

Input: a policy file at a path injected through
``managed_policy._PATHS_OVERRIDE`` (the system paths are never written),
with the admin uid set injected through ``_ADMIN_UIDS_OVERRIDE`` and the
isolation probe through ``_ISOLATION_PROBE``; tool calls through
``_DocToolExecutor.execute`` and ``cli_gate``-free entry points.

Pinned: a managed deny refuses in every mode, bypass included, before any
user allow or approval; a managed ask asks in every mode and refuses with
nobody to ask; ``only_allow`` bounds the bash auto-allow list; a disallowed
mode falls back (unattended) or is refused (attended, explicit, and from
/permissions); ``require_isolation`` forces the shell into a sandbox and
caps an unattended run on a host without one; the network ceiling
intersects the user's domains; the watcher cannot be switched off; a
user-owned or group-writable file only ever adds restrictions; the policy
file is unwritable through the write tools and the shell; the boundary
report and the audit destination carry the policy.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import managed_policy as MP
from delfin.agent import security_events
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

MODES = ("plan", "default", "acceptEdits", "bypassPermissions")


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
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "notes.txt").write_text("NOTES_BODY\n")
    etc = tmp_path / "etc" / "delfin"
    etc.mkdir(parents=True)
    etc.chmod(0o755)
    etc.parent.chmod(0o755)
    policy = etc / "managed_policy.json"
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(policy),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE",
                        frozenset({0, os.getuid()}))
    monkeypatch.setattr(MP, "_ISOLATION_PROBE", lambda backend="": True)
    MP.clear_cache()
    security_events.clear()
    yield {"ws": ws, "policy": policy, "tmp": tmp_path, "home": home}
    MP.clear_cache()


def _write_policy(path: Path, data, mode=0o644):
    path.write_text(data if isinstance(data, str) else json.dumps(data))
    path.chmod(mode)
    MP.clear_cache()


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
    return str(data.get("error") or "") if isinstance(data, dict) else ""


# ---------------------------------------------------------------------------
# No file: nothing changes
# ---------------------------------------------------------------------------

def test_without_a_policy_file_nothing_is_managed(env):
    pol = MP.load()
    assert not pol.active
    for mode in MODES:
        assert MP.coerce_mode(mode) == mode
    assert _perms(env, "bypassPermissions").mode == "bypassPermissions"


# ---------------------------------------------------------------------------
# Deny and ask
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
def test_a_managed_deny_refuses_a_write_in_every_mode(env, mode):
    _write_policy(env["policy"], {"permissions": {"deny": ["write(notes.txt)"]}})
    perms = _perms(env, mode, answer=True)
    perms.read_tracker[str((env["ws"] / "notes.txt").resolve())] = (
        (env["ws"] / "notes.txt").stat().st_mtime)
    err = _err(_run(perms, "write_file",
                    {"path": "notes.txt", "content": "CHANGED\n"}))
    assert "managed policy" in err
    assert (env["ws"] / "notes.txt").read_text() == "NOTES_BODY\n"
    assert perms.confirm_callback.asked == []


@pytest.mark.parametrize("mode", MODES)
def test_a_managed_deny_refuses_a_command_in_every_mode(env, mode):
    _write_policy(env["policy"], {"permissions": {"deny": ["bash(touch)"]}})
    perms = _perms(env, mode, answer=True)
    err = _err(_run(perms, "bash", {"command": "touch made.txt"}))
    assert "managed policy" in err
    assert not (env["ws"] / "made.txt").exists()


def test_a_user_allow_cannot_lift_a_managed_deny(env):
    _write_policy(env["policy"], {"permissions": {"deny": ["bash(touch)"]}})
    perms = _perms(env, "default", answer=True,
                   bash_auto_allow_patterns=(r"^touch\b",))
    err = _err(_run(perms, "bash", {"command": "touch made.txt"}))
    assert "managed policy" in err
    assert perms.confirm_callback.asked == []
    assert not (env["ws"] / "made.txt").exists()


def test_a_deny_matches_the_words_anywhere_in_a_chain(env):
    _write_policy(env["policy"], {"permissions": {"deny": ["bash(git push)"]}})
    perms = _perms(env, "bypassPermissions")
    for cmd in ("git push", "true && git push origin x",
                "sh -c 'git push'", "/usr/bin/git push"):
        assert "managed policy" in _err(_run(perms, "bash", {"command": cmd})), cmd


def test_a_read_deny_holds_for_the_shell_too(env):
    secret = env["ws"] / "secret"
    secret.mkdir()
    (secret / "k.txt").write_text("K\n")
    _write_policy(env["policy"], {"permissions": {"deny": ["read(secret)"]}})
    perms = _perms(env, "bypassPermissions")
    assert "managed policy" in _err(_run(perms, "read_file",
                                         {"path": "secret/k.txt"}))
    assert "managed policy" in _err(_run(perms, "bash",
                                         {"command": "cat secret/k.txt"}))


def test_a_tool_deny_names_the_tool(env):
    _write_policy(env["policy"], {"permissions": {"deny": ["web_search"]}})
    perms = _perms(env, "bypassPermissions")
    assert "managed policy" in _err(_run(perms, "web_search", {"query": "x"}))


@pytest.mark.parametrize("mode", ("default", "acceptEdits", "bypassPermissions"))
def test_a_managed_ask_asks_in_every_mode(env, mode):
    _write_policy(env["policy"], {"permissions": {"ask": ["bash(touch)"]}})
    perms = _perms(env, mode, answer=False)
    err = _err(_run(perms, "bash", {"command": "touch made.txt"}))
    assert err
    assert perms.confirm_callback.asked
    assert "MANAGED POLICY" in perms.confirm_callback.asked[0][2]
    assert not (env["ws"] / "made.txt").exists()


def test_a_managed_ask_with_nobody_to_ask_refuses(env):
    _write_policy(env["policy"], {"permissions": {"ask": ["bash(touch)"]}})
    perms = _perms(env, "bypassPermissions")
    err = _err(_run(perms, "bash", {"command": "touch made.txt"}))
    assert "managed policy" in err and "nobody" in err
    assert not (env["ws"] / "made.txt").exists()


def test_only_allow_bounds_the_auto_allow_list(env):
    _write_policy(env["policy"],
                  {"permissions": {"only_allow": ["bash(ls:*)"]}})
    perms = _perms(env, "default", bash_auto_allow_patterns=(r"^ls\b",
                                                             r"^cat\b"))
    assert perms.matches_bash_auto_allow("ls -la")
    assert not perms.matches_bash_auto_allow("cat notes.txt")
    assert not perms.matches_bash_auto_allow("ls; cat notes.txt")


def test_the_manual_example_is_a_policy_that_loads_clean(env):
    import re
    manual = (Path(__file__).resolve().parents[1] / "docs" / "USER_MANUAL.md"
              ).read_text(encoding="utf-8")
    # The user's permission-rule example comes earlier in the manual and
    # has the same opening; the managed one is in its own section.
    manual = manual.split("### For administrators: a managed policy", 1)[-1]
    block = re.search(r"```json\n(\{\n  \"permissions\": \{\n    \"deny\".*?)```",
                      manual, re.DOTALL)
    assert block, "the managed policy example is gone from the manual"
    _write_policy(env["policy"], block.group(1))
    pol = MP.load()
    assert pol.sources[0].trusted and not pol.errors, pol.errors
    assert len(pol.deny) == 3 and len(pol.ask) == 2
    assert pol.require_isolation and pol.disable_classifier_off
    assert MP.coerce_mode("bypassPermissions") == "default"


def test_a_broken_rule_is_reported(env):
    _write_policy(env["policy"],
                  {"permissions": {"deny": ["frob(x)", "bash(touch)"]}})
    pol = MP.load()
    assert pol.sources[0].trusted, pol.sources[0].note
    assert [r.text for r in pol.deny] == ["bash(touch)"]
    assert any("frob(x)" in e for e in pol.errors)


# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------

def test_a_disallowed_mode_falls_back_to_default(env):
    _write_policy(env["policy"], {"disallowed_modes": ["bypassPermissions"]})
    perms = _perms(env, "bypassPermissions")
    assert perms.mode == "default"
    perms.mode = "bypassPermissions"
    assert perms.mode == "default"
    perms.mode = "acceptEdits"
    assert perms.mode == "acceptEdits"
    # A mode adopted from inside a running call is the other route.
    with A._pin_permission_mode(perms):
        A.repin_permission_mode("bypassPermissions")
        assert A.effective_mode(perms) == "default"


def test_default_disallowed_falls_back_to_plan(env):
    _write_policy(env["policy"],
                  {"disallowed_modes": ["bypassPermissions", "default",
                                        "acceptEdits"]})
    assert MP.coerce_mode("bypassPermissions") == "plan"


def test_the_launch_decision(env):
    _write_policy(env["policy"], {"disallowed_modes": ["bypassPermissions"]})
    mode, note, refuse = MP.launch_decision(
        "bypassPermissions", explicit=True, attended=False)
    assert (mode, refuse) == ("default", False) and "managed" in note
    mode, note, refuse = MP.launch_decision(
        "bypassPermissions", explicit=True, attended=True)
    assert refuse and "bypassPermissions" in note
    # A persisted setting is not an explicit request: fall back, start.
    mode, note, refuse = MP.launch_decision(
        "bypassPermissions", explicit=False, attended=True)
    assert (mode, refuse) == ("default", False)
    assert MP.launch_decision("acceptEdits", explicit=True,
                              attended=True) == ("acceptEdits", "", False)


def test_permissions_command_refuses_a_disallowed_mode(env):
    from delfin.agent import repl_commands as RC
    _write_policy(env["policy"], {"disallowed_modes": ["acceptEdits"]})
    perms = _perms(env, "default")

    class Eng:
        kit_permissions = perms

        def set_kit_permission_mode(self, mode):
            perms.mode = mode
            return True

    out = RC._permissions(RC.ReplContext(engine=Eng()), "acceptEdits").output
    assert "managed" in out and perms.mode == "default"
    shown = RC._permissions(RC.ReplContext(engine=Eng()), "").output
    assert "managed policy" in shown and str(env["policy"]) in shown


# ---------------------------------------------------------------------------
# Isolation, network, watcher
# ---------------------------------------------------------------------------

def test_require_isolation_never_runs_a_plain_shell(env, monkeypatch):
    _write_policy(env["policy"], {"require_isolation": True})
    perms = _perms(env, "default")
    cmd = "echo hi"
    argv = A._bash_isolation_argv(cmd, env["ws"], perms, mode="off")
    assert argv != ["/bin/bash", "-c", cmd]
    assert argv != A._in_process_cage(["/bin/bash", "-c", cmd], env["ws"])


def test_require_isolation_caps_an_unattended_run_without_a_sandbox(env, monkeypatch):
    monkeypatch.setattr(MP, "_ISOLATION_PROBE", lambda backend="": False)
    _write_policy(env["policy"], {"require_isolation": True})
    assert MP.coerce_mode("bypassPermissions") == "default"
    assert MP.coerce_mode("acceptEdits") == "default"
    mode, note, refuse = MP.launch_decision("default", explicit=False,
                                            attended=True)
    assert refuse and "isolation" in note
    mode, note, refuse = MP.launch_decision("bypassPermissions",
                                            explicit=True, attended=False)
    assert (mode, refuse) == ("default", False) and "isolation" in note


def test_the_domain_ceiling_intersects_the_user_list(env, monkeypatch):
    from delfin.agent import egress_proxy as EP
    _write_policy(env["policy"],
                  {"network": {"allowed_domains_max": ["pypi.org",
                                                       "*.github.com"]}})
    proxy = EP.EgressProxy(allowed_domains=["pypi.org", "evil.example",
                                            "api.github.com"])
    try:
        assert proxy._check("pypi.org", 443) is None
        assert proxy._check("api.github.com", 443) is None
        assert "managed" in (proxy._check("evil.example", 443) or "")
    finally:
        proxy.close()
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "open"})
    assert EP.network_mode() == "proxy"


def test_web_fetch_outside_the_ceiling_is_refused(env):
    _write_policy(env["policy"],
                  {"network": {"allowed_domains_max": ["pypi.org"]}})
    perms = _perms(env, "bypassPermissions")
    err = _err(_run(perms, "web_fetch", {"url": "https://example.com/x"}))
    assert "managed policy" in err


def test_the_watcher_cannot_be_switched_off(env, monkeypatch):
    from delfin import user_settings as US
    monkeypatch.setattr(US, "load_settings", lambda *a, **k: {
        "agent": {"classifier": {"enabled": False}}})
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    assert A._classifier_settings() is None
    _write_policy(env["policy"], {"disable_classifier_off": True})
    assert A._classifier_settings() is not None


# ---------------------------------------------------------------------------
# Trust
# ---------------------------------------------------------------------------

def test_a_user_owned_file_is_not_trusted(env, monkeypatch):
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE", frozenset({0}))
    if os.getuid() == 0:
        pytest.skip("root owns every file here")
    _write_policy(env["policy"], {"permissions": {"deny": ["bash(touch)"]},
                                  "audit_log": str(env["tmp"] / "a.jsonl")})
    pol = MP.load()
    src = pol.sources[0]
    assert not src.trusted and "owned by" in src.note
    # It may still restrict ...
    assert [r.text for r in pol.deny] == ["bash(touch)"]
    # ... but nothing it names that only an administrator may set is used.
    assert MP.audit_destinations() == []


def test_a_group_writable_file_is_not_trusted(env):
    _write_policy(env["policy"], {"disallowed_modes": ["bypassPermissions"]},
                  mode=0o664)
    pol = MP.load()
    assert not pol.sources[0].trusted
    assert "writable" in pol.sources[0].note


def test_a_file_in_a_writable_directory_is_not_trusted(env):
    env["policy"].parent.chmod(0o777)
    try:
        _write_policy(env["policy"], {"disallowed_modes": []})
        assert not MP.load().sources[0].trusted
    finally:
        env["policy"].parent.chmod(0o755)


def test_an_untrusted_file_cannot_relax_a_trusted_one(env, monkeypatch):
    _write_policy(env["policy"], {"disallowed_modes": ["bypassPermissions"],
                                  "network": {"allowed_domains_max":
                                              ["pypi.org"]}})
    loose = env["tmp"] / "loose.json"
    loose.write_text(json.dumps({"disallowed_modes": [],
                                 "network": {"allowed_domains_max":
                                             ["pypi.org", "evil.example"]}}))
    monkeypatch.setenv(MP.ENV_VAR, str(loose))
    MP.clear_cache()
    pol = MP.load()
    assert len(pol.sources) == 2 and not pol.sources[1].trusted
    assert MP.coerce_mode("bypassPermissions") == "default"
    assert not MP.domain_within_max("evil.example")


def test_the_environment_variable_only_adds_restrictions(env, monkeypatch):
    extra = env["tmp"] / "extra.json"
    extra.write_text(json.dumps({"disallowed_modes": ["acceptEdits"],
                                 "audit_log": str(env["tmp"] / "x.jsonl")}))
    extra.chmod(0o644)
    monkeypatch.setenv(MP.ENV_VAR, str(extra))
    MP.clear_cache()
    pol = MP.load()
    assert not pol.sources[-1].trusted
    assert MP.coerce_mode("acceptEdits") == "default"
    assert MP.audit_destinations() == []


def test_a_broken_trusted_file_fails_closed(env):
    _write_policy(env["policy"], "{not json")
    pol = MP.load()
    assert pol.errors
    assert MP.coerce_mode("bypassPermissions") == "default"
    assert MP.coerce_mode("acceptEdits") == "default"


def test_user_settings_cannot_name_another_policy_file(env, monkeypatch):
    from delfin import user_settings as US
    other = env["tmp"] / "other.json"
    other.write_text("{}")
    monkeypatch.setattr(US, "load_settings", lambda *a, **k: {
        "agent": {"managed_policy": str(other),
                  "managed_policy_path": str(other)}})
    assert [str(p) for p, _ in MP.candidate_paths()] == [str(env["policy"])]


# ---------------------------------------------------------------------------
# The file itself
# ---------------------------------------------------------------------------

def test_the_policy_file_is_unwritable_through_every_route(env, monkeypatch):
    # The policy file placed INSIDE the workspace, so only the managed
    # guard -- not the workspace boundary -- can refuse the write.
    inside = env["ws"] / "managed_policy.json"
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(inside),))
    _write_policy(inside, {"disallowed_modes": ["bypassPermissions"]})
    before = inside.read_text()
    perms = _perms(env, "acceptEdits", answer=True)
    # Read first, so the read-before-write rule is not what refuses.
    assert not _err(_run(perms, "read_file", {"path": "managed_policy.json"}))
    assert "managed policy" in _err(_run(perms, "write_file", {
        "path": "managed_policy.json", "content": "{}"}))
    for cmd in ("echo '{}' > managed_policy.json",
                "cp /dev/null managed_policy.json"):
        assert "managed policy" in _err(_run(perms, "bash", {"command": cmd}))
    assert inside.read_text() == before


# ---------------------------------------------------------------------------
# Visibility and audit
# ---------------------------------------------------------------------------

def test_capabilities_name_the_policy(env):
    from delfin.agent import boundaries
    _write_policy(env["policy"], {"disallowed_modes": ["bypassPermissions"],
                                  "permissions": {"deny": ["bash(touch)"]}})
    caps = boundaries.capabilities(cwd=env["ws"],
                                   permission_mode="bypassPermissions")
    assert caps["mode"] == "default"
    mp = caps["managed_policy"]
    assert mp["active"] and mp["sources"][0]["path"] == str(env["policy"])
    text = boundaries.format_capabilities(caps)
    assert "managed" in text and str(env["policy"]) in text


def test_summary_lines_say_what_is_enforced(env):
    _write_policy(env["policy"], {"disallowed_modes": ["bypassPermissions"],
                                  "require_isolation": True,
                                  "disable_classifier_off": True,
                                  "permissions": {"deny": ["bash(touch)"]}})
    text = "\n".join(MP.summary_lines())
    for word in (str(env["policy"]), "bypassPermissions", "isolation",
                 "bash(touch)", "watcher"):
        assert word in text, word


def test_the_audit_record_reaches_the_managed_destination(env):
    from delfin.agent import audit_log
    dest = env["tmp"] / "audit" / "{user}.jsonl"
    dest.parent.mkdir()
    _write_policy(env["policy"], {"audit_log": str(dest)})
    [target] = MP.audit_destinations()
    audit_log.append({"tool": "bash", "command": "echo", "decision": "auto"},
                     log_path=env["home"] / "audit.log")
    lines = target.read_text().splitlines()
    assert json.loads(lines[-1])["command"] == "echo"
