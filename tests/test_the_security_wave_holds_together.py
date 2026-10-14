"""Defensive review of the hardening wave (PRs #198-#220): places where a
protection reported success without holding, or where two protections met so
that one lifted the other.

Each defect has one failing-then-fixed test, named after it. Interactions
that already hold are pinned too (``test_held_*``), so the review is visible:
what was checked and what proved it.

Everything runs through the gate functions and the executor -- command
strings to ``_run_permission_gate`` / ``_gate_mcp_tool`` and the redaction
and client builders -- under a temporary HOME. No attack command runs
against the real system; the real ``~/.delfin`` is never touched.
"""

from __future__ import annotations

import json
import os
import time

import pytest

from delfin.agent import api_client as A
from delfin.agent import output_guard as OG
from delfin.agent.api_client import KitToolPermissions, _doc_executor


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    (h / ".delfin").mkdir(parents=True)
    monkeypatch.setenv("HOME", str(h))
    # kit_settings.USER_SETTINGS_PATH is expanded at import, so it still
    # names the real home; point it at this temporary one.
    from delfin.agent import kit_settings as _ks
    monkeypatch.setattr(_ks, "USER_SETTINGS_PATH", h / ".delfin" / "settings.json")
    return h


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    return d


# ===========================================================================
# FINDING 1 -- a direct-built engine inherited a persisted bypassPermissions
# ===========================================================================
#
# launch_guard.resolve_posture and cli._headless_posture refuse to inherit a
# settings default of bypassPermissions: the threat model promises "every
# session ... never inherits it from your settings". But a caller that builds
# an AgentEngine WITHOUT going through those CLI doors -- the scheduler
# daemon's _default_engine_factory, the benchmark engine -- passed no
# permission_mode, and create_client then adopted the persisted bypass
# verbatim. The unattended scheduler started with more rights than the
# terminal and the headless `run`/`chat -p` doors with the same (no) flags.

def _client_mode(provider, home, ws, permission_mode=""):
    (home / ".delfin" / "settings.json").write_text(
        json.dumps({"kit": {"default_mode": "bypassPermissions"}}))
    c = A.create_client(backend="api", provider=provider, api_key="x",
                        model="m", permission_mode=permission_mode,
                        cwd=str(ws))
    return getattr(c._permissions, "mode", None)


@pytest.mark.parametrize("provider", ["kit", "ollama"])
def test_a_persisted_bypass_is_not_inherited_by_a_direct_client(
        provider, home, ws):
    """The defect: with no permission_mode passed, create_client adopted a
    persisted bypassPermissions. It must fall back to default instead; only
    an explicit argument reaches bypass."""
    assert _client_mode(provider, home, ws) == "default"


@pytest.mark.parametrize("provider", ["kit", "ollama"])
def test_an_explicit_bypass_argument_still_reaches_bypass(provider, home, ws):
    """The fix must not cost the one legitimate route: a caller that passes
    permission_mode=bypassPermissions (the CLI doors do, after their
    --unattended check) still gets it."""
    assert _client_mode(provider, home, ws,
                        permission_mode="bypassPermissions") == "bypassPermissions"


@pytest.mark.parametrize("provider", ["kit", "ollama"])
@pytest.mark.parametrize("persisted", ["plan", "default", "acceptEdits"])
def test_a_normal_persisted_default_is_still_inherited(provider, persisted,
                                                       home, ws):
    """Security must not cost function: every non-bypass persisted default
    is inherited exactly as before."""
    (home / ".delfin" / "settings.json").write_text(
        json.dumps({"kit": {"default_mode": persisted}}))
    c = A.create_client(backend="api", provider=provider, api_key="x",
                        model="m", permission_mode="", cwd=str(ws))
    assert getattr(c._permissions, "mode", None) == persisted


# The scheduler's own opt-in. After the fix above the persisted default no
# longer reaches the scheduler, so it gets an explicit, user-only switch
# mirroring `--permission-mode bypassPermissions --unattended`:
# agent.scheduler.permission_mode + agent.scheduler.unattended in the user
# settings file. The managed policy's disallowed_modes still wins.

_AGENT = {"provider": "ollama", "backend": "api", "model": "m"}


@pytest.fixture
def scheduler(home, monkeypatch):
    from delfin.agent import isolation_posture as IP
    from delfin.agent import managed_policy as MP
    from delfin.agent import scheduler_daemon as SD
    # No key resolution and no real endpoint: the factory only builds the
    # client whose mode is read.
    monkeypatch.setattr(SD._jm, "_resolve_provider_and_key",
                        lambda model, prov: ("ollama", "x"))
    # The isolation start check is pinned to "isolated" so the mode tests
    # do not depend on this host; its own test below unpins it.
    monkeypatch.setattr(IP, "isolation_status",
                        lambda perms=None: IP.IsolationStatus(True, layer="bwrap"))
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", ())
    MP.clear_cache()
    yield SD
    MP.clear_cache()


def _scheduled_mode(SD, ws, settings):
    engine = SD._default_engine_factory(str(ws), settings)
    return getattr(engine.client._permissions, "mode", None)


def test_the_scheduler_factory_does_not_inherit_bypass(home, ws, scheduler):
    """The real unattended entry: a persisted kit.default_mode of bypass
    and no scheduler opt-in -> default."""
    (home / ".delfin" / "settings.json").write_text(
        json.dumps({"kit": {"default_mode": "bypassPermissions"}}))
    assert _scheduled_mode(scheduler, ws, {"agent": dict(_AGENT)}) == "default"


def test_the_scheduler_opt_in_starts_a_scheduled_run_in_bypass(ws, scheduler):
    settings = {"agent": dict(_AGENT, scheduler={
        "permission_mode": "bypassPermissions", "unattended": True})}
    assert _scheduled_mode(scheduler, ws, settings) == "bypassPermissions"


def test_the_scheduler_opt_in_is_read_from_the_user_settings_file(
        ws, scheduler):
    """settings=None: the factory reads the user settings file itself."""
    from delfin import user_settings as US
    US.get_settings_path().write_text(json.dumps({"agent": dict(
        _AGENT, scheduler={"permission_mode": "bypassPermissions",
                           "unattended": True})}))
    assert _scheduled_mode(scheduler, ws, None) == "bypassPermissions"


@pytest.mark.parametrize("cfg", [
    {"permission_mode": "bypassPermissions"},
    {"permission_mode": "bypassPermissions", "unattended": False},
    {"permission_mode": "bypassPermissions", "unattended": "yes"},
    {"unattended": True},
    {"permission_mode": "bypass", "unattended": True},
    {},
])
def test_without_both_halves_of_the_opt_in_a_scheduled_run_is_default(
        cfg, ws, scheduler):
    settings = {"agent": dict(_AGENT, scheduler=cfg)}
    assert _scheduled_mode(scheduler, ws, settings) == "default"


@pytest.mark.parametrize("mode", ["plan", "acceptEdits"])
def test_the_scheduler_honours_a_named_non_bypass_mode(mode, ws, scheduler):
    settings = {"agent": dict(_AGENT, scheduler={"permission_mode": mode})}
    assert _scheduled_mode(scheduler, ws, settings) == mode


def test_a_managed_policy_disallowing_bypass_wins_over_the_opt_in(
        tmp_path, ws, scheduler, monkeypatch):
    from delfin.agent import managed_policy as MP
    etc = tmp_path / "etc" / "delfin"
    etc.mkdir(parents=True)
    for d in (etc, etc.parent):
        d.chmod(0o755)
    policy = etc / "managed_policy.json"
    policy.write_text(json.dumps({"disallowed_modes": ["bypassPermissions"]}))
    policy.chmod(0o644)
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(policy),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE", frozenset({0, os.getuid()}))
    MP.clear_cache()
    settings = {"agent": dict(_AGENT, scheduler={
        "permission_mode": "bypassPermissions", "unattended": True})}
    assert _scheduled_mode(scheduler, ws, settings) == "default"


def test_the_opt_in_is_lowered_where_nothing_isolates_a_command(
        ws, scheduler, monkeypatch):
    """The `delfin-agent run` start check: bypass with no isolation layer
    and nobody to confirm it runs in default."""
    from delfin.agent import isolation_posture as IP
    monkeypatch.setattr(IP, "isolation_status", lambda perms=None:
                        IP.IsolationStatus(False, reason="no layer here"))
    monkeypatch.setattr(IP, "_BYPASS_CONFIRMED", False)
    settings = {"agent": dict(_AGENT, scheduler={
        "permission_mode": "bypassPermissions", "unattended": True})}
    assert _scheduled_mode(scheduler, ws, settings) == "default"


# ===========================================================================
# FINDING 2 -- a large tool result skipped secret redaction entirely
# ===========================================================================
#
# _redact_tool_result returned any text over 400_000 chars unchanged -- a
# perf guard against the super-linear pattern scan. But it also skipped the
# exact-value scrub, so a credential DELFIN holds, printed in a big result (a
# read of a large log, a long stack trace), reached the model AND the tool
# trace on disk (which a bug report bundles) verbatim. The exact scrub is a
# few str.replace calls; only the pattern scan needs the cap.

SECRET = "kit-NOSHAPE-SECRET-abcdef0123456789xyz"


@pytest.fixture
def holds_secret(monkeypatch):
    monkeypatch.setattr("delfin.agent.web_tools._known_secret_values",
                        lambda: [SECRET])


def test_a_large_tool_result_still_has_the_held_credential_removed(
        holds_secret):
    """The defect: a >400 KB result carried the exact credential through."""
    big = ("x" * 500_000) + f" TOKEN={SECRET} rest"
    out = A._redact_tool_result(big)
    assert SECRET not in out
    assert "[redacted:credential]" in out


def test_a_small_tool_result_is_still_fully_scrubbed(holds_secret):
    """Below the cap, both the exact value and the pattern scan run, as
    before -- the exact value here, and a shaped token too."""
    small = f"log TOKEN={SECRET} and AKIA" + "A" * 16
    out = A._redact_tool_result(small)
    assert SECRET not in out


def test_the_exact_scrub_on_a_large_text_is_cheap():
    """Why only the pattern scan is capped: the exact pass is bounded by the
    text length, not super-linear. A guard against re-introducing the regex
    above the cap."""
    big = "a" * 4_000_000
    t0 = time.monotonic()
    OG.scrub_secrets(big, exact_only=True)
    assert time.monotonic() - t0 < 1.0


def test_scrub_secrets_exact_only_skips_the_pattern_scan(monkeypatch):
    """exact_only must not run the shape scan -- that is the whole point of
    the flag (it is what the size cap guards against)."""
    called = {"n": 0}

    def _spy(text, findings):
        called["n"] += 1
        return text

    monkeypatch.setattr(OG, "_redact_secrets", _spy)
    OG.scrub_secrets("some text AKIA" + "A" * 16, exact_only=True)
    assert called["n"] == 0
    OG.scrub_secrets("some text AKIA" + "A" * 16)
    assert called["n"] == 1


# ===========================================================================
# HELD interactions -- checked, and holding; the test is the proof
# ===========================================================================

def _gate(ws, cmd, mode="default", attended=False, **perm_kw):
    def _cb(tool, args, preview):
        return False
    perms = KitToolPermissions(
        mode=mode, workspace=str(ws),
        confirm_callback=_cb if attended else None, **perm_kw)
    return _doc_executor._run_permission_gate("bash", {"command": cmd}, perms)


def test_held_an_mcp_shell_is_not_admitted_by_the_caged_rule(home, ws,
                                                             monkeypatch):
    """Interaction 1/3: the caged-unasked rule is reachable only from the
    native shell dispatch (it passes caged_admission). An MCP shell routes
    through _gate_mcp_tool, which has no caged_admission list, so a command
    off the auto-allow list is asked about, never silently caged."""
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: True)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    perms = KitToolPermissions(mode="default", workspace=str(ws),
                               confirm_callback=None)
    # Native path admits it (caged).
    admitted: list = []
    assert _doc_executor._run_permission_gate(
        "bash", {"command": "./run.sh"}, perms,
        caged_admission=admitted) is None
    assert admitted == [True]
    # MCP shell path: no caged admission, so an off-list command with
    # nobody to ask is refused rather than run unasked.
    verdict = _doc_executor._gate_mcp_tool(
        "mcp__srv__bash", {"command": "./run.sh"}, perms)
    assert verdict is not None
    assert "auto-allow" in verdict or "not on" in verdict.lower()


def test_held_the_approval_answer_routes_are_refused_in_every_mode(home, ws):
    """Interaction 3: a shell command that answers an approval, or names an
    answer location, is refused in every mode -- the agent cannot approve
    its own request through the bash tool, caged or not."""
    (home / ".delfin" / "approval_codes").mkdir(parents=True, exist_ok=True)
    routes = [
        "delfin-agent approvals approve r1",
        "delfin-agent approvals deny r1",
        "cat ~/.delfin/approval_codes/r1.code",
    ]
    for mode in ("default", "acceptEdits", "bypassPermissions"):
        for cmd in routes:
            v = _gate(ws, cmd, mode=mode, attended=True,
                      session_read_dirs=(home / ".delfin",))
            assert v is not None and "blocked" in v.lower(), (mode, cmd, v)


def test_held_the_credential_store_and_keyring_reads_are_refused(home, ws):
    """Interaction 6: the credentials file is on the deny list and the
    commands that print a keyring entry are refused, in every mode."""
    (home / ".delfin" / "credentials.json").write_text("{}")
    cmds = [
        "cat ~/.delfin/credentials.json",
        "python3 -m keyring get svc KIT_TOOLBOX_API_KEY",
        "secret-tool lookup service delfin",
        "security find-generic-password -s delfin",
    ]
    for mode in ("default", "acceptEdits", "bypassPermissions"):
        for cmd in cmds:
            assert _gate(ws, cmd, mode=mode, attended=True) is not None, \
                (mode, cmd)


def test_held_a_persisted_bash_isolation_off_is_migrated_to_auto(home):
    """Interaction 2: a settings file that still carries
    agent.bash_isolation = "off" (set once, long ago) is migrated to "auto"
    so a locked or unattended session is confined again -- the one-time
    security update."""
    from delfin import user_settings as US
    out = US._apply_security_updates(
        {"agent": {"bash_isolation": "off"}})
    assert out["agent"]["bash_isolation"] == "auto"
    assert out.get("security_notices")
    # A later explicit "off" is respected (the update runs once).
    again = US._apply_security_updates(out | {"agent": {"bash_isolation": "off"}})
    assert again["agent"]["bash_isolation"] == "off"


def test_held_a_caged_unasked_command_never_gets_an_open_network(home, ws,
                                                                 monkeypatch):
    """Interaction 2/5: the caged-unasked admission is withheld when the
    network mode is "open" -- an unasked command must not reach an
    unrestricted network."""
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: ("bwrap", "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: True)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    monkeypatch.setattr("delfin.agent.egress_proxy.network_mode",
                        lambda: "open")
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    assert A._caged_command_may_run_unasked("./run.sh", perms, "default") is False


def test_held_network_open_is_capped_to_proxy_under_a_managed_ceiling(
        monkeypatch):
    """Interaction 5: a managed network.allowed_domains_max turns a user's
    sandbox_network.mode="open" into "proxy" so the ceiling can be
    enforced."""
    import delfin.agent.egress_proxy as EP
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "open"})
    monkeypatch.setattr("delfin.agent.managed_policy.network_ceiling_active",
                        lambda: True)
    assert EP.network_mode() == "proxy"
