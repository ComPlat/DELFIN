"""Safety invariants for the solo agent — the 'safe to use' contract.

These lock in the guarantees that make the solo coding agent safe to run
in production: a bash deny-list that blocks destructive commands, a
self-modification guard over the core files, auto-verification after code
edits, a bounded tool whitelist, capped tool rounds, and a cost soft
limit.  If any of these is accidentally removed or weakened, a test fails
— so the safe stand can't silently regress.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import KitToolPermissions


def _bwrap_functional() -> bool:
    """True only if bwrap is installed AND can actually sandbox. On unprivileged
    CI containers bwrap is often present but cannot create user namespaces, so
    it exits non-zero — `which("bwrap")` alone would wrongly let the test run."""
    import shutil
    import subprocess
    if shutil.which("bwrap") is None:
        return False
    try:
        r = subprocess.run(
            ["bwrap", "--ro-bind", "/", "/", "true"],
            capture_output=True, timeout=10,
        )
        return r.returncode == 0
    except Exception:
        return False


@pytest.fixture
def perms():
    return KitToolPermissions(workspace="/tmp")


# ---------------------------------------------------------------------------
# Bash deny-list blocks destructive commands, allows safe ones
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "rm -rf /",
    "rm -rf /home/user/data",
    "sudo apt install x",
    "mkfs.ext4 /dev/sda",
    "shutdown now",
    "dd if=/dev/zero of=/dev/sda",
    "git push --force origin main",
    "git reset --hard HEAD~5",
])
def test_dangerous_bash_is_denied(perms, cmd):
    assert perms.matches_bash_deny(cmd) is not None, f"NOT denied: {cmd}"


@pytest.mark.parametrize("cmd", [
    "ls -la",
    "python -m pytest tests/ -q",
    "git status",
    "grep -rn foo delfin/",
    "git push --force-with-lease origin feature",   # safe force variant
])
def test_safe_bash_is_allowed(perms, cmd):
    assert perms.matches_bash_deny(cmd) is None, f"wrongly denied: {cmd}"


# ---------------------------------------------------------------------------
# Self-modification guard protects the core agent files
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", [
    "delfin/agent/api_client.py",
    "delfin/agent/engine.py",
    "delfin/dashboard/tab_agent.py",
    "delfin/agent/kit_confirm.py",
    # The rest of the permission/sandbox layer + persisted config that
    # decides what runs unattended must be guarded too.
    "delfin/agent/sandbox.py",
    "delfin/agent/kit_settings.py",
    "delfin/user_settings.py",
    ".delfin/settings.json",
    "myproj/.delfin/settings.json",
])
def test_core_files_are_protected(perms, path):
    assert perms.matches_path_protected(path) is True, f"NOT protected: {path}"


def test_ordinary_file_is_not_protected(perms):
    assert perms.matches_path_protected("agent_workspace/task/notes.md") is False
    # A normal source file the user works on stays freely editable.
    assert perms.matches_path_protected("delfin/agent/model_routing.py") is False


# ---------------------------------------------------------------------------
# Auto-verify after code edits (solo + builder roles)
# ---------------------------------------------------------------------------

def _solo_engine(tmp_path):
    from unittest.mock import MagicMock, patch
    from delfin.agent.engine import AgentEngine
    # Stub the client (CI has no `claude` CLI); check_auto_verify inspects the
    # tool name/output, never the client, so a MagicMock keeps this in CI.
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        return AgentEngine(repo_dir=str(tmp_path), mode="solo")


def test_auto_verify_fires_after_python_edit(tmp_path):
    eng = _solo_engine(tmp_path)
    msg = eng.check_auto_verify("Edit", "edited /repo/foo.py")
    assert msg and "pytest" in msg


def test_auto_verify_silent_for_read(tmp_path):
    eng = _solo_engine(tmp_path)
    assert eng.check_auto_verify("Read", "read /repo/foo.py") is None


# ---------------------------------------------------------------------------
# Bounded tool whitelist
# ---------------------------------------------------------------------------

def test_solo_whitelist_is_bounded():
    from delfin.agent.engine import _SOLO_TOOLS
    assert {"Read", "Edit", "Write", "Bash", "Grep", "Glob"} <= _SOLO_TOOLS
    # No surprise tools beyond the documented code + research surface.
    assert _SOLO_TOOLS <= {
        "Read", "Grep", "Glob", "Bash", "Edit", "Write",
        "WebSearch", "WebFetch",
    }


# ---------------------------------------------------------------------------
# Bounded tool rounds + cost soft limit exist
# ---------------------------------------------------------------------------

def test_tool_rounds_are_bounded():
    from delfin.agent import model_profiles as mp
    # Every profile must cap tool rounds (no unbounded loops).
    for name in dir(mp):
        prof = getattr(mp, name)
        rounds = getattr(prof, "max_tool_rounds", None)
        if isinstance(rounds, int):
            assert 0 < rounds <= 100, f"{name}: implausible max_tool_rounds={rounds}"


def test_cost_soft_limit_default_present():
    from delfin.user_settings import DEFAULT_SETTINGS
    assert "cost_soft_limit_usd" in DEFAULT_SETTINGS["agent"]
    assert DEFAULT_SETTINGS["agent"]["cost_soft_limit_usd"] > 0


# ---------------------------------------------------------------------------
# Out-of-workspace access is refused with an actionable message
# (the Jerome case: project in another user's home, outside the sandbox)
# ---------------------------------------------------------------------------

def test_outside_workspace_path_rejected_with_grant_hint(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions
    ex = _DocToolExecutor()
    perms = KitToolPermissions(workspace=str(tmp_path))
    resolved, err = ex._resolve_in_workspace(
        "/home/other_user/project/experiments_all.xlsx", perms
    )
    assert resolved is None
    # Refused...
    assert "outside the allowed workspace roots" in err
    # ...AND tells the user how to fix it (grant the path), not a dead end.
    assert "--add-dir" in err or "extra_workspace_dirs" in err
    assert "GRANT" in err.upper()


def test_granting_a_dir_makes_it_accessible(tmp_path):
    from delfin.agent.api_client import KitToolPermissions
    ws = tmp_path / "ws"
    extra = tmp_path / "granted"
    ws.mkdir(); extra.mkdir()
    perms = KitToolPermissions(workspace=str(ws),
                               extra_workspace_dirs=(extra,))
    # A path under the granted dir resolves to a real root (not rejected).
    assert perms.find_root_for((extra / "file.xlsx").resolve()) is not None


# ---------------------------------------------------------------------------
# Reasoning models (gpt-5.x / o3 / o4) MUST receive tools — withholding them
# was why the agent "did nothing" / claimed no filesystem tool on Azure GPT-5.x
# ---------------------------------------------------------------------------

def test_tools_not_gated_off_for_reasoning_models():
    import re
    from pathlib import Path
    src = (Path(__file__).resolve().parent.parent
           / "delfin" / "agent" / "api_client.py").read_text(encoding="utf-8")
    # Find the line that advertises tools to the OpenAI request.
    idx = src.find('kwargs["tools"] = advertised_tools')
    assert idx > 0, "tool-advertisement line not found"
    # The guarding condition (the preceding `if ...:`) must NOT withhold
    # tools from reasoning models.
    guard = src[max(0, idx - 200):idx]
    assert "not is_reasoning" not in guard, (
        "tools are gated off for reasoning models — gpt-5.x would get no tools"
    )


def test_protected_path_detection_for_confirm_dialog():
    from delfin.agent.kit_confirm import _is_protected_path
    # Core self-mod-guarded files → protected (persist intentionally off).
    assert _is_protected_path("delfin/agent/engine.py") is True
    assert _is_protected_path("/x/software/delfin/delfin/dashboard/tab_agent.py") is True
    # Ordinary user files → NOT protected (so the dialog points to
    # acceptEdits instead of a misleading 'protected file' message).
    assert _is_protected_path("/home/u/calc/testrechnung/run.inp") is False
    assert _is_protected_path("experiments_all.xlsx") is False
    assert _is_protected_path("") is False


def test_gpt5_family_is_detected_as_reasoning():
    # Documents the classification the fix depends on: gpt-5.x needs
    # reasoning_effort AND tools (the bug was withholding the tools).
    import re
    for base in ("gpt-5", "gpt-5.1", "gpt-5.4", "gpt-5.5", "o3", "o4-mini"):
        is_reasoning = bool(re.match(r"^o\d", base) or base.startswith("gpt-5"))
        assert is_reasoning, f"{base} should be classified reasoning"


# ---------------------------------------------------------------------------
# Streaming fallback: proxy "can't stream" errors are detected (and only those)
# ---------------------------------------------------------------------------

def test_stream_unsupported_error_detection():
    from delfin.agent.api_client import _is_stream_unsupported_error
    prod = Exception(
        "Error code: 400 - {'detail': \"'async for' requires an object "
        "with __aiter__ method, got ModelResponse\"}"
    )
    assert _is_stream_unsupported_error(prod) is True
    assert _is_stream_unsupported_error(Exception("401 AuthenticationError")) is False
    assert _is_stream_unsupported_error(Exception("rate limit")) is False


# ---------------------------------------------------------------------------
# Greetings classify as simple → minimal reasoning budget (fast "Hallo")
# ---------------------------------------------------------------------------

def test_greetings_classify_simple_for_fast_turns():
    from delfin.agent.engine import AgentEngine, _COMPLEXITY_THINKING_MULT
    for t in ("Hallo", "hi", "Hallo, was kannst du?", "danke!", "ok"):
        assert AgentEngine.classify_task_complexity(t) == "simple", t
    # A greeting prefix on a REAL task must not downgrade it.
    assert AgentEngine.classify_task_complexity(
        "Hallo, bitte refactor delfin/agent/engine.py") == "complex"
    # simple budget maps to the OpenAI reasoning_effort "low" bucket (<16k).
    assert int(64000 * _COMPLEXITY_THINKING_MULT["simple"]) < 16000


# ---------------------------------------------------------------------------
# Bash filesystem isolation (Stufe 5): only granted roots are writable
# ---------------------------------------------------------------------------

def test_bash_isolation_off_is_plain_bash(tmp_path):
    from delfin.agent.api_client import _bash_isolation_argv, KitToolPermissions
    perms = KitToolPermissions(workspace=str(tmp_path))
    assert _bash_isolation_argv("echo hi", tmp_path, perms, mode="off") == [
        "/bin/bash", "-c", "echo hi"]


def test_bash_isolation_default_setting_is_off():
    from delfin.user_settings import DEFAULT_SETTINGS
    assert DEFAULT_SETTINGS["agent"]["bash_isolation"] == "off"


@pytest.mark.skipif(not _bwrap_functional(),
                    reason="bwrap not installed or not functional (e.g. unprivileged CI container)")
def test_bwrap_isolation_blocks_writes_outside_granted_roots(tmp_path):
    import subprocess
    from pathlib import Path
    from delfin.agent.api_client import _bash_isolation_argv, KitToolPermissions
    perms = KitToolPermissions(workspace=str(tmp_path))
    # Inside the granted workspace: writable.
    r1 = subprocess.run(
        _bash_isolation_argv(f"touch {tmp_path}/inside.txt", tmp_path, perms,
                             mode="bwrap"),
        capture_output=True, text=True, timeout=30)
    assert r1.returncode == 0 and (tmp_path / "inside.txt").exists()
    # Outside (the user's real home): read-only → blocked.
    outside = Path.home().resolve() / "delfin_isolation_probe.txt"
    r2 = subprocess.run(
        _bash_isolation_argv(f"touch {outside}", tmp_path, perms,
                             mode="bwrap"),
        capture_output=True, text=True, timeout=30)
    assert r2.returncode != 0 and not outside.exists()


@pytest.mark.skipif(not _bwrap_functional(),
                    reason="bwrap not installed or not functional (e.g. unprivileged CI container)")
def test_bwrap_isolation_grants_extra_dirs(tmp_path):
    import subprocess
    from delfin.agent.api_client import _bash_isolation_argv, KitToolPermissions
    ws = tmp_path / "ws"; ws.mkdir()
    extra = tmp_path / "granted"; extra.mkdir()
    perms = KitToolPermissions(workspace=str(ws),
                               extra_workspace_dirs=(extra,))
    r = subprocess.run(
        _bash_isolation_argv(f"touch {extra}/ok.txt", ws, perms,
                             mode="bwrap"),
        capture_output=True, text=True, timeout=30)
    assert r.returncode == 0 and (extra / "ok.txt").exists()
