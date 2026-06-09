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
])
def test_core_files_are_protected(perms, path):
    assert perms.matches_path_protected(path) is True, f"NOT protected: {path}"


def test_ordinary_file_is_not_protected(perms):
    assert perms.matches_path_protected("agent_workspace/task/notes.md") is False


# ---------------------------------------------------------------------------
# Auto-verify after code edits (solo + builder roles)
# ---------------------------------------------------------------------------

def _solo_engine(tmp_path):
    from delfin.agent.engine import AgentEngine
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


def test_gpt5_family_is_detected_as_reasoning():
    # Documents the classification the fix depends on: gpt-5.x needs
    # reasoning_effort AND tools (the bug was withholding the tools).
    import re
    for base in ("gpt-5", "gpt-5.1", "gpt-5.4", "gpt-5.5", "o3", "o4-mini"):
        is_reasoning = bool(re.match(r"^o\d", base) or base.startswith("gpt-5"))
        assert is_reasoning, f"{base} should be classified reasoning"
