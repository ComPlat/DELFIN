"""Run a single user query through the dashboard agent in plan-mode.

Safety guarantees (FOUR layers — the agent CANNOT damage the host):

1. ``--permission-mode plan`` — Bash, Edit, Write, NotebookEdit are
   hard-blocked by the CLI. The model emits tool-use blocks but the
   CLI denies their execution before they hit the OS.
2. Every mutating MCP tool listed by name in ``_hard_blocked_mcp`` is
   added to ``--disallowed-tools``, so even if a future tool ships
   with ``allow_mutate=True`` defaulted, plan-mode still blocks it.
3. ``--add-dir`` is restricted to the temporary sandbox only. The
   real DELFIN repo + the user's calc/ + archive/ are NOT in the
   writable surface; even read-only tools can't see them.
4. A pre-flight wall in ``run_agent_dryrun`` REFUSES to run if cwd
   sits inside the real DELFIN repo. A typo can't redirect writes at
   the user's actual source tree.

The agent still emits the tool_use blocks it WANTS to call — that's
exactly what we assert on.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class ToolCall:
    """One tool the model attempted (allowed or denied) in plan-mode."""
    name: str
    args: dict[str, Any]
    blocked: bool = False  # filled when a permission_denial event mentions it


@dataclass
class AgentDryRunResult:
    """Everything we captured from one plan-mode invocation."""
    user_query: str
    text_replies: list[str] = field(default_factory=list)
    tool_calls: list[ToolCall] = field(default_factory=list)
    permission_denials: list[str] = field(default_factory=list)
    cost_usd: float = 0.0
    stop_reason: str = ""
    timed_out: bool = False
    raw_events: list[dict] = field(default_factory=list)
    stderr_tail: str = ""

    @property
    def all_text(self) -> str:
        return "\n".join(self.text_replies)

    @property
    def tool_names(self) -> list[str]:
        return [c.name for c in self.tool_calls]


def _load_dashboard_prompt(sandbox_root: str | None = None) -> str:
    """Read the dashboard agent prompt verbatim with a synthetic live state.

    When ``sandbox_root`` is given the live-state block points at it so
    the agent's ACTION: lines reference paths it can actually see.
    """
    pack = (
        Path(__file__).resolve().parents[2]
        / "delfin" / "agent" / "pack" / "agents" / "dashboard_agent.md"
    )
    base = pack.read_text(encoding="utf-8")
    root = sandbox_root or "/tmp/dryrun_calc"
    live_state = (
        "\n--- Live state ---\n"
        "Branch: GUPPY (sandboxed dry-run)\n"
        "CONTROL.txt: functional=BP86, main_basisset=def2-SVP, PAL=8, maxcore=2000\n"
        "ORCA Builder: idle\n"
        f"Active calc folder: {root}/calc/test_a\n"
        f"Project root: {root}\n"
        "Recent jobs: (none)\n"
        "Permissions: plan-mode (every mutation will be blocked).\n"
        "Sandbox: this is an isolated tmp tree, not the real DELFIN repo. "
        "Don't write outside it.\n"
    )
    return base + live_state


def _ensure_ops_mcp_config() -> str:
    """Generate the delfin-ops MCP config file the CLI needs."""
    from delfin.ops_server.config import ensure_mcp_config
    return ensure_mcp_config()


def run_agent_dryrun(
    user_query: str,
    *,
    model: str = "haiku",
    timeout_s: int = 90,
    cwd: str | None = None,
    extra_disallowed_tools: list[str] | None = None,
    sandbox_root: str | None = None,
) -> AgentDryRunResult:
    """Drive the dashboard agent for one turn and return what it tried.

    Uses plan-mode so mutating tools are hard-blocked at the CLI layer.
    The model still emits the tool_use blocks it wanted to call — those
    are what we assert on.

    Parameters
    ----------
    user_query : str
        The user's message verbatim (German is fine).
    model : str
        ``haiku`` (cheap, fast — default), ``sonnet`` or ``opus``.
    timeout_s : int
        Hard wall-clock timeout. Plan-mode usually returns in 5-30 s.
    cwd : str, optional
        Working directory for the subprocess (mostly cosmetic).
    extra_disallowed_tools : list of str, optional
        Extra tool names to deny (on top of plan-mode's defaults).
    """
    claude_path = shutil.which("claude")
    if not claude_path:
        raise FileNotFoundError(
            "claude CLI not found on $PATH — cannot run dry-run tests"
        )

    # SAFETY WALL — never let the agent run with the real DELFIN repo as
    # cwd. Every dry-run test must be inside a tmp sandbox. Without this
    # wall a typo could send writes/reads at the user's actual code.
    effective_cwd = cwd or sandbox_root
    if effective_cwd:
        cwd_real = Path(effective_cwd).resolve()
        repo_real = Path(__file__).resolve().parents[2]  # delfin repo root
        try:
            cwd_real.relative_to(repo_real)
            raise RuntimeError(
                f"REFUSING dry-run: cwd '{cwd_real}' is inside the real "
                f"DELFIN repo '{repo_real}'. Use a tmp sandbox via "
                f"tests.agent_dryrun.sandbox.sandboxed()."
            )
        except ValueError:
            pass  # cwd is outside the repo — exactly what we want

    effective_root = sandbox_root or cwd
    system_prompt = _load_dashboard_prompt(sandbox_root=effective_root)
    mcp_config = _ensure_ops_mcp_config()

    cmd = [
        claude_path,
        "--print",
        "--permission-mode", "plan",
        "--output-format", "stream-json",
        "--verbose",
        "--include-hook-events",
        "--model", model,
        "--append-system-prompt", system_prompt,
        "--mcp-config", mcp_config,
        "--no-session-persistence",
    ]
    # When a sandbox is given, restrict --add-dir to it so the agent's
    # writable surface is exactly that one tmp tree (plan-mode already
    # hard-blocks mutations, but extra layered safety doesn't hurt).
    if sandbox_root:
        cmd.extend(["--add-dir", sandbox_root])

    # Defense-in-depth: even though plan-mode hard-blocks Bash/Edit/Write,
    # an MCP tool with allow_mutate=True default would slip past plan-mode
    # because the CLI doesn't see "mutating" as part of its permission
    # model. Pre-deny EVERY mutating MCP tool by name so the test machine
    # can't be touched even if a future tool ships with allow_mutate=True
    # by accident. Read-only and dry-run-default tools stay enabled.
    _hard_blocked_mcp = [
        "mcp__delfin-ops__cleanup",
        "mcp__delfin-ops__stop",
        "mcp__delfin-ops__pipeline_run",
        "mcp__delfin-ops__pipeline_prepare",
        "mcp__delfin-ops__run_orca_input",
        "mcp__delfin-ops__co2",
        "mcp__delfin-ops__tadf_xtb",
        "mcp__delfin-ops__hyperpol",
        "mcp__delfin-ops__submit_calculation",
        "mcp__delfin-ops__cancel_calculation",
        "mcp__delfin-ops__rename_calc_folder",
        "mcp__delfin-ops__create_calc_folder",
        "mcp__delfin-ops__move_calc_folder",
        "mcp__delfin-ops__move_to_archive",
        "mcp__delfin-ops__delete_calc_folder",
        "mcp__delfin-ops__kill_all_user_jobs",
        "mcp__delfin-ops__prepare_recalc",
        "mcp__delfin-ops__run_calc_option",
        "mcp__delfin-ops__index_new_pdf",
    ]
    blocked = list(_hard_blocked_mcp)
    if extra_disallowed_tools:
        blocked.extend(extra_disallowed_tools)
    cmd.extend(["--disallowed-tools", ",".join(blocked)])

    env = dict(os.environ)
    # Remove any auto-loaded API key so we go through the CLI's OAuth path
    env.pop("ANTHROPIC_API_KEY", None)

    proc = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=cwd,
        env=env,
    )
    try:
        stdout, stderr = proc.communicate(input=user_query, timeout=timeout_s)
    except subprocess.TimeoutExpired:
        proc.kill()
        stdout, stderr = proc.communicate()
        result = AgentDryRunResult(user_query=user_query, timed_out=True)
        result.stderr_tail = (stderr or "")[-500:]
        return result

    result = AgentDryRunResult(user_query=user_query)
    result.stderr_tail = (stderr or "")[-500:]

    for line in stdout.splitlines():
        if not line.strip():
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        result.raw_events.append(event)
        etype = event.get("type")

        if etype == "assistant":
            msg = event.get("message", {})
            for block in msg.get("content", []) or []:
                btype = block.get("type")
                if btype == "text":
                    text = (block.get("text") or "").strip()
                    if text:
                        result.text_replies.append(text)
                elif btype == "tool_use":
                    result.tool_calls.append(ToolCall(
                        name=block.get("name", ""),
                        args=block.get("input", {}) or {},
                    ))

        elif etype == "result":
            result.cost_usd = float(event.get("total_cost_usd", 0.0) or 0.0)
            result.stop_reason = event.get("stop_reason", "") or ""
            for denial in event.get("permission_denials", []) or []:
                if isinstance(denial, dict):
                    result.permission_denials.append(
                        denial.get("tool_name", "") or str(denial)
                    )
                else:
                    result.permission_denials.append(str(denial))

    # Mark calls that appear in permission_denials as blocked
    for call in result.tool_calls:
        if call.name in result.permission_denials:
            call.blocked = True

    return result


def has_tool_call(result: AgentDryRunResult, *substrings: str) -> bool:
    """Return True iff the agent attempted a tool whose name contains any
    of ``substrings`` (case-insensitive).

    "Attempted" includes calls that were later blocked by plan-mode —
    the test's purpose is to measure the agent's intent, not whether
    the CLI's permission layer let the call through. The fallback
    (Glob / Grep / etc.) the agent picks AFTER a block doesn't count
    against this check.
    """
    needle = [s.lower() for s in substrings]
    for call in result.tool_calls:
        nm = call.name.lower()
        if any(n in nm for n in needle):
            return True
    return False


def tool_call_args(
    result: AgentDryRunResult,
    name_substring: str,
) -> dict[str, Any] | None:
    """Return the args of the first tool call whose name contains the substring."""
    needle = name_substring.lower()
    for call in result.tool_calls:
        if needle in call.name.lower():
            return call.args
    return None


def _smoke():
    """Quick CLI smoke; running this module directly fires one trivial query."""
    t0 = time.time()
    res = run_agent_dryrun(
        "Was ist 1+1? Antworte in einer Zeile.",
        model="haiku",
        timeout_s=45,
    )
    dt = time.time() - t0
    print(f"--- {dt:.1f}s, ${res.cost_usd:.4f}, "
          f"{len(res.tool_calls)} tools, {len(res.text_replies)} replies")
    print(f"reply: {res.all_text[:200]}")
    print(f"tools: {res.tool_names}")
    if res.timed_out:
        print(f"TIMED OUT, stderr: {res.stderr_tail}")


if __name__ == "__main__":
    _smoke()
