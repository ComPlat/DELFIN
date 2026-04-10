"""DELFIN Agent Engine: orchestrates Claude conversations with role-based prompts."""

from __future__ import annotations

import json
import re
import subprocess
import threading
from pathlib import Path
from typing import Any, Callable

from .api_client import StreamEvent, create_client
from .prompt_loader import PromptLoader


# -- Legacy mode name migration ------------------------------------------------
_LEGACY_MODE_MAP = {
    "default": "quick",
    "high_risk": "reviewed",
    "runtime_cluster": "cluster",
    "release": "full",
}


def _migrate_mode(mode: str) -> str:
    """Map old mode names to new ones."""
    return _LEGACY_MODE_MAP.get(mode, mode)


# -- Mode auto-detection rules ------------------------------------------------
# Maps file patterns to the minimum mode that should be used.
# Order matters: first match wins (most specific first).
_ESCALATION_PATTERNS: list[tuple[str, str]] = [
    # cluster mode triggers
    (r"runtime_setup\.py", "cluster"),
    (r"qm_runtime\.py", "cluster"),
    (r"backend_slurm\.py", "cluster"),
    (r"backend_local\.py", "cluster"),
    (r"submit_templates/", "cluster"),
    (r"orca_recovery\.py", "cluster"),
    # reviewed mode triggers — infrastructure
    (r"delfin/cli\.py", "reviewed"),
    (r"pipeline\.py", "reviewed"),
    (r"parallel_classic_manually\.py", "reviewed"),
    (r"delfin/api\.py", "reviewed"),
    (r"config\.py", "reviewed"),
    # reviewed mode triggers — chemistry workflows
    (r"occupier\.py", "reviewed"),
    (r"occupier_auto\.py", "reviewed"),
    (r"esd_module\.py", "reviewed"),
    (r"esd_input_generator\.py", "reviewed"),
    (r"orca\.py", "reviewed"),
    (r"xtb_crest\.py", "reviewed"),
    (r"build_up_complex", "reviewed"),
    (r"stability_constant\.py", "reviewed"),
    (r"calculators\.py", "reviewed"),
    (r"ensemble_nmr\.py", "reviewed"),
    (r"smiles_converter\.py", "reviewed"),
    (r"tadf_xtb\.py", "reviewed"),
    (r"hyperpol\.py", "reviewed"),
]

# Mode escalation order (higher index = heavier mode)
_MODE_RANK = {
    "dashboard": 0,
    "research": 0,
    "solo": 0,
    "quick": 1,
    "reviewed": 2,
    "cluster": 3,
    "full": 4,
}

_DASHBOARD_KEYWORDS = (
    "/control", "/orca", "/calc", "/remote", "/job", "/submit", "/analyze",
    "dashboard", "control key", "orca builder", "calc browser", "remote archive",
    "job status", "set control", "submit job", "show jobs", "browse calculations",
    "open calc", "open calculation", "switch tab", "agent tab", "recalc",
    # UI manipulation (natural language → dashboard agent)
    "button", "btn", "widget", "style", "farbe", "color", "colour",
    "rot ", "grün", "blau", "gelb", "schwarz", "weiß",
    "red ", "green ", "blue ", "yellow ",
    "disable", "enable", "versteck", "zeig ",
    "send-btn", "send button", "stop button", "submit button",
    "input feld", "dropdown", "sichtbar", "unsichtbar",
    # CONTROL / ORCA setting commands (natural language)
    "setze ", "functional", "basis ", "basisset", "charge ", "mult ",
    "pal ", "maxcore", "dispersion", "solvent",
    "job name", "job-name", "ordner ",
    # Calc/analysis operations
    "tabelle", "table ", "report", "energien", "energie ",
    "konvergenz", "convergence", "archiv", "visuali",
)
_CHEMISTRY_KEYWORDS = (
    # Core QC concepts
    "dft", "functional", "basis set", "dispersion", "solvation", "solvent model",
    "redox", "nmr", "excited state", "uv-vis", "spin state", "thermochemistry",
    "electrochem", "orca method", "smiles", "metal complex", "ligand", "conformer",
    "geometry optimization", "oxidation state", "coordination", "reaction mechanism",
    # QC methods and programs
    "crest", "censo", "anmr", "gfn", "gfn2", "semiempirical",
    "cam-b3lyp", "pbe0", "tpss", "b97", "range-separated", "hybrid functional",
    "multireference", "ccsd", "coupled cluster", "casscf", "broken symmetry",
    # Properties and analysis
    "frequency", "vibration", "ir spectrum", "transition state", "irc",
    "homo", "lumo", "nbo", "population analysis", "spin-orbit",
    # DELFIN-specific workflows
    "occupier", "preoptimization", "preopt", "qm/mm", "molecular dynamics",
    "stability constant", "log k", "hyperpolarizability", "tadf",
)
_CODE_CHANGE_KEYWORDS = (
    "fix", "implement", "change", "patch", "refactor", "add", "update", "edit",
    "rewrite", "improve", "extend", "optimize", "debug", "repair", "modify",
)
_CODE_QUESTION_KEYWORDS = (
    "how does", "how do", "why does", "what does", "explain", "where is",
    "show me", "walk through", "understand", "question", "why is",
)
_FILE_OR_CODE_HINTS = (
    ".py", ".md", ".yaml", ".yml", ".json", "pytest", "test_", "stack trace",
    "traceback", "function", "class ", "module", "repo", "codebase", "diff",
)
_FULL_RISK_KEYWORDS = (
    "release", "milestone", "broad architecture", "major redesign",
    "multi entry point", "cross-cutting", "large refactor", "final confidence",
)
_REVIEWED_RISK_KEYWORDS = (
    "api semantics", "public behavior", "config semantics", "parser", "validation",
    "architecture", "cross-module", "user-facing behavior", "result semantics",
)
_CLUSTER_RISK_KEYWORDS = (
    "cluster", "slurm", "sbatch", "squeue", "scancel", "runtime", "scheduler",
    "scratch", "restart", "recovery", "backend_local", "backend_slurm", "submit template",
)

# -- Role-specific thinking budgets -------------------------------------------
# Builder needs the most thinking (implementation), review roles need moderate,
# planning roles need less.
# Per-role model routing: use the right model for each job.
# "auto" means use the user-selected model (no override).
# Role → model mapping.  "auto" = user's choice from the dropdown.
# Planning and review roles need quality models — haiku can't plan well.
_ROLE_MODEL_MAP: dict[str, str] = {
    "chief_agent": "sonnet",      # strategic decisions need quality
    "session_manager": "sonnet",  # planning is the most critical step
    "critic_agent": "sonnet",     # deep code review needs quality
    "runtime_agent": "haiku",     # runtime categorization
    "reviewer_agent": "sonnet",   # code review needs quality
    "builder_agent": "auto",      # user's choice (often sonnet/opus)
    "test_agent": "haiku",        # run tests, assert results
    "solo_agent": "auto",         # user's choice
    "dashboard_agent": "sonnet",  # data analysis + research needs quality
    "research_agent": "sonnet",   # chemistry method synthesis needs depth
}

_ROLE_THINKING_BUDGETS: dict[str, int] = {
    "chief_agent": 16000,         # strategic decisions need depth
    "session_manager": 32000,     # planning is critical — needs deep analysis
    "critic_agent": 16000,        # deep analysis
    "runtime_agent": 8000,        # categorization
    "builder_agent": 50000,       # complex implementation
    "reviewer_agent": 16000,      # code review needs depth
    "test_agent": 8000,           # test execution
    "research_agent": 16000,      # deep chemistry method analysis
    "solo_agent": 64000,          # scales: low=32k, medium=64k, high=128k
    "dashboard_agent": 64000,     # scales: low=32k, medium=64k, high=128k
}
_DEFAULT_THINKING_BUDGET = 10000

# Max output tokens per role.  Builder/Solo need long responses for
# complex implementations.  Review roles need less.
_ROLE_MAX_TOKENS: dict[str, int] = {
    "builder_agent": 32768,
    "solo_agent": 32768,
    "session_manager": 16384,
    "critic_agent": 16384,
    "reviewer_agent": 16384,
    "chief_agent": 8192,
    "research_agent": 8192,
    "test_agent": 8192,
    "runtime_agent": 8192,
    "dashboard_agent": 16384,    # analysis scripts + research reports need space
}
_DEFAULT_MAX_TOKENS = 8192

# Code-level tool whitelist per role.
# If a role emits a tool_use event for a tool NOT in its whitelist,
# the engine silently blocks it.  This prevents prompt-injection or
# model hallucination from bypassing role restrictions.
_READ_TOOLS = frozenset({"Read", "Grep", "Glob"})
_GIT_BASH = frozenset({"Read", "Grep", "Glob", "Bash"})  # Bash limited by prompt to git
_RESEARCH_TOOLS = frozenset({"Read", "Grep", "Glob", "Bash", "WebSearch", "WebFetch"})
_FULL_TOOLS = frozenset({"Read", "Grep", "Glob", "Bash", "Edit", "Write"})
_SOLO_TOOLS = frozenset({
    "Read", "Grep", "Glob", "Bash", "Edit", "Write",  # full code access
    "WebSearch", "WebFetch",                            # research capabilities
})
_DASHBOARD_TOOLS = frozenset({
    "Read", "Grep", "Glob",        # read code + data
    "Write", "Bash",               # agent_workspace only (CLI enforced via --add-dir)
    "WebSearch", "WebFetch",       # literature research
})

# MCP documentation server tools — allowed for ALL roles so every agent
# can look up ORCA manual sections, xTB docs, methodology, etc.
_DOC_TOOL_PREFIX = "mcp__delfin-docs__"
_ROLE_TOOL_WHITELIST: dict[str, frozenset[str]] = {
    "dashboard_agent": _DASHBOARD_TOOLS,        # full analysis + research, writes restricted to workspace
    "research_agent":  _RESEARCH_TOOLS,         # web search + code reading
    "chief_agent":     _GIT_BASH,
    "session_manager": _GIT_BASH,
    "critic_agent":    _GIT_BASH,
    "reviewer_agent":  _GIT_BASH,
    "runtime_agent":   _GIT_BASH,
    "test_agent":      _FULL_TOOLS,
    "builder_agent":   _FULL_TOOLS,
    "solo_agent":      _SOLO_TOOLS,
}


class AgentEngine:
    """Core orchestration engine for the DELFIN agent.

    Manages conversation state, role transitions, session persistence,
    and Claude communication via CLI or API backend.

    Parameters
    ----------
    repo_dir : Path
        Root of the DELFIN repository.
    backend : str
        ``"cli"`` (default, uses OAuth) or ``"api"`` (needs API key).
    api_key : str
        Only needed for ``"api"`` backend.
    model : str
        Model name or alias.
    mode : str
        DELFIN_AGENT_LITE mode (default: ``"quick"``).
    """

    def __init__(
        self,
        repo_dir: Path,
        backend: str = "cli",
        provider: str = "claude",
        api_key: str = "",
        model: str = "",
        mode: str = "quick",
        permission_mode: str = "",
        pack_dir: Path | None = None,
        mcp_config: str = "",
        allowed_tools: list[str] | None = None,
        extra_dirs: list[str] | None = None,
        agent_workspace_dir: str = "",
    ):
        self.repo_dir = Path(repo_dir)
        self.loader = PromptLoader(repo_dir=pack_dir)
        self.client = create_client(
            backend=backend, provider=provider, api_key=api_key,
            model=model, permission_mode=permission_mode,
            cwd=str(self.repo_dir), mcp_config=mcp_config,
            allowed_tools=allowed_tools, extra_dirs=extra_dirs,
        )
        self.backend = backend
        self.provider = provider
        self.mode = mode
        self._agent_workspace_dir = agent_workspace_dir
        self.route: list[str] = []
        self.mode_description: str = ""
        self.current_role_index: int = 0
        self.messages: list[dict[str, Any]] = []
        self.role_outputs: dict[str, str] = {}
        self.token_usage = {"input": 0, "output": 0}
        self.cost_usd: float = 0.0
        self.session_id: str = ""  # CLI session ID for conversation persistence
        self._stop_requested = False
        self._lock = threading.Lock()

        self._load_mode(mode)

    def _load_mode(self, mode: str) -> None:
        """Load a mode definition from the LITE manifest."""
        mode = _migrate_mode(mode)
        mode_data = self.loader.load_mode(mode)
        self.route = mode_data["route"]
        self.mode_description = mode_data["description"]
        self.mode = mode

    @property
    def current_role(self) -> str:
        """Return the current role ID."""
        if not self.route or self.current_role_index >= len(self.route):
            return ""
        return self.route[self.current_role_index]

    @property
    def is_cycle_complete(self) -> bool:
        """True if all roles in the route have been visited."""
        return self.current_role_index >= len(self.route)

    def get_status(self) -> dict[str, Any]:
        """Return current engine status for the UI."""
        return {
            "mode": self.mode,
            "backend": self.backend,
            "provider": self.provider,
            "role": self.current_role,
            "role_index": self.current_role_index,
            "role_total": len(self.route),
            "input_tokens": self.token_usage["input"],
            "output_tokens": self.token_usage["output"],
            "cost_usd": self.cost_usd,
            "cycle_complete": self.is_cycle_complete,
            "session_id": self.session_id,
        }

    def _build_current_system_prompt(self, memory_context: str = "") -> str:
        """Build the system prompt for the current role."""
        role = self.current_role
        if not role:
            role = self.route[0] if self.route else "builder_agent"

        return self.loader.build_system_prompt(
            role_id=role,
            mode_id=self.mode,
            mode_description=self.mode_description,
            route=self.route,
            role_index=self.current_role_index,
            prior_outputs=self.role_outputs if self.role_outputs else None,
            memory_context=memory_context,
        )

    def stream_response(
        self,
        user_message: str,
        memory_context: str = "",
        on_token: Callable[[str], None] | None = None,
        on_tool_use: Callable[[str, str], None] | None = None,
        on_tool_result: Callable[[str, str], None] | None = None,
        on_permission_denied: Callable[[str], None] | None = None,
        on_thinking: Callable[[str], None] | None = None,
        thinking_budget: int = 0,
        max_tokens: int = 0,
    ) -> str:
        """Send a user message and stream the response.

        Parameters
        ----------
        user_message : str
            The user's message text.
        memory_context : str
            Optional persistent memory to inject.
        on_token : callable, optional
            Called with each text chunk as it arrives.
        on_tool_use : callable, optional
            Called with (tool_name, tool_input_json) when agent uses a tool.
        on_tool_result : callable, optional
            Called with (tool_name, tool_output) when a tool returns its result.
        on_permission_denied : callable, optional
            Called with (description) when a tool call was blocked.
        on_thinking : callable, optional
            Called with thinking text chunks as the model reasons.
        thinking_budget : int
            Extended thinking budget in tokens (0 = default/auto).
        max_tokens : int
            Max output tokens (0 = use role default).

        Returns
        -------
        str
            The complete assistant response text.
        """
        self._stop_requested = False

        self.messages.append({"role": "user", "content": user_message})
        system_prompt = self._build_current_system_prompt(memory_context)

        # Resolve max_tokens: caller override > role default > global default
        effective_max = max_tokens or self.max_tokens_for_role(self.current_role)

        chunks: list[str] = []
        try:
            for event in self.client.stream_message(
                system=system_prompt,
                messages=self.messages,
                max_tokens=effective_max,
                session_id=self.session_id,
                thinking_budget=thinking_budget,
            ):
                if self._stop_requested:
                    break

                if event.type == "text_delta" and event.text:
                    chunks.append(event.text)
                    if on_token:
                        on_token(event.text)

                elif event.type == "thinking_delta" and event.text:
                    if on_thinking:
                        on_thinking(event.text)

                elif event.type == "tool_use":
                    # Code-level tool whitelist enforcement
                    role_id = self.route[self.current_role_index] if self.route else ""
                    allowed = _ROLE_TOOL_WHITELIST.get(role_id)
                    if allowed is not None and event.tool_name not in allowed:
                        # Allow MCP doc server tools for all roles
                        if not event.tool_name.startswith(_DOC_TOOL_PREFIX):
                            # Block unauthorized tool — don't call on_tool_use
                            continue
                    # Defense-in-depth: dashboard_agent Write restricted to workspace
                    if (role_id == "dashboard_agent"
                            and event.tool_name == "Write"
                            and self._agent_workspace_dir):
                        try:
                            _parsed = json.loads(event.tool_input)
                            _fp = _parsed.get("file_path", "")
                            if _fp:
                                _resolved = str(Path(_fp).resolve())
                                _ws = str(Path(self._agent_workspace_dir).resolve())
                                if not _resolved.startswith(_ws + "/") and _resolved != _ws:
                                    continue  # block write outside workspace
                        except Exception:
                            pass
                    if on_tool_use:
                        on_tool_use(event.tool_name, event.tool_input)

                elif event.type == "tool_result":
                    if on_tool_result and event.tool_output:
                        on_tool_result(event.tool_name, event.tool_output)

                elif event.type == "permission_denied":
                    if on_permission_denied:
                        on_permission_denied(event.tool_name)

                elif event.type == "session_init":
                    # Capture session ID from CLI for persistence
                    if event.text:
                        self.session_id = event.text

                elif event.type == "message_start":
                    # Input tokens tracked here (authoritative count
                    # including cache).  Do NOT also add in message_delta.
                    with self._lock:
                        self.token_usage["input"] += event.input_tokens

                elif event.type == "message_delta":
                    with self._lock:
                        # Only output tokens and cost from the final event.
                        # Input tokens already counted in message_start.
                        self.token_usage["output"] += event.output_tokens
                        self.cost_usd += event.cost_usd
                    # Capture session ID from result event
                    if event.text and not self.session_id:
                        self.session_id = event.text

        except Exception:
            if not chunks:
                self.messages.pop()
            raise

        full_response = "".join(chunks)
        if full_response:
            self.messages.append({"role": "assistant", "content": full_response})

        return full_response

    def advance_role(self) -> bool:
        """Advance to the next role in the route.

        Returns True if there is a next role, False if the cycle is complete.
        """
        role = self.current_role
        for msg in reversed(self.messages):
            if msg["role"] == "assistant":
                self.role_outputs[role] = msg["content"]
                break

        self.current_role_index += 1
        return not self.is_cycle_complete

    def request_stop(self) -> None:
        """Request the current streaming response to stop."""
        self._stop_requested = True

    def reset_cycle(self, mode: str | None = None) -> None:
        """Reset the engine for a new work cycle."""
        # Kill the persistent CLI process so a fresh one starts
        if hasattr(self.client, "kill"):
            self.client.kill()
        self.messages.clear()
        self.role_outputs.clear()
        self.current_role_index = 0
        self.token_usage = {"input": 0, "output": 0}
        self.cost_usd = 0.0
        self.session_id = ""  # New session for new cycle
        self._stop_requested = False
        if mode is not None and mode != self.mode:
            self._load_mode(mode)

    def export_state(self) -> dict:
        """Export engine state for session persistence."""
        return {
            "mode": self.mode,
            "role_index": self.current_role_index,
            "route": list(self.route),
            "role_outputs": dict(self.role_outputs),
            "engine_messages": list(self.messages),
            "token_usage": dict(self.token_usage),
            "cost_usd": self.cost_usd,
            "session_id": self.session_id,
        }

    def restore_state(self, data: dict) -> None:
        """Restore engine state from a saved session.

        The client and prompt loader remain as-is; only conversation
        state is restored so ``--resume`` can continue the CLI session.
        """
        if data.get("mode") and data["mode"] != self.mode:
            self._load_mode(data["mode"])
        self.current_role_index = data.get("role_index", 0)
        self.role_outputs = data.get("role_outputs", {})
        self.messages = data.get("engine_messages", [])
        self.token_usage = data.get("token_usage", {"input": 0, "output": 0})
        self.cost_usd = data.get("cost_usd", 0.0)
        self.session_id = data.get("session_id", "")

    def available_modes(self) -> list[str]:
        """Return list of available mode IDs."""
        return self.loader.available_modes()

    # -- test-failure retry ---------------------------------------------------

    def retry_from_builder(self) -> bool:
        """Rewind the cycle to the builder_agent after a test failure.

        Stores the test_agent output in ``role_outputs`` and sets the
        current role back to ``builder_agent``.  Returns True if the
        rewind succeeded, False if builder_agent is not in the route.
        """
        # Save test output first
        role = self.current_role
        for msg in reversed(self.messages):
            if msg["role"] == "assistant":
                self.role_outputs[role] = msg["content"]
                break

        try:
            builder_idx = self.route.index("builder_agent")
        except ValueError:
            return False
        self.current_role_index = builder_idx
        return True

    # -- mode auto-detection --------------------------------------------------

    @staticmethod
    def thinking_budget_for_role(role_id: str) -> int:
        """Return the recommended thinking budget for a role."""
        return _ROLE_THINKING_BUDGETS.get(role_id, _DEFAULT_THINKING_BUDGET)

    @staticmethod
    def model_for_role(role_id: str) -> str:
        """Return the recommended model for a role ('auto' = user's choice)."""
        return _ROLE_MODEL_MAP.get(role_id, "auto")

    @staticmethod
    def max_tokens_for_role(role_id: str) -> int:
        """Return the max output tokens for a role."""
        return _ROLE_MAX_TOKENS.get(role_id, _DEFAULT_MAX_TOKENS)

    def _load_cycle_memory(self, max_entries: int = 10) -> str:
        """Load recent cycle summaries for Session Manager context."""
        import json
        mem_path = Path.home() / "agent_workspace" / ".cycle_memory.jsonl"
        if not mem_path.exists():
            return ""
        try:
            lines = mem_path.read_text(encoding="utf-8").strip().splitlines()
            recent = lines[-max_entries:]
            entries = []
            for line in recent:
                try:
                    entry = json.loads(line)
                    entries.append(
                        f"- [{entry.get('timestamp', '?')[:10]}] "
                        f"{entry.get('mode', '?')}: "
                        f"{entry.get('task', '?')[:80]} "
                        f"→ {entry.get('verdict', '?')} "
                        f"(retries={entry.get('retries', 0)}, "
                        f"${entry.get('cost_usd', 0):.2f})"
                    )
                except (json.JSONDecodeError, KeyError):
                    continue
            if entries:
                return (
                    "\n\n--- Cycle Memory (recent tasks) ---\n"
                    + "\n".join(entries)
                )
        except Exception:
            pass
        return ""

    def _build_test_handoff(self, user_task: str, prior_summary: str) -> str:
        """Build a test handoff with extracted acceptance criteria."""
        sm_output = self.role_outputs.get("session_manager", "")
        criteria = self.extract_acceptance_criteria(sm_output)
        gates = self.extract_stage_gates(sm_output)
        checklist = ""
        if criteria:
            items = "\n".join(f"  - [ ] {c}" for c in criteria)
            checklist = f"\n\nAcceptance criteria checklist:\n{items}\n"
        gate_list = ""
        if gates:
            items = "\n".join(f"  - [ ] {g}" for g in gates)
            gate_list = f"\n\nStage gate checklist:\n{items}\n"
        return (
            f"Validate the implementation for this task:\n\n{user_task}\n\n"
            f"Prior agent outputs:\n{prior_summary}\n"
            f"{checklist}{gate_list}\n"
            f"IMPORTANT: Confirm acceptance criteria with the user BEFORE running tests. "
            f"Use QUESTION: to list the criteria and ask if anything should be added or skipped.\n\n"
            f"Run `python -m pytest tests/ -v --tb=short` and verify each "
            f"criterion and gate above as PASS / FAIL / UNTESTED."
        )

    @staticmethod
    def extract_plan_field(session_manager_output: str, field_name: str) -> str:
        """Extract a ``**Field:** value`` entry from the Session Manager plan."""
        pattern = rf"^\*\*{re.escape(field_name)}:\*\*\s*(.+)$"
        match = re.search(pattern, session_manager_output, flags=re.MULTILINE)
        return match.group(1).strip() if match else ""

    @staticmethod
    def _extract_list_section(
        session_manager_output: str,
        headings: tuple[str, ...],
    ) -> list[str]:
        """Extract bullet/numbered items from a named markdown section."""
        items: list[str] = []
        in_section = False
        heading_tokens = tuple(h.lower() for h in headings)
        for line in session_manager_output.splitlines():
            stripped = line.strip()
            lowered = stripped.lower()
            if any(lowered.startswith(token) for token in heading_tokens):
                in_section = True
                continue
            if in_section:
                if stripped.startswith("###") or stripped.startswith("## "):
                    break
                match = re.match(r"^(?:\d+\.\s*|-\s*)(.*)", stripped)
                if match and match.group(1).strip():
                    items.append(match.group(1).strip())
        return items

    def _build_locked_plan_contract(self) -> str:
        """Summarize the Session Manager plan for downstream handoffs."""
        sm_output = self.role_outputs.get("session_manager", "")
        if not sm_output:
            return ""

        task = self.extract_plan_field(sm_output, "Task")
        goal_lock = self._extract_list_section(sm_output, ("### Goal lock",))
        scope = self._extract_list_section(sm_output, ("### Scope",))
        out_of_scope = self._extract_list_section(sm_output, ("### Out of scope",))
        criteria = self.extract_acceptance_criteria(sm_output)
        gates = self.extract_stage_gates(sm_output)

        lines = ["Locked plan contract:"]
        if task:
            lines.append(f"- Task: {task}")
        if goal_lock:
            lines.append("- Goal lock:")
            lines.extend(f"  - {item}" for item in goal_lock)
        if scope:
            lines.append("- Scope:")
            lines.extend(f"  - {item}" for item in scope)
        if out_of_scope:
            lines.append("- Out of scope:")
            lines.extend(f"  - {item}" for item in out_of_scope)
        if criteria:
            lines.append("- Acceptance criteria:")
            lines.extend(f"  - {item}" for item in criteria)
        if gates:
            lines.append("- Stage gates:")
            lines.extend(f"  - {item}" for item in gates)
        lines.append(
            "- Downstream rule: do not redefine this contract silently. Raise QUESTION if it is wrong."
        )
        return "\n".join(lines)

    def build_handoff_message(self, user_task: str) -> str:
        """Build a context-rich handoff message for the current role.

        Includes the original user task and role-specific instructions so
        each agent knows exactly what to do without user intervention.
        """
        role = self.current_role
        prior_summary = ""
        # Role-aware prior output limits — SM plan is critical for Builder
        _PRIOR_LIMITS = {
            "session_manager": 8000,
            "critic_agent": 6000,
            "runtime_agent": 6000,
            "reviewer_agent": 4000,
            "research_agent": 4000,
        }
        if self.role_outputs:
            parts = []
            for rid, output in self.role_outputs.items():
                limit = _PRIOR_LIMITS.get(rid, 2000)
                text = output[:limit]
                if len(output) > limit:
                    text += "\n... [truncated]"
                parts.append(f"### {rid}\n{text}")
            prior_summary = "\n\n".join(parts)

        # Load cycle memory for Session Manager
        _cycle_memory_ctx = ""
        if role == "session_manager":
            _cycle_memory_ctx = self._load_cycle_memory()
        locked_contract = self._build_locked_plan_contract()
        contract_block = f"\n\n{locked_contract}\n" if locked_contract else "\n"

        # Role-specific handoff instructions
        instructions = {
            "session_manager": (
                f"The user wants:\n\n{user_task}\n\n"
                f"IMPORTANT: Ask the user clarifying questions BEFORE creating the plan. "
                f"Use QUESTION: to ask about scope, constraints, and success criteria. "
                f"Only write the PLAN after the user responds.\n\n"
                f"When you do write the plan: lock the real goal, define the success oracle, "
                f"and break non-trivial work into small stage gates with exit evidence.\n"
                f"If the task needs external research, add RESEARCH_NEEDED: [topic]. "
                f"If not, add SKIP_RESEARCH."
                + (_cycle_memory_ctx or "")
            ),
            "chief_agent": (
                f"The user wants:\n\n{user_task}\n\n"
                f"Provide strategic direction as a CHIEF DIRECTIVE."
            ),
            "research_agent": (
                f"Research the following for this task:\n\n{user_task}\n\n"
                f"{locked_contract}\n\n"
                f"Prior agent outputs:\n{prior_summary}\n\n"
                f"IMPORTANT: Confirm your research questions with the user BEFORE searching. "
                f"Use QUESTION: to list your 2-3 planned research questions and ask if "
                f"anything else should be investigated.\n\n"
                f"Focus on actionable findings for the Builder. "
                f"Use WebSearch and WebFetch. Max 5 searches."
            ),
            "critic_agent": (
                f"Review the plan and/or changes for this task:\n\n{user_task}\n\n"
                f"{locked_contract}\n\n"
                f"Prior agent outputs:\n{prior_summary}\n\n"
                f"IMPORTANT: After your analysis, present your top findings to the user "
                f"with QUESTION: and ask which the Builder should prioritize.\n\n"
                f"Produce a structured REVIEW. Focus on critical and major issues. "
                f"Explicitly flag goal drift, weak success proxies, and missing stage gates."
            ),
            "runtime_agent": (
                f"Review the runtime implications for this task:\n\n{user_task}\n\n"
                f"{locked_contract}\n\n"
                f"Prior agent outputs:\n{prior_summary}\n\n"
                f"Produce a structured RUNTIME REVIEW. Check local vs cluster behavior and "
                f"whether the runtime-facing success metric is actually sufficient."
            ),
            "builder_agent": (
                f"Implement the following task:\n\n{user_task}\n\n"
                f"{locked_contract}\n\n"
                f"Prior agent outputs:\n{prior_summary}\n\n"
                f"IMPORTANT: Confirm your approach with the user BEFORE writing code. "
                f"Use QUESTION: to summarize your planned approach and the key trade-off "
                f"or decision you need confirmed.\n\n"
                f"Follow the PLAN from Session Manager. Address all critical/major "
                f"findings from Critic/Runtime/Research. Work through stage gates in order. "
                f"Do not silently substitute an easier proxy metric for the locked goal. "
                f"Run tests after implementation."
            ),
            "reviewer_agent": (
                f"Review the implemented changes for this task:\n\n{user_task}\n\n"
                f"{locked_contract}\n\n"
                f"Prior agent outputs:\n{prior_summary}\n\n"
                f"Review the actual diff, not the hypothetical plan. "
                f"Check correctness, regressions, and whether the Builder stayed aligned "
                f"with the locked goal instead of solving an easier proxy problem.\n"
                f"If you find ambiguous changes, ask the user with QUESTION: before finalizing."
            ),
            "test_agent": (
                self._build_test_handoff(user_task, prior_summary)
                + contract_block
            ),
        }
        return instructions.get(role, f"Continue with the task:\n\n{user_task}")

    @staticmethod
    def suggest_mode(user_message: str, current_mode: str = "quick") -> str | None:
        """Suggest a higher mode based on file paths mentioned in the message.

        Returns the suggested mode ID if escalation is warranted,
        or None if the current mode is sufficient.
        """
        current_rank = _MODE_RANK.get(current_mode, 0)
        decision = AgentEngine.recommend_task_route(user_message, current_mode)
        suggested = decision.get("mode", current_mode)
        if _MODE_RANK.get(suggested, 0) > current_rank:
            return suggested
        return None

    @staticmethod
    def recommend_task_route(user_message: str, current_mode: str = "dashboard") -> dict[str, Any]:
        """Recommend the cheapest capable mode for a task.

        Returns a dict with:
        - ``mode``: recommended mode
        - ``task_class``: dashboard / chemistry / coding / general
        - ``intent``: operate / research / question / change
        - ``confidence``: low / medium / high
        - ``reasons``: short rationale list
        - ``risk_flags``: reviewed / cluster / full flags when present
        """
        text = user_message or ""
        lower = text.lower()
        stripped = lower.strip()

        def _contains_any(keywords: tuple[str, ...]) -> list[str]:
            return [kw for kw in keywords if kw in lower]

        dashboard_hits = _contains_any(_DASHBOARD_KEYWORDS)
        chemistry_hits = _contains_any(_CHEMISTRY_KEYWORDS)
        code_change_hits = _contains_any(_CODE_CHANGE_KEYWORDS)
        code_question_hits = _contains_any(_CODE_QUESTION_KEYWORDS)
        code_hint_hits = _contains_any(_FILE_OR_CODE_HINTS)
        reviewed_hits = _contains_any(_REVIEWED_RISK_KEYWORDS)
        cluster_hits = _contains_any(_CLUSTER_RISK_KEYWORDS)
        full_hits = _contains_any(_FULL_RISK_KEYWORDS)

        for pattern, mode in _ESCALATION_PATTERNS:
            if re.search(pattern, text):
                if mode == "cluster":
                    cluster_hits.append(pattern)
                elif mode == "reviewed":
                    reviewed_hits.append(pattern)

        dashboard_score = len(dashboard_hits) * 2 + (3 if stripped.startswith("/") else 0)
        chemistry_score = len(chemistry_hits) * 2
        code_score = len(code_change_hits) * 2 + len(code_hint_hits)
        if re.search(r"\bdelfin/[\w/.\-]+", lower):
            code_score += 3
        if "?" in text:
            chemistry_score += 1 if chemistry_hits else 0
            code_score += 1 if code_hint_hits or code_change_hits else 0

        task_class = "general"
        intent = "question"
        confidence = "low"
        mode = current_mode or "dashboard"
        reasons: list[str] = []

        if dashboard_score >= max(2, chemistry_score, code_score):
            task_class = "dashboard"
            intent = "operate"
            mode = "dashboard"
            confidence = "high" if dashboard_score >= 3 else "medium"
            if dashboard_hits:
                reasons.append(f"dashboard signals: {', '.join(dashboard_hits[:3])}")
            if stripped.startswith("/"):
                reasons.append("slash-command style request")
        elif chemistry_score >= 2:
            task_class = "chemistry"
            if chemistry_hits:
                reasons.append(f"chemistry terms: {', '.join(chemistry_hits[:3])}")
            if code_change_hits:
                # Chemistry-informed code change → needs research context before building
                intent = "change"
                mode = "reviewed"
                confidence = "high" if chemistry_score >= 4 else "medium"
                reasons.append("chemistry code change needs research context")
            else:
                intent = "research"
                mode = "research"
                confidence = "high" if chemistry_score >= 4 else "medium"
                reasons.append("chemistry/scientific method question")
        else:
            task_class = "coding"
            is_question = bool("?" in text or code_question_hits)
            is_change = bool(
                code_change_hits
                or re.search(r"\b(test|bug|feature|refactor|implement|fix|patch|change|edit|add|update)\b", lower)
            )
            if is_question:
                intent = "question"
                mode = "solo"
                confidence = "medium" if code_score or code_question_hits else "low"
                reasons.append("read-only codebase question")
            else:
                intent = "change"
                mode = "quick"
                confidence = "medium"
                if full_hits:
                    mode = "full"
                    confidence = "high"
                elif cluster_hits:
                    mode = "cluster"
                    confidence = "high"
                elif reviewed_hits:
                    mode = "reviewed"
                    confidence = "high"
                if code_change_hits:
                    reasons.append(f"code-change verbs: {', '.join(code_change_hits[:3])}")
                if code_hint_hits:
                    reasons.append(f"code targets: {', '.join(code_hint_hits[:3])}")

        if task_class == "coding" and current_mode == "dashboard" and mode == "quick" and not code_change_hits and code_hint_hits:
            # Questions about code from dashboard should use solo, not a full cycle.
            mode = "solo"
            intent = "question"
            confidence = "medium"
            reasons.append("de-escalated to single-agent code exploration")

        return {
            "mode": mode,
            "task_class": task_class,
            "intent": intent,
            "confidence": confidence,
            "reasons": reasons[:4],
            "risk_flags": {
                "reviewed": bool(reviewed_hits),
                "cluster": bool(cluster_hits),
                "full": bool(full_hits),
            },
        }

    # -- context compaction ---------------------------------------------------

    def compact_for_next_role(self) -> None:
        """Clear conversation messages but preserve structured role outputs.

        Called between roles to prevent context window overflow.
        The prior role outputs (stored in ``role_outputs``) are injected
        into the next role's system prompt, so the conversation history
        is no longer needed.

        NOTE: The CLI process MUST be restarted because each role uses a
        different system prompt (``--append-system-prompt`` is set at
        process startup).  Keeping it alive would use the wrong prompt.
        """
        self.messages.clear()
        if hasattr(self.client, "kill"):
            self.client.kill()
        self.session_id = ""

    # -- acceptance criteria extraction ---------------------------------------

    @staticmethod
    def extract_acceptance_criteria(session_manager_output: str) -> list[str]:
        """Parse acceptance criteria from the Session Manager's PLAN output.

        Looks for the ``### Acceptance criteria`` section and extracts
        numbered items.
        """
        criteria: list[str] = []
        in_section = False
        for line in session_manager_output.splitlines():
            stripped = line.strip()
            if stripped.lower().startswith("### acceptance criteria"):
                in_section = True
                continue
            if in_section:
                if stripped.startswith("###") or stripped.startswith("## "):
                    break  # Next section
                # Match numbered items: "1. ...", "2. ...", "- ..."
                m = re.match(r"^(?:\d+\.\s*|-\s*)(.*)", stripped)
                if m and m.group(1).strip():
                    criteria.append(m.group(1).strip())
        return criteria

    @staticmethod
    def extract_stage_gates(session_manager_output: str) -> list[str]:
        """Parse stage gates from the Session Manager's PLAN output."""
        return AgentEngine._extract_list_section(
            session_manager_output,
            ("### Stage gates",),
        )

    @staticmethod
    def extract_status_field(agent_output: str) -> str:
        """Extract the canonical ``**status:**`` verdict if present."""
        match = re.search(
            r"^\*\*status:\*\*\s*(approve_with_risks|approve|reject)\s*$",
            agent_output or "",
            flags=re.IGNORECASE | re.MULTILINE,
        )
        return match.group(1).lower() if match else ""

    @staticmethod
    def extract_named_verdict(agent_output: str, label: str) -> str:
        """Extract a non-status verdict line such as ``**Verdict:** PASS``."""
        pattern = rf"^\*\*{re.escape(label)}:\*\*\s*(.+?)\s*$"
        match = re.search(pattern, agent_output or "", flags=re.IGNORECASE | re.MULTILINE)
        return match.group(1).strip() if match else ""

    @staticmethod
    def extract_check_statuses(agent_output: str, heading: str) -> list[tuple[str, str]]:
        """Extract ``item -> status`` rows from a markdown checklist section."""
        results: list[tuple[str, str]] = []
        in_section = False
        target = heading.lower()
        for line in (agent_output or "").splitlines():
            stripped = line.strip()
            if stripped.lower().startswith(target):
                in_section = True
                continue
            if in_section:
                if stripped.startswith("###") or stripped.startswith("## ") or stripped.startswith("**"):
                    break
                match = re.match(
                    r"^(?:\d+\.\s*|-\s*)(.*?)\s+[-—]\s+(PASS|FAIL|UNTESTED|DONE|PARTIAL|BLOCKED)\b",
                    stripped,
                    flags=re.IGNORECASE,
                )
                if match:
                    results.append((match.group(1).strip(), match.group(2).upper()))
        return results

    @staticmethod
    def evaluate_role_gate(role_id: str, output: str) -> tuple[str, str, str]:
        """Return a communication-gate decision for a completed role.

        Returns ``(action, gate_type, message)`` where action is one of:
        - ``continue``: safe to auto-advance
        - ``pause``: stop and ask the user to review/approve
        """
        text = output or ""
        status = AgentEngine.extract_status_field(text)

        # Session Manager: validate plan completeness before routing work.
        # Conversational responses (greetings, clarifications) are not plan
        # attempts — let the pipeline pause for plan approval instead.
        if role_id == "session_manager":
            if not AgentEngine.is_conversational(role_id, text):
                errors = AgentEngine.validate_role_output(role_id, text)
                if errors:
                    return (
                        "pause",
                        "schema",
                        "Session Manager plan is incomplete: " + "; ".join(errors[:4]),
                    )

        if role_id in {"research_agent", "critic_agent", "runtime_agent"}:
            if status == "reject":
                return (
                    "pause",
                    "risk",
                    "reported `status: reject`; review the findings before continuing.",
                )
            # approve_with_risks: auto-continue — pausing wastes tokens
            # and forces agents to restart from scratch

        if role_id == "builder_agent":
            stage_statuses = AgentEngine.extract_check_statuses(text, "**Stage gate status:**")
            crit_statuses = AgentEngine.extract_check_statuses(text, "**Acceptance criteria status:**")
            blocked = [name for name, state in stage_statuses + crit_statuses if state == "BLOCKED"]
            partial = [name for name, state in stage_statuses + crit_statuses if state == "PARTIAL"]
            if status == "reject":
                return (
                    "pause",
                    "partial",
                    "reported `status: reject`; implementation did not close safely.",
                )
            if blocked:
                return (
                    "pause",
                    "partial",
                    "left blocked items: " + "; ".join(blocked[:4]),
                )
            if partial:
                return (
                    "pause",
                    "partial",
                    "left partial items: " + "; ".join(partial[:4]),
                )

        if role_id == "reviewer_agent":
            verdict = AgentEngine.extract_named_verdict(text, "Verdict").upper()
            goal_lock = AgentEngine.extract_named_verdict(text, "Goal-lock check").upper()
            if "ISSUES" in verdict:
                return (
                    "pause",
                    "review",
                    "reported `Verdict: ISSUES`; builder should address review findings.",
                )
            if "ISSUES" in goal_lock:
                return (
                    "pause",
                    "goal-lock",
                    "flagged goal-lock issues; review whether the builder solved the correct problem.",
                )

        return ("continue", "", "")

    @staticmethod
    def is_conversational(role_id: str, output: str) -> bool:
        """Check whether a role output is conversational rather than structured.

        Returns True if the output does not look like a structured plan or
        report attempt — i.e. the agent responded conversationally (greeting,
        clarification, waiting for input) rather than producing work output.
        """
        text = (output or "").strip()
        if not text:
            return False
        upper = text[:500].upper()
        if role_id == "session_manager":
            # If there's no plan heading and no template markers, it's conversational
            return "## PLAN" not in upper and "### GOAL LOCK" not in upper
        return False

    @staticmethod
    def validate_role_output(role_id: str, output: str) -> list[str]:
        """Validate a role output against the required structured contract."""
        text = (output or "").strip()
        if not text:
            return ["empty output"]

        # Conversational responses are valid — the agent is waiting for a real task
        if AgentEngine.is_conversational(role_id, text):
            return []

        upper_head = text[:300].upper()
        if "QUESTION:" in text:
            return []
        if role_id in {"critic_agent", "reviewer_agent", "runtime_agent", "research_agent"}:
            if upper_head.startswith("SKIP") or "\nSKIP" in upper_head:
                return []

        def _missing(token: str, label: str | None = None) -> str:
            return f"missing {label or token}"

        def _contains(token: str) -> bool:
            return token.lower() in text.lower()

        errors: list[str] = []

        if role_id == "session_manager":
            required_tokens = [
                "## PLAN",
                "**Task:**",
                "**Class:**",
                "**Risk:**",
                "**Mode:**",
                "### Affected files",
                "### Goal lock",
                "### Scope",
                "### Out of scope",
                "### Acceptance criteria",
                "### Stage gates",
                "### Execution plan",
                "### Known risks",
                "**confidence:**",
                "**reason:**",
            ]
            for token in required_tokens:
                if not _contains(token):
                    errors.append(_missing(token))
            if not AgentEngine.extract_acceptance_criteria(text):
                errors.append("missing numbered acceptance criteria")
            if not AgentEngine.extract_stage_gates(text):
                errors.append("missing numbered stage gates")
            goal_lock = AgentEngine._extract_list_section(text, ("### Goal lock",))
            if len(goal_lock) < 3:
                errors.append("goal lock must include primary goal, success metric/oracle, and wrong proxy")

        elif role_id == "critic_agent":
            required_tokens = [
                "## REVIEW",
                "**Overall assessment:**",
                "**Risk level:**",
                "### Critical findings",
                "### Major findings",
                "### Moderate findings",
                "### Goal-drift and measurement risks",
                "**confidence:**",
                "**reason:**",
                "**status:**",
                "**key findings:**",
                "**open risks:**",
                "**recommended next step:**",
            ]
            for token in required_tokens:
                if not _contains(token):
                    errors.append(_missing(token))

        elif role_id == "runtime_agent":
            required_tokens = [
                "## RUNTIME REVIEW",
                "**Overall assessment:**",
                "### Local implications",
                "### Cluster implications",
                "### Failure modes",
                "### Environment assumptions",
                "### Goal-drift and measurement risks",
                "**confidence:**",
                "**reason:**",
                "**status:**",
                "**key findings:**",
                "**open risks:**",
                "**recommended next step:**",
            ]
            for token in required_tokens:
                if not _contains(token):
                    errors.append(_missing(token))

        elif role_id == "reviewer_agent":
            required_tokens = [
                "## CODE REVIEW",
                "**Files reviewed:**",
                "**Findings:**",
                "**Goal-lock check:**",
                "**Verdict:**",
                "**Summary:**",
                "**confidence:**",
                "**reason:**",
            ]
            for token in required_tokens:
                if not _contains(token):
                    errors.append(_missing(token))

        elif role_id == "test_agent":
            required_tokens = [
                "## TEST REPORT",
                "**Test command:**",
                "**Result:**",
                "**Acceptance criteria verification:**",
                "**Stage gate verification:**",
                "**New tests written:**",
                "**Regression check:**",
                "**confidence:**",
                "**reason:**",
                "**status:**",
                "**key findings:**",
                "**open risks:**",
                "**recommended next step:**",
            ]
            for token in required_tokens:
                if not _contains(token):
                    errors.append(_missing(token))

        return errors

    # -- pipeline progress tracking -------------------------------------------

    def pipeline_status(self) -> list[dict[str, str]]:
        """Return a list of pipeline steps with their completion status.

        Each entry: ``{"role": "session_manager", "status": "done|active|pending"}``
        """
        steps = []
        for i, role_id in enumerate(self.route):
            if i < self.current_role_index:
                status = "done"
            elif i == self.current_role_index:
                status = "active"
            else:
                status = "pending"
            steps.append({"role": role_id, "status": status})
        return steps

    # -- git helpers ----------------------------------------------------------

    def git_diff_stat(self) -> str:
        """Run ``git diff --stat`` in the repo and return the output."""
        try:
            result = subprocess.run(
                ["git", "diff", "--stat"],
                cwd=str(self.repo_dir),
                capture_output=True, text=True, timeout=10,
            )
            return result.stdout.strip()
        except Exception:
            return ""
