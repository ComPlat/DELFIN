"""DELFIN Agent Engine: orchestrates Claude conversations with role-based prompts."""

from __future__ import annotations

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
    # reviewed mode triggers
    (r"delfin/cli\.py", "reviewed"),
    (r"pipeline\.py", "reviewed"),
    (r"parallel_classic_manually\.py", "reviewed"),
    (r"delfin/api\.py", "reviewed"),
    (r"config\.py", "reviewed"),
]

# Mode escalation order (higher index = heavier mode)
_MODE_RANK = {"quick": 0, "reviewed": 1, "cluster": 2, "full": 3}

# -- Role-specific thinking budgets -------------------------------------------
# Builder needs the most thinking (implementation), review roles need moderate,
# planning roles need less.
# Per-role model routing: use the right model for each job.
# "auto" means use the user-selected model (no override).
_ROLE_MODEL_MAP: dict[str, str] = {
    "chief_agent": "sonnet",
    "session_manager": "sonnet",
    "critic_agent": "haiku",
    "runtime_agent": "haiku",
    "reviewer_agent": "haiku",
    "builder_agent": "auto",      # user's choice (often sonnet/opus)
    "test_agent": "sonnet",
    "solo_agent": "auto",
    "dashboard_agent": "haiku",   # only parses slash commands, cheap
    "research_agent": "sonnet",
}

_ROLE_THINKING_BUDGETS: dict[str, int] = {
    "chief_agent": 8000,
    "session_manager": 8000,
    "critic_agent": 16000,
    "runtime_agent": 16000,
    "builder_agent": 50000,
    "reviewer_agent": 16000,
    "test_agent": 16000,
    "research_agent": 16000,
    "solo_agent": 50000,
    "dashboard_agent": 4000,
}
_DEFAULT_THINKING_BUDGET = 10000


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
        api_key: str = "",
        model: str = "",
        mode: str = "quick",
        permission_mode: str = "",
        pack_dir: Path | None = None,
    ):
        self.repo_dir = Path(repo_dir)
        self.loader = PromptLoader(repo_dir=pack_dir)
        self.client = create_client(
            backend=backend, api_key=api_key, model=model,
            permission_mode=permission_mode,
            cwd=str(self.repo_dir),
        )
        self.backend = backend
        self.mode = mode
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
        on_permission_denied: Callable[[str], None] | None = None,
        on_thinking: Callable[[str], None] | None = None,
        thinking_budget: int = 0,
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
        on_permission_denied : callable, optional
            Called with (description) when a tool call was blocked.
        on_thinking : callable, optional
            Called with thinking text chunks as the model reasons.
        thinking_budget : int
            Extended thinking budget in tokens (0 = default/auto).

        Returns
        -------
        str
            The complete assistant response text.
        """
        self._stop_requested = False

        self.messages.append({"role": "user", "content": user_message})
        system_prompt = self._build_current_system_prompt(memory_context)

        chunks: list[str] = []
        try:
            for event in self.client.stream_message(
                system=system_prompt,
                messages=self.messages,
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
                    if on_tool_use:
                        on_tool_use(event.tool_name, event.tool_input)

                elif event.type == "permission_denied":
                    if on_permission_denied:
                        on_permission_denied(event.tool_name)

                elif event.type == "session_init":
                    # Capture session ID from CLI for persistence
                    if event.text:
                        self.session_id = event.text

                elif event.type == "message_start":
                    with self._lock:
                        self.token_usage["input"] += event.input_tokens

                elif event.type == "message_delta":
                    with self._lock:
                        self.token_usage["input"] += event.input_tokens
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

    def _build_test_handoff(self, user_task: str, prior_summary: str) -> str:
        """Build a test handoff with extracted acceptance criteria."""
        sm_output = self.role_outputs.get("session_manager", "")
        criteria = self.extract_acceptance_criteria(sm_output)
        checklist = ""
        if criteria:
            items = "\n".join(f"  - [ ] {c}" for c in criteria)
            checklist = f"\n\nAcceptance criteria checklist:\n{items}\n"
        return (
            f"Validate the implementation for this task:\n\n{user_task}\n\n"
            f"Prior agent outputs:\n{prior_summary}\n"
            f"{checklist}\n"
            f"Run `python -m pytest tests/ -v --tb=short` and verify each "
            f"criterion above as PASS / FAIL / UNTESTED."
        )

    def build_handoff_message(self, user_task: str) -> str:
        """Build a context-rich handoff message for the current role.

        Includes the original user task and role-specific instructions so
        each agent knows exactly what to do without user intervention.
        """
        role = self.current_role
        prior_summary = ""
        if self.role_outputs:
            parts = []
            for rid, output in self.role_outputs.items():
                # Include more from critical roles
                limit = 4000 if rid in ("critic_agent", "runtime_agent") else 2000
                text = output[:limit]
                if len(output) > limit:
                    text += "\n... [truncated]"
                parts.append(f"### {rid}\n{text}")
            prior_summary = "\n\n".join(parts)

        # Role-specific handoff instructions
        instructions = {
            "session_manager": (
                f"The user wants:\n\n{user_task}\n\n"
                f"Create a structured PLAN in the mandatory format. "
                f"Start by running `git diff --stat` to see the current state."
            ),
            "chief_agent": (
                f"The user wants:\n\n{user_task}\n\n"
                f"Provide strategic direction as a CHIEF DIRECTIVE."
            ),
            "critic_agent": (
                f"Review the plan and/or changes for this task:\n\n{user_task}\n\n"
                f"Prior agent outputs:\n{prior_summary}\n\n"
                f"Produce a structured REVIEW. Focus on critical and major issues."
            ),
            "runtime_agent": (
                f"Review the runtime implications for this task:\n\n{user_task}\n\n"
                f"Prior agent outputs:\n{prior_summary}\n\n"
                f"Produce a structured RUNTIME REVIEW. Check local vs cluster behavior."
            ),
            "builder_agent": (
                f"Implement the following task:\n\n{user_task}\n\n"
                f"Prior agent outputs:\n{prior_summary}\n\n"
                f"Follow the PLAN from Session Manager. Address all critical/major "
                f"findings from Critic/Runtime. Run tests after implementation."
            ),
            "test_agent": self._build_test_handoff(user_task, prior_summary),
        }
        return instructions.get(role, f"Continue with the task:\n\n{user_task}")

    @staticmethod
    def suggest_mode(user_message: str, current_mode: str = "quick") -> str | None:
        """Suggest a higher mode based on file paths mentioned in the message.

        Returns the suggested mode ID if escalation is warranted,
        or None if the current mode is sufficient.
        """
        current_rank = _MODE_RANK.get(current_mode, 0)
        best_mode = current_mode

        for pattern, mode in _ESCALATION_PATTERNS:
            if re.search(pattern, user_message):
                if _MODE_RANK.get(mode, 0) > _MODE_RANK.get(best_mode, 0):
                    best_mode = mode

        if _MODE_RANK.get(best_mode, 0) > current_rank:
            return best_mode
        return None

    # -- context compaction ---------------------------------------------------

    def compact_for_next_role(self) -> None:
        """Clear conversation messages but preserve structured role outputs.

        Called between roles to prevent context window overflow.
        The prior role outputs (stored in ``role_outputs``) are injected
        into the next role's system prompt, so the conversation history
        is no longer needed.
        """
        self.messages.clear()
        # Kill the CLI process so a fresh one starts with clean context
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
