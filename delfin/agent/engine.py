"""DELFIN Agent Engine: orchestrates model conversations with role-based prompts."""

from __future__ import annotations

import dataclasses
import json
import re
import subprocess
import threading
import time as _time
from pathlib import Path
from typing import Any, Callable

from . import agent_metrics as _agent_metrics
from .api_client import StreamEvent, create_client
from .prompt_loader import AVAILABLE_ON_DEMAND, PromptLoader


# -- Legacy mode name migration ------------------------------------------------
def _load_settings() -> dict:
    """User settings, defensively. Routing must never break on a bad file."""
    try:
        from delfin import user_settings
        return user_settings.load_settings() or {}
    except Exception:
        return {}


# Every name a session could have been saved under, mapped to a mode the
# user can actually pick. Two generations of retirement are folded in here,
# and the older one used to stop halfway: `default` mapped to `quick`,
# `release` to `full` -- onto the very modes that were later retired. A
# session saved months ago therefore migrated onto a mode missing from the
# picker: the control showed a value it could not offer, and the engine
# started a seven-role pipeline the user could neither see nor switch off.
#
# The multi-role pipelines are gone because what genuinely differs between
# sessions is which folder they work in and what they may reach from there.
# Pipeline in particular was never a mode in any sense the rest of the
# system used: it routed to solo_agent, ran the same tools, and differed
# only by a page of instructions. A page of instructions is a procedure, so
# it is a skill (`pipeline-build`) -- reachable from the coding mode
# without having to enter a mode you can then be stuck in, and paid for
# only when invoked.
#
# The session store keeps the identical table; when the two disagreed, the
# same saved session migrated differently depending on which one read it.
_LEGACY_MODE_MAP = {
    # First generation of names.
    "default": "solo",
    "high_risk": "solo",
    "runtime_cluster": "solo",
    "release": "solo",
    # The multi-role pipelines themselves.
    "quick": "solo",
    "reviewed": "solo",
    "tdd": "solo",
    "cluster": "solo",
    "full": "solo",
    "pipeline": "solo",
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

# Task-complexity-based model routing for solo mode (Feature 1).
# Simple tasks → cheaper model, complex → user's choice (often Opus).
_SIMPLE_TASK_PATTERNS: tuple[str, ...] = (
    "git status", "git log", "git diff", "git branch",
    "was ist", "what is", "zeig mir", "show me",
    "lies", "read ", "cat ", "grep ",
    "wie heißt", "how many", "wieviele",
    "welche datei", "which file", "where is",
    "list ", "zeige ", "display",
)

_COMPLEX_TASK_PATTERNS: tuple[str, ...] = (
    "refactor", "implement", "rewrite", "migrate",
    "architecture", "redesign", "optimize",
    "fix.*bug.*in.*multiple", "across.*files",
    "pipeline", "workflow", "integration",
)

# Adaptive thinking budget multipliers by task complexity.
# simple is 0.2 (not 0.4) so the OpenAI reasoning_effort mapping lands on
# "low" (<16k) — a greeting on gpt-5.x with effort=high used to burn 30+ s
# of internal reasoning for "Hallo".
_COMPLEXITY_THINKING_MULT: dict[str, float] = {
    "simple": 0.2,    # greetings / one-liners → minimal budget
    "moderate": 1.0,  # default
    "complex": 1.5,   # complex multi-file tasks → 150%
}

# Greetings / acknowledgements — no reasoning needed at all.
_GREETING_PATTERNS: tuple[str, ...] = (
    "hallo", "hi", "hey", "hello", "moin", "servus", "guten morgen",
    "guten tag", "danke", "thanks", "thank you", "ok", "okay", "cool",
    "super", "perfekt", "nice", "👍",
)

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
    # Dashboard agent dispatches UI commands ("open tab X", "set
    # functional Y", "submit job") — no deep reasoning required. The
    # OpenAI/Azure path maps the budget to reasoning_effort buckets
    # (<16k = low, <64k = medium, ≥64k = high); keeping the base at
    # 4 k means even the high-effort dropdown setting (mult ×2.0 =
    # 8 k) still resolves to reasoning_effort=low so Azure GPT-5.x
    # responds in < 5 s instead of burning 100+ s on hidden reasoning.
    "dashboard_agent": 4000,
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
_OPS_TOOL_PREFIX = "mcp__delfin-ops__"
_MCP_TOOL_PREFIXES = (_DOC_TOOL_PREFIX, _OPS_TOOL_PREFIX)
# KIT-Toolbox coding-agent tools (write_file, edit_file, multi_edit, bash) are
# emitted with this prefix by api_client.py. They are gated by the
# KitToolPermissions layer (sandbox + denylist + confirm-callback), not by the
# role whitelist below — so we always let them through when KIT permissions
# have been wired up.
_KIT_CODING_PREFIX = "mcp__kit-coding__"
# Source-inspection tools the dashboard agent must NOT use (its prompt hard-rule:
# "no read_file/grep_file/list_files/glob_files on DELFIN source" — it researches
# via search_docs and drives the UI via ACTION: slash-commands). These bypass the
# role whitelist below because they arrive mcp__-namespaced, so they are blocked
# explicitly by bare name for dashboard_agent (bug 20260708-092217: a dashboard
# agent read the whole delfin/manta source tree).
_DASHBOARD_SOURCE_DENY = frozenset({
    "read_file", "grep_file", "list_files", "glob_files",
})
# Roles that drive a conversation and therefore receive the per-turn steering
# blocks (open tasks, finished jobs, budget, late answers). The pipeline roles
# run one scripted step each and have no list of their own to keep.
_STEERING_ROLES = ("solo_agent", "dashboard_agent", "office_agent")
# Of those blocks, the ones worth re-sending when they change mid-turn. The
# project pin never changes within a turn, and the backend notice is a
# session-opening statement — re-sending either would cost tokens to say
# what the system prompt already says.
_MID_TURN_STEERING_KEYS = frozenset({
    "context_status", "open_tasks", "unmet_delegation", "unmet_tasklist",
    "finished_jobs", "finished_subagents", "budget", "answered",
})
# Blocks backed by a store that drains: their content is new by construction,
# so comparing it against the last delivery would suppress a real event that
# happens to read the same as an earlier one.
#
# These are also the blocks EVERY role gets. The role gate below exists for
# the blocks that are recomputed from state on every turn -- a scripted
# pipeline role keeps no task list, no budget of its own and no project pin,
# so those would cost it tokens to say nothing. A drain-backed block is the
# opposite: it is empty unless something actually finished, it is consumed
# exactly once by its own store, and a builder or test role waiting on a
# cluster job needs the completion just as much as a solo one. Gating them
# meant eight of the eleven roles were never told a job had ended.
_DRAINED_STEERING_KEYS = ("finished_jobs", "finished_subagents", "answered")
# Ceiling on mid-turn steering deliveries. A block that flaps (a task list
# edited repeatedly) must not be able to spend a turn's context on itself.
_MAX_STEERING_REFRESHES = 6
# How late the open-work-from-previous-sessions notice may still arrive,
# measured in messages present when the prompt is built (turn N sees
# 2N-1). Seven is four user turns: enough that a session which opens with
# a couple of greetings still gets it, few enough that it cannot surface
# in the middle of a conversation it has nothing to do with.
_FOREIGN_TASKS_TURN_CAP = 7
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
    "office_agent":    _SOLO_TOOLS,          # documents + shell, no chemistry
}



def _granted_dir_note(resolved) -> str:
    """What a turn in flight needs to hear about a new writable root."""
    return (
        f"[permissions] The user just granted access to {resolved} "
        f"mid-turn. It is writable from your next tool call. If you "
        f"stopped short of something because that path was outside the "
        f"workspace, it is inside now — but do not go back and redo work "
        f"that already succeeded."
    )

class AgentEngine:
    """Core orchestration engine for the DELFIN agent.

    Manages conversation state, role transitions, session persistence,
    and model communication via CLI or API backend.

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
        read_only_dirs: list[str] | None = None,
        confirm_write_dirs: list[str] | None = None,
        agent_workspace_dir: str = "",
        effort: str = "",
        kit_confirm_callback=None,
    ):
        self.repo_dir = Path(repo_dir)
        # DELFIN-bias gate: chemistry/cluster escalation patterns and
        # DELFIN-specific tools are only applied when the engine is
        # rooted in a DELFIN tree. Generic projects get a clean agent.
        try:
            from .api_client import _is_delfin_workspace
            self._is_delfin_workspace: bool = _is_delfin_workspace(self.repo_dir)
        except Exception:
            self._is_delfin_workspace = False
        self.loader = PromptLoader(repo_dir=pack_dir)
        self.effort = (effort or "").strip().lower()
        self.client = create_client(
            backend=backend, provider=provider, api_key=api_key,
            model=model, permission_mode=permission_mode,
            cwd=str(self.repo_dir), mcp_config=mcp_config,
            allowed_tools=allowed_tools, extra_dirs=extra_dirs,
            read_only_dirs=read_only_dirs, confirm_write_dirs=confirm_write_dirs,
            effort=self.effort,
            kit_confirm_callback=kit_confirm_callback,
        )
        # A model switch has consequences the caller used to have to
        # remember: re-size the compaction budget to the new model's real
        # context window, and tell a turn already in flight that it is
        # now talking to something else. Installed at construction rather
        # than per turn, because the dashboard switches between turns as
        # well -- same probe shape as should_stop and turn_cost_cap.
        try:
            self.client.on_model_switched = self._on_model_switched
        except Exception:
            pass
        self.backend = backend
        self.provider = provider
        AgentEngine._active_provider = provider  # class-level for static methods
        self.loader._active_provider = provider
        # Orientation must name the directory the agent actually works in:
        # the permissions workspace when the client has one, else the launch
        # dir. Never the DELFIN source tree that hosts the prompt pack.
        try:
            _ws = getattr(self.kit_permissions, "workspace", None)
            self.loader.workspace_root = Path(_ws) if _ws else Path(self.repo_dir)
            self.loader.is_delfin_workspace = self._is_delfin_workspace
        except Exception:
            self.loader.workspace_root = None
        # Inject-once gating of stable prompt sections is only sound on the
        # persistent claude CLI process (its first system prompt stays alive
        # across turns). Codex CLI spawns per turn and the chat-API backends
        # rebuild every request, so the loader re-injects there (its default).
        self.loader.stateful_backend = (backend == "cli" and provider == "claude")
        self.mode = mode
        self._agent_workspace_dir = agent_workspace_dir
        self.route: list[str] = []
        self.mode_description: str = ""
        self.current_role_index: int = 0
        self.messages: list[dict[str, Any]] = []
        self.role_outputs: dict[str, str] = {}
        self.compaction_summaries: dict[str, str] = {}
        # Mechanical verification state (verdict tool + evidence ledger),
        # copied from the backend client after each turn and snapshotted
        # per role. Gates consult THESE first; prose regex is only the
        # fallback — a formatting slip must never flip a reject into an
        # auto-continue.
        self.role_verdicts: dict[str, dict] = {}
        self.role_test_evidence: dict[str, list] = {}
        self._last_structured_verdict: dict | None = None
        self._last_test_evidence: list = []
        self.token_usage = {"input": 0, "output": 0, "cached": 0}
        self.cost_usd: float = 0.0
        # What this session's turns DELEGATED, kept apart from what they
        # spent themselves. cost_usd above has only ever counted the
        # parent model's own tokens; a delegated run is billed separately
        # and is written to a telemetry file that carries no parent,
        # session or turn identifier, so nothing joined the two and a turn
        # that handed the work to five sub-agents reported a few cents.
        # See _meter_delegate_spend for where the join is made instead.
        self._delegate_spend = _agent_metrics.DelegateLedger()
        # Whether that number could be measured at all, per turn. cost_usd
        # is a sum of MEASURED spend only: a turn on a model with no
        # published rate adds 0.0, which no USD gate can tell from a cheap
        # turn -- so the ceiling is weakest exactly where the price is
        # least known. The counters are what let the gates say that
        # instead of passing silently. A NON_BILLING turn is counted
        # apart, because its zero is a fact and calling it unmeasured
        # would be false in the other direction.
        self._measured_cost_turns: int = 0
        self._non_billing_turns: int = 0
        self._unpriced_turns: int = 0
        self._unmeasured_budget_notice_shown: bool = False
        # Budget notices already raised in the durable attention inbox,
        # by level. The prompt block below reaches the MODEL and the
        # per-turn breaker writes to the CHAT STREAM, so a run that hit
        # its ceiling at 03:00 left an empty inbox and a doctor reporting
        # PASS -- the two surfaces a user checks in the morning.
        self._budget_attention_levels: set = set()
        # Exact system prompt of the most recent turn (for bug reports).
        # Kept for the running process only -- it carries the injected
        # memory, so it is not written into every session file.
        self.last_system_prompt: str = ""
        # Its size, which IS persisted: the context estimate has to count
        # the system prompt, and a resumed session that starts at zero
        # reads far emptier than the conversation really is.
        self._system_prompt_chars: int = 0
        # Auto-compact configuration. Threshold is a fraction of the
        # context window; once an estimated turn would exceed it, the
        # engine summarises older messages before sending. Set
        # ``auto_compact_pct`` to 0 to disable.
        #
        # Threshold raised 0.80 → 0.95 (2026-05-11): weak MoE models like
        # qwen3.5 397B-A17b lose long-context coherence around token ~3500
        # and the original 0.80 threshold was compacting away the "I just
        # built X" thread mid-task, causing the agent to re-discover and
        # re-build its own work. 0.95 keeps the full thread until the
        # context is almost full, trading a bit more API cost for far
        # fewer context-loss double-work cycles.
        self.context_window_tokens: int = 100_000
        # Last authoritative input-token count the provider reported for a full
        # request (from message_start). Used as the ground-truth floor of the
        # compaction estimate so system prompt + tool schemas are never missed.
        self._last_input_tokens: int = 0
        # Run start timestamp for the session-scoped run budget.
        import time as _time_mod
        self._run_started_at: float = _time_mod.time()
        # Chars removed from self.messages by in-place trims since the floor
        # above was last set. The floor is a PRE-trim provider count, so the
        # estimator must credit these removals or the trim/hard-clear stop
        # conditions can never be reached once the floor is the binding term.
        self._trimmed_chars_since_floor: int = 0
        # Whether the current turn delivered a message_start (authoritative
        # input count). Backends that only report usage in message_delta get
        # an accounting-only fallback there.
        self._saw_message_start: bool = False
        # Whether this turn's compaction floor has been taken.
        # Separate from _saw_message_start, which also tracks a
        # zero-token event and feeds the accounting fallback.
        self._floor_captured_this_turn: bool = False
        # Resolved per active model below (real window for Ollama/local/cloud)
        # so compaction matches each model's true context — and small local
        # models stop overflowing while large ones stop being throttled to 100k.
        self._active_capabilities = None
        # Background: the cold-cache KIT probe (~5s) must not block the first
        # turn. The window upgrades from the 100k default once the probe lands.
        self._refresh_context_window(background=True)
        self.auto_compact_pct: float = 0.95
        self.last_compaction_info: dict[str, int] | None = None
        # A6 — last cost_usd snapshot at outcome time; the next outcome's
        # delta is (current cost_usd - this). Reset together with cost_usd
        # in reset_cycle() so a fresh cycle starts at delta = full spend.
        self._last_outcome_cost: float = 0.0
        # Session id. The CLI backend fills this from its stream events; the
        # OpenAI / KIT / Ollama backends emit no such event, so without a value
        # here every API session collapses onto the same empty-id bucket —
        # task_list then falls back to "all workspace tasks", so a brand-new
        # session shows the entire history (bug 20260616-183359). Mint a fresh
        # id for those backends so tasks (and session save/load) scope per
        # session. CLI stays "" and gets its real id from the stream.
        self.session_id: str = self._fresh_session_id()
        self._sync_task_session()
        # One-shot: the open-foreign-tasks notice fires only on a session's
        # first prompt build (re-armed with the session in reset_cycle).
        self._foreign_tasks_shown: bool = False
        # True once a saved conversation has been loaded into self.messages.
        # The notice above used to infer this from len(self.messages) > 1,
        # which cannot tell a restored history from the second turn of a
        # live session — so a session that opened with a greeting burned
        # the one shot on the greeting and the parked work was never named.
        self._history_restored: bool = False
        # Whether the message that opened the CURRENT turn was a greeting
        # and nothing else. Set per turn in stream_response.
        self._turn_is_bare_greeting: bool = False
        import uuid as _uuid_trace
        self._trace_id: str = _uuid_trace.uuid4().hex[:12]   # stable trace key
        self._trace_pending: list = []   # (tool, input, t0) awaiting its result
        self._prompt_session_serial: int = 1
        self._stop_requested = False
        # A stop belongs to ONE turn. Without an owner the flag was engine
        # state that any later caller could reset, including a caller whose
        # own turn was about to be refused -- which erased the stop of the
        # turn still running. _turn_serial numbers the turns, _turn_id is
        # the turn holding the gate (0 = none), and _stop_owner_turn is the
        # turn the current stop was requested for.
        self._turn_serial: int = 0
        self._turn_id: int = 0
        self._stop_owner_turn: int | None = None
        # Set by the last turn that ended without producing any answer
        # text (see the empty-turn branch in stream_response); None after
        # any turn that answered.
        self.last_empty_turn: dict | None = None
        self._lock = threading.Lock()
        # Claim-grounding guard state (see _enforce_claim_grounding).
        # _claim_guard_active is the reentrancy latch: True only while the
        # single forced correction turn runs, so that turn can never spawn
        # another one. _claim_guard_spent records that this user turn has
        # used up its one correction — the budget, and nothing else.
        # _claim_guard_corrected is the VERDICT: True only after a re-scan
        # of the correction found the claims actually resolved by evidence
        # read during it. The two used to be one flag, assigned before the
        # correction turn was even attempted, so "answer corrected" was
        # reported for corrections that raised a new claim and for ones
        # that only added hedges. The ledger flag records whether the
        # backend client exposes an observed-files ledger at all — without
        # one, ungrounded-location detection cannot judge and stays off.
        self._claim_guard_active: bool = False
        self._claim_guard_spent: bool = False
        self._claim_guard_corrected: bool = False
        self._last_observed_files: set[str] = set()
        self._last_turn_tools: list[str] = []
        self._observed_ledger_available: bool = False
        # Executed-command ledger (session-cumulative, capped): what this
        # session actually RAN. Evidence input for the functional-claim
        # guard — "it works now" is grounded only when the artifact was
        # exercised. Filled live from tool calls (see _note_exec_command),
        # cleared on a new work cycle like the rest of the session state.
        self._exec_commands_session: list[str] = []
        # Commands issued but not yet known to have run. Committed to the
        # ledger above only once their result shows they did.
        self._exec_pending: list[tuple[str, str]] = []
        # Live state: a per-turn snippet (Dashboard widget state, calc folder,
        # etc.) appended to the system prompt — keeps it OUT of the user
        # message body so it doesn't accumulate in self.messages history.
        self._live_state: str = ""
        # What the per-turn steering blocks said when this turn's system
        # prompt was built, keyed by block. The system prompt is frozen for
        # the whole tool loop, so a block that becomes true DURING the turn
        # (a task created in round 2, the budget crossing its threshold)
        # could not reach the model at all. Comparing against this map is
        # what makes a mid-loop re-delivery carry only what CHANGED.
        self._steering_delivered: dict[str, str] = {}
        # Mid-turn deliveries used so far this turn (see _drain_turn_steering).
        self._steering_refreshes: int = 0
        # One turn at a time per conversation. See the gate in
        # stream_response: two workers on one engine interleave the history
        # and the per-turn cost state.
        self._turn_gate = threading.Lock()
        self._turn_in_flight: bool = False
        # Project-directory pin: the dir of the FIRST project write this
        # session. Re-injected every turn so the agent keeps writing there and
        # doesn't drift into sibling folders over a long session (a framework
        # guarantee, not a prompt nudge). Reset when the conversation resets.
        self._project_dir: str = ""

        # Context engineering features (all opt-in, default off)
        self._context_tracker = None  # type: Any  # ContextUsageTracker
        self._distiller = None  # type: Any  # ContextDistiller
        self._progressive_disclosure = False
        self._init_context_features()

        self._load_mode(mode)

    def _init_context_features(self) -> None:
        """Initialize context engineering features (tracker, distiller).

        Distiller is enabled by default for the solo role (where long
        sessions accumulate big repo-map + briefing + memory context that
        actually benefits from haiku-class compression). Dashboard runs
        on a much smaller prompt budget and the distiller's per-call
        cost ($~0.04) outweighs the savings, so it stays off there.

        Other roles inherit the conservative default off.
        """
        try:
            from .context_tracker import ContextUsageTracker
            self._context_tracker = ContextUsageTracker()
            self.loader._context_tracker = self._context_tracker
        except Exception:
            pass
        try:
            from .context_distiller import ContextDistiller
            _provider = (getattr(self, "provider", "") or "claude").lower()
            # Only where the distiller's compression path actually works.
            # Its API call builds an Anthropic client from the ambient
            # environment; on any other provider that call fails and the
            # lossy line-level fallback becomes the ONLY path -- silently,
            # with no marker in the prompt and nothing written to the
            # elision store. Measured on the shipped pack, that fallback
            # keeps 18-35% of a section and cuts in file order, so the
            # last sections get zero budget: live state and session
            # environment were deleted in full, every time it ran.
            _enable = (
                self.mode in ("solo", "quick", "reviewed", "tdd",
                              "cluster", "full")
                and _provider in ("claude", "anthropic"))
            self._distiller = ContextDistiller(
                enabled=_enable, provider=_provider)
        except Exception:
            pass

    def set_progressive_disclosure(self, enabled: bool) -> None:
        """Enable/disable progressive context disclosure."""
        self._progressive_disclosure = enabled
        self.loader._progressive_disclosure = enabled

    def set_context_distillation(self, enabled: bool) -> None:
        """Enable/disable cheap-model context distillation."""
        if self._distiller is not None:
            self._distiller.enabled = enabled

    # -- KIT-Toolbox coding-agent permissions --------------------------------

    @property
    def kit_permissions(self):
        """KitToolPermissions instance bound to the client (None if not KIT)."""
        return getattr(self.client, "_permissions", None)

    def run_gated_bash(self, command: str) -> str:
        """Run one shell command through the AGENT's gate, for `!cmd`.

        The whole point is that it is not a way around anything: the same
        deny-list, the same secret scan, the same auto-allow list and the
        same confirmation callback the model's own bash calls go through.
        A shell escape that skipped them would let someone do, from the
        agent's prompt, exactly what the agent may not — and it would look
        like a convenience rather than a hole.

        Returns the tool result string. Raises RuntimeError when this
        backend carries no permissions object, because on those the shell
        tool is unavailable to the model too.
        """
        perms = self.kit_permissions
        if perms is None:
            raise RuntimeError(
                "this backend has no permission gate, so it has no shell")
        from .api_client import _doc_executor
        return _doc_executor.execute("bash", {"command": str(command or "")},
                                     perms)

    # -- what a turn delegated ---------------------------------------------

    def delegate_spend(self) -> "_agent_metrics.DelegateLedger":
        """This session's delegated spend: per turn, per session, and the
        background part that belongs to no turn.

        Handed out as the ledger rather than as three numbers so a caller
        cannot read the per-turn bucket where it meant the session's.
        """
        ledger = getattr(self, "_delegate_spend", None)
        if not isinstance(ledger, _agent_metrics.DelegateLedger):
            ledger = _agent_metrics.DelegateLedger()
            self._delegate_spend = ledger
        return ledger

    def _meter_delegate_spend(self) -> bool:
        """Route every delegated run through a counter that knows which
        turn asked for it.

        The delegate's own telemetry file records tokens and duration but
        carries no parent, session or turn field, so nothing in that file
        can say which turn a delegate belonged to and a timestamp window
        over it would be a guess -- worse, a shared home directory means
        the guess would also pick up another session's delegates. The
        runner the tool executor calls is the one place where a delegate
        and the turn that started it are the same call, so the join is
        made HERE, in memory, and the file is not consulted at all.

        Re-installed at the start of every turn because the client
        re-binds ``subagent_runner`` on every ``set_permissions``. It
        unwraps an existing meter before wrapping, so re-installing can
        never stack two counters over one runner and bill each delegate
        twice.

        Returns True when a meter is in place. False means this backend
        has no permissions object or no runner, i.e. it cannot delegate.
        """
        perms = self.kit_permissions
        if perms is None:
            return False
        runner = getattr(perms, "subagent_runner", None)
        if runner is None:
            return False
        inner = getattr(runner, "_delegate_meter_inner", None) or runner

        def _metered(*args, **kwargs):
            # Read at ENTRY, not at return: a background delegate comes
            # back long after the turn that started it, and the turn
            # running by then is exactly the one it must NOT be charged
            # to.
            started = self._turn_id if self._turn_in_flight else 0
            payload = inner(*args, **kwargs)
            try:
                self._book_delegate_spend(payload, started_turn=started)
            except Exception:
                pass
            return payload

        _metered._delegate_meter_inner = inner
        try:
            perms.subagent_runner = _metered
        except Exception:
            return False
        return True

    def _book_delegate_spend(self, payload, *, started_turn: int) -> None:
        """Charge one finished delegate, to exactly one place.

        The turn is identified by the pair (a turn is in flight, its id),
        read at the runner's entry and again at its return. Turn ids
        strictly increase and the turn gate allows only one turn in flight
        at a time, so an unchanged pair proves the SAME turn has been
        running throughout -- which is what makes a foreground delegate
        this turn's, and a background one, which returns after its turn
        ended, nobody's.

        Every delegate is booked into the session total and into exactly
        one of the per-turn or the background bucket, so no delegate can
        be counted twice however often the numbers are read afterwards.
        """
        ledger = self.delegate_spend()
        one = _agent_metrics.DelegateSpend()
        if not one.add_payload(payload,
                               provider=str(getattr(self, "provider", "")
                                            or "")):
            return
        with self._lock:
            same_turn = bool(
                started_turn
                and getattr(self, "_turn_in_flight", False)
                and getattr(self, "_turn_id", 0) == started_turn)
            (ledger.turn if same_turn else ledger.background).merge(one)
            ledger.session.merge(one)

    def set_kit_confirm_callback(self, callback) -> bool:
        """Bind/unbind the KIT-Toolbox confirmation callback at runtime.

        Returns True when the callback was attached, False when the active
        client does not have KIT permissions (e.g. provider != 'kit').
        """
        perms = self.kit_permissions
        if perms is None:
            return False
        perms.confirm_callback = callback
        return True

    def _on_model_switched(self, previous: str, current: str) -> None:
        """Everything a model change implies, done by the change itself.

        The context window was the caller's job before -- "call after
        every ``client.switch_model(...)``", written in a docstring and
        enforced by nothing, so a third caller that forgot it left the
        compaction budget sized for the model the agent had stopped
        talking to. And nobody told the model: a turn in flight carries
        on under a different one, the user reads one answer written by
        two, and the agent never learns its own capabilities changed.
        """
        try:
            self._refresh_context_window()
        except Exception:
            pass
        self._tell_the_running_turn(
            f"[model] This session switched from '{previous}' to "
            f"'{current}' just now. Your capabilities and context window "
            f"may differ from the ones you started this turn with. Carry "
            f"on from where you are -- do not restart the work -- and say "
            f"in your answer that the model changed partway.")

    def _tell_the_running_turn(self, note: str) -> None:
        """Put a fact about the run in front of a turn already in flight.

        The permission mode, a newly granted directory, a command the user
        just allowed for good -- all of them change what the agent may do,
        all of them are written from the dashboard thread, and all of them
        took effect at the gate without the model ever hearing about it.
        The gap is not academic: an agent that was refused a write two
        rounds ago has already told the user it cannot do the thing, and
        will not try again on its own.

        Best-effort by design. A client without the rail behaves exactly
        as it did before, and a permission change must never fail because
        a notice could not be delivered.
        """
        try:
            self.client.push_run_note(note)
        except Exception:
            pass

    def set_kit_permission_mode(self, mode: str) -> bool:
        """Switch the KIT permission mode at runtime.

        Accepts 'plan', 'default', 'acceptEdits', 'bypassPermissions'.
        Returns True on success, False if no KIT permissions are bound or
        the mode is unknown.
        """
        perms = self.kit_permissions
        if perms is None:
            return False
        if mode not in {"plan", "default", "diff_approval", "acceptEdits", "bypassPermissions"}:
            return False
        previous = str(getattr(perms, "mode", "") or "")
        perms.mode = mode
        # A turn already in flight is working to the rules it was told at
        # the top of it, while the gate is already enforcing these. That
        # gap is the whole problem: switched to plan, the model keeps
        # writing and every write comes back refused for reasons it was
        # never given; switched to acceptEdits, it keeps proposing and
        # asking, and the user watches it ask permission it now has.
        # Delivered on the same rail as steering and finished jobs, so it
        # lands between tool rounds and again if the turn ends without
        # one.
        if previous and previous != mode:
            self._tell_the_running_turn(
                    f"[permissions] The user switched this session from "
                    f"'{previous}' to '{mode}' just now, mid-turn. That is "
                    f"already in force for every tool call from here on. "
                    f"Work to the new rules from your next step -- do not "
                    f"redo what is done, and say in your answer that the "
                    f"rules changed partway.")
        return True

    def add_kit_workspace_dir(self, path, *,
                              persist: bool = True,
                              scope: str = "user") -> tuple[bool, str]:
        """Allow the KIT-Toolbox agent to operate inside an additional directory.

        Returns (ok, message). The path must exist and be a directory; on
        success the directory is added to ``extra_workspace_dirs``. When
        ``persist`` is True (default) the path is also written to the
        persisted KIT settings (``~/.delfin/settings.json`` for ``scope='user'``,
        ``<repo>/.delfin/settings.json`` for ``scope='repo'``) so it survives
        across sessions.
        """
        perms = self.kit_permissions
        if perms is None:
            return False, "KIT permissions are not active (provider != 'kit')."
        try:
            resolved = perms.add_extra_dir(path)
        except ValueError as exc:
            return False, str(exc)
        if persist:
            try:
                from . import kit_settings as _kit_settings
                _kit_settings.persist_extra_dir(
                    resolved, scope=scope, repo_dir=self.repo_dir
                )
            except Exception as exc:
                self._tell_the_running_turn(_granted_dir_note(resolved))
                return True, f"added (in-memory only — persist failed: {exc}): {resolved}"
        self._tell_the_running_turn(_granted_dir_note(resolved))
        return True, f"added: {resolved}"

    def list_kit_workspace_dirs(self) -> list[str]:
        """Return all workspace roots the KIT agent may touch (as strings)."""
        perms = self.kit_permissions
        if perms is None:
            return []
        return [str(p) for p in perms.all_workspace_roots()]

    def remove_kit_workspace_dir(self, path, *,
                                 scope: str = "user") -> tuple[bool, str]:
        """Drop ``path`` from the persisted KIT settings (in-memory dirs stay
        until next reload — recreate the engine to clear them entirely)."""
        try:
            from . import kit_settings as _kit_settings
            _kit_settings.remove_extra_dir(
                path, scope=scope, repo_dir=self.repo_dir
            )
        except Exception as exc:
            return False, f"persist failed: {exc}"
        return True, f"removed from persisted settings: {path}"

    def persist_kit_pattern(self, pattern: str, *,
                            kind: str = "allow",
                            scope: str = "user") -> tuple[bool, str]:
        """Write a bash regex pattern to the persisted KIT settings and apply
        it to the live permissions. ``kind`` is 'allow' or 'deny'."""
        perms = self.kit_permissions
        if perms is None:
            return False, "KIT permissions are not active (provider != 'kit')."
        if kind not in {"allow", "deny"}:
            return False, f"kind must be 'allow' or 'deny', got {kind!r}"
        if not pattern:
            return False, "pattern must be non-empty"
        try:
            from . import kit_settings as _kit_settings
            _kit_settings.persist_pattern(
                pattern, kind=kind, scope=scope, repo_dir=self.repo_dir
            )
        except Exception as exc:
            return False, f"persist failed: {exc}"
        # Apply to live perms so the next call benefits without restart.
        if kind == "allow":
            if pattern not in perms.bash_auto_allow_patterns:
                perms.bash_auto_allow_patterns = (
                    perms.bash_auto_allow_patterns + (pattern,)
                )
            self._tell_the_running_turn(
                f"[permissions] The user just allowed {pattern!r} for good, "
                f"mid-turn. It no longer needs a confirmation. If you gave "
                f"up on something because that command was refused, it is "
                f"available now.")
        else:
            self._tell_the_running_turn(
                f"[permissions] The user just denied {pattern!r} for good, "
                f"mid-turn. Do not attempt it again; find another route or "
                f"say what you cannot do without it.")
            if pattern not in perms.bash_deny_patterns:
                perms.bash_deny_patterns = (
                    perms.bash_deny_patterns + (pattern,)
                )
        return True, f"persisted {kind}-pattern: {pattern}"

    def persist_kit_default_mode(self, mode: str, *,
                                 scope: str = "user") -> tuple[bool, str]:
        """Persist the default permission mode and apply it live."""
        if mode not in {"plan", "default", "diff_approval", "acceptEdits", "bypassPermissions"}:
            return False, f"unknown mode: {mode!r}"
        try:
            from . import kit_settings as _kit_settings
            _kit_settings.persist_default_mode(
                mode, scope=scope, repo_dir=self.repo_dir
            )
        except Exception as exc:
            return False, f"persist failed: {exc}"
        self.set_kit_permission_mode(mode)
        return True, f"persisted default_mode: {mode}"

    def kit_settings_snapshot(self) -> dict:
        """Return the merged KIT settings as a dict (for UI / debugging)."""
        try:
            from . import kit_settings as _kit_settings
            return _kit_settings.load(repo_dir=self.repo_dir).to_dict()
        except Exception:
            return {}

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
        """Return current engine status for the UI.

        ``cost_usd`` keeps its meaning exactly: the parent model's own
        measured spend, cumulative for the session. What it never
        included -- the runs this session handed to sub-agents -- is
        reported beside it under its own names, plus a total that says so
        in its own. Folding the delegated half into ``cost_usd`` would
        make every existing reader silently start meaning something else.
        """
        _dl = self.delegate_spend()
        return {
            "mode": self.mode,
            "backend": self.backend,
            "provider": self.provider,
            "role": self.current_role,
            "role_index": self.current_role_index,
            "role_total": len(self.route),
            "input_tokens": self.token_usage.get("input", 0),
            "output_tokens": self.token_usage.get("output", 0),
            "cached_tokens": self.token_usage.get("cached", 0),
            "cost_usd": self.cost_usd,
            # Cumulative for the session, like cost_usd above, so a caller
            # that already differences against a turn-start baseline gets
            # the turn's share of the delegated spend the same way.
            "delegated_cost_usd": _dl.session.cost_usd,
            "delegated_input_tokens": _dl.session.input_tokens,
            "delegated_output_tokens": _dl.session.output_tokens,
            "delegate_count": _dl.session.count,
            # Delegates whose model has no published rate. They are NOT in
            # delegated_cost_usd; a reader that shows the figure without
            # this count is reporting an unpriced run as a free one.
            "delegates_unpriced": _dl.session.unpriced,
            # The part of the delegated total that belongs to no single
            # turn: a backgrounded delegate finishes after the turn that
            # started it, so charging it to a turn would charge whichever
            # turn happened to be running when the thread came back.
            "background_delegated_cost_usd": _dl.background.cost_usd,
            "total_cost_usd": self.cost_usd + _dl.session.cost_usd,
            "cycle_complete": self.is_cycle_complete,
            "session_id": self.session_id,
        }

    def _refresh_context_window(self, *, background: bool = False) -> None:
        """Size the compaction budget to the ACTIVE model's real context.

        Resolves the model's true window + capabilities (Ollama ``/api/show``,
        KIT/vLLM ``/v1/models`` ``max_model_len``, else a curated static
        table) and stashes the capabilities for tool/vision decisions. Falls
        back to the legacy 100k on any error so behaviour never regresses.
        Call after construction and after every ``client.switch_model(...)``.

        With ``background=True`` the (one-time, ~5s on a cold KIT cache) live
        probe runs on a daemon thread, so engine construction never blocks the
        first turn: the window stays at its current value (the 100k default)
        until the probe lands, then upgrades for the next turn. The 24h disk
        cache makes every later construction instant regardless. A first turn
        that races the probe just uses the safe 100k window + name-regex tool/
        vision fallback for that one turn — strictly better than a 5s stall.
        """
        if background:
            import threading
            threading.Thread(
                target=self._refresh_context_window,
                daemon=True, name="caps-prewarm").start()
            return
        try:
            model = getattr(self.client, "model", "") or ""
            base_url = ""
            inner = getattr(self.client, "client", None)
            if inner is not None:
                try:
                    base_url = str(getattr(inner, "base_url", "") or "")
                except Exception:
                    base_url = ""
            from .model_capabilities import resolve as _resolve_caps
            # Pass the provider key so the KIT/vLLM /v1/models probe authenticates
            # (else it 401s → static window). Doing it HERE — at construction /
            # switch_model, off the per-turn hot path — also warms the (authed)
            # capability cache, so the first stream_message turn reads it instantly
            # instead of paying a ~5s synchronous probe before its first token.
            api_key = getattr(self.client, "_api_key", "") or ""
            caps = _resolve_caps(self.provider, model, base_url, api_key=api_key)
            if caps and caps.context_window > 0:
                self.context_window_tokens = int(caps.context_window)
                self._active_capabilities = caps
        except Exception:
            # Keep whatever window is already set (100k default).
            pass

    def _compaction_status_line(self) -> str:
        """What the last compaction actually did, including when it did
        nothing and why.

        ``last_compaction_info`` is written in five shapes and three of
        them carry a ``note`` that says the run is in trouble: the
        summariser produced nothing and a deterministic digest went in
        instead; nothing could be summarised at all and the context is
        still over budget; the system prompt alone exceeds the window, so
        no amount of trimming can reach it.

        All three were rendered as ``compacted 0 msgs, saved ~0 tokens``
        — the same sentence a harmless compaction writes, and shorter
        than the one a real one writes. The worse the situation, the
        quieter the line. The two shapes that count something other than
        messages (a sliding-window trim, a tool-result edit) printed
        ``compacted ? msgs`` for work they had genuinely done.
        """
        last = getattr(self, "last_compaction_info", None) or None
        if not last:
            return "- Last compaction: (none this session)"
        kind = str(last.get("kind", "") or "compaction")
        parts: list[str] = []
        if last.get("messages_compacted"):
            parts.append(f"{last['messages_compacted']} msg(s) compacted")
        if last.get("messages_trimmed"):
            parts.append(f"{last['messages_trimmed']} msg(s) trimmed")
        if last.get("tool_results_cleared"):
            parts.append(
                f"{last['tool_results_cleared']} tool result(s) cleared")
        if last.get("tokens_saved"):
            parts.append(f"~{last['tokens_saved']} tokens saved")
        if last.get("pinned_kept"):
            parts.append(f"{last['pinned_kept']} pinned kept")
        body = ", ".join(parts) if parts else "nothing was compacted"
        line = f"- Last compaction ({kind}): {body}"
        note = str(last.get("note", "") or "").strip()
        if note:
            # The note is the whole reason the record was written. It is
            # what the agent has to act on — usually by making the PROMPT
            # smaller or delegating, neither of which it would think to do
            # from a zero.
            line += f"\n  ⚠️ {note}"
        return line

    def _build_context_status_block(self) -> str:
        """Compact self-monitoring block injected into the solo prompt.

        Tells the agent how close it is to the auto-compaction trigger so it
        can proactively delegate to subagents (and conserve its own window)
        before compaction fires. Empty string when usage can't be estimated.
        """
        try:
            n_msgs = len(self.messages or [])
            window = int(self.context_window_tokens or 0) or 100_000
            compact_pct = float(self.auto_compact_pct or 0.0)
            tokens = int(self._estimate_context_tokens())
        except Exception:
            return ""
        if window <= 0:
            return ""
        pct = tokens / window * 100.0
        last_line = self._compaction_status_line()
        warn = " — WARNING: nearing auto-compaction" if pct >= 80.0 else ""
        # Prompt-cache health: a collapsing hit rate is the earliest signal
        # that something destabilised the cached prefix (and doubled prefill
        # cost) — surface it where it gets seen every turn.
        cache_line = ""
        try:
            _in = int(self.token_usage.get("input", 0) or 0)
            _cached = int(self.token_usage.get("cached", 0) or 0)
            if _in > 0 and _cached > 0:
                cache_line = (
                    f"\n- Prompt cache: {_cached:,} of {_in:,} input tokens "
                    f"served from cache ({_cached / _in * 100.0:.0f}%)")
        except Exception:
            cache_line = ""
        # Pinned messages are exempt from every compaction stage — surface
        # the count so the agent knows part of the window is reserved.
        pin_line = ""
        try:
            _n_pins = len(self.pinned_indices())
            if _n_pins:
                pin_line = (
                    f"\n- Pins: {_n_pins} pinned message(s) excluded "
                    f"from compaction")
        except Exception:
            pin_line = ""
        return (
            "# Context status (auto-injected each turn)\n"
            f"- Compaction trigger: {compact_pct*100:.0f}% of window "
            f"(gentle trim from 70%)\n"
            f"- Current usage: {n_msgs} msgs, ~{tokens:,} tokens "
            f"({pct:.1f}% of {window:,}){warn}\n"
            f"{last_line}"
            f"{cache_line}"
            f"{pin_line}"
        )

    def _build_open_tasks_block(self) -> str:
        """Per-turn reminder of OPEN tasks so the agent doesn't forget to mark
        finished steps ``completed`` (Jerome 2026-06-13: "sometimes forgets to
        actually check off done tasks"). Read-only; empty string when there are
        no open tasks. Mirrors the session filter of the ``task_list`` tool so
        the agent sees exactly the list it manages. The rendered list is capped
        so a long backlog can't bloat the prompt (the remainder count still
        shows).
        """
        try:
            perms = self.kit_permissions
            if perms is None:
                return ""
            # Plan mode is read-only: the agent must present a plan and stop for
            # approval, not self-drive a task list. Suppressing the reminder
            # removes the "keep working the list" pressure while planning
            # (part of the fix for bug 20260708-092217).
            if getattr(perms, "mode", "") == "plan":
                return ""
            from delfin.agent.agent_tasks import (
                OPEN_STATUSES, get_store, resolve_session_scope,
            )
            sid = getattr(perms, "task_session_id", "") or ""
            # The shared resolver, so this block and the user's ticker
            # describe the same store the same way for the same id.
            tasks = get_store(perms.workspace).list(
                session_id=resolve_session_scope(sid), with_seq=True)
        except Exception:
            return ""
        open_tasks = [
            t for t in (tasks or []) if t.get("status") in OPEN_STATUSES
        ]
        if not open_tasks:
            return ""
        in_prog = [t for t in open_tasks if t.get("status") == "in_progress"]
        pending = [t for t in open_tasks if t.get("status") == "pending"]
        blocked = [t for t in open_tasks if t.get("status") == "blocked"]
        lines = ["# Open tasks (auto-injected each turn — keep this list honest)"]
        shown = 0
        _CAP = 8
        def _label(t: dict) -> str:
            # Show the session-relative number to the user; keep the global id
            # in parentheses because task_update/task_get key off it.
            _seq = t.get("seq")
            _id = t.get("id")
            return f"task {_seq} (id {_id})" if _seq is not None else f"#{_id}"
        for t in in_prog:
            lines.append(f"- [in_progress] {_label(t)} {str(t.get('subject', ''))[:80]}")
            shown += 1
        # Blocked before pending: it is the bucket that needs a decision
        # rather than a next step, and burying it under the backlog is
        # how a task waiting on the user goes unmentioned for a day.
        for t in blocked:
            if shown >= _CAP:
                break
            reason = str(t.get("blocked_reason", ""))[:60]
            lines.append(
                f"- [blocked]     {_label(t)} {str(t.get('subject', ''))[:80]}"
                + (f" — waiting on {reason}" if reason else ""))
            shown += 1
        for t in pending:
            if shown >= _CAP:
                break
            lines.append(f"- [pending]     {_label(t)} {str(t.get('subject', ''))[:80]}")
            shown += 1
        remainder = len(open_tasks) - shown
        if remainder > 0:
            lines.append(f"- … +{remainder} more open task(s)")
        if in_prog:
            lines.append(
                "→ If any in_progress task above is actually finished, call "
                "task_update(task_id, status='completed') NOW — don't batch it "
                "to the end."
            )
        if blocked:
            lines.append(
                "→ A blocked task is not progress: either the block is "
                "resolved now, or say so to the user."
            )
        return "\n".join(lines)

    def _build_finished_jobs_block(self) -> str:
        """Event-driven completion notice for background work.

        Drains jobs that finished since the previous turn (persistent
        bash-job registry + agent-watched SLURM jobs) into a compact block,
        so the model reacts to a finished multi-hour calculation on its next
        turn instead of babysitting it with blocking status polls.
        Best-effort; empty string when nothing finished."""
        try:
            perms = self.kit_permissions
            ws = getattr(perms, "workspace", None) if perms else None
            if not ws:
                return ""
        except Exception:
            return ""
        events: list[str] = []
        rendered: list[dict] = []
        try:
            # The sweep, not the single-folder drain: a job belongs to the
            # workspace it was started in, which is not always the one the
            # session is in now.
            from delfin.agent.bash_jobs import drain_all_finished_events
            for ev in drain_all_finished_events(ws) or []:
                rc = ev.get("exit_code")
                state = "ok" if rc == 0 else (
                    f"exit {rc}" if rc is not None else "finished (exit unknown)")
                if ev.get("timed_out"):
                    state = "killed at its wall-clock cap"
                elif ev.get("children_running"):
                    state = f"{state}, children still running"
                tail = (ev.get("stderr_tail") if rc not in (0, None)
                        else ev.get("stdout_tail")) or ""
                tail = " ".join(str(tail).split())[:200]
                submitted = ev.get("watched_slurm_jobs") or []
                events.append(
                    f"- bash job {ev.get('job_id')} [{state}, "
                    f"{ev.get('runtime_s', 0):.0f}s] "
                    f"{str(ev.get('command', ''))[:100]}"
                    + (f" → {tail}" if tail else "")
                    + (f" (submitted cluster job(s) {', '.join(submitted)}; "
                       f"now watched)" if submitted else ""))
                rendered.append(ev)
        except Exception:
            pass
        try:
            from delfin.agent.job_monitor import check_agent_jobs
            for ev in check_agent_jobs(ws) or []:
                sig = ", ".join(ev.get("signatures") or [])
                degraded = str(ev.get("degraded") or "")
                events.append(
                    f"- {ev.get('kind', 'job')} {ev.get('job_id')} "
                    f"[{ev.get('state', '?')}] "
                    f"{str(ev.get('description', ''))[:80]}"
                    + (f" — signatures: {sig}" if sig else "")
                    + (f" — {degraded}" if degraded else ""))
        except Exception:
            pass
        try:
            # A registry cycle that ran without its cross-process lock is
            # reported here rather than nowhere. It belongs beside the
            # completions because those are what it may have damaged: the
            # lost update is an acknowledged flag or an exit code, which
            # shows up as a job id nothing can address or a completion
            # announced a second time.
            from delfin.agent.bash_jobs import take_lock_timeouts
            unlocked = take_lock_timeouts()
            if unlocked:
                events.append(
                    f"- NOTE: {len(unlocked)} job-registry write(s) went "
                    "through WITHOUT the cross-process lock (another process "
                    "held it past the wait). A record may have been "
                    "overwritten: an unknown job id or a completion you have "
                    "already seen has this as its cause, not the job.")
        except Exception:
            pass
        if not events:
            return ""
        block = "\n".join(
            ["# Background jobs finished since your last turn "
             "(act on these results now)"] + events)
        # Phase two: each event that actually reached the block being
        # returned can now become a permanent acknowledgement. Confirming
        # here rather than inside the drain is what makes a turn that dies
        # BEFORE this point re-deliver instead of losing the only word a
        # twelve-hour calculation ever gets.
        try:
            from delfin.agent.bash_jobs import confirm_finished_events
            confirm_finished_events(rendered)
        except Exception:
            pass
        return block

    # How much of a finished delegate's report the push carries. The report
    # IS the deliverable, so a 200-character tail (right for a job's stdout)
    # would deliver the fact and drop the substance; the full text can be
    # arbitrarily long, so the rest waits behind subagent_result.
    _SUBAGENT_REPORT_CHARS = 1200

    def _build_finished_subagents_block(self) -> str:
        """Event-driven completion notice for background DELEGATES.

        Background bash jobs have had a completion push since the events
        existed; background sub-agents had nothing. The only way the model
        learned a delegate had finished was to spend a whole tool round on
        ``subagent_result(sa_id)`` — against per-model round budgets of 10
        to 50, and against a tool result that ends "Continue other work
        meanwhile". A parent that ended its turn instead (the normal "I
        started it, I'll wait" ending) never took another turn: no timer,
        no callback and no autonomous turn exists to wake it, so the
        report was not late, it was never relayed at all.

        Drains the reports that ended since the previous turn, in the same
        shape as ``_build_finished_jobs_block`` and with the same
        exactly-once guarantee — the store claims each report for exactly
        one caller, so a completion cannot be announced twice. Empty
        string when nothing finished; best-effort, never raises."""
        try:
            from delfin.agent.subagents import drain_finished_subagents
            events = drain_finished_subagents() or []
        except Exception:
            return ""
        lines: list[str] = []
        for ev in events:
            if not isinstance(ev, dict):
                continue
            sa_id = str(ev.get("sa_id") or "?")
            head = (f"- sub-agent {ev.get('subagent_type') or '?'} [{sa_id}] "
                    f"{str(ev.get('description') or '')[:80]}")
            if str(ev.get("status")) == "died":
                lines.append(
                    f"{head} — ENDED WITHOUT A REPORT: "
                    f"{str(ev.get('error') or 'no reason recorded')[:200]}. "
                    "This work was not done; start it again or do it here.")
                continue
            report = str(ev.get("final_text") or "").strip()
            if len(report) > self._SUBAGENT_REPORT_CHARS:
                report = (report[:self._SUBAGENT_REPORT_CHARS]
                          + f"\n  […] full report: subagent_result(sa_id="
                            f"'{sa_id}')")
            if report:
                # The second route the same prose takes into this context.
                # A background delegate's report is pushed here rather than
                # returned through the tool result, so marking only the
                # tool path would leave the push unmarked — and the push is
                # the one that arrives unasked, between rounds, while the
                # parent is in the middle of something else.
                from .subagents import mark_delegate_text
                report = mark_delegate_text(report)
            err = str(ev.get("error") or "").strip()
            lines.append(head + (f" — run error: {err[:160]}" if err else ""))
            if report:
                lines.append(report)
            verification = ev.get("verification")
            notice = (str(verification.get("notice") or "").strip()
                      if isinstance(verification, dict) else "")
            if notice:
                lines.append(notice[:400])
            worktree = ev.get("worktree")
            warn = (str(worktree.get("warning") or "").strip()
                    if isinstance(worktree, dict) else "")
            if warn:
                lines.append(f"  ⚠️ {warn[:300]}")
        if not lines:
            return ""
        return "\n".join(
            ["# Background sub-agents that finished since your last turn "
             "(their reports are below — relay what matters to the user; "
             "calling subagent_result again only repeats what you can "
             "already read here)"] + lines)

    # Explicit delegation requests, in the languages the dashboard is used
    # in. Matching is deliberately narrow: only an explicit mention of
    # sub-agents or delegation counts, never a generic "split this up".
    _DELEGATION_REQUEST_RE = re.compile(
        r"(?i)\b(?:sub[\s_-]?agent(?:en|s)?|subagent(?:en|s)?|delegier\w*|"
        r"delegate|delegation)\b")

    def _delegation_was_requested(self) -> bool:
        """True when a user message in this session explicitly asked for
        sub-agent delegation. Never raises."""
        try:
            for m in self.messages:
                if str(m.get("role")) != "user":
                    continue
                if self._DELEGATION_REQUEST_RE.search(str(m.get("content", ""))):
                    return True
        except Exception:
            return False
        return False

    def _build_unmet_delegation_block(self) -> str:
        """Per-turn reminder while an explicit delegation request is open.

        A static prompt rule proved too weak: in a field run the user
        named five areas to delegate, the agent used no sub-agent and did
        not say why. A live reminder that disappears the moment a
        sub-agent runs (or the model states why it is unsuitable) puts
        the request where the model actually looks.
        """
        try:
            if getattr(self, "_delegation_satisfied", False):
                return ""
            if not self._delegation_was_requested():
                return ""
            used = any("subagent" in str(t or "")
                       for t in (getattr(self, "_session_tool_names", None) or ()))
            if used:
                self._delegation_satisfied = True
                return ""
            return (
                "# Open request: delegation\n"
                "The user explicitly asked for sub-agents and none has run "
                "yet. Delegate the next self-contained piece with "
                "`subagent` (independent pieces in ONE message run "
                "concurrently), or state in one sentence why this work "
                "cannot be split — silently doing everything yourself "
                "ignores the request."
            )
        except Exception:
            return ""

    # An explicit instruction to open a task list / roadmap up front.
    _TASKLIST_REQUEST_RE = re.compile(
        r"(?i)\b(?:task_create|task[\s_-]?list|taskliste|aufgabenliste|"
        r"roadmap)\b")

    def _build_unmet_tasklist_block(self) -> str:
        """Per-turn reminder while an explicit request for a task list is
        open.

        Same reasoning as the delegation reminder: an instruction stated
        once at the start of a long execution turn loses against the work
        in front of the model. Disappears as soon as a task exists.
        """
        try:
            if getattr(self, "_tasklist_satisfied", False):
                return ""
            asked = any(
                self._TASKLIST_REQUEST_RE.search(str(m.get("content", "")))
                for m in self.messages if str(m.get("role")) == "user")
            if not asked:
                return ""
            used = any("task_create" in str(t or "") for t in
                       (getattr(self, "_session_tool_names", None) or ()))
            if not used:
                # A tool-name ledger is per session; the task store is
                # not. After a resume, or after adopting a task from a
                # previous session, the list the user asked for already
                # exists and the reminder was still demanding it be
                # created. Ask the store, which is what "a task exists"
                # actually means.
                try:
                    from delfin.agent.agent_tasks import (
                        get_store, resolve_session_scope,
                    )
                    perms = self.kit_permissions
                    if perms is not None:
                        used = bool(get_store(perms.workspace).list(
                            session_id=resolve_session_scope(
                                getattr(perms, "task_session_id", ""))))
                except Exception:
                    pass
            if used:
                self._tasklist_satisfied = True
                return ""
            return (
                "# Open request: task list\n"
                "The user explicitly asked for the roadmap to be opened with "
                "`task_create` before execution starts, and no task exists "
                "yet. Create one task per planned step now, then begin task 1 "
                "with a real action — the list is the progress record the "
                "user reads, not optional bookkeeping."
            )
        except Exception:
            return ""

    def _build_open_foreign_tasks_block(self) -> str:
        """One-shot notice — first SUBSTANTIVE prompt build of a session —
        of open tasks left behind by PREVIOUS sessions of this workspace.

        The task store survives across sessions but every listing surface
        filters to the current session id, so a fresh session would
        otherwise never see the work a paused one parked (it dies
        silently). Deliberately a summary with explicit opt-in adoption
        (task_adopt / task_list(all_sessions=true)): session scoping is a
        leak fix (bug 20260616-183359) and lists are never auto-merged.
        Best-effort; empty string when nothing foreign is open.

        A greeting is not a substantive turn. The trigger used to be the
        session's first prompt build outright, so opening with "Hallo"
        handed the model a backlog it had no reason to act on — and the
        model, having been given a list and nothing else to do, read it
        out. The shot is not consumed there: the work is still parked, and
        the first real request is when naming it helps.
        """
        if getattr(self, "_foreign_tasks_shown", False):
            return ""
        if getattr(self, "_turn_is_bare_greeting", False):
            # Not consumed — deliberately. Skipping a greeting must not
            # cost the session the one notice it gets.
            return ""
        if (getattr(self, "_history_restored", False)
                or len(self.messages) > _FOREIGN_TASKS_TURN_CAP):
            # A restored conversation already carries its own context, and
            # a notice that arrives deep into a live session reads as a
            # non sequitur. Either way: never show it late.
            self._foreign_tasks_shown = True
            return ""
        try:
            perms = self.kit_permissions
            ws = getattr(perms, "workspace", None) if perms else None
            sid = str(getattr(perms, "task_session_id", "") or "") if perms else ""
            if not ws or not sid:
                # Empty session id → task_list already shows every
                # workspace task; nothing is invisible.
                return ""
            from delfin.agent.agent_tasks import open_foreign_tasks
            summary = open_foreign_tasks(ws, sid)
        except Exception:
            return ""
        self._foreign_tasks_shown = True
        count = int(summary.get("count", 0) or 0)
        if count <= 0:
            return ""
        oldest = int(summary.get("oldest_age_days", 0) or 0)
        # "oldest 0 days" is what the arithmetic gives for work parked an
        # hour ago, and it reads as a missing value rather than as recency.
        age = f"oldest {oldest} days" if oldest else "all from today"
        lines = [
            "# Open work from previous sessions",
            f"- {count} open task(s), {age} — adopt with "
            "task_adopt(id) or list via task_list(all_sessions=true); "
            "ignore if obsolete.",
        ]
        for t in summary.get("tasks", [])[:5]:
            lines.append(
                f"- id {t.get('id')} [{t.get('status')}] "
                f"{str(t.get('subject', ''))[:80]} ({t.get('age_days', 0)}d)")
        return "\n".join(lines)

    def _build_answered_attention_block(self) -> str:
        """Late answers to parked questions/confirms (attention inbox).

        A question/confirm that timed out unattended was parked in
        ~/.delfin/attention_inbox.jsonl; once the user resolves it, the
        answer is injected here exactly once — mirroring
        _build_finished_jobs_block."""
        try:
            perms = self.kit_permissions
            sid = str(getattr(perms, "task_session_id", "") or "") if perms else ""
        except Exception:
            sid = ""
        lines: list[str] = []
        try:
            from delfin.agent.attention import drain_resolved
            for ev in drain_resolved(sid) or []:
                ans = ev.get("answer")
                ans_txt = ", ".join(map(str, ans)) if isinstance(ans, list) \
                    else str(ans or "").strip()
                what = str(ev.get("kind", "question")).replace("_pending", "")
                lines.append(
                    f"- The user answered your earlier {what} "
                    f"\"{str(ev.get('title', ''))[:120]}\": {ans_txt or '(no text)'}")
        except Exception:
            return ""
        if not lines:
            return ""
        return "\n".join(
            ["# Answers to requests that previously timed out "
             "(act on them now — do not re-ask)"] + lines)

    _MUTATE_TOOLS_FOR_PIN = frozenset({
        "write_file", "edit_file", "multi_edit", "apply_patch", "notebook_edit",
    })

    def _maybe_pin_project_dir(self, tool_name: str, tool_input: Any) -> None:
        """Record the directory of the FIRST project write this session so the
        prompt can re-pin the agent there every turn (anti-drift). First write
        wins; never overwritten until the conversation resets."""
        if self._project_dir:
            return
        name = (tool_name.rsplit("__", 1)[-1]
                if isinstance(tool_name, str) and tool_name.startswith("mcp__")
                else tool_name)
        if name not in self._MUTATE_TOOLS_FOR_PIN:
            return
        path = ""
        try:
            data = json.loads(tool_input) if isinstance(tool_input, str) else tool_input
            if isinstance(data, dict):
                path = str(data.get("file_path") or data.get("path")
                           or data.get("notebook_path") or "")
        except Exception:
            m = re.search(r'"(?:file_path|path|notebook_path)"\s*:\s*"([^"]+)"',
                          tool_input if isinstance(tool_input, str) else "")
            path = m.group(1) if m else ""
        path = path.strip().rstrip("/")
        if not path:
            return
        # Directory of the first project write ("." == workspace root).
        self._project_dir = path.rsplit("/", 1)[0] if "/" in path else "."

    def _build_project_dir_block(self) -> str:
        """High-salience per-turn reminder pinning the project to its first-write
        directory — keeps a long session from drifting into sibling folders."""
        d = getattr(self, "_project_dir", "")
        if not d:
            return ""
        if d == ".":
            where, rule = ("the workspace ROOT",
                           "keep writing this project's files at the workspace root")
        else:
            where, rule = (f"`{d}/`", f"put EVERY further file of this project under `{d}/`")
        return (f"📁 PROJECT LOCATION (locked for this task): {where}. You already "
                f"started the project here — {rule}. Do NOT create a sibling "
                f"folder or switch to a different directory name mid-task; that "
                f"splits the project and breaks imports/tests.")

    def _scope_memory_domain(self, role: str) -> str:
        """Tell the memory store which kind of work this turn is.

        Recall is gated on the domain of the role that runs, and every
        stored memory carries the domain of the folder it was written for.
        Registering the folder here is what connects the two: the write
        path reaches the store through the permission workspace and has no
        other way to learn that this folder is where administrative work
        happens. Returns the domain; best-effort, never raises.
        """
        try:
            from .memory_store import (
                DOMAIN_OFFICE, domain_for_role, register_office_workspace,
            )
            domain = domain_for_role(role, getattr(self, "mode", "") or "")
            if domain == DOMAIN_OFFICE:
                workspace = getattr(self.kit_permissions, "workspace", None)
                register_office_workspace(workspace or self.repo_dir)
            return domain
        except Exception:
            return ""

    def _build_machine_grant_block(self) -> str:
        """What this machine actually granted the session, measured.

        Machine facts used to be written down once, on one host, and then
        recalled everywhere -- wrong the moment the session moved. The cap
        on background jobs already derives from the real grant, so the
        MECHANISM was honest; the model just never saw the number, so
        anything it sizes itself (``pytest -n``, a PAL value, an sbatch
        header) was a guess against the machine's total.

        First turn only. A grant does not change while a session runs, and
        the per-turn steering budget belongs to blocks that do.
        """
        try:
            if len(self.messages) > 1:
                return ""
            from .bash_jobs import _affinity_cpus, _available_cpus
            cpus = int(_available_cpus())
            if cpus < 1:
                return ""
            allocated = cpus < int(_affinity_cpus())
        except Exception:
            return ""
        where = " (SLURM allocation)" if allocated else ""
        return (f"This session may use {cpus} CPU(s){where}. Size -j / -n / "
                f"PAL to that, not to the machine's total.")

    def _steering_blocks(self, role: str) -> list[tuple[str, str]]:
        """The per-turn steering blocks for an interactive role, in prompt
        order, as ``(key, text)`` with empty blocks dropped.

        One builder for two consumers — the turn-start system prompt and the
        mid-loop refresh. Kept as a single list on purpose: two copies of
        this order would drift, and the second consumer would then silently
        re-send blocks the first had already delivered.
        """
        pairs: list[tuple[str, str]] = []
        # Context-usage snapshot and the project pin are solo-only (other
        # roles run in pipeline mode with their own budgets, and office is
        # pinned to its folder by construction); the open-tasks reminder
        # applies to every interactive role that drives the task list.
        if role == "solo_agent":
            pairs.append(("project_dir", self._build_project_dir_block()))
            pairs.append(("context_status", self._build_context_status_block()))
        pairs.append(("machine_grant", self._build_machine_grant_block()))
        pairs.append(("open_tasks", self._build_open_tasks_block()))
        pairs.append(("foreign_tasks", self._build_open_foreign_tasks_block()))
        pairs.append(("unmet_delegation", self._build_unmet_delegation_block()))
        pairs.append(("unmet_tasklist", self._build_unmet_tasklist_block()))
        pairs.append(("finished_jobs", self._build_finished_jobs_block()))
        pairs.append(("finished_subagents",
                      self._build_finished_subagents_block()))
        pairs.append(("budget", self._build_budget_block()))
        pairs.append(("answered", self._build_answered_attention_block()))
        # One-time backend capability notice: reduced surfaces (no tool
        # loop, no verify enforcement, ...) are stated at session start
        # instead of failing silently mid-task. Best-effort — a missing
        # attribute on a partially built engine must not break the
        # prompt build.
        try:
            from .backend_parity import degradation_notice as _deg_notice
            pairs.append(("backend_parity", _deg_notice(
                getattr(self, "backend", "") or "",
                getattr(self, "provider", "") or "",
                first_turn=len(self.messages) <= 1,
                has_permissions=self.kit_permissions is not None) or ""))
        except Exception:
            pass
        return [(key, text) for key, text in pairs if text]

    def _drained_steering_blocks(self) -> list[tuple[str, str]]:
        """Only the blocks whose store guarantees exactly-once delivery.

        Built separately from :meth:`_steering_blocks` because they are the
        ones no role is excluded from: what a drain returns has already
        happened, is gone from the store once read, and costs nothing when
        nothing happened. Building them via the full list would drag the
        role-gated blocks along with them.
        """
        builders = {
            "finished_jobs": self._build_finished_jobs_block,
            "finished_subagents": self._build_finished_subagents_block,
            "answered": self._build_answered_attention_block,
        }
        pairs: list[tuple[str, str]] = []
        for key in _DRAINED_STEERING_KEYS:
            build = builders.get(key)
            if build is None:
                continue
            try:
                text = build()
            except Exception:
                continue
            if text:
                pairs.append((key, text))
        return pairs

    def _drain_turn_steering(self) -> list[str]:
        """Steering blocks that became true *during* the running turn.

        The system prompt is built once and then frozen for the whole tool
        loop. Every block above therefore describes the world as it was
        before the first tool call: a task the agent creates in round 2 is
        absent from the reminder that exists to make it check the task off,
        and a run budget crossing its wind-down threshold in round 30 is
        not read until the user sends another message. Measured on this
        project's benchmark workspaces, every one of these blocks is empty
        at turn start — so for a single-turn task they never fire at all.

        Returns the blocks whose text has CHANGED since the last delivery,
        capped per turn. Drain-backed blocks (finished jobs, finished
        delegates, late answers) are consumed exactly once by their own
        store, so they are passed on as they come — for EVERY role, and
        outside the refresh cap: an event that a cap swallowed would not be
        delayed, it would be destroyed, since the drain that produced it
        has already forgotten it. Read-only and best-effort: this runs
        inside the tool loop and must never end a turn that is otherwise
        working.
        """
        role = self.current_role or (self.route[0] if self.route else "")
        # Off-switch, for measuring this mechanism against itself: a
        # benchmark arm can run with the refresh disabled and be compared
        # against one that has it, which is the only way to put a number on
        # what it is worth.
        import os as _os
        if _os.environ.get(
                "DELFIN_TURN_STEERING", "").strip() in ("0", "off", "false"):
            return []
        out: list[str] = []
        try:
            out.extend(text for _key, text in self._drained_steering_blocks())
        except Exception:
            pass
        if role in _STEERING_ROLES:
            out.extend(self._changed_steering_blocks(role))
        return [f"[live update — this changed since the turn started]\n{t}"
                for t in out]

    def _changed_steering_blocks(self, role: str) -> list[str]:
        """The role-gated blocks whose text moved since the last delivery.

        Capped per turn: these are recomputed from state, so a list that
        flaps could otherwise spend a whole turn's context restating
        itself. Nothing here is consumed by being read, so a block the cap
        holds back is delivered by the next turn's system prompt.
        """
        if self._steering_refreshes >= _MAX_STEERING_REFRESHES:
            return []
        try:
            pairs = self._steering_blocks(role)
        except Exception:
            return []
        out: list[str] = []
        for key, text in pairs:
            if key not in _MID_TURN_STEERING_KEYS:
                continue
            if key in _DRAINED_STEERING_KEYS:
                # Already delivered above, for every role.
                continue
            # A window-usage line that ticks up a percent every round would
            # be pure noise; it is worth a turn's attention only once it
            # warns about the compaction that is about to fire.
            if key == "context_status" and "WARNING" not in text:
                continue
            if self._steering_delivered.get(key) == text:
                continue
            self._steering_delivered[key] = text
            out.append(text)
        if out:
            self._steering_refreshes += 1
        return out

    #: How many recent user messages the module triggers may look at, and
    #: how much of them. The matcher is substring-only, so this costs a
    #: lowercase join and nothing else; the bound is there because the
    #: transcript is unbounded and a pasted logfile would otherwise decide
    #: which modules load.
    _TRIGGER_HISTORY_MSGS = 6
    _TRIGGER_HISTORY_CHARS = 4000

    def _recent_user_text(self) -> str:
        """The last few things the USER said, for the module triggers.

        ``build_system_prompt`` has taken ``conversation_text`` through
        four layers since the triggers were fixed, and no caller ever
        passed it — so the mechanism existed end to end and one hop in
        the middle dropped it. The triggers saw a single line, which is
        the case the sticky union was added to paper over.

        User turns only. Matching the ASSISTANT's own words would let the
        agent's phrasing decide which modules load next turn, and the
        modules shape that phrasing — a loop with nothing outside it.
        """
        out: list[str] = []
        budget = self._TRIGGER_HISTORY_CHARS
        for msg in reversed(self.messages[-(self._TRIGGER_HISTORY_MSGS * 2):]):
            try:
                if (msg or {}).get("role") != "user":
                    continue
                content = msg.get("content")
                if isinstance(content, list):      # image + text blocks
                    content = " ".join(
                        str(b.get("text", "")) for b in content
                        if isinstance(b, dict) and b.get("type") == "text")
                text = str(content or "")[:budget]
            except Exception:
                continue
            if not text:
                continue
            out.append(text)
            budget -= len(text)
            if budget <= 0:
                break
        return "\n".join(reversed(out))

    def _build_current_system_prompt(
        self,
        memory_context: str = "",
        task_text: str = "",
    ) -> str:
        """Build the system prompt for the current role."""
        role = self.current_role
        if not role:
            role = self.route[0] if self.route else "builder_agent"

        self._scope_memory_domain(role)

        # Inject session error context into memory (Feature 4)
        error_ctx = self.format_error_context()
        if error_ctx:
            memory_context = f"{memory_context}\n\n{error_ctx}" if memory_context else error_ctx

        # Solo only: prepend a live context-status snapshot so the agent
        # self-monitors window-usage and can delegate to subagents before
        # compaction fires (no value for other roles — they run in
        # pipeline mode with their own context budgets).
        live_state = self._live_state
        # The session's language, decided by the first message and stated
        # BEFORE the model writes its first word. Every other attempt at
        # this corrected afterwards: the end-of-turn guard replaces a
        # finished German answer, the mid-turn note reaches the model
        # after its opening sentence is already on the user's screen.
        # Neither can stop the draft, and a user filming the session sees
        # the draft. An instruction in the prompt is the only version of
        # this rule that acts before there is anything to correct.
        _lang_block = self._session_language_block()
        if _lang_block:
            live_state = (f"{_lang_block}\n\n{live_state}" if live_state
                          else _lang_block)
        # office_agent joins the interactive roles. It was left out when this
        # gate was written and never revisited, so the mode that works on
        # real records ran on prompt text alone -- on the model this project
        # has established prompt text does not bind. Measured across the
        # archived office runs: not one had a live-state block at all, while
        # one of them spent 1.6M input tokens.
        #
        # The blocks it gains are the role-agnostic ones: open tasks,
        # finished background jobs, the run budget, a late answer to a
        # question that timed out. The context snapshot and the project pin
        # stay solo-only -- office is pinned to its folder by construction,
        # so a pin block would state what the lock already guarantees.
        #
        # A scripted role gets the drain-backed blocks only. Those are the
        # ones the gate never had a reason to hold back: a finished cluster
        # job reached three roles of eleven, so a builder or test role that
        # had submitted one was told nothing and had to poll for it.
        pairs = (self._steering_blocks(role) if role in _STEERING_ROLES
                 else self._drained_steering_blocks())
        # Remember what each block said here, so the mid-loop refresh
        # (see _drain_turn_steering) re-sends a block only once it has
        # actually changed.
        self._steering_delivered = dict(pairs)
        self._steering_refreshes = 0
        extra_blocks: list[str] = [text for _, text in pairs]
        if extra_blocks:
            joined = "\n\n".join(extra_blocks)
            live_state = f"{joined}\n\n{live_state}" if live_state else joined

        try:
            _perm_mode = getattr(self.kit_permissions, "mode", "") or ""
        except Exception:
            _perm_mode = ""
        return self.loader.build_system_prompt(
            role_id=role,
            mode_id=self.mode,
            mode_description=self.mode_description,
            route=self.route,
            role_index=self.current_role_index,
            prior_outputs=self.role_outputs if self.role_outputs else None,
            memory_context=memory_context,
            task_text=task_text,
            conversation_text=self._recent_user_text(),
            session_key=f"engine-session-{self._prompt_session_serial}",
            live_state=live_state,
            model=getattr(self, "model", "") or "",
            # Plan is a permission profile now (not a mode): inject the plan
            # addendum whenever the active permission profile is "plan".
            permission_mode=_perm_mode,
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
        images: list[str] | None = None,
        on_notice: Callable[[str], None] | None = None,
        on_tool_result_meta: Callable[[str, dict], None] | None = None,
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
        on_notice : callable, optional
            Called with harness speech — retry banners, the stop and
            empty-turn notices, the cost ceiling, a blocking hook. When
            omitted these go to ``on_token``, which is what every caller
            saw before this parameter existed.
        on_tool_result_meta : callable, optional
            Called with (tool_name, meta) after each tool result, where
            *meta* carries ``chars``, ``truncated``, ``notes``, ``ok`` and
            ``error``. ``on_tool_result`` only ever gets a 2000-character
            head slice, so without this a blocked call and a successful
            one are indistinguishable to a UI.

        Returns
        -------
        str
            The complete assistant response text.
        """
        def _notice(text: str) -> None:
            """The harness talking about itself, kept apart from the answer.

            Answer deltas and harness notices shared ``on_token`` with no
            discriminator, and every consumer that treated a notice as
            answer text drew a wrong conclusion: a run whose entire
            recorded output was three retry banners was scored as a model
            that answered badly. A UI that prints ``on_token`` verbatim
            has the same problem one layer up — it prints machinery speech
            in the middle of a sentence.

            Falls back to ``on_token`` so callers that pass no
            ``on_notice`` behave exactly as they did before.
            """
            sink = on_notice or on_token
            if sink and text:
                try:
                    sink(text)
                except Exception:
                    pass

        # New user turn: re-arm the one-correction budget — unless this IS
        # a correction turn, which must not re-arm itself. Two shapes of
        # that: the engine's own nested call (guard active), and a
        # verification retry a UI layer drives as a fresh top-level call.
        # Both carry the "[Verify]" feedback prefix; re-arming on the
        # second would let one answer be corrected twice over, and would
        # drop the very evidence ledgers the retry is judged against.
        _is_guard_feedback = str(
            user_message or "").lstrip().startswith("[Verify]")
        if not self._claim_guard_active and not _is_guard_feedback:
            self._claim_guard_corrected = False
            self._claim_guard_spent = False
            # The ambiguity ledger is per TURN: the question is whether
            # THIS answer totals a column THIS turn was told it could not
            # read. Carrying it across turns would caveat a later, unrelated
            # figure. The correction turn is a nested call and keeps the
            # ledger, or the retry would lose the very evidence it is being
            # judged against.
            self._ambiguous_columns_turn = []
            # The figure guard needs it: a number the USER typed is
            # grounded, and flagging it would be the fastest way to teach
            # somebody to ignore the guard.
            self._last_user_message = user_message or ""
            # Per turn, like every other one-shot above it: the model is
            # told once that it is answering in the wrong language, and
            # the next turn gets to be told again if it happens again.
            self._language_note_sent = False
            # Per SESSION, and only the first time: what language this
            # conversation runs in. See _note_session_language.
            self._note_session_language(user_message or "")
            # The notes the USER reads follow the session too. They were
            # hardcoded German, so an English session got an English
            # answer with German warnings stapled underneath — the rule
            # and its own mechanism disagreeing in the one place the
            # disagreement is visible.
            try:
                from . import verify_guard as _vg_lang
                _vg_lang.set_caveat_language(
                    getattr(self, "_session_language", "") or "de")
            except Exception:
                pass
            # Same rule, same reason, and it was documented as per-turn
            # while never being cleared: once ANY tool truncated in a
            # session, every later answer carrying a two-digit count got
            # the caveat, however unrelated.
            self._truncated_tools_turn = []
            # The office figure ledger is per-turn for the same reason: a
            # total the tools produced two turns ago must not ground a
            # figure stated now. The turn gets its OWN token — the ledger
            # was one process-wide list, so a second engine's turn
            # boundary cleared this one's evidence mid-flight and a
            # sub-agent's document reads filled its cap. Both flagged a
            # correct, tool-computed total.
            try:
                from . import office as _office_ledger
                self._figure_ledger_token = _office_ledger.begin_figure_turn(
                    self.session_id)
            except Exception:
                self._figure_ledger_token = ""
            # The same question on the scientific side, and per-turn for
            # the same reason: is the energy in this answer one the tools
            # returned THIS turn. Kept apart from the ledger above because
            # the evidence has a different shape -- an ORCA output is text
            # a tool returned, not a result dict it filled in.
            try:
                from . import verify_guard as _vg_numbers
                _vg_numbers.reset_observed_numbers()
                _vg_numbers.reset_keyed_values()
            except Exception:
                pass

        # UserPromptSubmit hooks — fire BEFORE the message is appended so a
        # blocking hook can short-circuit the turn entirely. Stop hooks are
        # fired in the finally-block at the end of the turn.
        try:
            from . import hooks as _hooks_mod
            from .api_client import _session_hooks
            _kperms = getattr(self, "kit_permissions", None)
            _ws = _kperms.workspace if _kperms else None
            # Through the same funnel every other hook point uses. Loading
            # from the workspace alone never sees the permissions object,
            # so the two hook kinds fired from HERE — UserPromptSubmit and
            # Stop — went on running under a flag that says no ambient
            # configuration is read. A bounding flag that covers four of
            # six hook points is the silent non-delivery it exists to
            # prevent, one layer in.
            _hooks_cfg = (_session_hooks(_kperms) if _kperms
                          else _hooks_mod.HooksConfig())
            if not _hooks_cfg.is_empty():
                _ups = _hooks_mod.run_hooks(
                    "UserPromptSubmit", _hooks_cfg,
                    user_prompt=user_message, workspace=_ws,
                )
                _blk = _hooks_mod.first_block(_ups)
                if _blk is not None:
                    _reason = _blk.reason or _blk.stderr or "blocked by UserPromptSubmit hook"
                    _notice(f"\n[hook block] {_reason}\n")
                    return f"[hook block] {_reason}"
        except Exception:
            _hooks_cfg = None
            _ws = None

        # Run-budget hard gate: past 110% of the configured budget no new
        # turn starts (the 80-100% band already told the agent to wind
        # down). The session is intact — raise the budget or start a new
        # run to continue.
        try:
            _bfrac, _, _bdim = self._run_budget_detail()
            if _bfrac >= 1.1:
                _msg = (
                    "🛑 Run budget exhausted (>110%). No further turns will "
                    "run in this session. State is saved — raise "
                    f"{_bdim or 'agent.run_budget_usd'} or start a new run "
                    "to continue.")
                self._emit_budget_attention(
                    "hard_stop",
                    "Run stopped: budget ceiling reached",
                    _msg)
                _notice(_msg)
                return _msg
        except Exception:
            pass

        # One conversation, one turn at a time. Nothing used to refuse a
        # second turn while the first was still running: a Stop clears the
        # widget flag without waiting for the worker, so the next send
        # started a second thread on THIS engine. Both appended to this
        # message list, both mutated the per-turn cost state, and both ran
        # tools in the same workspace. The guard belongs here, on the object
        # that owns the conversation, not on a flag a Stop has already
        # cleared.
        with self._turn_gate:
            if self._turn_in_flight:
                _busy = (
                    "A turn is already running on this session. Stop it "
                    "first, or send this once it has finished — starting a "
                    "second turn would interleave both into one history.")
                _notice(_busy)
                return _busy
            self._turn_in_flight = True
            self._turn_serial += 1
            self._turn_id = self._turn_serial
            # The stop reset lives HERE, behind the gate, and nowhere
            # earlier. It used to be the first statement of the method,
            # which meant a turn the gate was about to refuse still wiped
            # the flag on its way out: Stop during a long tool call, send
            # the next message, and that message's turn -- refused, having
            # done nothing -- deleted the stop the running turn was about
            # to read. The running turn then carried on at its next round
            # boundary and billed for a plan the user had abandoned.
            # Whoever owns the turn owns its stop; a caller that never got
            # a turn owns nothing.
            self._stop_requested = False
            self._stop_owner_turn = None
        _turn_id = self._turn_id

        # This turn's delegated spend starts at zero, and every delegated
        # run from here is routed through a counter that knows which turn
        # asked for it. Installed per turn rather than once, because the
        # client re-binds the runner whenever permissions change; the
        # installer unwraps an existing meter first, so repeating it
        # cannot stack two counters over one runner.
        try:
            self.delegate_spend().turn = _agent_metrics.DelegateSpend()
            self._meter_delegate_spend()
        except Exception:
            pass

        # Read once per turn, before the prompt is built: the steering
        # blocks are assembled from engine state and never see the message
        # that triggered them, so anything that has to know what was asked
        # has to be captured here.
        self._turn_is_bare_greeting = self.is_bare_greeting(user_message)

        self.messages.append(self._build_user_message(user_message, images))
        # Sanitize message history: ensure proper user/assistant alternation.
        # Concurrent stop/send can leave consecutive user messages.
        self._sanitize_messages()

        # Mid-conversation compaction. Fires for solo/dashboard on the
        # legacy message-count threshold, and for any role when the
        # token-budget threshold is exceeded.
        self._compact_history()

        system_prompt = self._build_current_system_prompt(
            memory_context,
            task_text=user_message,
        )

        # Context distillation for expensive calls (Feature 5)
        if self._distiller and self._distiller.should_distill(
            system_prompt, self.messages,
        ):
            system_prompt = self._distiller.distill(
                system_prompt, task_text=user_message,
            )

        # Retain the exact system prompt that drove this turn so a bug
        # report can show *why* the agent behaved as it did (the injected
        # playbook + memory determine behaviour). Kept on the engine, not
        # in the message list, since the API takes it as a separate arg.
        self.last_system_prompt = system_prompt
        self._system_prompt_chars = len(system_prompt or "")

        # Let the tool loop re-read the steering blocks between rounds. The
        # system prompt above is frozen from here until the turn ends, so
        # without this callback anything that becomes true during the turn
        # waits for the next user message. Best-effort: a client that does
        # not look for the attribute is unaffected by carrying it.
        try:
            self.client.steering_provider = self._drain_turn_steering
        except Exception:
            pass

        # Let the tool loop see a Stop. The engine only checks its own flag
        # between STREAM EVENTS, and the loop that runs tools between those
        # events never looked at it -- the name did not occur in the client
        # module at all. During a long command or a stalled stream, Stop set
        # a flag nobody read.
        try:
            self.client.should_stop = lambda: bool(self._stop_requested)
        except Exception:
            pass

        # ...and the turn's cost ceiling, for the same reason: the check
        # below fires on message_delta, which every client emits once, at
        # the end. Checked in the loop it bounds the spend instead of
        # reporting it.
        try:
            self.client.turn_cost_cap = lambda: float(
                getattr(self, "_cost_cap_value", 0.0) or 0.0)
        except Exception:
            pass

        # Resolve max_tokens: caller override > role default > global default
        effective_max = max_tokens or self.max_tokens_for_role(self.current_role)

        # Sync the current role onto the KitToolPermissions so the
        # api_client can gate the advertised tool surface per-role
        # (dashboard_agent gets read-only; coding roles get the full set).
        try:
            _kp = getattr(self, "kit_permissions", None)
            if _kp is not None:
                _kp.agent_role = self.current_role or ""
        except Exception:
            pass

        chunks: list[str] = []
        # Events still dispatched after a stop has been seen (see the loop
        # below): enough for a client's closing notice, far too few for a
        # client that ignores the stop to keep generating.
        _stop_drain = 3
        # Reasoning characters this turn. A model that streams its whole
        # answer on the reasoning channel produces zero chunks and zero
        # exceptions, and the empty-turn branch below needs to be able to
        # say so instead of reporting "no answer" with no explanation.
        _thinking_chars = 0
        # Per-turn timing: capture time-to-first-token + tool count so a slow
        # turn can be diagnosed after the fact (backend stall vs generation vs
        # tool rounds). Recorded in the finally below so errored/partial turns
        # are captured too. See turn_metrics.
        _turn_t0 = _time.monotonic()
        # Cleared per turn: a diagnostic left over from an earlier empty
        # turn must not be read as a report about this one.
        self.last_empty_turn = None
        _usage_before = dict(self.token_usage)
        _turn_ttft: float | None = None
        # Time of the FIRST stream event of any kind — message_start and
        # thinking deltas included, neither of which stamps ttft. Without
        # it a live connection to a silent model and a transport that
        # delivered nothing at all wrote the same record, and the log is
        # the only evidence left afterwards. See turn_metrics.silence_kind.
        _turn_first_event: float | None = None
        # The exception this turn died on, kept for the record below. The
        # except clause's own name is unbound by the time `finally` runs,
        # and `record()` has had an `error` parameter from the start that
        # nothing ever filled -- so a crash and a silent backend logged
        # identically.
        _turn_error = ""
        _turn_tool_calls = 0
        # Tool NAMES this turn — evidence input for the claim-grounding
        # guard (a lookup tool used this turn grounds quantity claims).
        _turn_tool_names: list[str] = []
        # Mid-turn crash-checkpoint throttle: first write fires 10 tool
        # results or 60s into the tool loop, whichever comes first — short
        # turns never touch disk. See session_store.save_turn_checkpoint.
        _ckpt_events = 0
        _ckpt_last = _turn_t0
        self._saw_message_start = False
        self._floor_captured_this_turn = False
        # Per-turn runaway circuit-breaker: snapshot the cost at turn start so a
        # single turn's tool-loop can't run away forever. Resets every turn — it
        # is NOT a cumulative session budget. (_MAX_TOOL_ROUNDS already bounds
        # rounds; this is the extra, very-high cost ceiling.)
        self._turn_start_cost = self.cost_usd
        self._cost_cap_hit = False
        self._cost_cap_value = self._cost_hard_cap()   # read once per turn
        try:
            # Pass no_tools only to a client that accepts it. Adding a
            # parameter must not break a caller built against the older
            # signature -- test doubles and any out-of-tree backend. I made
            # exactly this mistake once already today with _fetch_bytes;
            # here the cost would have been every mocked engine test.
            _stream_kwargs = {}
            try:
                import inspect as _inspect
                if "no_tools" in _inspect.signature(
                        self.client.stream_message).parameters:
                    _stream_kwargs["no_tools"] = self.is_bare_greeting(
                        user_message)
            except (TypeError, ValueError):
                pass
            for event in self.client.stream_message(
                system=system_prompt,
                # Wire view: identical to self.messages except private
                # bookkeeping keys (e.g. _pinned) are stripped — strict
                # backends reject unknown message fields.
                messages=self._wire_messages(),
                max_tokens=effective_max,
                session_id=self.session_id,
                thinking_budget=thinking_budget,
                **_stream_kwargs,
            ):
                # ANY event, before any filtering: this stamp answers
                # "did the transport deliver anything at all?", which is
                # a different question from "did a token arrive?".
                if _turn_first_event is None:
                    _turn_first_event = _time.monotonic()

                if self._stop_requested and event.type != "message_delta":
                    # A stop ends the turn here -- but not before the
                    # events the client emits BECAUSE of the stop. This
                    # used to be an unconditional break, and what a client
                    # sends the moment it sees the flag is its closing
                    # notice followed by the message_delta carrying this
                    # turn's tokens and cost. Breaking on the notice threw
                    # away both: a stopped turn reported $0.00 however many
                    # rounds it had already billed, the notice never
                    # reached the user, and the run-budget gate was blind
                    # to the whole turn.
                    #
                    # The notice is forwarded to the caller but NOT added
                    # to ``chunks``: it is machinery talking, not the
                    # model's answer, and the stop branch below still needs
                    # to see an empty response to take the unanswered
                    # message back out of the history.
                    #
                    # The budget is the reason this cannot become a way to
                    # keep streaming through a stop: a client that ignores
                    # the flag gets a handful of events, not a free run.
                    # A notice counts here as well: the sentence that says
                    # WHY the turn ended is emitted as one, and dropping
                    # it at the stop boundary would end the turn without a
                    # word -- which is the failure this whole drain exists
                    # to prevent.
                    if (event.type in ("text_delta", "notice") and event.text
                            and _stop_drain > 0):
                        _stop_drain -= 1
                        if event.type == "notice":
                            _notice(event.text)
                        elif on_token:
                            try:
                                on_token(event.text)
                            except Exception:
                                pass
                        continue
                    break

                if event.type == "text_delta" and event.text:
                    if _turn_ttft is None:
                        _turn_ttft = _time.monotonic()
                    chunks.append(event.text)
                    if on_token:
                        on_token(event.text)

                elif event.type == "notice" and event.text:
                    # The harness talking about itself: a retry banner, a
                    # stop, a cost ceiling. Shown to the user and to
                    # nothing else -- it is not the model's answer, and
                    # every place that treated it as one drew a wrong
                    # conclusion from it. A run whose entire recorded
                    # output was three retry banners was scored as a
                    # model that answered badly, at a quality number, and
                    # written into the file baselines compare against;
                    # the same text stamped time-to-first-token, so the
                    # sickest turn produced the healthiest-looking
                    # record. Not appended to `chunks`, and _turn_ttft is
                    # deliberately left alone.
                    _notice(event.text)

                elif event.type == "thinking_delta" and event.text:
                    # A reasoning token IS the backend starting to produce.
                    # Left unstamped, a turn that thought for two minutes
                    # and then answered recorded ttft_ms=None, which is the
                    # same record a backend that never said anything
                    # writes -- and on an endpoint that reports no usage
                    # the output-token corroboration is absent too, so a
                    # thinking-only turn could be counted as never
                    # started. Widens what ttft measures from "first
                    # visible word" to "first token of any kind", which is
                    # the question the health report actually asks.
                    if _turn_ttft is None:
                        _turn_ttft = _time.monotonic()
                    _thinking_chars += len(event.text)
                    if on_thinking:
                        on_thinking(event.text)

                elif event.type == "tool_use":
                    # The first thing the user reads is the sentence before
                    # this call, and the end-of-turn guard cannot reach it.
                    self._nudge_language_if_wrong(chunks)
                    # Code-level tool whitelist enforcement
                    role_id = self.route[self.current_role_index] if self.route else ""
                    allowed = _ROLE_TOOL_WHITELIST.get(role_id)
                    if allowed is not None and event.tool_name not in allowed:
                        # Allow MCP doc + ops server tools for all roles, plus
                        # the KIT-Toolbox coding tools (already permission-gated).
                        if not (
                            event.tool_name.startswith(_MCP_TOOL_PREFIXES)
                            or event.tool_name.startswith(_KIT_CODING_PREFIX)
                        ):
                            # Block unauthorized tool — don't call on_tool_use
                            continue
                    # Dashboard scope: the mcp__ exemption above lets KIT-backend
                    # tools bypass the role whitelist, but the dashboard agent
                    # must never inspect DELFIN source. Re-block the forbidden
                    # source-reading tools by their bare (de-namespaced) name.
                    if role_id == "dashboard_agent":
                        _bare = event.tool_name
                        for _pfx in (*_MCP_TOOL_PREFIXES, _KIT_CODING_PREFIX):
                            if _bare.startswith(_pfx):
                                _bare = _bare[len(_pfx):]
                                break
                        if _bare in _DASHBOARD_SOURCE_DENY:
                            continue  # dashboard scope: no source inspection
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
                    # Trace: stash this call until its result arrives.
                    if _turn_ttft is None:
                        _turn_ttft = _time.monotonic()
                    _turn_tool_calls += 1
                    _turn_tool_names.append(event.tool_name)
                    # Evidence for the functional-claim guard: record WHAT
                    # was executed, not just that some tool ran. Held here
                    # and committed when the RESULT arrives -- recording it
                    # at the call meant a command that was denied by a
                    # gate, blocked by a hook, or died with a traceback
                    # counted as "run", and the guard that exists to catch
                    # "it works now" then found its evidence. The ledger is
                    # about what happened, not what was attempted.
                    self._exec_pending.append(
                        (event.tool_name, event.tool_input))
                    if len(self._exec_pending) > 64:
                        del self._exec_pending[:-64]
                    self._trace_pending.append(
                        (event.tool_name, event.tool_input, _time.monotonic()))
                    self._maybe_pin_project_dir(
                        event.tool_name, event.tool_input)
                    if on_tool_use:
                        on_tool_use(event.tool_name, event.tool_input)

                elif event.type == "tool_result":
                    if on_tool_result and event.tool_output:
                        on_tool_result(event.tool_name, event.tool_output)
                    # Commit the held command only if it actually ran.
                    self._commit_exec_command(
                        event.tool_name, event.tool_output or "")
                    # Evidence for the ambiguous-column guard: which columns
                    # a reader said it could not decide. Taken from the
                    # reader's own note rather than re-derived, so the two
                    # cannot disagree about what is in question.
                    #
                    # tool_output is a HEAD slice, and read_document writes
                    # its notes after the grid — on a 200-row sheet the note
                    # is ~10 kB in and never arrived. The notes therefore
                    # ride along as their own field.
                    _out = event.tool_output or ""
                    _notes = str(getattr(event, "output_notes", "") or "")
                    self._note_ambiguous_columns(
                        _out + ("\n" + _notes if _notes else ""))
                    # Every number this result carried, so a quantity in
                    # the answer can be checked against what came back
                    # rather than against the fact that something ran.
                    try:
                        from . import verify_guard as _vg_numbers
                        _vg_numbers.record_tool_numbers(
                            _out,
                            truncated=bool(
                                getattr(event, "output_truncated", False)))
                        # The same numbers, keyed by record and field, so
                        # two of them can be compared rather than merely
                        # counted as present. The path comes from the call
                        # that produced this result — a value read out of
                        # <folder>/DELFIN_Data.json and one read back out
                        # of a table row keyed on <folder> describe the
                        # same thing, and that is the only way to see it.
                        _src = ""
                        for _pn, _pi, _ in reversed(self._trace_pending):
                            if _pn == event.tool_name:
                                _src = str(_pi)
                                break
                        _vg_numbers.record_keyed_values(_out, source=_src)
                    except Exception:
                        pass
                    # And whether the model saw the whole result. A count
                    # taken from output that was cut short is an estimate
                    # wearing the clothes of a measurement. The backend
                    # reports the cut as a field, because the marker itself
                    # is written past the slice this event carries. Backends
                    # that pass the result through whole (the CLI path) set
                    # no field and still have the marker in the text.
                    if getattr(event, "output_truncated", False) or (
                            "truncated," in _out and "chars" in _out):
                        self._note_truncated_tool(event.tool_name or "a tool")
                    # Whether it worked is decided by the result, not
                    # asserted. ok=True was hardcoded here, and the only
                    # writer of ok=False is the permission_denied event --
                    # which ONLY the CLI backend emits. On every other
                    # backend a gate refusal arrives as an ordinary
                    # tool_result carrying {"error": ...} and was traced
                    # green. Five consumers read this flag: /trace printed
                    # a checkmark for the blocked call, the aggregates
                    # reported an error rate of zero forever, `/agents
                    # tools` showed a permanent 0% error column, and the
                    # live panel rendered a green tick. A user looking for
                    # what went wrong found a clean list.
                    _failed, _reason = self._tool_result_failed(_out)
                    self._record_tool_trace(
                        event.tool_name, _out,
                        ok=not _failed, error=_reason)
                    # The same verdict, offered to the caller. on_tool_result
                    # hands over a 2000-character head slice and nothing
                    # else, so a refusal and a success look identical to a
                    # UI -- the very confusion the trace flag above was
                    # fixed for, one layer further out.
                    if on_tool_result_meta:
                        try:
                            on_tool_result_meta(event.tool_name or "", {
                                "chars": int(getattr(
                                    event, "output_chars", 0) or len(_out)),
                                "truncated": bool(getattr(
                                    event, "output_truncated", False)),
                                "notes": _notes,
                                "ok": not _failed,
                                "error": _reason or "",
                            })
                        except Exception:
                            pass
                    # Crash insurance: full session saves happen only at
                    # turn boundaries, but one turn can run hundreds of
                    # tool rounds — persist a cheap checkpoint (throttled,
                    # best-effort) so a SIGKILL mid-loop costs the last few
                    # rounds, not the whole turn. Cleared at turn end.
                    _ckpt_events += 1
                    _ckpt_now = _time.monotonic()
                    if _ckpt_events >= 10 or (_ckpt_now - _ckpt_last) >= 60.0:
                        _ckpt_events = 0
                        _ckpt_last = _ckpt_now
                        try:
                            from . import session_store as _ss_ckpt
                            _ss_ckpt.save_turn_checkpoint(
                                self.session_id or "", {
                                    "user_message": user_message[:2000],
                                    "partial_response": "".join(chunks)[-20000:],
                                    "tool_calls": _turn_tool_calls,
                                    "ts": _time.time(),
                                })
                        except Exception:
                            pass

                elif event.type == "permission_denied":
                    if on_permission_denied:
                        on_permission_denied(event.tool_name)
                    self._record_tool_trace(
                        event.tool_name, "", ok=False, error="permission denied")

                elif event.type == "session_init":
                    # Capture session ID from CLI for persistence
                    if event.text:
                        self.session_id = event.text

                elif event.type == "message_start":
                    # Input tokens tracked here (authoritative count
                    # including cache).  Do NOT also add in message_delta.
                    with self._lock:
                        self.token_usage["input"] += event.input_tokens
                        self._saw_message_start = True
                        # Snapshot the real per-request input count so the
                        # compaction budget can use it as a ground-truth
                        # floor (it includes the system prompt + tool
                        # schemas that self.messages omits) -- but only
                        # from the FIRST request of the turn.
                        #
                        # The tool-loop backends emit one message_start per
                        # ROUND, and each later round carries that turn's
                        # accumulated tool results. Those results never
                        # enter self.messages and are gone by the time the
                        # next request is built, so keeping round N made the
                        # floor the intra-turn peak. Measured on a real
                        # loop: true next-request size 1,155 tokens, floor
                        # kept 17,634 -- 15x. The next turn's compaction
                        # then saw pressure that did not exist and trimmed
                        # every older message on a conversation occupying
                        # 2% of the window, every turn, for the rest of the
                        # session; above the summary threshold it also
                        # bought an LLM summarisation call to destroy nine
                        # messages worth 1.5% of the window. The same
                        # estimate feeds the self-monitoring block, which
                        # told the agent it was at 89% and should wind down
                        # while it held one token of history.
                        #
                        # The docstring already said "the last real count
                        # for THIS context". Round N's prompt is this
                        # context plus ephemeral tool results, so the code
                        # contradicted its own contract.
                        if event.input_tokens and not self._floor_captured_this_turn:
                            self._floor_captured_this_turn = True
                            self._last_input_tokens = int(event.input_tokens)
                            # Fresh provider count -> trims applied since the
                            # previous floor are already reflected in it.
                            self._trimmed_chars_since_floor = 0
                        # Prompt tokens served from the endpoint cache (Anthropic
                        # reports it here). Visibility into how much of the input
                        # was free — drives the caching/efficiency work.
                        self.token_usage["cached"] = self.token_usage.get(
                            "cached", 0) + (getattr(event, "cached_tokens", 0) or 0)

                elif event.type == "message_delta":
                    with self._lock:
                        # Only output tokens and cost from the final event.
                        # Input tokens already counted in message_start.
                        self.token_usage["output"] += event.output_tokens
                        self.cost_usd += event.cost_usd
                        # Accounting fallback for backends that never emit
                        # message_start: count the turn's input here so
                        # /context and the self-monitoring block show real
                        # numbers. NOT used as the compaction floor — this
                        # value is cumulative across tool rounds and would
                        # overestimate a single request.
                        if (not self._saw_message_start
                                and getattr(event, "input_tokens", 0)):
                            self.token_usage["input"] += event.input_tokens
                        # OpenAI/vLLM report cached prompt tokens here.
                        self.token_usage["cached"] = self.token_usage.get(
                            "cached", 0) + (getattr(event, "cached_tokens", 0) or 0)
                    # Per-turn runaway breaker: if THIS turn's spend crosses the
                    # (very high) hard cap, stop so the loop can't run forever.
                    _cap = getattr(self, "_cost_cap_value", 0.0)
                    if (_cap > 0 and not self._cost_cap_hit
                            and (self.cost_usd - self._turn_start_cost) >= _cap):
                        self._cost_cap_hit = True
                        self._stop_requested = True
                        # This stop belongs to THIS turn, like any other.
                        self._stop_owner_turn = _turn_id
                        _note = (
                            f"\n\n🛑 Cost circuit-breaker: this turn hit the "
                            f"${_cap:.0f} per-turn hard cap and was stopped so "
                            f"the loop can't run away. (Per turn, NOT a session "
                            f"budget — raise via agent.cost_hard_limit_usd.)")
                        # The stream is the only place this was ever said.
                        # One event per broken turn, so a run that trips it
                        # repeatedly is visible as repetition, not as one line.
                        self._emit_budget_attention(
                            f"cost_breaker:{_turn_id}",
                            "Turn stopped by the cost circuit-breaker",
                            _note.strip())
                        chunks.append(_note)
                        _notice(_note)
                    # Capture session ID from result event
                    if event.text and not self.session_id:
                        self.session_id = event.text
                    if self._stop_requested:
                        # The accounting above is what the stop check at the
                        # top of the loop was let past for. Nothing further
                        # from this stream may be dispatched.
                        break

        except Exception as _turn_exc:
            # Take the plain form first, so the turn log gets an error even
            # if the classification below fails.
            _turn_error = f"{type(_turn_exc).__name__}: {_turn_exc}"
            if not chunks:
                self.messages.pop()
            else:
                # Text streamed and then the turn died. Re-raising here
                # skips the append below, so the history ended on a user
                # message with no answer -- and the alternation sanitiser
                # resolves two consecutive user messages by keeping the
                # NEWEST. The next thing the user typed therefore replaced
                # their original question, which the dashboard could not
                # show because it had already printed the partial text.
                # Commit what arrived, marked as unfinished so the model
                # does not read its own truncated output as an answer.
                self.messages.append({
                    "role": "assistant",
                    "content": "".join(chunks)
                    + "\n\n[this turn ended with an error before it "
                      "finished; the text above is incomplete]",
                })
            # The FAIL side of the learning loop: errored turns previously
            # recorded NO outcome at all, so provider profiles only ever
            # learned from successes. Classify + record before re-raising.
            try:
                _etype = type(_turn_exc).__name__
                try:
                    from .api_client import (_is_context_length_error,
                                             _is_transient_api_error)
                    if _is_context_length_error(_turn_exc):
                        _etype = "context_overflow"
                    elif _is_transient_api_error(_turn_exc):
                        _etype = "transient_api"
                except Exception:
                    pass
                # The same classification the profiles get. It was already
                # being computed here while the turn log was told nothing.
                _turn_error = f"{_etype}: {_turn_exc}"
                self.record_cycle_outcome(
                    "FAIL", user_message, error_type=_etype,
                    start_time=_turn_t0)
            except Exception:
                pass
            raise
        finally:
            # Release the one-turn gate on EVERY exit path. A failed turn
            # that left it set would make the engine refuse work forever.
            with self._turn_gate:
                self._turn_in_flight = False
            # Stop hooks — fire after the stream finishes (success or
            # exception). Failures are silenced so a misconfigured hook
            # never bubbles up.
            try:
                if _hooks_cfg is not None and not _hooks_cfg.is_empty():
                    from . import hooks as _hm2
                    _hm2.run_hooks("Stop", _hooks_cfg, workspace=_ws)
            except Exception:
                pass
            # Per-turn timing telemetry (best-effort; never affects the turn).
            # Captures total + time-to-first-token so a slow turn is diagnosable
            # after the fact (backend stall vs generation vs tool rounds).
            #
            # The failure fields travel with it: the exception the turn died
            # on, the first-EVENT stamp that says whether the transport ever
            # delivered, and the client's retry count. Without all three,
            # every way of producing nothing wrote the same record.
            try:
                from . import turn_metrics as _tm
                _tm.record(
                    self.trace_session(),
                    provider=getattr(self, "provider", "") or "",
                    model=getattr(self.client, "model", "") or "",
                    role=self.current_role or "",
                    total_ms=int((_time.monotonic() - _turn_t0) * 1000),
                    ttft_ms=(int((_turn_ttft - _turn_t0) * 1000)
                             if _turn_ttft is not None else None),
                    first_event_ms=(int((_turn_first_event - _turn_t0) * 1000)
                                    if _turn_first_event is not None else None),
                    output_chars=sum(len(c) for c in chunks),
                    tool_calls=_turn_tool_calls,
                    stopped=self._stop_requested,
                    error=_turn_error,
                    retries=self._client_retry_count(),
                    # Per-TURN deltas, not the running totals: a cache hit
                    # rate is only meaningful per request.
                    input_tokens=max(
                        0, int(self.token_usage.get("input", 0))
                        - int(_usage_before.get("input", 0))),
                    output_tokens=max(
                        0, int(self.token_usage.get("output", 0))
                        - int(_usage_before.get("output", 0))),
                    cached_tokens=max(
                        0, int(self.token_usage.get("cached", 0))
                        - int(_usage_before.get("cached", 0))),
                )
            except Exception:
                pass
            # What the turn delegated, said once, by the turn that paid
            # for it. The figure every other surface shows for a turn is
            # the parent model's own tokens, so a turn that handed the
            # work to five sub-agents read as one of the cheapest of the
            # session. Nothing is printed when nothing was delegated: a
            # "delegated $0.00" line after every turn is noise that trains
            # a reader to skip the line the real number appears on.
            try:
                _deleg_line = _agent_metrics.format_turn_delegation(
                    float(self.cost_usd) - float(self._turn_start_cost),
                    self.delegate_spend().turn)
                if _deleg_line:
                    _notice("\n" + _deleg_line)
            except Exception:
                pass
            # Whether this turn's spend could be measured at all. The USD
            # gates read a bare self.cost_usd, to which an unpriced turn
            # contributes 0.0 -- indistinguishable from a cheap one. See
            # _usd_budget_enforced.
            try:
                self._note_turn_price_state(
                    float(self.cost_usd) - float(self._turn_start_cost))
            except Exception:
                pass
            # Turn is over (normally or with a surfaced error) — the
            # mid-turn crash checkpoint is obsolete. Only an unclean
            # process death (which skips this finally) leaves one behind
            # for the next session load to recover from.
            try:
                from . import session_store as _ss_ckpt
                _ss_ckpt.clear_turn_checkpoint(self.session_id or "")
            except Exception:
                pass

        # Copy the per-turn verification mechanics from the backend client
        # (report_verdict tool + test-evidence ledger) and snapshot them for
        # the role that just ran, so the acceptance gate can still consult
        # the TEST role's turn after later roles have overwritten the
        # per-turn state. Best-effort: CLI backends have neither attribute.
        try:
            self._last_structured_verdict = getattr(
                self.client, "_last_structured_verdict", None)
            self._last_test_evidence = list(
                getattr(self.client, "_test_evidence", None) or [])
            # Observed-files ledger (session-cumulative) for the code-claim
            # citation check — files the model actually read/grepped. The
            # availability flag distinguishes "backend keeps no ledger"
            # (claim grounding cannot judge -> stays off) from "ledger kept
            # but empty" (zero evidence -> enforcement may fire).
            #
            # UNION, never replace. restore_state deliberately leaves the
            # client alone and the client mints its ledger fresh, so a
            # replace here threw away everything the resumed session had
            # read before the save -- on the FIRST resumed turn, before any
            # guard had read it. The grounding guard then flagged every
            # file:line restated from the restored history as unsupported,
            # which is precisely the false correction the evidence export
            # exists to prevent.
            _obs_ledger = getattr(self.client, "_observed_files_session", None)
            _restored_flag = bool(
                getattr(self, "_observed_ledger_available", False))
            self._observed_ledger_available = (
                _obs_ledger is not None or _restored_flag)
            try:
                _carried = set(getattr(self, "_last_observed_files", None) or ())
                self._last_observed_files = _carried | set(_obs_ledger or ())
            except TypeError:
                self._last_observed_files = set(
                    getattr(self, "_last_observed_files", None) or ())
                self._observed_ledger_available = _restored_flag
            self._last_turn_tools = list(_turn_tool_names)
            # Session-wide tool names: an open delegation request is
            # satisfied by a sub-agent run in ANY turn of the session.
            if not hasattr(self, "_session_tool_names"):
                self._session_tool_names: set[str] = set()
            self._session_tool_names.update(
                str(t) for t in _turn_tool_names)
            _v_role = self.current_role
            if _v_role:
                if isinstance(self._last_structured_verdict, dict):
                    self.role_verdicts[_v_role] = self._last_structured_verdict
                if self._last_test_evidence:
                    self.role_test_evidence.setdefault(_v_role, []).extend(
                        self._last_test_evidence)
        except Exception:
            pass

        full_response = "".join(chunks)
        # True once the empty-turn branch below has written a diagnostic
        # into full_response — it is a report about the turn, not an answer
        # the model gave, so no guard may treat it as one.
        _empty_turn = False
        # Repair corrupted output (harmony tool-channel leaks + glitch
        # tokens from gpt-5.x via the OpenAI-compatible endpoint) BEFORE it
        # enters the conversation context — otherwise the garbage is replayed
        # to the model on every subsequent turn.
        if full_response:
            try:
                from delfin.agent.text_sanitize import sanitize_agent_text
                full_response = sanitize_agent_text(full_response).text
            except Exception:
                pass
        # Output-guard stage: redact credential material from the FINAL
        # answer before it enters the transcript. Known limitation: the
        # individual tokens were already streamed out via on_token, so a
        # secret may have been displayed once live — this stage protects the
        # stored history and every downstream consumer (context replay,
        # session files, subagent reports, memory distillation). Guarding
        # the live stream as well would require buffering the whole
        # response before emitting anything, which is deliberately not done.
        _guard_note = ""
        if full_response:
            try:
                from delfin.agent.output_guard import run_output_guards
                _guard = run_output_guards(full_response)
                if _guard.changed:
                    full_response = _guard.text
                    _redactions = [f for f in _guard.findings
                                   if f.get("check") == "secret_redaction"]
                    if _redactions:
                        _kinds = sorted({f.get("detail", "")
                                         for f in _redactions})
                        _guard_note = (
                            f"\n[output-guard] {len(_redactions)} finding(s): "
                            f"{', '.join(_kinds)} — sensitive content "
                            f"redacted.")
            except Exception:
                pass
        if full_response:
            self.messages.append({"role": "assistant", "content": full_response})

            # Context usage tracking (Feature 4)
            self._record_context_usage(full_response)

            # Progressive disclosure: NEED_CONTEXT detection (Feature 1)
            if self._progressive_disclosure:
                full_response = self._handle_need_context(
                    full_response, user_message, memory_context,
                    on_token=on_token,
                    on_tool_use=on_tool_use,
                    on_tool_result=on_tool_result,
                    on_permission_denied=on_permission_denied,
                    on_thinking=on_thinking,
                    thinking_budget=thinking_budget,
                    max_tokens=max_tokens,
                )
        elif self._stop_requested:
            # Stop before any answer text — during thinking, or before the
            # first token. The message appended at the top of this turn has
            # nothing after it, and two consecutive user messages break the
            # API on the next turn, so it comes back out. That much was
            # already right.
            #
            # Saying so was missing, and the silence was worse than saying
            # nothing: the stop notice tells the user "the rounds completed
            # so far are kept; send a message to continue from here", which
            # is true of a stop inside the tool loop and false here. Nothing
            # was kept, and the question they asked is no longer in the
            # history to continue FROM. Measured before this: the turn
            # returned '', on_token was never called, and the history came
            # back empty -- so the only sentence the user got was the one
            # that was wrong.
            #
            # Same repair as the empty-turn branch below, which already
            # names what it took out.
            if self.messages and self.messages[-1].get("role") == "user":
                self.messages.pop()
                full_response = (
                    "[stopped] Stopped before any answer text. Nothing was "
                    "added to the history and your message is not in it "
                    "either — send it again to pick the question back up."
                    + self._hand_back_undelivered_steer()
                )
                _notice(full_response)
        else:
            # No text, no stop, no exception. This branch did not exist,
            # and its absence destroyed the question: the user message
            # appended above stayed in the history with no answer after
            # it, and on the next turn the alternation sanitiser resolves
            # two consecutive user messages by keeping the NEWEST -- so
            # the follow-up silently overwrote the original task and the
            # model answered a question nobody had asked.
            #
            # It needs no crash to happen. A reasoning model that streams
            # its whole answer on the reasoning channel produces thinking
            # deltas only, which never reach ``chunks``; and a CLI turn
            # whose process dies without a diagnostic used to end the same
            # way. The exception path already handled this hazard; the
            # quiet path never did.
            #
            # Same repair as the exception path: take the unanswered
            # message back out, and say what was observed instead of
            # returning "" for the caller to interpret as an answer.
            _empty_turn = True
            if self.messages and self.messages[-1].get("role") == "user":
                self.messages.pop()
            _empty_elapsed = _time.monotonic() - _turn_t0
            full_response = (
                f"[empty turn] The backend ended this turn without any "
                f"answer text — {_thinking_chars} characters of reasoning, "
                f"{_turn_tool_calls} tool call(s), {_empty_elapsed:.1f}s. "
                f"Nothing was added to the history, so your message is "
                f"unchanged: send it again, or switch model if this "
                f"repeats (a model that answers on the reasoning channel "
                f"only produces exactly this)."
            )
            # Named, so the turn log can tell an empty turn apart from a
            # normal finish instead of recording a successful cycle that
            # answered nothing.
            self.last_empty_turn = {
                "thinking_chars": _thinking_chars,
                "tool_calls": _turn_tool_calls,
                "duration_s": round(_empty_elapsed, 3),
                "model": str(getattr(self.client, "model", "") or ""),
                "role": self.current_role or "",
            }
            try:
                self.record_cycle_outcome(
                    "FAIL", user_message, error_type="empty_turn",
                    start_time=_turn_t0)
            except Exception:
                pass
            _notice(full_response)

        # Claim-grounding enforcement — every mode funnels through this
        # method, so the guard runs here (not in any UI layer): a final
        # answer asserting a specific code location or a measured quantity
        # without matching observed evidence gets exactly ONE forced
        # correction turn; a still-ungrounded correction gets a visible
        # caveat. The nested correction turn (guard active) skips this
        # block, which structurally rules out a loop. Best-effort: a guard
        # failure must never break the turn.
        # _turn_describes_intent no longer gates ENTRY. It exempts the
        # claim scanners inside, which is what it was measured for — a
        # plan names files that do not exist yet and grounding them turns
        # every plan into a false alarm. It never had anything to say
        # about the LANGUAGE a plan is written in, and gating entry on it
        # meant a persisted `default_mode: plan` switched off the whole
        # family. Measured on one sandbox, one line apart: with the
        # setting an English question came back in German, without it in
        # English.
        if (full_response and not _empty_turn and not self._stop_requested
                and not self._claim_guard_active):
            try:
                full_response = self._enforce_claim_grounding(
                    full_response,
                    memory_context=memory_context,
                    on_token=on_token,
                    on_tool_use=on_tool_use,
                    on_tool_result=on_tool_result,
                    on_permission_denied=on_permission_denied,
                    on_thinking=on_thinking,
                    thinking_budget=thinking_budget,
                    max_tokens=max_tokens,
                    on_notice=on_notice,
                    on_tool_result_meta=on_tool_result_meta,
                )
            except Exception:
                pass

        return full_response + _guard_note

    def _note_session_language(self, user_message: str) -> None:
        """Fix the session's language from the FIRST message that says one.

        The rule the user asked for, in their words: the language you start
        writing in is the one the session runs in. It is stronger than
        "the language of the latest message" for the case that actually
        hurt — a long unattended run where the answer arrives minutes
        later — because it can be STATED UP FRONT rather than checked
        afterwards.

        Set once. A later message in another language does not move it;
        that is the point of a session language, and it is what stops an
        English run drifting German halfway through.
        """
        if getattr(self, "_session_language", ""):
            return
        text = str(user_message or "")
        if text.lstrip().startswith("[Verify]"):
            return                      # the harness talking, always English
        try:
            from . import verify_guard as _vg
            found = _vg.detect_language(text)
        except Exception:
            found = ""
        if found:
            self._session_language = found

    def _session_language_block(self) -> str:
        """The one line that goes in front of the model, or ""."""
        want = str(getattr(self, "_session_language", "") or "")
        try:
            from . import verify_guard as _vg
            name = _vg._LANGUAGE_NAMES.get(want, "")
        except Exception:
            name = ""
        if not name:
            return ""
        return (
            f"SESSION LANGUAGE: {name}. Every answer in this session is "
            f"written in {name} — the first message set it, and a later "
            f"message in another language does not change it. This covers "
            f"what you SAY. What goes into code is English either way: "
            f"comments, docstrings, identifiers, log and error strings."
        )

    def _nudge_language_if_wrong(self, chunks) -> None:
        """Tell the model it is in the wrong language, mid-turn.

        The end-of-turn guard corrects the ANSWER. What a user reads FIRST
        is the sentence before the first tool call, and that one is on
        screen before anything has looked at it — measured in a field
        report where an English question got a German opening line, four
        tool rounds, and no final answer at all: the whole turn the user
        saw was the one sentence no guard could reach, and the report
        asked why it was still both languages.

        Nothing here can un-print that sentence. What it can do is stop
        the rest of the turn from continuing in it. The note rides the
        rail built for facts that changed under the model mid-run
        (``push_run_note`` → drained between tool rounds), so the model is
        told before it writes its next word rather than after its last.

        Once per turn. A second telling on every round would be nagging
        about a thing already said, and the model has to be allowed to act
        on the first one.
        """
        if getattr(self, "_language_note_sent", False):
            return
        try:
            asked = str(getattr(self, "_last_user_message", "") or "")
            # A correction turn's body is always English; reading one back
            # would demand English of every German session.
            if asked.lstrip().startswith("[Verify]"):
                return
            said = "".join(chunks or ())
            from . import verify_guard as _vg
            # Same source of truth as the end-of-turn guard.
            session_lang = str(getattr(self, "_session_language", "") or "")
            if session_lang:
                spoken = _vg.detect_language(said)
                want = session_lang if spoken and spoken != session_lang else ""
            else:
                want = _vg.scan_for_language_mismatch(said, asked)
            if not want:
                return
            note = _vg.language_mismatch_feedback(want).replace(
                "[Verify] ", "", 1)
            if not note:
                return
            self._language_note_sent = True
            self.client.push_run_note(note)
        except Exception:
            # A nudge that raises would cost the turn it was meant to
            # improve. It is advice, not a gate.
            pass

    def _turn_describes_intent(self) -> bool:
        """True when this turn proposes future work rather than asserting
        present state.

        A plan names the files it INTENDS to create; they do not exist yet
        by definition. Grounding those names against the workspace turns
        every plan into a false alarm and burns a correction turn on it
        (observed in the field). Plan mode, and any turn that submitted a
        plan, are therefore exempt — the guard's subject is claims about
        what IS, not proposals about what will be.
        """
        try:
            if str(getattr(self.kit_permissions, "mode", "") or "") == "plan":
                return True
            return any("exit_plan_mode" in str(t or "")
                       for t in (getattr(self, "_last_turn_tools", None) or ()))
        except Exception:
            return False

    # Result text that means the command did not run, or ran and failed.
    # A functional claim ("it works now") must not be able to cite any of
    # these as the run that proves it.
    # Deliberately narrow, and narrower than the first attempt at it.
    #
    # That attempt matched "traceback", "permission denied" and "no such
    # file or directory" as substrings anywhere in the output. A pytest
    # run with one failing test contains a traceback, so a run that
    # demonstrably happened was discarded and the guard then told an
    # accurate report that the file had never been exercised. Ordinary
    # stdout says "3 files skipped (Permission denied)" all the time.
    # A marker that fires on real output makes the guard louder on
    # CORRECT answers, which is the failure mode that teaches a model to
    # write around a guard instead of satisfying it.
    #
    # What is left are statements that the command itself did not run:
    # the shell could not find it, a gate refused it, a hook blocked it.
    # Whether what ran then FAILED is a different question, answered by
    # the exit code below.
    _EXEC_FAILURE_MARKERS = (
        "command not found", "blocked by hook",
        "is not available to the", "permission denied: cannot execute",
    )
    _EXIT_CODE_RE = re.compile(r"exit(?:\s+code)?[:= ]\s*([1-9]\d*)\b",
                               re.IGNORECASE)

    def _commit_exec_command(self, tool_name: str, tool_output: str) -> None:
        """Record a held command, but only if its result shows it ran.

        The ledger used to be written at the CALL. Everything that can go
        wrong between deciding to run something and it succeeding -- a
        permission gate, a hook block, a missing interpreter, a traceback,
        a non-zero exit -- was therefore invisible to it, and the guard
        that asks "was this artifact ever exercised?" found evidence for a
        command that never produced output. Best-effort; bookkeeping must
        never break a turn.
        """
        try:
            pending = self._exec_pending
            if not pending:
                return
            # Results arrive in call order, so the OLDEST unmatched call of
            # this name is the one this result belongs to. Taking the
            # newest was guaranteed wrong for a batched round: two bash
            # calls in one round, the first one's traceback arrives first
            # and pops the second one's entry -- so the crashed command
            # ended up in the ledger and the successful one was dropped.
            for idx, (pname, _) in enumerate(pending):
                if pname == tool_name:
                    name, tool_input = pending.pop(idx)
                    break
            else:
                return
            out = (tool_output or "").strip()
            if out.lstrip().startswith('{"error"'):
                return
            low = out.lower()
            if any(marker in low for marker in self._EXEC_FAILURE_MARKERS):
                return
            # A non-zero exit reported in the last few lines, where a
            # runner puts it. Searching the whole body would match a
            # program that merely PRINTS about exit codes.
            if self._EXIT_CODE_RE.search("\n".join(out.splitlines()[-3:])):
                return
            self._note_exec_command(name, tool_input)
        except Exception:
            pass

    def _note_exec_command(self, tool_name: str, tool_input: str) -> None:
        """Record an executed command in the session ledger.

        Only execution tools contribute (verify_guard.extract_exec_command
        decides); read-only job inspectors and every other tool are ignored.
        Bounded and de-duplicated so a long session cannot grow it without
        limit. Best-effort — bookkeeping must never break a turn."""
        try:
            from . import verify_guard as _vg
            cmd = _vg.extract_exec_command(tool_name, tool_input)
            if not cmd:
                return
            ledger = self._exec_commands_session
            if cmd in ledger:
                return
            ledger.append(cmd)
            if len(ledger) > 200:
                del ledger[:len(ledger) - 200]
        except Exception:
            pass

    def _note_ambiguous_columns(self, tool_output: str) -> None:
        """Record columns a document reader reported as undecidable.

        Per TURN, not per session: the question is whether THIS answer
        totals a column THIS turn could not read. Best-effort — the
        bookkeeping must never break the turn it observes."""
        try:
            from . import verify_guard as _vg
            found = _vg.extract_ambiguous_columns(tool_output)
            if not found:
                return
            ledger = getattr(self, "_ambiguous_columns_turn", None)
            if ledger is None:
                ledger = self._ambiguous_columns_turn = []
            for name in found:
                if name not in ledger:
                    ledger.append(name)
            if len(ledger) > 40:
                del ledger[:len(ledger) - 40]
        except Exception:
            pass

    def _note_truncated_tool(self, tool_name: str) -> None:
        """Record that a tool result was cut short this turn.

        Per TURN, like the ambiguous-column ledger and for the same
        reason: the question is whether THIS answer counts from output
        THIS turn could not see in full."""
        try:
            ledger = getattr(self, "_truncated_tools_turn", None)
            if ledger is None:
                ledger = self._truncated_tools_turn = []
            name = str(tool_name or "a tool")
            if name not in ledger:
                ledger.append(name)
            if len(ledger) > 20:
                del ledger[:len(ledger) - 20]
        except Exception:
            pass

    def _scan_truncated_counts(self, text: str) -> list:
        """Counts this answer states although its source was cut short."""
        try:
            from . import verify_guard as _vg
            return _vg.scan_for_counts_over_truncated_output(
                text, list(getattr(self, "_truncated_tools_turn", ()) or ()))
        except Exception:
            return []

    def _scan_ambiguous_column_totals(self, text: str) -> list:
        """Columns this answer totals although the reader could not read them."""
        try:
            from . import verify_guard as _vg
            return _vg.scan_for_totals_over_ambiguous_columns(
                text, list(getattr(self, "_ambiguous_columns_turn", ()) or ()))
        except Exception:
            return []

    def _append_self_consistency_caveat(
        self, text: str,
        on_token: Callable[[str], None] | None = None,
        *, source: str | None = None,
    ) -> str:
        """Mark an answer that contradicts its own list.

        The other two guards compare the answer to EVIDENCE -- a ledger of
        what was read, a note about what could not be. Neither compares it
        to itself, which is how "ich habe 31 PDF-Dateien verifiziert"
        survived above a list of 29 and a table restating 31: three
        inconsistent counts in one message, after both layers had run.

        Caveat and not correction, and naming BOTH numbers rather than
        picking one. A retry cannot fix a counting error the model just
        made twice, and the framework does not know which figure is the
        right one -- only that they disagree.

        ``source`` is the text to SCAN when it differs from the text to
        append to: at the single exit the earlier caveats are already on
        the answer, and a guard must never read another guard's output as
        if the model had written it."""
        try:
            from . import verify_guard as _vg
            pairs = _vg.scan_for_count_vs_enumeration(
                source if source is not None else text)
            caveat = _vg.count_vs_enumeration_caveat(pairs)
        except Exception:
            return text
        if not caveat:
            return text
        if on_token:
            try:
                on_token(caveat)
            except Exception:
                pass
        return text + caveat

    def _append_truncated_count_caveat(
        self, text: str,
        on_token: Callable[[str], None] | None = None,
        *, source: str | None = None,
    ) -> str:
        """Mark a count whose only source was cut short.

        Same consequence as the ambiguous-column caveat and for the same
        reason: the number is not demonstrably wrong, it is unfounded. A
        correction turn cannot help either -- the full listing is still not
        in context, so a retry would count the same truncated output
        again. In the field this produced "31 PDF-Dateien verifiziert"
        followed by a list of 29, after the framework's own verify nudge
        had already fired once.

        ``source`` is the text to SCAN when it differs from the text to
        append to — see _append_self_consistency_caveat."""
        try:
            counts = self._scan_truncated_counts(
                source if source is not None else text)
            if not counts:
                return text
            from . import verify_guard as _vg
            caveat = _vg.truncated_output_caveat(
                counts, list(getattr(self, "_truncated_tools_turn", ()) or ()))
        except Exception:
            return text
        if not caveat:
            return text
        if on_token:
            try:
                on_token(caveat)
            except Exception:
                pass
        return text + caveat

    def _append_ambiguous_column_caveat(
        self, text: str, columns: list,
        on_token: Callable[[str], None] | None = None,
    ) -> str:
        """Mark a figure that rests on a guessed reading, in the answer itself.

        A caveat rather than a forced correction: the number is not wrong,
        it is unfounded, and the reader needs it labelled rather than
        withheld. A retry would also be free to make the same assumption
        again — the tool result and the role prompt both already said not
        to, and neither bound. Appending puts the warning where the figure
        is, and into the transcript, so a later turn sees it marked."""
        if not columns:
            return text
        try:
            from . import verify_guard as _vg
            caveat = _vg.ambiguous_column_caveat(columns)
        except Exception:
            return text
        if not caveat:
            return text
        if on_token:
            try:
                on_token(caveat)
            except Exception:
                pass
        try:
            if (self.messages
                    and self.messages[-1].get("role") == "assistant"
                    and isinstance(self.messages[-1].get("content"), str)):
                self.messages[-1]["content"] += caveat
        except Exception:
            pass
        return text + caveat

    def _scan_functional_claims(self, text: str) -> list:
        """Scan a finished answer for claims that something now works which
        this session never exercised. Returns the flag list; never raises."""
        try:
            from . import verify_guard as _vg
            cmds = getattr(self, "_exec_commands_session", None)
            return _vg.scan_for_unexercised_functional_claims(
                text,
                exec_commands=list(cmds or ()),
                exec_ledger_available=cmds is not None,
                tools_used=set(getattr(self, "_last_turn_tools", None) or ()),
            )
        except Exception:
            return []

    def _append_functional_caveat(
        self, text: str, flags: list,
        on_token: Callable[[str], None] | None = None,
    ) -> str:
        """Append the visible functional-claim caveat to a finished answer.

        This class gets a caveat, never a forced correction turn: the model
        cannot verify browser or interactive behavior headlessly either, so
        a retry would only invite a second confident assertion. The caveat
        also enters the transcript, so later turns see the claim marked
        unconfirmed instead of standing bare."""
        if not flags:
            return text
        try:
            from . import verify_guard as _vg
            caveat = _vg.functional_claim_caveat(flags)
        except Exception:
            return text
        if not caveat:
            return text
        if on_token:
            try:
                on_token(caveat)
            except Exception:
                pass
        try:
            if (self.messages
                    and self.messages[-1].get("role") == "assistant"
                    and isinstance(self.messages[-1].get("content"), str)):
                self.messages[-1]["content"] += caveat
        except Exception:
            pass
        return text + caveat

    def _append_answer_caveats(
        self, text: str, *, functional: list, ambiguous: list,
        on_token: Callable[[str], None] | None = None,
        scan_source: str | None = None,
    ) -> str:
        """Apply the whole caveat chain ONCE, at the guard's single exit.

        Each of the four exits used to append its own subset, so a single
        location or quantity flag silently suppressed the truncated-count
        and the count-vs-own-enumeration guards — the two least related to
        it. ``ambiguous`` was even computed and then dropped on two of
        them. An answer is caveated for what it says, not for which branch
        of the guard happened to return it.

        The later scanners read the model's own text, never the caveats
        added ahead of them in this chain — ``scan_source`` names it when
        a guard note is already appended to ``text``."""
        source = scan_source if scan_source is not None else text
        out = self._append_functional_caveat(text, functional, on_token)
        out = self._append_ambiguous_column_caveat(out, ambiguous, on_token)
        out = self._append_truncated_count_caveat(out, on_token, source=source)
        out = self._append_self_consistency_caveat(out, on_token, source=source)
        out = self._append_figure_coverage_caveat(out, on_token, source=source)
        return out

    def _append_figure_coverage_caveat(
        self, text: str, on_token: Callable[[str], None] | None,
        *, source: str,
    ) -> str:
        """Name a money figure or a count the tools never produced.

        The office path had no mechanism at all: a total that silently
        dropped a row read exactly like a total that did not, and the
        quantity scanner only knows the chemistry units — eV, kcal/mol —
        so a currency amount and a document count passed with zero
        evidence acts behind them.

        Deliberately quiet. Measured against twelve answers a careful
        person would write about the fixture workbooks: none flagged. It
        reads the model's own text, never the caveats appended ahead of
        it in this chain.
        """
        try:
            from . import office as _office_ledger
            note = _office_ledger.figure_coverage_caveat(
                source,
                user_text=getattr(self, "_last_user_message", "") or "",
                # This turn's ledger, named: the answer is checked against
                # what ITS tools produced, never against whatever else the
                # process happened to record.
                token=getattr(self, "_figure_ledger_token", "") or "",
            )
        except Exception:
            return text
        if not note:
            return text
        if on_token:
            try:
                on_token(note)
            except Exception:
                pass
        return text + note

    def _scan_claim_grounding(
        self, text: str, turn_tools: list[str] | None,
    ) -> tuple[list, list]:
        """Scan a finished answer for ungrounded code-location and
        physical-quantity claims against this session's observed evidence.
        Returns (location_flags, quantity_flags); never raises."""
        try:
            from . import verify_guard as _vg
            observed = getattr(self, "_last_observed_files", None) or set()
            ledger = bool(getattr(self, "_observed_ledger_available", False))
            loc = _vg.scan_for_ungrounded_location_claims(
                text, observed_files=observed, ledger_available=ledger)
            qty = _vg.scan_for_unsourced_quantities(
                text, observed_files=observed,
                evidence_tools_used=set(turn_tools or ()),
                numbers=_vg.observed_numbers())
            return loc, qty
        except Exception:
            return [], []

    def _enforce_claim_grounding(
        self,
        response_text: str,
        *,
        memory_context: str = "",
        on_token: Callable[[str], None] | None = None,
        on_tool_use: Callable[[str, str], None] | None = None,
        on_tool_result: Callable[[str, str], None] | None = None,
        on_permission_denied: Callable[[str], None] | None = None,
        on_thinking: Callable[[str], None] | None = None,
        thinking_budget: int = 0,
        max_tokens: int = 0,
        on_notice: Callable[[str], None] | None = None,
        on_tool_result_meta: Callable[[str, dict], None] | None = None,
    ) -> str:
        """Enforce evidence grounding on a finished answer (all modes).

        Ungrounded specific claims (code locations / measured quantities,
        per verify_guard) trigger exactly one forced correction turn that
        asks the model to verify with tools or restate with explicit
        uncertainty. One retry only — if the correction is still
        ungrounded, a visible caveat is appended instead of failing or
        looping.

        Functional claims ("it works now") run through the SAME gate but
        take the other consequence: a visible caveat naming what was never
        exercised, and no correction turn — see _append_functional_caveat.

        Whether the correction WORKED is decided by re-scanning it against
        the observed-files ledger, never by the fact that a correction ran:
        a retry that only rephrased the claim with hedges leaves it exactly
        as unverified as it was, and is reported as such.

        Returns the (possibly extended) answer text."""
        from . import verify_guard as _vg
        loc, qty = self._scan_claim_grounding(
            response_text, getattr(self, "_last_turn_tools", None))
        func = self._scan_functional_claims(response_text)
        # A figure over a column the reader could not read joins the caveat
        # class rather than the correction class: the number is not wrong,
        # it is unfounded, and a retry is free to guess the same way again.
        ambiguous = self._scan_ambiguous_column_totals(response_text)
        # Two values for one field of one record. This one does not judge
        # the ANSWER at all — it judges what the turn's tools returned, so
        # it fires even when the answer only reports aggregates, which is
        # exactly how the field case hid: the chat gave min/max/mean and
        # the contradiction lived between a stored file and a table the
        # agent had just written.
        conflicts = _vg.scan_for_conflicting_figures()
        # A plan proposes; it does not assert. Every scanner about what IS
        # stands down for it — see _turn_describes_intent — and only the
        # language question survives, because a plan is still an answer to
        # somebody who wrote in a particular language.
        if self._turn_describes_intent():
            loc, qty, conflicts, func, ambiguous = [], [], [], [], []
        # The language the user wrote in. Measured to be necessary: with
        # the prompt rule delivered and nothing German in the prompt, the
        # backend answered an English question in German four times in
        # four runs. Returns "" whenever either side is too short or too
        # technical to say, so it is silent far more often than it fires.
        # The engine already records what the user wrote, at the top of
        # the turn. A "[Verify] …" body means a correction turn has since
        # overwritten it — those are always English, and reading one back
        # would demand English of every German session, so it is ignored.
        # ONE source of truth. The session language is set by the first
        # message and stated in the prompt; judging the answer against the
        # LATEST message instead would let the two disagree, and a rule
        # whose two halves disagree is the defect this project keeps
        # finding. A session that never got one (an older restore) falls
        # back to the per-message question rather than to silence.
        _session_lang = str(getattr(self, "_session_language", "") or "")
        if _session_lang:
            _said = _vg.detect_language(response_text)
            wrong_language = (_session_lang
                              if _said and _said != _session_lang else "")
        else:
            _asked = str(getattr(self, "_last_user_message", "") or "")
            wrong_language = _vg.scan_for_language_mismatch(
                response_text,
                "" if _asked.lstrip().startswith("[Verify]") else _asked)
        if not loc and not qty and not conflicts and not wrong_language:
            return self._append_answer_caveats(
                response_text, functional=func, ambiguous=ambiguous,
                on_token=on_token)
        if self._claim_guard_spent:
            # The single correction for this user turn is spent (e.g. a
            # nested continuation re-entered here) — annotate, never loop.
            return self._append_answer_caveats(
                response_text + _vg.grounding_caveat(loc, qty)
                + _vg.conflicting_figure_caveat(conflicts)
                + _vg.language_mismatch_caveat(wrong_language),
                functional=func, ambiguous=ambiguous, on_token=on_token)
        parts: list[str] = []
        if loc:
            parts.append(_vg.location_claim_feedback(loc))
        if qty:
            parts.append(_vg.quantity_claim_feedback(qty))
        if conflicts:
            parts.append(_vg.conflicting_figure_feedback(conflicts)
                         .replace("[Verify] ", "", 1))
        if wrong_language:
            parts.append(_vg.language_mismatch_feedback(wrong_language)
                         .replace("[Verify] ", "", 1))
        feedback = "[Verify] " + " ".join(parts)
        # What the session had observed BEFORE the retry. The retry either
        # reads something new or it does not, and that is the only thing
        # that can turn an unverified claim into a verified one.
        observed_before = set(getattr(self, "_last_observed_files", None) or ())
        self._claim_guard_spent = True
        self._claim_guard_active = True
        if on_token:
            try:
                on_token("\n\n")
            except Exception:
                pass
        try:
            correction = self.stream_response(
                user_message=feedback,
                memory_context=memory_context,
                on_token=on_token,
                on_tool_use=on_tool_use,
                on_tool_result=on_tool_result,
                on_permission_denied=on_permission_denied,
                on_thinking=on_thinking,
                thinking_budget=thinking_budget,
                max_tokens=max_tokens,
                on_notice=on_notice,
                on_tool_result_meta=on_tool_result_meta,
            )
        except Exception:
            correction = ""
        finally:
            self._claim_guard_active = False
        if not correction:
            return self._append_answer_caveats(
                response_text + _vg.grounding_caveat(loc, qty)
                + _vg.conflicting_figure_caveat(conflicts)
                + _vg.language_mismatch_caveat(wrong_language),
                functional=func, ambiguous=ambiguous, on_token=on_token)
        # Appending is right when the correction ADDS something — the
        # claim and its correction both belong on the page. A language
        # correction replaces instead: the reader asked in English, and
        # showing them the German first and the English underneath gives
        # them the thing they did not ask for plus the thing they did.
        # Only when language was the SOLE reason, so a grounding fix is
        # never swallowed along with it.
        if wrong_language and not (loc or qty or conflicts):
            combined = correction
        else:
            combined = response_text + "\n\n" + correction
        # Re-scan the correction: the recursive turn refreshed the
        # observed-files snapshot and _last_turn_tools.
        loc2, qty2 = self._scan_claim_grounding(
            correction, getattr(self, "_last_turn_tools", None))
        # The correction may restate the functional claim — scan it too and
        # merge (order-stable, de-duplicated: the flags are frozen).
        func = list(dict.fromkeys(
            func + self._scan_functional_claims(correction)))
        ambiguous = self._scan_ambiguous_column_totals(combined)
        new_files = set(
            getattr(self, "_last_observed_files", None) or ()) - observed_before
        if loc2 or qty2:
            # The retry produced its own ungrounded claims.
            caveat = _vg.grounding_caveat(loc2, qty2)
        elif not new_files:
            # It read nothing. The scanners are silent because the wording
            # changed, not because anything was checked — so the ORIGINAL
            # claims are named, exactly as unverified as before. The
            # feedback offers "restate as unverified" as a way out, and
            # taking it must not be indistinguishable from verifying.
            caveat = _vg.grounding_caveat(loc, qty)
        else:
            caveat = ""
            self._claim_guard_corrected = True
        # A conflict is resolved by NAMING a number, not by acknowledging
        # the question. The field case answered its own comparison with
        # "deutet auf eine unterschiedliche Formel-Konvention hin" and
        # delivered its values unchanged — prose that would satisfy any
        # check made of words. Requiring one of the two figures to appear
        # in the correction cannot be met that way.
        unresolved = [c for c in conflicts
                      if not _vg.conflict_is_addressed(c, correction)]
        if unresolved:
            caveat += _vg.conflicting_figure_caveat(unresolved)
        # Judged on the CORRECTION, so switching resolves it and restating
        # in the same language does not. Silence counts as resolved: a
        # short correction cannot be judged, and a caveat built on a
        # verdict the detector refused to give would be this guard doing
        # what the project keeps catching it doing.
        if wrong_language and _vg.detect_language(correction) not in (
                "", wrong_language):
            caveat += _vg.language_mismatch_caveat(wrong_language)
        marker = _vg.verification_marker(new_files) if caveat == "" else ""
        note = caveat or marker
        if note:
            if on_token:
                try:
                    on_token(note)
                except Exception:
                    pass
            # Record it in the transcript too, so later turns see the claim
            # marked rather than standing bare.
            try:
                if (self.messages
                        and self.messages[-1].get("role") == "assistant"
                        and isinstance(self.messages[-1].get("content"), str)):
                    self.messages[-1]["content"] += note
            except Exception:
                pass
        return self._append_answer_caveats(
            combined + note, functional=func, ambiguous=ambiguous,
            on_token=on_token, scan_source=combined)

    def trace_session(self) -> str:
        """Stable key for this engine's tool-call trace — the backend session
        id when present, else a per-engine uuid (so OpenAI/KIT/Ollama turns,
        which have no server session id, still get a stable trace file)."""
        return self.session_id or self._trace_id

    @staticmethod
    def _tool_result_failed(output: str) -> tuple[bool, str]:
        """Whether a tool result reports a failure, and the reason.

        Every executor returns a refusal or an error as a JSON object with
        an ``error`` key -- that is the one shape the whole tool layer
        agrees on, which is why the grounding ledgers key on it too. Read
        it here rather than trusting the caller's assertion.

        Deliberately NOT a judgement about the content: a command that ran
        and exited non-zero is a successful tool call reporting a failed
        command, and conflating the two would make the error rate measure
        the user's code instead of the agent's tooling.
        """
        text = (output or "").lstrip()
        if not text.startswith('{"error"'):
            return False, ""
        try:
            reason = str((json.loads(output) or {}).get("error", ""))
        except Exception:
            reason = text[:200]
        return True, reason[:300]

    def _record_tool_trace(
        self, name: str, output: str = "", *, ok: bool = True, error: str = "",
    ) -> None:
        """Pair the latest tool_use with this result and append a trace entry."""
        try:
            t0 = None
            inp = ""
            for i, (pn, pi, pt) in enumerate(self._trace_pending):
                if pn == name:
                    inp, t0 = pi, pt
                    self._trace_pending.pop(i)
                    break
            else:
                if self._trace_pending:
                    _, inp, t0 = self._trace_pending.pop(0)
            dur = int((_time.monotonic() - t0) * 1000) if t0 else 0
            from .tool_trace import record as _rec
            _rec(self.trace_session(), tool=name, tool_input=inp,
                 output=output, duration_ms=dur, ok=ok, error=error)
        except Exception:
            pass

    def _build_user_message(
        self, text: str, images: list[str] | None,
    ) -> dict[str, Any]:
        """Build the user message dict.

        When images are attached AND the active model is vision-capable on an
        OpenAI-compatible backend (Ollama / KIT / OpenAI), the pixels are sent
        as multimodal ``image_url`` content so the model actually SEES them.
        Otherwise a plain-text message (the caller adds a note about the
        attached files). The CLI backend uses a different image format → text fallback.
        """
        if not images:
            return {"role": "user", "content": text}
        try:
            if self.provider not in ("openai", "kit", "ollama"):
                return {"role": "user", "content": text}
            from .image_input import (
                load_image, model_supports_vision, to_openai_content)
            model = getattr(self.client, "model", "") or ""
            caps = getattr(self, "_active_capabilities", None)
            if not model_supports_vision(model, caps):
                return {"role": "user", "content": text}
            loaded = []
            for p in images:
                try:
                    loaded.append(load_image(p))
                except Exception:
                    continue
            if not loaded:
                return {"role": "user", "content": text}
            content = to_openai_content(text, loaded, model=model, caps=caps)
            return {"role": "user", "content": content}
        except Exception:
            return {"role": "user", "content": text}

    def _sanitize_messages(self) -> None:
        """Ensure message history has proper user/assistant alternation.

        Concurrent stop/send can leave consecutive user messages or other
        structural issues.  This removes duplicates and ensures the
        conversation can be sent to the API without errors.

        Pinned messages (``msg["_pinned"]``) are never merged away: when a
        pinned message would be replaced by a newer same-role message, a
        short opposite-role filler is inserted instead so BOTH survive
        verbatim and alternation still holds for strict chat templates.
        """
        if len(self.messages) < 2:
            return
        cleaned: list[dict] = [self.messages[0]]
        for msg in self.messages[1:]:
            role = msg.get("role")
            if role == cleaned[-1].get("role") and role in ("user", "assistant"):
                if cleaned[-1].get("_pinned"):
                    # The older message is pinned — keep it verbatim and
                    # restore alternation with a filler instead of merging.
                    filler_role = "assistant" if role == "user" else "user"
                    cleaned.append({
                        "role": filler_role,
                        "content": "[alternation filler - pinned message "
                                   "above kept verbatim]",
                    })
                    cleaned.append(msg)
                else:
                    # Consecutive same-role messages: keep the latest
                    # (also keeps a pinned NEWER message automatically).
                    cleaned[-1] = msg
            else:
                cleaned.append(msg)
        if cleaned != self.messages:
            self.messages[:] = cleaned

    # -- Pinned context regions --------------------------------------------
    #
    # A pin (``msg["_pinned"] = True``) protects one message from EVERY
    # destructive compaction stage: the sliding-window trim, the hard-clear
    # of old tool results, and the full summary compaction (the message is
    # carried verbatim into the kept set instead of being summarised). The
    # flag lives on the message dict itself so it survives export_state /
    # save_session / load_session / restore_state unchanged; the backend
    # request path strips private keys (see _wire_messages) so providers
    # never see it.

    def pin_message(self, index: int) -> bool:
        """Pin the message at ``index`` so compaction never alters it.

        ``index`` addresses ``self.messages``; negative values resolve
        Python-style (``-1`` = newest). Returns True when the pin was set,
        False for an out-of-range or malformed index/message.
        """
        msgs = getattr(self, "messages", None) or []
        try:
            idx = int(index)
        except (TypeError, ValueError):
            return False
        if idx < 0:
            idx += len(msgs)
        if not (0 <= idx < len(msgs)) or not isinstance(msgs[idx], dict):
            return False
        msgs[idx]["_pinned"] = True
        return True

    def unpin_message(self, index: int) -> bool:
        """Remove the pin at ``index``.

        Returns True when a pin was actually removed; False when the index
        is out of range or the message was not pinned.
        """
        msgs = getattr(self, "messages", None) or []
        try:
            idx = int(index)
        except (TypeError, ValueError):
            return False
        if idx < 0:
            idx += len(msgs)
        if not (0 <= idx < len(msgs)) or not isinstance(msgs[idx], dict):
            return False
        return msgs[idx].pop("_pinned", None) is not None

    def pinned_indices(self) -> list[int]:
        """Indices of all currently pinned messages, ascending."""
        msgs = getattr(self, "messages", None) or []
        return [
            i for i, m in enumerate(msgs)
            if isinstance(m, dict) and m.get("_pinned")
        ]

    def _wire_messages(self) -> list[dict[str, Any]]:
        """The message list as sent to the backend.

        Private (underscore-prefixed) bookkeeping keys such as ``_pinned``
        must never reach a provider: the Anthropic backend forwards the
        dicts verbatim and strict APIs reject unknown fields. Returns
        ``self.messages`` itself when no message carries a private key (the
        common case — preserves list identity for backends that bind the
        live list); otherwise a shallow per-message copy with private keys
        stripped (same length and indices, so ``live:<i>`` refs stay valid).
        """
        msgs = self.messages
        if not any(
            isinstance(m, dict) and any(str(k).startswith("_") for k in m)
            for m in msgs
        ):
            return msgs
        out: list[dict[str, Any]] = []
        for m in msgs:
            if isinstance(m, dict) and any(str(k).startswith("_") for k in m):
                out.append({k: v for k, v in m.items()
                            if not str(k).startswith("_")})
            else:
                out.append(m)
        return out

    def _elide_original(
        self, index: int, role: str, content: str, *, reason: str,
    ) -> str:
        """Best-effort lossless-elision write before a destructive trim.

        Appends the FULL original ``content`` to this session's elided
        store (``<sessions_dir>/<session_id>.elided.jsonl``) and returns
        the short ref id for the replacement marker, so the model can page
        the dropped text back in via ``history_get('elided:<ref>')``.
        Returns ``""`` on any failure or when no session id exists — the
        caller then keeps its plain marker and compaction proceeds
        unchanged (an elision-store failure must never break compaction).
        """
        try:
            sid = str(getattr(self, "session_id", "") or "")
            if not sid:
                return ""
            from delfin.agent.session_store import append_elided_record
            ref = append_elided_record(
                sid, index=index, role=role, content=content, reason=reason,
            )
            return ref or ""
        except Exception:
            return ""

    # -- Mid-conversation compaction (Feature 3) ---------------------------

    # The legacy message-count trigger. It is NOT a gate any more: pressure
    # decides, and the only count that still matters is the structural floor
    # (_KEEP_RECENT) below which there is nothing between the summary block
    # and the kept tail. Retained as the "a conversation this short is not
    # worth summarising" reference point.
    # Retired. It gated _compact_history until pressure was asked
    # first; since the count can only matter once the budget is
    # already over the line, its whole remaining effect was to block
    # compaction exactly when it was needed. The floor that survives
    # is _KEEP_RECENT: there must be something to summarise beyond
    # the messages that are kept.
    _KEEP_RECENT = 4            # keep last 4 messages intact
    # Header that marks a compaction summary message — used both to build the
    # block and to recognise a PRIOR summary on re-compaction so it isn't
    # re-truncated like a fresh user goal (which compounded loss).
    _SUMMARY_BLOCK_PREFIX = "[Conversation summary"

    def _estimate_context_tokens(self) -> int:
        """Estimate the tokens in the NEXT request, not just self.messages.

        The system prompt (repo map, memory, playbook, live-state, tasks) and
        the advertised tool schemas ride on every request but are not in
        self.messages. Counting only the messages undercounts the budget, so a
        request could satisfy ``auto_compact_pct`` yet still blow the provider
        window (400 context_length_exceeded). So fold in the last built system
        prompt, and never estimate below the last real tokenizer count the
        provider reported for this context (that count already includes the
        tool schemas, so no separate schema constant is needed once a turn has
        run — and the first turn is far too short to compact anyway).
        """
        total_chars = 0
        for m in self.messages:
            content = m.get("content", "")
            if isinstance(content, str):
                total_chars += len(content)
            elif isinstance(content, list):
                for part in content:
                    if not isinstance(part, dict):
                        continue
                    if part.get("text"):
                        total_chars += len(str(part["text"]))
                    elif (part.get("type") in ("image", "image_url")
                          or "image_url" in part or "source" in part):
                        # An image is invisible to char/4 but real to the model;
                        # count a fixed ~1500-token allowance instead of 0.
                        total_chars += 6000
                    else:
                        total_chars += len(str(part))
        est = total_chars // 4
        # The live prompt text is authoritative while a process is running;
        # the persisted size stands in only after a resume, when the text is
        # gone (it is not written into session files -- it carries the
        # injected memory) but its weight in the window is unchanged.
        _sp_live = len(getattr(self, "last_system_prompt", "") or "")
        est += (_sp_live or int(getattr(self, "_system_prompt_chars", 0) or 0)) // 4
        # The provider floor is a PRE-trim snapshot; credit chars that
        # in-place trims removed since it was taken, or the trim loops'
        # stop conditions ("estimate <= budget") could never be reached
        # and the paid LLM summary would always fire at the 95% cliff.
        floor = int(getattr(self, "_last_input_tokens", 0) or 0)
        floor -= int(getattr(self, "_trimmed_chars_since_floor", 0) or 0) // 4
        return max(est, max(0, floor))

    def _should_auto_compact(self) -> bool:
        """True if estimated context exceeds the configured fraction."""
        if not self.auto_compact_pct or not self.context_window_tokens:
            return False
        budget = int(self.context_window_tokens * float(self.auto_compact_pct))
        return self._estimate_context_tokens() > budget

    # Sliding-window threshold (proactive trim long before the auto-compact
    # cliff at auto_compact_pct=0.95). Default 0.70 of the window — every
    # message past that triggers a gentle in-place trim of the oldest
    # tool_result first. The full auto-compact still fires at 0.95.
    _SLIDING_WINDOW_PCT = 0.70

    def _should_slide(self) -> bool:
        """Trigger the gentler sliding-window trim BEFORE auto_compact."""
        if not self.context_window_tokens:
            return False
        budget = int(self.context_window_tokens * self._SLIDING_WINDOW_PCT)
        return self._estimate_context_tokens() > budget

    def _irreducible_tokens(self) -> int:
        """The part of the estimate that trimming messages cannot change.

        The system prompt, which is rebuilt every turn and does not shrink
        when history does. Kept separate because the trim loops compare
        against a budget covering the whole estimate, and a constant on the
        wrong side of that comparison is what made them unstoppable.
        """
        live = len(getattr(self, "last_system_prompt", "") or "")
        return (live or int(getattr(self, "_system_prompt_chars", 0) or 0)) // 4

    def _trimming_cannot_reach(self, budget: int, kind: str) -> bool:
        """True when the target is out of reach before anything is trimmed.

        Measured on a 32k window with a 90 kB prompt: budget 22400, system
        prompt 22500 on its own. The loops' break condition -- estimate
        under budget -- was therefore unsatisfiable, so they ran to the end
        of the eligible prefix and cut every message in it down to a stub,
        finished still over budget, and did it again the next turn on
        whatever had grown since. The agent's own answers are the only
        place its conclusions and file lists survive.

        Shredding the history does not fix a prompt that is too large for
        the window; it only costs the session its memory of itself. The
        remedy is a smaller prompt, and which part to give up is the user's
        decision, so this is recorded rather than done quietly.

        The record carries ``archived_at`` like every other compaction
        record. Without it the post-turn notice — which keys on exactly
        that field to detect "something happened this turn" — never fired,
        so the one diagnosis that needs the USER to act was the only one
        never shown to them.
        """
        irreducible = self._irreducible_tokens()
        if irreducible < budget:
            return False
        self.last_compaction_info = {
            "kind": kind,
            "trimmed": 0,
            "messages_compacted": 0,
            "tokens_saved": 0,
            "archived_at": _time.time(),
            "note": (
                f"the system prompt alone is ~{irreducible} tokens against a "
                f"{budget}-token budget for this window, so trimming the "
                f"conversation cannot reach it. History left intact — the "
                f"prompt is what has to get smaller."
            ),
        }
        return True

    # Machine-generated turns carry the "user" role in ``self.messages``.
    #
    # Enumerating every append site shows this history only ever holds the
    # roles "user" and "assistant": tool calls and their results live in a
    # per-request list inside the backend client and are discarded when the
    # request returns. What DOES reach the history in place of a tool result
    # is a synthetic user turn -- command output fed back to the model,
    # verification feedback, injected context sections. Those are the bulk,
    # and treating them as GOALS (which "role == user" did) meant the trim
    # stages skipped them and cut the agent's own answers instead: the only
    # place its conclusions, file lists and line numbers survive.
    _MACHINE_TURN_PREFIXES = (
        "[Command results]",
        "[Verify]",
        "[System]",
    )

    @classmethod
    def _is_machine_turn(cls, msg: dict) -> bool:
        """True for a synthetic user turn — this history's tool output."""
        if not isinstance(msg, dict) or msg.get("role") != "user":
            return False
        content = msg.get("content")
        if not isinstance(content, str):
            return False
        return content.lstrip().startswith(cls._MACHINE_TURN_PREFIXES)

    def _trim_candidates(self, msgs: list[dict]) -> list[tuple[int, dict]]:
        """(index, message) pairs a destructive stage may touch, cheapest
        loss first: machine turns before the agent's own answers, oldest
        first within each group. A real user goal and a pinned message are
        in neither group and are therefore never returned.
        """
        machine: list[tuple[int, dict]] = []
        own: list[tuple[int, dict]] = []
        for idx, msg in enumerate(msgs):
            if not isinstance(msg, dict) or msg.get("_pinned"):
                continue
            if self._is_machine_turn(msg):
                machine.append((idx, msg))
            elif msg.get("role", "") != "user":
                own.append((idx, msg))
        return machine + own

    def _shorten_oldest_non_goal_messages(self) -> int:
        """In-place trim: progressively shorten the OLDEST messages that are
        not user goals until we're back under the sliding-window threshold.
        Machine turns go first, then the agent's own answers; real user
        goals and pinned messages are never touched.

        Returns the number of messages that were trimmed (not removed —
        just shortened with an ``... [trimmed by sliding window] ...``
        middle marker so the agent still sees the head + tail).
        """
        budget = int(self.context_window_tokens * self._SLIDING_WINDOW_PCT)
        if self._trimming_cannot_reach(budget, "sliding_window"):
            return 0
        trimmed = 0
        # Don't touch the last KEEP_RECENT messages — they're conversation
        # state the agent is actively reasoning over.
        protected_from = max(0, len(self.messages) - self._KEEP_RECENT)
        for idx, msg in self._trim_candidates(self.messages[:protected_from]):
            if self._estimate_context_tokens() <= budget:
                break
            role = msg.get("role", "")
            content = msg.get("content", "")
            if not isinstance(content, str):
                continue
            # Only trim large messages — small ones (short goals, short
            # acknowledgements) survive untouched.
            if len(content) < 2000:
                continue
            # Lossless elision: persist the FULL original before trimming so
            # the dropped middle stays retrievable. Best-effort — an empty
            # ref just means the marker carries no retrieval hint.
            ref = self._elide_original(
                idx, role, content, reason="sliding_window")
            hint = f" — retrievable via history_get('elided:{ref}')" if ref else ""
            # Keep the head + tail of the oldest long non-goal message.
            head_len = 600
            tail_len = 400
            new_content = (
                content[:head_len]
                + f"\n... [trimmed by sliding window, "
                + f"{len(content) - head_len - tail_len} chars dropped"
                + f"{hint}] ...\n"
                + content[-tail_len:]
            )
            msg["content"] = new_content
            self._trimmed_chars_since_floor = (
                getattr(self, "_trimmed_chars_since_floor", 0)
                + max(0, len(content) - len(new_content)))
            trimmed += 1
        return trimmed

    def _stub_oldest_non_goal_messages(self, old_msgs: list[dict]) -> int:
        """Deterministic context-editing pass: stub the body of old, bulky
        NON-GOAL messages down to a short marker so the lossy +
        non-deterministic + paid LLM summary becomes a genuine LAST RESORT
        rather than the default at the 95% cliff.

        Structure-preserving: machine turns first, then the agent's own
        answers, and only string content (real user goals are the GOALS;
        list-content tool_use blocks and their tool_call pairing are left
        intact), so the message list stays valid for the backend. Returns
        how many messages were cleared.

        Runs AFTER the 0.70 sliding-window trim, so it only bites when a
        session is genuinely at the compaction cliff; then, if it frees enough,
        ``_compact_history`` skips the summary entirely.
        """
        pct = float(self.auto_compact_pct or 0.95)
        budget = int(self.context_window_tokens * pct)
        # Same rule as the sliding window: a budget the prompt already
        # exceeds cannot be reached by stubbing messages, and running the
        # loop anyway empties the history for nothing.
        if self._trimming_cannot_reach(budget, "hard_clear"):
            return 0
        cleared = 0
        for idx, msg in self._trim_candidates(old_msgs):
            if self._estimate_context_tokens() <= budget:
                break
            content = msg.get("content", "")
            if not isinstance(content, str) or len(content) < 800:
                continue
            if content.startswith("[cleared:"):
                continue  # already stubbed on an earlier compaction
            # Lossless elision: keep the full original retrievable. old_msgs
            # is a prefix slice of self.messages, so idx is the live index.
            ref = self._elide_original(
                idx, msg.get("role", ""), content, reason="hard_clear")
            hint = f" — retrievable via history_get('elided:{ref}')" if ref else ""
            head = content[:200]
            msg["content"] = (
                f"[cleared: {len(content)} chars of an older message body "
                f"elided to save context{hint}]\n{head}"
            )
            self._trimmed_chars_since_floor = (
                getattr(self, "_trimmed_chars_since_floor", 0)
                + max(0, len(content) - len(msg["content"])))
            cleared += 1
        return cleared

    def _llm_summarize_old_messages(self, old_msgs: list[dict]) -> str:
        """Call a cheap model to summarise the older half of the
        conversation into a high-fidelity recap. Falls back to ``""`` on
        any failure so the caller can drop back to extractive mode.

        The summarisation prompt asks for: user goals, key decisions,
        files touched, and pending items; drops chit-chat. The cost is
        ~1k tokens of input + a few hundred output; for an Opus session
        that's <$0.01.

        Backend handling: API-style clients (APIClient, OpenAIClient)
        are stateless per stream_message call — safe to invoke for an
        out-of-band summary. CLI-style clients (CLIClient, CodexCLIClient)
        use a *persistent* subprocess with the user's session_id; sending
        the summarisation prompt there would pollute the live transcript.
        For CLI backends we therefore decline LLM-mode and the caller
        falls back to extractive summarisation.
        """
        if not old_msgs or not getattr(self, "client", None):
            return ""
        if getattr(self, "backend", "") == "cli":
            # Persistent subprocess — running the summariser would append
            # the side-conversation to the user's live session. Skip.
            return ""

        rendered: list[str] = []
        carried_summary = ""   # a prior compaction's recap, carried forward verbatim
        for msg in old_msgs:
            role = msg.get("role", "")
            content = msg.get("content", "")
            if isinstance(content, str) and content.lstrip().startswith(
                    self._SUMMARY_BLOCK_PREFIX):
                # A prior summary block is ALREADY a faithful recap of many
                # messages. Re-summarising it compounds loss across repeated
                # compactions on a long session, so carry it forward verbatim
                # (header stripped) and summarise only the newer messages.
                body = (content.split("\n", 1)[1].strip()
                        if "\n" in content else content.strip())
                if body:
                    carried_summary = body
                continue
            if not isinstance(content, str):
                # Flatten tool/structured content for the summariser
                try:
                    content = " ".join(
                        str(part.get("text", part))
                        for part in content
                        if isinstance(part, (dict, str))
                    )
                except Exception:
                    content = str(content)
            if not content:
                continue
            # Cap each message individually so a single huge tool_result
            # can't dominate the summariser's input.
            if len(content) > 4000:
                content = content[:2000] + "\n... [truncated for summary] ...\n" + content[-1500:]
            rendered.append(f"[{role}]\n{content}")
        if not rendered:
            # Nothing new to summarise beyond a prior recap — carry it forward.
            return carried_summary

        system_prompt = (
            "You are a conversation-summarisation assistant. Produce a "
            "compact but faithful recap of the dialog below that the next "
            "turn can use as drop-in context. Required structure:\n"
            "\n"
            "## Goal\n"
            "One sentence — what the user was working toward.\n\n"
            "## Key facts & decisions\n"
            "- 3-7 bullets covering decisions made, files changed, "
            "constraints discovered, errors hit + fixes, current state.\n"
            "- Cite file paths + line numbers when they appeared.\n\n"
            "## Open items / next step\n"
            "- 1-3 bullets — what is still TODO or unresolved.\n\n"
            "Be terse. Skip pleasantries, repeated greetings, and dead-end "
            "explorations the user already moved past. Preserve specific "
            "numbers, identifiers, error messages, and exact filenames."
        )
        joined = "\n\n".join(rendered)
        if len(joined) > 60_000:
            joined = joined[:30_000] + "\n... [middle truncated] ...\n" + joined[-25_000:]
        user_msg = (
            "Summarise the conversation below according to the structure "
            "in the system prompt.\n\n=== CONVERSATION START ===\n"
            f"{joined}\n=== CONVERSATION END ==="
        )

        try:
            text_parts: list[str] = []
            stream = self.client.stream_message(
                messages=[{"role": "user", "content": user_msg}],
                system=system_prompt,
                max_tokens=1500,
            )
            for event in stream:
                if getattr(event, "type", "") == "text_delta":
                    chunk = getattr(event, "text", "") or ""
                    if chunk:
                        text_parts.append(chunk)
                # Stop once we have a reasonable summary length so a
                # runaway model can't waste tokens.
                if sum(len(p) for p in text_parts) > 5000:
                    break
        except Exception:
            return ""
        summary = "".join(text_parts).strip()
        # Keep the earlier recap AHEAD of the new one so the oldest context
        # isn't re-compressed (and eroded) on every compaction.
        if carried_summary:
            summary = (carried_summary + "\n\n" + summary).strip() if summary \
                else carried_summary
        return summary

    def _deterministic_digest(self, msgs: list[dict]) -> str:
        """Last-resort recap when both summarisers came back empty.

        No model call, no settings, no failure mode: the role and the
        opening of each message that is about to be replaced. It is a poor
        summary, which is why it is last; it is still enormously better
        than the silent ``return`` that used to sit here and left the
        session over budget with no record that anything had been tried.
        """
        lines: list[str] = []
        for msg in msgs:
            if not isinstance(msg, dict):
                continue
            content = msg.get("content", "")
            if not isinstance(content, str):
                try:
                    content = " ".join(
                        str(p.get("text", p)) if isinstance(p, dict) else str(p)
                        for p in content
                    )
                except Exception:
                    content = str(content)
            content = " ".join(content.split())
            if not content:
                continue
            role = str(msg.get("role", "?"))
            lines.append(f"[{role}] {content[:300]}")
        if not lines:
            return ""
        return (
            f"No summary could be produced for these {len(lines)} "
            f"message(s); their openings follow verbatim.\n"
            + "\n".join(lines)
        )

    def _compact_history(self, *, force: bool = False) -> None:
        """Summarize older messages, keeping recent ones intact.

        Three-stage strategy for multi-day session resilience:
        1. ``_shorten_oldest_non_goal_messages`` — gentle in-place trim of
           the oldest long non-goal content when usage crosses 70% of the
           window. Real user goals survive verbatim; machine turns and
           assistant payloads get head + tail with a middle marker.
        2. Full ``_compact_history`` summarisation when usage crosses
           ``auto_compact_pct`` (default 95%).
        3. LLM-quality summary preferred (API backends), extractive
           fallback (CLI backend, summary failure, opt-out via
           ``agent.llm_compaction: false``).

        ``force`` runs the full summarisation regardless of pressure — the
        manual ``/compact`` command, where the user has already decided.
        The structural floor and every protection still apply.

        Pinned messages (``msg["_pinned"]``) are exempt from every stage:
        the trims skip them and the summary path carries them verbatim
        into the kept set. Trimmed/cleared bodies are first persisted to
        the session's elided store so ``history_get('elided:<ref>')`` can
        recover them.
        """
        # Stage 1: sliding-window in-place trim. Fires every turn at 70%
        # so the cliff at 95% is much rarer. Cheap — no LLM call.
        try:
            if self._should_slide():
                _trimmed = self._shorten_oldest_non_goal_messages()
                if _trimmed:
                    # Record the trim event into last_compaction_info so
                    # /context shows what happened.
                    self.last_compaction_info = {
                        "kind": "sliding_window",
                        "messages_trimmed": _trimmed,
                        "tokens_after": self._estimate_context_tokens(),
                        "archived_at": _time.time(),
                    }
        except Exception:
            pass

        # Full compaction is driven by CONTEXT PRESSURE (token budget), never
        # by raw message count alone: a dozen short messages can sit at ~15%
        # of the window, and summarising there throws away the live working
        # context — forcing the agent to re-discover its own work (Jerome
        # 2026-06-13: "compacted at 15% full → agent confused"). This was the
        # legacy ``msg_count >= 12 OR tokens > 95%`` trigger; the OR meant the
        # 0.95 token threshold (raised precisely to avoid early compaction)
        # was undermined by the message count.
        #
        # Pressure is now asked FIRST. The count gate used to return before
        # anything looked at the budget, and since it can only have an effect
        # once pressure is already over the line, its whole effect was to
        # block compaction exactly when it was needed: eleven messages at 99%
        # of the window compacted nothing and the request went out over the
        # window. What remains of the count is the structural floor it was
        # described as — there must be something between the summary block
        # and the kept tail, or there is nothing to summarise.
        if not force and not self._should_auto_compact():
            return
        if len(self.messages) <= self._KEEP_RECENT:
            return

        old_msgs = self.messages[:-self._KEEP_RECENT]
        recent = self.messages[-self._KEEP_RECENT:]
        # Pinned messages are carried into the kept set VERBATIM — they are
        # excluded from summarisation entirely (summarising them would both
        # duplicate and endanger their exact wording).
        pinned_old = [
            m for m in old_msgs if isinstance(m, dict) and m.get("_pinned")
        ]
        compactable = [
            m for m in old_msgs
            if not (isinstance(m, dict) and m.get("_pinned"))
        ]

        # Stage 2a: deterministic context-editing BEFORE any LLM summary.
        # Structure-preserving hard-clear of old bulky tool output — free,
        # deterministic, lossless for user goals + assistant reasoning. If it
        # brings us back under the auto-compact budget, skip the lossy/paid
        # summary and keep the (edited) history intact.
        try:
            _cleared = self._stub_oldest_non_goal_messages(old_msgs)
            if _cleared and not force and not self._should_auto_compact():
                self.last_compaction_info = {
                    "kind": "context_edit",
                    "tool_results_cleared": _cleared,
                    "tokens_after": self._estimate_context_tokens(),
                    "archived_at": _time.time(),
                }
                return
        except Exception:
            pass

        # 1) Try LLM summarisation (only when hard-clear wasn't enough)
        summary = ""
        try:
            from delfin.user_settings import load_settings as _load_settings
            _agent_cfg = (_load_settings() or {}).get("agent", {}) or {}
            use_llm = bool(_agent_cfg.get("llm_compaction", True))
        except Exception:
            use_llm = True
        if use_llm:
            summary = self._llm_summarize_old_messages(compactable)

        # 2) Fall back to extractive summarisation if LLM unavailable
        if not summary:
            summary_parts: list[str] = []
            # Selective extractive: user messages (= GOALS) survive in
            # near-full form because they're cheap and load-bearing.
            # Assistant text gets summarised. Tool-result-shaped content
            # gets aggressively dropped — the answers it produced are
            # usually already mentioned in the assistant text we keep.
            for msg in compactable:
                role = msg.get("role", "")
                content = msg.get("content", "")
                if self._is_machine_turn(msg):
                    # This history's tool output: a synthetic user turn.
                    # The branch that used to stand here tested
                    # ``role == "tool"``, a role no append site ever
                    # produces, so it selected nothing while the bulk it
                    # was written for sailed through the "user = GOAL"
                    # branch at full 400-char weight.
                    text = content if isinstance(content, str) else ""
                    if text.strip().startswith('{"error"'):
                        summary_parts.append(f"Tool (error): {text[:120]}")
                    elif "error" in text[:400].lower():
                        summary_parts.append(f"Tool (error): {text[:200]}")
                    continue
                if role == "user":
                    text = content if isinstance(content, str) else ""
                    # A prior compaction's summary block is ALREADY a
                    # compressed recap of many messages. On a second
                    # compaction it lands here as a "user" message — and
                    # truncating it to 400 chars (as if it were a one-line
                    # goal) silently drops everything deeper than ~400 chars,
                    # compounding loss across repeated compactions on long
                    # sessions. Carry it forward near-whole instead (bounded,
                    # header stripped so we don't nest "[summary]" markers).
                    if text.lstrip().startswith(self._SUMMARY_BLOCK_PREFIX):
                        body = text.split("\n", 1)[1].strip() if "\n" in text else text.strip()
                        if len(body) > 3000:
                            body = body[:3000] + "\n... [older summary detail elided] ..."
                        if body:
                            summary_parts.append(body)
                        continue
                    # Goals deserve more than one line — keep up to
                    # 400 chars of the first sentence/paragraph so the
                    # post-compact agent can read the original intent.
                    keep = text.strip()[:400]
                    if keep:
                        summary_parts.append(f"User: {keep}")
                elif role == "assistant":
                    extracted = self._summarize_output_for_handoff(content, limit=300)
                    if extracted:
                        summary_parts.append(f"Assistant: {extracted}")
            summary = "\n".join(summary_parts)
        summary_kind = "summary"
        if not summary:
            # Both summarisers came back empty. Returning here left the
            # session over budget with NO record that compaction had been
            # attempted at all, so the next turn repeated it and /context
            # still reported the compaction before last. A deterministic
            # digest always has something to say as long as there is
            # anything to compact, and the record names which path ran.
            summary = self._deterministic_digest(compactable)
            summary_kind = "deterministic_digest"
        if not summary:
            self.last_compaction_info = {
                "kind": "summary_unavailable",
                "messages_compacted": 0,
                "tokens_saved": 0,
                "archived_at": _time.time(),
                "note": (
                    f"nothing could be summarised out of "
                    f"{len(compactable)} compactable message(s); the "
                    f"history is unchanged and the context is still over "
                    f"budget."
                ),
            }
            return

        n_compacted = len(compactable)
        tokens_before = self._estimate_context_tokens()

        # Archive the pre-compaction transcript so the user can scroll back
        # to it later via /session archive ls — nothing is ever lost in a
        # long session, only summarised in-place for the next turn. Pinned
        # messages stay in the live list, so only the summarised ones are
        # archived here.
        try:
            from delfin.agent.session_store import archive_pre_compaction_transcript
            archive_pre_compaction_transcript(
                getattr(self, "session_id", "") or "",
                compactable,
                info={
                    "messages_compacted": n_compacted,
                    "tokens_before": tokens_before,
                },
            )
        except Exception:
            pass

        self.messages = [
            {"role": "user", "content": f"[Conversation summary — older messages compacted]\n{summary}"},
            {"role": "assistant", "content": "Understood. I have the context from our earlier conversation."},
        ] + pinned_old + recent
        # `recent` starts with an assistant turn whenever compaction fires right
        # after a user message was appended, so the assistant_ack above lands
        # next to it → two consecutive assistant messages. Strict chat templates
        # (vLLM/KIT/Ollama) reject a non-alternating list; sanitize BEFORE the
        # request goes out rather than relying on the next turn's cleanup.
        self._sanitize_messages()
        # The context just shrank a lot; drop the stale pre-compaction input
        # floor so the estimate reflects the compacted size immediately (the
        # next request's message_start repopulates the real count).
        self._last_input_tokens = 0
        self._trimmed_chars_since_floor = 0

        tokens_after = self._estimate_context_tokens()
        self.last_compaction_info = {
            "kind": summary_kind,
            "messages_compacted": n_compacted,
            "pinned_kept": len(pinned_old),
            "tokens_before": tokens_before,
            "tokens_after": tokens_after,
            "tokens_saved": max(0, tokens_before - tokens_after),
            "archived_at": _time.time(),
            "forced": bool(force),
        }
        if summary_kind == "deterministic_digest":
            self.last_compaction_info["note"] = (
                "the summariser produced nothing, so the older messages "
                "were replaced by a deterministic digest of their openings "
                "— less faithful than a summary. The full transcript is in "
                "the archive."
            )

        # CLI backend: by default tear down the persistent process so the
        # next stream_message starts fresh and the compacted history
        # actually reaches the API. Opt-out via the user setting
        # ``agent.compact_resets_cli`` keeps the CLI alive and treats
        # compaction as local-RAM-only (cheaper restart, but the CLI's
        # internal session memory still holds the full untruncated
        # history, so API-level cost reduction won't materialise).
        if hasattr(self.client, "kill") and self.backend == "cli":
            try:
                from delfin.user_settings import load_settings as _load_settings
                _agent_cfg = (_load_settings() or {}).get("agent", {}) or {}
                _resets = bool(_agent_cfg.get("compact_resets_cli", True))
            except Exception:
                _resets = True
            if _resets:
                self.client.kill()

    # -- Context usage tracking (Feature 4) --------------------------------

    def _record_context_usage(self, response_text: str) -> None:
        """Record which injected sections the agent actually referenced."""
        if self._context_tracker is None:
            return
        injected = getattr(self.loader, "_last_injected_sections", [])
        if not injected:
            return
        try:
            self._context_tracker.record_usage(
                sections_injected=injected,
                response_text=response_text,
                role_id=self.current_role,
                provider=self.provider,
            )
        except Exception:
            pass  # never let tracking break the main flow

    # -- Progressive disclosure: NEED_CONTEXT handling (Feature 1) ---------

    _NEED_CONTEXT_RE = re.compile(r"NEED_CONTEXT:\s*(\w+)", re.IGNORECASE)
    _MAX_NEED_CONTEXT_ROUNDS = 3

    def _handle_need_context(
        self,
        response_text: str,
        original_user_message: str,
        memory_context: str = "",
        **stream_kwargs: Any,
    ) -> str:
        """Detect NEED_CONTEXT requests and inject the requested sections.

        Returns the final response text (may include multiple rounds).
        """
        for _ in range(self._MAX_NEED_CONTEXT_ROUNDS):
            matches = self._NEED_CONTEXT_RE.findall(response_text)
            if not matches:
                break

            # Load requested sections
            injected_parts: list[str] = []
            for section_name in matches:
                section_name = section_name.lower().strip()
                if section_name not in AVAILABLE_ON_DEMAND:
                    continue
                content = self._load_on_demand_section(
                    section_name, original_user_message, memory_context,
                )
                if content:
                    injected_parts.append(
                        f"[System] Requested context '{section_name}':\n{content}"
                    )

            if not injected_parts:
                break

            # Inject as a user message and continue the conversation
            context_msg = "\n\n".join(injected_parts)
            self.messages.append({"role": "user", "content": context_msg})

            # Stream continuation
            chunks: list[str] = []
            try:
                system_prompt = self._build_current_system_prompt(
                    memory_context, task_text=original_user_message,
                )
                for event in self.client.stream_message(
                    system=system_prompt,
                    messages=self._wire_messages(),
                    max_tokens=stream_kwargs.get("max_tokens") or self.max_tokens_for_role(self.current_role),
                    session_id=self.session_id,
                    thinking_budget=stream_kwargs.get("thinking_budget", 0),
                ):
                    if self._stop_requested:
                        break
                    if event.type == "text_delta" and event.text:
                        chunks.append(event.text)
                        on_token = stream_kwargs.get("on_token")
                        if on_token:
                            on_token(event.text)
                    elif event.type == "message_start":
                        with self._lock:
                            self.token_usage["input"] += event.input_tokens
                            self.token_usage["cached"] = self.token_usage.get(
                                "cached", 0) + (getattr(event, "cached_tokens", 0) or 0)
                    elif event.type == "message_delta":
                        with self._lock:
                            self.token_usage["output"] += event.output_tokens
                            self.cost_usd += event.cost_usd
                            self.token_usage["cached"] = self.token_usage.get(
                                "cached", 0) + (getattr(event, "cached_tokens", 0) or 0)
            except Exception:
                break

            continuation = "".join(chunks)
            if continuation:
                self.messages.append({"role": "assistant", "content": continuation})
                response_text = continuation
            else:
                break

        return response_text

    def _load_on_demand_section(
        self,
        section_name: str,
        task_text: str,
        memory_context: str = "",
    ) -> str:
        """Load a specific context section on demand."""
        if section_name == "playbook":
            return self.loader._load_relevant_playbook_context(task_text)
        elif section_name == "repo_map":
            return self.loader._load_repo_map_context(task_text)
        elif section_name == "profile":
            return self.loader._load_profile_context(self.mode)
        elif section_name == "memory":
            return memory_context
        elif section_name == "prior_outputs":
            if not self.role_outputs:
                return ""
            parts = []
            for rid, output in self.role_outputs.items():
                parts.append(f"### {rid}\n{output[:3000]}")
            return "\n\n".join(parts)
        return ""

    # -- Auto-verification after code edits (Feature 2) --------------------

    _AUTO_VERIFY_ROLES = {"solo_agent", "builder_agent"}
    _CODE_EDIT_TOOLS = {"Edit", "Write"}

    def check_auto_verify(self, tool_name: str, tool_output: str) -> str | None:
        """Check if auto-verification should run after a tool call.

        Returns a verification message to inject, or None.
        Called by the dashboard worker after each tool result.
        """
        role = self.current_role
        if role not in self._AUTO_VERIFY_ROLES:
            return None
        if tool_name not in self._CODE_EDIT_TOOLS:
            return None
        # Check if the edited file is a Python file
        # (tool_output from Edit/Write contains the file path)
        if ".py" not in tool_output and ".py" not in str(getattr(self, "_last_edit_path", "")):
            return None
        return (
            "[System] A Python file was just modified. "
            "Run `python -m pytest tests/ -x -q --tb=short` to verify. "
            "If tests fail, fix the issue before continuing."
        )

    # -- Within-session error learning (Feature 4) -------------------------

    def __init_error_memory(self) -> None:
        """Initialize error memory if not present."""
        if not hasattr(self, "_session_errors"):
            self._session_errors: list[dict[str, str]] = []

    def record_session_error(
        self,
        tool_name: str,
        error_text: str,
    ) -> None:
        """Record a tool error within the current session."""
        self.__init_error_memory()
        # Deduplicate
        key = f"{tool_name}:{error_text[:100]}"
        for existing in self._session_errors:
            if existing.get("key") == key:
                existing["count"] = str(int(existing.get("count", "1")) + 1)
                return
        self._session_errors.append({
            "key": key,
            "tool": tool_name,
            "error": error_text[:300],
            "count": "1",
        })
        # Keep bounded
        if len(self._session_errors) > 20:
            self._session_errors = self._session_errors[-20:]

    def format_error_context(self) -> str:
        """Format session errors as context for the system prompt."""
        self.__init_error_memory()
        if not self._session_errors:
            return ""
        lines = ["[Session errors — do NOT repeat these:]"]
        for err in self._session_errors[-5:]:
            count = err.get("count", "1")
            lines.append(f"- {err['tool']}: {err['error'][:150]} (x{count})")
        return "\n".join(lines)

    def advance_role(self) -> bool:
        """Advance to the next role in the route.

        Returns True if there is a next role, False if the cycle is complete.
        """
        role = self.current_role
        for msg in reversed(self.messages):
            if msg["role"] == "assistant":
                self.role_outputs[role] = msg["content"]
                self.compaction_summaries[role] = self._summarize_output_for_handoff(
                    msg["content"]
                )
                break

        self.current_role_index += 1
        return not self.is_cycle_complete

    # -- Parallel agent execution (Feature 5) ------------------------------

    # Roles that can run in parallel (same pipeline position).
    _PARALLEL_GROUPS: list[set[str]] = [
        {"critic_agent", "runtime_agent"},  # both review the plan independently
    ]

    def can_run_parallel(self) -> list[str] | None:
        """Check if the current + next role can run in parallel.

        Returns a list of role IDs that can run concurrently, or None.
        """
        if self.current_role_index >= len(self.route) - 1:
            return None
        current = self.route[self.current_role_index]
        next_role = self.route[self.current_role_index + 1]
        for group in self._PARALLEL_GROUPS:
            if current in group and next_role in group:
                return [current, next_role]
        return None

    def advance_past_parallel(self, outputs: dict[str, str]) -> bool:
        """Skip past parallel roles, storing their outputs.

        Call this after running parallel roles instead of advance_role().
        ``outputs`` maps role_id → output text.
        """
        for rid, output in outputs.items():
            self.role_outputs[rid] = output
            self.compaction_summaries[rid] = self._summarize_output_for_handoff(output)
        # Advance past all parallel roles
        self.current_role_index += len(outputs)
        return not self.is_cycle_complete

    @staticmethod
    def _summarize_output_for_handoff(text: str, limit: int = 1200) -> str:
        """Compress a role output keeping structured content (headers, lists, decisions).

        Default 1200 chars keeps enough detail for downstream agents to act
        without re-reading the full output.  The ``build_handoff_message``
        method selects between summary and full text per-role.
        """
        if not text:
            return ""
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        selected: list[str] = []
        total = 0
        for line in lines:
            lowered = line.lower()
            is_structural = (
                lowered.startswith("**status:")
                or lowered.startswith("**question:")
                or lowered.startswith("**summary:")
                or lowered.startswith("**goal:")
                or lowered.startswith("**plan:")
                or lowered.startswith("## ")
                or lowered.startswith("### ")
                or lowered.startswith("- ")
                or re.match(r"^\d+\.\s", line)
            )
            if is_structural:
                selected.append(line)
                total += len(line)
            if total >= limit:
                break
        # If nothing structural found, take the first lines
        if not selected:
            for line in lines[:8]:
                selected.append(line)
                total += len(line)
                if total >= limit:
                    break
        summary = "\n".join(selected)
        if len(summary) > limit:
            summary = summary[: limit - 3].rstrip() + "..."
        return summary

    def _capture_current_role_summary(self) -> None:
        """Persist the latest assistant output before message compaction."""
        role = self.current_role
        if not role:
            return
        for msg in reversed(self.messages):
            if msg.get("role") == "assistant" and msg.get("content"):
                content = msg["content"]
                self.role_outputs.setdefault(role, content)
                self.compaction_summaries.setdefault(
                    role,
                    self._summarize_output_for_handoff(content),
                )
                return

    def request_stop(self) -> None:
        """Request the current streaming response to stop.

        The stop is stamped with the turn it is aimed at — the one holding
        the gate, or the last one to have held it. Without that stamp the
        flag was anonymous engine state, and anonymous state gets cleared
        by whoever comes next (see ``clear_stop``).
        """
        self._stop_requested = True
        self._stop_owner_turn = getattr(self, "_turn_id", 0)

    def clear_stop(self) -> None:
        """Cancel a previous stop so a NEW turn can run.

        The flag used to be cleared in exactly one place -- the top of
        ``stream_response`` -- which works only if that method is reached.
        A caller that checks the flag BEFORE calling it can never clear
        it again: the code that resets the flag sits behind the check the
        flag itself fails. That is not a hypothetical; it silenced the
        dashboard permanently after any Stop, and every re-send reported
        "no output" without ever contacting the backend.

        A stop belongs to the turn it interrupted. Starting a new turn is
        the moment it stops applying, so the owner of that decision is
        whoever begins the turn -- and ONLY once that turn has actually
        begun. Clearing a stop whose turn is still in flight is not
        starting fresh, it is disarming a running turn's only brake: the
        stopped turn polls this flag between rounds, so a caller that
        cleared it here made the abandoned plan resume and bill for
        itself. Refused, and the caller's own turn is refused by the gate
        a moment later anyway.
        """
        if (getattr(self, "_stop_requested", False)
                and getattr(self, "_turn_in_flight", False)
                and getattr(self, "_stop_owner_turn", None)
                == getattr(self, "_turn_id", 0)):
            return
        self._stop_requested = False
        self._stop_owner_turn = None

    def _hand_back_undelivered_steer(self) -> str:
        """Return what the user typed mid-run and never got delivered.

        The client's tool loop has its own careful version of this in its
        stop branch. TRACED live, 2026-08-14, and it never runs:

            0.42s  turn starts
            4.16s  _stop_was_requested() -> False      (round 1, top)
           13.13s  TOOL DISPATCHED + stop requested
           13.14s  turn returned

        Ten milliseconds. The engine checks the stop flag BETWEEN STREAM
        EVENTS and breaks out there, and the client's loop is a generator —
        a consumer that stops iterating simply abandons it. Its stop
        branch, its drains and its hand-back never execute. The engine
        wins that race in the ordinary case, not an edge one, so the
        client-side handling is unreachable whenever a stop lands during a
        tool round.

        Measured consequence, three runs out of three: the typed message
        stayed in the client's queue. It would then be injected into some
        later, unrelated turn — out of order with whatever the user typed
        AFTER pressing Stop, which is exactly the harm the client's branch
        was written to prevent.

        Whoever ends the turn owns the hand-back. Run notes are left where
        they are on purpose: they are facts about the session that are
        still true, the client object outlives this turn, and the next
        turn drains them at its first round.
        """
        try:
            pending = self.client._drain_steer()
        except Exception:
            return ""
        # Only a list of real strings is handed back. A client whose drain
        # returns something else -- an older backend, a stub, a mock -- must
        # not get its repr rendered into the user's answer; that is how a
        # MagicMock ended up in the transcript once already.
        if not isinstance(pending, (list, tuple)):
            return ""
        pending = [p for p in pending if isinstance(p, str) and p.strip()]
        if not pending:
            return ""
        quoted = "; ".join(f"“{p}”" for p in pending[:3])
        more = " …" if len(pending) > 3 else ""
        return (" Not delivered, because the turn stopped first: "
                + quoted + more + " — send it again if it still applies.")

    def steer(self, text: str) -> bool:
        """Inject a user message into the RUNNING tool loop (mid-loop steering).
        The model reacts to it on the next round, without ending the turn.
        Returns True if the active client supports it (API/KIT/Ollama engine)."""
        client = getattr(self, "client", None)
        if client is not None and hasattr(client, "push_steer"):
            try:
                client.push_steer(text)
                return True
            except Exception:
                return False
        return False

    def set_live_state(self, text: str) -> None:
        """Update the per-turn live-state snippet.

        The string is appended to the system prompt (cache-friendly,
        regenerated only when needed) instead of being baked into the
        user message body — so old turns in ``self.messages`` don't
        carry stale dashboard state forever.
        """
        self._live_state = text or ""

    def _fresh_session_id(self) -> str:
        """A new session id for backends that supply none.

        A backend that announces its own id on the stream returns "" here
        and is filled in when the event arrives; everyone else gets a
        minted UUID, so each session is a distinct, stable bucket for task
        scoping and session save/load.

        Asked of the CLIENT, not of the backend string. ``backend="cli"``
        with ``provider="kit"`` builds an OpenAIClient — create_client
        routes on the provider — and that client announces nothing. Keying
        on the string left the id empty forever, and an empty task session
        id means UNSCOPED: a fresh session's first turn reported every
        earlier session's open tasks as its own unfinished work.
        """
        # `is True`, not truthiness: an arbitrary stand-in answers every
        # attribute lookup with something truthy, and the failure would be
        # silent — an empty id is not an error, it is the unscoped bucket.
        # Anything that does not declare the capability exactly gets an id,
        # which is the safe direction.
        if getattr(self.client, "supplies_session_id", False) is True:
            return ""
        import uuid as _u
        return _u.uuid4().hex

    def _sync_task_session(self) -> None:
        """Bind the agent's task list to the current session id so task_create/
        task_list scope to THIS session (empty id → task_list shows all
        workspace tasks). Best-effort; no-op when there are no KIT permissions."""
        try:
            kp = self.kit_permissions
            if kp is not None:
                kp.task_session_id = self.session_id or ""
        except Exception:
            pass

    def _run_budget(self) -> tuple[float, float]:
        """(usd, seconds) run-scoped budget for this SESSION — the cumulative
        ceiling that makes unattended runs deployable (the per-turn cap below
        only stops a single runaway turn; turns compose into unbounded run
        cost without this). 0 disables either dimension. An explicit
        ``self.run_budget_usd``/``self.run_budget_s`` attribute (set e.g. by
        the scheduler daemon per entry) overrides the settings values."""
        usd = float(getattr(self, "run_budget_usd", 0.0) or 0.0)
        secs = float(getattr(self, "run_budget_s", 0.0) or 0.0)
        if usd <= 0 or secs <= 0:
            try:
                from delfin.user_settings import load_settings as _ls
                ag = (_ls() or {}).get("agent", {}) or {}
                if usd <= 0:
                    usd = max(0.0, float(ag.get("run_budget_usd", 0.0) or 0.0))
                if secs <= 0:
                    secs = max(0.0, float(ag.get("run_budget_s", 0.0) or 0.0))
            except Exception:
                pass
        return usd, secs

    def _run_budget_turns(self) -> int:
        """Turn ceiling for this run — 0 disables it.

        The third dimension exists because the first one cannot always be
        measured. A turn is countable whatever the model costs, so this
        is the ceiling that still holds when the USD one is inert (see
        :meth:`_usd_budget_enforced`). Same precedence as the other two:
        an explicit attribute set per scheduler entry beats the setting.
        """
        turns = int(getattr(self, "run_budget_turns", 0) or 0)
        if turns <= 0:
            try:
                from delfin.user_settings import load_settings as _ls
                ag = (_ls() or {}).get("agent", {}) or {}
                turns = int(ag.get("run_budget_turns", 0) or 0)
            except Exception:
                turns = 0
        return max(0, turns)

    def _turns_taken(self) -> int:
        """Turns this run has actually run, priced or not.

        Deliberately the SUM of the three price states and not the
        unpriced count: a run that starts on an unpriced model and
        switches to a priced one would otherwise stop spending its
        ceiling half way through.
        """
        return (int(getattr(self, "_measured_cost_turns", 0) or 0)
                + int(getattr(self, "_non_billing_turns", 0) or 0)
                + int(getattr(self, "_unpriced_turns", 0) or 0))

    def _client_retry_count(self) -> int | None:
        """How often the client re-issued this turn's request, if it says.

        ``None`` means the backend does not report retries — which is not
        a count of zero. A client that retries internally (transient 5xx /
        timeout, 1.5 + 3 + 6 s of sleep) exposes ``last_turn_retries``;
        without it a turn that fought for twenty seconds and a quiet one
        record the same thing.
        """
        try:
            raw = getattr(self.client, "last_turn_retries", None)
            return int(raw) if raw is not None else None
        except (TypeError, ValueError):
            return None

    def _note_turn_price_state(self, cost_delta: float = 0.0) -> None:
        """Count the turn that just ran by what its cost could mean.

        Three ways the number means something, resolved in the order the
        dashboard's own cost line already uses:

        * NON_BILLING — the provider charges no USD for the call at all
          (KIT quota, local hardware). A measured zero.
        * The backend measured it itself: a reported spend above zero, or
          the CLI backend, where the provider states the turn's cost
          (often zero, on a subscription) instead of leaving it to a rate
          table. An observation outranks the table either way.
        * PRICED — a published rate applied to counted tokens.

        Everything else is UNKNOWN: no rate, nothing observed, and a 0.0
        in ``self.cost_usd`` that every USD gate reads as thrift.
        """
        from . import pricing as _pricing
        state = _pricing.resolve(
            str(getattr(self.client, "model", "") or ""),
            str(getattr(self, "provider", "") or ""),
        ).state
        if state == _pricing.NON_BILLING:
            self._non_billing_turns += 1
        elif (state == _pricing.PRICED
                or float(cost_delta or 0.0) > 0
                or getattr(self, "backend", "") == "cli"):
            self._measured_cost_turns += 1
        else:
            self._unpriced_turns += 1

    def _usd_budget_enforced(self) -> bool:
        """Whether the configured USD ceiling actually bounds this run.

        False in one case that matters: a USD budget is set and at least
        one turn ran on a model with no published rate. Such a turn adds
        0.0 to ``self.cost_usd``, so the fraction the gates compare is a
        sum over the turns that happened to be measurable — it can sit at
        0% through a run of any size, and the refusal above 110% can never
        fire. The tier aliases the model picker ships (``opus``,
        ``sonnet``, ``haiku``) are precisely the unpriced ids, so this is
        the default case, not an exotic one.

        A quota-funded or locally served run is NOT this case: its zero is
        a measurement, its ceiling holds, and calling it unmeasured would
        be a false statement in the other direction.
        """
        usd, _secs = self._run_budget()
        return usd > 0 and int(getattr(self, "_unpriced_turns", 0) or 0) == 0

    def _unmeasured_budget_block(self) -> str:
        """Say ONCE that no monetary ceiling is in force, and what is.

        Silent unless a USD budget is configured AND an unpriced turn has
        run: with no budget there is nothing to mislead anyone about, and
        with every turn priced the ordinary wind-down block does the job.
        Turn counts and token totals are always measured, so they are what
        this offers in place of the ceiling that cannot be enforced.
        """
        usd, secs = self._run_budget()
        unpriced = int(getattr(self, "_unpriced_turns", 0) or 0)
        if (usd <= 0 or unpriced == 0
                or getattr(self, "_unmeasured_budget_notice_shown", False)):
            return ""
        self._unmeasured_budget_notice_shown = True
        turn_cap = self._run_budget_turns()
        # The prompt reaches the model and nobody else. A run budget is
        # what makes an unattended run deployable, so the one fact that
        # it is not in force has to reach whoever is not watching.
        self._emit_budget_attention(
            "unenforceable",
            "Run budget: the USD ceiling is not in force",
            f"{unpriced} turn(s) ran on a model with no published rate, so "
            f"the ${usd:.2f} ceiling cannot be enforced and neither can the "
            f"per-turn cost breaker."
            + (f" The {turn_cap}-turn ceiling IS in force."
               if turn_cap > 0 else
               " Set agent.run_budget_turns or agent.run_budget_s for a "
               "ceiling that can be."))
        turns = (unpriced + int(getattr(self, "_measured_cost_turns", 0) or 0)
                 + int(getattr(self, "_non_billing_turns", 0) or 0))
        tok_in = int(self.token_usage.get("input", 0) or 0)
        tok_out = int(self.token_usage.get("output", 0) or 0)
        lines = [
            "# Run budget NOT enforceable in USD (auto-injected, once)",
            f"- {unpriced} of {turns} turn(s) so far ran on a model with no "
            f"published rate, so their cost was never measured. The "
            f"${usd:.2f} ceiling is not in force, and neither is the "
            f"per-turn cost breaker.",
            f"- Measured instead: {turns} turn(s), {tok_in:,} input + "
            f"{tok_out:,} output tokens. Judge the size of this run by "
            f"those, not by a dollar figure.",
        ]
        if turn_cap > 0:
            lines.append(f"- The {turn_cap} turn(s) ceiling IS measured and "
                         "still ends this run.")
        if secs > 0:
            lines.append(f"- The {secs:.0f}s wall-clock ceiling IS measured "
                         "and still ends this run.")
        if turn_cap <= 0 and secs <= 0:
            lines.append("- For a ceiling that can be enforced, set "
                         "agent.run_budget_turns or agent.run_budget_s, or "
                         "pick a model id with a published rate.")
        return "\n".join(lines)

    def _run_budget_status(self) -> tuple[float, bool]:
        """(worst fraction spent, exhausted). Fraction is max over the
        enabled dimensions; (0.0, False) when no budget is configured.

        The USD term is a sum over measurable turns only — see
        :meth:`_usd_budget_enforced` for what it cannot see.
        """
        frac, exhausted, _dim = self._run_budget_detail()
        return frac, exhausted

    def _run_budget_detail(self) -> tuple[float, bool, str]:
        """As :meth:`_run_budget_status`, plus the setting that is worst.

        A refusal that names ``run_budget_usd`` when the wall clock or the
        turn count stopped the run sends its reader to the wrong dial.
        """
        usd, secs = self._run_budget()
        turns_cap = self._run_budget_turns()
        frac, dim = 0.0, ""
        if usd > 0:
            share = float(self.cost_usd) / usd
            if share > frac:
                frac, dim = share, "agent.run_budget_usd"
        if secs > 0:
            import time as _t
            started = float(getattr(self, "_run_started_at", 0.0) or 0.0)
            if started > 0:
                share = (_t.time() - started) / secs
                if share > frac:
                    frac, dim = share, "agent.run_budget_s"
        if turns_cap > 0:
            share = float(self._turns_taken()) / turns_cap
            if share > frac:
                frac, dim = share, "agent.run_budget_turns"
        return frac, frac >= 1.0, dim

    def _emit_budget_attention(self, level: str, title: str,
                               detail: str) -> None:
        """Raise ONE durable budget event per level for this run.

        ``budget_warning`` was declared, labelled and ranked in the
        attention inbox with nothing anywhere that emitted it. Everything
        the budget machinery says goes either into the system prompt (the
        model reads it, nobody else) or into the chat stream (gone with
        the scrollback), so an unattended run that stopped on its ceiling
        overnight showed a silent desktop, an empty inbox and a green
        doctor. Best-effort: a notice must never break the turn it
        describes.
        """
        try:
            if level in self._budget_attention_levels:
                return
            self._budget_attention_levels.add(level)
            from delfin.agent.attention import emit_attention
            emit_attention(
                "budget_warning",
                session_id=str(getattr(self, "session_id", "") or ""),
                title=title[:200],
                detail=detail[:400],
                workspace=str(getattr(self, "repo_dir", "") or ""),
            )
        except Exception:
            pass

    def _build_budget_block(self) -> str:
        """Per-turn budget status for the prompt: silent when unconfigured
        or comfortably below the wind-down threshold; from 80% it instructs
        the agent to wrap up gracefully instead of being cut off mid-work.

        The unmeasured-cost notice comes first and independently of the
        threshold, because the threshold is exactly what an unpriced run
        can never reach: its fraction stays at 0% however long it runs.
        """
        parts: list[str] = []
        try:
            parts.append(self._unmeasured_budget_block())
        except Exception:
            pass
        try:
            frac, exhausted = self._run_budget_status()
        except Exception:
            frac, exhausted = 0.0, False
        if frac >= 0.8:
            pct = min(999, int(frac * 100))
            usd, secs = self._run_budget()
            turn_cap = self._run_budget_turns()
            ceiling = ", ".join(
                [p for p in (f"${usd:.2f}" if usd > 0 else "",
                             f"{secs:.0f}s" if secs > 0 else "",
                             f"{turn_cap} turns" if turn_cap > 0 else "")
                 if p])
            if exhausted:
                parts.append(
                    "# Run budget EXHAUSTED\n"
                    f"- {pct}% of this run's budget is spent. Finish your "
                    "CURRENT sentence of work only: save state, mark task "
                    "statuses honestly, and summarise what remains. No new "
                    "work."
                )
                self._emit_budget_attention(
                    "exhausted",
                    "Run budget exhausted — the agent is wrapping up",
                    f"{pct}% of the run budget ({ceiling}) is spent. The "
                    f"agent finishes the current step and starts no new "
                    f"work. Raise agent.run_budget_usd / run_budget_s or "
                    f"start a new run to continue.")
            else:
                parts.append(
                    "# Run budget wind-down (auto-injected)\n"
                    f"- {pct}% of this run's budget is spent. Start wrapping "
                    "up NOW: complete or checkpoint the current step, "
                    "commit/save state, mark tasks, and summarise open work. "
                    "If a follow-up run makes sense, note exactly where to "
                    "resume. Do not start new large work items."
                )
                self._emit_budget_attention(
                    "wind_down",
                    "Run budget at 80% — the agent is winding down",
                    f"{pct}% of the run budget ({ceiling}) is spent. The "
                    f"agent is checkpointing and will stop soon.")
        return "\n\n".join(p for p in parts if p)

    def _cost_hard_cap(self) -> float:
        """Per-turn hard cost ceiling in USD — a runaway circuit-breaker, NOT a
        cumulative session budget. The default is deliberately very high so it
        only ever trips on a loop gone wrong, never on normal use. 0 disables
        it. Configurable via settings ``agent.cost_hard_limit_usd``."""
        try:
            from delfin.user_settings import load_settings as _ls
            v = ((_ls() or {}).get("agent", {}) or {}).get(
                "cost_hard_limit_usd", 50.0)
            return max(0.0, float(v))
        except Exception:
            return 50.0

    def reset_cycle(
        self,
        mode: str | None = None,
        preserve_messages: bool = False,
    ) -> None:
        """Reset the engine for a new work cycle.

        Parameters
        ----------
        mode : str, optional
            If given and different from current, switch to this mode
            (reloads role definitions and system prompts).
        preserve_messages : bool
            If True, KEEP the message history. Used for mode-switches
            mid-conversation (dashboard → solo) so the new agent sees
            the user's prompt without forcing a re-paste. Defaults to
            False for backward compatibility (true reset).
        """
        # Kill the persistent CLI process so a fresh one starts
        if hasattr(self.client, "kill"):
            self.client.kill()
        if not preserve_messages:
            self.messages.clear()
            # Every session-scoped field, from the one declaration that
            # export and restore also read. Hand-clearing a subset here was
            # the third face of the same drift: three of the eight ledgers
            # were listed, so a brand-new conversation started holding the
            # previous one's tool-name set and a delegation flag already
            # satisfied — which suppressed the delegation reminder for work
            # that had never delegated anything.
            self._reset_session_fields()
            self.session_id = self._fresh_session_id()  # fresh session for new cycle
            self._sync_task_session()
            self._foreign_tasks_shown = False  # new session → notice re-armed
        else:
            # Mode-switch with continuity: keep messages so the new
            # agent sees the user's prompt history, but reset
            # role-tracking since the new mode has a different role
            # route. role_outputs and compaction_summaries are
            # role-keyed and stale for the new mode — drop them.
            self._reset_session_fields(role_scoped_only=True)
        self.loader.reset_session_prompt_state(
            f"engine-session-{self._prompt_session_serial}"
        )
        self._prompt_session_serial += 1
        self._stop_requested = False
        if mode is not None and mode != self.mode:
            self._load_mode(mode)

    @dataclasses.dataclass(frozen=True)
    class _SessionField:
        """One piece of session-scoped state that has to survive a resume.

        ``dump`` makes it JSON-safe, ``load`` takes it back including the
        None case (so a session saved before the field existed loads with
        the same default a fresh session starts from), and ``reset`` gives
        the value a brand-new conversation starts from. All three live
        together because all three were maintained separately and all
        three had drifted.

        ``section`` says where the value rides in the saved dict:
        ``"state"`` is a top-level key, ``"evidence"`` is inside the
        evidence sub-dict the guards read.

        ``role_scoped`` marks the values a mode-switch invalidates — they
        are keyed by role, and the new mode has a different role route.
        """
        attr: str
        key: str
        section: str
        dump: "Callable[[Any], Any]"
        load: "Callable[[Any], Any]"
        reset: "Callable[[], Any]"
        role_scoped: bool = False

    @dataclasses.dataclass(frozen=True)
    class RestoreReport:
        """What a restore actually brought back.

        ``_restore_evidence`` swallowed every per-field exception into an
        empty default with no counter, no log and no return value, and
        ``restore_state`` returned None — so a restore that recovered two
        fields out of twenty was indistinguishable from a complete one,
        and the session went on to judge as if it had the record.
        """
        schema_version: int
        restored: tuple[str, ...] = ()
        missing: tuple[str, ...] = ()
        failed: tuple[str, ...] = ()
        migrations: tuple[str, ...] = ()

        @property
        def complete(self) -> bool:
            return not self.missing and not self.failed

        def summary(self) -> str:
            parts = [
                f"schema v{self.schema_version}",
                f"{len(self.restored)} field(s) restored",
            ]
            if self.missing:
                parts.append("missing: " + ", ".join(sorted(self.missing)))
            if self.failed:
                parts.append("failed: " + ", ".join(sorted(self.failed)))
            if self.migrations:
                parts.append("migrated: " + ", ".join(self.migrations))
            return "; ".join(parts)

    # The single declaration. Export, restore, reset and the coverage test
    # all read it, because the failure here was never one missing field --
    # it was a hand-written list beside a set of attributes maintained by
    # need, with nothing noticing when they drifted. They had drifted three
    # times: twice in the export (the location scanner ran disabled while
    # holding a populated ledger, and a restored PASS faced no test-evidence
    # veto) and once in the reset (a new conversation inherited the previous
    # one's tool-name set and a satisfied delegation flag).
    _SESSION_FIELDS: "tuple[_SessionField, ...]" = (
        # The language this conversation runs in, set by its first message.
        # Carried, because a resumed session is the SAME conversation: a
        # run picked up tomorrow that forgot which language it was started
        # in would answer the resumed half differently from the half
        # already on the page, which is the drift the pin exists to stop.
        _SessionField(
            "_session_language", "session_language", "state",
            str, str, lambda: ""),
        # -- the ledgers the guards judge against -------------------------
        _SessionField(
            "_last_observed_files", "observed_files", "evidence",
            lambda v: sorted(v or ()), lambda v: set(v or ()), set),
        # Whether a ledger EXISTS is not the same question as whether it is
        # empty -- an empty one means nothing was read, which is a fact.
        # Kept as its own value for that reason, and exported for the same.
        _SessionField(
            "_observed_ledger_available", "observed_ledger_available",
            "evidence", bool, bool, lambda: False),
        _SessionField(
            "_exec_commands_session", "exec_commands", "evidence",
            lambda v: list(v or ()), lambda v: list(v or ()), list),
        _SessionField(
            "_session_tool_names", "session_tool_names", "evidence",
            lambda v: sorted(v or ()), lambda v: set(v or ()), set),
        _SessionField(
            "_delegation_satisfied", "delegation_satisfied", "evidence",
            bool, bool, lambda: False),
        # The ninth and tenth ledger, found by the enumeration below the
        # moment it started walking the engine instead of comparing the
        # declaration to a copy of itself. Both have exactly the shape of
        # the eight above and neither survived a resume or a new cycle:
        # the task-list reminder re-fired for a list that already existed,
        # and the "do NOT repeat these" error memory came back empty.
        _SessionField(
            "_tasklist_satisfied", "tasklist_satisfied", "evidence",
            bool, bool, lambda: False),
        _SessionField(
            "_session_errors", "session_errors", "evidence",
            lambda v: [dict(x) for x in (v or ()) if isinstance(x, dict)],
            lambda v: [dict(x) for x in (v or ()) if isinstance(x, dict)]
            if isinstance(v, list) else [], list),
        _SessionField(
            "role_verdicts", "role_verdicts", "evidence",
            lambda v: dict(v or {}), lambda v: dict(v or {}), dict,
            role_scoped=True),
        _SessionField(
            "role_test_evidence", "role_test_evidence", "evidence",
            lambda v: {k: list(x) for k, x in (v or {}).items()},
            lambda v: {k: list(x) for k, x in (v or {}).items()}
            if isinstance(v, dict) else {}, dict, role_scoped=True),
        _SessionField(
            "_trimmed_chars_since_floor", "trimmed_chars_since_floor",
            "evidence", lambda v: int(v or 0), lambda v: int(v or 0),
            lambda: 0),
        # -- top-level session state --------------------------------------
        _SessionField(
            "current_role_index", "role_index", "state",
            lambda v: int(v or 0), lambda v: int(v or 0), lambda: 0,
            role_scoped=True),
        _SessionField(
            "role_outputs", "role_outputs", "state",
            lambda v: dict(v or {}), lambda v: dict(v or {}), dict,
            role_scoped=True),
        _SessionField(
            "compaction_summaries", "compaction_summaries", "state",
            lambda v: dict(v or {}), lambda v: dict(v or {}), dict,
            role_scoped=True),
        _SessionField(
            "token_usage", "token_usage", "state",
            lambda v: dict(v or {}),
            # Sessions saved before the cached-token counter existed have
            # no such key; the estimate must not read None for it.
            lambda v: {"input": 0, "output": 0, "cached": 0, **(v or {})},
            lambda: {"input": 0, "output": 0, "cached": 0}),
        _SessionField(
            "cost_usd", "cost_usd", "state",
            lambda v: float(v or 0.0), lambda v: float(v or 0.0),
            lambda: 0.0),
        # How much of that total could be measured at all. Carried with it
        # for the same reason the run clock is carried as elapsed: a resume
        # that dropped these would restore a partial spend figure and
        # present it as a ceiling the run is under, when the turns it
        # cannot see are the ones that were never priced.
        _SessionField(
            "_measured_cost_turns", "measured_cost_turns", "state",
            lambda v: int(v or 0), lambda v: int(v or 0), lambda: 0),
        _SessionField(
            "_non_billing_turns", "non_billing_turns", "state",
            lambda v: int(v or 0), lambda v: int(v or 0), lambda: 0),
        _SessionField(
            "_unpriced_turns", "unpriced_turns", "state",
            lambda v: int(v or 0), lambda v: int(v or 0), lambda: 0),
        _SessionField(
            "_project_dir", "project_dir", "state",
            lambda v: str(v or ""), lambda v: str(v or ""), lambda: ""),
        _SessionField(
            "_last_input_tokens", "last_input_tokens", "state",
            lambda v: int(v or 0), lambda v: int(v or 0), lambda: 0),
        _SessionField(
            "_system_prompt_chars", "system_prompt_chars", "state",
            lambda v: int(v or 0), lambda v: int(v or 0), lambda: 0),
        # The run clock, carried as ELAPSED rather than as a timestamp: the
        # wall-clock run budget is measured from _run_started_at, and a
        # resume that reset it to now handed the run its whole time budget
        # again -- so resuming was a way to launder a time budget on an
        # unattended run. Saving the elapsed seconds and setting the start
        # to "now minus that" continues the same clock.
        _SessionField(
            "_run_started_at", "run_elapsed_s", "state",
            lambda v: max(0.0, _time.time() - float(v or 0.0)),
            lambda v: _time.time() - max(0.0, float(v or 0.0)),
            _time.time),
        # The cost baseline the next outcome's delta is measured from. Zero
        # after a resume meant the first outcome record booked the WHOLE
        # session's spend as that one turn, which the activity view then
        # summed on top of the records already written for it.
        _SessionField(
            "_last_outcome_cost", "last_outcome_cost", "state",
            lambda v: float(v or 0.0), lambda v: float(v or 0.0),
            lambda: 0.0),
        # What this session's turns delegated. Carried for the same reason
        # the price-state counters are: a resumed run's spend figure means
        # nothing without the half of it that was delegated, and dropping
        # it would report a resumed session as cheaper than it was --
        # under a budget gate that reads the same figure. The per-TURN
        # bucket inside it is deliberately not carried; DelegateLedger's
        # loader zeroes it, because a resumed session has no turn in
        # flight for a previous process's delegates to be charged to.
        _SessionField(
            "_delegate_spend", "delegate_spend", "state",
            lambda v: (v if isinstance(v, _agent_metrics.DelegateLedger)
                       else _agent_metrics.DelegateLedger()).as_dict(),
            _agent_metrics.DelegateLedger.from_dict,
            _agent_metrics.DelegateLedger),
    )

    def _dump_field(self, spec: "AgentEngine._SessionField") -> Any:
        """getattr default + a swallowed failure: engines built for a
        targeted test never ran the full __init__, and export must not be
        the thing that breaks."""
        try:
            return spec.dump(getattr(self, spec.attr, None))
        except Exception:
            try:
                return spec.dump(spec.reset())
            except Exception:
                return None

    def _reset_session_fields(self, *, role_scoped_only: bool = False) -> None:
        """Set every declared field back to what a fresh session starts
        from. ``role_scoped_only`` is the mode-switch case: the message
        history and the evidence for it stay, only the role-keyed values go.
        """
        for spec in self._SESSION_FIELDS:
            if role_scoped_only and not spec.role_scoped:
                continue
            try:
                fresh = spec.reset()
            except Exception:
                continue
            cur = getattr(self, spec.attr, None)
            # Clear containers IN PLACE where possible: other components
            # hold references to these, and a rebind would leave them
            # pointing at the previous cycle's contents.
            if (isinstance(cur, (dict, list, set))
                    and type(cur) is type(fresh) and not fresh):
                cur.clear()
            else:
                setattr(self, spec.attr, fresh)

    def export_state(self) -> dict:
        """Export engine state for session persistence.

        Every declared field is written from the declaration, so a call
        site cannot drop one by listing the keys it happens to know.
        """
        # getattr throughout: engines built for a targeted test never ran
        # the full __init__, and export must not be the thing that breaks.
        out: dict = {
            "mode": getattr(self, "mode", ""),
            "route": list(getattr(self, "route", None) or ()),
            "engine_messages": list(getattr(self, "messages", None) or ()),
            "session_id": getattr(self, "session_id", "") or "",
        }
        for spec in self._SESSION_FIELDS:
            if spec.section == "state":
                out[spec.key] = self._dump_field(spec)
        # One sub-dict, written and read as a unit, so a new ledger is
        # added in one place instead of four.
        out["evidence"] = self._export_evidence()
        return out

    def _export_evidence(self) -> dict:
        """The ledgers the guards judge against.

        A resumed session got the whole conversation and none of the
        evidence -- while the "a ledger exists" flags came back True, which
        is the ENFORCING branch. So the first answer after a resume that
        restated a file:line from restored history was corrected as
        unsupported, and every "it works now" got an unverified caveat for
        commands that had run before the save.

        _trimmed_chars_since_floor is the sharpest of them: restore_state
        brings back _last_input_tokens, a PRE-trim provider snapshot, while
        the credit that offsets it resets to zero. A resumed session then
        reads its context as larger than it is and compacts early.
        """
        out: dict = {}
        for spec in self._SESSION_FIELDS:
            if spec.section == "evidence":
                out[spec.key] = self._dump_field(spec)
        # When it was taken. Nothing reads it yet, and that is deliberate:
        # a resumed session asserting "the tests pass now" is believed on a
        # ledger recorded against a tree that may have moved since, and the
        # answer to that is to say where the evidence CAME FROM. Discarding
        # stale evidence is what produced the false "unverified" caveats
        # that made exporting it necessary in the first place.
        out["saved_at"] = int(_time.time())
        return out

    def _load_declared_fields(
        self, data: dict, *, section: str,
    ) -> tuple[list[str], list[str], list[str]]:
        """Read one section back. Missing or malformed means the fresh
        default, never an exception: a session saved before a field existed
        must still load.

        Returns (restored, missing, failed) attribute names — the counters
        that used to not exist, which is why a nearly empty restore looked
        exactly like a complete one.
        """
        if section == "evidence":
            src = data.get("evidence")
            if not isinstance(src, dict):
                src = {}
        else:
            src = data
        restored: list[str] = []
        missing: list[str] = []
        failed: list[str] = []
        for spec in self._SESSION_FIELDS:
            if spec.section != section:
                continue
            if spec.key not in src:
                missing.append(spec.attr)
                try:
                    setattr(self, spec.attr, spec.reset())
                except Exception:
                    failed.append(spec.attr)
                continue
            try:
                setattr(self, spec.attr, spec.load(src.get(spec.key)))
                restored.append(spec.attr)
            except Exception:
                failed.append(spec.attr)
                try:
                    setattr(self, spec.attr, spec.reset())
                except Exception:
                    pass
        return restored, missing, failed

    def _restore_evidence(self, data: dict) -> "AgentEngine.RestoreReport":
        """Read the guard ledgers back, driven by the same declaration as
        the export so the two cannot drift apart."""
        restored, missing, failed = self._load_declared_fields(
            data, section="evidence")
        return AgentEngine.RestoreReport(
            schema_version=int(data.get("schema_version") or 1),
            restored=tuple(restored), missing=tuple(missing),
            failed=tuple(failed),
        )

    def _seed_client_observed_files(self) -> None:
        """Hand the restored read-ledger to the client.

        ``restore_state`` deliberately leaves the client alone, and the
        client mints ``_observed_files_session`` fresh on its first turn.
        Without this the restored set was overwritten by an empty one
        before any guard read it, and the grounding check then flagged
        every file:line the answer restated from the restored history.
        Seeded only when the save says a ledger really existed, so a
        backend that keeps none does not acquire a fictitious one.
        """
        if not getattr(self, "_observed_ledger_available", False):
            return
        client = getattr(self, "client", None)
        if client is None:
            return
        try:
            files = set(getattr(self, "_last_observed_files", None) or ())
            current = getattr(client, "_observed_files_session", None)
            if isinstance(current, set):
                current |= files
            else:
                client._observed_files_session = files
        except Exception:
            pass

    def restore_state(self, data: dict) -> "AgentEngine.RestoreReport":
        """Restore engine state from a saved session.

        The client and prompt loader remain as-is; only conversation
        state is restored so ``--resume`` can continue the CLI session.

        Returns a :class:`RestoreReport` naming what came back and what did
        not. Raises ``session_store.SessionSchemaError`` when the file was
        written by a NEWER schema than this code understands — a stated
        refusal, because reading it with today's loaders would silently
        drop whatever the newer version added.
        """
        from . import session_store as _ss

        try:
            version = int(data.get("schema_version") or 1)
        except (TypeError, ValueError):
            version = 1
        if version > _ss.SESSION_SCHEMA_VERSION:
            raise _ss.SessionSchemaError(
                f"session was written with schema v{version}; this build "
                f"understands v{_ss.SESSION_SCHEMA_VERSION}. Refusing to "
                f"load it rather than silently dropping what it carries."
            )
        migrations: tuple[str, ...] = ()
        if version < _ss.SESSION_SCHEMA_VERSION:
            data, migrations = _ss.migrate_session_data(data, version)

        restored: list[str] = []
        missing: list[str] = []
        drifted: list[str] = []

        if data.get("mode"):
            restored.append("mode")
            if data["mode"] != self.mode:
                self._load_mode(data["mode"])
        else:
            missing.append("mode")

        # ``route`` was exported and read by nobody: the mode reload above
        # re-derives it, so the saved value looked carried and was not.
        # It is the only record of the route the session ACTUALLY ran, and
        # a pack edited between two sessions resumes on a different one
        # with nothing said. The current pack still decides what runs —
        # reinstating a route the pack no longer defines would resume into
        # roles that are not there — but the difference is named.
        saved_route = list(data.get("route") or ())
        if saved_route and saved_route != list(getattr(self, "route", ()) or ()):
            drifted.append(
                f"route (session ran {'→'.join(str(r) for r in saved_route)}; "
                f"this pack yields "
                f"{'→'.join(str(r) for r in (self.route or ())) or 'nothing'})"
            )
        elif saved_route:
            restored.append("route")

        if "engine_messages" in data:
            self.messages = list(data.get("engine_messages") or [])
            restored.append("engine_messages")
            # Stated, not inferred. The cross-session task notice needs to
            # know a conversation came off disk; a message count cannot
            # tell that apart from a session's second turn.
            self._history_restored = bool(self.messages)
        else:
            self.messages = []
            missing.append("engine_messages")

        self.session_id = str(data.get("session_id", "") or "")
        if self.session_id:
            restored.append("session_id")
        else:
            missing.append("session_id")
        # Bind the task list to the restored id. The dashboard patched this
        # from outside afterwards; the headless path did not, so a headless
        # resume filtered tasks on a throwaway id and saw zero open tasks.
        self._sync_task_session()

        s_restored, s_missing, s_failed = self._load_declared_fields(
            data, section="state")
        ev = self._restore_evidence(data)
        # The restored read-ledger has to reach the client before the first
        # turn's post-stream update, or it is discarded unread.
        self._seed_client_observed_files()

        report = AgentEngine.RestoreReport(
            schema_version=version,
            restored=tuple(restored + s_restored + list(ev.restored)),
            missing=tuple(missing + s_missing + list(ev.missing)),
            failed=tuple(drifted + s_failed + list(ev.failed)),
            migrations=migrations,
        )

        # Crash recovery: a surviving mid-turn checkpoint means the
        # previous process died inside a turn, so that turn's work was
        # never committed. Inject the recovery note using the house
        # user-note + assistant-ack pair (same shape as the compaction
        # summary) so message alternation survives _sanitize_messages.
        try:
            from . import session_store as _ss_ckpt
            _note = _ss_ckpt.consume_crash_recovery_note(self.session_id)
            if _note:
                self.messages.append({"role": "user", "content": _note})
                self.messages.append({
                    "role": "assistant",
                    "content": "Understood. I will verify the workspace "
                               "state before continuing.",
                })
        except Exception:
            pass
        self.last_restore_report = report
        return report

    def available_modes(self) -> list[str]:
        """Return list of available mode IDs."""
        return self.loader.available_modes()

    def record_cycle_outcome(
        self,
        verdict: str,
        user_task: str,
        error_type: str | None = None,
        denied_commands: list[str] | None = None,
        start_time: float | None = None,
    ) -> dict[str, str]:
        """Record a cycle outcome and update the provider profile.

        Returns a dict of profile changes made (for transparency logging).
        """
        from delfin.agent.outcome_tracker import CycleOutcome, append_outcome
        from delfin.agent.provider_profile import update_from_outcome

        duration = 0.0
        if start_time:
            import time
            duration = time.monotonic() - start_time

        task_class = ""
        try:
            route_info = self.recommend_task_route(
                user_task, self.mode,
                is_delfin_workspace=self._is_delfin_workspace,
            )
            task_class = route_info.get("task_class", "")
        except Exception:
            pass

        # A6 — honest delta. ``self.cost_usd`` is the cumulative engine cost
        # across every turn of this session. The Δ for THIS turn is the new
        # total minus the total at the previous outcome — that's the value
        # the activity tab should sum across cycles to get real spend.
        prev_total = float(getattr(self, "_last_outcome_cost", 0.0) or 0.0)
        cost_delta = max(0.0, float(self.cost_usd) - prev_total)
        self._last_outcome_cost = float(self.cost_usd)

        outcome = CycleOutcome(
            task=user_task[:200],
            provider=self.provider,
            model=getattr(self.client, "model", ""),
            mode=self.mode,
            verdict=verdict,
            cost_usd=self.cost_usd,
            cost_usd_delta=round(cost_delta, 6),
            duration_s=round(duration, 1),
            retries=0,
            denied_commands=denied_commands or [],
            error_type=error_type,
            task_class=task_class,
            timestamp="",
        )
        try:
            append_outcome(outcome)
        except Exception:
            pass
        try:
            return update_from_outcome(self.provider, outcome)
        except Exception:
            return {}

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
    def thinking_budget_for_role(
        role_id: str,
        task_class: str = "",
        task_text: str = "",
    ) -> int:
        """Return the recommended thinking budget for a role.

        Applies both provider-specific and task-complexity multipliers.
        Simple tasks get 40% of base, complex tasks get 150%.
        """
        base = _ROLE_THINKING_BUDGETS.get(role_id, _DEFAULT_THINKING_BUDGET)

        # Task complexity multiplier (Feature 7)
        complexity = AgentEngine.classify_task_complexity(task_text)
        complexity_mult = _COMPLEXITY_THINKING_MULT.get(complexity, 1.0)

        try:
            from delfin.agent.provider_profile import (
                load_provider_profile,
                load_task_profile,
            )
            _prov = getattr(AgentEngine, "_active_provider", "claude")
            _profile = load_provider_profile(_prov)
            _task_profile = load_task_profile(_prov, task_class)
            provider_mult = _task_profile.get(
                "thinking_budget_mult",
                _profile.get("thinking_budget_mult", 1.0),
            )
            # Clamp multiplier to safe range
            provider_mult = max(0.0, min(3.0, provider_mult))
            return int(base * provider_mult * complexity_mult)
        except Exception:
            return int(base * complexity_mult)

    @staticmethod
    def is_greeting(text: str) -> bool:
        """True for short greeting/acknowledgement messages ("Hallo",
        "danke!").  Used for the fast-path: minimal reasoning budget AND
        a no-tools hint so a greeting never triggers tool calls or
        confirm prompts."""
        lower = (text or "").lower().strip()
        return bool(lower) and len(lower) < 40 and any(
            lower == g or lower.startswith(g + " ")
            or lower.startswith(g + "!") or lower.startswith(g + ",")
            for g in _GREETING_PATTERNS
        )

    @staticmethod
    def is_bare_greeting(text: str) -> bool:
        """True only when the message is a greeting and NOTHING else.

        Deliberately stricter than is_greeting, which is used to pick a
        reasoning budget and is happy with a greeting PREFIX. That is far
        too loose for the decision here: "hi, lies bitte buchungen.csv" and
        "hallo, kannst du die tabelle prüfen?" both satisfy it, and
        suppressing the tool surface on either would leave the agent unable
        to do the work it was just asked to do.

        So every word has to be a greeting word, the agent's name, or
        punctuation, and a question mark disqualifies outright: a question
        is a request even when it is polite.
        """
        lower = (text or "").strip().lower()
        if not lower or len(lower) > 40 or "?" in lower:
            return False
        cleaned = re.sub(r"[^\w\s]+", " ", lower, flags=re.UNICODE)
        # Multi-word greetings first, so "guten morgen" is not judged as
        # two unknown words.
        for phrase in sorted(_GREETING_PATTERNS, key=len, reverse=True):
            if " " in phrase:
                cleaned = cleaned.replace(phrase, " ")
        words = [w for w in cleaned.split() if w]
        if not words:
            return bool(lower)      # punctuation or an emoji alone
        allowed = set(_GREETING_PATTERNS) | {"delfin", "du", "dir", "mal"}
        return all(w in allowed for w in words)

    @staticmethod
    def classify_task_complexity(text: str) -> str:
        """Classify task complexity: simple / moderate / complex."""
        lower = (text or "").lower().strip()
        if not lower:
            return "moderate"
        # Simple: short questions, read-only operations
        # Greetings / acknowledgements: short message that *starts* with a
        # greeting token → no reasoning needed (fast, cheap turn).
        if AgentEngine.is_greeting(lower):
            return "simple"
        if len(lower) < 60 and any(p in lower for p in _SIMPLE_TASK_PATTERNS):
            return "simple"
        # Complex: multi-file refactors, architecture changes
        if any(re.search(p, lower) for p in _COMPLEX_TASK_PATTERNS):
            return "complex"
        # Long messages with multiple file paths → likely complex
        file_refs = re.findall(r"[\w./]+\.py\b", lower)
        if len(file_refs) >= 3:
            return "complex"
        return "moderate"

    @staticmethod
    def model_for_task(role_id: str, task_text: str = "") -> str:
        """Intelligent model routing based on role AND task complexity.

        Simple tasks (git status, read file) → haiku (cheap, fast).
        Moderate tasks → sonnet or user's choice.
        Complex tasks → user's choice (often opus).
        """
        role_model = _ROLE_MODEL_MAP.get(role_id, "auto")
        if role_model != "auto":
            return role_model
        # Solo/builder: route by complexity
        if role_id in ("solo_agent", "dashboard_agent") and task_text:
            complexity = AgentEngine.classify_task_complexity(task_text)
            if complexity == "simple":
                return "haiku"
        return "auto"

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

        # Roles that act on prior outputs need more detail than read-only roles.
        # Builder and Test need the full SM plan + critic findings.
        # Research/Reviewer/Runtime only need compact summaries.
        _FULL_PRIOR_ROLES = {"builder_agent", "test_agent"}
        _PRIOR_CHAR_LIMITS = {
            "session_manager": 4000,
            "critic_agent": 3000,
            "research_agent": 1500,
            "runtime_agent": 1500,
            "reviewer_agent": 1500,
        }

        if role in _FULL_PRIOR_ROLES and self.role_outputs:
            # Builder/Test get truncated-but-substantial prior outputs
            parts = ["--- Prior Role Outputs ---"]
            for rid, output in self.role_outputs.items():
                limit = _PRIOR_CHAR_LIMITS.get(rid, 3000)
                text = output[:limit]
                if len(output) > limit:
                    text += "\n... [truncated]"
                parts.append(f"### {rid}\n{text}")
            prior_summary = "\n\n".join(parts)
        elif self.compaction_summaries:
            # Other roles get compact summaries (saves tokens)
            parts = ["Prior outputs (compact):"]
            for rid, summary in self.compaction_summaries.items():
                if summary:
                    parts.append(f"### {rid}\n{summary}")
            if len(parts) > 1:
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
                f"{prior_summary}\n\n"
                f"IMPORTANT: Confirm your research questions with the user BEFORE searching. "
                f"Use QUESTION: to list your 2-3 planned research questions and ask if "
                f"anything else should be investigated.\n\n"
                f"Focus on actionable findings for the Builder. "
                f"Use WebSearch and WebFetch. Max 5 searches."
            ),
            "critic_agent": (
                f"Review the plan and/or changes for this task:\n\n{user_task}\n\n"
                f"{locked_contract}\n\n"
                f"{prior_summary}\n\n"
                f"IMPORTANT: After your analysis, present your top findings to the user "
                f"with QUESTION: and ask which the Builder should prioritize.\n\n"
                f"Produce a structured REVIEW. Focus on critical and major issues. "
                f"Explicitly flag goal drift, weak success proxies, and missing stage gates."
            ),
            "runtime_agent": (
                f"Review the runtime implications for this task:\n\n{user_task}\n\n"
                f"{locked_contract}\n\n"
                f"{prior_summary}\n\n"
                f"Produce a structured RUNTIME REVIEW. Check local vs cluster behavior and "
                f"whether the runtime-facing success metric is actually sufficient."
            ),
            "builder_agent": (
                f"Implement the following task:\n\n{user_task}\n\n"
                f"{locked_contract}\n\n"
                f"{prior_summary}\n\n"
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
                f"{prior_summary}\n\n"
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

    # `suggest_mode` lived here. It ranked modes on a ladder (quick <
    # reviewed < cluster < full) and proposed a climb up it. There is no
    # ladder any more: the three modes that remain differ by WHICH FOLDER
    # the session works in, which is not something to escalate into. It
    # returned `cluster` and `full` to a user who cannot select either.
    #
    # Nothing in the package called it -- checked across delfin/ -- which
    # is exactly why it survived the retirement untouched. Deleted rather
    # than repaired: a suggestion function over a set with no order has
    # nothing left to suggest, and leaving it would hand the next caller a
    # mode that does not exist.

    @staticmethod
    def recommend_task_route(
        user_message: str, current_mode: str = "dashboard",
        *, is_delfin_workspace: bool = True,
    ) -> dict[str, Any]:
        """Recommend the cheapest capable mode for a task.

        Returns a dict with:
        - ``mode``: recommended mode
        - ``task_class``: dashboard / chemistry / coding / general
        - ``intent``: operate / research / question / change
        - ``confidence``: low / medium / high
        - ``reasons``: short rationale list
        - ``risk_flags``: reviewed / cluster / full flags when present

        When ``is_delfin_workspace`` is False, chemistry / DELFIN-file
        escalation patterns are skipped so generic projects don't get
        flagged as "review-needed" just because they mention "orca" or
        contain a file called ``calculators.py``.
        """
        text = user_message or ""
        lower = text.lower()
        stripped = lower.strip()

        def _contains_any(keywords: tuple[str, ...]) -> list[str]:
            return [kw for kw in keywords if kw in lower]

        dashboard_hits = (
            _contains_any(_DASHBOARD_KEYWORDS) if is_delfin_workspace else []
        )
        chemistry_hits = (
            _contains_any(_CHEMISTRY_KEYWORDS) if is_delfin_workspace else []
        )
        code_change_hits = _contains_any(_CODE_CHANGE_KEYWORDS)
        code_question_hits = _contains_any(_CODE_QUESTION_KEYWORDS)
        code_hint_hits = _contains_any(_FILE_OR_CODE_HINTS)
        reviewed_hits = _contains_any(_REVIEWED_RISK_KEYWORDS)
        cluster_hits = _contains_any(_CLUSTER_RISK_KEYWORDS)
        full_hits = _contains_any(_FULL_RISK_KEYWORDS)

        if is_delfin_workspace:
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

        # Adaptive routing: prefer a higher-success mode when this
        # installation's recorded rate for the task class is poor.
        #
        # OFF unless asked for, because the number it reads is not
        # trustworthy and the rewrite it performs is invisible. Measured
        # on this machine: outcome_history.jsonl holds 55 records, all
        # PASS, in modes solo and dashboard -- not one `quick` record
        # exists -- while provider_profile_state.json holds kit solo
        # 0.005 and kit coding 0.007. Two accumulators for one fact, with
        # no reconciliation: update_from_outcome never reads the history,
        # and the two writes sit in separate try/except blocks so the
        # profile mutates even when the audit write fails. The EMA also
        # seeds at 0.5, below the threshold, so a class recovering from a
        # single FAIL escalates until several consecutive passes.
        #
        # It had already corrupted two tests, both carrying written-up
        # workarounds saying the assertion measures the machine rather
        # than the code. A mechanism that rewrites a decision has to be
        # visible and resettable first; /profile now shows and resets it.
        try:
            _adaptive = bool(((_load_settings().get("agent") or {})
                              .get("routing") or {}).get(
                                  "adaptive_escalation", False))
        except Exception:
            _adaptive = False
        try:
            from delfin.agent.provider_profile import load_provider_profile
            # The provider actually in use. This fell back to "claude",
            # so a KIT session was routed by numbers recorded for
            # Anthropic -- which is why the live reason string said 17%
            # while KIT's own rate for the same class was 0.7%.
            _prov = getattr(AgentEngine, "_active_provider", "") or "claude"
            _profile = load_provider_profile(_prov) if _adaptive else {}
            _task_perf = _profile.get("task_performance", {}).get(task_class, {})
            _task_rate = _task_perf.get("success_rate")
            if (
                _adaptive
                and isinstance(_task_rate, (int, float))
                and mode == "quick"
                and task_class in ("chemistry", "coding")
                and _task_rate < 0.75
            ):
                reasons.append(
                    f"adaptive: escalate {task_class} task from quick to reviewed "
                    f"({_prov} recorded {_task_rate:.0%} success for this class)"
                )
                mode = "reviewed"
        except Exception:
            pass

        return {
            # The scoring above still reasons in terms of the retired
            # pipelines, and that reasoning is worth keeping -- "this looks
            # like a release" is a real signal, and it survives in
            # risk_flags and reasons. What must not survive is a `mode` the
            # user cannot select: this value is a public part of the
            # contract, and naming `quick` or `full` here invites the next
            # caller to act on a mode that no longer exists.
            "mode": _migrate_mode(mode),
            "task_class": task_class,
            "intent": intent,
            "confidence": confidence,
            "reasons": reasons[:5],
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
        self._capture_current_role_summary()
        self.messages.clear()
        if hasattr(self.client, "kill"):
            self.client.kill()
        # Only the CLI needs a fresh session here (its process is restarted so
        # the next role gets a new --append-system-prompt). API backends keep
        # the same session across role compaction so tasks stay scoped to it.
        if (getattr(self, "backend", "") or "").lower() == "cli":
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
        """Extract the canonical ``**status:**`` verdict if present.

        Deliberately tolerant: trailing text ("reject — see findings") and
        word variants ("rejected"/"approved") still count. The old
        end-anchored regex extracted '' for those, and an empty status made
        a REJECTING critic auto-continue — the gate failed OPEN on exactly
        the format slips weak models make most.
        """
        match = re.search(
            r"^\*\*status:\*\*\s*"
            r"(approve_with_risks|approved|approve|rejected|reject)\b",
            agent_output or "",
            flags=re.IGNORECASE | re.MULTILINE,
        )
        if not match:
            return ""
        value = match.group(1).lower()
        return {"approved": "approve", "rejected": "reject"}.get(value, value)

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
    def evaluate_role_gate(
        role_id: str,
        output: str,
        structured_verdict: dict | None = None,
    ) -> tuple[str, str, str]:
        """Return a communication-gate decision for a completed role.

        Returns ``(action, gate_type, message)`` where action is one of:
        - ``continue``: safe to auto-advance
        - ``pause``: stop and ask the user to review/approve

        ``structured_verdict`` is the parsed ``report_verdict`` tool call
        from the role's turn, when one was made. It wins over prose: a tool
        argument cannot suffer the formatting slips that used to make a
        rejecting critic auto-continue.
        """
        text = output or ""
        status = ""
        if isinstance(structured_verdict, dict):
            sv = str(structured_verdict.get("status", "")).strip().lower()
            if sv in ("approve", "approve_with_risks", "reject"):
                status = sv
        if not status:
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
            if status == "reject":
                return (
                    "pause",
                    "review",
                    "reported `status: reject`; review the findings before continuing.",
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
        # Solo and dashboard agents have no structured format — skip validation
        if role_id in ("solo_agent", "dashboard_agent"):
            return []

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
