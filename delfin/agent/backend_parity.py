"""Backend capability parity for the DELFIN agent.

The full agent surface is carried by the OpenAI-compatible API client
(``OpenAIClient`` in ``api_client.py``): the agentic tool loop, the
sandboxed file/bash executors, subagent orchestration, task tracking,
verify enforcement, the undo journal, web tools and MCP integration all
live in its ``stream_message`` loop and the shared ``_DocToolExecutor``.
The other client backends only relay model output:

- ``CLIClient`` delegates the whole surface to an external agent CLI
  process (that process brings its own tools, permission handling and
  MCP support; DELFIN's gates and journals never see its actions).
- ``APIClient`` (Anthropic SDK) streams text/thinking only — its
  request carries no tools parameter at all.
- ``CodexCLIClient`` spawns a per-turn external process whose sandbox
  flags are mapped from the DELFIN permission profile; nothing else of
  the DELFIN tool surface is forwarded.

This module states those differences in one queryable place so a
session can be told ONCE, at start, what its backend cannot do —
instead of each gap surfacing only when a task silently lacks it.

All functions are pure, take no environment input, and never raise.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Capability keys + human-readable labels
# ---------------------------------------------------------------------------

CAP_TOOL_LOOP = "tool_loop"
CAP_FILE_TOOLS = "file_tools"
CAP_BASH = "bash"
CAP_SUBAGENTS = "subagents"
CAP_TASK_TOOLS = "task_tools"
CAP_MEMORY_RECALL = "memory_recall"
CAP_VERIFY = "verify_enforcement"
CAP_UNDO_JOURNAL = "undo_journal"
CAP_WEB_TOOLS = "web_tools"
CAP_MCP = "mcp"

CAPABILITY_LABELS: dict[str, str] = {
    CAP_TOOL_LOOP: "agentic tool loop (gated tool execution)",
    CAP_FILE_TOOLS: "file tools (read_file/write_file/edit_file/multi_edit)",
    CAP_BASH: "shell execution (bash with deny-list and confirm gating)",
    CAP_SUBAGENTS: "subagent delegation and orchestration",
    CAP_TASK_TOOLS: "task tracking (task_create/task_update/task_list)",
    CAP_MEMORY_RECALL: "memory recall injection",
    CAP_VERIFY: "verify enforcement (auto-verify + verdict/evidence ledger)",
    CAP_UNDO_JOURNAL: "undo journal (change capture + undo_changes)",
    CAP_WEB_TOOLS: "web tools (web_search/web_fetch)",
    CAP_MCP: "MCP servers (tools/resources/prompts)",
}

#: Support levels used in the matrix.
YES = "yes"
CONDITIONAL = "conditional"
NO = "no"

#: Canonical backend identifiers (one per client class in api_client.py).
BACKEND_OPENAI = "openai"            # OpenAIClient (providers kit/ollama/openai)
BACKEND_CLI = "cli"                  # CLIClient (external agent CLI process)
BACKEND_ANTHROPIC_API = "anthropic_api"  # APIClient (Anthropic SDK streaming)
BACKEND_CODEX_CLI = "codex_cli"      # CodexCLIClient (per-turn codex exec)

#: The backend every other one is compared against.
REFERENCE_BACKEND = BACKEND_OPENAI

_BACKEND_LABELS: dict[str, str] = {
    BACKEND_OPENAI: "OpenAI-compatible API",
    BACKEND_CLI: "external agent CLI",
    BACKEND_ANTHROPIC_API: "Anthropic API (streaming only)",
    BACKEND_CODEX_CLI: "Codex CLI (external process)",
}

# ---------------------------------------------------------------------------
# The matrix. Every cell is (support, reason) and was derived from code
# inspection of api_client.py / engine.py / prompt_loader.py — the reason
# names the deciding symbol so the claim stays checkable as lines move.
# ---------------------------------------------------------------------------

_MATRIX: dict[str, dict[str, tuple[str, str]]] = {
    # OpenAIClient.stream_message: tool-round loop, executor dispatch,
    # auto-verify gate, MCP discovery. Conditions inside the same method:
    # tools are suppressed for models whose resolved capabilities report
    # supports_tools=False, weak-model profiles trim to a core tool set,
    # and every coding tool needs a KitToolPermissions policy (wired by
    # the kit/ollama factory paths in create_client; the plain openai
    # provider path constructs the client without one).
    BACKEND_OPENAI: {
        CAP_TOOL_LOOP: (CONDITIONAL, (
            "full agentic loop in OpenAIClient.stream_message; suppressed "
            "for models whose resolved capabilities report no native tool "
            "calling (model_capabilities.supports_tools) — those run "
            "chat-only; weak-model profiles trim to a core tool set"
        )),
        CAP_FILE_TOOLS: (CONDITIONAL, (
            "read tools (read_file/grep_file/list_files) are always "
            "advertised; write_file/edit_file/multi_edit/apply_patch need "
            "a workspace permissions policy — the kit and ollama factory "
            "paths wire one, the plain openai provider path does not"
        )),
        CAP_BASH: (CONDITIONAL, (
            "needs a workspace permissions policy; commands pass the "
            "deny-list / auto-allow / confirm gate and the secret and "
            "egress scans in _DocToolExecutor"
        )),
        CAP_SUBAGENTS: (CONDITIONAL, (
            "needs a permissions policy; the runner is bound via "
            "permissions.subagent_runner in "
            "OpenAIClient._attach_subagent_runner"
        )),
        CAP_TASK_TOOLS: (CONDITIONAL, (
            "need a permissions policy; the task executors return an "
            "error when permissions are missing"
        )),
        CAP_MEMORY_RECALL: (YES, (
            "the engine injects project/external/episodic memory into "
            "every prompt build; the remember/forget write tools are "
            "executable in the tool loop"
        )),
        CAP_VERIFY: (CONDITIONAL, (
            "auto-verify runs when Python files were edited and "
            "agent.auto_verify is not off; the structured verdict and "
            "test-evidence ledger are recorded per turn"
        )),
        CAP_UNDO_JOURNAL: (CONDITIONAL, (
            "pre-images are captured for changes made through the DELFIN "
            "executors; the undo_changes tool needs a workspace "
            "permissions policy"
        )),
        CAP_WEB_TOOLS: (CONDITIONAL, (
            "web_search/web_fetch are part of the coding tool surface — "
            "advertised only when a permissions policy is configured"
        )),
        CAP_MCP: (YES, (
            "MCP tools/resources/prompts are discovered and gated inside "
            "the tool loop; workspace-scoped server config needs a "
            "permissions policy"
        )),
    },
    # CLIClient: spawns a persistent external agent CLI process and only
    # relays its stream events. The external process executes its own
    # tools under its own permission mode; DELFIN forwards a tool
    # restriction (--allowedTools), extra dirs (--add-dir) and an MCP
    # config (--mcp-config) but is not in the execution path.
    BACKEND_CLI: {
        CAP_TOOL_LOOP: (CONDITIONAL, (
            "runs inside the external agent CLI process, which brings its "
            "own tool loop; DELFIN's loop, gates and per-turn budgets do "
            "not apply"
        )),
        CAP_FILE_TOOLS: (CONDITIONAL, (
            "the external CLI's own file tools (optionally restricted via "
            "--allowedTools); DELFIN's sandboxed executors and change "
            "previews are bypassed"
        )),
        CAP_BASH: (CONDITIONAL, (
            "the external CLI's own shell tool under its own permission "
            "mode; DELFIN's deny-list and secret/egress scans are not "
            "applied"
        )),
        CAP_SUBAGENTS: (CONDITIONAL, (
            "only what the external CLI itself provides; DELFIN's "
            "subagent/orchestrate tools are not wired to this client"
        )),
        CAP_TASK_TOOLS: (NO, (
            "the DELFIN task-store tools are dispatched only by the "
            "OpenAI-compatible tool loop"
        )),
        CAP_MEMORY_RECALL: (CONDITIONAL, (
            "memory is injected into the system prompt once at process "
            "spawn (--append-system-prompt); the remember/forget write "
            "tools are unavailable"
        )),
        CAP_VERIFY: (NO, (
            "auto-verify and the verdict/evidence ledger exist only in "
            "the OpenAI-compatible loop; the engine reads neither from "
            "CLI clients"
        )),
        CAP_UNDO_JOURNAL: (NO, (
            "change pre-images are captured only by the DELFIN executors; "
            "edits made by the external CLI never reach the journal"
        )),
        CAP_WEB_TOOLS: (CONDITIONAL, (
            "only the external CLI's own web features, if enabled there; "
            "DELFIN's web_search/web_fetch executors are unreachable"
        )),
        CAP_MCP: (CONDITIONAL, (
            "a server config is forwarded via --mcp-config when provided; "
            "the external CLI handles MCP itself, outside DELFIN's "
            "registry and gates"
        )),
    },
    # APIClient: direct Anthropic SDK streaming. Its request kwargs carry
    # model/max_tokens/system/messages (+thinking) and never a tools
    # parameter, so no tool can be called at all on this backend.
    BACKEND_ANTHROPIC_API: {
        CAP_TOOL_LOOP: (NO, (
            "APIClient.stream_message sends no tools parameter and has no "
            "executor loop — text/thinking streaming only"
        )),
        CAP_FILE_TOOLS: (NO, "no tools are advertised on this backend"),
        CAP_BASH: (NO, "no tools are advertised on this backend"),
        CAP_SUBAGENTS: (NO, "no tools are advertised on this backend"),
        CAP_TASK_TOOLS: (NO, "no tools are advertised on this backend"),
        CAP_MEMORY_RECALL: (YES, (
            "the engine rebuilds the system prompt (with memory blocks) "
            "every turn; only the remember/forget write tools are "
            "unavailable"
        )),
        CAP_VERIFY: (NO, (
            "auto-verify and the verdict/evidence ledger exist only in "
            "the OpenAI-compatible loop"
        )),
        CAP_UNDO_JOURNAL: (NO, (
            "nothing edits files on this backend and the journal is never "
            "written"
        )),
        CAP_WEB_TOOLS: (NO, "no tools are advertised on this backend"),
        CAP_MCP: (NO, "no MCP discovery or forwarding on this backend"),
    },
    # CodexCLIClient: per-turn external process. The DELFIN permission
    # profile maps to sandbox flags; nothing else of the tool surface is
    # forwarded and no MCP config is passed.
    BACKEND_CODEX_CLI: {
        CAP_TOOL_LOOP: (CONDITIONAL, (
            "delegated to a per-turn external codex process; DELFIN's "
            "loop and gates do not apply"
        )),
        CAP_FILE_TOOLS: (CONDITIONAL, (
            "the external process edits files under sandbox modes mapped "
            "from the DELFIN permission profile; DELFIN's executors and "
            "change previews are bypassed"
        )),
        CAP_BASH: (CONDITIONAL, (
            "the external process's own shell under its sandbox flags; "
            "DELFIN's bash gates are not applied"
        )),
        CAP_SUBAGENTS: (NO, (
            "no delegation wiring exists on this client"
        )),
        CAP_TASK_TOOLS: (NO, (
            "the DELFIN task-store tools are dispatched only by the "
            "OpenAI-compatible tool loop"
        )),
        CAP_MEMORY_RECALL: (YES, (
            "the engine's system prompt (with memory blocks) is embedded "
            "in every per-turn prompt; the remember/forget write tools "
            "are unavailable"
        )),
        CAP_VERIFY: (NO, (
            "auto-verify and the verdict/evidence ledger exist only in "
            "the OpenAI-compatible loop"
        )),
        CAP_UNDO_JOURNAL: (NO, (
            "change pre-images are captured only by the DELFIN executors; "
            "edits made by the external process never reach the journal"
        )),
        CAP_WEB_TOOLS: (NO, (
            "DELFIN's web tools are unreachable and nothing equivalent is "
            "forwarded to the external process"
        )),
        CAP_MCP: (NO, "no MCP config is forwarded to the external process"),
    },
}

#: Capabilities on the reference backend that exist only once a workspace
#: permissions policy is configured (see the reference-row reasons above).
#: Used by the degradation notice for the plain openai provider path,
#: which constructs OpenAIClient without a policy.
_PERMISSION_GATED_CAPS: tuple[str, ...] = (
    CAP_FILE_TOOLS, CAP_BASH, CAP_SUBAGENTS, CAP_TASK_TOOLS,
    CAP_UNDO_JOURNAL, CAP_WEB_TOOLS,
)

#: Aliases accepted by capability_gaps / degradation_notice in addition to
#: the create_client (backend, provider) vocabulary.
_BACKEND_ALIASES: dict[str, str] = {
    "api": BACKEND_ANTHROPIC_API,
    "anthropic": BACKEND_ANTHROPIC_API,
    "codex": BACKEND_CODEX_CLI,
    # Provider names accepted in the single-argument form: both factory
    # paths construct OpenAIClient with a permissions policy; the default
    # provider yields the CLI client.
    "kit": BACKEND_OPENAI,
    "ollama": BACKEND_OPENAI,
    "claude": BACKEND_CLI,
}


def canonical_backend(backend: str = "", provider: str = "") -> str:
    """Map create_client's (backend, provider) pair to a matrix key.

    Mirrors the dispatch order of ``create_client`` in api_client.py:
    provider kit/ollama -> OpenAIClient; provider openai -> CodexCLIClient
    for backend 'cli', else OpenAIClient; backend 'api' -> APIClient;
    everything else -> CLIClient. Canonical names and a few aliases are
    accepted directly. Never raises; unknown input falls back to the
    factory default ('cli').
    """
    b = str(backend or "").strip().lower()
    p = str(provider or "").strip().lower()
    if p in ("kit", "ollama"):
        return BACKEND_OPENAI
    if p == "openai":
        return BACKEND_CODEX_CLI if b == "cli" else BACKEND_OPENAI
    if b in _MATRIX:
        # Already canonical ('openai', 'cli', 'anthropic_api', 'codex_cli').
        # 'cli' is also the factory default, so this covers both readings.
        return b
    if b in _BACKEND_ALIASES:
        return _BACKEND_ALIASES[b]
    return BACKEND_CLI


def parity_matrix() -> dict:
    """The full backend/capability matrix as plain nested dicts.

    Shape: ``{backend: {capability: {"support", "reason", "label"}}}``
    with support one of ``"yes" | "conditional" | "no"``. The returned
    structure is a fresh copy — callers may mutate it freely.
    """
    return {
        backend: {
            cap: {
                "support": support,
                "reason": reason,
                "label": CAPABILITY_LABELS.get(cap, cap),
            }
            for cap, (support, reason) in row.items()
        }
        for backend, row in _MATRIX.items()
    }


def capability_gaps(backend: str) -> list[str]:
    """Human-readable names of capabilities MISSING on *backend*.

    Missing means support ``"no"`` while the reference backend (the
    OpenAI-compatible API client) has the capability at ``"yes"`` or
    ``"conditional"``. Conditional cells are NOT gaps — there the
    capability exists in a reduced or delegated form (see the matrix
    reasons). Accepts canonical names, create_client backend values and
    the aliases in ``_BACKEND_ALIASES``. Never raises.
    """
    canon = canonical_backend(backend)
    row = _MATRIX.get(canon, {})
    ref = _MATRIX.get(REFERENCE_BACKEND, {})
    gaps: list[str] = []
    for cap in CAPABILITY_LABELS:
        support = row.get(cap, (NO, ""))[0]
        ref_support = ref.get(cap, (NO, ""))[0]
        if support == NO and ref_support != NO:
            gaps.append(CAPABILITY_LABELS[cap])
    return gaps


def _delegated_labels(canon: str) -> list[str]:
    """Labels of capabilities that exist only in delegated/reduced form."""
    row = _MATRIX.get(canon, {})
    return [
        CAPABILITY_LABELS[cap]
        for cap in CAPABILITY_LABELS
        if row.get(cap, (NO, ""))[0] == CONDITIONAL
    ]


def degradation_notice(
    backend: str,
    provider: str = "",
    *,
    first_turn: bool = True,
    has_permissions: bool | None = None,
) -> str:
    """One-time session-start notice about missing backend capabilities.

    Returns a prompt block naming what this backend/config cannot do and
    how to get the full surface, or ``""`` when there is nothing to say
    (full surface) or when *first_turn* is False (the notice must appear
    exactly once, never per turn). Pure: the caller supplies backend,
    provider and — for the OpenAI-compatible backend — whether a
    workspace permissions policy is configured (``has_permissions``;
    None means infer from the provider: the kit/ollama factory paths
    wire one, the plain openai path does not). Never raises.
    """
    if not first_turn:
        return ""
    canon = canonical_backend(backend, provider)
    gaps = capability_gaps(canon)
    perm_gaps: list[str] = []
    if canon == BACKEND_OPENAI:
        no_perms = (
            has_permissions is False
            or (has_permissions is None
                and str(provider or "").strip().lower() == "openai")
        )
        if no_perms:
            perm_gaps = [CAPABILITY_LABELS[c] for c in _PERMISSION_GATED_CAPS]
    if not gaps and not perm_gaps:
        return ""

    label = _BACKEND_LABELS.get(canon, canon)
    lines: list[str] = [
        "# Backend capability notice (one-time)",
        f"This session runs on the {label} backend.",
    ]
    if gaps:
        lines.append(
            "Unavailable here (present on the full OpenAI-compatible API "
            "backend): " + "; ".join(gaps) + "."
        )
    if perm_gaps:
        lines.append(
            "No workspace permissions policy is configured on this "
            "client, so these are also unavailable: "
            + "; ".join(perm_gaps) + "."
        )
    delegated = _delegated_labels(canon)
    if delegated and canon != BACKEND_OPENAI:
        lines.append(
            "Delegated to the external process (DELFIN gating and "
            "journaling do not apply): " + "; ".join(delegated) + "."
        )
    lines.append(
        "Full surface: start the session on the OpenAI-compatible API "
        "backend with a workspace permissions policy (provider 'kit' or "
        "'ollama' in create_client)."
    )
    lines.append(
        "Tell the user about these limits once at the start of your next "
        "reply — do not repeat the notice afterwards, and if a request "
        "needs a missing capability, say plainly that this backend does "
        "not provide it."
    )
    return "\n".join(lines)


__all__ = [
    "CAPABILITY_LABELS",
    "REFERENCE_BACKEND",
    "canonical_backend",
    "capability_gaps",
    "degradation_notice",
    "parity_matrix",
]
