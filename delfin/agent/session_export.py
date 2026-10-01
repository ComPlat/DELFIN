"""Export a saved agent session as a replayable Jupyter notebook.

Provenance principle (J3): every result must trace back to the steps that
produced it, and those steps must be repeatable without the agent. A
session file carries the *why* (the chat); the tool trace carries the
*what* (each call with its arguments). This module joins the two into one
artefact: for every chemical step a Markdown cell (what and why) and a
code cell holding the native DELFIN call, so the notebook replays the
chemistry without the agent.

Internal agent tooling (tasks, memory, subagents, permissions) is left
out — a notebook documents chemistry, not the agent's own bookkeeping.
All exported text is scrubbed with ``memory_store._without_secrets``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

# Tool-name fragments that mark a chemical step worth replaying. The MCP
# surface prefixes chemistry tools with the server name; DELFIN's own
# importable API lives in delfin.api and delfin.smiles_converter.
_CHEMISTRY_TOOLS: dict[str, str] = {
    # trace tool name (suffix after last "__") -> native Python call
    "smiles_to_xyz": "from delfin.smiles_converter import smiles_to_xyz",
    "smiles_to_xyz_quick": "from delfin.smiles_converter import smiles_to_xyz_quick",
    "build_orca_input": "from delfin.stability_constant import build_orca_input",
    "submit_calculation": "from delfin import api",
    "qm_run": "from delfin import api",
    "pipeline_run": "from delfin import api",
    "run_orca_input": "from delfin import api",
    "extract_energy_table": "from delfin.api import extract_energy_table",
    "extract_imaginary_frequencies": "from delfin.api import extract_imaginary_frequencies",
    "extract_orbital_energies": "from delfin.api import extract_orbital_energies",
    "extract_excited_states": "from delfin.api import extract_excited_states",
    "extract_calc_summary_table": "from delfin import api",
    "extract_delfin_json": "from delfin import api",
    "extract_vibrational_modes": "from delfin import api",
    "extract_dipole": "from delfin import api",
    "extract_mulliken_charges": "from delfin import api",
    "extract_loewdin_charges": "from delfin import api",
    "extract_optimization_trajectory": "from delfin import api",
    "extract_scf_convergence": "from delfin import api",
    "extract_somf": "from delfin import api",
}

# Fragments of tool names that are the agent's own bookkeeping. Those
# steps cannot be replayed by a human and say nothing about the
# chemistry, so they never reach the notebook.
_EXCLUDED_FRAGMENTS = (
    "task",          # task_create / task_list / task_update
    "remember",      # memory writes
    "forget",
    "memory",
    "subagent",
    "permission",
    "skill",
    "session_export",
    "ask_user",
    "push_notification",
    "schedule",
    "cron",
    "search_docs",
    "search_calcs",
    "read_section",
    "list_docs",
    "list_sections",
    "list_tools",
    "describe_tool",
    "explain_delfin",
    "watch_job",
    "audit",
    "watchpost",
)


def _scrub(text: str) -> str:
    """Every exported text goes through the live secret redaction."""
    from .memory_store import _without_secrets
    return _without_secrets(text or "")


def _short_name(tool: str) -> str:
    """The suffix after the last ``__`` of an MCP-style tool name."""
    return str(tool or "").rsplit("__", 1)[-1]


def _is_chemistry(tool: str) -> bool:
    short = _short_name(tool)
    if short in _CHEMISTRY_TOOLS:
        return True
    if any(f in str(tool).lower() for f in _EXCLUDED_FRAGMENTS):
        return False
    # Other MCP chemistry tools (parsing/explainer surface) still count:
    # anything served by delfin-ops whose name starts with a known verb.
    return short.startswith(("extract_", "parse_", "find_", "check_"))


def _is_excluded(tool: str) -> bool:
    low = str(tool or "").lower()
    return any(f in low for f in _EXCLUDED_FRAGMENTS)


def _py_name(arg_key: str) -> str:
    return arg_key.replace("-", "_")


def _call_line(short: str, args: dict[str, Any]) -> str:
    """A best-effort native DELFIN call reconstructed from the arguments."""
    if not args:
        return f"# call replay not possible: no recorded arguments\n# original tool: {short}"
    # Map trace args onto the importable function; unknown keys are
    # passed through as keywords so a human sees exactly what was set.
    parts = ", ".join(
        f"{_py_name(k)}={v!r}" for k, v in list(args.items())[:8])
    return f"{short}({parts})"


def _cells_for_step(entry: dict, md) -> list:
    """Markdown (what and why) + code (native call) for one trace entry."""
    import nbformat  # local: only the export path needs it

    short = _short_name(entry.get("tool", ""))
    args = _args_of(entry)
    ok = bool(entry.get("ok", True))
    error = _scrub(str(entry.get("error") or ""))

    why = ""
    # best-effort: reuse nothing from the chat here; the pairing of
    # chat context to steps happens in export_session via _context_for.
    title = f"### Step: `{short}`"
    status = "succeeded" if ok else "**failed**"
    lines = [title, "", f"Tool call {status}." if ok
             else f"Tool call {status}: `{error}`"]
    cells = [md(_scrub("\n".join(lines)))]
    code_src = _scrub(_native_call(short, args))
    if not ok:
        code_src = (
            f"# This call failed in the session: {error}\n" + code_src)
    code = nbformat.v4.new_code_cell(code_src)
    return [cells[0], code]


def _native_call(short: str, args: dict[str, Any]) -> str:
    import_line = _CHEMISTRY_TOOLS.get(short)
    call = _call_line(short, args)
    if import_line:
        return f"{import_line}\n{call}"
    return f"from delfin import api  # native surface for `{short}`\n{call}"


def _args_of(entry: dict) -> dict:
    from . import tool_trace
    return tool_trace.call_args(entry)


def _chat_context(session: dict) -> list[str]:
    """The chat turns as Markdown blocks (what was asked, what was done)."""
    msgs = session.get("chat_messages") or session.get("engine_messages") or []
    blocks = []
    for m in msgs:
        role = str(m.get("role") or m.get("role_label") or "").strip()
        content = _scrub(str(m.get("content") or "").strip())
        if not content or role not in ("user", "assistant"):
            continue
        label = "User" if role == "user" else "Assistant"
        blocks.append(f"**{label}:** {content}")
    return blocks


def export_session(session: dict, *, trace_root: "Path | str" = ""):
    """Build a notebook (nbformat node) from a stored session + its trace.

    ``trace_root`` reads the tool trace from a directory other than the
    real ``~/.delfin`` — tests and archives use it; production passes
    nothing and gets this machine's own store.
    """
    import nbformat

    nb = nbformat.v4.new_notebook()
    md = nbformat.v4.new_markdown_cell

    session_id = str(session.get("session_id") or "unknown-session")
    title = _scrub(str(session.get("title") or session_id))
    workspace = _scrub(str(session.get("workspace") or ""))

    header = [
        f"# Session: {title}",
        "",
        f"Session ID: `{session_id}`",
    ]
    if workspace:
        header.append(f"Workspace: `{workspace}`")
    header += [
        "",
        "Every chemical step below is replayable without the agent: the",
        "Markdown says what and why, the code cell is the native DELFIN",
        "call with the arguments recorded in the session's tool trace.",
        "Text is scrubbed with `memory_store._without_secrets`.",
    ]
    cells = [md("\n".join(header))]

    # What was asked / answered, before the steps: the *why* of the work.
    for block in _chat_context(session):
        cells.append(md(block))

    from . import tool_trace
    entries = tool_trace.read(session_id, root=trace_root)
    n_steps = 0
    for entry in entries:
        tool = str(entry.get("tool") or "")
        if _is_excluded(tool) or not _is_chemistry(tool):
            continue
        cells.extend(_cells_for_step(entry, md))
        n_steps += 1

    if n_steps == 0 and entries:
        cells.append(md(
            "_No chemical tool calls found in this session's trace._"))
    nb.cells = cells
    return nb


def export_session_to_file(
    session: dict, out_path: "Path | str", *, trace_root: "Path | str" = "",
) -> Path:
    """Export a session and write the notebook to ``out_path``."""
    import nbformat

    nb = export_session(session, trace_root=trace_root)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, out)
    return out
