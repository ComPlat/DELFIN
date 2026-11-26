"""Which tools does the agent have but never call?

A read-only diagnostic over two real stores:

* ``delfin.agent.tool_trace`` — the JSONL of every tool call the agent
  made, per session, with ``ts`` / ``tool`` / ``ok``. This is the ONLY
  record of what was actually called (a call that failed is still a
  call; a refusal is a failure like any other).
* ``delfin.agent.api_client`` — the tool catalogue the client
  advertises (``_DOC_TOOLS_OPENAI``) filtered by ``advertisable_tools``
  with a default ``ToolSurfaceContext``. The comparison set comes from
  THAT catalogue, never from a list typed out here — a typed list would
  drift the moment a tool is added, and the point of the report is to
  measure the drift.

Finding baked into the design: the last time this was measured, 33 of
72 advertised tools had never been called in 2509 runs. That number is
not hardcoded anywhere; the report recomputes it so the drift is visible
each time it runs.

Usage: ``python -m delfin.agent.report_tools [--since-hours H]
[--trace-root DIR]``
"""

from __future__ import annotations

import argparse
import time
from typing import Any

from . import tool_trace
from .api_client import ToolSurfaceContext, _DOC_TOOLS_OPENAI, advertisable_tools


def _advertised_names(catalogue: list[dict[str, Any]] | None = None) -> list[str]:
    """Names of the tools the client actually advertises.

    Uses ``advertisable_tools`` with a default context (the surface a
    full session sees), not the raw catalogue: a tool held back from
    every session could never be called and is a catalogue question,
    not an agent-behaviour question.
    """
    cat = _DOC_TOOLS_OPENAI if catalogue is None else catalogue
    return [
        t.get("function", {}).get("name", "")
        for t in advertisable_tools(cat, ToolSurfaceContext())
        if t.get("function", {}).get("name")
    ]


def collect(
    *, since_s: float | None = None,
    root: Any = None, catalogue: list[dict[str, Any]] | None = None,
) -> dict:
    """Count tool calls in the trace store and diff against the advertised set.

    Returns a dict (never raises on missing data — an empty store yields
    an empty table, not an error):

    ``tools``: tool name -> {calls, failures, failure_share}, sorted by
    call count descending (most used first); names are the bare tool
    name, with any ``server__`` MCP prefix stripped so a traced
    ``mcp__delfin-ops__list_tools`` compares against the advertised
    ``list_tools``.
    ``never_used``: advertised tools with zero calls, sorted.
    ``advertised_only_calls``: tools called but never advertised — a
    back-compat check, reported rather than hidden.
    ``n_advertised`` / ``n_used_advertised`` / ``n_never_used`` /
    ``total_calls`` / ``sessions``: the summary counts.
    ``window_s``: the window actually used (None = all history).
    """
    now = time.time()
    names: list[str] = tool_trace.sessions(root=root) if root else tool_trace.sessions()
    counts: dict[str, dict[str, Any]] = {}
    total = 0
    for session in names:
        for entry in tool_trace.read(session, root=root or ""):
            if not isinstance(entry, dict):
                continue
            if since_s is not None:
                try:
                    if float(entry.get("ts") or 0) < now - float(since_s):
                        continue
                except (TypeError, ValueError):
                    continue
            name = str(entry.get("tool") or "")
            if not name:
                continue
            # mcp__server__tool -> tool, so traced names compare against
            # the flat advertised surface.
            base = name.rsplit("__", 1)[-1] if name.startswith("mcp__") else name
            row = counts.setdefault(base, {"calls": 0, "failures": 0})
            row["calls"] += 1
            if entry.get("ok") is False:
                row["failures"] += 1
            total += 1
    for row in counts.values():
        calls = row["calls"]
        row["failure_share"] = (row["failures"] / calls) if calls else 0.0

    ordered = dict(sorted(counts.items(), key=lambda kv: (-kv[1]["calls"], kv[0])))
    advertised = set(_advertised_names(catalogue))
    never_used = sorted(advertised - set(ordered))
    called_not_advertised = sorted(set(ordered) - advertised)
    return {
        "tools": ordered,
        "never_used": never_used,
        "called_not_advertised": called_not_advertised,
        "n_advertised": len(advertised),
        "n_used_advertised": len(advertised & set(ordered)),
        "n_never_used": len(never_used),
        "total_calls": total,
        "sessions": len(names),
        "window_s": since_s,
    }


def format_text(data: dict) -> str:
    """Render collect()'s output: usage table first, then the never-used list."""
    tools = data.get("tools") or {}
    lines: list[str] = []
    window = data.get("window_s")
    scope = (f"last {window / 3600:g}h" if window else "all history")
    lines.append(
        f"Tool usage across {data.get('sessions', 0)} traced session(s), {scope}:"
    )
    if not tools:
        lines.append("No tool calls recorded.")
    else:
        header = f"{'tool':<36}{'calls':>8}{'fails':>8}{'fail%':>8}"
        lines.append(header)
        lines.append("-" * len(header))
        for name, row in tools.items():
            lines.append(
                f"{name[:36]:<36}{row['calls']:>8}{row['failures']:>8}"
                f"{100 * row['failure_share']:>7.1f}%"
            )
        lines.append("-" * len(header))
    lines.append(f"Total calls: {data.get('total_calls', 0)}")
    lines.append("")
    n_adv = data.get("n_advertised", 0)
    n_never = data.get("n_never_used", 0)
    lines.append(
        f"Advertised tools never called: {n_never} of {n_adv}"
        f" ({n_adv and 100 * n_never / n_adv:.1f}%)"
    )
    for name in data.get("never_used") or []:
        lines.append(f"  - {name}")
    stray = data.get("called_not_advertised") or []
    if stray:
        lines.append("")
        lines.append(
            "Called but not in the advertised catalogue (stale trace or "
            "renamed tool): " + ", ".join(stray)
        )
    return "\n".join(lines)


def format_roles(report: dict[str, dict[str, Any]]) -> str:
    """Render ``role_tool_surface_report``: what each role costs to offer.

    The listing above says what a session is offered. This says what that
    costs, and what the role scoping buys: every row is compared with the
    unrestricted baseline ``""``, so "office_agent saves N tools and M
    schema tokens" is a measurement rather than an assertion.
    """
    base = report.get("", {})
    base_count = int(base.get("count", 0) or 0)
    base_tokens = int(base.get("total_tokens", 0) or 0)
    header = f"{'role':<22}{'tools':>7}{'tokens':>9}{'saves':>8}{'tok saved':>11}"
    lines = [
        "Tool surface per role, with the unrestricted baseline first:",
        header,
        "-" * len(header),
    ]
    # The baseline first, then the scoped roles cheapest-surface first: the
    # question a reader has is "which role offers least", not alphabetical.
    ordered = [""] + sorted(
        (r for r in report if r != ""),
        key=lambda r: (int(report[r].get("total_tokens", 0) or 0), r),
    )
    for role in ordered:
        row = report[role]
        count = int(row.get("count", 0) or 0)
        tokens = int(row.get("total_tokens", 0) or 0)
        if role == "":
            lines.append(
                f"{'(unrestricted)':<22}{count:>7}{tokens:>9}{'--':>8}{'--':>11}")
            continue
        lines.append(
            f"{role:<22}{count:>7}{tokens:>9}"
            f"{base_count - count:>8}{base_tokens - tokens:>11}")
    return "\n".join(lines)


def collect_roles(
    roles: list[str] | None = None,
) -> dict[str, dict[str, Any]]:
    """The per-role advertised surface and its schema-token cost.

    A thin, documented door onto ``api_client.role_tool_surface_report`` so
    this module's entry point reaches it by CALL rather than by import: the
    report is the product's own answer to "what does role scoping save",
    and until now nothing a person could run printed it.

    ``defer=False`` is passed explicitly because the question here is what
    each role can CALL, not what ships on a request; see that function's
    docstring for the measured reason the default had to stay put.

    THE ROLE LIST IS BROADER THAN THE REPORT'S DEFAULT ONE. With ``roles``
    left as ``None`` the report enumerates only ``_ROLE_EXEC_ALLOWLIST`` and
    ``_ROLE_EXEC_DENYLIST`` -- 7 roles. But ``_role_allowlist``, which is
    what actually refuses a tool at the execution call site, also answers
    for every subagent preset carrying a ``tools:`` line: 17 more, so 19
    restricted roles in all. A table built from the report's default
    therefore showed 8 of them and left a reader to conclude the other 11
    were unrestricted. Naming them here rather than widening the report's
    default keeps ``test_the_office_role_appears_in_the_surface_report``
    and every other existing caller measuring what they measured before.
    """
    from .api_client import (_ROLE_EXEC_ALLOWLIST, _ROLE_EXEC_DENYLIST,
                             _role_allowlist, role_tool_surface_report)
    from .subagents import SUBAGENT_PRESETS

    if roles is None:
        roles = ["", *sorted(
            set(_ROLE_EXEC_ALLOWLIST) | set(_ROLE_EXEC_DENYLIST)
            | {r for r in SUBAGENT_PRESETS if _role_allowlist(r) is not None})]
    return role_tool_surface_report(roles, defer=False)



def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Which tools does the agent never use?")
    parser.add_argument(
        "--since-hours", type=float, default=None,
        help="Only count calls newer than this many hours (default: all)")
    parser.add_argument(
        "--trace-root", default=None,
        help="Read traces from this directory instead of the local store")
    parser.add_argument(
        "--roles", action="store_true",
        help="Print the advertised tool surface and schema-token cost per "
             "role, against the unrestricted baseline, instead of the usage "
             "table")
    args = parser.parse_args(argv)
    if args.roles:
        print(format_roles(collect_roles()))
        return 0
    data = collect(
        since_s=(args.since_hours * 3600.0
                 if args.since_hours is not None else None),
        root=args.trace_root,
    )
    print(format_text(data))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
