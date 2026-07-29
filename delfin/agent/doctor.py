"""One-surface environment health report for the DELFIN agent stack.

Every prerequisite the agent depends on already has *some* check buried
in its own module — the doc-index gate in the API client, the masked
credential listing, ``shutil.which`` probes in the tool adapters, the
scheduler pid file, the optimize-check gate.  Users hit each of them
only when a task fails on it.  ``run_doctor`` aggregates those existing
probes (reusing the real module logic, not re-implementing it) into a
single report the CLI and the dashboard can render before anything
breaks.

Contract:

- ``run_doctor`` NEVER raises.  A probe that explodes becomes a FAIL
  row with the exception message; the report is always complete.
- Credential checks report only *which* providers are configured —
  no value (masked or otherwise) ever appears in a detail line.
- ``fast=True`` (the default) touches no network and starts no MCP
  server process: it only lists what is configured.
"""

from __future__ import annotations

import importlib
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Callable

PASS = "PASS"
WARN = "WARN"
FAIL = "FAIL"

_ICONS = {PASS: "✅", WARN: "⚠️", FAIL: "❌"}

# Provider label → credential name consumed by the engine (see
# ``credentials._WELL_KNOWN_KEYS``).  Labels are what the report shows;
# values are looked up but NEVER printed.
_PROVIDER_KEYS: tuple[tuple[str, str], ...] = (
    ("KIT", "KIT_TOOLBOX_API_KEY"),
    ("Anthropic", "ANTHROPIC_API_KEY"),
    ("OpenAI", "OPENAI_API_KEY"),
)

_CHEM_BINARIES: tuple[str, ...] = ("xtb", "orca")
_PYTHON_DEPS: tuple[str, ...] = ("rdkit", "openbabel")

_DISK_WARN_GB = 1.0


def _row(check: str, status: str, detail: str, fix: str = "") -> dict:
    return {"check": check, "status": status, "detail": detail, "fix": fix}


def _tilde(path: Path | str) -> str:
    return str(path).replace(str(Path.home()), "~")


# ---------------------------------------------------------------------------
# Individual checks — each returns a list of result rows.  Underlying
# logic is REUSED from the module that owns it; the doctor only decides
# PASS/WARN/FAIL and phrases the fix.
# ---------------------------------------------------------------------------


def _check_doc_index(ctx: dict) -> list[dict]:
    """Manual/doc search index — same path the search_docs gate loads."""
    import json

    from delfin.doc_server.indexer import get_default_index_path

    idx_path = get_default_index_path()
    if not idx_path.exists():
        return [_row(
            "doc index", WARN,
            f"no index at {_tilde(idx_path)} — search_docs unavailable",
            "run delfin-docs-index to build the manual search index",
        )]
    try:
        data = json.loads(idx_path.read_text(encoding="utf-8"))
    except Exception as exc:
        return [_row(
            "doc index", FAIL,
            f"index at {_tilde(idx_path)} is unreadable: {exc}",
            "rebuild it with delfin-docs-index",
        )]
    if isinstance(data, list) and data:
        data = data[0]
    docs = data.get("documents", {}) if isinstance(data, dict) else {}
    if not docs:
        return [_row(
            "doc index", WARN,
            f"index at {_tilde(idx_path)} contains no documents",
            "re-run delfin-docs-index over the literature/ directory",
        )]
    return [_row(
        "doc index", PASS,
        f"{len(docs)} document(s) indexed at {_tilde(idx_path)}",
    )]


def _check_credentials(ctx: dict) -> list[dict]:
    """Which providers have keys — presence only, values never shown."""
    from .credentials import load_credential

    have = [label for label, name in _PROVIDER_KEYS
            if bool(load_credential(name))]
    missing = [label for label, _ in _PROVIDER_KEYS if label not in have]
    if not have:
        return [_row(
            "credentials", FAIL,
            "no provider keys configured "
            f"(missing: {', '.join(missing)})",
            "python -m delfin.agent.cli credentials set <NAME> "
            "(e.g. ANTHROPIC_API_KEY)",
        )]
    detail = f"set: {', '.join(have)}"
    if missing:
        detail += f" | missing: {', '.join(missing)}"
    return [_row("credentials", PASS, detail)]


def _check_binaries(ctx: dict) -> list[dict]:
    """Chemistry binaries on PATH — reuses the tools install-policy table."""
    out: list[dict] = []
    for name in _CHEM_BINARIES:
        found = shutil.which(name)
        if found:
            out.append(_row(f"binary: {name}", PASS, f"found at {found}"))
            continue
        fix = f"install {name} and add it to PATH"
        try:
            from delfin.tools._environment import tool_info
            info = tool_info(name, "binary")
            hint = (info.install_hint or "").split(". ")[0].strip()
            if hint:
                fix = hint
            if info.source:
                fix += f" ({info.source})"
        except Exception:
            pass
        out.append(_row(f"binary: {name}", WARN, "not found on PATH", fix))
    return out


def _check_python_deps(ctx: dict) -> list[dict]:
    """Optional chemistry Python deps — import probe, importorskip style."""
    out: list[dict] = []
    for mod in _PYTHON_DEPS:
        try:
            importlib.import_module(mod)
        except Exception as exc:
            fix = f"pip install {mod}  (or via conda)"
            try:
                from delfin.tools._environment import tool_info
                hint = tool_info(mod, "python").install_hint
                if hint:
                    fix = hint
            except Exception:
                pass
            out.append(_row(
                f"python: {mod}", WARN,
                f"not importable ({type(exc).__name__})", fix,
            ))
        else:
            out.append(_row(f"python: {mod}", PASS, "importable"))
    return out


def _check_mcp(ctx: dict) -> list[dict]:
    """MCP servers — configured list; reachability only when fast=False."""
    from .mcp_client import _load_configs

    workspace = ctx.get("workspace")
    configs = _load_configs(Path(workspace) if workspace else None)
    if not configs:
        return [_row(
            "mcp servers", WARN, "no MCP servers configured",
            "add servers to ~/.delfin/mcp_servers.json "
            "(builtin delfin-tools was disabled)",
        )]
    names = ", ".join(sorted(configs))
    if ctx.get("fast", True):
        return [_row(
            "mcp servers", PASS,
            f"{len(configs)} configured: {names} (not probed; fast mode)",
        )]
    # Slow path: actually start each server and list its tools.
    from .mcp_client import MCPRegistry
    reg = MCPRegistry()
    unreachable: list[str] = []
    try:
        reg.load(Path(workspace) if workspace else None)
        for name, srv in reg.servers.items():
            try:
                srv.list_tools()
            except Exception:
                unreachable.append(name)
    finally:
        try:
            reg.shutdown()
        except Exception:
            pass
    if unreachable:
        return [_row(
            "mcp servers", WARN,
            f"{len(configs)} configured, unreachable: "
            f"{', '.join(sorted(unreachable))}",
            "check the server command/URL in mcp_servers.json",
        )]
    return [_row(
        "mcp servers", PASS, f"{len(configs)} configured + reachable: {names}",
    )]


def _check_scheduler(ctx: dict) -> list[dict]:
    """Scheduler daemon pid file + liveness — reuses daemon_status."""
    from .scheduler_daemon import daemon_status

    st = daemon_status()
    if st.get("running"):
        return [_row(
            "scheduler daemon", PASS,
            f"running (PID {st.get('pid')}), "
            f"{st.get('entries', 0)} scheduled entries",
        )]
    if st.get("entries", 0) > 0:
        return [_row(
            "scheduler daemon", WARN,
            f"not running but {st['entries']} scheduled entries exist "
            "— they will not fire",
            "python -m delfin.agent.cli scheduler start",
        )]
    return [_row(
        "scheduler daemon", PASS, "not running (no scheduled entries)",
    )]


def _check_attention(ctx: dict) -> list[dict]:
    """Attention inbox — pending events awaiting a user answer."""
    from .attention import list_pending

    pending = list_pending()
    if pending:
        return [_row(
            "attention inbox", WARN,
            f"{len(pending)} pending event(s) awaiting your answer",
            "open the dashboard Agent tab and answer the pending events",
        )]
    return [_row("attention inbox", PASS, "no pending events")]


def _check_benchmark(ctx: dict) -> list[dict]:
    """Benchmark tasks + ground truth — one-line optimize_check summary."""
    from .optimize_check import run_checks

    issues = run_checks()
    errors = [i for i in issues if i.severity == "error"]
    warns = [i for i in issues if i.severity == "warn"]
    fix = "python -m delfin.agent.optimize_check for the full list"
    if errors:
        return [_row(
            "benchmark truth", FAIL,
            f"{len(errors)} error(s), {len(warns)} warning(s) in "
            "benchmark tasks / ground truth", fix,
        )]
    if warns:
        return [_row(
            "benchmark truth", WARN,
            f"{len(warns)} warning(s) in benchmark tasks / ground truth",
            fix,
        )]
    return [_row(
        "benchmark truth", PASS, "benchmark tasks + ground truth OK",
    )]


def _check_memory_store(ctx: dict) -> list[dict]:
    """~/.delfin writable — memory/credential/session stores live there."""
    delfin_dir = Path.home() / ".delfin"
    if not delfin_dir.exists():
        if os.access(Path.home(), os.W_OK):
            return [_row(
                "memory store", PASS,
                f"{_tilde(delfin_dir)} will be created on first use",
            )]
        return [_row(
            "memory store", FAIL,
            f"{_tilde(Path.home())} is not writable — cannot create "
            f"{_tilde(delfin_dir)}",
            "fix the home-directory permissions",
        )]
    probe = delfin_dir / ".doctor_probe"
    try:
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()
    except OSError as exc:
        return [_row(
            "memory store", FAIL,
            f"{_tilde(delfin_dir)} is not writable: {exc}",
            f"chmod u+rwx {_tilde(delfin_dir)}",
        )]
    return [_row(
        "memory store", PASS, f"{_tilde(delfin_dir)} writable",
    )]


def _check_disk(ctx: dict) -> list[dict]:
    """Free disk space in the workspace — calculations need headroom."""
    workspace = Path(ctx.get("workspace") or Path.cwd())
    usage = shutil.disk_usage(workspace)
    free_gb = usage.free / (1024 ** 3)
    if free_gb < _DISK_WARN_GB:
        return [_row(
            "disk space", WARN,
            f"{free_gb:.2f} GB free at {_tilde(workspace)} "
            f"(< {_DISK_WARN_GB:g} GB)",
            "free up disk space before starting calculations",
        )]
    return [_row(
        "disk space", PASS,
        f"{free_gb:.1f} GB free at {_tilde(workspace)}",
    )]


# Names are looked up on the module at run time so tests can monkeypatch
# a single ``_check_*`` function without rebuilding this table.
_CHECK_ATTRS: tuple[tuple[str, str], ...] = (
    ("doc index", "_check_doc_index"),
    ("credentials", "_check_credentials"),
    ("chemistry binaries", "_check_binaries"),
    ("python deps", "_check_python_deps"),
    ("mcp servers", "_check_mcp"),
    ("scheduler daemon", "_check_scheduler"),
    ("attention inbox", "_check_attention"),
    ("benchmark truth", "_check_benchmark"),
    ("memory store", "_check_memory_store"),
    ("disk space", "_check_disk"),
)


def _normalise(row: Any, group: str) -> dict:
    """Coerce whatever a (possibly monkeypatched) check returned."""
    if not isinstance(row, dict):
        return _row(group, FAIL, f"check returned {type(row).__name__}")
    status = row.get("status", FAIL)
    if status not in (PASS, WARN, FAIL):
        status = FAIL
    return _row(
        str(row.get("check", group)), status,
        str(row.get("detail", "")), str(row.get("fix", "")),
    )


def run_doctor(
    workspace: str | Path | None = None,
    *,
    settings: dict | None = None,
    fast: bool = True,
) -> list[dict]:
    """Run every health check; never raises.

    Returns a list of ``{check, status: PASS|WARN|FAIL, detail, fix}``
    rows in a stable order.  A probe that raises becomes a single FAIL
    row for its group so one broken subsystem cannot hide the rest of
    the report.
    """
    ctx: dict = {
        "workspace": str(workspace) if workspace else "",
        "settings": settings or {},
        "fast": bool(fast),
    }
    module = sys.modules[__name__]
    results: list[dict] = []
    for group, attr in _CHECK_ATTRS:
        fn: Callable[[dict], list[dict]] | None = getattr(
            module, attr, None)
        if fn is None:
            results.append(_row(group, FAIL, f"check {attr} missing"))
            continue
        try:
            rows = fn(ctx)
        except Exception as exc:  # noqa: BLE001 — report, never raise
            results.append(_row(
                group, FAIL,
                f"check crashed: {type(exc).__name__}: {exc}",
            ))
            continue
        if not isinstance(rows, list):
            rows = [rows]
        if not rows:
            results.append(_row(group, FAIL, "check returned no result"))
            continue
        results.extend(_normalise(r, group) for r in rows)
    return results


def format_doctor(results: list[dict]) -> str:
    """Aligned per-check report + fix hints + one summary line."""
    if not results:
        return "No checks ran.\n0 pass, 0 warn, 0 fail"
    width = max(len(str(r.get("check", ""))) for r in results)
    lines: list[str] = []
    counts = {PASS: 0, WARN: 0, FAIL: 0}
    for r in results:
        status = r.get("status", FAIL)
        counts[status if status in counts else FAIL] += 1
        icon = _ICONS.get(status, "❌")
        name = str(r.get("check", ""))
        lines.append(f"{icon} {name:<{width}}  {r.get('detail', '')}")
        fix = str(r.get("fix", "") or "").strip()
        if status != PASS and fix:
            lines.append(f"   {'':<{width}}  fix: {fix}")
    lines.append("")
    lines.append(
        f"{counts[PASS]} pass, {counts[WARN]} warn, {counts[FAIL]} fail"
    )
    return "\n".join(lines)


__all__ = ["run_doctor", "format_doctor", "PASS", "WARN", "FAIL"]
