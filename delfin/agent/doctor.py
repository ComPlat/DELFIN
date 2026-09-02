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
    failed = data.get("failed_documents", []) if isinstance(data, dict) else []
    failed = [f for f in failed if isinstance(f, dict)]
    if not docs:
        detail = f"index at {_tilde(idx_path)} contains no documents"
        if failed:
            detail += f" ({len(failed)} could not be extracted)"
        return [_row(
            "doc index", WARN, detail,
            "re-run delfin-docs-index over the literature/ directory",
        )]
    rows = [_row(
        "doc index", PASS,
        f"{len(docs)} document(s) indexed at {_tilde(idx_path)}",
    )]
    if failed:
        # A file that yielded no text is not a document with nothing in
        # it — it is a file that was never indexed, and every search
        # against it returns nothing forever. Named here with the reason
        # because "indexed" was the only word the user got.
        names = ", ".join(
            f"{f.get('doc_id') or f.get('path', '?')}: {f.get('reason', '?')}"
            for f in failed[:4])
        rows.append(_row(
            "doc index", WARN,
            f"{len(failed)} document(s) yielded no text — {names}",
            "install pypdf for PDF text, or OCR a scanned manual before "
            "indexing it",
        ))
    stale = _stale_documents(data)
    if stale:
        rows.append(_row(
            "doc index", WARN,
            f"{len(stale)} indexed document(s) changed or vanished since the "
            f"index was built ({data.get('built_at', 'unknown time')}): "
            + ", ".join(stale[:4]),
            "re-run delfin-docs-index",
        ))
    return rows


def _stale_documents(index: dict) -> list[str]:
    """doc_ids whose source changed or disappeared after the index was built.

    ``built_at`` was written by both indexers and read by nobody: no
    consumer ever compared it to a source mtime, so a section deleted from
    the source last month still answered today, byte-identical to a fresh
    hit.
    """
    try:
        from delfin.doc_server.indexer import stale_documents
    except Exception:       # pragma: no cover - doc_server not importable
        return []
    try:
        return [s["doc_id"] for s in stale_documents(index)]
    except Exception:       # pragma: no cover
        return []


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
            "delfin-agent credentials set <NAME> "
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


def _uncontained_mcp(configs: dict) -> int:
    """How many configured servers no namespace is built around."""
    try:
        from .mcp_isolation import parse_isolation
        return sum(1 for cfg in configs.values()
                   if not cfg.get("url") and parse_isolation(cfg) is None)
    except Exception:
        return 0


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
    # Stated, not warned about: running a server uncontained is the default
    # and an ordinary choice. What is not ordinary is believing the shell's
    # sandbox covers it, and the doctor is where that belief gets checked.
    loose = _uncontained_mcp(configs)
    if loose:
        names += f" — {loose} without declared roots (outside the shell's isolation)"
    if ctx.get("fast", True):
        return [_row(
            "mcp servers", PASS,
            f"{len(configs)} configured: {names} (not probed; fast mode)",
        )]
    # Slow path: actually start each server and list its tools. The verdict
    # comes from ``unreachable_servers``, which reads ``last_error`` after
    # the call — this loop used to catch an exception that ``list_tools``
    # never raises (it fails closed and returns ``[]``), so a server with a
    # missing binary was reported as "configured + reachable".
    from .mcp_client import MCPRegistry, unreachable_servers
    reg = MCPRegistry()
    unreachable: list[str] = []
    try:
        reg.load(Path(workspace) if workspace else None)
        unreachable = unreachable_servers(reg)
    finally:
        try:
            reg.shutdown()
        except Exception:
            pass
    if unreachable:
        return [_row(
            "mcp servers", WARN,
            f"{len(configs)} configured, unreachable: "
            f"{'; '.join(sorted(unreachable))}",
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
            "delfin-agent scheduler start",
        )]
    return [_row(
        "scheduler daemon", PASS, "not running (no scheduled entries)",
    )]


def _check_attention(ctx: dict) -> list[dict]:
    """Attention inbox — what is waiting, in both directions.

    Two states, not one: events waiting for the user, and answers the
    user has already given that no session has picked up. The second was
    invisible everywhere (every surface filtered on ``pending``), so the
    user's answer could sit in the file while the report said all clear.
    """
    from .attention import _BLOCKING_KINDS, list_pending, list_undelivered

    out: list[dict] = []
    pending = list_pending()
    blocking = [ev for ev in pending if ev.get("kind") in _BLOCKING_KINDS]
    if blocking:
        out.append(_row(
            "attention inbox", WARN,
            f"{len(blocking)} event(s) blocking the agent"
            + (f", {len(pending) - len(blocking)} notice(s)"
               if len(pending) > len(blocking) else ""),
            "/attention in the dashboard Agent tab, then "
            "/attention answer <id> <text>",
        ))
    elif pending:
        out.append(_row(
            "attention inbox", PASS,
            f"{len(pending)} unread notice(s), nothing blocking the agent",
        ))
    else:
        out.append(_row("attention inbox", PASS, "no pending events"))

    waiting = list_undelivered()
    if waiting:
        out.append(_row(
            "attention answers", WARN,
            f"{len(waiting)} answered event(s) not yet delivered to a "
            "session — they reach the agent on its next turn",
            "start or continue an agent session in that workspace",
        ))
    return out


def _check_attention_transports(ctx: dict) -> list[dict]:
    """Can the agent reach the user out-of-band at all?

    Everything else in this report checks what the agent can DO. This
    checks whether anyone would be told when it stops and waits: with no
    desktop notifier, no webhook and no hook command, the fan-out for a
    blocking question is a no-op, and every other surface still reports
    healthy while an unattended run sits blocked until morning.
    """
    from .attention import transport_status

    st = transport_status(ctx.get("settings") or None)
    usable = list(st.get("usable") or [])
    detail = st.get("detail") or {}
    if usable:
        return [_row(
            "attention transports", PASS,
            f"{len(usable)} usable: {', '.join(usable)}",
        )]
    reasons = "; ".join(
        f"{name}: {detail[name]}" for name in ("desktop", "webhook", "hook")
        if detail.get(name))
    return [_row(
        "attention transports", WARN,
        "no usable transport — a blocked or finished run reaches you only "
        "in the inbox (" + (reasons or "nothing configured") + ")",
        "set agent.attention.notify_command to your own notifier, or "
        "agent.job_monitor.webhook_url to an https endpoint",
    )]


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


def _check_bash_isolation(ctx: dict) -> list[dict]:
    """Report whether shell commands run in a filesystem namespace.

    Write-target gating refuses paths outside the workspace, but it reads
    the command text — a subprocess started by that command is beyond it.
    Only namespace isolation contains that, and it is opt-in because it
    can disturb cluster workflows. Surfacing the state (and whether the
    machine even supports it) lets the user decide instead of assuming.
    """
    mode = "auto"
    try:
        from delfin.user_settings import load_settings
        mode = str(((load_settings() or {}).get("agent") or {})
                   .get("bash_isolation", "auto") or "auto").strip().lower()
    except Exception:
        pass
    try:
        from delfin.agent.api_client import _bwrap_functional
        usable = bool(_bwrap_functional())
    except Exception:
        usable = False

    fix = ("set agent.bash_isolation = \"bwrap\" to contain shell writes "
           "in every permission mode")
    if mode == "bwrap":
        if usable:
            return [{"check": "bash isolation", "status": "PASS",
                     "detail": "bwrap namespace active for every command"}]
        return [{"check": "bash isolation", "status": "FAIL",
                 "detail": "configured as bwrap but bwrap does not work here",
                 "fix": "install bwrap, or set agent.bash_isolation "
                        "= \"auto\" to fall back to the write gate"}]
    if mode == "off":
        return [{"check": "bash isolation", "status": "WARN",
                 "detail": "explicitly off — only the write-target gate "
                           "protects paths outside the workspace",
                 "fix": fix}]
    # "auto": isolated only in the unattended permission mode.
    detail = ("auto — isolated in bypassPermissions only"
              + ("" if usable else "; bwrap unusable here, so never isolated"))
    return [{"check": "bash isolation", "status": "WARN",
             "detail": detail, "fix": fix if usable else
             "install bwrap for real containment; the write-target gate "
             "stays active either way"}]


def _check_document_backends(ctx: dict) -> list[dict]:
    """Spreadsheet / PDF / Word support — one row per backend.

    A missing backend is a WARN, not a FAIL: the agent works without it,
    the affected tools are simply not advertised. It is reported because
    the alternative is discovering the gap halfway through a document
    task, when the tool the model was told to use is not there.
    """
    out: list[dict] = []
    for kind, module, dist, capability in (
        ("spreadsheets", "openpyxl", "openpyxl",
         "read_document / edit_sheet on .xlsx"),
        ("PDF", "pypdf", "pypdf", "read_document / fill_pdf_form"),
        ("Word", "docx", "python-docx", "read_document on .docx"),
        # A separate dependency from pypdf: taking PDF pages apart is not
        # the same library as laying text out on one.
        ("PDF writing", "reportlab", "reportlab", "create_pdf"),
        # LibreOffice formats. Reading only, which is why the capability
        # names reading: writing ODF is refused rather than approximated.
        ("OpenDocument", "odf", "odfpy",
         "read_document / compare_tables on .ods and .odt"),
    ):
        label = f"documents: {kind}"
        try:
            importlib.import_module(module)
        except Exception as exc:
            out.append(_row(
                label, WARN,
                f"{dist} not importable ({type(exc).__name__}) — "
                f"{capability} unavailable",
                "pip install 'delfin-complat[office]'",
            ))
        else:
            out.append(_row(label, PASS, f"{dist} available"))

    # OCR gets its own row because its failure mode is not a missing
    # import: pytesseract installs cleanly and then fails at the first
    # call because the program it drives is not on the machine. The row
    # names the component that is actually missing.
    try:
        from .office import ocr_availability
        status = ocr_availability()
    except Exception as exc:  # noqa: BLE001 — a probe may not crash the report
        out.append(_row(
            "documents: OCR", WARN,
            f"could not be determined ({type(exc).__name__})",
            "check the office module",
        ))
        return out
    if status["available"]:
        out.append(_row(
            "documents: OCR", PASS,
            f"{status['engine']} available — read_document(ocr=true) can "
            "read scanned pages"))
    else:
        out.append(_row(
            "documents: OCR", WARN,
            ("; ".join(status["detail"]) or "no OCR engine found")
            + " — scanned PDFs can be detected but not read",
            status["next_step"],
        ))
    return out


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
    ("attention transports", "_check_attention_transports"),
    ("benchmark truth", "_check_benchmark"),
    ("memory store", "_check_memory_store"),
    ("disk space", "_check_disk"),
    ("bash isolation", "_check_bash_isolation"),
    ("document backends", "_check_document_backends"),
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
