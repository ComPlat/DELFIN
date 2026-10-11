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
import re
import shutil
import subprocess
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


def _row(check: str, status: str, detail: str, fix: str = "",
         *, command: str = "", setting: tuple | None = None) -> dict:
    """One report row.

    ``fix`` is prose for a person. ``command`` and ``setting`` are the
    machine-actionable form, and a check that has one DECLARES it rather
    than leaving it to be parsed back out of the prose: a remedy read out
    of a sentence is a remedy a wording change breaks, and this one is
    offered to the user for approval and then executed.

    Most prerequisites have no such form -- installing a system package,
    logging a credential helper in, reordering a library path. Those
    carry ``fix`` alone, and saying so is part of the answer: the agent
    must not improvise around them.
    """
    row = {"check": check, "status": status, "detail": detail, "fix": fix}
    if command:
        row["command"] = command
    if setting:
        row["setting"] = list(setting)
    return row


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
    """Chemistry binaries — found the way DELFIN itself finds them.

    DELFIN does not require xtb/orca on PATH: it also resolves them from
    the qm_tools directories (``~/.delfin/qm_tools/bin``), ``*_BINARY``
    environment variables and known system locations
    (``qm_runtime.resolve_tool``). A bare ``shutil.which`` here reported
    every qm_tools install as "not found on PATH" while DELFIN happily
    ran the same binary — so this check reuses the real resolver instead
    of re-implementing a PATH-only subset of it.
    """
    out: list[dict] = []
    for name in _CHEM_BINARIES:
        found = None
        try:
            from delfin import qm_runtime
            found = qm_runtime.find_tool_executable(name)
        except Exception:
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


def _check_test_runner(ctx: dict) -> list[dict]:
    """Can the agent run the suite in the interpreter it would use.

    The one prerequisite with a remedy DELFIN can carry out, and the
    reason the proposal mechanism exists: pytest was declared only in an
    extra the default install does not select, `python -m pytest` with no
    pytest writes no report, and the agent read "no report file produced"
    as something to work around. The field reports show what it did --
    a venv of its own, or a wrapper script in the home directory.

    Asked of THIS interpreter, which is the one `run_tests` launches.
    """
    missing = importlib.util.find_spec("pytest") is None
    if not missing:
        return [_row("test runner", PASS, f"pytest available to {sys.executable}")]
    return [_row(
        "test runner", WARN,
        f"pytest is not installed in {sys.executable}",
        "install the test extra; until then the agent cannot run the suite "
        "and must say so rather than build a runner of its own",
        command=f"{sys.executable} -m pip install 'delfin-complat[test]'",
    )]


def _uncontained_mcp(configs: dict, workspace=None) -> int:
    """How many configured servers no namespace is built around.

    Asks the registry's own decision rather than re-deriving it: a
    built-in whose roots come from the settings is contained, and a
    doctor that counted it as loose would be reporting on a rule it had
    reimplemented.
    """
    try:
        from .mcp_client import _isolation_for
        return sum(1 for name, cfg in configs.items()
                   if not cfg.get("url")
                   and _isolation_for(name, cfg, workspace) is None)
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
    loose = _uncontained_mcp(configs, Path(workspace) if workspace else None)
    # The switch that contains DELFIN's own servers exists
    # (agent.mcp_isolation = "builtin"), opt-in by its author's decision
    # until it has run against a real session. A row that only says
    # "without declared roots" leaves the user to find the switch; this
    # names it, as a proposal `/fix` can apply with approval -- and only
    # where the loose servers are the built-ins the switch covers.
    offer_fix = ""
    offer_setting = None
    if loose:
        names += f" — {loose} without declared roots (outside the shell's isolation)"
        try:
            from .mcp_isolation import builtin_isolation_enabled
            from .mcp_client import _BUILTIN_SERVERS
            # Declared roots live at the top of the entry ("roots" /
            # "read_roots"), and "isolation": "off" is the escape hatch;
            # a built-in with neither is what the switch would contain.
            covered = any(n in _BUILTIN_SERVERS and not cfg.get("url")
                          and not cfg.get("roots") and not cfg.get("read_roots")
                          and str(cfg.get("isolation", "") or "").lower() != "off"
                          for n, cfg in configs.items())
            if covered and not builtin_isolation_enabled():
                offer_fix = ("set agent.mcp_isolation = \"builtin\" to run "
                             "DELFIN's own servers inside derived roots (the "
                             "workspace, the office folder and its state "
                             "directories); a third-party server needs an "
                             "explicit \"isolation\" entry in its config")
                offer_setting = ("agent.mcp_isolation", "builtin")
        except Exception:
            pass
    if ctx.get("fast", True):
        return [_row(
            "mcp servers", PASS,
            f"{len(configs)} configured: {names} (not probed; fast mode)",
            offer_fix, setting=offer_setting,
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


def _check_git_tooling(ctx: dict) -> list[dict]:
    """git and gh: present, authenticated, and able to sign a commit.

    Four questions, four rows, because they fail independently and the
    remedy differs for each. None of them was asked anywhere in the
    product. git's absence was noticed only as a side effect of
    ``_check_push`` -- reported under the name "git remote", with no fix
    string, and only on a checkout that has an origin remote. gh was never
    probed at all, while the write gate REQUIRES the pull-request route
    (a contributor's push to the default branch is refused and the agent
    is told to open a PR), so the one tool the sanctioned path needs was
    the one nothing checked. An absent or logged-out gh surfaced as bash
    exit 127 or gh's own stderr, with nothing mapping it to an action.

    Read-only throughout: ``--version``, ``config --get`` and ``gh auth
    status`` publish nothing. WARN rather than FAIL, the module's
    convention for a missing prerequisite, and every row carries a fix --
    a row without one only tells the reader that something is wrong.
    """
    import shutil as _sh

    cwd = str(ctx.get("workspace") or ".")
    rows: list[dict] = []

    git = _sh.which("git")
    if not git:
        # Identity cannot be asked about without the binary, and gh is a
        # separate tool: name this one and go on.
        rows.append(_row(
            "git installed", WARN, "git is not on PATH",
            "install git; the agent cannot branch, commit or push without it"))
    else:
        rows.append(_row("git installed", PASS, git))
        missing = []
        for field in ("user.name", "user.email"):
            try:
                done = subprocess.run(["git", "config", "--get", field],
                                      capture_output=True, text=True,
                                      timeout=10, cwd=cwd)
                if done.returncode != 0 or not (done.stdout or "").strip():
                    missing.append(field)
            except (OSError, subprocess.SubprocessError):
                missing.append(field)
        if missing:
            rows.append(_row(
                "git identity", WARN,
                "not configured: " + ", ".join(missing),
                "git config --global user.name \"...\" and "
                "git config --global user.email \"...\" -- a commit cannot "
                "be made without them, and the failure arrives at commit "
                "time rather than now"))
        else:
            rows.append(_row("git identity", PASS, "user.name and user.email set"))

    if ctx.get("skip_gh"):
        # The push gate asks for git alone: a plain ``git push`` needs no
        # gh, and ``gh auth status`` validates the token over the network
        # (up to 20 s) for rows the caller would discard.
        return rows

    gh = _sh.which("gh")
    if not gh:
        rows.append(_row(
            "gh installed", WARN, "the GitHub CLI is not on PATH",
            "install it (https://cli.github.com) -- the write gate routes "
            "changes to the default branch through a pull request, and that "
            "is the tool for it; without gh the agent can only hand the user "
            "a compare URL"))
        return rows

    rows.append(_row("gh installed", PASS, gh))
    try:
        done = subprocess.run(["gh", "auth", "status"],
                              capture_output=True, text=True, timeout=20,
                              cwd=cwd)
    except subprocess.TimeoutExpired:
        rows.append(_row(
            "gh authenticated", WARN, "gh auth status did not answer in 20 s",
            "check network access to github from this host"))
        return rows
    except (OSError, subprocess.SubprocessError) as exc:
        rows.append(_row("gh authenticated", WARN, f"could not ask gh: {exc}",
                         "run 'gh auth status' by hand"))
        return rows
    if done.returncode == 0:
        rows.append(_row("gh authenticated", PASS, "a login is active"))
    else:
        # Presence and authentication are different questions, and gh
        # answers the second one only when asked.
        rows.append(_row(
            "gh authenticated", WARN, "no active GitHub login",
            "gh auth login -- pull requests cannot be opened until this "
            "succeeds"))
    return rows


def _remote_url(workspace: str, remote: str) -> str | None:
    """The push URL ``remote`` names in ``workspace``, or None.

    A configured remote name resolves through ``git remote get-url
    --push`` (pushurl and insteadOf applied; reads configuration, runs
    nothing). Anything else is taken as the URL itself, the way
    ``git push <url>`` takes it; a relative local path is made absolute
    against ``workspace`` because the probe runs elsewhere. None: the
    word is neither a configured remote nor a URL or path.
    """
    try:
        done = subprocess.run(["git", "remote", "get-url", "--push", remote],
                              capture_output=True, text=True, timeout=10,
                              cwd=workspace, stdin=subprocess.DEVNULL)
        if done.returncode == 0 and done.stdout.strip():
            remote = done.stdout.strip()
        elif not any(c in remote for c in ":/\\") and remote != ".":
            # A bare word that is not a configured remote: git would fail
            # with "does not appear to be a git repository".
            return None
    except (OSError, subprocess.SubprocessError):
        pass
    if "://" in remote:
        return remote
    colon, slash = remote.find(":"), remote.find("/")
    if colon > 0 and (slash < 0 or colon < slash):
        return remote  # scp-like host:path
    path = Path(remote).expanduser()
    return str(path if path.is_absolute() else (Path(workspace) / path))


def _check_push(ctx: dict) -> list[dict]:
    """Can this host reach the git remote at all.

    Read-only: ``git ls-remote`` asks the remote what refs it has and
    publishes nothing. It is the same question a push answers the hard
    way, and answering it here costs one call instead of thirteen --
    measured 2026-09-28, when two sessions spent eleven and thirteen
    attempts varying transports and proxies against a host that could
    not reach github at all. Nothing told them the first message was
    final, so they kept rephrasing it.

    This does not push and does not make pushing easier. It says, before
    the work, whether the last step will be the user's.
    """
    import tempfile

    from .push_diagnosis import diagnose

    workspace = str(ctx.get("workspace") or ".")
    remote = str(ctx.get("remote") or "origin")
    url = _remote_url(workspace, remote)
    if url is None:
        return [_row("git remote", WARN,
                     f"this checkout has no remote named '{remote}'",
                     f"git remote add {remote} <url>, or push to a remote "
                     "that exists (git remote -v lists them)")]
    # The probe runs in this process, outside any sandbox the agent's own
    # commands run in, and the checkout's .git/config is a file the agent
    # can edit: remote.<name>.uploadpack, core.sshCommand and a local
    # credential.helper are commands git would run for ls-remote (a
    # probe in the checkout executed an uploadpack written there). So the
    # URL is read from the checkout and asked from an empty directory
    # with discovery stopped above it -- only the user's global and
    # system configuration apply. No controlling terminal and no prompt:
    # a passphrase or password question becomes a failure, not a hang
    # that takes the user's terminal.
    env = {k: v for k, v in os.environ.items()
           if k not in ("GIT_DIR", "GIT_WORK_TREE", "GIT_COMMON_DIR",
                        "GIT_INDEX_FILE", "GIT_CONFIG_PARAMETERS",
                        "GIT_CONFIG_COUNT")}
    env["GIT_TERMINAL_PROMPT"] = "0"
    try:
        with tempfile.TemporaryDirectory(prefix="delfin-ls-remote-") as away:
            env["GIT_CEILING_DIRECTORIES"] = str(Path(away).parent)
            done = subprocess.run(
                ["git", "-c", "protocol.ext.allow=never",
                 "ls-remote", "--exit-code", url, "HEAD"],
                capture_output=True, text=True, timeout=25, cwd=away,
                env=env, stdin=subprocess.DEVNULL, start_new_session=True)
    except FileNotFoundError:
        return [_row("git remote", WARN, "git is not installed",
                     "install git; the git tooling check above says the "
                     "same thing with the remedy")]
    except subprocess.TimeoutExpired:
        return [_row(
            "git remote", WARN,
            "the remote did not answer within 25 s",
            "push from a node with outbound access, or ask the user")]
    except OSError as exc:
        return [_row("git remote", WARN, f"could not be asked: {exc}",
                     "check that git runs, then push from a node with "
                     "outbound access, or ask the user")]

    rows: list[dict] = []
    # A credential helper inherited from somebody else's session can never
    # work and is worth naming even when the remote answers: the socket in
    # it belongs to another uid, so every push through it fails with
    # "Missing or invalid credentials" and an EACCES nobody reads as
    # "wrong owner". Reported, not rewritten -- this is the user's git
    # configuration.
    try:
        helper = subprocess.run(
            ["git", "config", "--get", "credential.helper"],
            capture_output=True, text=True, timeout=10,
            cwd=str(ctx.get("workspace") or ".")).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        helper = ""
    if helper:
        sock = re.search(r"(/run/user/(\d+)/\S*\.sock)", helper)
        if sock and os.getuid() != int(sock.group(2)):
            rows.append(_row(
                "git credential helper", WARN,
                f"points at another account's socket (uid {sock.group(2)}, "
                f"this is {os.getuid()})",
                "git config --unset credential.helper in this checkout, "
                "or let the user push from their own session"))

    if done.returncode == 0:
        return rows + [_row("git remote", PASS, "reachable, and it answers")]
    if done.returncode == 2:
        # --exit-code: the remote answered and has no HEAD -- an empty
        # repository before its first push. Reachable, not a failure.
        return rows + [_row("git remote", PASS,
                            "reachable, and it has no branches yet")]

    text = (done.stderr or "") + "\n" + (done.stdout or "")
    said = " ".join(text.split())[:160]
    found = diagnose(text)
    if found is None:
        return rows + [_row(
            "git remote", WARN, "unreachable: " + said,
            "read the message; a push from here will fail the same way")]
    # git's own words go with the cause: a remedy that says "read the line
    # git printed" is read where that line is not shown (the push gate's
    # refusal, the report). git strips credentials from URLs it prints.
    return rows + [_row("git remote", WARN, f"{found.cause} (git: {said})",
                        found.remedy)]


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


def _isolation_mechanism() -> str:
    """Which mechanism can hold a command on THIS host, or "".

    The product tries three, in this order: bubblewrap, then Landlock on
    a Linux kernel new enough for it, then Seatbelt on macOS. Asking only
    about the first one is how a host that is protected gets told to
    install bubblewrap, and how a host running under Landlock reads as
    unprotected. Order and names follow ``_bash_isolation_argv``, so the
    answer is what will actually run rather than what is installed.
    """
    try:
        from delfin.agent.api_client import (
            _bwrap_functional, _landlock_functional, _seatbelt_functional)
    except Exception:
        return ""
    for name, probe in (("bwrap", _bwrap_functional),
                        ("Landlock", _landlock_functional),
                        ("Seatbelt", _seatbelt_functional)):
        try:
            if probe():
                return name
        except Exception:
            continue
    return ""


#: Pairs that must load in ONE interpreter, and the subsystem each
#: belongs to. The first element is imported first on purpose: the point
#: is the ORDER, because the dynamic loader resolves a shared library once
#: per process and the first resolution wins for everything after it.
#:
#: pymupdf is the agent's PDF reader (office.py) and drags in libmupdf,
#: which needs libstdc++. sqlite3 is pulled by stk via atomlite, and
#: needs a newer C++ ABI through libicu. Where a directory on
#: LD_LIBRARY_PATH ships an older libstdc++ than the interpreter's own
#: libraries need, importing the first makes the second fatal -- and the
#: optional-import fallbacks downstream turn that into a capability that
#: is silently absent rather than an error.
_IMPORT_ORDER_PAIRS: tuple[tuple[str, str, str], ...] = (
    ("pymupdf", "sqlite3", "PDF reading and the structure database"),
)


def _check_loader_path(ctx: dict) -> list[dict]:
    """Can the native stack load in one interpreter, in both orders.

    Asked in a SUBPROCESS, with this process's environment, because the
    loader caches LD_LIBRARY_PATH at process start: by the time any check
    could run in-process the answer is already fixed, and importing the
    pair here would poison the very interpreter doing the asking.

    Reports the directory that won, when it can, rather than only that
    something broke -- a row saying "an import failed" sends the reader
    looking in Python, and the fault is three layers below it.

    Universal: nothing here names ORCA or any site. It imports two
    modules DELFIN itself needs and reports what the loader did. On a host
    whose library order is sound, both orders succeed and the row passes.
    """
    import shutil as _sh
    rows: list[dict] = []
    py = sys.executable or _sh.which("python3") or "python3"
    for first, second, what in _IMPORT_ORDER_PAIRS:
        code = (f"import {first}\n"
                f"import {second}\n"
                "print('ok')")
        try:
            done = subprocess.run([py, "-c", code], capture_output=True,
                                  text=True, timeout=120)
        except (OSError, subprocess.SubprocessError) as exc:
            rows.append(_row(
                "loader path", WARN,
                f"could not ask the interpreter: {exc}",
                f"run: {py} -c 'import {first}, {second}'"))
            continue
        if done.returncode == 0:
            rows.append(_row("loader path", PASS,
                             f"{first} + {second} load together"))
            continue
        err = (done.stderr or "").strip().splitlines()
        detail = err[-1][:300] if err else f"{first} then {second} failed"
        # The loader names the file it took and the version it lacked;
        # pull the PATH out so the fix is a directory and not a guess.
        # Matched as a path rather than by splitting on ":" -- the line
        # opens with "ImportError:", and splitting took that as the
        # library.
        culprit = ""
        for line in err:
            if "not found" not in line:
                continue
            m = re.search(r"(/[^\s:]*\.so(?:\.\d+)*)", line)
            if m:
                culprit = m.group(1)
                break
        fix = (
            "a directory on LD_LIBRARY_PATH provides an older shared "
            "library than this interpreter's own libraries need, and the "
            "loader takes the first match. Put the interpreter's lib "
            "directory first, or drop that entry for Python processes")
        if culprit:
            fix = f"the loader took {culprit} first -- " + fix
        rows.append(_row(
            f"loader path ({what})", WARN, detail, fix))
    return rows


def _check_bash_isolation(ctx: dict) -> list[dict]:
    """Report whether shell commands run in a filesystem namespace.

    Write-target gating refuses paths outside the workspace, but it reads
    the command text — a subprocess started by that command is beyond it.
    Only namespace isolation contains that, and it is opt-in because it
    can disturb cluster workflows. Surfacing the state (and whether the
    machine even supports it) lets the user decide instead of assuming.

    Which mechanism does it is part of the state: a cluster login node
    rarely allows the user namespace bubblewrap needs, and the same
    kernel usually offers Landlock, which the product uses in its place.
    """
    mode = "auto"
    try:
        from delfin.user_settings import load_settings
        mode = str(((load_settings() or {}).get("agent") or {})
                   .get("bash_isolation", "auto") or "auto").strip().lower()
    except Exception:
        pass
    held_by = _isolation_mechanism()

    # For the "off" row only: auto already contains every mode wherever
    # this host can, so proposing bwrap there would propose the state the
    # host is already in.
    fix = ("set agent.bash_isolation = \"bwrap\" to contain shell writes "
           "in every permission mode")
    no_mechanism = ("install bubblewrap, or run on a kernel with Landlock "
                    "(Linux 5.13+); the write-target gate stays active "
                    "either way")
    if mode == "bwrap":
        if held_by:
            return [{"check": "bash isolation", "status": "PASS",
                     "detail": f"{held_by} active for every command"}]
        # Not a silent downgrade: with nothing able to hold the command,
        # the product refuses to run it rather than run it unisolated.
        return [{"check": "bash isolation", "status": "FAIL",
                 "detail": "isolation is switched on and nothing here can "
                           "provide it — shell commands are refused",
                 "fix": "install bubblewrap, or run on a kernel with "
                        "Landlock (Linux 5.13+); or set "
                        "agent.bash_isolation = \"auto\" to fall back to "
                        "the write gate"}]
    if mode == "off":
        return [{"check": "bash isolation", "status": "WARN",
                 "detail": "explicitly off — only the write-target gate "
                           "protects paths outside the workspace",
                 "fix": fix}]
    # "auto", the shipped default. It isolates wherever this host can
    # actually hold a command -- in EVERY permission mode, not only the
    # unattended one.
    #
    # This row said "isolated in bypassPermissions only" and offered the
    # bwrap setting as its fix. That was true of an older resolver: auto
    # used to wall the unattended mode and a locked scope and nothing
    # else. `_bash_isolation_argv` changed -- an approval is given on the
    # command TEXT, and the text is not the act, so an interpreter or a
    # symlink walks past it whether or not somebody is watching -- and
    # this row did not. Measured on a host with bwrap: all four of
    # default, acceptEdits, plan and bypassPermissions come back walled,
    # and none of them can read outside the workspace.
    #
    # Understating protection is not the harmless direction. It tells the
    # user attended sessions are unguarded when they are not, and sends
    # them to change a setting that is already in force.
    if held_by:
        return [{"check": "bash isolation", "status": "PASS",
                 "detail": (f"auto — {held_by} active in every permission "
                            "mode")}]
    # Nothing here can hold a command. Then auto is honest about what is
    # left: the command still runs (refusing every shell command in an
    # attended session would be secure and useless), and what protects
    # paths outside the workspace is the write-target gate reading the
    # command text, plus the socket guard. Saying which is the point of
    # the row.
    return [_row(
        "bash isolation", WARN,
        ("auto — nothing here can isolate a command; the write-target "
         "gate and the socket guard are what hold"),
        no_mechanism,
    )]


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
    ("test runner", "_check_test_runner"),
    ("mcp servers", "_check_mcp"),
    ("scheduler daemon", "_check_scheduler"),
    ("git tooling", "_check_git_tooling"),
    ("loader path", "_check_loader_path"),
    ("git remote", "_check_push"),
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
    # The actionable form is carried through, VALIDATED. This function
    # coerces whatever a check returned, including a monkeypatched one,
    # and what it returns is offered to the user for approval -- so a
    # command is accepted only as a plain one-line string and a setting
    # only as a (dotted key, value) pair. Anything else is dropped and
    # the row keeps its prose.
    command = row.get("command")
    if not isinstance(command, str) or "\n" in command:
        command = ""
    setting = row.get("setting")
    if (isinstance(setting, (list, tuple)) and len(setting) == 2
            and isinstance(setting[0], str) and setting[0].strip()):
        setting = (setting[0], setting[1])
    else:
        setting = None
    return _row(
        str(row.get("check", group)), status,
        str(row.get("detail", "")), str(row.get("fix", "")),
        command=command.strip(), setting=setting,
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


def ready_for_push(workspace: str = ".", *, remote: str | None = None,
                   gh: bool = True) -> list[dict]:
    """Can this checkout push / open a pull request, before either is tried.

    The gate that guards ``git push`` / ``gh pr create`` asks this and
    refuses, with the remedy, when any row is not PASS -- instead of
    letting the command fail and diagnosing afterwards. Five questions,
    five rows, because they fail independently and the remedy differs:
    is git on PATH, is an identity set, does the remote answer, is gh on
    PATH, is gh logged in.

    Every row comes from the two doctor checks ``_check_git_tooling``
    and ``_check_push``, so there is exactly one set of facts about this
    checkout and a second copy cannot drift from it: the gate asks the
    same probes the report shows. Both checks are read-only (--version,
    config --get, git ls-remote, gh auth status) and never push. A test
    shadows ``shutil.which`` and ``subprocess.run`` -- the doctor's own
    no-network trick -- so no live tool and no network is touched.

    Never raises: every probe is caught inside the two checks and
    degrades to a WARN row with a fix. Every row carries prose ``fix``
    only -- installing git or gh is a system-package change and the
    login is the user's (``! gh auth login``), so the ``_row`` contract
    (a check that has an actionable command DECLARES it) and the module
    rule that an agent must not improvise around a system package both
    say: no ``command`` on a readiness row.

    ``remote`` is the remote the push names (a name or a URL; None means
    origin). ``gh=False`` leaves out the two gh rows and their network
    call, for a ``git push`` that does not use gh.
    """
    ctx = {"workspace": workspace, "remote": remote, "skip_gh": not gh}
    return _check_git_tooling(ctx) + _check_push(ctx)


__all__ = ["run_doctor", "format_doctor", "ready_for_push",
           "PASS", "WARN", "FAIL"]
