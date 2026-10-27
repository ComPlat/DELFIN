"""The boundary a session would run under, read before it starts.

Two read-only reports behind ``delfin-agent capabilities`` and
``delfin-agent gate explain``.

``capabilities(cwd, ...) -> dict``
    The effective permission mode and where it came from; writable,
    read-only and readable roots; the deny and protected path globs and
    the masked home directories (names only, never contents); the shell
    deny patterns and auto-allow list; the isolation the product would
    apply to a command on this host, read off the argv it would build
    (``_bash_isolation_argv``) rather than off what is installed; the
    availability of bubblewrap, Landlock, Seatbelt, the socket guard, the
    process cage and the process guard; who can be asked; the git role
    resolved for the workspace's repository (``git_role``); the session
    tool lists.

``explain(tool, arguments, cwd, ...) -> dict``
    ``verdict`` (``allow`` / ``ask`` / ``refuse``), ``reason_code``,
    ``reason``, the ``questions`` the gate would put, and ``if_approved``
    (``allow`` / ``refuse``: what a yes to every question leads to). It
    calls the executor's own entry point, ``_DocToolExecutor.execute``,
    under ``api_client.gate_dry_run``: every gate runs unchanged, and the
    tool returns after its last gate instead of acting. Questions go to a
    recorder that answers no (a second pass answers yes for
    ``if_approved``). Only ``bash`` and ``write_file`` read the dry-run
    flag, so only they are explained; any other tool is refused with
    ValueError rather than executed.

Nothing here writes into the user's stores. During ``explain`` every
writable sink of ``state_paths.USER_STATE_SINKS`` and
``USER_STATE_RESOLVERS`` points into a temporary directory removed
afterwards (the two settings files stay real, and are only read), and
``security_events.record`` is replaced by a list: no panel entry, no
attention-inbox line, no audit or failure-log record, no approval
request. Session PreToolUse hooks are not run (they are commands); the
number configured is reported instead. Neither function changes a
permission: the policy object is built fresh, used and dropped.

The policy is built by the same functions the client uses
(``api_client._kit_sandbox`` for the kit provider,
``_workspace_sandbox`` for openai and ollama), from the grants
``launch_guard.inspect_launch_dir`` accepts, with the mode
``launch_guard.resolve_posture`` gives. Other providers do not run
DELFIN's executor gate and are refused with ValueError.

Rejected: re-implementing the gate as a pure predicate. The executor's
decision is spread over the entry point, the policy gate and the tool's
own path gates; a second copy is a second answer, and two answers drift.
"""

from __future__ import annotations

import contextlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

ASKERS = ("terminal", "approvals", "nobody")
MODES = ("plan", "default", "acceptEdits", "bypassPermissions")
GATED_PROVIDERS = ("kit", "openai", "ollama")

#: Containment announcements made once per process. They say what holds
#: a command, not why one was allowed or refused, so they never become a
#: reason code (and a once-per-process event cannot make two otherwise
#: identical decisions look different).
#: ``network`` / ``network_missing`` are api_client._announce_child_network:
#: what holds the network of the commands, recorded by the first command
#: of a process.
_ANNOUNCEMENTS = frozenset({"isolation", "isolation_missing",
                            "handover_root", "network", "network_missing"})

#: Events that record the answer to a question rather than its cause.
_ANSWERS = frozenset({"denied_by_user", "approval_timeout"})

#: The two settings files: read while the gate decides, never written.
_READ_ONLY_SINKS = frozenset({
    ("delfin.agent.hooks_editor", "_USER_SETTINGS"),
    ("delfin.agent.kit_settings", "USER_SETTINGS_PATH"),
})

#: Refusal text -> reason code, for refusals that record no security
#: event. First match wins; matched case-insensitively.
_TEXT_CODES: tuple[tuple[str, str], ...] = (
    ("managed policy", "managed_policy"),
    ("plan mode", "plan_mode"),
    ("was excluded from this session", "session_tool_list"),
    ("is not available to the", "role"),
    ("blocked_by_hook", "hook"),
    ("without a prior read_file", "read_first"),
    ("modified since last read_file", "read_first"),
    ("requires permissions to be configured", "no_gate"),
    ("no approval dialog", "nobody_to_ask"),
    ("no one to ask", "nobody_to_ask"),
    ("not on the auto-allow", "not_auto_allowed"),
    ("escapes workspace sandbox", "outside_ws"),
    ("outside the allowed workspace roots", "outside_ws"),
    ("read-only", "read_only_write"),
    ("blocked:", "self_preservation"),
    ("is required", "malformed"),
)


# ---------------------------------------------------------------------------
# The policy a session would get
# ---------------------------------------------------------------------------

def _settings_provider() -> str:
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent", {}) or {}
        return str(agent.get("provider", "") or "")
    except Exception:
        return ""


def default_asker(supervised: bool = False) -> str:
    """Who a session started from here could ask, by the rule ``chat``
    applies: a terminal on stdin, else the approvals directory with
    ``--supervised``, else nobody."""
    try:
        if sys.stdin is not None and sys.stdin.isatty():
            return "terminal"
    except Exception:
        pass
    return "approvals" if supervised else "nobody"


def resolve_mode(cwd: Path, permission_mode: str = "") -> tuple[str, str]:
    """(mode, where it came from). An explicit mode is taken as given --
    including bypassPermissions, which an interactive ``chat`` reaches
    only with ``--unattended`` beside it; without one, the persisted
    default as ``chat`` would read it, else plan."""
    from . import launch_guard
    from .cli import _persisted_default_mode
    persisted, source = _persisted_default_mode(Path(cwd))
    return launch_guard.resolve_posture(
        flag_mode=permission_mode or "", persisted_mode=persisted,
        unattended_opt_in=True,
        settings_path=source or "~/.delfin/settings.json")


def session_permissions(
    cwd: Any,
    *,
    permission_mode: str = "",
    add_dirs: Iterable[str] = (),
    read_dirs: Iterable[str] = (),
    allowed_tools: Iterable[str] = (),
    disallowed_tools: Iterable[str] = (),
    asker: str = "nobody",
    provider: str = "",
    confirm_callback: Optional[Callable[[str, dict, str], bool]] = None,
):
    """``(KitToolPermissions, info)`` for a session started in ``cwd``.

    ``confirm_callback`` is bound only when ``asker`` is not ``nobody``.
    ``info`` carries the launch report, the mode and its source, and the
    provider. Raises ValueError for a provider without the executor gate
    or an unknown asker; propagates the constructor's refusal of a home
    directory or a system root as a workspace.
    """
    from . import api_client as ac
    from . import launch_guard
    if asker not in ASKERS:
        raise ValueError(f"asker must be one of {ASKERS}, got {asker!r}")
    provider = (provider or _settings_provider() or "kit").strip().lower()
    if provider not in GATED_PROVIDERS:
        raise ValueError(
            f"provider {provider!r} does not run DELFIN's executor gate "
            f"(only {', '.join(GATED_PROVIDERS)} do), so there is nothing "
            "here to explain")
    workspace = Path(cwd).expanduser().resolve()
    report = launch_guard.inspect_launch_dir(
        workspace, add_dirs=tuple(add_dirs or ()),
        read_dirs=tuple(read_dirs or ()))
    mode, mode_source = resolve_mode(workspace, permission_mode)
    callback = confirm_callback if asker != "nobody" else None
    granted = [str(p) for p in report.granted_dirs]
    readable = [str(p) for p in report.read_dirs]
    if provider == "kit":
        perms = ac._kit_sandbox(str(workspace), mode, callback, granted,
                                readable, None)
    else:
        perms = ac._workspace_sandbox(str(workspace), mode, callback,
                                      readable, None)
    perms.session_allowed_tools = ac._normalise_tool_names(
        list(allowed_tools or ()))
    perms.session_denied_tools = ac._normalise_tool_names(
        list(disallowed_tools or ()))
    info = {"launch": report, "mode": perms.mode, "mode_source": mode_source,
            "provider": provider, "asker": asker}
    return perms, info


def assume_read(perms, tool: str, arguments: dict) -> None:
    """Mark an existing write_file target as read this session, so the
    verdict is about the permission and not about the read-before-write
    precondition a session meets by reading first."""
    if tool != "write_file":
        return
    raw = str(arguments.get("path") or arguments.get("file_path") or "")
    if not raw:
        return
    p = Path(raw).expanduser()
    if not p.is_absolute():
        p = Path(perms.workspace) / p
    try:
        p = p.resolve()
        if p.is_file():
            perms.read_tracker[str(p)] = p.stat().st_mtime
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Classification: one function for the explained and the executed result
# ---------------------------------------------------------------------------

def _load(result: str) -> Any:
    try:
        return json.loads(result)
    except (TypeError, ValueError):
        return None


def _refusal_argv(argv: list) -> bool:
    return argv[:2] == ["/bin/sh", "-c"] and argv[3:4] == ["refused"]


def classify(tool: str, result: str, questions: list, events: list,
             executed: Optional[bool] = None) -> dict:
    """The verdict a tool result stands for.

    Input: the result string of ``execute``, the questions put, the
    security events recorded (``{"kind", "blocked"}``). Output:
    ``verdict``, ``reason_code``, ``reason``, ``questions``, ``reached``
    (the tool got past its last gate). Used for the dry run and for a
    real run alike, so the two are compared on one reading.
    """
    from .api_client import GATE_DRY_RUN_KEY
    obj = _load(result)
    text = str(result or "")
    reached = False
    reason = ""
    isolation_refused = False
    if isinstance(obj, dict) and obj.get(GATE_DRY_RUN_KEY):
        reached = True
        if _refusal_argv(list(obj.get("argv_head") or [])):
            isolation_refused = True
    elif tool == "bash" and isinstance(obj, dict) and "exit_code" in obj:
        reached = True
        err = str(obj.get("stderr") or "")
        if (obj.get("exit_code") == 126 and err.startswith("refused: ")
                and "The command was not run." in err):
            isolation_refused = True
    elif tool == "write_file" and not (isinstance(obj, dict)
                                       and obj.get("error")):
        reached = not text.lstrip().startswith('{"error"')
    if isinstance(obj, dict) and obj.get("error"):
        reason = str(obj.get("error"))
        if obj.get("reason"):
            reason += ": " + str(obj.get("reason"))
    if isolation_refused:
        reached = False
        reason = reason or ("no isolation on this host can hold a command "
                            "for this session, so it is not run")
    kinds = [e for e in events if e.get("kind") not in _ANNOUNCEMENTS]
    if questions:
        # Why it asks, not what the recorder answered.
        why = [e for e in kinds if e.get("kind") not in _ANSWERS]
        verdict = "ask"
        code = why[0]["kind"] if why else "confirm_required"
    elif reached:
        verdict = "allow"
        code = kinds[0]["kind"] if kinds else "allowed"
    else:
        verdict = "refuse"
        blocked = [e for e in kinds if e.get("blocked")]
        if isolation_refused:
            code = "isolation_unavailable"
        elif blocked:
            code = blocked[0]["kind"]
        else:
            low = reason.lower() or text.lower()
            code = next((c for needle, c in _TEXT_CODES if needle in low),
                        "refused")
    return {"verdict": verdict, "reason_code": code,
            "reason": reason[:2000], "questions": list(questions),
            "reached": reached}


# ---------------------------------------------------------------------------
# explain
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def _events_into(events: list):
    """``security_events.record`` appends ``{"kind", "blocked"}`` to
    ``events`` instead of reaching the panel and the attention inbox; and
    building an isolation argv names the network mode without starting
    the egress proxy (a socket and a thread) or asking SLURM for its
    controllers. Neither changes what is allowed: the network flags are
    not read by any gate."""
    from . import api_client as ac
    from . import security_events
    real_record = security_events.record
    real_net = ac._sandbox_network_args

    def _capture(kind, tool, detail, *, blocked=True):
        events.append({"kind": str(kind), "blocked": bool(blocked)})

    def _net_without_proxy():
        try:
            from . import egress_proxy
            return ["--net", egress_proxy.network_mode()], ""
        except Exception:
            return ["--net", "none"], ""

    # The network notice is spent once per process. The dry run must leave
    # it for the session that really runs a command here.
    announced = ac._CHILD_NETWORK_ANNOUNCED
    security_events.record = _capture
    ac._sandbox_network_args = _net_without_proxy
    try:
        yield
    finally:
        security_events.record = real_record
        ac._sandbox_network_args = real_net
        ac._CHILD_NETWORK_ANNOUNCED = announced


@contextlib.contextmanager
def _nothing_recorded(events: list):
    """Every writable user-state sink in a temporary directory, security
    events into ``events``, the dry-run flag on. Removed on exit."""
    from . import api_client as ac
    from .state_paths import RedirectedUserState
    scratch = tempfile.mkdtemp(prefix="delfin-gate-explain-")
    try:
        with _events_into(events):
            with RedirectedUserState(scratch, keep=_READ_ONLY_SINKS):
                with ac.gate_dry_run():
                    yield
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _pass(tool: str, arguments: dict, session_kw: dict,
          answer: bool) -> dict:
    from .api_client import _DocToolExecutor
    questions: list = []

    def _recorder(name, args, preview):
        questions.append(str(preview))
        return answer

    events: list = []
    with _nothing_recorded(events):
        perms, info = session_permissions(confirm_callback=_recorder,
                                          **session_kw)
        perms.skip_hook_discovery = True
        assume_read(perms, tool, arguments)
        result = _DocToolExecutor().execute(tool, dict(arguments), perms)
    out = classify(tool, result, questions, events)
    out["_info"] = info
    return out


def _hook_count(workspace: Path) -> int:
    try:
        from . import hooks
        return len(hooks.load_hooks(workspace).for_event("PreToolUse"))
    except Exception:
        return 0


def explain(tool: str, arguments: dict, *, cwd: Any = ".",
            permission_mode: str = "", add_dirs: Iterable[str] = (),
            read_dirs: Iterable[str] = (), allowed_tools: Iterable[str] = (),
            disallowed_tools: Iterable[str] = (), asker: str = "terminal",
            provider: str = "") -> dict:
    """What the gate would do with one call. See the module docstring."""
    from .api_client import GATE_DRY_RUN_TOOLS
    if tool not in GATE_DRY_RUN_TOOLS:
        raise ValueError(
            f"only {', '.join(sorted(GATE_DRY_RUN_TOOLS))} can be explained "
            f"without being run; {tool!r} cannot")
    session_kw = dict(
        cwd=cwd, permission_mode=permission_mode,
        add_dirs=tuple(add_dirs or ()), read_dirs=tuple(read_dirs or ()),
        allowed_tools=tuple(allowed_tools or ()),
        disallowed_tools=tuple(disallowed_tools or ()),
        asker=asker, provider=provider)
    first = _pass(tool, arguments, session_kw, answer=False)
    info = first.pop("_info")
    if first["verdict"] == "ask":
        second = _pass(tool, arguments, session_kw, answer=True)
        if_approved = "allow" if second["reached"] else "refuse"
    else:
        if_approved = first["verdict"]
    workspace = Path(cwd).expanduser().resolve()
    first.pop("reached", None)
    return {
        "tool": tool,
        "arguments": {k: (v if k != "content" else f"<{len(str(v))} chars>")
                      for k, v in arguments.items()},
        "workspace": str(workspace),
        "mode": info["mode"],
        "mode_source": info["mode_source"],
        "provider": info["provider"],
        "asker": asker,
        **first,
        "if_approved": if_approved,
        "hooks_not_run": _hook_count(workspace),
    }


# ---------------------------------------------------------------------------
# capabilities
# ---------------------------------------------------------------------------

def _describe_argv(argv: list) -> dict:
    """What holds a command started with ``argv``."""
    from . import landlock_exec
    helper = str(Path(landlock_exec.__file__).resolve())
    argv = [str(a) for a in argv]
    if _refusal_argv(argv):
        return {"command_isolation": "refused",
                "detail": argv[-1] if argv else ""}
    in_cage = bool(argv) and argv[0] == "bwrap" and "--unshare-pid" in argv
    from .container_exec import RUNTIMES
    runtime = Path(argv[0]).name if argv else ""
    if runtime in RUNTIMES and argv[1:2] in (["run"], ["exec"]):
        fs = f"container ({runtime})"
    elif argv and argv[0] == "bwrap" and argv[1:4] == ["--ro-bind", "/", "/"]:
        fs = "bwrap"
    elif helper in argv and "--write" in argv:
        fs = "Landlock"
    elif argv and "sandbox-exec" in argv[0]:
        fs = "Seatbelt"
    else:
        fs = "none"
    return {"command_isolation": fs, "process_cage_applied": in_cage,
            "socket_guard_applied": helper in argv}


def _sandbox_layers(perms) -> dict:
    from . import api_client as ac
    out: dict = {}

    def _probe(name, fn):
        try:
            out[name] = bool(fn())
        except Exception:
            out[name] = False

    _probe("bwrap", ac._bwrap_functional)
    _probe("landlock", ac._landlock_functional)
    _probe("seatbelt", ac._seatbelt_functional)
    try:
        from . import landlock_exec
        out["landlock_abi"] = int(landlock_exec.abi_version())
    except Exception:
        out["landlock_abi"] = 0
    try:
        from . import socket_guard
        out["socket_guard"] = bool(socket_guard.available())
    except Exception:
        out["socket_guard"] = False
    cage_off = os.environ.get(ac._PROCESS_CAGE_ENV, "").strip().lower() == "off"
    _probe("process_cage", lambda: (not cage_off)
           and ac._process_cage_functional())
    guard_off = os.environ.get("DELFIN_PROCESS_GUARD", "").strip().lower() == "off"
    supported = sys.platform.startswith("linux") or sys.platform == "darwin"
    out["process_guard"] = {
        "enabled": bool(supported and not guard_off),
        "applied_to": "chat and run processes, at start",
        "off_by": "DELFIN_PROCESS_GUARD=off" if guard_off else "",
    }
    setting = "auto"
    try:
        from delfin.user_settings import load_settings
        setting = str(((load_settings() or {}).get("agent") or {})
                      .get("bash_isolation", "auto") or "auto")
    except Exception:
        pass
    out["isolation_setting"] = ac._BASH_ISOLATION_OVERRIDE or setting
    if out["isolation_setting"].strip().lower() == "container":
        try:
            from . import container_exec
            c = container_exec.choose()
            out["container"] = {"runtime": c.runtime, "image": c.image,
                                "unavailable_because": c.reason}
        except Exception as exc:
            out["container"] = {"runtime": "", "image": "",
                                "unavailable_because": str(exc)}
    try:
        from . import egress_proxy
        out["network"] = egress_proxy.network_mode()
    except Exception:
        out["network"] = ""
    try:
        argv = ac._bash_isolation_argv("true", perms.workspace, perms)
    except Exception:
        argv = []
    out.update(_describe_argv(argv))
    return out


def _git(workspace: Path) -> dict:
    from . import git_role as GR
    key = GR.repo_key(workspace)
    role, source = GR.role_for_keys([key] if key else [],
                                    fallback=GR.global_role)
    return {"repository_key": key or "", "role": role, "source": source}


def _managed_summary() -> dict:
    """managed_policy.summary(), or an inactive one when it cannot load."""
    try:
        from . import managed_policy
        return managed_policy.summary()
    except Exception as exc:
        return {"active": False, "sources": [], "errors": [str(exc)]}


def _rules_summary(workspace: Path) -> dict:
    """The permission rules in force, managed and the user's and the
    project's together, each with its origin and source file
    (permission_rules.rule_records), in evaluation order."""
    try:
        from . import permission_rules as PR
        rs = PR.load(workspace)
        notice = PR.scope_notice(rs)
        return {"rules": PR.rule_records(rs),
                "repository": rs.repo_key or "",
                "allow_off": [r.text for r in rs.allow_off],
                "notice": notice or "",
                "errors": PR.error_lines(rs)}
    except Exception as exc:
        return {"rules": [], "errors": [f"cannot be read: {exc}"]}


#: Home-relative paths the bwrap sandbox hides (``masked_home_dirs`` in
#: the report). Path names, not secrets. A separate literal on purpose:
#: code scanning classifies the sandbox's list by its variable name as a
#: secret and follows every alias or import of that object, so reading it
#: here made ``capabilities --json`` a clear-text secret leak
#: (py/clear-text-logging-sensitive-data); an alias was still followed.
#: tests/test_gate_patterns_finish_in_linear_time.py pins this tuple
#: equal, in order, to the sandbox's list.
MASKED_HOME_ENTRIES: tuple[str, ...] = (
    ".ssh", ".aws", ".gnupg", ".kube", ".docker", ".netrc", ".pgpass",
    ".azure", ".gcp", ".config/gcloud", ".config/gh",
    ".git-credentials", ".npmrc", ".pypirc",
    ".anthropic", ".openai", ".claude",
    ".delfin/credentials.json",
    ".delfin/approval_codes",
    ".delfin/trust_store",
    # Which sessions the phone may reach (delfin.messenger.activation).
    ".delfin/messenger_sessions",
)


def capabilities(*, cwd: Any = ".", permission_mode: str = "",
                 add_dirs: Iterable[str] = (), read_dirs: Iterable[str] = (),
                 allowed_tools: Iterable[str] = (),
                 disallowed_tools: Iterable[str] = (), asker: str = "nobody",
                 provider: str = "") -> dict:
    """The boundary a session started with these options would have."""
    with _events_into([]):
        perms, info = session_permissions(
            cwd, permission_mode=permission_mode, add_dirs=add_dirs,
            read_dirs=read_dirs, allowed_tools=allowed_tools,
            disallowed_tools=disallowed_tools, asker=asker,
            provider=provider, confirm_callback=lambda *a: False)
        sandbox = _sandbox_layers(perms)
    report = info["launch"]
    writable = [str(p) for p in perms.all_workspace_roots()]
    readable = [str(p) for p in perms.all_readable_roots()]
    workspace = Path(perms.workspace)
    return {
        "workspace": str(workspace),
        "provider": info["provider"],
        "mode": perms.mode,
        "mode_source": info["mode_source"],
        "asker": asker,
        "scope_locked": bool(perms.scope_locked),
        "writable_roots": writable,
        "read_only_roots": [str(p) for p in perms.read_only_workspace_dirs],
        "confirm_write_roots": [str(p) for p in perms.confirm_write_dirs],
        "readable_roots": readable,
        "deny_paths": list(perms.path_deny_globs),
        "protected_paths": list(perms.path_protected_globs),
        "masked_home_dirs": list(MASKED_HOME_ENTRIES),
        "write_scope_globs": list(perms.write_allow_globs),
        "bash": {"deny_patterns": list(perms.bash_deny_patterns),
                 "auto_allow_patterns": list(perms.bash_auto_allow_patterns)},
        "tools": {"allowed": sorted(perms.session_allowed_tools),
                  "denied": sorted(perms.session_denied_tools)},
        "sandbox": sandbox,
        "git": _git(workspace),
        "pre_tool_hooks": _hook_count(workspace),
        "managed_policy": _managed_summary(),
        "permission_rules": _rules_summary(workspace),
        "launch": {"refused": bool(report.refused),
                   "findings": [{"level": f.level, "code": f.code,
                                 "message": f.message}
                                for f in report.findings]},
    }


# ---------------------------------------------------------------------------
# Text renderings
# ---------------------------------------------------------------------------

_ASKER_WORDS = {
    "terminal": "the person at this terminal",
    "approvals": "whoever answers the approvals directory",
    "nobody": "nobody (what would be asked is refused or decided by mode)",
}


def format_capabilities(caps: dict) -> str:
    sb = caps["sandbox"]
    yes = lambda b: "yes" if b else "no"  # noqa: E731
    lines = [
        f"workspace   {caps['workspace']}",
        f"provider    {caps['provider']}",
        f"mode        {caps['mode']}  [{caps['mode_source']}]",
        f"asks        {_ASKER_WORDS.get(caps['asker'], caps['asker'])}",
    ]
    if caps["scope_locked"]:
        lines.append("scope       locked to the workspace")
    for label, key in (("writable", "writable_roots"),
                       ("read-only", "read_only_roots"),
                       ("readable", "readable_roots")):
        for i, root in enumerate(caps[key]):
            lines.append(f"{label if i == 0 else '':<11} {root}")
    lines.append(f"deny paths  {', '.join(caps['deny_paths'])}")
    lines.append(f"masked      ~/{', ~/'.join(caps['masked_home_dirs'])}")
    lines.append(f"shell       {len(caps['bash']['deny_patterns'])} deny "
                 f"patterns, {len(caps['bash']['auto_allow_patterns'])} "
                 "auto-allowed")
    tools = caps["tools"]
    lines.append("tools       "
                 + (f"only {', '.join(tools['allowed'])}" if tools["allowed"]
                    else "no allow list")
                 + (f"; refused {', '.join(tools['denied'])}"
                    if tools["denied"] else ""))
    held = [w for w, on in (("socket guard", sb.get("socket_guard_applied")),
                            ("process cage", sb.get("process_cage_applied")))
            if on]
    lines.append(f"isolation   {sb['command_isolation']} for a command"
                 + (f", with {' and '.join(held)}" if held else "")
                 + f"  (setting: {sb['isolation_setting']}, network: "
                 f"{sb.get('network') or '?'})")
    box = sb.get("container")
    if box:
        lines.append("container   " + (
            f"{box['runtime']}, image {box['image']}" if box["runtime"]
            else f"not used: {box['unavailable_because']}"))
    lines.append(f"layers      bwrap {yes(sb['bwrap'])} · Landlock "
                 f"{yes(sb['landlock'])}"
                 + (f" (ABI {sb['landlock_abi']})" if sb["landlock_abi"]
                    else "")
                 + f" · Seatbelt {yes(sb['seatbelt'])} · socket guard "
                 f"{yes(sb['socket_guard'])} · process cage "
                 f"{yes(sb['process_cage'])} · process guard "
                 f"{yes(sb['process_guard']['enabled'])}")
    git = caps["git"]
    lines.append(f"git role    {git['role']}  [{git['source']}"
                 + (f": {git['repository_key']}" if git["repository_key"]
                    else "") + "]")
    if caps["pre_tool_hooks"]:
        lines.append(f"hooks       {caps['pre_tool_hooks']} PreToolUse")
    managed = caps.get("managed_policy") or {}
    if managed.get("active") or managed.get("errors"):
        from . import managed_policy
        # Its deny/ask/only_allow are listed with the user's rules below.
        lines.extend(managed_policy.summary_lines(include_rules=False))
    rules = caps.get("permission_rules") or {}
    if rules.get("repository"):
        lines.append(f"{'repository':<11} {rules['repository']}")
    for i, r in enumerate(rules.get("rules") or ()):
        scope = r.get("scope") or ""
        where = (f", {scope}" if r["list"] == "allow" and scope else "")
        lines.append(f"{'rules' if i == 0 else '':<11} {r['list']:<10} "
                     f"{r['rule']}  [{r['origin']}: {r['source']}{where}]")
    for err in rules.get("errors") or ():
        lines.append(f"rules       {err}")
    if rules.get("notice"):
        lines.append(f"rules       note: {rules['notice']}")
    if caps["launch"]["refused"]:
        lines.append("launch      REFUSED: " + "; ".join(
            f["message"] for f in caps["launch"]["findings"]
            if f["level"] == "refuse"))
    return "\n".join(lines)


def format_explain(said: dict) -> str:
    target = (said["arguments"].get("command")
              or said["arguments"].get("path") or "")
    lines = [f"{said['verdict'].upper()}  ({said['reason_code']})  "
             f"{said['tool']}: {target}",
             f"mode {said['mode']} · asks "
             f"{_ASKER_WORDS.get(said['asker'], said['asker'])}"]
    for q in said["questions"]:
        lines.append(f"question: {q}")
    if said["verdict"] == "ask":
        lines.append(f"on yes: {said['if_approved']}")
    if said["reason"]:
        lines.append(f"reason: {said['reason']}")
    if said["hooks_not_run"]:
        lines.append(f"note: {said['hooks_not_run']} PreToolUse hook(s) "
                     "configured; explain does not run them")
    return "\n".join(lines)


__all__ = ["ASKERS", "MODES", "GATED_PROVIDERS", "default_asker",
           "resolve_mode", "session_permissions", "assume_read", "classify",
           "explain", "capabilities", "format_capabilities",
           "format_explain"]
