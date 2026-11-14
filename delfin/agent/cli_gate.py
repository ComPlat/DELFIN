"""DELFIN's file boundary for Claude Code running as the agent backend.

With ``backend="cli"`` the tool loop runs inside the ``claude`` process,
so DELFIN's read and write gates never saw its tools. Under Bypass that
process was started with ``--dangerously-skip-permissions`` and could
read -- and write -- anywhere the account could, while the same session
on any other backend asked before reading outside its folders and
refused writing outside them. Same chip, two different promises.

This module closes that by putting DELFIN's own gate in front of the
CLI's tools, as a PreToolUse hook:

* ``CLIClient`` writes a small state file (workspace, extra and read-only
  roots, the mode, whether anyone can be asked) and hands ``claude`` a
  ``--settings`` hook that runs ``python -m delfin.agent.cli_gate <token>``
  before every tool call.
* The hook rebuilds ``KitToolPermissions`` from that state and asks the
  SAME functions the native tools go through -- ``_check_read_access``,
  ``_gate_write_path``, for a command the native bash gate
  ``_run_permission_gate`` followed by ``_gate_bash_read_paths`` and the
  write-target gate, for WebFetch / WebSearch ``_run_permission_gate``
  (egress host decision, payload check), for an
  ``mcp__`` tool ``_mcp_call_refusal``. Every tool call gets a decision
  (``main``) except the program's own session tools (``_CLI_OWN_TOOLS``);
  one the gate does not know is asked about where a person can be asked
  and refused otherwise. Plan mode is refused here, as natively, instead
  of being left to the program's version.
* A question becomes a request in the approvals directory
  (``file_confirm``). The ``Relay`` in the DELFIN process picks it up and
  puts it to the dashboard panel or the terminal prompt -- the user sees
  the dialog they see on every other backend.
* What a click grants (a directory opened for reads, a refusal) is
  written back to the state file, so the next tool call -- a new hook
  process -- knows it.

A hook that cannot decide denies. "Could not check" is not permission.
"""

from __future__ import annotations

import fcntl
import json
import os
import re
import sys
import threading
import uuid
from collections.abc import Callable
from contextlib import contextmanager
from pathlib import Path
from typing import Any

#: Seconds the hook waits for an answer. A little longer than the dialog
#: itself (300 s), so an expired dialog is seen as expiry, not as a lost
#: answer.
HOOK_WAIT_S = 330.0
#: What ``claude`` is told the hook may take, above HOOK_WAIT_S.
HOOK_TIMEOUT_S = 360

#: Claude Code tools that READ a path, and the argument that names it.
_READ_TOOLS = {
    "Read": "file_path", "NotebookRead": "notebook_path",
    "Grep": "path", "Glob": "path", "LS": "path",
}
#: Claude Code tools that WRITE a path.
_WRITE_TOOLS = {
    "Write": "file_path", "Edit": "file_path", "MultiEdit": "file_path",
    "NotebookEdit": "notebook_path",
}
#: The external CLI program's web tools, as the native tools that decide
#: them (``_run_permission_gate``).
_WEB_TOOLS = {"WebFetch": "web_fetch", "WebSearch": "web_search"}
#: Tools that act on the external CLI program's own session only: its
#: task list, loading a deferred tool's schema, the output or the stop of
#: a background command the gate admitted, a skill's text, a question to
#: the person through the program's own channel, entering and leaving its
#: plan mode. No decision is printed for them, so the program's own mode
#: applies, as before: an "allow" for ExitPlanMode would answer the
#: plan-approval question in the person's place. Names measured on the
#: program's tool list (2.1.322: Agent, Bash, Edit, ListAgents, Read,
#: ReportFindings, ScheduleWakeup, Skill, ToolSearch, Workflow, Write and
#: the deferred CronCreate, CronDelete, CronList, EnterWorktree,
#: ExitWorktree, Monitor, NotebookEdit, SendMessage, TaskCreate, TaskGet,
#: TaskList, TaskStop, TaskUpdate, WebFetch, WebSearch) and the earlier
#: spellings (TodoWrite, BashOutput, KillShell, ...).
_CLI_OWN_TOOLS = frozenset({
    "ToolSearch", "DeferredToolPlaceholder",
    "TodoWrite", "TodoRead", "TaskCreate", "TaskGet", "TaskList",
    "TaskUpdate", "TaskOutput", "TaskStop", "BashOutput", "KillShell",
    "KillBash", "Monitor", "CronList", "ListAgents", "ReportFindings",
    "Skill", "SlashCommand", "AskUserQuestion", "ExitPlanMode",
    "EnterPlanMode",
})
#: Tools that start or steer another agent, or a later turn, inside the
#: external CLI program. Decided like a tool the gate does not know
#: (``_unknown_tool``); named so the refusal says what the tool does.
_AGENT_TOOLS = frozenset({
    "Agent", "Task", "Workflow", "SendMessage", "ScheduleWakeup",
    "CronCreate", "CronDelete",
})


def state_dir() -> Path:
    return Path.home() / ".delfin" / "cli_gate"


def _state_path(token: str) -> Path:
    safe = "".join(c for c in str(token) if c.isalnum() or c in "-_")[:64]
    if not safe:
        raise ValueError("empty gate token")
    return state_dir() / f"{safe}.json"


@contextmanager
def _locked(token: str):
    path = _state_path(token)
    lock = path.with_suffix(".lock")
    fd = os.open(lock, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield path
    finally:
        try:
            fcntl.flock(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)


def _write_json(path: Path, data: dict) -> None:
    tmp = path.with_name(f".{path.name}.{uuid.uuid4().hex[:6]}")
    tmp.write_text(json.dumps(data), encoding="utf-8")
    os.chmod(tmp, 0o600)
    tmp.replace(path)


def new_token() -> str:
    return uuid.uuid4().hex[:16]


def write_state(token: str, *, workspace: str, mode: str, can_ask: bool,
                extra_dirs=(), read_only_dirs=()) -> Path:
    """Create or refresh the state the hook reads. Keeps what earlier
    clicks granted (session read dirs, refusals). In bypass the egress
    proxy of this process is attached (``_attach_proxy``)."""
    d = state_dir()
    d.mkdir(parents=True, exist_ok=True)
    os.chmod(d, 0o700)
    with _locked(token) as path:
        old = _read_state(path)
        state = {
            "workspace": str(workspace),
            "mode": str(mode or "default"),
            "can_ask": bool(can_ask),
            "extra_dirs": [str(p) for p in extra_dirs or ()],
            "read_only_dirs": [str(p) for p in read_only_dirs or ()],
            "session_read_dirs": old.get("session_read_dirs", []),
            "denied_paths": old.get("denied_paths", []),
            "user_named_urls": old.get("user_named_urls", []),
            "web_seen_urls": old.get("web_seen_urls", []),
            "egress_answers": old.get("egress_answers", {}),
        }
        state.update(_attach_proxy(state, token))
        _write_json(path, state)
    return path


#: Permissions the egress proxy keys a CLI session's token on, kept alive
#: here: the proxy forgets a token when its permissions object goes.
_PROXY_PERMS: dict = {}


def _attach_proxy(state: dict, token: str) -> dict:
    """The proxy fields of the state: in bypass, where the hook cages each
    command, the socket and token file of this process's egress proxy
    (the hook's own process ends before the command runs). Empty in any
    other mode, or when the network is not set to go through the proxy;
    the cage then has the network the setting names, or none."""
    if _native_mode(state.get("mode")) != "bypassPermissions":
        return {}
    try:
        from . import egress_proxy as _ep
        if _ep.network_mode() != "proxy":
            return {}
        perms = _PROXY_PERMS.get(token)
        if perms is None:
            perms = _PROXY_PERMS[token] = _permissions(state, token)
        proxy = _ep.get_proxy()
        return {"proxy_socket": proxy.socket_path,
                "proxy_token_file": proxy.token_file_for(perms)}
    except Exception:
        return {}


def _native_mode(mode) -> str:
    """The native name of a CLI permission mode."""
    from . import api_client as ac
    return str(ac._map_kit_permission_mode(str(mode or "default")))


def _read_state(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def hook_settings(token: str) -> dict:
    """The ``hooks`` block for ``claude --settings``."""
    command = f"{sys.executable} -m delfin.agent.cli_gate {token}"
    return {"PreToolUse": [{
        "matcher": "*",
        "hooks": [{"type": "command", "command": command,
                   "timeout": HOOK_TIMEOUT_S}],
    }]}


def _permissions(state: dict, token: str):
    from . import api_client as ac
    from .file_confirm import FileConfirmBroker
    try:
        from . import kit_settings
        persisted = ac._persisted_read_dirs(kit_settings.load())
    except Exception:
        persisted = ()
    session = []
    for d in state.get("session_read_dirs", []):
        try:
            session.append(Path(d))
        except Exception:
            continue
    # No code: this hook runs inside the external agent's process tree and
    # has no human channel; DELFIN's Relay answers it from the person's
    # dialog (approval_code, "What the code does not cover").
    broker = (FileConfirmBroker(session_id=f"cli-gate-{token}",
                                timeout_s=HOOK_WAIT_S, require_code=False)
              if state.get("can_ask") else None)
    perms = ac.KitToolPermissions(
        workspace=Path(state["workspace"]),
        mode=ac._map_kit_permission_mode(state.get("mode") or "default"),
        confirm_callback=broker.callback if broker else None,
        extra_workspace_dirs=tuple(Path(p) for p in state.get("extra_dirs", [])),
        read_only_workspace_dirs=tuple(
            Path(p) for p in state.get("read_only_dirs", [])),
        session_read_dirs=tuple(persisted) + tuple(session),
    )
    try:
        perms.denied_paths.update(state.get("denied_paths", []))
    except Exception:
        pass
    # Where a URL came from, for the egress decision (_web_egress_refusal):
    # what the person typed and what a web result showed, recorded by the
    # DELFIN process (note_user_message, note_tool_result).
    for field in ("user_named_urls", "web_seen_urls"):
        try:
            getattr(perms, field).update(
                str(u) for u in state.get(field, []) if isinstance(u, str))
        except Exception:
            pass
    return perms


#: At most this many URLs of each kind are kept in the state file.
_URLS_KEPT = 2000


def _add_urls(token: str, field: str, urls) -> None:
    """Merge ``urls`` into the state's ``field`` (first seen first, capped
    at ``_URLS_KEPT``). Nothing happens without a state file."""
    urls = [u for u in urls if isinstance(u, str) and u]
    if not urls:
        return
    try:
        with _locked(token) as path:
            state = _read_state(path)
            if not state.get("workspace"):
                return
            have = [u for u in state.get(field, []) if isinstance(u, str)]
            known = set(have)
            for u in urls:
                if len(have) >= _URLS_KEPT:
                    break
                if u not in known:
                    known.add(u)
                    have.append(u)
            state[field] = have
            _write_json(path, state)
    except Exception:
        return


def note_user_message(token: str, content) -> None:
    """The URLs the person typed in ``content`` (a user message the DELFIN
    process sends to the external CLI program), for the hook's egress
    decision: a WebFetch of a URL the person wrote proceeds, as natively
    (``_note_user_named_urls``, under the same rule as ``_grant_push_from``:
    harness text and machine messages name nothing). Never raises."""
    try:
        from . import api_client as ac
        if isinstance(content, ac._HarnessText):
            return
        text = ac._message_text(content)
        if text.lstrip().startswith(ac._MACHINE_MESSAGE_PREFIXES):
            return
        holder = type("_Urls", (), {})()
        holder.user_named_urls = set()
        ac._note_user_named_urls(holder, text)
        _add_urls(token, "user_named_urls", sorted(holder.user_named_urls))
    except Exception:
        return


def note_tool_result(token: str, tool: str, text: str) -> None:
    """The URLs a WebSearch or WebFetch result of the external CLI program
    showed, for the hook's egress decision: a search result can be opened
    (``web_seen_urls``, as ``_note_web_result`` keeps them natively; at
    most ``_WEB_SEEN_PER_PAGE`` per result). Other tools add nothing.
    Never raises."""
    if tool not in _WEB_TOOLS:
        return
    try:
        from . import api_client as ac
        found = []
        for m in ac._URL_IN_CONTENT_RE.finditer(str(text or "")):
            u = m.group(0).rstrip(".,;:)'\">}]")
            if u:
                found.append(u)
            if len(found) >= ac._WEB_SEEN_PER_PAGE:
                break
        _add_urls(token, "web_seen_urls", found)
    except Exception:
        return


def _abs(raw: str, cwd: str, workspace: str) -> Path:
    p = Path(str(raw)).expanduser()
    if not p.is_absolute():
        p = Path(cwd or workspace) / p
    return p


def _unjson(err: str | None) -> str | None:
    """The gates return either prose or a JSON ``{"error": ...}`` (a
    blocking PreToolUse hook as ``{"error": "blocked_by_hook", "reason"}``)."""
    if err is None:
        return None
    try:
        data = json.loads(err)
        if isinstance(data, dict) and data.get("error"):
            if data["error"] == "blocked_by_hook" and data.get("reason"):
                return f"blocked by a PreToolUse hook: {data['reason']}"
            return str(data["error"])
    except Exception:
        pass
    return str(err)


#: The CLI's tool names as the native tools the permission rules name.
_NATIVE_NAMES = {
    "Read": "read_file", "NotebookRead": "notebook_read", "Grep": "grep_file",
    "Glob": "list_files", "LS": "list_files", "Write": "write_file",
    "Edit": "edit_file", "MultiEdit": "multi_edit",
    "NotebookEdit": "notebook_edit", "Bash": "bash", "WebFetch": "web_fetch",
    "WebSearch": "web_search",
}


def _native_call(tool: str, args: dict, cwd: str, workspace: str) -> tuple:
    """(native tool name, arguments) for a CLI tool call; MCP names as
    they are."""
    name = _NATIVE_NAMES.get(tool, tool)
    args = dict(args or {})
    key = _READ_TOOLS.get(tool) or _WRITE_TOOLS.get(tool)
    if key and args.get(key):
        args = {**args, "path": str(_abs(args[key], cwd, workspace))}
    if name == "bash":
        args = {"command": str(args.get("command") or ""), "cwd": cwd}
    return name, args


def decide(token: str, event: dict) -> str | None:
    """Why this tool call is refused, or None when it may run.

    ``event`` is the CLI's PreToolUse input: ``tool_name``,
    ``tool_input``, ``cwd``.
    """
    return decide_call(token, event)[0]


def decide_call(token: str, event: dict) -> tuple[str | None, str | None]:
    """``(refusal, caged command)`` for this tool call.

    The refusal as ``decide``; the caged command is the shell text the CLI
    runs instead of the model's, for a Bash call in bypass (``_bypass_cage``),
    else None.

    By tool: a read and a write as before (plan mode refuses a write); a
    command through ``_decide_bash``; WebFetch / WebSearch through the
    native network gate (``_decide_web``); an ``mcp__`` tool through
    ``_mcp_call_refusal``, which applies the permission rules itself; the
    program's own session tools (``_CLI_OWN_TOOLS``) unchecked; any other
    tool through ``_unknown_tool``.
    """
    from . import api_client as ac
    with _locked(token) as path:
        state = _read_state(path)
    if not state.get("workspace"):
        return ("refused: DELFIN's gate state for this session is missing, "
                "so this call cannot be checked. Restart the session."), None
    perms = _permissions(state, token)
    ex = ac._doc_executor
    tool = str(event.get("tool_name") or "")
    args = event.get("tool_input") or {}
    if not isinstance(args, dict):
        args = {}
    cwd = str(event.get("cwd") or state["workspace"])
    err: str | None = None
    caged: str | None = None
    rule_call = None

    if tool.startswith("mcp__"):
        # The MCP path, as the streaming loop decides an mcp__ call of a
        # native session: managed and user permission rules, the
        # PreToolUse hooks, plan mode, the side-effect question.
        err = _unjson(ac._mcp_call_refusal(tool, dict(args), perms))
        _keep_state(token, perms, ac)
        return err, None

    # The permission rules (permission_rules.py) -- the administrator's
    # managed rules first, then the user's and the project's, one matcher
    # -- on the native name of the call: deny refuses, ask asks through
    # the same relay, and with nobody to ask an ask refuses. An allow rule
    # is not consulted here: what passes is passed by the gates below. An
    # ask rule the person approved answers the gates' own question for
    # this call, as natively (_RULE_CALL "asked", _rule_answers): one
    # question, not two.
    native, native_args = _native_call(tool, args, cwd, state["workspace"])
    hit = ac._permission_rule_hit(perms, native, native_args)
    if hit is not None and hit[0] == "deny":
        return ac._refuse_by_rule(native, hit[1]), None
    # Plan mode, after a deny rule and before an ask rule's question, as in
    # _execute_pinned: nobody is asked about a call plan mode refuses.
    if perms.mode == "plan":
        plan_err = _plan_gate(tool, args)
        if plan_err is not None:
            ac._record_gate_refusal("plan_mode", tool, "")
            return plan_err, None
    if hit is not None:
        refusal = ex._ask_by_rule(perms, native, native_args, hit[1])
        if refusal is not None:
            return refusal, None
        rule_call = ac._RULE_CALL.set({"name": native, "args": native_args,
                                       "asked": True})
    try:
        if tool in _READ_TOOLS:
            raw = args.get(_READ_TOOLS[tool]) or cwd
            err = ex._check_read_access(perms, _abs(raw, cwd,
                                                    state["workspace"]))
        elif tool in _WRITE_TOOLS:
            raw = args.get(_WRITE_TOOLS[tool])
            if raw:
                # The native tool and its own arguments (content,
                # old_string, new_string, edits share their names): the
                # dialog previews the change. With the path alone a new
                # file was "(no changes)".
                err = ex._gate_write_path(
                    str(_abs(raw, cwd, state["workspace"])), perms,
                    _NATIVE_NAMES.get(tool, "write_file"),
                    {**args, "path": str(raw)})
        elif tool == "Bash":
            err, caged = _decide_bash(state, token, perms,
                                      str(args.get("command") or ""), cwd, ex)
        elif tool in _WEB_TOOLS:
            err = _decide_web(state, perms, native, native_args, ex)
        elif tool in _CLI_OWN_TOOLS:
            pass
        else:
            err = _unknown_tool(perms, tool, args, ex)
    finally:
        if rule_call is not None:
            ac._RULE_CALL.reset(rule_call)

    _keep_state(token, perms, ac)
    err = _unjson(err)
    return err, (caged if err is None else None)


def _keep_state(token: str, perms, ac) -> None:
    """Keep what this call granted or refused for the next one: directories
    opened for reads, refused paths, the answers to the unknown-host
    question (``EgressProxy.host_answers``)."""
    with _locked(token) as path:
        fresh = _read_state(path)
        if not fresh.get("workspace"):
            return
        grants = set(fresh.get("session_read_dirs", []))
        grants.update(str(p) for p in perms.session_read_dirs)
        try:
            from . import kit_settings
            saved = {str(p) for p in ac._persisted_read_dirs(
                kit_settings.load())}
        except Exception:
            saved = set()
        fresh["session_read_dirs"] = sorted(grants - saved)
        fresh["denied_paths"] = sorted(
            set(fresh.get("denied_paths", []))
            | set(getattr(perms, "denied_paths", set()) or ()))
        answers = _host_answers(perms)
        if answers:
            merged = dict(fresh.get("egress_answers") or {})
            merged.update(answers)
            fresh["egress_answers"] = merged
        _write_json(path, fresh)


def _live_proxy():
    """This process's egress proxy when it decides hosts itself (not the
    stub that only routes a cage to another process's proxy), else None.
    Never starts one."""
    try:
        from . import egress_proxy as _ep
        proxy = _ep._PROXY
        if isinstance(proxy, _ep.EgressProxy) and not proxy._closed:
            return proxy
    except Exception:
        pass
    return None


def _host_answers(perms) -> dict:
    proxy = _live_proxy()
    if proxy is None:
        return {}
    try:
        return proxy.host_answers(perms)
    except Exception:
        return {}


def _same_dir(a: str, b: str) -> bool:
    try:
        return Path(a).resolve(strict=False) == Path(b).resolve(strict=False)
    except Exception:
        return False


def _decide_bash(state: dict, token: str, perms, cmd: str, cwd: str,
                 ex) -> tuple[str | None, str | None]:
    """``(refusal, caged command)`` for a Bash call, in the native order.

    1. The approval guard (a pending approval is answered by the person).
    2. Plan mode (a read-only form; ``_plan_gate`` refused every other):
       the native bash gate as default mode with nobody to ask decides it
       (deny list, secret paths, tree walks; it passes only on the
       auto-allow list).
       Every other mode: the native bash gate ``_run_permission_gate``
       with this session's permissions -- deny patterns, the git-role push
       gate, the egress scan, the irreversible-needs-a-person rule, the
       auto-allow list and read-only shapes, and the question where a
       person can be asked. The executor's cage admission is not offered
       (no ``caged_admission``): outside bypass the program's shell runs
       outside DELFIN's cage.
    3. The read-path and write-target gates (as ``_execute_bash`` runs
       them after the gate).
    4. Bypass: DELFIN's cage as the command the program runs.
    """
    from . import api_client as ac
    # The same refusal as the native bash gate: a pending approval is
    # answered by the person (approval_guard), in every mode.
    from . import approval_guard as _approval_guard
    _answers = _approval_guard.answer_command(ac._prose_blanked(cmd))
    if _answers:
        ac._record_security_event("approval_self_answer", "Bash",
                                  cmd[:120], blocked=True)
        return _approval_guard.refusal(_answers), None
    _operator = _approval_guard.operator_command(ac._prose_blanked(cmd))
    if _operator:
        ac._record_security_event("operator_action_refused", "Bash",
                                  cmd[:120], blocked=True)
        return _approval_guard.operator_refusal(_operator), None
    # The native gate reads ``cwd`` as the model's own argument: absent
    # when the command runs in the workspace (and then not shown in the
    # question).
    gate_args = {"command": cmd,
                 "cwd": "" if _same_dir(cwd, state["workspace"]) else cwd}
    if perms.mode == "plan":
        # Only a read-only form reaches this (_plan_gate).
        view = _permissions({**state, "mode": "default", "can_ask": False},
                            token)
        err = _unjson(ex._run_permission_gate("bash", gate_args, view))
    else:
        err = _unjson(ex._run_permission_gate("bash", gate_args, perms))
    if err is None:
        err = _unjson(ex._gate_bash_read_paths(cmd, perms, cwd))
    targets: list = []
    if err is None:
        # git apply / patch / checkout -- / restore name their files
        # in the patch or the pathspec, not as a redirection: the same
        # reader the executor's bash gate uses, and the same refusal
        # for a patch that cannot be read from the text.
        patch_paths, patch_err = ex._patch_write_targets(cmd, cwd, perms)
        if patch_err is not None:
            err = _unjson(patch_err)
        else:
            targets = ac._bash_write_targets(cmd) + list(patch_paths or [])
    if err is None:
        for target in targets:
            p = _abs(target, cwd, state["workspace"])
            if ac._is_ephemeral_sink(p, perms.workspace):
                continue
            gate = ex._gate_write_path(str(p), perms, "bash",
                                       {"command": cmd})
            if gate is not None:
                err = (f"blocked: this command would write to "
                       f"'{target}', which is outside what you may "
                       f"modify. {_unjson(gate)}")
                break
    if err is None and perms.mode == "bypassPermissions":
        return _bypass_cage(state, perms, cmd, cwd)
    return err, None


#: Programs a plan-mode command may run: each prints what it reads and
#: changes nothing, with the options ``_plan_segment_read_only`` refuses.
_PLAN_READ_PROGRAMS = frozenset({
    "ls", "pwd", "cat", "head", "tail", "wc", "file", "stat", "which",
    "basename", "dirname", "realpath", "readlink", "echo", "grep", "egrep",
    "fgrep", "rg", "diff", "cmp", "comm", "cut", "tr", "uniq", "nl", "sort",
    "du", "df", "whoami", "uname", "id", "true", "false", "md5sum",
    "sha1sum", "sha256sum", "sha512sum", "cksum", "find", "cd", "git",
})
#: git sub-commands that read the repository.
_PLAN_GIT_READS = frozenset({
    "status", "diff", "log", "show", "blame", "ls-files", "ls-tree",
    "rev-parse", "describe", "shortlog", "cat-file", "rev-list",
    "merge-base", "name-rev", "for-each-ref", "count-objects",
    "check-ignore", "grep", "branch", "remote", "tag", "stash",
})
_PLAN_GIT_BRANCH_FLAGS = frozenset({
    "-a", "--all", "-r", "--remotes", "-v", "-vv", "--verbose",
    "--show-current", "--list", "-l"})
#: Shell text a plan-mode command may not contain: substitution,
#: redirection, sub-shells and groups, escapes, line breaks.
_PLAN_REFUSED_CHARS = frozenset("`$<>(){}\\\n\r")
#: Programs whose options can write or run a program: with an unquoted
#: glob in the command, a planted file name such as ``--output=x`` or
#: ``--pre=./y`` would become one of their options.
_PLAN_GLOB_UNSAFE = frozenset({"find", "sort", "rg", "git", "file", "uniq"})


def _unquoted_glob(cmd: str) -> bool:
    """Whether ``cmd`` has a ``*``, ``?`` or ``[`` outside quotes (one the
    shell expands into file names)."""
    quote = ""
    for c in cmd:
        if quote:
            if c == quote:
                quote = ""
        elif c in "'\"":
            quote = c
        elif c in "*?[":
            return True
    return False


def _plan_segment_read_only(words: list, globbed: bool = False) -> bool:
    """Whether one simple command (its words) only reads. ``globbed``: the
    command has an unquoted glob (``_unquoted_glob``)."""
    if not words:
        return False
    prog, rest = words[0], words[1:]
    if prog not in _PLAN_READ_PROGRAMS:
        return False
    if globbed and prog in _PLAN_GLOB_UNSAFE:
        return False
    if prog == "find":
        return not any(w in ("-exec", "-execdir", "-ok", "-okdir", "-delete",
                             "-fls") or w.startswith("-fprint") for w in rest)
    if prog == "sort":
        # -o FILE / --output writes; --compress-program runs a program.
        return not any(re.match(r"-[A-Za-z]*o", w) or w.startswith("--o")
                       or w.startswith("--com") for w in rest
                       if w.startswith("-"))
    if prog == "uniq":
        return len([w for w in rest if not w.startswith("-")]) <= 1
    if prog == "rg":
        # --pre runs a program on every file; --hostname-bin runs one.
        return not any(w.startswith(("--pre", "--hostname-bin"))
                       for w in rest)
    if prog == "file":
        return not any(w in ("-C", "--compile") for w in rest)
    if prog != "git":
        return True
    i = 0
    while i < len(rest) and rest[i].startswith("-"):
        if rest[i] == "--no-pager":
            i += 1
        elif rest[i] == "-C" and i + 1 < len(rest):
            i += 2
        else:
            return False            # -c, --git-dir, --exec-path, ...
    if i >= len(rest) or rest[i] not in _PLAN_GIT_READS:
        return False
    sub, sargs = rest[i], rest[i + 1:]
    # --output writes a file; --ext-diff runs a program; -O /
    # --open-files-in-pager opens the files in one.
    if any(w.startswith(("--ou", "--ext", "--op")) for w in sargs):
        return False
    if sub == "grep" and any(re.match(r"-[A-Za-z]*O", w) for w in sargs):
        return False
    if sub == "branch":
        return all(w in _PLAN_GIT_BRANCH_FLAGS for w in sargs)
    if sub == "tag":
        return not sargs or sargs[0] in ("-l", "--list")
    if sub == "stash":
        return bool(sargs) and sargs[0] in ("list", "show")
    if sub == "remote":
        from .sandbox import _read_only_form
        return _read_only_form("git", "remote", sargs)
    return True


def _plan_read_only(cmd: str) -> bool:
    """Whether a shell command only reads, so plan mode may run it.

    An allow list, not a filter: every simple command of a pipeline or
    list (``|``, ``&&``, ``||``, ``;``) is one of ``_PLAN_READ_PROGRAMS``,
    spelled by its bare name, with no option that writes or runs a
    program; the text has no substitution, redirection, sub-shell,
    escape, background ``&`` or line break; a program whose options can
    write (``_PLAN_GLOB_UNSAFE``) is not run with an unquoted glob. Anything
    this reader cannot place is not read-only.

    Rejected: the auto-allow list (it admits touch, mkdir, cp, mv,
    git commit, python and pytest, which change files); refusing every
    command as the native plan gate does by name (the external CLI
    program's plan mode runs read-only commands, and its sessions look at
    the repository with git log and grep; the program has no list_files
    equivalent for every such question).
    """
    import shlex
    cmd = str(cmd or "").strip()
    if not cmd or len(cmd) > 4000 or any(c in _PLAN_REFUSED_CHARS for c in cmd):
        return False
    try:
        lexer = shlex.shlex(cmd, posix=True, punctuation_chars=";&|")
        lexer.whitespace_split = True
        tokens = list(lexer)
    except ValueError:
        return False
    globbed = _unquoted_glob(cmd)
    words: list = []
    for tok in tokens + [";"]:
        if tok in ("|", "||", "&&", ";"):
            if not _plan_segment_read_only(words, globbed):
                return False
            words = []
        elif tok and all(c in ";&|" for c in tok):
            return False            # &, |&, ;; and the like
        else:
            words.append(tok)
    return True


def _plan_refusal(tool: str, what: str) -> str:
    """The plan-mode refusal for a CLI tool, worded as the native one."""
    return (f"plan mode (read-only) — '{tool}' {what}, so DELFIN refused it "
            "before it ran. To look at something use Read, or a read-only "
            "command; describe the change in your plan instead. It runs "
            "once the plan is accepted and the mode changed.")


def _what_unknown(tool: str) -> str:
    return ("starts another agent inside the external CLI program"
            if tool in _AGENT_TOOLS else
            "is a tool of the external CLI program DELFIN's gate has no rule "
            "for")


def _plan_gate(tool: str, args: dict) -> str | None:
    """The plan-mode refusal for a CLI tool call, or None to decide it on.

    Plan mode runs what changes nothing, as the native plan gate
    (``_execute_pinned``, ``_PLAN_READONLY_TOOLS``): the read tools, a
    command of a read-only form (``_plan_read_only``) and the program's
    own session tools pass on to their gates; a write tool, any other
    command, WebFetch / WebSearch (a request out of the machine, refused
    in plan mode natively) and any other tool are refused. Before this the
    hook printed an explicit allow for writes and commands in plan mode,
    and only the program's own plan mode refused them (measured, 2.1.322).
    """
    if tool in _READ_TOOLS or tool in _CLI_OWN_TOOLS:
        return None
    if tool in _WRITE_TOOLS:
        return _plan_refusal(tool, "changes files")
    if tool == "Bash":
        if _plan_read_only(str(args.get("command") or "")):
            return None
        return _plan_refusal(
            "Bash", "runs a command that is not one of the read-only forms "
            "plan mode runs (ls, cat, head, grep, rg, find without "
            "-exec/-delete, git status/diff/log/show, ... with no "
            "redirection or substitution)")
    if tool in _WEB_TOOLS:
        return _plan_refusal(tool, "sends a request out of this machine")
    return _plan_refusal(tool, _what_unknown(tool))


def _decide_web(state: dict, perms, native: str, native_args: dict,
                ex) -> str | None:
    """A WebFetch / WebSearch call through the native network gate.

    ``_run_permission_gate(native, args, perms)`` as the native executor
    calls it for web_fetch / web_search: plan mode, the tainted-turn rule,
    the locked-scope question, ``_web_egress_refusal`` (a URL the person
    typed or a web result showed proceeds; else the egress proxy's host
    decision -- a known host passes, an unknown one is asked once per
    session where a person can be asked and refused otherwise; then the
    payload check of the URL or the query). The person's earlier answers
    to the unknown-host question are restored first and kept afterwards
    (``_keep_state``), since every call is a new process here.

    The tainted-turn rule runs but has no input on this backend: the
    state records no untrusted text the program's own tools read
    (``turn_untrusted_sources`` stays empty), so it never refuses here.
    """
    try:
        from . import egress_proxy as _ep
        proxy = _ep.get_proxy()
        if isinstance(proxy, _ep.EgressProxy):
            proxy.restore_host_answers(perms, state.get("egress_answers"))
    except Exception:
        pass
    return _unjson(ex._run_permission_gate(native, dict(native_args), perms))


def _unknown_tool(perms, tool: str, args: dict, ex) -> str | None:
    """A tool the gate has no rule for: an agent tool (``_AGENT_TOOLS``)
    or any name it does not know.

    It may do anything the program lets it do, so it is not let through
    unchecked: refused in bypass (which asks nobody) and where nobody can
    be asked (and in plan mode, by ``_plan_gate``); put to the person
    otherwise, the call shown as it is. Rejected: leaving it to the
    program's own mode, which under bypass runs it without anyone deciding.
    """
    from . import api_client as ac
    what = _what_unknown(tool)
    if perms.mode == "bypassPermissions" or perms.confirm_callback is None:
        why = ("this session runs in bypass, which asks nobody"
               if perms.mode == "bypassPermissions"
               else "nobody can be asked in this session")
        ac._record_gate_refusal("unknown_tool", tool, why)
        return (f"refused: '{tool}' {what}, and {why}. Do the work with the "
                "tools DELFIN decides (Read, Edit, Write, Bash, WebFetch, "
                "WebSearch, the MCP tools), or tell the user what it is for.")
    try:
        shown = json.dumps(args, ensure_ascii=False, default=str)[:600]
    except Exception:
        shown = str(args)[:600]
    preview = f"{tool} {shown}\n(this tool {what})"
    try:
        ok = bool(ac._ask_confirm(perms, tool, dict(args), preview))
    except Exception as exc:
        return f"refused: the question could not be asked ({exc})"
    if ok:
        return None
    if ex._confirm_timed_out(perms):
        return (f"approval for '{tool}' TIMED OUT -- the user is away, this "
                "is NOT a denial. Do not retry it now; say what is waiting "
                "for them.")
    return (f"user denied '{tool}'. Do NOT retry it or reach the same effect "
            "another way." + ex._said_why(perms))


def _bypass_cage(state: dict, perms, cmd: str,
                 cwd: str) -> tuple[str | None, str | None]:
    """``(refusal, caged command)`` for a Bash call in bypass that the
    native gate admitted (``_decide_bash``).

    Bypass asks nothing, so the command gets what a native bypass command
    gets: the whole native gate (deny list, push gate, egress scan, the
    action classifier's second look -- ``_run_permission_gate``, run by
    ``_decide_bash`` in every mode), then DELFIN's cage. The CLI runs its shell outside DELFIN's process, so the
    cage reaches it as text: the command the CLI runs becomes the cage's
    argv, joined for a shell (``updatedInput``). Isolation is forced
    (``mode="bwrap"``: bubblewrap, else Landlock or Seatbelt); with no
    layer the command is refused, as ``_bash_isolation_argv`` refuses a
    native bypass command there. The network goes through the egress proxy
    of the DELFIN process (``_attach_proxy``), the one decision point for
    every route.

    Rejected: refusing the CLI backend bypass outright -- a person who
    picks it would get a session that asks for every command, the opposite
    of what bypass is for, while a cage gives the same promise as on every
    other backend. Rejected: the CLI's own sandbox -- a second policy with
    its own rules and its own network helper (absent on this server), not
    the one the security panel and the audit log describe.
    """
    import shlex
    from . import api_client as ac
    if state.get("proxy_socket") and state.get("proxy_token_file"):
        from . import egress_proxy as _ep
        _ep.use_external(state["proxy_socket"], state["proxy_token_file"])
    argv = ac._bash_isolation_argv(cmd, cwd, perms, mode="bwrap")
    if argv[:2] == ["/bin/sh", "-c"] and argv[3:4] == ["refused"]:
        return str(argv[4]), None
    return None, shlex.join(argv)


def main(argv: list[str]) -> int:
    """Hook entry point: PreToolUse JSON on stdin, a decision on stdout.

    A call DELFIN refuses is denied. Every other call it checked --
    asking the user where that is the rule -- is allowed, so the CLI does
    not ask a second time: under ``-p`` it has nobody to ask and refused
    every write the user had just approved. The CLI's own deny rules (the
    secret files, _cli_secret_deny_rules) still hold over an allow
    (measured with claude 2.1.309). Only the program's own session tools
    (``_CLI_OWN_TOOLS``) get no decision and are left to the CLI.
    """
    token = argv[1] if len(argv) > 1 else ""
    tool = ""
    event: dict = {}
    caged: str | None = None
    try:
        event = json.loads(sys.stdin.read() or "{}")
        tool = str(event.get("tool_name") or "")
        reason, caged = decide_call(token, event)
    except Exception as exc:
        reason, caged = (f"refused: DELFIN's gate could not check this call "
                         f"({type(exc).__name__}: {exc})."), None
    if reason:
        print(json.dumps({"hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": reason,
        }}))
    elif tool not in _CLI_OWN_TOOLS:
        out = {"hookEventName": "PreToolUse",
               "permissionDecision": "allow",
               "permissionDecisionReason": "DELFIN's gate let this call through"}
        if caged is not None:
            # The command runs in DELFIN's cage: the CLI runs this text
            # instead of the model's (measured: the CLI runs a hook's
            # updatedInput with an allow).
            tool_input = event.get("tool_input")
            out["updatedInput"] = {
                **(tool_input if isinstance(tool_input, dict) else {}),
                "command": caged}
            out["permissionDecisionReason"] = ("DELFIN's gate let this call "
                                               "through, in DELFIN's cage")
        print(json.dumps({"hookSpecificOutput": out}))
    return 0


class Relay:
    """Carries the hook's questions to the dialog of the DELFIN process.

    Polls the approvals directory for requests from this session's hook,
    puts each to ``callback`` (the dashboard broker or the terminal one)
    on its own thread, and writes the answer back. An expired dialog is
    left unanswered: the hook then sees expiry, not a refusal.
    """

    def __init__(self, token: str,
                 callback: Callable[[str, dict, str], bool],
                 poll_s: float = 0.3) -> None:
        self.session_id = f"cli-gate-{token}"
        self.callback = callback
        self.poll_s = poll_s
        self._seen: set[str] = set()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True,
                                        name="delfin-cli-gate-relay")

    def start(self) -> Relay:
        self._thread.start()
        return self

    def stop(self) -> None:
        self._stop.set()

    def _run(self) -> None:
        from . import file_confirm
        while not self._stop.is_set():
            try:
                for record in file_confirm.pending():
                    rid = str(record.get("id") or "")
                    if (record.get("session_id") != self.session_id
                            or rid in self._seen):
                        continue
                    self._seen.add(rid)
                    threading.Thread(target=self._ask, args=(record,),
                                     daemon=True).start()
            except Exception:
                pass
            self._stop.wait(self.poll_s)

    def _ask(self, record: dict) -> None:
        from . import file_confirm
        rid = str(record.get("id") or "")
        args = {"path": record.get("path", "")} if record.get("path") else {}
        if record.get("command"):
            # So the dialog can offer to change it, and a grant can match it.
            args["command"] = str(record.get("command"))
        try:
            ok = bool(self.callback(record.get("tool", ""), args,
                                    record.get("preview", "")))
        except Exception:
            ok = False
        owner: Any = getattr(self.callback, "__self__", None)
        if not ok and bool(getattr(owner, "last_timed_out", False)):
            return
        reason = "" if ok else str(getattr(owner, "last_refusal_reason", "")
                                   or "")
        changed = "" if ok else str(getattr(owner, "last_changed_to", "") or "")
        file_confirm._write_answer(rid, ok, by="delfin-cli-gate",
                                   reason=reason, changed_to=changed)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
