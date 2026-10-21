"""Who may answer a pending approval: the person, never an agent session.

An approval waits as a file. ``file_confirm`` puts a headless session's
question in the approvals directory, ``terminal_confirm`` publishes a
terminal or dashboard question in its own directory, and either is
answered by an ``<id>.answer.json`` beside it -- written by
``delfin-agent approvals approve|deny|answer``. The CLI backend's hook
(``cli_gate``) keeps the session's roots and mode in a state file. The
answer checks (owner, mode, age, id) establish that the answer came from
this account, which is also the account the agent's shell runs as. So
none of them tells the person from the session that is waiting, or from
another session of the same account.

This module is the one place that names those locations and those
commands. Its callers:

* the permission gate (``api_client._run_permission_gate``) refuses a
  shell command that answers an approval or names an answer location, in
  every mode, before anything is asked;
* the write gate (``api_client._DocToolExecutor._gate_write_path``)
  refuses a write into an answer location from every write tool and
  every shell write target it reads;
* the bubblewrap cage binds an answer location read-only when it lies
  under a writable root;
* the CLI backend hook refuses the same commands for its Bash tool;
* ``cmd_approvals`` and the answer functions refuse when the process
  carries :data:`ENV_MARKER`, which every agent child process gets from
  ``api_client._scrubbed_bash_env``.

Reading the queue (``approvals ls``, ``show``, ``watch``) stays allowed:
showing the request to the person is what the session should do.

Rejected: a confirm dialog for these routes (the person would be asked to
approve the agent approving itself, and an unattended session has nobody
to ask); refusing ``~/.delfin`` as a workspace root (it moves every
answer-room path test-wide, since the suite redirects the rooms into
``tmp_path``, and does not cover a room reached through the shell);
detecting an outside agent by its own environment variables (a person's
terminal can inherit them from the editor that hosts that agent).

Limits: the shell route is read from the command text. With no
filesystem isolation, a command that builds both the command words and
the answer path at run time (``a=appro; ...``), and also strips the
marker from its environment, is not recognised. Where bubblewrap,
Landlock or Seatbelt holds the command, the rooms are outside every
writable root (or bound read-only) and such a write fails.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Optional

#: Set in every agent child process (shell tool, background job, test
#: run, external agent CLI). Its presence means "an agent session started
#: this process"; the value is the DELFIN process id, for the record only.
ENV_MARKER = "DELFIN_AGENT_SESSION"

#: The answer locations under ``~/.delfin``, by name. The module
#: attributes below may be redirected (tests, a benchmark attempt); both
#: the redirected and the canonical location are protected.
_ROOM_NAMES = ("approvals", "terminal_confirmations", "cli_gate",
               "approval_codes")

#: The subcommands that answer. ``ls``, ``show`` and ``watch`` read.
_ANSWER_VERBS = ("approve", "deny", "answer")

# After quotes and backslashes are removed: ``approvals`` followed by an
# answering verb, flags allowed between. ``del'fin'-agent "approvals"
# appro\ve`` reads as the shell will run it.
_ANSWER_COMMAND_RE = re.compile(
    r"(?<![\w.-])approvals\s+(?:-{1,2}[\w-]+(?:=\S+)?\s+)*"
    r"(" + "|".join(_ANSWER_VERBS) + r")(?![\w-])", re.IGNORECASE)

# ``service stop|restart|end`` under any program name, after the same
# unquoting: stop and restart end the agent service and every session in
# it, and when the socket does not answer they do so by signal, past the
# service's watch-only rule for an agent's processes; end ends sessions
# of the user's. ``status`` and ``logs`` read.
_SERVICE_STOP_RE = re.compile(
    r"(?<![\w.-])service\s+(?:-{1,2}[\w-]+(?:=\S+)?\s+)*"
    r"(stop|restart|end)(?![\w-])", re.IGNORECASE)

# The Python functions that write an answer. Names only these functions
# carry match anywhere; the generic ``answer`` matches only together with
# an import of one of the answering modules, so ``grep -n "answer("
# delfin/agent/file_confirm.py`` (reading the code) is not refused.
_ANSWER_ONLY_NAMES_RE = re.compile(
    r"\b(?:answer_waiting|write_choice_answer|cmd_approvals|_write_answer)\b")
_ANSWER_MODULE_RE = re.compile(
    r"\b(?:file_confirm|approval_answers|terminal_confirm|cli_gate)\b")
_ANSWER_CALL_RE = re.compile(r"\banswer\s*\(|\bimport\b[^\n;]*\banswer\b")
_IMPORT_RE = re.compile(r"\bimport\b|__import__|importlib")

# The agent service's socket and run files, where a window answers a
# session's approvals: ``/…/delfin-<uid>/agent.sock`` and
# ``.delfin/agent_service/run``, joined the same three ways. A literal
# path only, like the rooms above; what keeps a command away from the
# socket itself is the cage's session doors (api_client._session_doors).
_SERVICE_ROOM_RE = re.compile(
    r"(?<=/)delfin-\d+[\s/,]+agent\.sock\b"
    r"|\.delfin[\s/,]+agent_service[\s/,]+run(?![\w-])",
    re.IGNORECASE)

# Names only an answer file or an answer location carry: ``<id>.answer
# .json`` (and its partial), and a room joined under ``.delfin`` by a
# slash, a Path ``/`` or an ``os.path.join`` comma. A grep for the bare
# word ``terminal_confirmations`` in the code stays allowed.
_ANSWER_FILE_RE = re.compile(
    r"[\w-]\.answer\.(?:json|partial)\b"
    r"|\.delfin[\s/,]+(?:approvals|terminal_confirmations|cli_gate"
    r"|approval_codes)(?![\w-])",
    re.IGNORECASE)


def answer_rooms() -> tuple[Path, ...]:
    """Every directory an answer or the CLI gate state lives in, resolved.

    Output: the redirected module locations and the canonical
    ``~/.delfin/<name>`` for each, de-duplicated. Never raises; a
    location that cannot be computed is left out.
    """
    found: list[Path] = []

    def _add(p) -> None:
        try:
            q = Path(p).expanduser().resolve()
        except Exception:
            return
        if q not in found:
            found.append(q)

    try:
        from . import file_confirm as _fc
        _add(_fc.requests_dir())
    except Exception:
        pass
    try:
        from . import terminal_confirm as _tc
        _add(_tc._PENDING_DIR)
    except Exception:
        pass
    try:
        from . import cli_gate as _cg
        _add(_cg.state_dir())
    except Exception:
        pass
    try:
        # Not an answer room but a code room: the deny globs and the cages
        # keep it unreadable; here it is kept unwritable too.
        from . import approval_code as _ac
        _add(_ac.codes_dir())
    except Exception:
        pass
    try:
        home = Path.home()
        for name in _ROOM_NAMES:
            _add(home / ".delfin" / name)
    except Exception:
        pass
    try:
        # The agent service (delfin-agent tui): a window answers a session's
        # approvals through its socket, found by its run files.
        from .service import rundir as _rd
        _add(_rd.runtime_dir())
        _add(_rd.run_file().parent)
    except Exception:
        pass
    return tuple(found)


def is_answer_room_path(path, base=None) -> bool:
    """True when *path* is an answer location or lies inside one.

    Input: a path as a tool argument names it; relative paths are taken
    against *base*. Both the path as written and its resolved form are
    compared, so a symlink inside the workspace that points into a room
    is caught, and so is a room reached through a symlinked parent.
    Unresolvable input fails closed only when its written form is inside
    a room.
    """
    try:
        p = Path(str(path)).expanduser()
        if not p.is_absolute() and base is not None:
            p = Path(base) / p
        forms = [Path(os.path.abspath(p))]
        try:
            forms.append(p.resolve())
        except Exception:
            pass
    except Exception:
        return False
    for room in answer_rooms():
        for f in forms:
            if f == room or room in f.parents:
                return True
    return False


def _unquoted(cmd: str) -> str:
    """The command with quote characters, backslashes and the characters
    of ``$(...)``, ``${...}`` and backticks removed: ``"$(echo
    approvals)" approve`` reads as the words it runs."""
    return re.sub(r"[`$(){}]", "", _dequoted(cmd))


def _dequoted(cmd: str) -> str:
    """The command with quote characters and backslashes removed."""
    return re.sub(r"[\"'\\]", "", cmd or "")


def answer_command(cmd: str) -> Optional[str]:
    """What in *cmd* answers an approval, or None.

    Input: a shell command, commit-message prose already blanked by the
    caller. Output: a short description of the match. Four shapes:
    the ``approvals approve|deny|answer`` subcommand under any program
    name, ``service stop|restart|end`` (the agent service holds the
    sessions whose approvals a window answers), the Python functions that
    write an answer, and any naming of an answer file or answer location.
    """
    if _ANSWER_COMMAND_RE.search(_unquoted(cmd)):
        return "the command answers a pending approval"
    stop = _SERVICE_STOP_RE.search(_unquoted(cmd))
    if stop and stop.group(1).lower() == "end":
        return ("the command ends sessions of the agent service, which are "
                "the user's, with the approvals they wait on")
    if stop:
        return ("the command stops the agent service, which holds the "
                "user's sessions and the approvals they wait on")
    text = _dequoted(cmd)
    if _ANSWER_ONLY_NAMES_RE.search(text) or (
            _ANSWER_MODULE_RE.search(text) and _ANSWER_CALL_RE.search(text)
            and _IMPORT_RE.search(text)):
        return "the command calls the function that answers a pending approval"
    if _ANSWER_FILE_RE.search(text) or _SERVICE_ROOM_RE.search(text):
        return ("the command names the location where approvals are "
                "answered")
    return None


# -- the operator's actions --------------------------------------------------
#
# ``delfin-agent ops send|escalate|pause`` and the top-level ``delfin-agent
# pause|resume`` are the person's (ops_interface): a message sent with
# them reads as the supervisor's, a pause or a resume decides over another
# session. ops_interface and cmd_pause refuse when the process carries
# ENV_MARKER; a command can strip the marker (``env -u``), so the command
# text is refused as well. The words are read as the shell will run them
# (api_client._command_words: shell bodies, wrappers, assignments,
# quotes), and the program may stand anywhere in a simple command, which
# also covers runners such as ``micromamba run -n env delfin-agent``.

#: The ops subcommands that act. sessions/session/approvals/events/wait read.
_OPS_ACTS = frozenset({"send", "escalate", "pause"})
#: The top-level subcommands that decide over another session.
_TOP_ACTS = frozenset({"pause", "resume"})
_OPERATOR_MODULE_RE = re.compile(r"\b(?:ops_interface|session_pause)\b")
_OPERATOR_CALL_RE = re.compile(r"\b(?:send|escalate|pause|resume)\s*\(")
_CLI_MODULE_RE = re.compile(r"\bdelfin\.agent\.cli\b|\bdelfin\.agent\s+import\s+cli\b")
_CLI_ACT_RE = re.compile(
    r"\bops\b\W{0,8}(?:send|escalate|pause)\b|\bmain\s*\(\s*\[\W*(?:pause|resume)\b")


def _agent_args(words: list) -> Optional[list]:
    """The words after a delfin-agent program in one simple command, or
    None when none is there: ``delfin-agent`` under any directory,
    ``python -m delfin.agent.cli``, or a program held in a variable."""
    for i, w in enumerate(words):
        base = os.path.basename(w)
        if base == "delfin-agent":
            return words[i + 1:]
        if (w == "-m" and i + 1 < len(words)
                and words[i + 1] in ("delfin.agent.cli", "delfin.agent")):
            return words[i + 2:]
    if words and words[0].startswith("$"):
        return words[1:]
    return None


def _acting_subcommand(args: list, xargs: bool) -> str:
    """``ops <verb>`` or ``pause``/``resume`` when *args* run one, else "".
    With xargs on the line, a missing verb or argument comes from stdin."""
    rest = [a for a in args if not a.startswith("-")]
    if not rest:
        return "delfin-agent with its words from stdin" if xargs else ""
    if rest[0] == "ops":
        if len(rest) > 1 and rest[1] in _OPS_ACTS:
            return f"delfin-agent ops {rest[1]}"
        if len(rest) == 1 and xargs:
            return "delfin-agent ops with its verb from stdin"
        return ""
    if rest[0] in _TOP_ACTS and (len(rest) > 1 or xargs):
        return f"delfin-agent {rest[0]}"
    # Accepting a changed MCP tool description is the person's decision.
    if rest[0] == "mcp-pins" and (
            (len(rest) > 1 and rest[1] == "accept") or (len(rest) == 1
                                                        and xargs)):
        return "delfin-agent mcp-pins accept"
    return ""


def operator_command(cmd: str) -> Optional[str]:
    """What in *cmd* acts as the operator, or None.

    Input: a shell command, commit-message prose already blanked by the
    caller. Output: the subcommand it runs. Shapes: ``delfin-agent ops
    send|escalate|pause`` and ``delfin-agent pause|resume`` in every
    spelling ``_command_words`` reads (and with xargs supplying the verb),
    and the Python functions that do the same (``ops_interface.send``,
    ``session_pause.pause``, ``cli.main(["ops", "pause", ...])``)."""
    text = cmd or ""
    try:
        from .api_client import _command_words
        commands = _command_words(text)
    except Exception:
        commands = [part.split() for part in re.split(r"[;&|\n]", text)]
    xargs = bool(re.search(r"(?<![\w.-])xargs(?![\w-])", text))
    for words in commands:
        args = _agent_args(list(words))
        if args is None:
            continue
        found = _acting_subcommand(args, xargs)
        if found:
            return found
    flat = _dequoted(text)
    if _IMPORT_RE.search(flat) or "-m" in flat.split():
        if _OPERATOR_MODULE_RE.search(flat) and _OPERATOR_CALL_RE.search(flat):
            return "the Python function that acts as the operator"
        if _CLI_MODULE_RE.search(flat) and _CLI_ACT_RE.search(flat):
            return "the Python entry point of delfin-agent with an operator action"
    return None


def operator_refusal(what: str) -> str:
    """The message the model reads. Starts with ``blocked:``."""
    return (
        f"blocked: `{what}` acts as the operator. Sending as the operator, "
        "passing an approval on and pausing or resuming a session are the "
        "person's actions, in every permission mode -- not this session's. "
        "To tell another session something, use the session_message tool; "
        "to ask for a pause, tell the user. Reading (`delfin-agent ops "
        "sessions`, `ops events`, `ops approvals`) stays allowed.")


def refusal(what: str) -> str:
    """The message the model reads. Starts with ``blocked:``."""
    return (
        f"blocked: {what}. Only the user answers an approval -- not this "
        "session, not another session, whatever the permission mode. "
        "Show the request to the user (`delfin-agent approvals ls` and "
        "`delfin-agent approvals show <id>` may be read) and ask them to "
        "answer it themselves.")


def inside_agent_session(env=None) -> bool:
    """True when this process was started by an agent session."""
    env = os.environ if env is None else env
    return bool(str(env.get(ENV_MARKER, "") or "").strip())


__all__ = ["ENV_MARKER", "answer_rooms", "is_answer_room_path",
           "answer_command", "refusal", "inside_agent_session",
           "operator_command", "operator_refusal"]
