"""The allowed form for a command refused by a bash deny pattern.

``suggest(command, reason)`` takes the refused command line and the deny
pattern that refused it and returns a short note naming the way forward, or
``""`` when there is none.

Semantics:

* The note is chosen by the PATTERN that refused the command, never by the
  shape of the command. Keyed on the command, the first cut answered
  ``sudo cat /etc/shadow 2>&1`` (refused by ``\\bsudo\\b``) with "Run
  `sudo cat /etc/shadow` instead" and ``rm -rf ~`` with a ``mv`` of the
  home directory.
* A backticked command in a note is re-checked by the caller
  (``_DocToolExecutor._run_permission_gate`` drops a note whose command the
  deny list refuses again).
* A path is echoed into a note only when it is a plain relative path (no
  ``~``, ``$``, ``..``, quotes or leading ``/``). Anything else gets a note
  without a command.

Patterns handled: the default ``rm`` deletes, and the shapes an installation
may add to ``agent.bash_deny_patterns`` (``||``, a stderr redirect, ``;``,
pytest, git fetch). A pattern with a route in ``api_client._DENY_HINTS`` is
answered there; the caller asks this module only when that table is silent.
"""

from __future__ import annotations

import re

#: A plain relative path: no shell syntax, no home, variable or root anchor.
_PLAIN_RELATIVE = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_./+-]*")

#: The stderr redirect forms a reduced command drops.
_STDERR_REDIRECT = re.compile(r"\s*2>(?:&\d+|\s*/dev/null)")


def _rm_target(command: str) -> str:
    """The single operand of a command that is one plain ``rm``, or ""."""
    words = command.split()
    if not words or words[0] != "rm":
        return ""
    operands = [w for w in words[1:] if not w.startswith("-")]
    if len(operands) != 1:
        return ""
    target = operands[0].rstrip("/")
    if (not _PLAIN_RELATIVE.fullmatch(target)
            or ".." in target.split("/")):
        return ""
    return target


def _hint_delete(command: str) -> str:
    target = _rm_target(command)
    if target:
        return (f"To set it aside without destroying it: "
                f"`mv {target} {target}.bak`.")
    return ("Name the files you mean to remove in your answer; the user "
            "deletes them.")


def _hint_stderr(command: str) -> str:
    reduced = re.sub(r"\s+", " ", _STDERR_REDIRECT.sub("", command)).strip()
    if not reduced or reduced == command:
        return ""
    return f"Drop the stderr redirect: `{reduced}`."


def suggest(command: str, reason: str = "") -> str:
    """The note for ``command`` refused by deny pattern ``reason``, or ""."""
    cmd = (command or "").strip()
    pattern = reason or ""
    if not cmd or not pattern:
        return ""
    if pattern.startswith((r"\brm\b", r"\brm\s")):
        return _hint_delete(cmd)
    if "pytest" in pattern:
        return ("pytest is refused by this installation's deny list; ask the "
                "user how tests are run here.")
    if "fetch" in pattern:
        return ("This installation denies fetching from a remote in the "
                "shell; say what you need from it and the user fetches.")
    if r"\|\|" in pattern or pattern == "||":
        return ("Run the fallback as its own bash call and read the first "
                "call's exit code.")
    if "2>" in pattern:
        return _hint_stderr(cmd)
    if pattern.strip() in (";", r"\;", r"\s*;\s*"):
        return "Run each part of the ';' chain as its own bash call."
    return ""
