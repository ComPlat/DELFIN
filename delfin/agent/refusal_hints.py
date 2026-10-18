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


def _checkout_shape(command: str) -> str:
    """Which checkout shape a git line is: "paths" (``git checkout --
    paths`` or ``git checkout <ref> -- paths``), "branch" (a plain branch
    switch), or "" (anything else, e.g. git log that merely mentions the
    word). The split on `` -- `` is what separates a destructive checkout
    from a read; the word "checkout" alone never does."""
    words = command.split()
    if not words or words[0] != "git" or "checkout" not in words[1:2]:
        if not (words and words[0] == "git" and len(words) > 1
                and words[1] == "checkout"):
            return ""
    if " -- " in command:
        return "paths"
    if "--" in words[2:]:
        return "paths"
    return "branch"


def _hint_checkout(command: str) -> str:
    shape = _checkout_shape(command)
    if shape == "paths":
        return ("A checkout with -- discards or imports files around the "
                "write gate: commit your uncommitted work first, and bring "
                "another branch's file in by merge, not by checking it "
                "out.")
    return ""


def _hint_amend() -> str:
    return ("A commit others may already have is never amended; make a new "
            "commit with the correction instead.")


def _hint_revert_merge() -> str:
    return ("A merge commit is not reverted in place; undo the specific "
            "change in a new commit so the branch's history stays intact.")


def _hint_script_write() -> str:
    return ("A script that writes a file bypasses the write gate; make the "
            "change with edit_file so the write scope holds.")


def _hint_repeat_refusal() -> str:
    return ("A refused action is not retried in another spelling; ask the "
            "user how to proceed instead.")


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
    # Package-B2 additions: keyed on the deny PATTERN (the gate's reason),
    # never on the command's shape -- the shape only refines the checkout
    # case, after the pattern named checkout at all.
    pattern_text = pattern.lower()
    if "checkout" in pattern_text or "restore" in pattern_text:
        return _hint_checkout(cmd)
    if "amend" in pattern_text:
        return _hint_amend()
    if "revert" in pattern_text:
        return _hint_revert_merge()
    if ("script" in pattern_text
            and ("> " in cmd or "-c " in cmd)):
        return _hint_script_write()
    if ("repeat" in pattern_text or "again" in pattern_text):
        return _hint_repeat_refusal()
    if "push" in pattern_text and "force" in cmd:
        return _hint_repeat_refusal()
    return ""
