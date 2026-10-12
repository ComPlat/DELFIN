"""Command forms that are approved by rule instead of by a dialog.

Most approval dialogs of a working session are for a handful of commands
that cannot do harm: checking whether a patch applies, hashing a file,
asking a tool for its version. This module names those forms. The gate
asks it only AFTER the deny list, the secret scan and the refusal ledger
have passed, and the path gates still run afterwards.

A form is an ALLOW LIST, not a filter. A command matches only when it is
exactly one of the spellings below -- same words, same count, nothing
else. The first version of this module took a command apart and vetoed
the options it knew to be dangerous; ``git grep --open-files-in-p=sh x``
went through, because git accepts any unique prefix of a long option and
a list of bad names cannot hold them all. So there are no options here at
all, and no shape for ``git show`` / ``git grep`` (the gate's own
read-only git rule already decides those).

  git apply --check <path>          reads the patch, writes nothing
  git hash-object <path>            prints a hash, stores no object
  git switch -c pr/<name> origin/main   a WRITE, approved by rule: a new
                                    local branch from the remote main
  gh --version | gh version

<path> is one word, relative or absolute, that resolves (links followed)
to something inside the workspace.

"Harmless" means: no option, so no write and no program named in the
command. git still runs what the REPOSITORY configures -- hash-object
and apply --check pass the file through its clean filter, switch runs
the post-checkout hook, the smudge filter and core.fsmonitor -- exactly
as the already auto-allowed git status / git add do. A configuration the
model planted (plain files ``HEAD``, ``config``, ``objects/``, ``refs/``
or a ``.git`` gitfile in the workspace make it a repository) therefore
runs a program here without a dialog. What holds that program is the
shell sandbox every command runs in (bubblewrap, or Landlock: writes
only beneath the workspace roots, network through the egress proxy).
Rejected: requiring ``--no-filters`` on hash-object. apply has no such
option, and git status, reached without a dialog already, runs the same
filter, so it would close nothing.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

# Every character a shape command may contain. Anything the shell would
# interpret -- quotes, a backslash, ``$``, a glob, ``~``, a redirect, a
# separator, a newline, a tab -- is absent, so the words the shell passes
# to the program are exactly the words split here.
_PLAIN = re.compile(r"[A-Za-z0-9_./:@+=,-]+(?: [A-Za-z0-9_./:@+=,-]+)*")

_PR_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}")


def _words(command: str):
    """The command's words, or None when the shell could read it otherwise."""
    if not isinstance(command, str) or not _PLAIN.fullmatch(command):
        return None
    return command.split(" ")


def _inside(arg: str, workspace: str) -> bool:
    """True when the path word resolves, links followed, into the workspace."""
    if not arg or arg.startswith("-"):
        return False
    try:
        root = Path(workspace).resolve(strict=False)
        p = Path(arg)
        if not p.is_absolute():
            p = root / p
        p = p.resolve(strict=False)
    except Exception:
        return False
    return p == root or root in p.parents


def _git_apply_check(w, workspace: str) -> bool:
    return (len(w) == 4 and w[:3] == ["git", "apply", "--check"]
            and _inside(w[3], workspace))


def _git_hash_object(w, workspace: str) -> bool:
    return (len(w) == 3 and w[:2] == ["git", "hash-object"]
            and _inside(w[2], workspace))


def _git_switch_pr(w, workspace: str) -> bool:
    if len(w) != 5 or w[:3] != ["git", "switch", "-c"] or w[4] != "origin/main":
        return False
    if not w[3].startswith("pr/"):
        return False
    name = w[3][3:]
    return bool(_PR_NAME.fullmatch(name)) and ".." not in name \
        and not name.endswith((".", ".lock"))


def _gh_version(w, workspace: str) -> bool:
    return w in (["gh", "--version"], ["gh", "version"])


_SHAPES = (
    ("git_apply_check", _git_apply_check),
    ("git_hash_object", _git_hash_object),
    ("git_switch_pr", _git_switch_pr),
    ("gh_version", _gh_version),
)

#: Ordered shape names.
SHAPE_NAMES: tuple[str, ...] = tuple(name for name, _ in _SHAPES)

#: The spelling each shape approves, as the settings box shows it.
SHAPE_SPELLINGS: dict[str, str] = {
    "git_apply_check": "git apply --check <path>",
    "git_hash_object": "git hash-object <path>",
    "git_switch_pr": "git switch -c pr/<name> origin/main",
    "gh_version": "gh --version / gh version",
}


def matches(command: str, workspace: str) -> str:
    """The name of the shape *command* is, or ``""``."""
    return matches_with(command, workspace, None)


def matches_with(command: str, workspace: str,
                 enabled: Iterable[str] | None) -> str:
    """As :func:`matches`, but only the shapes in *enabled* may fire.
    ``None`` means every shape; an empty collection means none."""
    allowed = set(SHAPE_NAMES) if enabled is None else set(enabled)
    words = _words(command)
    if words is None:
        return ""
    for name, fn in _SHAPES:
        if name in allowed and fn(words, workspace):
            return name
    return ""


# ── stored settings contract ────────────────────────────────────────────────
#
# The gate (api_client.py, protected) consults ``agent.read_only_shapes`` from
# the user settings payload before the confirm dialog; the settings tab
# (tab_settings.py, owned) renders and persists that same block. The keys, the
# defaults and the normalisation live HERE so both sides read one definition:
#
#     agent.read_only_shapes = {
#         "enabled": bool,          # master switch, default True
#         "enabled_shapes": [name], # active shape names, default = SHAPE_NAMES
#     }
#
# Normalisation on load: unknown shape names are dropped, the two known keys
# are always present, and a partial save can never corrupt the block.

SETTINGS_KEY = "read_only_shapes"


def read_only_shapes_defaults() -> dict:
    """The canonical on-disk settings for the read-only-shape gate."""
    return {
        "enabled": True,
        "enabled_shapes": list(SHAPE_NAMES),
    }


def read_only_shapes_from_payload(payload: dict) -> dict:
    """Return the normalised ``read_only_shapes`` block from *payload*.

    *payload* is a user-settings dict; the block lives under
    ``payload['agent']['read_only_shapes']``. Missing payload, missing block or
    malformed values fall back to the defaults (master switch ON, every shape
    enabled). Unknown shape names are dropped so a stale name cannot silently
    widen the gate, and an empty list switches every shape off.
    """
    unknown = set(SHAPE_NAMES)
    agent = (payload or {}).get("agent") or {}
    block = agent.get(SETTINGS_KEY) or {}
    if not isinstance(block, dict):
        return read_only_shapes_defaults()
    enabled = block.get("enabled", True)
    if not isinstance(enabled, bool):
        enabled = True
    wanted = block.get("enabled_shapes")
    if not isinstance(wanted, list):
        shapes = list(SHAPE_NAMES)
    else:
        # An empty list is the user unticking every shape: none, not all.
        shapes = [n for n in wanted if n in unknown]
    return {"enabled": enabled, "enabled_shapes": shapes}


def read_only_shapes_into_payload(payload: dict, settings: dict) -> dict:
    """Write the normalised *settings* block back into *payload* (mutating)
    without disturbing other ``agent`` keys, and return *payload*."""
    payload.setdefault("agent", {})
    payload["agent"][SETTINGS_KEY] = read_only_shapes_from_payload(payload)
    stored = payload["agent"][SETTINGS_KEY]
    stored["enabled"] = bool(settings.get("enabled", True))
    stored["enabled_shapes"] = [
        n for n in settings.get("enabled_shapes", SHAPE_NAMES)
        if n in set(SHAPE_NAMES)
    ]
    return payload
