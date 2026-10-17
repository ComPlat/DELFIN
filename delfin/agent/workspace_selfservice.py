"""Self-service shapes: commands that may run without a question inside the
agent's own workspace, and nowhere else.

The gate's bash auto-allow list covers read-only commands and a handful of
shape-free write forms. Everything else asks a person. A working session
reaches, again and again, for a small set of commands that only write inside
the workspace and only touch known, pinned inputs: creating a venv, installing
*hashed* requirements into it, switching to an own branch, asking a known tool
for its version. Those are the shapes this module approves by rule instead of
by a dialog, but ONLY when the caller has already confined execution to the
own workspace (the sandbox/cage). When it is off, or the command is not
exactly one of these shapes, classify returns None and the gate decides as it
does today -- an ask or a deny is never turned into an allow here except for
the exact rated shapes below, inside the exact workspace given.

This is a *structure* matcher, not a word list. It reuses the gate's own
tokenizer (``_split_shell_segments``) and write-target extractor
(``_bash_write_targets``) from api_client -- there is deliberately no second
parser here. A command that does not parse into exactly the tokens of one
rated shape is refused by returning None.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Callable, Literal, Optional

# The single verdict this classifier may emit. None means "the gate decides".
AUTO = "auto"
Verdict = Optional[Literal["auto"]]


# ---------------------------------------------------------------------------
# Primitive reuse -- the gate's own tokenizer and write-target extractor. The
# caller (the gate in api_client) has already imported this module; importing
# these two pure functions back from it is the intended wiring, and keeps the
# single parser that has authority over what a command writes.
# ---------------------------------------------------------------------------
def _segments(cmd: str) -> list[str]:
    from delfin.agent import api_client as _ac
    return list(_ac._split_shell_segments(cmd))


def _write_targets(cmd: str) -> list[str]:
    from delfin.agent import api_client as _ac
    return list(_ac._bash_write_targets(cmd))


# ---------------------------------------------------------------------------
# Workspace containment. A write target is only ever allowed when it resolves,
# symlinks followed, strictly inside the given workspace.
# ---------------------------------------------------------------------------
def _inside(path: str, workspace: str) -> bool:
    """True when `path` (absolute or relative) resolves inside `workspace`.

    A relative path is resolved against the WORKSPACE, not the process cwd:
    the gate runs the command with the workspace as its cwd, so
    ``.venv/bin/pip`` means ``<workspace>/.venv/bin/pip``.

    ``resolve(strict=False)`` folds ``..`` lexically and follows whatever
    symlinks exist on disk, so a symlink inside the workspace that points at
    a directory outside it resolves to the outside path and is rejected --
    exactly the escape the sandbox would also stop, made cheap to test.
    """
    try:
        root = Path(workspace).resolve(strict=False)
        p = Path(os.path.expanduser(path))
        # Join BEFORE resolving: Path.resolve() folds a relative path against
        # the process cwd, which is not the workspace -- measured:
        # _inside(".venv", wd) was False because ".venv" resolved against the
        # gate's cwd and became absolute before the workspace join.
        if not p.is_absolute():
            p = root / p
        p = p.resolve(strict=False)
    except (OSError, ValueError):
        return False
    return root in p.parents  # strictly inside; the workspace root itself is not a target


def _all_targets_inside(cmd: str, workspace: str) -> bool:
    """False when the write-target extractor yields a target outside the
    workspace, or yields no targets at all for a shape that is pinned by its
    write-target check (S1 -- the pip command line must produce at least the
    requirements path it names, or the extractor missed something and the
    write gate gets to see the command instead)."""
    targets = _write_targets(cmd)
    if not targets:
        return False
    return all(_inside(t, workspace) for t in targets)


def _no_target_outside(cmd: str, workspace: str) -> bool:
    """True when the extractor yields NO target that resolves outside the
    workspace. Used by the shapes whose pinned command forms produce no
    targets at all (git switch, <tool> --version): the exact word forms are
    pinned by the shape matcher itself, so anything the extractor WOULD see
    (a redirect, a dest option) appearing as an outside target still refuses
    here, while a target-free pinned form is allowed."""
    return all(_inside(t, workspace) for t in _write_targets(cmd))


# ---------------------------------------------------------------------------
# Token helpers for the exact shapes.
# ---------------------------------------------------------------------------
_ENV_ASSIGN_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")


def _plain_words(cmd: str) -> Optional[list[str]]:
    """shlex the single segment into words, or None when the command is not
    exactly one chain-free segment whose first word is not an env assignment.

    An env assignment (``FOO=bar pip …``) is a shape the gate must still
    decide on -- a model could smuggle a flag through it -- so any command
    opening with one returns None.
    """
    if len(_segments(cmd)) != 1:
        return None
    import shlex
    try:
        words = shlex.split(cmd)
    except ValueError:
        return None
    if not words or _ENV_ASSIGN_RE.match(words[0]):
        return None
    return words


def _venv_prefix_of(words: list[str]) -> Optional[str]:
    """The virtualenv prefix when the command runs through an in-workspace
    venv: ``.venv/bin/pip …`` or ``.venv/bin/python -m pip …``. Returns the
    prefix string (e.g. ``.venv``) or None when the interpreter is not a
    venv-prefixed path.
    """
    head = words[0]
    if "/bin/pip" in head or "/bin/python" in head or "/bin/python3" in head:
        marker = "/bin/"
        idx = head.rfind(marker)
        return head[:idx]
    return None


# ---------------------------------------------------------------------------
# Shape S1 -- pip install --require-hashes into an in-workspace venv.
#
#   <venv>/bin/pip install -r <file> --require-hashes
#   <venv>/bin/python -m pip install -r <file> --require-hashes
#
# Approved only when every condition holds:
#   * exactly one segment (no &&, ||, ;, |, newline) and a venv-prefixed
#     interpreter -- the venv itself must resolve inside the workspace;
#   * a requirements file ``-r <file>`` that itself resolves inside the
#     workspace;
#   * ``--require-hashes`` present, so a remote index cannot substitute
#     unvetted wheels;
#   * NO index override (``--index-url``/``-i``/``--extra-index-url``), NO
#     editable VCS source (``-e``/``--editable``), NO ``--user`` install;
#   * every write target the gate's own extractor finds is inside the
#     workspace.
# ---------------------------------------------------------------------------
_INDEX_OVERRIDE = ("--index-url", "--extra-index-url", "--extra-index", "-i")
_EDITABLE = ("-e", "--editable")
_DANGEROUS_PIP_OPTS = _INDEX_OVERRIDE + _EDITABLE + ("--user", "--prefix", "--root")


def _matches_dangerous_opt(w: str) -> bool:
    """True when the word is a dangerous pip option in ANY spelling pip
    accepts: standalone (``-e``), attached value (``-e.``), or equals form
    (``--index-url=https://...``). pip's own parser accepts all of these
    (measured in s25's adversarial cases F1/F2), so containment cannot
    match only the standalone word."""
    for opt in _DANGEROUS_PIP_OPTS:
        if w == opt:
            return True
        if w.startswith(opt + "="):
            return True
        # attached value: only for short options (-e. / -e/tmp/evil); for
        # long options an attached form without '=' is not one option.
        if not opt.startswith("--") and w.startswith(opt) and len(w) > len(opt):
            return True
    return False


def _pip_install_r(cmd: str, workspace: str) -> bool:
    words = _plain_words(cmd)
    if words is None:
        return False
    prefix = _venv_prefix_of(words)
    if prefix is None or not _inside(prefix, workspace):
        return False
    # The interpreter form: <venv>/bin/python -m pip install …
    idx = 0
    if words[0].endswith("/python") or words[0].endswith("/python3"):
        if len(words) < 4 or words[1:3] != ["-m", "pip"]:
            return False
        idx = 3
    if idx == 0 and not words[0].endswith("/pip"):
        return False
    if idx == 0:
        idx = 1
    if idx < len(words) and words[idx] != "install":
        return False
    idx += 1  # past "install"
    req_values: list[str] = []
    consumed_flag_word = False
    i = idx
    while i < len(words):
        w = words[i]
        if _matches_dangerous_opt(w):
            return False
        if w in ("-r", "--requirement"):
            # -r consumes the next word as its value. If there is no next
            # word, or the next word looks like another option (starts with
            # '-'), the option is malformed -- pip would abort or swallow
            # it (s25 F6/F7); the pinned form cannot be verified -> None.
            if i + 1 >= len(words) or words[i + 1].startswith("-"):
                return False
            req_values.append(words[i + 1])
            i += 2
            continue
        if w == "--require-hashes":
            consumed_flag_word = True
        i += 1
    if not req_values or any(not _inside(r, workspace) for r in req_values):
        return False
    if not consumed_flag_word:
        return False
    # The pinned word form produces no targets in the gate's extractor
    # (measured: pip install with -r and --require-hashes yields []); the
    # extractor still gets a vote on anything it DOES see -- a redirect or
    # a dest option the word pinning missed must not point outside.
    return _no_target_outside(cmd, workspace)


# ---------------------------------------------------------------------------
# Shape S2 -- git switch / checkout to an own, non-foreign branch.
#
#   git switch <branch>
#   git switch -c <branch>
#   git checkout -b <branch>
#
# git cannot tell which branch "belongs to" this session, so the branch name
# is pinned to a SHAPE, not matched against a session name. A branch is
# approved when it is a plain local refname that is not ``main``/``master``/
# ``HEAD``, not an ``origin/*`` remote-tracking ref, and not another working
# session's branch (``agent/sNN-*``, bare ``sNN-*``). A switch through any
# other spelling (``git checkout <name>``, ``git switch <name> <start>``)
# stays with the gate.
# ---------------------------------------------------------------------------
_SAFE_BRANCH_RE = re.compile(r"^(?!\.)[A-Za-z0-9][A-Za-z0-9._-]{0,200}$")
_BAD_BRANCH = {"main", "master", "HEAD"}


def _is_own_branch_switch(words: list[str]) -> Optional[str]:
    """The safe branch name when `words` is one of the accepted switch forms
    and the branch is an own, non-foreign branch; else None."""
    if words[:2] == ["git", "switch"]:
        if len(words) == 3:
            branch = words[2]
        elif len(words) == 4 and words[2] == "-c":
            branch = words[3]
        else:
            return None
    elif words[:3] == ["git", "checkout", "-b"] and len(words) == 4:
        branch = words[3]
    else:
        return None
    if branch in _BAD_BRANCH:
        return None
    if branch.startswith("origin/"):
        return None
    # A refname component may not be '..' (git check-ref-format); a branch
    # containing '..' would resolve to another ref or fail -- s25 F3 pins
    # it to None.
    if ".." in branch:
        return None
    if not _SAFE_BRANCH_RE.match(branch):
        return None
    if branch.startswith("agent/") or re.match(r"^s\d+[-_]", branch):
        return None  # another working session's branch
    return branch


def _git_switch_own(cmd: str, workspace: str) -> bool:
    words = _plain_words(cmd)
    if words is None:
        return False
    branch = _is_own_branch_switch(words)
    if branch is None:
        return False
    # The switch and any write it triggers must stay in the workspace.
    return _no_target_outside(cmd, workspace)


# ---------------------------------------------------------------------------
# Shape S3 -- <known-tool> --version. No other words, no redirects, no options
# that could read or write anything.
#
# The tools on the allow-list are exactly those the gate already treats as
# read-only elsewhere (auto-allow logs / read-only-shapes). A tool that is not
# on the list returns None -- unknown tools are a hole and stay with the gate.
# ---------------------------------------------------------------------------
_KNOWN_VERSION_TOOLS = frozenset(
    ["gh", "git", "python", "python3", "pip", "pip3", "bash", "sh", "zsh",
     "make", "cmake", "ninja", "meson", "tar", "unzip", "cmake3"]
)


def _tool_version(cmd: str) -> bool:
    words = _plain_words(cmd)
    if words is None:
        return False
    if len(words) != 2:
        return False
    tool, flag = words
    # Only real version flags: `gh version` is not gh's version command
    # (s25 F4) and `-V` is not bash's. A bare `git version` is git's.
    if flag in ("version", "-V") and tool != "git":
        return False
    if flag not in ("--version", "-V", "version"):
        return False
    if tool not in _KNOWN_VERSION_TOOLS:
        return False
    # A bare `git version` is fine; `git --version` prints a version string.
    if tool == "git" and flag == "version":
        return True
    targets = _write_targets(cmd)
    return not targets  # reads nothing, writes nothing


# ---------------------------------------------------------------------------
# The classifier. Returns "auto" for exactly the rated shapes confined to the
# given workspace; None for everything else (the gate decides as today).
#
# The sandbox-off guard is NOT this function's: it takes no sandbox state.
# The caller (the gate in api_client) consults classify ONLY after it has
# established that the command runs inside the sandbox/cage confined to the
# own workspace and that the deny/secret/path gates have passed.
# ---------------------------------------------------------------------------
def classify(cmd: str, workspace: str) -> Verdict:
    """Verdict for a bash command confined to the own ``workspace``.

    ``auto`` for the exact, unparse-safe self-service shapes (pip install
    ``--require-hashes`` from an in-workspace venv; ``git switch`` to an own
    branch; ``<known-tool> --version``); ``None`` for everything else --
    an ask or deny of the existing gate is never made an allow here.
    """
    if not cmd or not str(cmd).strip():
        return None
    cmd_s = str(cmd)
    if _tool_version(cmd_s):
        return AUTO
    if _git_switch_own(cmd_s, workspace):
        return AUTO
    if _pip_install_r(cmd_s, workspace):
        return AUTO
    return None


# ── stored settings contract ────────────────────────────────────────────────
#
# The gate (api_client.py, protected) consults ``agent.workspace_selfservice``
# from the user settings payload at the caged bash admission; the person can
# switch the shapes off there. The keys, the defaults and the normalisation
# live HERE so both sides read one definition:
#
#     agent.workspace_selfservice = {
#         "enabled": bool,   # master switch, default True
#     }
#
# Normalisation on load: a malformed value falls back to the default, and a
# partial save can never corrupt the block.

SETTINGS_KEY = "workspace_selfservice"


def selfservice_defaults() -> dict:
    """The canonical on-disk settings for the workspace self-service gate."""
    return {"enabled": True}


def settings_from_payload(payload: dict) -> dict:
    """Return the normalised ``workspace_selfservice`` block from *payload*.

    *payload* is a user-settings dict; the block lives under
    ``payload['agent']['workspace_selfservice']``. Missing payload, missing
    block or malformed values fall back to the default (master switch ON).
    """
    agent = (payload or {}).get("agent") or {}
    block = agent.get(SETTINGS_KEY) or {}
    if not isinstance(block, dict):
        return selfservice_defaults()
    enabled = block.get("enabled", True)
    if not isinstance(enabled, bool):
        enabled = True
    return {"enabled": enabled}


def settings_into_payload(payload: dict, settings: dict) -> dict:
    """Write the normalised *settings* block back into *payload* (mutating)
    without disturbing other ``agent`` keys, and return *payload*."""
    payload.setdefault("agent", {})
    stored = dict(selfservice_defaults())
    stored["enabled"] = bool(settings.get("enabled", True))
    payload["agent"][SETTINGS_KEY] = stored
    return payload


__all__ = ["classify", "AUTO", "Verdict",
           "SETTINGS_KEY", "selfservice_defaults",
           "settings_from_payload", "settings_into_payload"]
