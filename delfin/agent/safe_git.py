"""Git as DELFIN's own code runs it on a repository it does not control.

DELFIN reads repository state -- branch, status, recent commits, a diff --
from outside the sandbox: for the launch banner, the system prompt, the
evidence fingerprint, project introspection, held-branch recovery and
worktrees. The repository is the session's workspace, which a command in
the sandbox can write. Git runs programs named in the repository's own
configuration and attributes (fsmonitor, clean/smudge/process filters,
hooks, pager, textconv, external diff, ssh command, credential helper,
...). Run plainly, such a read would execute whatever the repository names,
with the user's rights and outside every sandbox layer.

``run`` is the one way DELFIN's code runs git on such a repository:

- fixed ``-c`` overrides after git's own name and before the subcommand,
  so no caller argument can weaken them: no fsmonitor, no hooks, no pager,
  no editor, no askpass, no external diff, no signature program;
- every program-naming key defined in the repository's own config (scope
  ``local`` or ``worktree``, includes followed) is overridden with an
  empty value; a filter is also marked not required, so a blanked filter
  passes the content through instead of failing the command. The user's
  global and the system configuration are trusted and left alone;
- a key whose name cannot be overridden on the command line (a
  subsection containing ``=``) refuses the call instead of running it;
- ``diff``/``log``/``show`` get ``--no-textconv --no-ext-diff``,
  ``blame``/``grep`` ``--no-textconv``; ``status``/``diff`` ignore submodules (a submodule is
  another repository with its own configuration);
- ``--no-optional-locks``: a read never takes the session's index lock;
- the environment carries no inherited ``GIT_*`` variable, no terminal
  prompt and no optional locks.

Never raises: a failure is a ``CompletedProcess`` with a non-zero return
code and the reason in ``stderr``, the shape every caller already handles.

Rejected: reimplementing ``status`` from the index stat data (measured
equivalent to ``status`` with the overrides above, at the cost of a second
implementation to keep in step); dropping the system configuration
(``GIT_CONFIG_NOSYSTEM``), which on shared clusters carries
``safe.directory`` and would stop git from reading shared repositories.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from typing import Optional, Sequence

#: Keys whose value git runs as a program (or, for include paths, reads as
#: more configuration), matched on the lower-cased key: git compares
#: section and variable names case-insensitively. Subsection names
#: (``<driver>``, ``<url>``, ``<format>``) are case-sensitive in git; the
#: match lower-cases them too, which only widens it. The bash gate
#: (api_client) uses the same expression for ``git config`` writes.
PROGRAM_KEY_RE = re.compile(
    r"(?:core\.(?:hookspath|fsmonitor|sshcommand|pager|editor|askpass|"
    r"gitproxy|alternaterefscommand)"
    r"|sequence\.editor|gpg\.(?:.+\.)?program|credential\.(?:.+\.)?helper"
    r"|diff\.external|diff\..+\.(?:textconv|command)"
    r"|filter\..+\.(?:clean|smudge|process)|merge\..+\.driver"
    r"|(?:difftool|mergetool)\..+\.cmd|interactive\.difffilter"
    r"|remote\..+\.(?:uploadpack|receivepack)|submodule\..+\.update"
    r"|pager\..+|alias\..+|include\.path|includeif\..+\.path)")

#: Overrides that hold whatever the repository or the caller says.
FIXED_OVERRIDES: tuple[str, ...] = (
    "core.fsmonitor=false",
    "core.hooksPath=/dev/null",
    "core.pager=cat",
    "core.editor=:",
    "sequence.editor=:",
    "core.askPass=",
    "core.alternateRefsCommand=",
    "diff.external=",
    "log.showSignature=false",
    "protocol.ext.allow=never",
    "submodule.recurse=false",
    "status.submoduleSummary=false",
    "diff.ignoreSubmodules=all",
    "gc.auto=0",
    "maintenance.auto=false",
)

#: Repository scopes whose program keys are neutralised.
_UNTRUSTED_SCOPES = frozenset({"local", "worktree"})

#: Per subcommand, the flags that keep it from running a textconv or an
#: external diff program (git grep has no --no-ext-diff).
_NO_PROGRAM_FLAGS = {
    "diff": ("--no-textconv", "--no-ext-diff"),
    "log": ("--no-textconv", "--no-ext-diff"),
    "show": ("--no-textconv", "--no-ext-diff"),
    "blame": ("--no-textconv",),
    "grep": ("--no-textconv",),
}
_SUBMODULE_AWARE = frozenset({"status", "diff"})

_DEFAULT_TIMEOUT_S = 30.0


class Refused(Exception):
    """The repository names a program in a way an override cannot reach."""


def _env(extra: Optional[dict] = None) -> dict:
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update({
        "GIT_TERMINAL_PROMPT": "0",
        "GIT_OPTIONAL_LOCKS": "0",
        "GIT_PAGER": "cat",
        "GIT_EDITOR": ":",
        "GIT_ASKPASS": "",
        "SSH_ASKPASS": "",
    })
    if extra:
        env.update(extra)
    return env


def _base(cwd: str) -> list[str]:
    argv = ["git", "-C", str(cwd), "--no-pager", "--no-optional-locks"]
    for kv in FIXED_OVERRIDES:
        argv += ["-c", kv]
    return argv


def repository_program_overrides(cwd: str, *, timeout: float = 10.0) -> list[str]:
    """``key=`` overrides for every program-naming key the repository at
    *cwd* defines in its own (local or worktree) configuration, includes
    followed. Raises ``Refused`` for a key the command line cannot name."""
    try:
        proc = subprocess.run(
            _base(cwd) + ["config", "--null", "--list", "--show-scope",
                          "--includes"],
            capture_output=True, timeout=timeout, env=_env(), check=False)
    except (OSError, subprocess.TimeoutExpired):
        return []
    if proc.returncode != 0:
        return []
    out: list[str] = []
    seen: set[str] = set()
    # --null --show-scope: "<scope>\0<key>\n<value>\0" per entry.
    parts = proc.stdout.split(b"\0")
    for i in range(0, len(parts) - 1, 2):
        scope = parts[i].decode("utf-8", "replace").strip()
        entry = parts[i + 1].decode("utf-8", "replace")
        key = entry.split("\n", 1)[0].strip()
        if scope not in _UNTRUSTED_SCOPES or not key:
            continue
        if not PROGRAM_KEY_RE.fullmatch(key.lower()):
            continue
        if "=" in key:
            raise Refused(f"the repository's configuration names a program "
                          f"under a key DELFIN cannot override: {key[:80]!r}")
        if key in seen:
            continue
        seen.add(key)
        out.append(f"{key}=")
        low = key.lower()
        if low.startswith("filter.") and low.rsplit(".", 1)[-1] in (
                "clean", "smudge", "process"):
            name = key[len("filter."):].rsplit(".", 1)[0]
            req = f"filter.{name}.required=false"
            if req not in seen:
                seen.add(req)
                out.append(req)
    return out


def argv_for(args: Sequence[str], cwd: str, *, timeout: float = 10.0) -> list[str]:
    """The full command line ``run`` executes. Raises ``Refused``."""
    args = [str(a) for a in args]
    argv = _base(cwd)
    for kv in repository_program_overrides(cwd, timeout=timeout):
        argv += ["-c", kv]
    if not args:
        return argv
    sub, rest = args[0], args[1:]
    argv.append(sub)
    argv += list(_NO_PROGRAM_FLAGS.get(sub, ()))
    if sub in _SUBMODULE_AWARE:
        argv.append("--ignore-submodules=all")
    return argv + rest


def run(args: Sequence[str], cwd, *, timeout: float = _DEFAULT_TIMEOUT_S,
        text: bool = True, input=None, env: Optional[dict] = None
        ) -> subprocess.CompletedProcess:
    """Run ``git <args>`` in *cwd* as described in the module docstring."""
    try:
        argv = argv_for(args, str(cwd), timeout=min(timeout, 10.0))
    except Refused as exc:
        return subprocess.CompletedProcess(list(args), 128, "" if text else b"",
                                           str(exc) if text else str(exc).encode())
    try:
        return subprocess.run(argv, capture_output=True, text=text, input=input,
                              timeout=timeout, env=_env(env), check=False)
    except FileNotFoundError as exc:
        msg = f"git is not installed: {exc}"
        return subprocess.CompletedProcess(argv, 127, "" if text else b"",
                                           msg if text else msg.encode())
    except (OSError, subprocess.TimeoutExpired) as exc:
        return subprocess.CompletedProcess(argv, 1, "" if text else b"",
                                           str(exc) if text else str(exc).encode())


def out(cwd, *args: str, timeout: float = _DEFAULT_TIMEOUT_S) -> str:
    """stdout of a successful run, "" otherwise."""
    r = run(args, cwd, timeout=timeout)
    return r.stdout if r.returncode == 0 else ""


def is_repository(path) -> bool:
    try:
        return Path(path).exists() and run(
            ["rev-parse", "--git-dir"], path, timeout=5).returncode == 0
    except Exception:
        return False
