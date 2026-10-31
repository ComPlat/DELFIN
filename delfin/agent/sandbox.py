"""Runner for the dashboard agent's approved CLI-backend ``Bash`` calls.

Input: a command the person approved in the dashboard (Layer 4) and the
repository directory. Output: a ``RunResult``. In every mode the command
first passes the shared-checkout rule (``shared_checkout.refusal``) and then
runs as ``bash -c`` with the scrubbed shell environment (no provider keys),
in its own session and process group, ended at exit and at the timeout
(``contained_run.run``).

Layers, as the code applies them:

  1. Command allow-list (every mode except ``off``) -- each pipeline
     segment is judged by ``_check_argv``: a tool with a sub-command list
     (``_ALLOW_SUBCOMMAND``: ``git``, ``pip``) is judged by that list in
     every spelling (``pip3``, ``pip3.N``, ``<python> -m pip``), and a
     sub-command with a writing form passes only in its read forms
     (``_read_only_form``: ``git config``, ``git remote``, ``pip config``);
     ``env``, ``nohup`` and ``xargs`` (``_RUNS_ITS_ARGUMENTS``) are judged by
     the command they run (``env -S`` is refused); every other first word
     must be in ``_ALLOW_FIRST_TOKEN``, which does not contain ``pip``. No
     segment may match ``_DENY_SUBSTRINGS``, ``unsafe_install_refusal``
     must pass, and a tree walk rooted at a whole file system is refused.
     This is a text check, not a sandbox: the list admits ``python``,
     ``awk``, ``sed``, ``find`` and ``tee``, each of which can run any
     program or write any file the account can, and a deny substring is
     matched as text only.
  2. Confinement, by ``DELFIN_AGENT_SANDBOX`` (``detect_config``); what is
     applied on this host is ``confinement()``, the same value the runner
     builds its argv from and the dashboard prints
     (``mode=..., fs=..., net=...``):

     * ``bwrap``: read-only root file system; ``/tmp`` and ``/root`` are
       fresh tmpfs; ``$HOME`` stays visible read-only, with the entries of
       ``_HOME_SECRET_DIRS`` masked (tmpfs for a directory, ``/dev/null``
       for a file) and the approval-code directory masked; the repository
       bound read-write; PID, IPC and UTS namespaces unshared; network
       unshared unless ``DELFIN_AGENT_SANDBOX_NETWORK=1``.
     * ``firejail``: ``--noprofile`` with a private ``/tmp`` and ``/dev``,
       the repository whitelisted, ``/etc`` and ``/usr`` read-only,
       ``--noroot``, all capabilities dropped; ``--net=none`` unless
       ``DELFIN_AGENT_SANDBOX_NETWORK=1``.
     * ``allowlist``: no namespace. Where the kernel has Landlock the
       command runs under ``landlock_exec`` (writes only in the repository
       and a private temp directory, credential paths hidden), and the
       network is denied unless ``DELFIN_AGENT_SANDBOX_NETWORK=1``
       (``net=on``). The deny is ``net=off`` where the seccomp socket
       guard is available (every connect judged, UDP included),
       ``net=off-tcp`` where only Landlock ABI >= 4 holds it (TCP connects
       refused, UDP not), and ``net=open`` with a notice where neither
       exists. A restricted network is passed with ``--strict 1``: if the
       helper finds a layer missing that the probe reported, it refuses
       the command (exit 126) instead of running it with the network open.
       Without Landlock the command runs with the file system and network
       of the account (``fs=none, net=open``, with a notice); only Layer 1
       applies. If the Landlock argv cannot be built, the command is
       refused, not run unconfined.
     * ``off``: neither Layer 1 nor confinement.

  3. Audit log -- every command, blocked or run, appended to
     ``$XDG_CACHE_HOME/delfin/agent-audit.jsonl`` (created 0600); a command
     that ran is recorded with exit code, output sizes and the applied
     ``fs`` and ``net``.
  4. Approval -- the caller's responsibility: the dashboard shows the full
     command for approval, then prints the applied confinement and its
     notice when it executes it; see ``dashboard/tab_agent.py``.

Env vars (read on every call so tests can flip them with monkeypatch):

  ``DELFIN_AGENT_SANDBOX``           one of {auto,off,allowlist,bwrap,firejail}.
                                     Default: auto, which resolves to the
                                     first that works of bwrap -> firejail
                                     (installed) -> allowlist. ``off``
                                     disables both the allow-list and the
                                     confinement (debug only).
  ``DELFIN_AGENT_SANDBOX_NETWORK``   ``1`` keeps the network inside bwrap,
                                     firejail and the Landlock allowlist
                                     mode. Default: denied there, as far as
                                     the host can hold it (see above).
  ``DELFIN_AGENT_SANDBOX_TIMEOUT``   per-command timeout in seconds. Default 60.

The agent's own shell (``api_client._bash_isolation_argv``) is a different
path with different layers; see docs/SECURITY_MODEL.md.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Layer 1 — allow-list

# First token of each pipeline segment must be in this set OR (for the few
# commands with sub-command grammar) match an entry in ``_ALLOW_SUBCOMMAND``.
# A tool with a sub-command list is never also in this set: _check_segment
# read this set first, so ``pip`` here made ``pip install`` pass and the pip
# list unread.
#
# The limit of this layer: it reads the first word, and an interpreter on
# it runs anything. ``python -c``, a script file, ``awk 'BEGIN{system()}'``,
# ``find -exec`` and ``sed e`` are not contained by it; the OS confinement
# of the runner (bwrap, firejail, Landlock -- see ``confinement``) is what
# limits them, and where there is none the dashboard says so.
_ALLOW_FIRST_TOKEN = frozenset({
    # Read-only inspection
    "ls", "cat", "head", "tail", "less", "more", "wc", "grep", "egrep",
    "fgrep", "rg", "find", "stat", "file", "tree", "du", "df", "echo",
    "pwd", "env", "which", "whoami", "uname", "date", "id", "hostname",
    "basename", "dirname", "realpath", "readlink", "sort", "uniq", "tr",
    "cut", "awk", "sed", "diff", "cmp", "tee", "xargs",
    # Python tooling
    "python", "python3", "pytest", "ruff", "black", "mypy",
    "flake8", "isort",
    # DELFIN entry-points
    "delfin", "delfin-build", "delfin-voila", "delfin-json",
    "delfin_ESD", "delfin_IR",
    # Chemistry / QM tools commonly invoked from DELFIN
    "xtb", "orca", "crest", "obabel", "x2t", "t2x", "molden", "multiwfn",
    "censo", "anmr",
    # Process / shell built-ins that are safe
    "ps", "true", "false", "test", "[", "yes", "sleep", "wait",
    "nohup",
})

# Per-binary subcommand allow-list. Used when the binary itself is
# multi-purpose (e.g. ``git``, ``pip``) — we allow only read-only or harmless
# subcommands here. Other subcommands fall through to deny. ``pip`` also
# covers ``pip3``, ``pip3.N`` and ``<python> -m pip`` (``_check_argv``).
_ALLOW_SUBCOMMAND = {
    "git": frozenset({
        "status", "log", "diff", "show", "branch", "remote", "rev-parse",
        "rev-list", "stash", "ls-files", "ls-tree", "blame", "describe",
        "tag", "shortlog", "fsck", "cat-file", "for-each-ref", "name-rev",
        "config",  # read forms only (_SUBCOMMAND_READ_ONLY)
    }),
    "pip": frozenset({"list", "show", "freeze", "check", "config"}),
}

#: Sub-commands on the lists above whose writing forms rewrite configuration
#: that a later git or pip OUTSIDE the sandbox trusts: ``core.pager``,
#: ``core.hooksPath``, ``core.fsmonitor`` and an alias run code, a remote URL
#: receives the next push, ``global.index-url`` serves the next install.
#: For these only the read forms pass (``_read_only_form``).
_GIT_CONFIG_READ_FLAGS = frozenset({
    "--get", "--get-all", "--get-regexp", "--get-urlmatch", "--get-color",
    "--get-colorbool", "--list", "-l"})
_GIT_CONFIG_WRITE_WORDS = frozenset({
    "--add", "--unset", "--unset-all", "--replace-all", "--rename-section",
    "--remove-section", "--edit", "-e", "set", "unset", "edit",
    "rename-section", "remove-section"})
#: git config options whose value is the next word.
_GIT_CONFIG_VALUE_OPTS = frozenset({
    "-f", "--file", "--blob", "--type", "--default", "--comment"})

# Always-deny substrings (case-insensitive). These win even if the
# first-token check passes — designed to catch obvious destructive or
# escalation patterns and unsafe redirects.
_DENY_SUBSTRINGS = (
    " --no-verify",
    " --force ", " --force\t",
    " push --force", " push -f ",
    " reset --hard",
    " clean -fd", " clean -ffd",
    "rm -rf /", "rm -rf ~", "rm -rf .",
    "rm -fr /", "rm -fr ~", "rm -fr .",
    "mkfs", "dd if=/", "dd of=/dev/sd",
    "> /dev/sd", "> /dev/nvme",
    ":(){ :|:& };:",
    "sudo ", "sudo\t", " su -", " su\t-",
    "chmod 777", "chmod -r 777", "chmod -r 7",
    "git config --global",
    "ssh-keygen", "ssh-add",
    "curl http",  # discourage remote-fetch + pipe-to-shell
    "wget http",
    " | sh", " | bash", " |sh", " |bash",
    "/etc/passwd", "/etc/shadow",
    "~/.ssh", "~/.aws", "~/.gnupg",
    "/.ssh/id_", "/.ssh/known_hosts",
    "/.aws/credentials",
)


# Commands that recurse through a directory tree. Harmless on a job folder,
# but rooted at a whole file system they cost one metadata RPC per file --
# hundreds of thousands of them on a populated HOME -- which degrades a shared
# machine for every user on it. A site's operations team will read that as
# abuse and block the account, which is what prompted this guard.
#
# The rule derives its roots at runtime from the current user's HOME and the
# machine's mount table, so it holds for any user on any cluster without
# knowing site-specific paths.
_TREE_WALK_COMMANDS = frozenset({"du", "ncdu", "find", "tree", "rg"})
#: Walk only when told to recurse: grep -r/-R, ls -R.
_RECURSIVE_ON_FLAG = frozenset({"grep", "egrep", "fgrep", "ls"})


def _is_tree_walk_root(resolved: Path) -> bool:
    """True if a recursive walk rooted here would sweep a whole file system.

    Two universal signals, no hardcoded site paths:
      * the path is the user's HOME or an ancestor of it, and
      * the path is a mount point (``/scratch``, ``/work``, project file
        systems, ``/`` itself) -- a job directory never is.
    """
    if str(resolved) == os.sep:
        return True
    try:
        home = Path.home().resolve()
        if resolved == home or resolved in home.parents:
            return True
    except Exception:
        pass
    try:
        return os.path.ismount(str(resolved))
    except OSError:
        return False


def _recurses(base: str, toks: list[str]) -> bool:
    """Whether this command walks a tree at all."""
    if base in _TREE_WALK_COMMANDS:
        return True
    if base not in _RECURSIVE_ON_FLAG:
        return False
    letter = "R" if base == "ls" else "rR"
    for t in toks[1:]:
        if t in ("--recursive", "--dereference-recursive"):
            return True
        if t.startswith("-") and not t.startswith("--") \
                and any(c in t[1:] for c in letter):
            return True
    return False


def _check_tree_walk(toks: list[str], cwd=None) -> tuple[bool, str]:
    """Deny recursive walks rooted at a whole file system; subdirs are fine.

    The one implementation of this rule. The approval runner of the CLI
    backend reaches it through is_allowed, and the gate of every API and
    terminal session through tree_walk_refusal: it used to live on the
    CLI path only, so the sessions that run most of the work walked past
    it (measured 2026-09-22: `find ~ -name '*.py'`, three times, and walks
    over the archive, all executed).
    """
    base = os.path.basename(toks[0])
    if not _recurses(base, toks):
        return True, "ok"

    targets = [t for t in toks[1:] if not t.startswith("-")]
    if not targets or (base == "rg" and len(targets) == 1) \
            or (base in ("grep", "egrep", "fgrep") and len(targets) == 1):
        # No path given (rg/grep: only the pattern): the walk is rooted
        # at the working directory.
        targets.append(".")

    for target in targets:
        try:
            path = Path(os.path.expanduser(target))
            if not path.is_absolute() and cwd is not None:
                path = Path(cwd) / path       # where the COMMAND runs
            resolved = path.resolve()
        except Exception:
            continue
        if _is_tree_walk_root(resolved):
            return False, (
                f"{base} would walk all of {str(resolved)!r}; on a shared file "
                f"system that overloads the metadata servers. For disk usage "
                f"use delfin.quota.home_usage(); otherwise point {base} at a "
                f"specific subdirectory."
            )
    return True, "ok"


def tree_walk_refusal(cmd: str, cwd=None) -> Optional[str]:
    """Why *cmd* may not run as a walk over a whole file system, or None.

    Checks every pipeline segment, whatever the command's other merits --
    a walk rooted at the home, one of its ancestors or a mount point is
    refused in every mode, before anything is asked. Relative paths are
    read from *cwd*, the directory the command runs in (not this
    process's), and a `cd` inside the command moves it: `cd ~ && find .`
    is a walk over the home.
    """
    try:
        segments = _split_pipeline(cmd.strip())
    except ValueError:
        segments = [s for s in re.split(r"\|\||&&|[|;&\n]", cmd) if s.strip()]
    for seg in segments:
        try:
            toks = shlex.split(seg, posix=True)
        except ValueError:
            toks = seg.split()
        while toks and "=" in toks[0] and toks[0].split("=", 1)[0].isidentifier():
            toks = toks[1:]
        # Wrappers that run the next word: the walk is that word's.
        while toks and os.path.basename(toks[0]) in ("nice", "ionice", "time",
                                                      "timeout", "command",
                                                      "exec"):
            toks = toks[1:]
            # their own options and values: -n 10, -c2, a duration like 5s
            while toks and (toks[0].startswith("-")
                            or re.fullmatch(r"\d+[smhd]?", toks[0])):
                toks = toks[1:]
        if not toks:
            continue
        if toks[0] == "cd":
            dest = toks[1] if len(toks) > 1 else "~"
            try:
                p = Path(os.path.expanduser(dest))
                if not p.is_absolute() and cwd is not None:
                    p = Path(cwd) / p
                cwd = str(p.resolve())
            except Exception:
                pass
            continue
        ok, reason = _check_tree_walk(toks, cwd)
        if not ok:
            return reason
    return None


# ---------------------------------------------------------------------------
# Installs that write outside the session venv

#: pip / uv long options that move an install out of the interpreter's own
#: environment. pip parses with optparse, which accepts any unambiguous
#: prefix (``--targ`` is ``--target``), so a token is matched as a prefix
#: of these names, not only as the full name.
_OUTSIDE_LONG = ("user", "target", "prefix", "root", "system",
                 "break-system-packages", "python")
#: Full option names that are a prefix of one above and mean something else.
_OUTSIDE_LONG_EXACT_SAFE = frozenset({"pre"})
#: Short options whose value follows in the same token (``-rreq.txt``): the
#: rest of a cluster after one of these is a value, not more options.
_SHORT_WITH_VALUE = frozenset("rceif")
#: Environment forms of the same options.
_OUTSIDE_ENV_NAMES = frozenset({"PIP_USER", "PIP_TARGET", "PIP_PREFIX",
                                "PIP_ROOT", "PIP_PYTHON", "UV_PYTHON",
                                "UV_SYSTEM_PYTHON"})
_OUTSIDE_ENV_RE = re.compile(
    r"\b(?:" + "|".join(sorted(_OUTSIDE_ENV_NAMES)) + r")\s*=")
#: A command substitution inside a word (``"$(pip ...)"``, backticks): the
#: tokenizer keeps it inside the quoted token, so it is read separately.
_SUBSTITUTION_RE = re.compile(r"\$\(([^()]*)\)|`([^`]*)`")
_PIP_NAME_RE = re.compile(r"pip(?:3(?:\.\d+)?)?")
_PYTHON_NAME_RE = re.compile(r"python(?:3(?:\.\d+)?)?")
_WRAPPERS = frozenset({"nice", "ionice", "nohup", "time", "timeout", "command",
                       "exec", "stdbuf", "sudo", "doas"})
_SHELLS = frozenset({"sh", "bash", "zsh", "dash", "ksh"})


def _same_env_dir(path: str, session_python: str) -> bool:
    """Whether *path* is an executable in the session interpreter's bin dir.

    The directory is compared resolved, the file name is not: a venv's
    ``python`` is a symlink to the base interpreter, so resolving it would
    make the system Python equal to the session venv.
    """
    if not session_python or "/" not in path:
        return False
    try:
        return (os.path.realpath(os.path.dirname(os.path.abspath(path)))
                == os.path.realpath(os.path.dirname(session_python)))
    except Exception:
        return False


def _in_workspace(path: str, workspace: Optional[str]) -> bool:
    """Whether the executable *path* lies inside *workspace*.

    A venv the session made in its workspace (``.venv/bin/pip``) installs
    into the workspace, which is where the session may write; the
    workspace-local venv tools are auto-allowed for that reason. A
    relative path is read from the workspace.
    """
    if not workspace or "/" not in path:
        return False
    try:
        root = Path(os.path.realpath(workspace))
        p = Path(os.path.expanduser(path))
        if not p.is_absolute():
            p = root / p
        return Path(os.path.realpath(p.parent)).is_relative_to(root)
    except Exception:
        return False


def _is_session_python(path: str, session_python: str) -> bool:
    if not session_python:
        return False
    if path == session_python:
        return True
    return (_PYTHON_NAME_RE.fullmatch(os.path.basename(path)) is not None
            and _same_env_dir(path, session_python))


def _segments(cmd: str) -> list[list[str]]:
    """*cmd* as argv lists, split at shell operators, quotes removed."""
    try:
        lex = shlex.shlex(cmd, posix=True, punctuation_chars="();<>|&`")
        lex.whitespace_split = True
        toks = list(lex)
    except ValueError:
        toks = [t for t in re.split(r"(\s+|[;&|()])", cmd) if t.strip()]
    out: list[list[str]] = [[]]
    for t in toks:
        if t and all(c in ";&|()<>`" for c in t):
            out.append([])
        else:
            out[-1].append(t)
    return [seg for seg in out if seg]


def _strip_prefix(argv: list[str], env_names: set) -> list[str]:
    """argv without leading assignments and wrappers (``env``, ``nice`` ...)."""
    while argv:
        head = argv[0]
        if "=" in head and head.split("=", 1)[0].isidentifier():
            env_names.add(head.split("=", 1)[0])
            argv = argv[1:]
            continue
        base = os.path.basename(head)
        if base == "env" or base in _WRAPPERS:
            argv = argv[1:]
            while argv and (argv[0].startswith("-")
                            or re.fullmatch(r"\d+[smhd]?", argv[0])):
                opt = argv[0]
                argv = argv[1:]
                if base == "env" and opt in ("-u", "--unset", "-C",
                                             "--chdir") and argv:
                    argv = argv[1:]
            continue
        return argv
    return argv


def _pip_args_refusal(args: list[str], session_python: str, *,
                      uv: bool) -> Optional[str]:
    """Why these pip / uv-pip arguments install outside the session venv."""
    if "config" in args and "set" in args:
        at = args.index("set") + 1
        key = args[at] if at < len(args) else ""
        if key.rsplit(".", 1)[-1].replace("_", "-") in _OUTSIDE_LONG:
            return f"pip config set {key} makes every later install leave the venv"
        return None
    if "install" not in args and not (uv and "sync" in args):
        return None
    for i, tok in enumerate(args):
        if tok == "--":
            break
        nxt = args[i + 1] if i + 1 < len(args) else ""
        if tok.startswith("--") and len(tok) > 2:
            name, eq, value = tok[2:].partition("=")
            if not eq:
                value = nxt
            for opt in _OUTSIDE_LONG:
                if name == opt or (len(name) >= 3 and opt.startswith(name)
                                   and name not in _OUTSIDE_LONG_EXACT_SAFE):
                    if opt == "python" and _is_session_python(value, session_python):
                        break
                    return f"--{opt} installs outside the session venv"
        elif tok.startswith("-") and len(tok) > 1:
            for j, ch in enumerate(tok[1:], start=1):
                if ch == "t":
                    return "-t (--target) installs outside the session venv"
                if uv and ch == "p":
                    if _is_session_python(tok[j + 1:] or nxt, session_python):
                        break
                    return "-p (--python) installs into another interpreter"
                if ch in _SHORT_WITH_VALUE:
                    break
    return None


def _python_module(rest: list[str]) -> tuple[str, list[str]]:
    """The ``-m`` module of a python argv and the arguments after it.

    Interpreter options may precede ``-m`` (``-I -m pip``), be clustered
    with it (``-Im pip``) or carry the module in the same token (``-mpip``).
    """
    i = 0
    while i < len(rest):
        tok = rest[i]
        if not tok.startswith("-") or tok == "-" or tok.startswith("--"):
            return "", []
        for j, ch in enumerate(tok[1:], start=1):
            if ch == "m":
                if tok[j + 1:]:
                    return tok[j + 1:], rest[i + 1:]
                if i + 1 < len(rest):
                    return rest[i + 1], rest[i + 2:]
                return "", []
            if ch == "c":
                return "", []
            if ch in "XW":
                if not tok[j + 1:]:
                    i += 1          # the value is the next token
                break
        i += 1
    return "", []


def _segment_refusal(argv: list[str], session_python: str, raw: str,
                     depth: int, workspace: Optional[str]) -> Optional[str]:
    env_names: set = set()
    argv = _strip_prefix(argv, env_names)
    if not argv:
        return None
    head = argv[0]
    base = os.path.basename(head)
    rest = argv[1:]

    if base in _SHELLS and "-c" in rest and depth < 3:
        # The payload of `sh -c` is a command line of its own.
        at = rest.index("-c") + 1
        return unsafe_install_refusal(rest[at] if at < len(rest) else "",
                                      session_python, workspace=workspace,
                                      _depth=depth + 1)
    if base == "eval" and depth < 3:
        return unsafe_install_refusal(" ".join(rest), session_python,
                                      workspace=workspace, _depth=depth + 1)

    pip_args: Optional[list[str]] = None
    uv = False
    if _PIP_NAME_RE.fullmatch(base):
        if "/" in head and "install" in rest \
                and not _same_env_dir(head, session_python) \
                and not _in_workspace(head, workspace):
            return f"{head} is another interpreter's pip, not the session venv's"
        pip_args = rest
    elif _PYTHON_NAME_RE.fullmatch(base):
        module, margs = _python_module(rest)
        if module != "pip":
            return None
        if "install" in margs and not _is_session_python(head, session_python) \
                and not (_PYTHON_NAME_RE.fullmatch(base)
                         and _in_workspace(head, workspace)):
            return f"{head} -m pip installs into {head}, not into the session venv"
        pip_args = margs
    elif base == "uv" and rest[:1] == ["pip"]:
        pip_args, uv = rest[1:], True
    elif base == "uv" and rest[:2] == ["tool", "install"]:
        return "uv tool install puts the package into the home directory"
    elif base == "pipx" and rest[:1] and rest[0] in ("install", "inject"):
        return f"pipx {rest[0]} puts the package into the home directory"
    elif base in ("conda", "mamba", "micromamba") and "install" in rest:
        for i, tok in enumerate(rest):
            if (tok in ("-n", "--name") and rest[i + 1:i + 2] == ["base"]) \
                    or tok in ("-nbase", "--name=base"):
                return f"{base} install into base changes the shared environment"
        return None
    if pip_args is None:
        return None
    if ("install" in pip_args or (uv and "sync" in pip_args)) and (
            env_names & _OUTSIDE_ENV_NAMES or _OUTSIDE_ENV_RE.search(raw)):
        return "a PIP_/UV_ variable moves the install outside the session venv"
    return _pip_args_refusal(pip_args, session_python, uv=uv)


def unsafe_install_refusal(cmd: str, session_python: str, *,
                           workspace: Optional[str] = None,
                           _depth: int = 0) -> Optional[str]:
    """Why *cmd* installs Python packages outside the session venv, or None.

    *session_python* is the interpreter DELFIN installs into
    (:func:`delfin.installer.session_python`); '' means it is unknown, and
    then every ``python -m pip install`` and every pip named by path is
    refused. Read per shell segment after quote removal, with leading
    assignments and wrappers (``env``, ``nice``, ``sudo`` ...) dropped and
    ``sh -c`` / ``eval`` payloads read as commands. Refused:

      * pip / ``uv pip`` install or sync with ``--user``, ``--target``,
        ``--prefix``, ``--root``, ``--system``, ``--break-system-packages``
        or ``--python`` naming another interpreter, each also as an
        unambiguous prefix (optparse resolves ``--targ``), and ``-t`` /
        uv ``-p`` also inside a short-option cluster (``-qt``);
      * the same through ``PIP_USER``/``PIP_TARGET``/``PIP_PREFIX``/
        ``PIP_ROOT``/``PIP_PYTHON``/``UV_PYTHON``/``UV_SYSTEM_PYTHON``,
        and ``pip config set <scope>.user|target|prefix|root|python``;
      * ``<python> -m pip install`` (interpreter options before ``-m``
        and ``-mpip`` included) unless *python* is the session interpreter
        or a ``python*`` in its bin directory, and ``/path/to/pip install``
        from another bin directory. Both are allowed for an interpreter
        inside *workspace* (a venv the session made there writes only the
        workspace). A bare ``pip`` is the one on PATH and is judged by its
        options only;
      * ``pipx install|inject`` and ``uv tool install`` (both write the
        home directory) and conda/mamba/micromamba install into ``base``.

    Not covered: a Python program that calls pip itself (``python -c``, a
    script file). The command text does not show that act; the write roots
    of the OS sandbox are what stop it.
    """
    if not cmd or not cmd.strip():
        return None
    if _depth < 3:
        for m in _SUBSTITUTION_RE.finditer(cmd):
            reason = unsafe_install_refusal(m.group(1) or m.group(2) or "",
                                            session_python, workspace=workspace,
                                            _depth=_depth + 1)
            if reason:
                return reason
    for seg in _segments(cmd):
        reason = _segment_refusal(seg, session_python, cmd, _depth, workspace)
        if reason:
            return reason
    return None


def _session_python() -> str:
    """:func:`delfin.installer.session_python`, or '' when it cannot load."""
    try:
        from delfin import installer as _installer
        return _installer.session_python()
    except Exception:
        return ""


@dataclass
class AllowResult:
    allowed: bool
    reason: str


def is_allowed(cmd: str) -> AllowResult:
    """Layer 1 — allow-list check on a shell command string.

    Splits the command into pipeline segments at ``|``, ``||``, ``&&``,
    ``;``, ``&`` and verifies each segment's first token. Also rejects on
    any always-deny substring.
    """
    cmd_stripped = cmd.strip()
    if not cmd_stripped:
        return AllowResult(False, "empty command")

    # Always-deny substring scan (case-insensitive, with surrounding spaces
    # so word-boundary patterns like " --force " don't false-positive on a
    # filename like ``--force-tag.txt``).
    haystack = " " + cmd_stripped.lower() + " "
    for needle in _DENY_SUBSTRINGS:
        if needle in haystack:
            return AllowResult(False, f"matched deny pattern: {needle.strip()!r}")

    # The same install rule as the gate of the API and terminal sessions.
    install = unsafe_install_refusal(cmd_stripped, _session_python())
    if install:
        return AllowResult(False, f"unsafe install target: {install}")

    # Pipeline-aware first-token check
    try:
        segments = _split_pipeline(cmd_stripped)
    except ValueError as exc:
        return AllowResult(False, f"unparseable command: {exc}")
    for seg in segments:
        ok, reason = _check_segment(seg)
        if not ok:
            return AllowResult(False, reason)
    return AllowResult(True, "ok")


def _split_pipeline(cmd: str) -> list[str]:
    """Split *cmd* at unquoted ``|``, ``||``, ``&&``, ``;``, ``&``."""
    toks = shlex.split(cmd, posix=True)
    out: list[list[str]] = [[]]
    for t in toks:
        if t in ("|", "||", "&&", ";", "&"):
            out.append([])
        else:
            out[-1].append(t)
    return [" ".join(seg) for seg in out if seg]


def _check_segment(seg: str) -> tuple[bool, str]:
    toks = shlex.split(seg, posix=True)
    if not toks:
        return False, "empty pipeline segment"
    return _check_argv(toks)


#: Allowed commands that run the command given as their arguments.
_RUNS_ITS_ARGUMENTS = frozenset({"env", "nohup", "xargs"})
#: xargs options whose value is the next word when the option ends the
#: word (GNU: -a -d -E -I -L -n -P -s; BSD adds -J -R -S).
_XARGS_VALUE_SHORT = frozenset("adEILnPsJRS")
_XARGS_VALUE_LONG = frozenset({"--arg-file", "--delimiter", "--max-args",
                               "--max-procs", "--max-chars",
                               "--process-slot-var"})


def _wrapped_command(base: str, args: list[str]) -> tuple[list[str], str]:
    """The command *base* (env, nohup, xargs) runs, and a refusal reason.

    Returns ``(argv, "")``; ``argv`` is empty when the wrapper runs no
    command of its own (``env`` alone prints the environment) and
    ``["echo"]`` for a bare ``xargs``. A word this parser cannot place is
    taken as the command, so a misreading refuses rather than skips."""
    i = 0
    if base == "nohup":
        if args[:1] == ["--"]:
            i = 1
        return args[i:], ""
    if base == "env":
        while i < len(args):
            a = args[i]
            if a.startswith("-S") or a.startswith("--split-string") \
                    or (a.startswith("-") and not a.startswith("--")
                        and "S" in a[1:]):
                return [], ("env -S splits a string into a command line; "
                            "write the command out")
            if a == "--":
                i += 1
                break
            if a in ("-u", "--unset", "-C", "--chdir"):
                i += 2
                continue
            if a.startswith("-") or ("=" in a and not a.startswith("=")):
                i += 1
                continue
            break
        return args[i:], ""
    # xargs
    while i < len(args):
        a = args[i]
        if a == "--":
            i += 1
            break
        if a.startswith("--"):
            i += 2 if a in _XARGS_VALUE_LONG else 1
            continue
        if a.startswith("-") and len(a) > 1:
            for j, ch in enumerate(a[1:], start=1):
                if ch in _XARGS_VALUE_SHORT:
                    if j == len(a) - 1:
                        i += 1          # the value is the next word
                    break
                if ch in "eil":         # optional value, attached only
                    break
            i += 1
            continue
        break
    return (args[i:] or ["echo"]), ""


def _read_only_form(tool: str, sub: str, args: list[str]) -> bool:
    """Whether ``<tool> <sub> <args>`` only reads (see
    ``_GIT_CONFIG_READ_FLAGS``). Sub-commands without a writing form
    are read-only by being on the list."""
    if (tool, sub) == ("pip", "config"):
        return bool(args) and args[0] in ("list", "get", "debug")
    if (tool, sub) == ("git", "remote"):
        rest = [a for a in args if a not in ("-v", "--verbose")]
        return not rest or rest[0] in ("show", "get-url")
    if (tool, sub) == ("git", "config"):
        if any(a.split("=", 1)[0] in _GIT_CONFIG_WRITE_WORDS for a in args):
            return False
        if args and args[0] in ("get", "list"):
            return True
        if any(a.split("=", 1)[0] in _GIT_CONFIG_READ_FLAGS for a in args):
            return True
        positional, skip = [], False
        for a in args:
            if skip:
                skip = False
            elif a in _GIT_CONFIG_VALUE_OPTS:
                skip = True
            elif not a.startswith("-"):
                positional.append(a)
        return len(positional) == 1     # ``git config <name>`` reads it
    return True


def _check_argv(toks: list[str], depth: int = 0) -> tuple[bool, str]:
    """Layer 1 for one segment's argv.

    Order: a tool with a sub-command list is judged by that list (and by
    ``_read_only_form``), whatever else names it -- ``pip``, ``pip3``,
    ``pip3.N``, ``<python> -m pip`` and ``-m pip.<sub>``; a wrapper
    (``_RUNS_ITS_ARGUMENTS``) is judged by the command it runs; then the
    first-token set and the tree-walk rule."""
    while toks and "=" in toks[0] and toks[0].split("=", 1)[0].isidentifier():
        toks = toks[1:]
    if not toks:
        return False, "command was just env-vars"
    base = os.path.basename(toks[0])
    tool, sub_argv = base, toks[1:]
    if _PIP_NAME_RE.fullmatch(base):
        tool = "pip"
    elif _PYTHON_NAME_RE.fullmatch(base):
        module, margs = _python_module(toks[1:])
        if module == "pip" or module.startswith("pip."):
            tool, sub_argv = "pip", margs
    if tool in _ALLOW_SUBCOMMAND:
        sub = sub_argv[0] if sub_argv else ""
        if sub not in _ALLOW_SUBCOMMAND[tool]:
            return False, f"{tool} subcommand not in allow-list: {sub!r}"
        if not _read_only_form(tool, sub, sub_argv[1:]):
            return False, (f"{tool} {sub}: only its read forms are allowed "
                           "(the written setting is used outside the sandbox)")
        return True, "ok"
    if base in _RUNS_ITS_ARGUMENTS:
        inner, reason = _wrapped_command(base, toks[1:])
        if reason:
            return False, reason
        if inner:
            if depth >= 3:
                return False, f"{base} nested too deep to read"
            return _check_argv(inner, depth + 1)
        return True, "ok"
    if base in _ALLOW_FIRST_TOKEN:
        return _check_tree_walk(toks)
    return False, f"first token not in allow-list: {base!r}"


# ---------------------------------------------------------------------------
# Layer 2 — sandbox

@dataclass
class SandboxConfig:
    mode: str           # "off" | "allowlist" | "bwrap" | "firejail"
    allow_network: bool
    timeout_s: int


def detect_config(env: Optional[dict] = None) -> SandboxConfig:
    """Read env vars and decide which sandbox mode to use."""
    e = env if env is not None else os.environ
    requested = (e.get("DELFIN_AGENT_SANDBOX") or "auto").strip().lower()
    allow_net = e.get("DELFIN_AGENT_SANDBOX_NETWORK", "0").strip() == "1"
    try:
        timeout_s = int(e.get("DELFIN_AGENT_SANDBOX_TIMEOUT", "60"))
    except ValueError:
        timeout_s = 60

    if requested == "off":
        return SandboxConfig("off", allow_net, timeout_s)
    if requested == "allowlist":
        return SandboxConfig("allowlist", allow_net, timeout_s)
    if requested == "bwrap":
        return SandboxConfig("bwrap", allow_net, timeout_s)
    if requested == "firejail":
        return SandboxConfig("firejail", allow_net, timeout_s)

    # auto (and any unknown value) — pick best available. Installed is not
    # working: a bwrap whose user namespace is refused made every approved
    # command fail.
    if shutil.which("bwrap") and _bwrap_works():
        return SandboxConfig("bwrap", allow_net, timeout_s)
    if shutil.which("firejail"):
        return SandboxConfig("firejail", allow_net, timeout_s)
    return SandboxConfig("allowlist", allow_net, timeout_s)


def _bwrap_works() -> bool:
    try:
        from .api_client import _bwrap_functional
        return bool(_bwrap_functional())
    except Exception:
        return False


@dataclass(frozen=True)
class Confinement:
    """What the runner applies to a command in one mode on this host.

    ``filesystem``: "bwrap", "firejail", "landlock" or "none".
    ``network``: "off" (no connect leaves the command), "off-tcp" (TCP
    connects refused by Landlock ABI 4, UDP not restricted), "on" (opened
    by DELFIN_AGENT_SANDBOX_NETWORK=1) or "open" (not restricted although
    the configuration asked for it, or a mode with no sandbox).
    ``landlock_net``: the ``--net`` value handed to ``landlock_exec``; ""
    when the command does not run under the Landlock helper.
    ``notice``: one sentence when something is not in force; "" otherwise.
    """
    mode: str
    filesystem: str
    network: str
    landlock_net: str
    notice: str

    @property
    def network_restricted(self) -> bool:
        return self.network in ("off", "off-tcp")

    def label(self) -> str:
        return f"mode={self.mode}, fs={self.filesystem}, net={self.network}"


def confinement(cfg: SandboxConfig) -> Confinement:
    """The confinement :func:`run_agent_command` applies for *cfg* here.

    The one answer to "what holds this command": the runner builds its
    argv from it (``_unsandboxed_argv``) and the dashboard prints it, so a
    state is shown only where the runner applies it. Probes, in the order
    of ``landlock_exec._run_guarded``: Landlock usable
    (``api_client._landlock_functional``), seccomp user notification
    (``socket_guard.available``: every connect judged, UDP included),
    Landlock ABI >= 4 with the seccomp filter (TCP connects only).

    Rejected: refusing every command where the network cannot be held.
    The command was approved by the person at the screen, and the hosts
    without bubblewrap are the ones where this mode is the only one; it
    runs, and the state and the notice say what is missing.
    Never raises; a failed probe counts as "not available".
    """
    mode = cfg.mode
    if mode in ("bwrap", "firejail"):
        return Confinement(mode, mode, "on" if cfg.allow_network else "off",
                           "", "")
    if mode != "allowlist":
        return Confinement(mode, "none", "open", "", (
            "sandbox off: no allow-list, no filesystem or network "
            "confinement"))
    try:
        from .api_client import _landlock_functional
        landlock = bool(_landlock_functional())
    except Exception:
        landlock = False
    if not landlock:
        return Confinement(mode, "none", "open", "", (
            "no bubblewrap, firejail or Landlock on this host: the command "
            "runs unconfined (allow-list and scrubbed environment only); "
            "the filesystem and the network are NOT restricted"))
    if cfg.allow_network:
        return Confinement(mode, "landlock", "on", "open", "")
    try:
        from . import socket_guard as _sg
        guard = bool(_sg.available())
        filt = bool(_sg.filter_supported())
    except Exception:
        guard = filt = False
    try:
        from . import landlock_exec as _ll
        abi = int(_ll.abi_version())
    except Exception:
        abi = 0
    if guard:
        return Confinement(mode, "landlock", "off", "none", "")
    if abi >= 4 and filt:
        return Confinement(mode, "landlock", "off-tcp", "none", (
            "network partly restricted: TCP connects are refused (Landlock), "
            "UDP is NOT restricted on this host (no seccomp user "
            "notification)"))
    return Confinement(mode, "landlock", "open", "open", (
        "network NOT restricted on this host: no seccomp user notification "
        "and no Landlock TCP rule (Linux 6.7+); the filesystem is confined "
        "(Landlock)"))


def _unsandboxed_argv(cmd: str, repo_dir: Path, mode: str, *,
                      allow_network: bool = False) -> list[str]:
    """``bash -c cmd`` for the modes without bubblewrap or firejail.

    The allowlist mode ran commands with nothing around them: python, pip,
    awk and env are on the list, and every key was in the environment. Where
    the kernel has Landlock the command now runs under it (writes only in the
    repository and a private temp directory, credential folders unreadable),
    as the agent's own shell does on such a host. "off" stays off.

    Filesystem and network are the ones :func:`confinement` names, so the
    state the dashboard prints is the state applied. A restricted network
    is passed with ``--strict 1``: should the helper find a layer missing
    that the probe here reported, it refuses the command (exit 126)
    instead of running it with the network open; unrestricted, ``--strict
    0`` as before. Raises when the Landlock argv cannot be built: running
    ``bash -c`` then would be a confinement the dashboard did not show."""
    base = ["bash", "-c", cmd]
    conf = confinement(SandboxConfig(mode, allow_network, 0))
    if conf.filesystem != "landlock":
        return base
    from .api_client import (_landlock_hide_paths, _private_tmp_dir,
                             _session_scratch_dir)
    from . import landlock_exec as _ll
    repo = str(Path(repo_dir).resolve())
    argv = [sys.executable, "-I", str(Path(_ll.__file__).resolve()),
            "--write", repo]
    # The session scratch dir: writable, and the tmpdir the helper maps.
    try:
        scratch = _session_scratch_dir()
        argv += ["--write", scratch]
        tmpdir = scratch
    except Exception:
        tmpdir = _private_tmp_dir()
    for h in _landlock_hide_paths([repo, _private_tmp_dir()]):
        argv += ["--hide", h]
    strict = "1" if conf.network_restricted else "0"
    return argv + ["--tmpdir", tmpdir,
                   "--net", conf.landlock_net,
                   "--strict", strict, "--"] + base


# Home-relative directories that are tmpfs'd inside the sandbox so the agent
# cannot read or exfiltrate them. Everything else under $HOME stays bound
# read-only so language envs (micromamba, conda, venvs) remain usable for
# legitimate commands like ``python -m pytest`` or ``delfin --recalc``.
_HOME_SECRET_DIRS = (
    ".ssh", ".aws", ".gnupg", ".kube", ".docker", ".netrc", ".pgpass",
    ".azure", ".gcp", ".config/gcloud", ".config/gh",
    # NOT .config/git: it is git's own config and ignore file, no
    # credential, and hiding it made git fatal ("cannot use
    # ~/.config/git/ignore as an exclude file"). Credentials live in
    # .git-credentials and a helper, both denied separately.
    ".git-credentials", ".npmrc", ".pypirc",
    ".anthropic", ".openai", ".claude",
    # The framework's own credential store. The other providers' folders
    # were masked and this one was not, which is the wrong way round: it
    # is the file holding the key the running agent is using.
    ".delfin/credentials.json",
    # One-time approval codes of runs with no terminal (approval_code):
    # the person's `delfin-agent approvals` reads them, no agent may.
    ".delfin/approval_codes",
    # The account key and the MCP pins (trust_store): read by DELFIN's
    # own processes, never by a command an agent runs.
    ".delfin/trust_store",
    ".delfin/messenger_sessions",
    # The user's schedules (delfin.agent.schedules): what runs unattended
    # is the user's to write and confirm, never an agent's.
    ".delfin/user_schedules",
)


#: Paths a run hides from its own cage: absolute paths joined by
#: ``os.pathsep``, set by whoever starts the session (a swarm launcher).
RUN_HIDE_ENV = "DELFIN_CAGE_HIDE"


def run_hidden_paths() -> list[str]:
    """The existing absolute paths named by ``DELFIN_CAGE_HIDE``, resolved.

    Input: the environment of the session that builds the cage. Output:
    each entry that is absolute and exists, resolved, once. Semantics:
    hiding only removes access -- every cage layer masks these AFTER its
    binds like the home secret directories, and the file tools refuse them
    inside a readable root. Relative, empty and missing entries are
    ignored (a mount point must exist). Read per cage build, from the
    session's own environment, so a command inside the cage cannot change
    what its cage hides. Never raises."""
    out: list[str] = []
    for raw in str(os.environ.get(RUN_HIDE_ENV) or "").split(os.pathsep):
        raw = raw.strip()
        if not raw or not os.path.isabs(raw):
            continue
        try:
            p = Path(raw).resolve()
            if p.exists() and str(p) not in out:
                out.append(str(p))
        except OSError:
            continue
    return out


def _bwrap_hide_args(paths) -> list[str]:
    """bubblewrap masks for ``paths``: tmpfs over a directory, /dev/null
    over a file -- the same two forms the home secret paths use."""
    args: list[str] = []
    for p in paths:
        try:
            if os.path.isdir(p):
                args += ["--tmpfs", p]
            elif os.path.isfile(p):
                args += ["--ro-bind", "/dev/null", p]
        except OSError:
            continue
    return args


def runtime_masked_dirs() -> list[str]:
    """The masked directories that appear at run time, created now and
    resolved: a cage built now masks them instead of skipping a path that
    does not exist yet, also where the location is redirected (tests, a
    benchmark attempt). A cage masks these AFTER its binds, so a granted
    parent cannot re-expose them. Never raises."""
    out: list[str] = []
    try:
        from . import approval_code as _ac
        d = _ac.ensure_codes_dir()
        if d is not None:
            out.append(str(d.resolve()))
    except Exception:
        pass
    try:
        from . import trust_store as _ts
        d = _ts.ensure_store_dir()
        if d is not None and str(d.resolve()) not in out:
            out.append(str(d.resolve()))
    except Exception:
        pass
    # Not the messenger's activation records: they live in
    # ~/.delfin/messenger_sessions, masked with the other home secrets
    # (_HOME_SECRET_DIRS) and refused by the path deny globs, like the trust
    # store. As a runtime mask the directory was created on every cage build,
    # and a workspace containing it (a workspace at $HOME) lost its write
    # root under Landlock, which cannot carve a hole into an allowed tree.
    return out


def _bwrap_argv(cmd: str, repo_dir: Path, allow_network: bool) -> list[str]:
    """Build a bwrap argv that runs *cmd* in a read-only namespace.

    The repo dir is bound read-write; common credential directories under
    ``$HOME`` are hidden behind tmpfs mounts; ``/root`` and ``/tmp`` are
    fresh tmpfs; ``--unshare-net`` blocks outbound traffic by default.
    """
    repo = str(repo_dir.resolve())
    # Resolve $HOME fully: on HPC clusters /home is often a symlink
    # (e.g. /home -> /pfs) and bwrap cannot create mountpoints through
    # symlinked path components ("Can't mkdir parents").
    home = str(Path.home().resolve())
    args: list[str] = [
        "bwrap",
        "--ro-bind", "/", "/",
        "--dev", "/dev",
        "--proc", "/proc",
        "--tmpfs", "/tmp",
        "--tmpfs", "/root",
    ]
    # Mask each secret path individually if it exists on the host. tmpfs
    # only works on directories; for single files we mask with /dev/null
    # (read-only bind of an empty source). Missing paths are skipped.
    late_masks = runtime_masked_dirs()
    for rel in _HOME_SECRET_DIRS:
        p = Path(home) / rel
        try:
            is_dir = p.is_dir()
            is_file = p.is_file()
        except OSError:
            continue
        if is_dir:
            args += ["--tmpfs", str(p)]
        elif is_file:
            args += ["--ro-bind", "/dev/null", str(p)]
    args += ["--bind", repo, repo]
    # The session scratch dir: bound read-write, TMPDIR/TMP/TEMP point at
    # it (in-cage path == host path, "/" is bound at itself), so a file
    # one command writes survives to the next and the file tools reach it.
    try:
        from .api_client import _session_scratch_dir
        scratch = _session_scratch_dir()
        args += ["--bind-try", scratch, scratch]
        args += ["--setenv", "TMPDIR", scratch,
                 "--setenv", "TMP", scratch,
                 "--setenv", "TEMP", scratch]
    except Exception:
        pass
    for d in late_masks:
        args += ["--tmpfs", d]
    # The run's own hidden paths, after every bind like the masks above.
    args += _bwrap_hide_args(run_hidden_paths())
    args += [
        "--chdir", repo,
        "--unshare-pid",
        "--unshare-ipc",
        "--unshare-uts",
        "--die-with-parent",
        "--new-session",
    ]
    if not allow_network:
        args.append("--unshare-net")
    args += ["--", "bash", "-c", cmd]
    return args


#: Characters firejail expands in a path option: file globs and the
#: ``${HOME}``-style macros.
_FIREJAIL_PATTERN_CHARS = frozenset("*?[]$")


def _firejail_argv(cmd: str, repo_dir: Path, allow_network: bool) -> list[str]:
    """Build a firejail argv with comparable restrictions.

    Raises ``ValueError`` for a repository path holding a glob or macro
    character (``_FIREJAIL_PATTERN_CHARS``) or a control character:
    firejail expands those in ``--whitelist=`` and ``--chdir=``, so the
    whitelist would name every path the pattern matches."""
    repo = str(repo_dir.resolve())
    for ch in repo:
        if ch in _FIREJAIL_PATTERN_CHARS or ord(ch) < 0x20 or ord(ch) == 0x7f:
            raise ValueError(
                f"firejail cannot whitelist the repository {repo!r}: it "
                f"contains {ch!r}, which firejail expands as a pattern. Use a "
                "repository path without * ? [ ] $ or control characters")
    args: list[str] = [
        "firejail",
        "--quiet",
        "--noprofile",
        "--private-tmp",
        "--private-dev",
        f"--whitelist={repo}",
        "--read-only=/etc",
        "--read-only=/usr",
        "--noroot",
        "--caps.drop=all",
    ]
    if not allow_network:
        args.append("--net=none")
    for p in run_hidden_paths():
        if any(ch in _FIREJAIL_PATTERN_CHARS or ord(ch) < 0x20 for ch in p):
            raise ValueError(
                f"firejail cannot blacklist {p!r}: it contains a pattern or "
                "control character")
        args.append(f"--blacklist={p}")
    args += [f"--chdir={repo}", "bash", "-c", cmd]
    return args


def _command_env(mode: str) -> dict:
    """The scrubbed environment for a command run under sandbox ``mode``.

    ``_scrubbed_bash_env`` points TMPDIR/TMP/TEMP at the session scratch dir,
    which is correct for the layers that bind that dir read-write into their
    cage (bwrap, Landlock, Seatbelt, and no lock): a file one command writes
    to $TMPDIR survives to the next and the file tools reach it. firejail is
    the exception: it runs ``--private-tmp`` and does not bind the scratch
    dir, so inside its cage a TMPDIR naming the host scratch would point at a
    path that does not exist and mktemp plus every TMPDIR/TMP/TEMP user would
    fail. For firejail keep the variables at /tmp (firejail's private
    writable tmp), restoring the behaviour that predates the scratch dir.
    """
    from .api_client import _scrubbed_bash_env
    env = _scrubbed_bash_env()
    if mode == "firejail":
        for key in ("TMPDIR", "TMP", "TEMP"):
            env[key] = "/tmp"
    return env


# ---------------------------------------------------------------------------
# Layer 3 — audit log

def _audit_log_path() -> Path:
    cache_root = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(cache_root) / "delfin" / "agent-audit.jsonl"


def _audit(record: dict) -> None:
    """Append one JSON record to the per-user audit log (mode 0600)."""
    path = _audit_log_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except OSError:
        return  # audit must never break the agent
    line = json.dumps(record, ensure_ascii=False, default=str) + "\n"
    try:
        # O_APPEND is atomic for small writes; create with 0600.
        fd = os.open(
            str(path),
            os.O_WRONLY | os.O_CREAT | os.O_APPEND,
            0o600,
        )
        try:
            os.write(fd, line.encode("utf-8"))
        finally:
            os.close(fd)
    except OSError:
        return


# ---------------------------------------------------------------------------
# Public entry point

@dataclass
class RunResult:
    blocked: bool
    block_reason: Optional[str]
    returncode: int
    stdout: str
    stderr: str
    mode: str               # the sandbox mode actually used
    timed_out: bool = False
    elapsed_s: float = 0.0
    extra: dict = field(default_factory=dict)


def run_agent_command(
    cmd: str,
    repo_dir: Path,
    *,
    config: Optional[SandboxConfig] = None,
) -> RunResult:
    """Run *cmd* through the layered defense and return its result.

    The caller is responsible for Layer 4 (interactive approval) — this
    function trusts that the user has already approved *cmd* in the UI.
    """
    cfg = config or detect_config()
    t0 = time.monotonic()

    # The shared checkout stays on main, in every mode ('off' included):
    # the same rule the permission gate of the API and terminal sessions
    # applies (shared_checkout.refusal).
    try:
        from .shared_checkout import refusal as _shared_refusal
        _shared = _shared_refusal(cmd, repo_dir)
    except Exception:
        _shared = ""
    if _shared:
        _audit({"ts": time.time(), "cwd": str(repo_dir), "mode": cfg.mode,
                "cmd": cmd, "blocked": True, "block_reason": _shared,
                "elapsed_s": round(time.monotonic() - t0, 4)})
        return RunResult(blocked=True, block_reason=_shared, returncode=-1,
                         stdout="", stderr="", mode=cfg.mode,
                         elapsed_s=time.monotonic() - t0)

    # Layer 1 — allow-list (skipped only in 'off' mode)
    if cfg.mode != "off":
        decision = is_allowed(cmd)
        if not decision.allowed:
            elapsed = time.monotonic() - t0
            _audit({
                "ts": time.time(),
                "cwd": str(repo_dir),
                "mode": cfg.mode,
                "cmd": cmd,
                "blocked": True,
                "block_reason": decision.reason,
                "elapsed_s": round(elapsed, 4),
            })
            return RunResult(
                blocked=True,
                block_reason=decision.reason,
                returncode=-1,
                stdout="",
                stderr="",
                mode=cfg.mode,
                elapsed_s=elapsed,
            )

    # Layer 2 — choose argv based on sandbox mode
    conf = confinement(cfg)
    if cfg.mode == "bwrap":
        argv = _bwrap_argv(cmd, repo_dir, cfg.allow_network)
        shell = False
    elif cfg.mode == "firejail":
        try:
            argv = _firejail_argv(cmd, repo_dir, cfg.allow_network)
        except ValueError as exc:
            elapsed = time.monotonic() - t0
            _audit({
                "ts": time.time(),
                "cwd": str(repo_dir),
                "mode": cfg.mode,
                "cmd": cmd,
                "blocked": True,
                "block_reason": str(exc),
                "elapsed_s": round(elapsed, 4),
            })
            return RunResult(
                blocked=True,
                block_reason=str(exc),
                returncode=-1,
                stdout="",
                stderr="",
                mode=cfg.mode,
                elapsed_s=elapsed,
            )
        shell = False
    else:
        # 'off' or 'allowlist' or 'auto-fell-back-to-allowlist'. Routed
        # through bash -c explicitly so the argv stays inspectable; the
        # allowlist mode runs under Landlock where the kernel has it, with
        # the network confinement() names.
        try:
            argv = _unsandboxed_argv(cmd, repo_dir, cfg.mode,
                                     allow_network=cfg.allow_network)
        except Exception as exc:
            reason = (f"the Landlock confinement could not be built ({exc}); "
                      "the command was not run")
            _audit({"ts": time.time(), "cwd": str(repo_dir), "mode": cfg.mode,
                    "cmd": cmd, "blocked": True, "block_reason": reason,
                    "elapsed_s": round(time.monotonic() - t0, 4)})
            return RunResult(blocked=True, block_reason=reason, returncode=-1,
                             stdout="", stderr="", mode=cfg.mode,
                             elapsed_s=time.monotonic() - t0)
        shell = False

    try:
        # Every mode: the shell's scrubbed environment (no provider keys),
        # its own session and process group, no terminal, the whole group
        # ended at exit and at the timeout. TMPDIR points at the session
        # scratch dir except under firejail (which does not bind it): see
        # ``_command_env`` (which applies _scrubbed_bash_env itself).
        from . import contained_run as _contained
        proc = _contained.run(
            argv,
            cwd=str(repo_dir),
            env=_command_env(cfg.mode),
            timeout=cfg.timeout_s,
            budget=True,
        )
        elapsed = time.monotonic() - t0
        _audit({
            "ts": time.time(),
            "cwd": str(repo_dir),
            "mode": cfg.mode,
            "cmd": cmd,
            "blocked": False,
            "exit": proc.returncode,
            "fs": conf.filesystem,
            "net": conf.network,
            "stdout_bytes": len(proc.stdout),
            "stderr_bytes": len(proc.stderr),
            "elapsed_s": round(elapsed, 4),
        })
        return RunResult(
            blocked=False,
            block_reason=None,
            returncode=proc.returncode,
            stdout=proc.stdout,
            stderr=proc.stderr,
            mode=cfg.mode,
            elapsed_s=elapsed,
            extra={"confinement": conf.label()},
        )
    except subprocess.TimeoutExpired as exc:
        elapsed = time.monotonic() - t0
        _audit({
            "ts": time.time(),
            "cwd": str(repo_dir),
            "mode": cfg.mode,
            "cmd": cmd,
            "blocked": False,
            "exit": -1,
            "timed_out": True,
            "elapsed_s": round(elapsed, 4),
        })
        return RunResult(
            blocked=False,
            block_reason=None,
            returncode=-1,
            stdout=(exc.stdout or b"").decode("utf-8", errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or ""),
            stderr=(exc.stderr or b"").decode("utf-8", errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or ""),
            mode=cfg.mode,
            timed_out=True,
            elapsed_s=elapsed,
        )
