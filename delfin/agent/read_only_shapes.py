"""Recognize read-only command shapes approved by rule.

Many approval dialogs of the last waves were read-only commands the bash
auto-allow list (api_client.py:_DEFAULT_BASH_AUTO_ALLOW) does not know:
``git show <ref>:<path>``, ``git grep``, ``git apply --check``,
``git hash-object <rel>``, ``gh --version``,
and ``git switch -c pr/<name> origin/main``.

This module classifies such commands as *shapes*. ``matches`` returns the
shape's name when the command is exactly one of the recognised read-only
forms and nothing more; it returns ``""`` for everything else.

Design rules (reviewer-hardened, wave 15):

* Strictly NARROWER than the git backstop, never floor-equal to it. The
  backstop ``_git_reads_only`` (api_client.py:1071) approves ``--textconv``,
  bare ``-o`` and ``--exec`` on read-only subcommands; this module must
  positively reject them by name.
* A command is parsed with :mod:`shlex` into ONE simple command. Compounds
  (``;`` ``&&`` ``||`` ``|``), substitutions (``$( ... )``, backticks, ``$``),
  redirections (``<`` ``>``) and newlines are rejected outright: a shape is
  a single argv that runs exactly one program.
* Every path argument must resolve inside the workspace.
* A per-shape switch can disable a shape; ``matches_with`` honours that. A
  shape ships green only when every escape variant is tested both ways.
"""

from __future__ import annotations

import shlex
from pathlib import Path
from typing import Iterable, Set

#: Canonical option names that make a reading program write a file or run a
#: program. Every shape rejects any of these anywhere in its argv (some
#: shapes additionally exempt a specifically-safe option, e.g. ``-c`` for
#: ``git switch -c``, which creates a local branch and execs nothing).
#:
#: Derived from the task and from reviewer findings (wave 15): the existing
#: gate regex only catches the ``-O`` spelling for the multi-program escapes
#: and misses ``--textconv`` (runs a configured external filter) and bare
#: ``-o`` (e.g. ``git format-patch -o DIR`` writes). ``_git_reads_only``
#: approves all of them on read-only subcommands, so the backstop is NOT a
#: safe floor — this set is.
_WRITE_OR_EXEC_OPTS: frozenset[str] = frozenset({
    # -- git / file readers that write or run ------------------------------
    "--output",           # write the diff/hunk to a file
    "-o",                 # short form of --output (path-carrying)
    "--ext-diff",         # run an external diff program
    "--textconv",         # run the configured textconv filter (executes)
    "--exec",             # run a program as the outer command
    "--exec-path",        # re-point where git looks for sub-programs
    "--upload-pack",      # specify a program to run on the remote side
    "--open-files-in-pager",  # run a program on each match
    "-c",                 # git global config: git -c core.pager=x <cmd>
    "--config-env",       # config value from an env var
    # -- git hash-object ----------------------------------------------------
    "--write",            # git hash-object -w writes an object into the repo
    # -- git apply / apply --check mutators ---------------------------------
    "--index",            # apply updates the index (writes)
    "--cached",           # apply to the index only (writes)
    "--3way",             # 3-way merge during apply (writes)
    "--reject",           # write .rej files on conflict
    # -- git grep pattern-file read (reviewer-s28) ---------------------------
    # `git grep -f FILE` / `--file=FILE` reads the pattern list from FILE,
    # which may be ANY path on the host. A shape guarantees every path inside
    # the workspace, so the pattern-file option is a read of an unvalidated
    # path and is vetoed. No shape uses -f/--file as a legitimate option:
    # git show / apply --check / hash-object / switch -c / gh --version /
    # gate probe take no -f. The ATTACHED short form `-fFILE` is not caught by
    # _token_is_output_or_run (it only bundles -o and -O); _git_grep rejects
    # that spelling explicitly.
    "-f",                 # git grep -f FILE (read pattern file)
    "--file",             # git grep --file=FILE (long form)
})


def _token_is_output_or_run(tok: str) -> bool:
    """True when *tok* names a write/execute option, in any spelling.

    Matches a long option exactly (``--output``) or with a ``=`` value
    (``--output=/tmp/x``), the short ``-o``, the pager short bundle
    ``-O``/``-O<prog>`` (capital O, git grep), ``-c<name>=<value>`` the git
    config one-token form, ``--3way``, and the apply mutators.
    """
    if tok in _WRITE_OR_EXEC_OPTS:
        return True
    if tok.startswith("-") and "=" in tok and not tok.startswith("--"):
        name = tok.split("=", 1)[0]
        return name in _WRITE_OR_EXEC_OPTS
    # --index=..., --cached=... etc are not real, but fail closed anyway.
    # Capital-O pager bundle: git grep -O sh (two tokens) / -Osh (one).
    if tok.startswith("-O") and len(tok) >= 2:
        return True
    # Long options carrying an = value: --output=/tmp/x, --config-env=X.
    if tok.startswith("--") and "=" in tok:
        name = tok.split("=", 1)[0]
        return name in _WRITE_OR_EXEC_OPTS
    return False


def _no_escape(args, exempt: Iterable[str] = ()) -> bool:
    """True when no argv token is a write/execute option, apart from those in
    *exempt* (a per-shape safe reuse, e.g. ``-c`` for ``git switch -c``)."""
    exempt_set = set(exempt)
    return not any(
        tok in exempt_set and False
        for tok in args
    ) and not any(
        tok not in exempt_set and _token_is_output_or_run(tok)
        for tok in args
    )


def _compound_or_substitution(cmd: str) -> bool:
    """True when *cmd* is more than one simple command or hides a program.

    Conservative by design: a control/metacharacter ANYWHERE in the raw text
    (quotes included) fails closed. A shape is a single argv running exactly
    one program with no shell interpretation of its arguments, so even a
    literal ``;`` or ``$`` in a quoted pattern means we do not have that
    guarantee and the command is not a shape.
    """
    if "\n" in cmd:
        return True
    return (
        ";" in cmd or "&&" in cmd or "||" in cmd
        or "|" in cmd or ">" in cmd or "<" in cmd or "&" in cmd
        or "$(" in cmd or "`" in cmd or "$" in cmd
    )


def _resolve_inside(arg: str, workspace: str) -> bool:
    """True when *arg*, resolved against the WORKSPACE root (not the process
    cwd — a shape runs in the workspace), stays inside it."""
    try:
        root = Path(workspace).resolve(strict=False)
        p = Path(arg)
        if not p.is_absolute():
            p = root / p
        p = p.resolve(strict=False)
    except Exception:
        return False
    return p == root or root in p.parents


def _repo_relative(arg: str) -> bool:
    """True when *arg* is a git-tree path that stays inside the repo.

    ``git show <rev>:<path>`` and ``git grep <pat> <rev>:<path>`` name a path
    inside the repository's tree, not on the host filesystem. Such a path is
    repo-relative: not absolute and without a ``..`` escape.
    """
    if not arg or arg.startswith("/"):
        return False
    parts = arg.split("/")
    return ".." not in parts


def _split_one(cmd: str):
    """Return the parsed argv, or None when *cmd* is not a single clean
    command."""
    if _compound_or_substitution(cmd):
        return None
    try:
        argv = shlex.split(cmd)
    except ValueError:
        return None
    if not argv:
        return None
    return argv


_SHAPES: "list[tuple[str, object]]" = []


def _shape(name: str):
    def deco(fn):
        _SHAPES.append((name, fn))
        return fn
    return deco


@_shape("git_show_path")
def _git_show_path(argv, workspace: str) -> bool:
    # git [global-opts...] show <rev-or-rev:path>...
    if argv[0] != "git" or "show" not in argv:
        return False
    try:
        i = argv.index("show")
    except ValueError:
        return False
    rest = argv[i + 1:]
    if not rest:
        return False
    # global options before the subcommand: --no-pager / --paginate only
    for g in argv[1:i]:
        if g not in ("--no-pager", "--paginate"):
            return False
    if not _no_escape(argv):
        return False
    for a in rest:
        if a.startswith("-"):
            if a in ("--stat", "--name-only", "--no-color", "-p",
                     "--name-status", "--format"):
                if a == "--format":
                    # --format can carry an external command only via the
                    # write set; keep it but it is read-only in show.
                    continue
                continue
            return False
        if ":" in a:
            # rev:path -> the path must stay inside the repo
            _rev, sep, path = a.partition(":")
            if not sep or not _repo_relative(path):
                return False
        elif not _resolve_inside(a, workspace):
            # a bare positional arg (a rev, or a `rev::path`-less path) must
            # still resolve inside the workspace: an absolute or `..`-escaping
            # token would otherwise ride along unvalidated (reviewer-s27).
            return False
    return True


@_shape("git_grep")
def _git_grep(argv, workspace: str) -> bool:
    # git [-C dir] grep [options] <pattern> [--] [<path>...]
    if argv[0] != "git" or "grep" not in argv:
        return False
    try:
        i = argv.index("grep")
    except ValueError:
        return False
    rest = argv[i + 1:]
    if not rest:
        return False
    # global options before the subcommand: --no-pager / --paginate only
    # (reject -C dir: it points git at another repository)
    for g in argv[1:i]:
        if g not in ("--no-pager", "--paginate"):
            return False
    if not _no_escape(argv):
        return False
    # A short attached pattern-file option `-fFILE` is not in the token rule
    # (only -o and -O are bundled there); git grep has no other lowercase
    # -f* flag, so reject any argv token that is -f followed by a value.
    if any(a.startswith("-f") and len(a) > 2 and not a.startswith("--")
           for a in argv):
        return False
    # first non-option is the pattern; the rest may be revs or paths
    pat = None
    for idx, a in enumerate(rest):
        if a == "--":
            break
        if a.startswith("-"):
            continue
        pat = a
        rest_after = rest[idx + 1:]
        break
    if pat is None:
        # pattern supplied, e.g. -e PAT; no positional path required
        return True
    # optional -- then paths; each path must resolve inside the workspace
    if rest_after and rest_after[0] == "--":
        rest_after = rest_after[1:]
    for p in rest_after:
        if p.startswith("-"):
            continue
        if ":" in p:
            _rev, sep, path = p.partition(":")
            if not sep or not _repo_relative(path):
                return False
        else:
            if not _resolve_inside(p, workspace):
                return False
    return True


#: The closed set of options `git apply --check` may carry and still be the
#: read-only shape. Everything read-only and in-memory; see finding (s28 #3).
_APPLY_CHECK_READONLY_OPTS: frozenset[str] = frozenset({
    "--check",            # the marker itself (read-only)
    "--numstat",          # print a report, writes nothing
    "--summary",          # print a summary, writes nothing
    "-R", "--reverse",    # apply the diff reversed (in-memory)
    "--unidiff-zero",     # zero-context matching
    "--directory",        # re-interpret patched paths (in-memory)
    "--ignore-space-change", "--ignore-whitespace",  # matching tolerance
    "-z",                 # NUL-terminated report output
    "--verbose",
})


@_shape("git_apply_check")
def _git_apply_check(argv, workspace: str) -> bool:
    # git apply --check <path>   (apply --check reads, never writes)
    if len(argv) < 4 or argv[:2] != ["git", "apply"]:
        return False
    if "--check" not in argv:
        return False
    # global write/exec veto, then mutator options even alongside --check
    if not _no_escape(argv):
        return False
    if any(a in ("--index", "--cached", "--3way", "-3", "--reject")
           for a in argv[2:]):
        return False
    # A shape is "exactly one known read-only form": any option outside the
    # closed whitelist of genuinely read-only apply options makes the command
    # not that form (reviewer-s28 finding 3: `--tool=x` matched although it is
    # not a git apply option). Fail closed — such a command still goes through
    # the normal approval dialog.
    for a in argv[2:]:
        if a.startswith("-"):
            name = a.split("=", 1)[0] if "=" in a else a
            if name not in _APPLY_CHECK_READONLY_OPTS:
                return False
    # the positional patch path(s): exactly one, inside the workspace
    pos = [a for a in argv[2:] if not a.startswith("-")]
    if len(pos) != 1:
        return False
    return _resolve_inside(pos[0], workspace)


@_shape("git_hash_object")
def _git_hash_object(argv, workspace: str) -> bool:
    # git hash-object <path>   (no -w: hashes without writing an object)
    if len(argv) < 3 or argv[:2] != ["git", "hash-object"]:
        return False
    if not _no_escape(argv, exempt=()):
        return False
    if "-w" in argv or "--write" in argv:
        return False
    pos = [a for a in argv[2:] if not a.startswith("-")]
    if len(pos) != 1:
        return False
    return _resolve_inside(pos[0], workspace)


@_shape("git_switch_pr")
def _git_switch_pr(argv, workspace: str) -> bool:
    # git switch -c pr/<name> origin/main   creates a LOCAL branch -- this
    # is a WRITE that is approved by rule (it writes a branch ref, moves
    # HEAD and the work tree), NOT a read. It runs no external program and
    # never touches a host path outside the workspace. The <name> is
    # strictly bounded: [A-Za-z0-9._-]+ with no leading dash and no "..",
    # so pr/.. cannot escape or smuggle an option.
    if argv[:2] != ["git", "switch"]:
        return False
    rest = argv[2:]
    if "-c" not in rest:
        return False
    # exactly: -c pr/<name> origin/main, nothing before or after
    i = rest.index("-c")
    if i != 0:
        return False
    after = rest[i + 1:]
    if len(after) != 2:
        return False
    branch, origin_ref = after
    if not branch.startswith("pr/"):
        return False
    name = branch[3:]
    if not name or any(c not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-" for c in name):
        return False
    # no leading dash and no empty/traversal component for the name itself
    if name.startswith("-") or ".." in name:
        return False
    if origin_ref != "origin/main":
        return False
    # no other tokens anywhere
    if len(rest) != 3:
        return False
    return True


@_shape("gh_version")
def _gh_version(argv, workspace: str) -> bool:
    # gh --version  /  gh version
    if len(argv) == 2 and argv[0] == "gh" and argv[1] in ("--version", "version"):
        return True
    return False


#: Ordered shape names (registration order).
SHAPE_NAMES: tuple[str, ...] = tuple(name for name, _ in _SHAPES)


def matches(command: str, workspace: str) -> str:
    """Return the name of the read-only shape *command* matches, or ``""``."""
    return matches_with(command, workspace, None)


def matches_with(command: str, workspace: str,
                 enabled: Iterable[str] | None) -> str:
    """As :func:`matches`, but only shapes *enabled* may fire. ``None`` (or
    an omitted ``enabled``) means every shape is eligible."""
    allowed = set(enabled) if enabled is not None else set(SHAPE_NAMES)
    argv = _split_one(command)
    if argv is None:
        return ""
    for name, fn in _SHAPES:
        if name not in allowed:
            continue
        try:
            if fn(argv, workspace):
                return name
        except Exception:
            continue
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
    widen the gate.
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
    if not isinstance(wanted, list) or not wanted:
        shapes = list(SHAPE_NAMES)
    else:
        shapes = [n for n in wanted if n in unknown]
        if not shapes:
            shapes = list(SHAPE_NAMES)
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
    ] or list(SHAPE_NAMES)
    return payload
