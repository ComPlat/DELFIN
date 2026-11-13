"""What removes the way back: git history and other people's jobs.

Two pure readings of a shell command, used by the bash gate
(``api_client._run_permission_gate``):

:func:`git_history_refusal`
    Input: a command line (commit-message prose already blanked). Output:
    "" when every git invocation in it leaves the history recoverable,
    else the refusal text with the route that works.

    Recoverable means: what a branch pointed at before can be found again
    through the reflog (``git reflog``, ``<branch>@{1}``). Moving a branch
    (``branch -f``, a soft or mixed ``reset``, ``rebase``, ``commit
    --amend``, ``checkout -B``, ``switch -C``, ``update-ref`` of a branch)
    writes a reflog entry and stays free. Refused are the acts that remove
    that record or act where it does not exist:

    * a push that replaces or deletes refs on the remote -- the remote's
      old state is in no reflog of this repository: ``--force-with-lease``
      (``--force`` and ``-f`` are on the deny list), ``--force-if-includes``,
      a ``+`` refspec, ``--mirror``, ``--prune``, ``--delete``/``-d``, a
      ``:dst`` refspec, a short-option bundle carrying ``f`` or ``d``, and
      any abbreviation git accepts of those long options;
    * ``reflog expire|delete|drop``, ``gc --prune=<date>`` and ``prune``
      with an expiry shorter than a day (``now``, ``all``) or none;
    * the reflog switched off or shortened: ``core.logAllRefUpdates``
      false, ``gc.reflogExpire*`` / ``gc.pruneExpire`` below a day, set
      for one command (``-c``, ``--config-env``, ``GIT_CONFIG_*``) or for
      good (``git config``); a forced push configured
      (``remote.<r>.push`` with ``+``, ``remote.<r>.mirror``);
    * moving a tag (``tag -f``, ``update-ref refs/tags/...``): git keeps
      no reflog for tags, so the old target is gone;
    * ``update-ref --stdin``: deletes and moves read from a file the gate
      does not see.

    Rejected: a backup ref (``refs/delfin-backup/...``) written before each
    ref move. With the reflog protected it records nothing the reflog
    does not already hold for 90 days (30 for unreachable entries, git's
    defaults), it would have to be written by the gate before the command
    runs -- a write into ``.git`` on every move, refused by the cage's
    read-only git metadata on some hosts -- and it is pushed by
    ``--mirror`` and listed by ``--all`` everywhere afterwards. Tags,
    where it would add something, are refused instead.

:func:`foreign_job_cancels` / :func:`submitted_job_ids`
    A cancel (``scancel``, ``qdel``, ``bkill``) whose every target is a
    literal job id this session submitted is the session undoing its own
    work. Anything else -- another id, an id from a variable or a
    substitution, a cancel by filter with no id, ``xargs scancel``,
    ``bkill 0``, an option this module does not know -- may end a job
    somebody else is waiting on, and is returned as foreign. The ledger
    is filled from the output of a command that ran a submit program
    (``sbatch``, ``qsub``, ``bsub``), never from a file the agent can
    write.

The words of each simple command come from ``api_client._command_words``
(segments, ``sh -c`` bodies, wrappers such as ``env``, ``timeout`` and
``xargs`` stepped over), so the same spellings are read as by the other
command rules.
"""
from __future__ import annotations

import re
from typing import Iterable, Optional

# -- git ----------------------------------------------------------------------

#: Push options that replace or delete refs on the remote. git accepts any
#: unambiguous prefix of a long option, so a token is matched by prefix.
_PUSH_DESTRUCTIVE_LONG = ("force", "force-with-lease", "force-if-includes",
                          "mirror", "prune", "delete")
#: Push options that carry a value in the next word.
_PUSH_VALUE_OPTS = {"-o", "--push-option", "--repo", "--receive-pack",
                    "--exec"}

_FALSE = {"false", "no", "off", "0", ""}
_TRUE = {"true", "yes", "on", "1", "always"}

#: An expiry git keeps for at least a day: "never", or N units of a day or
#: longer (``2.weeks.ago``, ``90.days``).
_KEEPING_EXPIRY = re.compile(
    r"never|(?:[1-9]\d*)[. ](?:days?|weeks?|months?|years?)(?:[. ]ago)?",
    re.IGNORECASE)

_PROTECTED_KEY_RE = re.compile(
    r"core\.logallrefupdates|gc\.(?:[^=\s'\"]+\.)?reflogexpire(?:unreachable)?"
    r"|gc\.pruneexpire|remote\.[^=\s'\"]+\.(?:push|mirror)", re.IGNORECASE)

_ENV_CONFIG_RE = re.compile(
    r"\bGIT_CONFIG_(?:PARAMETERS|KEY_\d+)\s*=\s*['\"]?[^;&|\n]*?"
    r"(core\.logallrefupdates|gc\.[\w.*-]*reflogexpire\w*|gc\.pruneexpire"
    r"|remote\.[\w.-]+\.(?:push|mirror))", re.IGNORECASE)

_REMOTE_ALTERNATIVE = (
    " Push the work under a NEW branch name instead "
    "(`git push <remote> HEAD:refs/heads/<new-name>`), open a pull request "
    "from it if the work should land, and leave replacing or deleting the "
    "existing branch to the user; say so in your report.")

_REFLOG_ALTERNATIVE = (
    " The reflog is what makes a moved branch, a rebase or an amended "
    "commit recoverable (`git reflog`, `git reset --keep <branch>@{1}`). "
    "Keep it on and leave garbage collection to git's own schedule (a "
    "plain `git gc`, or nothing); if the repository really needs a "
    "clean-up, name it in your report for the user.")

_TAG_ALTERNATIVE = (
    " Git keeps no reflog for tags, so the old target would be lost. "
    "Create a new tag name instead and leave moving the existing tag to "
    "the user.")


def _keeps(value: Optional[str]) -> bool:
    return value is not None and bool(_KEEPING_EXPIRY.fullmatch(value.strip()))


def _config_refusal(key: str, value: Optional[str]) -> str:
    """Why setting *key* to *value* removes the way back, or ""."""
    k = key.strip().lower()
    v = None if value is None else value.strip().strip("'\"")
    if k == "core.logallrefupdates":
        if v is None or v.lower() in _TRUE:
            return ""
        return (f"setting {key}={v} switches the reflog off, so a later ref "
                "move cannot be undone." + _REFLOG_ALTERNATIVE)
    if re.fullmatch(r"gc\.(?:.+\.)?reflogexpire(?:unreachable)?|gc\.pruneexpire",
                    k):
        if _keeps(v):
            return ""
        return (f"setting {key}={v} lets git drop reflog entries or "
                "unreachable commits at once." + _REFLOG_ALTERNATIVE)
    if re.fullmatch(r"remote\..+\.push", k) and (v or "").startswith("+"):
        return (f"setting {key}={v} makes every push a forced push."
                + _REMOTE_ALTERNATIVE)
    if re.fullmatch(r"remote\..+\.mirror", k) and (
            v is None or v.lower() not in _FALSE):
        return (f"setting {key}={v} makes every push a mirror push, which "
                "deletes remote branches." + _REMOTE_ALTERNATIVE)
    return ""


def _split_kv(text: str) -> tuple[str, Optional[str]]:
    if "=" in text:
        k, v = text.split("=", 1)
        return k, v
    return text, None


def _global_config(rest: list[str]) -> tuple[str, int]:
    """(refusal, index of the subcommand) for git's own options in *rest*."""
    i = 0
    while i < len(rest):
        arg = rest[i]
        if not arg.startswith("-"):
            return "", i
        if arg == "-c" and i + 1 < len(rest):
            why = _config_refusal(*_split_kv(rest[i + 1]))
            if why:
                return why, i
            i += 2
            continue
        if arg.startswith("-c") and len(arg) > 2 and "=" in arg:
            why = _config_refusal(*_split_kv(arg[2:]))
            if why:
                return why, i
            i += 1
            continue
        if arg == "--config-env" or arg.startswith("--config-env="):
            spec = (arg.split("=", 1)[1] if "=" in arg
                    else (rest[i + 1] if i + 1 < len(rest) else ""))
            key = spec.split("=", 1)[0]
            if _PROTECTED_KEY_RE.fullmatch(key):
                return (f"--config-env sets {key} from the environment, "
                        "where the gate cannot read its value."
                        + _REFLOG_ALTERNATIVE), i
            i += 1 if "=" in arg else 2
            continue
        if arg in ("-C", "--git-dir", "--work-tree", "--namespace",
                   "--exec-path") and "=" not in arg:
            i += 2
            continue
        i += 1
    return "", len(rest)


def _push_refusal(args: list[str]) -> str:
    positional: list[str] = []
    i = 0
    while i < len(args):
        a = args[i]
        if a == "--":
            positional.extend(args[i + 1:])
            break
        if a.startswith("--"):
            name = a[2:].split("=", 1)[0]
            if name.startswith("no-"):
                i += 1
                continue
            if name and any(full.startswith(name)
                            for full in _PUSH_DESTRUCTIVE_LONG):
                hit = next(f for f in _PUSH_DESTRUCTIVE_LONG
                           if f.startswith(name))
                return (f"`git push --{name}` ({hit}) replaces or deletes "
                        "refs on the remote, where no reflog of this "
                        "repository can bring the old state back."
                        + _REMOTE_ALTERNATIVE)
            if a in _PUSH_VALUE_OPTS:
                i += 2
                continue
            i += 1
            continue
        if a.startswith("-") and len(a) > 1:
            for ch in a[1:]:
                if ch in "fd":
                    what = "forces" if ch == "f" else "deletes"
                    return (f"`git push {a}` {what} a ref on the remote, "
                            "where no reflog of this repository can bring "
                            "the old state back." + _REMOTE_ALTERNATIVE)
                if ch == "o":
                    if a.endswith("o"):
                        i += 1          # the value is the next word
                    break
            i += 1
            continue
        positional.append(a)
        i += 1
    for spec in positional[1:]:
        if spec.startswith("+"):
            return (f"the refspec `{spec}` forces the update on the remote, "
                    "where no reflog of this repository can bring the old "
                    "state back." + _REMOTE_ALTERNATIVE)
        if spec.startswith(":"):
            return (f"the refspec `{spec}` deletes a ref on the remote."
                    + _REMOTE_ALTERNATIVE)
    return ""


def _config_cmd_refusal(args: list[str]) -> str:
    value_opts = {"-f", "--file", "--blob", "--type", "--default",
                  "--comment", "--value"}
    reading = {"--get", "--get-all", "--get-regexp", "--get-urlmatch",
               "--list", "-l", "--unset", "--unset-all", "--edit", "-e",
               "--remove-section", "--rename-section", "--show-origin",
               "--show-scope", "--name-only"}
    positional: list[str] = []
    i = 0
    while i < len(args):
        a = args[i]
        if a in reading:
            return ""
        if a.startswith("-"):
            i += 2 if (a in value_opts and "=" not in a) else 1
            continue
        positional.append(a)
        i += 1
    if positional and positional[0] == "set":
        positional = positional[1:]
    elif positional and positional[0] in ("get", "list", "unset", "edit",
                                          "rename-section", "remove-section"):
        return ""
    if len(positional) >= 2:
        return _config_refusal(positional[0], positional[1])
    return ""


def _option_value(args: list[str], name: str) -> Optional[str]:
    for i, a in enumerate(args):
        if a.startswith(name + "="):
            return a.split("=", 1)[1]
        if a == name and i + 1 < len(args):
            return args[i + 1]
    return None


def _subcommand_refusal(sub: str, args: list[str]) -> str:
    if sub == "push":
        return _push_refusal(args)
    if sub == "config":
        return _config_cmd_refusal(args)
    if sub == "reflog":
        action = next((a for a in args if not a.startswith("-")), "")
        if action in ("expire", "delete", "drop"):
            return (f"`git reflog {action}` removes reflog entries."
                    + _REFLOG_ALTERNATIVE)
        return ""
    if sub == "gc":
        for a in args:
            if a.startswith("--prune="):
                if not _keeps(a.split("=", 1)[1]):
                    return (f"`git gc {a}` deletes unreachable commits at "
                            "once." + _REFLOG_ALTERNATIVE)
        return ""
    if sub == "prune":
        if any(a in ("-n", "--dry-run") for a in args):
            return ""
        if _keeps(_option_value(args, "--expire")):
            return ""
        return ("`git prune` deletes every unreachable object, with no "
                "grace period unless --expire names one." + _REFLOG_ALTERNATIVE)
    if sub == "tag":
        for a in args:
            if a == "--force" or (a.startswith("-") and not a.startswith("--")
                                  and "f" in a[1:].split("m")[0]
                                  .split("F")[0].split("u")[0]):
                return f"`git tag {a}` moves an existing tag." + _TAG_ALTERNATIVE
        return ""
    if sub == "update-ref":
        if "--stdin" in args:
            return ("`git update-ref --stdin` moves or deletes refs listed "
                    "in input the gate does not see. Move one branch at a "
                    "time instead (`git branch -f <branch> <commit>`); the "
                    "reflog keeps each move recoverable.")
        positional = [a for a in args if not a.startswith("-")]
        if (positional and positional[0].startswith("refs/tags/")
                and "--create-reflog" not in args):
            return (f"`git update-ref {positional[0]}` moves a tag."
                    + _TAG_ALTERNATIVE)
        return ""
    return ""


def git_history_refusal(command: str) -> str:
    """The refusal for a git act that removes the way back, or ""."""
    if "git" not in (command or "").lower():
        return ""
    m = _ENV_CONFIG_RE.search(command or "")
    if m:
        return (f"the command sets `{m.group(1)}` through the environment "
                "(GIT_CONFIG_*), where it can switch the reflog off or "
                "shorten it." + _REFLOG_ALTERNATIVE)
    from .api_client import _command_words
    for words in _command_words(command):
        if not words or words[0] != "git":
            continue
        rest = words[1:]
        why, at = _global_config(rest)
        if why:
            return why
        if at >= len(rest):
            continue
        why = _subcommand_refusal(rest[at], rest[at + 1:])
        if why:
            return why
    return ""


# -- jobs ---------------------------------------------------------------------

_CANCEL_PROGRAMS = ("scancel", "qdel", "bkill")
_SUBMIT_PROGRAMS = ("sbatch", "qsub", "bsub")

#: scancel options that carry a value in the next word. With at least one
#: job id named, every scancel filter narrows the set further; without one
#: a filter selects jobs by name, partition, state ... -- not provably own.
_SCANCEL_VALUE_OPTS = {"-A", "--account", "-n", "--name", "--jobname",
                       "-p", "--partition", "-q", "--qos", "-R",
                       "--reservation", "-s", "--signal", "-t", "--state",
                       "-w", "--nodelist", "-M", "--clusters", "--wckey"}
_SCANCEL_FLAGS = {"-b", "--batch", "-f", "--full", "-Q", "--quiet", "-v",
                  "--verbose", "-i", "--interactive", "-H", "--hurry",
                  "--me", "--ctld", "--sibling"}
_QDEL_VALUE_OPTS = {"-W"}
_QDEL_FLAGS = {"-x"}
_BKILL_VALUE_OPTS = {"-s"}
_BKILL_FLAGS = {"-b", "-l"}

_JOB_ID_SHAPES = {
    "scancel": re.compile(r"(\d+)(?:_(?:\d+|\[[\d,%-]+\]))?(?:\.\w+)?"),
    "qdel": re.compile(r"(\d+)(?:\[\d*\])?(?:\.[\w.-]+)?"),
    "bkill": re.compile(r"([1-9]\d*)(?:\[[\d,-]+\])?"),
}


def _cancel_targets(program: str, args: list[str]) -> Optional[list[str]]:
    """The base job ids *args* name, or None when the set is not literal."""
    value_opts, flags = {
        "scancel": (_SCANCEL_VALUE_OPTS, _SCANCEL_FLAGS),
        "qdel": (_QDEL_VALUE_OPTS, _QDEL_FLAGS),
        "bkill": (_BKILL_VALUE_OPTS, _BKILL_FLAGS),
    }[program]
    ids: list[str] = []
    i = 0
    while i < len(args):
        a = args[i]
        if a.startswith("-") and len(a) > 1:
            name = a.split("=", 1)[0]
            if name in value_opts:
                i += 1 if "=" in a else 2
                continue
            if a in flags:
                i += 1
                continue
            return None                 # unknown option: not provable
        m = _JOB_ID_SHAPES[program].fullmatch(a)
        if not m:
            return None                 # a variable, a glob, a name, "all"
        ids.append(m.group(1))
        i += 1
    return ids or None


def foreign_job_cancels(command: str, own: Iterable[str]) -> list[str]:
    """The cancels in *command* not provably of this session's own jobs.

    Output: one short description per such cancel (``"scancel 4712"``);
    empty when there is no cancel or every target is in *own*."""
    low = (command or "").lower()
    if not any(p in low for p in _CANCEL_PROGRAMS):
        return []
    from .api_client import _command_words
    mine = {str(j).strip() for j in own or ()}
    out: list[str] = []
    for words in _command_words(command):
        if not words or words[0] not in _CANCEL_PROGRAMS:
            continue
        ids = _cancel_targets(words[0], words[1:])
        shown = " ".join(words)[:120]
        if ids is None or any(j not in mine for j in ids):
            out.append(shown)
    return out


_SLURM_SUBMITTED = re.compile(r"Submitted batch job\s+(\d+)")
_SLURM_PARSABLE = re.compile(r"^(\d+)(?:;[\w.-]+)?$")
_PBS_ID = re.compile(r"^(\d+)(?:\.[\w.-]+)?$")
_LSF_SUBMITTED = re.compile(r"Job <(\d+)> is submitted")


def submitted_job_ids(command: str, stdout: str) -> list[str]:
    """The base ids of the jobs *command* submitted, read from its output.

    Only a command whose every program is a submit program (or ``cd``;
    wrappers such as ``env`` and ``timeout`` are stepped over by
    ``_command_words``) counts: any other program in the line -- an
    ``echo``, a ``printf``, a ``cat`` of an old log, a script, a shell
    body -- could print a submit line for an id this session never
    submitted, so such a line records nothing. A missed own id only
    makes a later cancel of it asked or refused (fail closed)."""
    try:
        from .api_client import _command_words
        programs = {w[0] for w in _command_words(command or "") if w}
    except Exception:
        return []
    if (not programs & set(_SUBMIT_PROGRAMS)
            or programs - set(_SUBMIT_PROGRAMS) - {"cd"}):
        return []
    found: list[str] = []
    text = stdout or ""
    if "sbatch" in programs:
        found += _SLURM_SUBMITTED.findall(text)
        if "--parsable" in (command or ""):
            found += [m.group(1) for line in text.splitlines()
                      if (m := _SLURM_PARSABLE.match(line.strip()))]
    if "qsub" in programs:
        found += [m.group(1) for line in text.splitlines()
                  if (m := _PBS_ID.match(line.strip()))]
    if "bsub" in programs:
        found += _LSF_SUBMITTED.findall(text)
    out: list[str] = []
    for jid in found:
        if jid not in out:
            out.append(jid)
    return out
