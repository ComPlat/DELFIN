"""What an approval request is about, in words the person reads.

Input: one call the gate is about to put to a person -- a tool name, its
arguments, the directory a command runs in, the session's permissions.
Output: ``classify`` returns a :class:`Consequence` of one of three levels

* ``reversible``: nothing is destroyed or sent, or a copy is kept (the undo
  journal, git, the trash);
* ``recoverable``: undone with effort -- the replaced commits stay in git's
  reflog;
* ``irreversible``: it deletes or overwrites what nothing keeps, publishes
  (a push, a pull request, a release), sends data or a message off the
  machine, cancels a job, ends a process, or acts in the lab;

with the facts that make it so, one concrete sentence each (which file and
how big, which branch on which remote and how many commits, which host and
how many bytes, which job), and the reversible way to the same end where
there is one (the trash instead of a delete, a new branch instead of a
forced push).

``statement`` writes the facts as the block an approval opens with
(``with_statement`` puts it at the head of a preview); ``split`` reads that
block back off a preview as ``(template, values)`` pairs. The block is the
transport: every surface receives ``(tool, args, preview)`` -- the terminal
broker, the dashboard panel, the service and the phone bridge -- and a new
argument would have to pass five brokers and a socket protocol. The English
template is also the key of its German text (``messenger.lang.DE``), so the
phone says the same sentence in the person's language while the terminal and
the dashboard, which are English, show it as written. A line that matches no
template is shown as it stands.

Semantics the surfaces apply to an ``irreversible`` statement: no grant that
outlives the one action is offered (``always``, ``commands like this``,
``stop asking about writes``, ``Allow + Permanent``), the default answer is
no, and an allow pattern, an allow rule or the cage does not answer it in
place of the person (api_client).

Rejected: a second list of dangerous commands. The facts come from what
already decides: ``unrestorable`` (what nothing could restore), the gate's
own push parsing (``api_client._push_targets``), auto_guard's names for
tools that send. Rejected: classifying from the preview text on each
surface -- the preview is prose for a person, and the gate is the one place
that knows the directory a command runs in.

Not seen: what an interpreter or a script does inside (``python -c``,
``./clean.sh``), a target held in a variable or a substitution, a write
tool's overwrite (the change journal keeps it), and an MCP tool whose name
does not say that it sends. Such a request is asked as before, without a
statement.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import string
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Optional

REVERSIBLE = "reversible"
RECOVERABLE = "recoverable"
IRREVERSIBLE = "irreversible"
_RANK = {REVERSIBLE: 0, RECOVERABLE: 1, IRREVERSIBLE: 2}

#: The first line of a statement; it decides the level ``split`` reads.
HEAD_IRREVERSIBLE = "⛔ This cannot be undone."
HEAD_RECOVERABLE = "↺ This can be undone only with effort."
BULLET = "• "

# One sentence each; ``{name}`` is a value, never English prose, so the
# German text (messenger/lang.DE) can be filled with the same values.
DELETE = ("Deletes {path} ({size}). Nothing keeps a copy: the undo journal "
          "cannot hold it, git does not have it, and this session did not "
          "write it.")
OVERWRITE = ("Overwrites {path} ({size}). Nothing keeps a copy of what is in "
             "it now: the undo journal cannot hold it, git does not have it, "
             "and this session did not write it.")
MORE_FILES = "And {n} more file(s) that nothing keeps."
TRASH = ("Instead: move it to the trash, which keeps it {days} days: "
         "mv {path} {trash}/")
PUSH = ("Publishes branch {branch} to {remote} ({where}): {n} commit(s) it "
        "does not have yet. Whoever can read it can fetch them at once; "
        "deleting them later does not take back the copies already fetched.")
PUSH_NEW = ("Publishes the new branch {branch} on {remote} ({where}). Whoever "
            "can read it can fetch it at once; deleting it later does not "
            "take back the copies already fetched.")
FORCE = ("Replaces the history of branch {branch} on {remote} ({where}): "
         "{n} commit(s) there that your branch does not contain are removed "
         "from it.")
FORCE_UNKNOWN = ("Replaces the history of branch {branch} on {remote} "
                 "({where}): commits there that your branch does not contain "
                 "are removed from it.")
NEW_BRANCH = ("Instead: push to a new branch and leave {branch} as it is: "
              "git push {remote} HEAD:{branch}-2")
PR_CREATE = ("Opens a pull request for branch {branch}: its commits and its "
             "description become visible to everyone who can read the "
             "repository.")
PR_MERGE = "Merges pull request {pr} into its base branch on the forge."
RELEASE = ("Publishes release {tag}: whoever can read the repository can "
           "download it at once.")
REPO_DELETE = ("Deletes the repository {repo} on the forge, with its issues "
               "and pull requests.")
API_WRITE = ("Sends a {method} request to the forge's API ({path}): it "
             "changes data there.")
HTTP_WRITE = ("Sends a {method} request to {host}: what it changes or deletes "
              "there cannot be taken back from here.")
UPLOAD = ("Sends {n} file(s), {size} ({files}), to {host}: the data leaves "
          "this machine, and the copy there is out of DELFIN's reach.")
UPLOAD_DATA = ("Sends data to {host}: it leaves this machine, and the copy "
               "there is out of DELFIN's reach.")
MAIL = "Sends an e-mail to {to}: once sent it cannot be called back."
MESSAGE = "Sends a message to {to}: once delivered it cannot be called back."
NOTIFY = ("Sends a notification to the paired phone: once delivered it "
          "cannot be called back.")
WEBHOOK = ("Calls {url}: what the other side does with the request cannot be "
           "taken back.")
CANCEL = ("Cancels job {jobs}: a running calculation stops at once, and its "
          "work since the last checkpoint is lost.")
CANCEL_ON = ("Cancels job {jobs} on {cluster}: a running calculation stops at "
             "once, and its work since the last checkpoint is lost.")
KILL = ("Ends process {pid} ({name}): what it holds in memory and has not "
        "written yet is lost.")
KILL_MATCH = ("Ends every process whose name matches {pattern}: what they hold "
              "in memory and have not written yet is lost.")
DISCARD = ("Discards the uncommitted changes to {n} file(s) ({files}): git "
           "keeps no copy of them.")
SAVE_PATCH = "Instead: save them first: git diff > before-discard.patch"
CLEAN = ("Deletes {n} untracked file(s) ({files}): git has never had a copy "
         "of them.")
MCP = ("Calls {tool} of the MCP server {server}: what it sends, creates or "
       "deletes there is out of DELFIN's reach.")
LAB = ("Sets {device}.{channel} to {value} in the lab: a physical action, "
       "which setting it back does not undo.")
REFLOG = ("Rewrites the history of branch {branch}: the commits it replaces "
          "stay in git's reflog for about 90 days (git reflog lists them, "
          "git reset brings one back).")

#: Every text a statement is made of, the two heads first. The phone looks
#: each one up in messenger/lang.DE (a test holds every one to an entry).
TEMPLATES: tuple[str, ...] = (
    HEAD_IRREVERSIBLE, HEAD_RECOVERABLE, DELETE, OVERWRITE, MORE_FILES, TRASH,
    PUSH, PUSH_NEW, FORCE, FORCE_UNKNOWN, NEW_BRANCH, PR_CREATE, PR_MERGE,
    RELEASE, REPO_DELETE, API_WRITE, HTTP_WRITE, UPLOAD, UPLOAD_DATA, MAIL,
    MESSAGE, NOTIFY, WEBHOOK, CANCEL, CANCEL_ON, KILL, KILL_MATCH, DISCARD,
    SAVE_PATCH, CLEAN, MCP, LAB, REFLOG,
)

#: A command line lists at most this many names; the count says the rest.
_SHOWN = 3
#: Bounds of the walk that sizes an upload's directory.
_WALK_MAX = 5_000
_WALK_S = 1.0


@dataclass(frozen=True)
class Fact:
    """One sentence: its template and the values filled in."""
    template: str
    values: tuple = ()

    @property
    def params(self) -> dict:
        return dict(self.values)

    def text(self, translate: Optional[Callable[[str], str]] = None) -> str:
        tpl = translate(self.template) if translate else self.template
        try:
            return tpl.format(**self.params)
        except (KeyError, IndexError, ValueError):
            return self.template.format(**self.params) if self.values \
                else self.template


@dataclass(frozen=True)
class Consequence:
    level: str = REVERSIBLE
    facts: tuple = ()
    alternative: Optional[Fact] = None

    @property
    def irreversible(self) -> bool:
        return self.level == IRREVERSIBLE


@dataclass
class _Acc:
    level: str = REVERSIBLE
    facts: list = field(default_factory=list)
    alternative: Optional[Fact] = None

    def add(self, level: str, template: str, alternative: Optional[Fact] = None,
            **values) -> None:
        if _RANK[level] > _RANK[self.level]:
            self.level = level
        self.facts.append(_fact(template, **values))
        if alternative is not None and self.alternative is None:
            self.alternative = alternative

    def done(self) -> Consequence:
        if self.level == REVERSIBLE:
            return Consequence()
        return Consequence(self.level, tuple(self.facts[:4]), self.alternative)


def _clean(value: Any) -> str:
    """A value as one line of plain text: control characters shown as '?'."""
    text = str(value)
    return "".join("?" if (ord(c) < 32 or ord(c) == 127) else c for c in text)


def _fact(template: str, **values) -> Fact:
    return Fact(template, tuple(sorted((k, _clean(v)) for k, v in values.items())))


# --------------------------------------------------------------------------
# The statement: written at the head of a preview, read back off it
# --------------------------------------------------------------------------

def statement(cons: Consequence) -> str:
    """The block an approval of ``cons`` opens with; "" when reversible."""
    if cons is None or cons.level == REVERSIBLE or not cons.facts:
        return ""
    head = HEAD_IRREVERSIBLE if cons.irreversible else HEAD_RECOVERABLE
    facts = list(cons.facts) + ([cons.alternative] if cons.alternative else [])
    return "\n".join([head] + [BULLET + f.text() for f in facts])


def with_statement(preview: str, cons: Consequence) -> str:
    """``preview`` with the statement of ``cons`` at its head (once)."""
    block = statement(cons)
    if not block or split(preview)[0] is not None:
        return preview
    return f"{block}\n\n{preview}" if preview else block


def _pattern(template: str) -> re.Pattern:
    parts, seen = [], set()
    for literal, name, _spec, _conv in string.Formatter().parse(template):
        parts.append(re.escape(literal))
        if name:
            parts.append(f"(?P={name})" if name in seen else f"(?P<{name}>.+?)")
            seen.add(name)
    return re.compile("^" + "".join(parts) + "$", re.S)


_PATTERNS = [(t, _pattern(t)) for t in TEMPLATES[2:]]


@dataclass(frozen=True)
class Statement:
    """A statement read back off a preview."""
    level: str
    lines: tuple            # (template, values dict) or (raw line, None)

    @property
    def irreversible(self) -> bool:
        return self.level == IRREVERSIBLE

    @property
    def head(self) -> str:
        return HEAD_IRREVERSIBLE if self.irreversible else HEAD_RECOVERABLE

    def sentences(self, translate: Optional[Callable[[str], str]] = None
                  ) -> list[str]:
        """The bullet sentences, each translated before it is filled."""
        out = []
        for tpl, values in self.lines:
            if values is None:
                out.append(tpl)
                continue
            out.append(Fact(tpl, tuple(sorted(values.items()))).text(translate))
        return out

    def render(self, translate: Optional[Callable[[str], str]] = None
               ) -> list[str]:
        head = translate(self.head) if translate else self.head
        return [head] + [BULLET + s for s in self.sentences(translate)]


def _parse_line(line: str):
    for tpl, rx in _PATTERNS:
        m = rx.match(line)
        if m:
            return tpl, m.groupdict()
    return line, None


def split(preview: str) -> tuple[Optional[Statement], str]:
    """``(statement, the rest of the preview)``; ``(None, preview)`` when
    the preview does not open with one. Never raises."""
    text = str(preview or "")
    first, _, rest = text.partition("\n")
    level = {HEAD_IRREVERSIBLE: IRREVERSIBLE,
             HEAD_RECOVERABLE: RECOVERABLE}.get(first.strip())
    if level is None:
        return None, text
    lines = rest.split("\n")
    items, i = [], 0
    while i < len(lines) and lines[i].startswith(BULLET):
        items.append(_parse_line(lines[i][len(BULLET):]))
        i += 1
    return Statement(level, tuple(items)), "\n".join(lines[i:]).lstrip("\n")


def is_irreversible(preview: str) -> bool:
    found = split(preview)[0]
    return bool(found and found.irreversible)


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------

_CACHE: dict = {}
_CACHE_LOCK = threading.Lock()
_CACHE_S = 5.0


def classify(tool: str, args: Optional[dict] = None, *, cwd: Any = None,
             perms: Any = None) -> Consequence:
    """The consequence of ``tool(args)`` run in ``cwd`` (default: the
    workspace of ``perms``). Never raises; a failure reads as reversible,
    which shows no statement and changes no answer.

    The gate asks twice for one call (whether a grant may answer it, then
    the dialog), so a result is kept for a few seconds per call."""
    args = args if isinstance(args, dict) else {}
    try:
        key = (str(tool), json.dumps(args, sort_keys=True, default=str),
               str(cwd), id(perms))
    except Exception:
        key = None
    now = time.monotonic()
    if key is not None:
        with _CACHE_LOCK:
            hit = _CACHE.get(key)
        if hit is not None and now - hit[0] < _CACHE_S:
            return hit[1]
    try:
        cons = _classify(str(tool or ""), args, cwd, perms)
    except Exception:
        cons = Consequence()
    if key is not None:
        with _CACHE_LOCK:
            if len(_CACHE) > 64:
                _CACHE.clear()
            _CACHE[key] = (now, cons)
    return cons


def _workdir(cwd: Any, perms: Any, args: dict) -> Path:
    base = Path(str(getattr(perms, "workspace", "") or os.getcwd()))
    if cwd:
        base = base / Path(str(cwd)).expanduser()
    elif args.get("cwd"):
        base = base / Path(str(args.get("cwd"))).expanduser()
    return base


def _classify(tool: str, args: dict, cwd: Any, perms: Any) -> Consequence:
    bare = tool.rsplit("__", 1)[-1] if tool.startswith("mcp__kit-coding__") \
        else tool
    acc = _Acc()
    if bare in ("bash", "bash_background"):
        _bash(acc, str(args.get("command") or ""), _workdir(cwd, perms, args),
              perms)
    elif bare == "cluster_cancel":
        jobs = args.get("job_ids") or []
        jobs = jobs if isinstance(jobs, list) else [jobs]
        if jobs:
            acc.add(IRREVERSIBLE, CANCEL_ON, jobs=_names(jobs),
                    cluster=args.get("cluster") or "?")
    elif bare == "cluster_upload":
        _tool_upload(acc, args, perms)
    elif bare in ("session_message", "subagent_message"):
        to = (args.get("to") or args.get("session_id") or args.get("target")
              or args.get("session") or "?")
        acc.add(IRREVERSIBLE, MESSAGE, to=to)
    elif bare == "push_notification":
        acc.add(IRREVERSIBLE, NOTIFY)
    elif bare == "remote_trigger":
        acc.add(IRREVERSIBLE, WEBHOOK, url=args.get("url") or args.get("name")
                or "?")
    elif tool.startswith("mcp__") and not tool.startswith("mcp__kit-coding__"):
        parts = tool.split("__")
        name = parts[-1] if len(parts) >= 3 else tool
        try:
            from .auto_guard import _MCP_SENDING
            sends = bool(_MCP_SENDING.search(name))
        except Exception:
            sends = False
        if sends:
            acc.add(IRREVERSIBLE, MCP, tool=name,
                    server=parts[1] if len(parts) >= 3 else "?")
    return acc.done()


def _names(items, shown: int = _SHOWN) -> str:
    items = [str(i) for i in items]
    head = ", ".join(items[:shown])
    return head + (", …" if len(items) > shown else "")


def _size_text(size: int) -> str:
    try:
        from .unrestorable import _size_text as _st
        return _st(int(size))
    except Exception:
        return f"{int(size)} B"


def _show(path: Path, base: Path) -> str:
    """``path`` relative to ``base`` when it lies under it."""
    try:
        return os.path.relpath(str(path), str(base)) \
            if str(Path(path).resolve()).startswith(str(Path(base).resolve()) + os.sep) \
            else str(path)
    except Exception:
        return str(path)


def _files_and_size(paths: list[Path]) -> tuple[list[str], int]:
    """Regular files under ``paths`` (bounded walk) and their total size."""
    names, total, seen = [], 0, 0
    deadline = time.monotonic() + _WALK_S
    for p in paths:
        try:
            if p.is_file():
                names.append(p.name)
                total += p.stat().st_size
                continue
            for dirpath, _dirs, files in os.walk(p, followlinks=False):
                for f in files:
                    seen += 1
                    if seen > _WALK_MAX or time.monotonic() > deadline:
                        return names, total
                    fp = Path(dirpath) / f
                    if fp.is_file() and not fp.is_symlink():
                        names.append(f)
                        total += fp.stat().st_size
        except OSError:
            continue
    return names, total


def _tool_upload(acc: _Acc, args: dict, perms: Any) -> None:
    src = str(args.get("src") or "")
    base = Path(str(getattr(perms, "workspace", "") or "."))
    names, size = _files_and_size([base / src]) if src else ([], 0)
    host = str(args.get("cluster") or "?")
    if names:
        acc.add(IRREVERSIBLE, UPLOAD, n=len(names), size=_size_text(size),
                files=_names(names), host=host)
    else:
        acc.add(IRREVERSIBLE, UPLOAD_DATA, host=host)


# ---- the shell -------------------------------------------------------------

def _segments(cmd: str) -> list[list[str]]:
    """Each simple command of ``cmd`` as words, wrappers dropped."""
    try:
        from .api_client import _split_shell_segments, _with_shell_bodies
        segs = _split_shell_segments(_with_shell_bodies(cmd))
    except Exception:
        segs = re.split(r"\|\||&&|[;|&\n]", cmd)
    from .unrestorable import _strip_prefix
    out = []
    for seg in segs:
        try:
            words = shlex.split(seg, comments=False)
        except ValueError:
            continue
        kept, skip = [], False
        for w in words:
            if skip:
                skip = False
            elif re.fullmatch(r"\d*(?:[<>]{1,2}|&>|>\||<<<?)", w):
                skip = True                     # the redirection's file
            elif not re.match(r"^\d*[<>]", w):
                kept.append(w)
        words = _strip_prefix(kept)
        if words:
            out.append(words)
    return out


def _bash(acc: _Acc, cmd: str, cwd: Path, perms: Any) -> None:
    if not cmd.strip():
        return
    _destroyed_files(acc, cmd, cwd, perms)
    for words in _segments(cmd):
        prog = os.path.basename(words[0])
        rest = words[1:]
        try:
            if prog == "git":
                _git(acc, rest, cwd)
            elif prog == "gh":
                _gh(acc, rest, cwd)
            elif prog in ("scp", "sftp", "rsync"):
                _remote_copy(acc, prog, rest, cwd)
            elif prog in ("curl", "wget"):
                _http(acc, prog, rest, cwd)
            elif prog in ("mail", "mailx", "sendmail", "mutt", "ssmtp", "msmtp"):
                to = [w for w in _operands(rest, {"-s", "-c", "-b", "-a", "-r",
                                                  "-F", "-f", "-S"})]
                acc.add(IRREVERSIBLE, MAIL, to=_names(to) or "?")
            elif prog in ("kill", "pkill", "killall"):
                _kill(acc, prog, rest)
            elif prog in ("scancel", "qdel", "bkill"):
                jobs = _operands(rest, {"-n", "--name", "-p", "--partition",
                                        "-A", "--account", "-M", "--clusters",
                                        "-s", "--signal", "-t", "--state",
                                        "-w", "--nodelist", "-q", "--qos"})
                if jobs:
                    acc.add(IRREVERSIBLE, CANCEL, jobs=_names(jobs))
        except Exception:
            continue


def _operands(args: list[str], takes_value: set) -> list[str]:
    from .unrestorable import _operands as _ops
    return _ops(args, takes_value)


def _destroyed_files(acc: _Acc, cmd: str, cwd: Path, perms: Any) -> None:
    if perms is None:
        # Without a session the directory is the workspace (a workspace
        # under /tmp is not scratch then, as for a real session there).
        from types import SimpleNamespace
        perms = SimpleNamespace(workspace=str(cwd),
                                all_workspace_roots=lambda: [str(cwd)])
    try:
        from . import unrestorable as U
        found = U.destroyed(cmd, cwd)
        bad = U.unrestorable(found, perms) if found else []
    except Exception:
        return
    if not bad:
        return
    verb, path, size = bad[0]
    template = OVERWRITE if verb in ("overwrite", "truncate") else DELETE
    try:
        root = U._root_of(path, perms) if perms is not None else Path(cwd)
    except Exception:
        root = Path(cwd)
    day = Path(root) / U.TRASH / datetime.now(timezone.utc).strftime("%Y-%m-%d")
    shown = _show(path, cwd)
    acc.add(IRREVERSIBLE, template, path=shown, size=size,
            alternative=_fact(TRASH, days=U.TRASH_DAYS,
                              path=shlex.quote(shown), trash=_show(day, cwd)))
    if len(bad) > 1:
        acc.add(IRREVERSIBLE, MORE_FILES, n=len(bad) - 1)


def _git_run(cwd: Path, *args: str) -> Optional[str]:
    try:
        from .api_client import _git_out
        return _git_out(cwd, *args)
    except Exception:
        return None


def _current_branch(cwd: Path) -> str:
    name = _git_run(cwd, "rev-parse", "--abbrev-ref", "HEAD") or ""
    return name if name and name != "HEAD" else "HEAD"


def _where(url: str) -> str:
    """A remote's location without credentials: host and path."""
    url = str(url or "").strip()
    m = re.match(r"^[a-z][a-z0-9+.-]*://(?:[^@/]*@)?([^/]+)(/.*)?$", url, re.I)
    if m:
        return m.group(1) + (m.group(2) or "")
    m = re.match(r"^(?:[^@/:]+@)?([A-Za-z0-9._-]+):(.+)$", url)
    if m and not url.startswith(("/", ".")):
        return f"{m.group(1)}:{m.group(2)}"
    return url or "?"


def _git(acc: _Acc, words: list[str], cwd: Path) -> None:
    i = 0
    while i < len(words) and words[i].startswith("-"):
        if words[i] in ("-C", "-c") and i + 1 < len(words):
            if words[i] == "-C":
                cwd = cwd / Path(words[i + 1]).expanduser()
            i += 2
            continue
        i += 1
    if i >= len(words):
        return
    sub, rest = words[i], words[i + 1:]
    if sub == "push":
        _push(acc, rest, cwd)
    elif sub in ("checkout", "restore"):
        _discard(acc, sub, rest, cwd)
    elif sub == "clean":
        _clean_untracked(acc, rest, cwd)
    elif (sub == "rebase" and not {"--abort", "--continue", "--skip", "--quit",
                                    "--show-current-patch"} & set(rest)) \
            or (sub == "commit" and "--amend" in rest) \
            or (sub == "reset" and rest and "--hard" not in rest
                and any(not w.startswith("-") for w in rest)):
        acc.add(RECOVERABLE, REFLOG, branch=_current_branch(cwd))


def _push(acc: _Acc, words: list[str], cwd: Path) -> None:
    try:
        from .api_client import _PUSH_VALUE_OPTS as value_opts
    except Exception:
        value_opts = {"-o", "--push-option", "--receive-pack", "--exec", "--repo"}
    force = False
    positional: list[str] = []
    i = 0
    while i < len(words):
        w = words[i]
        if w in value_opts:
            i += 2
            continue
        if w.startswith("--"):
            if w == "--dry-run":
                return
            if w in ("--force", "--mirror") or w.startswith("--force-with-lease"):
                force = True
        elif w.startswith("-") and len(w) > 1:
            if "f" in w[1:]:
                force = True
            if "n" in w[1:]:
                return                          # --dry-run, short form
        else:
            positional.append(w)
        i += 1
    current = _current_branch(cwd)
    if positional:
        remote = positional[0]
    else:
        remote = (_git_run(cwd, "config", f"branch.{current}.remote")
                  or "origin")
    url = _git_run(cwd, "remote", "get-url", remote)
    where = _where(url if url else remote)
    specs = positional[1:] or ["HEAD"]
    for spec in specs[:_SHOWN]:
        plus = spec.startswith("+")
        src, _, dst = spec.lstrip("+").partition(":")
        dst = dst or src
        if not src and dst:
            continue                            # a delete: refused by the deny list
        src = src or "HEAD"
        branch = current if dst == "HEAD" else dst
        if branch.startswith("refs/heads/"):
            branch = branch[len("refs/heads/"):]
        tracking = f"refs/remotes/{remote}/{branch}"
        known = _git_run(cwd, "rev-parse", "--verify", "-q", tracking)
        if force or plus:
            lost = (_git_run(cwd, "rev-list", "--count", f"{src}..{tracking}")
                    if known else None)
            alt = _fact(NEW_BRANCH, branch=branch, remote=remote)
            if lost is not None and lost.isdigit() and int(lost) == 0:
                pass                            # nothing known is lost
            elif lost is not None and lost.isdigit():
                acc.add(IRREVERSIBLE, FORCE, alternative=alt, branch=branch,
                        remote=remote, where=where, n=lost)
                continue
            else:
                acc.add(IRREVERSIBLE, FORCE_UNKNOWN, alternative=alt,
                        branch=branch, remote=remote, where=where)
                continue
        if known:
            n = _git_run(cwd, "rev-list", "--count", f"{tracking}..{src}") or "?"
            acc.add(IRREVERSIBLE, PUSH, branch=branch, remote=remote,
                    where=where, n=n)
        else:
            acc.add(IRREVERSIBLE, PUSH_NEW, branch=branch, remote=remote,
                    where=where)


def _discard(acc: _Acc, sub: str, words: list[str], cwd: Path) -> None:
    if sub == "checkout":
        if "--" in words:
            paths = words[words.index("--") + 1:]
        elif words == ["."]:
            paths = ["."]
        else:
            return
    else:
        if "--staged" in words and "--worktree" not in words and "-W" not in words:
            return                              # only unstages
        paths = _operands(words, {"-s", "--source"})
    if not paths:
        return
    changed = _git_run(cwd, "diff", "--name-only", "--", *paths)
    files = [f for f in (changed or "").splitlines() if f.strip()]
    if files:
        acc.add(IRREVERSIBLE, DISCARD, n=len(files), files=_names(files),
                alternative=_fact(SAVE_PATCH))


def _clean_untracked(acc: _Acc, words: list[str], cwd: Path) -> None:
    forced = any(w == "--force" or (w.startswith("-") and not w.startswith("--")
                                    and "f" in w) for w in words)
    if not forced:
        return
    dry = []
    for w in words:
        if w == "--force":
            continue
        if w.startswith("-") and not w.startswith("--"):
            w = w.replace("f", "")
            if w == "-":
                continue
        dry.append(w)
    listed = _git_run(cwd, "clean", "-n", *dry) or ""
    files = [ln.split(" ", 2)[-1] for ln in listed.splitlines()
             if ln.startswith("Would remove ")]
    if files:
        acc.add(IRREVERSIBLE, CLEAN, n=len(files), files=_names(files))


def _gh(acc: _Acc, words: list[str], cwd: Path) -> None:
    if len(words) < 2 and not (words and words[0] == "api"):
        return
    head = words[:2]
    if head == ["pr", "create"]:
        branch = _current_branch(cwd)
        for i, w in enumerate(words):
            if w in ("--head", "-H") and i + 1 < len(words):
                branch = words[i + 1]
        acc.add(IRREVERSIBLE, PR_CREATE, branch=branch)
    elif head == ["pr", "merge"]:
        ops = _operands(words[2:], {"-t", "--subject", "-b", "--body",
                                    "--match-head-commit", "-A", "--author-email"})
        acc.add(IRREVERSIBLE, PR_MERGE, pr=ops[0] if ops else _current_branch(cwd))
    elif head == ["release", "create"]:
        ops = _operands(words[2:], {"-t", "--title", "-n", "--notes",
                                    "-F", "--notes-file", "--target"})
        acc.add(IRREVERSIBLE, RELEASE, tag=ops[0] if ops else "?")
    elif head == ["repo", "delete"]:
        ops = _operands(words[2:], set())
        acc.add(IRREVERSIBLE, REPO_DELETE, repo=ops[0] if ops else "?")
    elif words[0] == "api":
        rest = words[1:]
        method = ""
        for i, w in enumerate(rest):
            if w in ("-X", "--method") and i + 1 < len(rest):
                method = rest[i + 1].upper()
            elif w.startswith("--method="):
                method = w.split("=", 1)[1].upper()
            elif w.startswith("-X") and len(w) > 2:
                method = w[2:].upper()
        if not method and any(w in ("-f", "-F", "--field", "--raw-field",
                                    "--input") for w in rest):
            method = "POST"
        if method in ("POST", "PUT", "PATCH", "DELETE"):
            ops = _operands(rest, {"-X", "--method", "-f", "-F", "--field",
                                   "--raw-field", "--input", "-H", "--header",
                                   "-q", "--jq", "-t", "--template",
                                   "--hostname", "--cache"})
            acc.add(IRREVERSIBLE, API_WRITE, method=method,
                    path=ops[0] if ops else "?")


_REMOTE_OPERAND = re.compile(r"^(?:[^@/:\s]+@)?([A-Za-z0-9][A-Za-z0-9._-]*):")


def _remote_copy(acc: _Acc, prog: str, words: list[str], cwd: Path) -> None:
    takes = {"-P", "-i", "-o", "-F", "-c", "-l", "-S", "-J", "-e", "--rsh",
             "-p", "--port", "--exclude", "--include", "--filter", "-f"}
    ops = _operands(words, takes)
    if prog == "sftp":
        host = next((m.group(1) for m in map(_REMOTE_OPERAND.match, ops) if m),
                    ops[0].split("@")[-1] if ops else "?")
        acc.add(IRREVERSIBLE, UPLOAD_DATA, host=host)
        return
    if len(ops) < 2:
        return
    dest = ops[-1]
    m = _REMOTE_OPERAND.match(dest) if not dest.startswith(("/", ".")) else None
    if dest.startswith("rsync://"):
        host = _where(dest).split("/", 1)[0]
    elif m:
        host = m.group(1)
    else:
        return                                  # a download: lands here
    local = [cwd / Path(s).expanduser() for s in ops[:-1]
             if not _REMOTE_OPERAND.match(s)]
    names, size = _files_and_size(local)
    if names:
        acc.add(IRREVERSIBLE, UPLOAD, n=len(names), size=_size_text(size),
                files=_names(names), host=host)
    else:
        acc.add(IRREVERSIBLE, UPLOAD_DATA, host=host)


def _http(acc: _Acc, prog: str, words: list[str], cwd: Path) -> None:
    url, method, files, inline = "", "", [], False
    i = 0
    while i < len(words):
        w = words[i]
        nxt = words[i + 1] if i + 1 < len(words) else ""
        if w in ("-T", "--upload-file"):
            files.append(nxt)
            method = method or "PUT"
            i += 2
            continue
        if w in ("-F", "--form", "-d", "--data", "--data-binary", "--data-raw",
                 "--data-urlencode", "--json", "--post-data", "--post-file",
                 "--body-data", "--body-file"):
            value = nxt.split("=", 1)[-1] if w in ("-F", "--form") else nxt
            if w in ("--post-file", "--body-file"):
                files.append(value)
            elif value.startswith(("@", "<")):
                files.append(value[1:])
            else:
                inline = True
            method = method or "POST"
            i += 2
            continue
        if w.startswith(("--post-data=", "--post-file=", "--body-data=",
                         "--body-file=", "--data=", "--json=")):
            key, value = w.split("=", 1)
            if key in ("--post-file", "--body-file") or value.startswith("@"):
                files.append(value.lstrip("@"))
            else:
                inline = True
            method = method or "POST"
        elif w in ("-X", "--request", "--method"):
            method = nxt.upper()
            i += 2
            continue
        elif w.startswith("-X") and len(w) > 2:
            method = w[2:].upper()
        elif w == "--url":
            url = nxt
            i += 2
            continue
        elif re.match(r"^[a-z][a-z0-9+.-]*://", w, re.I) and not url:
            url = w
        i += 1
    if method not in ("POST", "PUT", "PATCH", "DELETE"):
        return
    host = _where(url).split("/", 1)[0] if url else "?"
    local = [cwd / Path(f).expanduser() for f in files if f and f != "-"]
    names, size = _files_and_size(local)
    if names:
        acc.add(IRREVERSIBLE, UPLOAD, n=len(names), size=_size_text(size),
                files=_names(names), host=host)
    elif inline or files:
        acc.add(IRREVERSIBLE, UPLOAD_DATA, host=host)
    else:
        acc.add(IRREVERSIBLE, HTTP_WRITE, method=method, host=host)


def _kill(acc: _Acc, prog: str, words: list[str]) -> None:
    if prog == "kill":
        if any(w in ("-l", "-L", "--list", "-0") for w in words):
            return
        if "-s" in words:
            sig = words[words.index("-s") + 1] if words.index("-s") + 1 < len(words) else ""
            if sig in ("0",):
                return
        pids = [w for w in _operands(words, {"-s", "-n"}) if w.isdigit()]
        for pid in pids[:_SHOWN]:
            try:
                name = Path(f"/proc/{pid}/comm").read_text().strip() or "?"
            except OSError:
                name = "?"
            acc.add(IRREVERSIBLE, KILL, pid=pid, name=name)
        return
    ops = _operands(words, {"-s", "--signal", "-u", "--user", "-g", "-G",
                            "-P", "-t", "-U", "-n"})
    if ops:
        acc.add(IRREVERSIBLE, KILL_MATCH, pattern=ops[-1])


# --------------------------------------------------------------------------
# The lab: a hazardous write, as the kernel files it for a person
# --------------------------------------------------------------------------

def lab_consequence(device: str, channel: str, value: Any,
                    unit: Optional[str] = None) -> Consequence:
    """A hazardous lab write (delfin/lab manifest level ``hazardous``)."""
    shown = f"{value} {unit}".strip() if unit else str(value)
    return Consequence(IRREVERSIBLE, (_fact(LAB, device=device,
                                            channel=channel, value=shown),))


__all__ = [
    "REVERSIBLE", "RECOVERABLE", "IRREVERSIBLE", "Consequence", "Statement",
    "classify", "statement", "with_statement", "split", "is_irreversible",
    "lab_consequence", "TEMPLATES", "HEAD_IRREVERSIBLE", "HEAD_RECOVERABLE",
]
