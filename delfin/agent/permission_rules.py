"""Structured permission rules: ``permissions.deny`` / ``ask`` / ``allow``.

Input: the ``permissions`` block of the user's settings file
(``~/.delfin_settings.json``, written only by the user) and of the
project file ``<workspace>/.delfin/settings.json``::

    "permissions": {"deny": [...], "ask": [...], "allow": [...]}

and the administrator's managed policy (``managed_policy``), whose
``deny``/``ask`` are parsed by ``parse_rule`` here, marked ``managed``,
and placed ahead of the user's and the project's in the ``RuleSet``; its
``only_allow`` bounds every allow (``within_only_allow``). This module is
the one parser and matcher for all of them.

Rule forms:

* ``bash(git push)`` / ``bash(pytest:*)`` -- a command prefix, as words.
* ``read(<glob>)`` / ``write(<glob>)`` -- a path glob. Absolute or ``~``
  globs are absolute; anything else is relative to the workspace. ``*``
  is one path segment, ``**`` any number.
* ``web_fetch(domain:example.org)`` -- that host and its subdomains.
* ``mcp__server__tool``, ``mcp__server__*`` or ``mcp__server`` -- MCP tools.
* a bare tool name (``web_search``) -- every call of that tool.

Output: a ``RuleSet`` and three questions asked of it -- ``friction`` (the
deny or ask rule a call meets), ``allows`` (the allow rule that answers one
soft question the gate would otherwise put to the user) and
``denies_path`` (a deny rule on a path a tool resolved itself).

Semantics, in the order they are applied by the gate (api_client):

* deny > ask > allow. Deny and ask only add friction and hold in every
  mode, ``bypassPermissions`` included; an ask with nobody to ask is a
  refusal.
* Deny and ask match generously: a bash rule matches its words anywhere in
  any segment of a chain, substitution or ``sh -c`` body (``echo git push``
  is a false positive -- friction only); read and write rules match every
  path-like token of a shell command; a glob without a wildcard covers the
  path and everything below it.
* Allow matches narrowly: a bash rule covers one simple command whose words
  start with the rule's (``:*``) or equal them -- a command with any of
  ``; & | < > ` $ ( ) { }``, a newline or a backslash is never covered; a
  path rule covers a path only when the path itself matches the glob.
* Allow never reaches a hard refusal: the gate consults it only at the
  places where it would otherwise ask (see ``api_client._rule_answers``).
* A project file may add deny and ask rules; its allow list (and any
  ``repos``/``global_allow``) is ignored and reported. A repository must
  not grant itself permissions.
* Scope. The user's file may hold ``"repos": {"<key>": {"allow": [...],
  "ask": [...], "deny": [...]}}``: an entry's rules apply only in the
  workspace whose repository it names. A key is the ``agent.git_roles``
  key (``github.com/owner/repo`` or a remote URL, compared with
  ``git_role.repo_key`` of the workspace: origin's push URL, else the
  only remote's, else the top-level path), or an absolute path / ``~``
  path glob matched against the checkout's top level and the workspace.
  The top-level ``deny``/``ask`` apply everywhere. The top-level
  ``allow`` applies everywhere too while ``global_allow`` is unset, with
  a start notice naming each rule (``scope_notice``); ``true`` keeps it
  without the notice, ``false`` (or a value that is not a bool) turns it
  off. Rejected: turning an unscoped allow off outright (a working setup
  would start asking with nothing said), and scoping it to the first
  repository it answered in (state the user never wrote).
* A rule that cannot be parsed grants nothing. A broken deny or ask rule
  is an error in ``RuleSet.errors``, shown at startup and by
  ``/permissions``.

Rejected: regexes as the rule language (the classifier's free-text rules
and ``bash_auto_allow_patterns`` already exist and are both easy to write
too wide); evaluating rules by filtering the advertised tool list
(advertising is not enforcement); reading the file through
``user_settings.load_settings``, which merges and WRITES defaults back --
this module only reads.

Cost: two ``stat`` calls per gated tool call; the files are parsed again
only when their size or mtime changes.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

LISTS = ("deny", "ask", "allow")
#: Keys of the user's block that scope the allow list (not lists).
SCOPE_KEYS = ("repos", "global_allow")
PROJECT_RELPATH = Path(".delfin") / "settings.json"

_KIND_ALIASES = {
    "bash": "bash", "shell": "bash",
    "read": "read",
    "write": "write", "edit": "write",
    "web_fetch": "web_fetch", "webfetch": "web_fetch",
}
_RULE_RE = re.compile(r"^([A-Za-z_]+)\((.*)\)$", re.DOTALL)
_TOOL_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_\-]*$")
_MCP_RE = re.compile(r"^mcp__[A-Za-z0-9_\-]+(?:__(?:[A-Za-z0-9_\-]+|\*))?$")
_HOST_RE = re.compile(r"^[a-z0-9](?:[a-z0-9\-]*[a-z0-9])?(?:\.[a-z0-9](?:[a-z0-9\-]*[a-z0-9])?)*$")
#: A command holding any of these is never covered by an allow rule.
_NOT_SIMPLE = set(";&|<>`$(){}\n\r\\")
_BASH_TOOLS = frozenset({"bash", "bash_background"})
#: Tools whose path arguments are write targets.
WRITE_TOOLS = frozenset({
    "write_file", "edit_file", "multi_edit", "apply_patch", "notebook_edit",
    "make_plot", "edit_sheet", "fill_pdf_form", "fill_docx_template",
    "create_docx", "fill_series", "merge_pdfs", "split_pdf", "create_pdf",
    "draft_email",
})
_PATH_KEYS = ("path", "file_path", "filename", "file", "target", "notebook_path",
              "directory", "dir", "cwd", "root", "source", "src", "dest",
              "destination", "output_path", "input_path", "pdf_path",
              "template_path", "output", "paths", "files", "sources")


@dataclass(frozen=True)
class Rule:
    """One parsed rule. ``value``: words (bash), a glob (read/write), a host
    (web_fetch) or a tool name (tool)."""
    text: str
    kind: str
    value: Any
    wildcard: bool = False
    list: str = ""
    source: str = ""
    #: True for a rule of the administrator's managed policy.
    managed: bool = False
    #: The ``permissions.repos`` key a user rule was written under, as
    #: written; "" for a rule of the top-level lists.
    scope: str = ""


@dataclass
class RuleSet:
    #: Managed deny/ask rules first, then the user's, then the project's.
    deny: tuple = ()
    ask: tuple = ()
    allow: tuple = ()
    #: (source path, list, rule text, reason) for every entry not in force.
    errors: list = field(default_factory=list)
    user_path: str = ""
    project_path: str = ""
    #: The managed ``only_allow`` lists, one tuple of rules per policy file
    #: that sets one; a grant must fit every tuple. () = no bound.
    only_allow: tuple = ()
    #: The ``managed_policy.Policy`` merged in, or None.
    policy: Any = None
    #: The workspace's repository key (``git_role.repo_key``: normalised
    #: remote URL, else top-level path); "" outside a git checkout or
    #: when it cannot be told.
    repo_key: str = ""
    #: The resolved top-level path of the checkout ("" outside one).
    repo_path: str = ""
    #: ``permissions.global_allow`` as written: None (unset: the global
    #: allow list applies, with a start notice), True, or False.
    global_allow: Any = None
    #: Global allow rules not in force because ``global_allow`` is off.
    allow_off: tuple = ()
    #: User rules of ``permissions.repos`` entries that do not match this
    #: workspace, for the ``/permissions`` view.
    other_repos: tuple = ()

    def is_empty(self) -> bool:
        return not (self.deny or self.ask or self.allow)


def parse_rule(text: Any) -> "Rule | str":
    """A ``Rule``, or the reason ``text`` is not one."""
    if not isinstance(text, str):
        return "a rule is a string"
    raw = text.strip()
    if not raw:
        return "empty rule"
    m = _RULE_RE.match(raw)
    if m:
        kind = _KIND_ALIASES.get(m.group(1).lower())
        body = m.group(2).strip()
        if kind is None:
            return (f"unknown rule kind {m.group(1)!r} (bash, read, write, "
                    "web_fetch, or a bare tool name)")
        if not body:
            return f"{kind}() needs a value"
        if kind == "bash":
            wildcard = body.endswith(":*")
            if wildcard:
                body = body[:-2].rstrip()
            if not body or any(c in _NOT_SIMPLE for c in body):
                return "a bash rule names the words of one command"
            try:
                words = tuple(shlex.split(body))
            except ValueError as exc:
                return f"bash rule does not parse: {exc}"
            if not words:
                return "a bash rule names the words of one command"
            return Rule(raw, "bash", words, wildcard)
        if kind in ("read", "write"):
            if "\0" in body:
                return "a path glob holds no NUL"
            return Rule(raw, kind, body)
        # web_fetch
        if not body.lower().startswith("domain:"):
            return "web_fetch rules are web_fetch(domain:<host>)"
        host = body[len("domain:"):].strip().lower().rstrip(".")
        # ``*.example.org``: the subdomains only, not the domain itself.
        sub_only = host.startswith("*.")
        if sub_only:
            host = host[2:]
        if not host.isascii() or not _HOST_RE.match(host):
            return f"{host!r} is not a host name"
        return Rule(raw, "web_fetch", host, sub_only)
    if raw.startswith("mcp__"):
        if not _MCP_RE.match(raw):
            return "an MCP rule is mcp__<server>, mcp__<server>__<tool> or mcp__<server>__*"
        return Rule(raw, "tool", raw)
    if _TOOL_RE.match(raw):
        return Rule(raw, "tool", raw)
    return "not a rule (bash(...), read(...), write(...), web_fetch(domain:...), mcp__..., or a tool name)"


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

_CACHE: dict = {}
_CACHE_LOCK = threading.Lock()


def user_settings_path() -> Path:
    from delfin import user_settings
    return Path(user_settings.get_settings_path())


def project_settings_path(workspace: Any) -> Optional[Path]:
    try:
        return Path(workspace).expanduser().resolve() / PROJECT_RELPATH
    except Exception:
        return None


def _stamp(path: Optional[Path]) -> tuple:
    if path is None:
        return ("", None)
    try:
        st = path.stat()
        return (str(path), (st.st_mtime_ns, st.st_size, st.st_ino))
    except OSError:
        return (str(path), None)


def _read_block(path: Optional[Path], errors: list) -> dict:
    """The ``permissions`` block of ``path``; {} when there is none."""
    if path is None:
        return {}
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return {}
    except OSError as exc:
        errors.append((str(path), "", "", f"cannot be read: {exc}"))
        return {}
    try:
        data = json.loads(text)
    except ValueError as exc:
        errors.append((str(path), "", "",
                       f"is not valid JSON ({exc}); no rule from it is in force"))
        return {}
    if not isinstance(data, dict):
        return {}
    block = data.get("permissions")
    if block is None:
        return {}
    if not isinstance(block, dict):
        errors.append((str(path), "", "",
                       "'permissions' is not an object; no rule from it is in force"))
        return {}
    return block


def _parse_lists(block: dict, source: str, *, project: bool,
                 errors: list, scope: str = "") -> dict:
    out = {k: [] for k in LISTS}
    for key, entries in block.items():
        if key in SCOPE_KEYS and not scope:
            if project:
                errors.append((source, key, "",
                               "ignored: a project file cannot grant "
                               "permissions, only add deny and ask rules"))
            continue
        if key not in LISTS:
            errors.append((source, str(key), "",
                           f"unknown list {key!r} (deny, ask, allow"
                           + ("" if project or scope else
                              "; repos, global_allow") + ")"))
            continue
        if entries is None:
            continue
        if not isinstance(entries, list):
            errors.append((source, key, "", f"'{key}' is not a list"))
            continue
        for entry in entries:
            if project and key == "allow":
                errors.append((source, key, str(entry),
                               "ignored: a project file cannot grant "
                               "permissions, only add deny and ask rules"))
                continue
            rule = parse_rule(entry)
            if isinstance(rule, str):
                errors.append((source, key, str(entry), rule))
                continue
            out[key].append(Rule(rule.text, rule.kind, rule.value,
                                 rule.wildcard, key, source, scope=scope))
    return out


# ---------------------------------------------------------------------------
# The repository a workspace belongs to (permissions.repos)
# ---------------------------------------------------------------------------

_REPO_CACHE: dict = {}


def _git_marker(start: Path) -> tuple:
    """A stamp of the nearest ``.git`` at or above ``start`` (for a
    directory its inode and ``config``'s stamp; for a worktree's ``.git``
    file its own); () outside a checkout. A remote changed in that config
    or a new ``git init`` changes it; a remote changed in a linked
    worktree's common config is seen on the next process start."""
    for d in (start, *start.parents):
        marker = d / ".git"
        try:
            st = marker.stat()
        except OSError:
            continue
        if marker.is_dir():
            # Not the directory's mtime: every git command touches it, and
            # the remote lives in config.
            return (str(marker), st.st_ino) + _stamp(marker / "config")
        return (str(marker), st.st_mtime_ns, st.st_ino)
    return ()


def repo_identity(workspace: Any) -> tuple:
    """(repo key, top-level path, resolved workspace) for ``workspace``.

    The key is ``git_role.repo_key`` -- the key ``agent.git_roles`` uses:
    origin's (else the only remote's) push URL normalised, else the
    top-level path. ("", "", ws) outside a checkout or when git cannot
    tell. Cached per workspace and ``.git`` stamp: the git subprocesses
    run once per checkout and process, not per gated call."""
    if workspace is None:
        return ("", "", "")
    try:
        ws = Path(workspace).expanduser().resolve()
    except Exception:
        return ("", "", "")
    marker = _git_marker(ws)
    if not marker:
        return ("", "", str(ws))
    key = (str(ws), marker)
    with _CACHE_LOCK:
        hit = _REPO_CACHE.get(key)
    if hit is not None:
        return hit
    try:
        from . import git_role
        top = git_role.toplevel(ws) or ""
        rkey = (git_role.repo_key(ws) or "") if top else ""
    except Exception:
        top, rkey = "", ""
    out = (rkey, top, str(ws))
    with _CACHE_LOCK:
        if len(_REPO_CACHE) > 64:
            _REPO_CACHE.clear()
        _REPO_CACHE[key] = out
    return out


def repo_entry_matches(entry: str, identity: tuple) -> "bool | None":
    """True when the ``permissions.repos`` key ``entry`` names the
    repository of ``identity`` (``repo_identity``); None when ``entry`` is
    not a key.

    A key starting with ``/`` or ``~`` is a path: with ``*``/``?`` a glob
    over the top-level path or the workspace, otherwise that path and
    everything below it -- the checkout's location, whatever its remote.
    Any other key is a repository key as ``agent.git_roles`` takes it
    (``github.com/owner/repo`` or a remote URL) and must equal the
    workspace's ``repo_key``."""
    text = str(entry or "").strip()
    rkey, top, ws = identity
    if text.startswith(("/", "~")):
        pattern = os.path.expanduser(text)
        if not pattern.startswith("/"):
            return None
        wild = any(ch in pattern for ch in "*?")
        return any(path_matches(pattern, p, None, subtree=not wild)
                   for p in (top, ws) if p)
    try:
        from . import git_role
        nk = git_role.entry_key(text)
    except Exception:
        nk = None
    if not nk:
        return None
    return bool(rkey) and nk == rkey


def _repo_lists(block: dict, source: str, identity: tuple,
                errors: list) -> tuple:
    """(lists matching this workspace, rules of the other entries) from
    the user's ``permissions.repos``."""
    here = {k: [] for k in LISTS}
    other: list = []
    repos = block.get("repos")
    if repos is None:
        return here, other
    if not isinstance(repos, dict):
        errors.append((source, "repos", "",
                       "'repos' is not an object; no rule from it is in force"))
        return here, other
    for entry, sub in repos.items():
        where = f"{source} repos[{entry}]"
        match = repo_entry_matches(entry, identity)
        if match is None:
            errors.append((where, "repos", "",
                           "not a repository key (github.com/owner/repo, a "
                           "remote URL, or an absolute path or path glob)"))
            continue
        if not isinstance(sub, dict):
            errors.append((where, "repos", "",
                           "an entry is {\"allow\": [...], \"ask\": [...], "
                           "\"deny\": [...]}; no rule from it is in force"))
            continue
        got = _parse_lists(sub, source, project=False, errors=errors,
                           scope=str(entry))
        for k in LISTS:
            if match:
                here[k].extend(got[k])
            else:
                other.extend(got[k])
    return here, other


def _managed_policy():
    """The managed policy in force, or None when there is none."""
    try:
        from . import managed_policy as _mp
        pol = _mp.load()
    except Exception:
        return None
    return pol if (pol.active or pol.errors) else None


def load(workspace: Any = None) -> RuleSet:
    """The rules in force for ``workspace``: the administrator's managed
    deny and ask rules ahead of the user's and the project's, and the
    managed ``only_allow`` bound on every allow. Never raises."""
    local = _load_local(workspace)
    pol = _managed_policy()
    if pol is None:
        return local
    key = ("combined", id(local), id(pol))
    with _CACHE_LOCK:
        hit = _CACHE.get(key)
    # The objects are kept in the entry: an id is only a key while the
    # object it named is alive.
    if hit is not None and hit[0] is local and hit[1] is pol:
        return hit[2]
    # The policy's own errors stay in ``policy.errors`` (shown by
    # ``managed_policy.summary_lines``), not in ``errors``: one report each.
    rs = RuleSet(tuple(pol.deny) + local.deny, tuple(pol.ask) + local.ask,
                 local.allow, list(local.errors), local.user_path,
                 local.project_path,
                 tuple(pol.only_allow), pol,
                 repo_key=local.repo_key, repo_path=local.repo_path,
                 global_allow=local.global_allow, allow_off=local.allow_off,
                 other_repos=local.other_repos)
    with _CACHE_LOCK:
        if len(_CACHE) > 64:
            _CACHE.clear()
        _CACHE[key] = (local, pol, rs)
    return rs


def _load_local(workspace: Any = None) -> RuleSet:
    """The user's and the project's rules for ``workspace``. Never raises:
    a failure is an error entry and grants nothing."""
    try:
        upath = user_settings_path()
    except Exception:
        upath = None
    ppath = project_settings_path(workspace) if workspace is not None else None
    identity = repo_identity(workspace)
    key = (_stamp(upath), _stamp(ppath), identity)
    with _CACHE_LOCK:
        hit = _CACHE.get(key)
    if hit is not None:
        return hit
    errors: list = []
    lists = {k: [] for k in LISTS}
    gallow, allow_off, other = None, [], []
    for path, project in ((upath, False), (ppath, True)):
        block = _read_block(path, errors)
        if not block:
            continue
        got = _parse_lists(block, str(path), project=project, errors=errors)
        if project:
            for k in LISTS:
                lists[k].extend(got[k])
            continue
        gallow = block.get("global_allow")
        if gallow is not None and not isinstance(gallow, bool):
            errors.append((str(path), "global_allow", repr(gallow),
                           "is neither true nor false; the global allow "
                           "list is off until it is"))
            gallow = False
        here, other = _repo_lists(block, str(path), identity, errors)
        # The repository's own rules after the global deny/ask (every
        # deny before any ask either way) and before the global allows,
        # so the scoped rule is the one reported.
        lists["deny"].extend(got["deny"] + here["deny"])
        lists["ask"].extend(got["ask"] + here["ask"])
        lists["allow"].extend(here["allow"])
        if gallow is False:
            allow_off = got["allow"]
        else:
            lists["allow"].extend(got["allow"])
    rs = RuleSet(tuple(lists["deny"]), tuple(lists["ask"]),
                 tuple(lists["allow"]), errors,
                 str(upath or ""), str(ppath or ""),
                 repo_key=identity[0], repo_path=identity[1],
                 global_allow=gallow, allow_off=tuple(allow_off),
                 other_repos=tuple(other))
    with _CACHE_LOCK:
        if len(_CACHE) > 64:
            _CACHE.clear()
        _CACHE[key] = rs
    return rs


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def bare_tool_name(name: str) -> str:
    text = str(name or "")
    return text.rsplit("__", 1)[-1] if text.startswith("mcp__") else text


def _tool_matches(rule: Rule, name: str, *, strict: bool) -> bool:
    value = str(rule.value)
    if value.startswith("mcp__"):
        if value.endswith("__*"):
            return name.startswith(value[:-1])
        if value.count("__") == 1:          # mcp__server
            return name.startswith(value + "__")
        return name == value
    if name == value:
        return True
    if strict:
        # A bare allow names the native tool, or that tool in the native
        # coding namespace; never a same-named tool of another server.
        return name == f"mcp__kit-coding__{value}"
    return bare_tool_name(name) == value


def _glob_regex(glob: str) -> str:
    out, i = [], 0
    while i < len(glob):
        c = glob[i]
        if glob.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
            continue
        if glob.startswith("**", i):
            out.append(".*")
            i += 2
            continue
        if c == "*":
            out.append("[^/]*")
        elif c == "?":
            out.append("[^/]")
        else:
            out.append(re.escape(c))
        i += 1
    return "".join(out)


def _anchor(glob: str, workspace: Any) -> str:
    """``glob`` as an absolute pattern, its literal head resolved."""
    g = glob.strip()
    if g.startswith("~"):
        g = os.path.expanduser(g)
    if not g.startswith("/"):
        base = Path(workspace or ".").expanduser()
        try:
            base = base.resolve()
        except Exception:
            pass
        g = str(base).rstrip("/") + "/" + (g[2:] if g.startswith("./") else g)
    parts = g.split("/")
    head = []
    for p in parts:
        if any(ch in p for ch in "*?"):
            break
        head.append(p)
    rest = parts[len(head):]
    literal = "/".join(head) or "/"
    try:
        literal = str(Path(literal).resolve(strict=False))
    except Exception:
        pass
    return literal if not rest else literal.rstrip("/") + "/" + "/".join(rest)


def path_matches(glob: str, path: Any, workspace: Any, *,
                 subtree: bool) -> bool:
    """True when ``path`` (absolute) matches ``glob``. ``subtree``: a glob
    without a wildcard also covers everything below it (deny/ask)."""
    try:
        target = str(path)
        pattern = _anchor(glob, workspace)
    except Exception:
        return False
    if not any(ch in pattern for ch in "*?"):
        pattern = pattern.rstrip("/") or "/"
        if target == pattern:
            return True
        return subtree and target.startswith(pattern.rstrip("/") + "/")
    return re.fullmatch(_glob_regex(pattern), target) is not None


def _host_matches(rule_host: str, host: str, sub_only: bool = False) -> bool:
    host = (host or "").lower().rstrip(".")
    if sub_only:
        return bool(host) and host.endswith("." + rule_host)
    return bool(host) and (host == rule_host or host.endswith("." + rule_host))


def url_host(url: Any) -> str:
    text = str(url or "").strip()
    m = re.match(r"^[A-Za-z][A-Za-z0-9+.\-]*://(?:[^@/\s]*@)?(\[[^\]]*\]|[^:/\s?#]+)", text)
    return m.group(1).lower() if m else ""


def _words(segment: str) -> list:
    try:
        return shlex.split(segment, comments=False, posix=True)
    except ValueError:
        return [w.strip("'\"") for w in segment.split()]


def _segments(cmd: str, depth: int = 0) -> list:
    """Word lists of every command a shell line may run, over-split on
    purpose: operators, substitutions and grouping split it, and a word
    with whitespace (an ``sh -c`` body, an ``eval`` string) is split again."""
    out = []
    pieces = re.split(r"\|\||&&|\$\(|<\(|>\(|[;|&\n\r()`{}]", str(cmd or ""))
    for piece in pieces:
        words = _words(piece)
        if not words:
            continue
        out.append(words)
        if depth < 3:
            for w in words:
                if any(ch.isspace() for ch in w):
                    out.extend(_segments(w, depth + 1))
    return out


def _bash_rule_hits(rule: Rule, cmd: str) -> bool:
    """Deny/ask: the rule's first word (by basename) anywhere in a segment,
    its other words following in order."""
    head, rest = rule.value[0], list(rule.value[1:])
    head_base = head.rsplit("/", 1)[-1]
    for words in _segments(cmd):
        for i, w in enumerate(words):
            w = w.lstrip("\\")
            if w != head and w.rsplit("/", 1)[-1] != head_base:
                continue
            want = list(rest)
            for later in words[i + 1:]:
                if want and later == want[0]:
                    want.pop(0)
            if not want:
                return True
    return False


def simple_command_words(cmd: str) -> Optional[list]:
    """The words of ``cmd`` when it is one simple command; else None."""
    text = str(cmd or "").strip()
    if not text or any(c in _NOT_SIMPLE for c in text):
        return None
    try:
        words = shlex.split(text, posix=True)
    except ValueError:
        return None
    if not words or "=" in words[0]:
        return None
    return words


#: Commands that run another command named in their words. An allow rule
#: for one of them would cover whatever program follows.
_WRAPPERS = frozenset({
    "env", "xargs", "timeout", "nice", "nohup", "time", "command", "builtin",
    "exec", "eval", "sudo", "doas", "su", "stdbuf", "chroot", "setsid",
    "flock", "watch", "parallel", "find", "ssh", "script", "unbuffer",
    "taskset", "ionice", "strace", "ltrace", "gdb", "busybox",
})
#: Interpreters, and the flags with which their words carry program text.
_INTERPRETER_RE = re.compile(
    r"^(?:python[0-9.]*|py|pypy[0-9.]*|sh|bash|zsh|dash|ksh|fish|csh|tcsh|"
    r"perl[0-9.]*|ruby|node|nodejs|deno|bun|Rscript|R|julia|php|lua[0-9.]*|"
    r"tclsh|wish|awk|gawk|mawk|osascript|pwsh|powershell)$")
_INLINE_FLAGS = frozenset({"-c", "-e", "-E", "--eval", "--command", "-r",
                           "--exec", "-x", "-p", "-n", "-i"})


def carries_a_program(words: list) -> bool:
    """True when the command runs a program its words carry: a wrapper
    (``env``, ``xargs``, ``timeout`` ...), or an interpreter given inline
    code (``python3 -c``, ``sh -c``, ``perl -e``) or read from stdin."""
    head = words[0].rsplit("/", 1)[-1]
    if head in _WRAPPERS:
        return True
    if _INTERPRETER_RE.match(head):
        rest = words[1:]
        if not rest or rest[0] == "-":
            return True
        if head in ("awk", "gawk", "mawk"):
            return True
        return any(w in _INLINE_FLAGS or (w.startswith("-") and not
                   w.startswith("--") and "c" in w[1:] and head in (
                       "sh", "bash", "zsh", "dash", "ksh"))
                   for w in rest)
    return False


def _bash_allow_hits(rule: Rule, cmd: str) -> bool:
    words = simple_command_words(cmd)
    if words is None or carries_a_program(words):
        return False
    want = list(rule.value)
    if rule.wildcard:
        return words[:len(want)] == want
    return words == want


_ABS_TOKEN_RE = re.compile(r"(?<![A-Za-z0-9_])(?:~|\$HOME|/)[^\s;|&<>()`'\"]+")
_REL_TOKEN_RE = re.compile(r"(?<![^\s;|&<>()=])(?![-~$/])[^\s;|&<>()`'\"=]+")


def bash_path_tokens(cmd: str, workspace: Any, cwd: Any = None) -> list:
    """Absolute paths for every path-like token of a shell command.
    Generous on purpose: used for deny and ask only."""
    cleaned = re.sub(r"['\"]", " ", str(cmd or ""))
    base = Path(workspace or ".")
    if cwd:
        c = Path(str(cwd)).expanduser()
        base = c if c.is_absolute() else base / c
    out = []
    for tok in _ABS_TOKEN_RE.findall(cleaned):
        out.append(Path(os.path.expanduser(tok.replace("$HOME", "~"))))
    for tok in _REL_TOKEN_RE.findall(cleaned):
        if "://" in tok or tok.startswith("-"):
            continue
        out.append(base / tok)
    return out


def _arg_paths(args: dict, workspace: Any) -> list:
    out = []
    base = Path(workspace or ".")
    for key in _PATH_KEYS:
        v = (args or {}).get(key)
        values = v if isinstance(v, (list, tuple)) else [v]
        for item in values:
            if isinstance(item, str) and item.strip() and "\0" not in item:
                p = Path(os.path.expanduser(item.strip()))
                out.append(p if p.is_absolute() else base / p)
    return out


def _patch_paths(args: dict, workspace: Any) -> list:
    diff = str((args or {}).get("diff") or "")
    out = []
    for m in re.finditer(r"^(?:\+\+\+|---) (?:[ab]/)?(\S+)", diff, re.MULTILINE):
        if m.group(1) != "/dev/null":
            out.append(Path(workspace or ".") / m.group(1))
    return out


def _forms(paths: Iterable[Path]) -> list:
    """Each path as written (absolute) and resolved: a deny matches either."""
    out = []
    for p in paths:
        s = os.path.normpath(str(p))
        out.append(s)
        try:
            r = str(Path(s).resolve(strict=False))
            if r != s:
                out.append(r)
        except Exception:
            pass
    return out


def _subjects(name: str, args: dict, workspace: Any) -> dict:
    bare = bare_tool_name(name)
    args = args if isinstance(args, dict) else {}
    read, write, hosts, cmd = [], [], [], None
    if bare in _BASH_TOOLS:
        cmd = str(args.get("command") or "")
        toks = bash_path_tokens(cmd, workspace, args.get("cwd"))
        read.extend(toks)
        write.extend(toks)
        hosts.extend(h for h in (url_host(t) for t in re.findall(r"\S+://\S+", cmd)) if h)
    else:
        paths = _arg_paths(args, workspace)
        if bare == "apply_patch":
            paths += _patch_paths(args, workspace)
        read.extend(paths)
        if bare in WRITE_TOOLS:
            write.extend(paths)
    if bare == "web_fetch":
        h = url_host(args.get("url"))
        if h:
            hosts.append(h)
    if bare == "verify_source":
        # SRC1 (operator decision 2026-11-01): the network policy binds
        # this tool like web_fetch. Its target host is NOT user-supplied
        # -- it is the fixed API host of the action, derived in ONE place
        # (delfin.agent.source_apis.host_for_action), so a new action
        # cannot silently bypass the rule layer. An unknown action gives
        # no host, so nothing can talk the rule layer into one.
        h = _source_api_host(args.get("action"))
        if h:
            hosts.append(h)
    return {"bash": cmd, "read": _forms(read), "write": _forms(write),
            "hosts": hosts}


def _hits(rule: Rule, name: str, subj: dict, workspace: Any) -> bool:
    if rule.kind == "tool":
        return _tool_matches(rule, name, strict=False)
    if rule.kind == "bash":
        return subj["bash"] is not None and _bash_rule_hits(rule, subj["bash"])
    if rule.kind in ("read", "write"):
        pool = subj["read"] if rule.kind == "read" else subj["write"]
        return any(path_matches(rule.value, p, workspace, subtree=True)
                   for p in pool)
    if rule.kind == "web_fetch":
        return any(_host_matches(rule.value, h) for h in subj["hosts"])
    return False


def _source_api_host(action) -> Optional[str]:
    """The fixed API host a verify_source action dials (or None).

    Reads the mapping from delfin.agent.source_apis -- one derivation
    place; the rule layer never keeps its own host list that could
    drift from what the tool actually dials.
    """
    try:
        from .source_apis import host_for_action
        return host_for_action(action)
    except Exception:
        # source_apis unavailable -> no host -> the rule layer treats
        # the call as having no web subject; the tool itself will fail
        # its own import and the call goes nowhere.
        return None


def _ceiling_rule(rs: RuleSet, name: str, args: dict) -> Optional[Rule]:
    """The managed network ceiling as a deny rule, when ``web_fetch``
    (or a verify_source action, SRC1 2026-11-01) names a host outside
    it (``managed_policy.domain_within_max``)."""
    pol = rs.policy
    bare = bare_tool_name(name)
    if (pol is None or not getattr(pol, "allowed_domains_max", ())
            or bare not in ("web_fetch", "verify_source")):
        return None
    try:
        from . import managed_policy as _mp
        if bare == "verify_source":
            host = _source_api_host((args or {}).get("action"))
        else:
            host = url_host((args or {}).get("url"))
        if _mp.domain_within_max(host, pol):
            return None
        src = pol.sources[0].path if pol.sources else ""
    except Exception:
        src = ""
    return Rule("network.allowed_domains_max", "web_fetch", "", list="deny",
                source=src, managed=True)


def friction(rs: RuleSet, name: str, args: dict, workspace: Any, *,
             managed_only: bool = False) -> "tuple[str, Rule] | None":
    """("deny", rule) or ("ask", rule) for the first rule the call meets.

    Order: every deny before any ask (a deny anywhere outranks a
    question), and within each list the managed rules before the user's
    and the project's; the managed network ceiling is a deny.
    ``managed_only``: the user's and the project's rules are skipped (an
    entry point where only the managed policy is decided)."""
    args = args if isinstance(args, dict) else {}
    ceiling = _ceiling_rule(rs, name, args)
    if not (rs.deny or rs.ask or ceiling):
        return None
    try:
        subj = _subjects(name, args, workspace)
    except Exception:
        subj = {"bash": None, "read": [], "write": [], "hosts": []}
    for which, rules in (("deny", rs.deny), ("ask", rs.ask)):
        for rule in rules:
            if managed_only and not rule.managed:
                continue
            try:
                if _hits(rule, name, subj, workspace):
                    return which, rule
            except Exception:
                if rule.managed:
                    # An administrator's rule that cannot be evaluated
                    # holds: deny refuses, ask asks.
                    return which, rule
                continue
        if which == "deny" and ceiling is not None:
            return "deny", ceiling
    return None


def denies_path(rs: RuleSet, kind: str, path: Any, workspace: Any, *,
                managed_only: bool = False) -> Optional[Rule]:
    """The deny rule on ``path`` for a ``kind`` ("read"/"write") access a
    tool resolved itself. A read rule covers writes too. Managed rules
    come first (``rs.deny`` is in that order)."""
    kinds = ("read",) if kind == "read" else ("read", "write")
    forms = _forms([Path(str(path))])
    for rule in rs.deny:
        if managed_only and not rule.managed:
            continue
        if rule.kind in kinds and any(
                path_matches(rule.value, f, workspace, subtree=True)
                for f in forms):
            return rule
    return None


def allows(rs: RuleSet, site: str, name: str, args: dict, workspace: Any,
           target: Any = None) -> Optional[Rule]:
    """The allow rule answering the question at ``site`` for this call.

    ``site``: "bash" (the ordinary command dialog), "read" (a read outside
    the roots, ``target`` the resolved path), "write" (the default-mode
    write dialog, ``target`` the resolved path), "tool" (a side-effecting
    MCP tool). A call answers only the site of its own kind: an allowed
    command answers no read question it raises, and the reverse.

    A managed ``only_allow`` bounds the answer: an allow rule for a call
    outside it answers nothing (``within_only_allow``)."""
    if not rs.allow:
        return None
    rule = _allow_match(rs.allow, site, name, args, workspace, target)
    if rule is None:
        return None
    if not within_only_allow(rs, site, name, args, workspace, target):
        return None
    return rule


def _only_allow_segment_fits(rules: tuple, name: str, words: list) -> bool:
    for r in rules:
        if r.kind == "tool" and _tool_matches(r, name, strict=True):
            return True
        if r.kind == "bash" and (words[:len(r.value)] == list(r.value)
                                 if r.wildcard else words == list(r.value)):
            return True
    return False


#: Chain operators split a command into segments for ``only_allow``; any
#: other shell syntax keeps it outside every grant.
_ONLY_ALLOW_SPLIT = re.compile(r"\|\||&&|[;|&]")
_ONLY_ALLOW_UNSAFE = set("<>`$(){}\n\r\\")


def bash_within_only_allow(rs: RuleSet, cmd: str,
                           name: str = "bash") -> bool:
    """False when the managed ``only_allow`` leaves ``cmd`` outside every
    grant. True when there is no bound.

    Narrow: each segment of a ``; && || | &`` chain must start with the
    words of an ``only_allow`` bash rule (``:*``) or equal them, under
    every policy file that sets a list; a command with a substitution, a
    redirection, grouping, a newline or a backslash, or a segment that
    runs a program its words carry (``env``, ``xargs``, ``sh -c`` ...;
    ``carries_a_program``) fits none."""
    if not rs.only_allow:
        return True
    text = str(cmd or "")
    if not text.strip() or any(c in _ONLY_ALLOW_UNSAFE for c in text):
        return False
    segments = [s.strip() for s in _ONLY_ALLOW_SPLIT.split(text) if s.strip()]
    if not segments:
        return False
    parsed = []
    for seg in segments:
        try:
            words = shlex.split(seg, posix=True)
        except ValueError:
            return False
        if not words or "=" in words[0] or carries_a_program(words):
            return False
        parsed.append(words)
    return all(_only_allow_segment_fits(group, name, words)
               for group in rs.only_allow for words in parsed)


def within_only_allow(rs: RuleSet, site: str, name: str, args: dict,
                      workspace: Any, target: Any = None) -> bool:
    """True when the call fits the managed ``only_allow`` at ``site``
    (or there is no bound). ``site`` as for ``allows``; a bash call is
    judged by ``bash_within_only_allow`` whatever the site."""
    if not rs.only_allow:
        return True
    args = args if isinstance(args, dict) else {}
    if bare_tool_name(name) in _BASH_TOOLS:
        return bash_within_only_allow(rs, str(args.get("command") or ""),
                                      name)
    return all(_allow_match(group, site, name, args, workspace, target)
               is not None for group in rs.only_allow)


def _allow_match(rules: tuple, site: str, name: str, args: dict,
                 workspace: Any, target: Any = None) -> Optional[Rule]:
    """The first of ``rules`` that answers ``site`` for this call (the
    narrow allow matching), or None."""
    rs = RuleSet(allow=tuple(rules))
    bare = bare_tool_name(name)
    args = args if isinstance(args, dict) else {}
    if site == "bash":
        if bare not in _BASH_TOOLS:
            return None
        cmd = str(args.get("command") or "")
        words = simple_command_words(cmd)
        if words is None or carries_a_program(words):
            return None
        for rule in rs.allow:
            if rule.kind == "bash" and _bash_allow_hits(rule, cmd):
                return rule
            if rule.kind == "tool" and _tool_matches(rule, name, strict=True):
                return rule
        return None
    if site in ("read", "write"):
        if bare in _BASH_TOOLS:
            return None
        if site == "write" and bare not in WRITE_TOOLS:
            return None
        if site == "read" and bare in WRITE_TOOLS:
            return None
        if target is None:
            return None
        try:
            resolved = str(Path(str(target)).resolve(strict=False))
        except Exception:
            return None
        for rule in rs.allow:
            if rule.kind == site and path_matches(
                    rule.value, resolved, workspace, subtree=False):
                return rule
            if rule.kind == "tool" and _tool_matches(rule, name, strict=True):
                return rule
        return None
    if site == "tool":
        for rule in rs.allow:
            if rule.kind == "tool" and _tool_matches(rule, name, strict=True):
                return rule
            if (rule.kind == "web_fetch" and bare == "web_fetch"
                    and _host_matches(rule.value, url_host(args.get("url")),
                                       rule.wildcard)):
                return rule
        return None
    return None


# ---------------------------------------------------------------------------
# What the user is shown
# ---------------------------------------------------------------------------

def error_lines(rs: RuleSet) -> list:
    """One line per entry not in force; deny/ask errors say ERROR."""
    lines = []
    for source, which, text, reason in rs.errors:
        where = f"{which} " if which else ""
        shown = f"{text!r} " if text else ""
        tag = ("ERROR" if which in ("deny", "ask") or not which else
               "ignored")
        lines.append(f"{tag}: {where}{shown}in {source}: {reason}")
    return lines


def rule_origin(rs: RuleSet, rule: Rule) -> str:
    """"managed", "project" or "user": the layer a rule comes from."""
    if rule.managed:
        return "managed"
    return "project" if rule.source == rs.project_path else "user"


EVERY_REPOSITORY = "every repository"


def rule_scope(rs: RuleSet, rule: Rule) -> str:
    """Where a rule applies: its ``permissions.repos`` key, "this project"
    for a project file's rule, else "every repository"."""
    if rule.scope:
        return rule.scope
    if not rule.managed and rule.source and rule.source == rs.project_path:
        return "this project"
    return EVERY_REPOSITORY


def _repo_line(rs: RuleSet) -> str:
    if rs.repo_key:
        where = f"repository {rs.repo_key}"
        if rs.repo_path and rs.repo_path != rs.repo_key:
            where += f" at {rs.repo_path}"
        return where
    return "repository: none (not a git checkout; only path keys can match)"


def scope_notice(rs: RuleSet) -> Optional[str]:
    """The start notice for global allow rules, or None.

    Shown while ``permissions.global_allow`` is unset and the user's
    top-level ``allow`` list holds a rule: those rules answer in every
    repository. ``true`` keeps them without the notice, ``false`` turns
    them off; ``permissions.repos`` scopes a rule to one repository."""
    if rs.global_allow is not None:
        return None
    glob = [r.text for r in rs.allow
            if not r.managed and not r.scope
            and rule_scope(rs, r) == EVERY_REPOSITORY]
    if not glob:
        return None
    key = rs.repo_key or "github.com/owner/repo"
    return ("allow rules in " + (rs.user_path or "~/.delfin_settings.json")
            + f" apply in {EVERY_REPOSITORY}: " + ", ".join(glob)
            + f". Scope them with permissions.repos {{\"{key}\": "
            "{\"allow\": [...]}}, or set permissions.global_allow to true "
            "(keep, no notice) or false (off).")


def rule_records(rs: RuleSet) -> list:
    """Every rule in force as ``{list, rule, origin, source, scope}``, in
    the order the gate evaluates them; ``only_allow`` entries carry
    ``list: "only_allow"``. ``scope`` is ``rule_scope``. For the boundary
    report."""
    out = []
    for which in LISTS:
        for r in getattr(rs, which):
            out.append({"list": which, "rule": r.text,
                        "origin": rule_origin(rs, r), "source": r.source,
                        "scope": rule_scope(rs, r)})
    for group in rs.only_allow:
        for r in group:
            out.append({"list": "only_allow", "rule": r.text,
                        "origin": "managed", "source": r.source,
                        "scope": EVERY_REPOSITORY})
    return out


def describe(rs: RuleSet) -> str:
    """The read-only view ``/permissions`` prints: the managed policy's
    rules and the user's and the project's, each with its source, and
    for every allow the repository it applies in."""
    lines = [f"permission rules (edit {rs.user_path or '~/.delfin_settings.json'}"
             ", key \"permissions\"; the agent cannot write it)",
             "this workspace: " + _repo_line(rs)]
    if rs.project_path:
        lines.append(f"project rules: {rs.project_path} (deny and ask only)")
    pol = rs.policy
    for src in (getattr(pol, "sources", ()) or ()):
        lines.append(f"managed rules: {src.path} (administrator; "
                     + ("trusted" if src.trusted else
                        "NOT trusted, restrictions only") + ")")
    for which in LISTS:
        rules = getattr(rs, which)
        lines.append(f"  {which}:" + ("" if rules else " (none)"))
        for r in rules:
            origin = rule_origin(rs, r)
            tag = f"managed {r.source}" if origin == "managed" else origin
            if r.scope:
                tag += f", repository {r.scope}"
            elif which == "allow" and origin == "user":
                tag += f", {EVERY_REPOSITORY}"
            lines.append(f"    {r.text}  [{tag}]")
    if rs.allow_off:
        lines.append("  allow not in force (global_allow is false): "
                     + ", ".join(r.text for r in rs.allow_off))
    if rs.other_repos:
        lines.append("  rules for other repositories (not this repository):")
        lines.extend(f"    {r.list} {r.text}  [repository {r.scope}]"
                     for r in rs.other_repos)
    notice = scope_notice(rs)
    if notice:
        lines.append(f"  note: {notice}")
    for group in rs.only_allow:
        src = group[0].source if group else "managed policy"
        lines.append("  only_allow (an allow outside it grants nothing): "
                     + (", ".join(r.text for r in group) or "nothing")
                     + f"  [managed {src}]")
    errs = error_lines(rs)
    for err in (getattr(pol, "errors", ()) or ()):
        errs.append(f"ERROR: managed policy: {err}")
    if errs:
        lines.append("  not in force:")
        lines.extend(f"    {e}" for e in errs)
    lines.append("  deny > ask > allow, managed rules before the user's; deny "
                 "and ask hold in every mode; allow only skips a question and "
                 "never lifts a refusal.")
    return "\n".join(lines)
