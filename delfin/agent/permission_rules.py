"""Structured permission rules: ``permissions.deny`` / ``ask`` / ``allow``.

Input: the ``permissions`` block of the user's settings file
(``~/.delfin_settings.json``, written only by the user) and of the
project file ``<workspace>/.delfin/settings.json``::

    "permissions": {"deny": [...], "ask": [...], "allow": [...]}

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
* A project file may add deny and ask rules; its allow list is ignored and
  reported. A repository must not grant itself permissions.
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


@dataclass
class RuleSet:
    deny: tuple = ()
    ask: tuple = ()
    allow: tuple = ()
    #: (source path, list, rule text, reason) for every entry not in force.
    errors: list = field(default_factory=list)
    user_path: str = ""
    project_path: str = ""

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
        if not host.isascii() or not _HOST_RE.match(host):
            return f"{host!r} is not a host name"
        return Rule(raw, "web_fetch", host)
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
                 errors: list) -> dict:
    out = {k: [] for k in LISTS}
    for key, entries in block.items():
        if key not in LISTS:
            errors.append((source, str(key), "",
                           f"unknown list {key!r} (deny, ask, allow)"))
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
                                 rule.wildcard, key, source))
    return out


def load(workspace: Any = None) -> RuleSet:
    """The rules in force for ``workspace``. Never raises: a failure is an
    error entry and grants nothing."""
    try:
        upath = user_settings_path()
    except Exception:
        upath = None
    ppath = project_settings_path(workspace) if workspace is not None else None
    key = (_stamp(upath), _stamp(ppath))
    with _CACHE_LOCK:
        hit = _CACHE.get(key)
    if hit is not None:
        return hit
    errors: list = []
    lists = {k: [] for k in LISTS}
    for path, project in ((upath, False), (ppath, True)):
        block = _read_block(path, errors)
        if block:
            got = _parse_lists(block, str(path), project=project, errors=errors)
            for k in LISTS:
                lists[k].extend(got[k])
    rs = RuleSet(tuple(lists["deny"]), tuple(lists["ask"]),
                 tuple(lists["allow"]), errors,
                 str(upath or ""), str(ppath or ""))
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


def _host_matches(rule_host: str, host: str) -> bool:
    host = (host or "").lower().rstrip(".")
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


def _bash_allow_hits(rule: Rule, cmd: str) -> bool:
    words = simple_command_words(cmd)
    if words is None:
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


def friction(rs: RuleSet, name: str, args: dict, workspace: Any
             ) -> "tuple[str, Rule] | None":
    """("deny", rule) or ("ask", rule) for the first rule the call meets."""
    if not (rs.deny or rs.ask):
        return None
    try:
        subj = _subjects(name, args, workspace)
    except Exception:
        subj = {"bash": None, "read": [], "write": [], "hosts": []}
    for which, rules in (("deny", rs.deny), ("ask", rs.ask)):
        for rule in rules:
            try:
                if _hits(rule, name, subj, workspace):
                    return which, rule
            except Exception:
                continue
    return None


def denies_path(rs: RuleSet, kind: str, path: Any, workspace: Any
                ) -> Optional[Rule]:
    """The deny rule on ``path`` for a ``kind`` ("read"/"write") access a
    tool resolved itself. A read rule covers writes too."""
    kinds = ("read",) if kind == "read" else ("read", "write")
    forms = _forms([Path(str(path))])
    for rule in rs.deny:
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
    command answers no read question it raises, and the reverse."""
    if not rs.allow:
        return None
    bare = bare_tool_name(name)
    args = args if isinstance(args, dict) else {}
    if site == "bash":
        if bare not in _BASH_TOOLS:
            return None
        cmd = str(args.get("command") or "")
        for rule in rs.allow:
            if rule.kind == "bash" and _bash_allow_hits(rule, cmd):
                return rule
            if (rule.kind == "tool" and _tool_matches(rule, name, strict=True)
                    and simple_command_words(cmd) is not None):
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
                    and _host_matches(rule.value, url_host(args.get("url")))):
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


def describe(rs: RuleSet) -> str:
    """The read-only view ``/permissions`` prints."""
    lines = [f"permission rules (edit {rs.user_path or '~/.delfin_settings.json'}"
             ", key \"permissions\"; the agent cannot write it)"]
    if rs.project_path:
        lines.append(f"project rules: {rs.project_path} (deny and ask only)")
    for which in LISTS:
        rules = getattr(rs, which)
        lines.append(f"  {which}:" + ("" if rules else " (none)"))
        for r in rules:
            src = "project" if r.source == rs.project_path else "user"
            lines.append(f"    {r.text}  [{src}]")
    errs = error_lines(rs)
    if errs:
        lines.append("  not in force:")
        lines.extend(f"    {e}" for e in errs)
    lines.append("  deny > ask > allow; deny and ask hold in every mode; allow "
                 "only skips a question and never lifts a refusal.")
    return "\n".join(lines)
