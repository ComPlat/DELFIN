"""The administrator's managed policy: restrictions no user layer can loosen.

Input: JSON files at fixed system paths, readable by every user and
writable only by an administrator::

    /etc/delfin/managed_policy.json                          (Linux)
    /Library/Application Support/DELFIN/managed_policy.json  (macOS)
    <installed delfin package>/managed_policy.json           (any host)

plus, optionally, the file named by ``DELFIN_MANAGED_POLICY``. Schema::

    {
      "permissions": {"deny": [...], "ask": [...], "only_allow": [...]},
      "disallowed_modes": ["bypassPermissions"],
      "require_isolation": true,
      "container_backend": "bwrap",
      "network": {"allowed_domains_max": ["pypi.org", "*.github.com"]},
      "disable_classifier_off": true,
      "audit_log": "/var/log/delfin/{user}.jsonl"
    }

Rule syntax is that of the user's permission rules: ``bash(git push)``,
``bash(pytest:*)``, ``read(<glob>)``, ``write(<glob>)``,
``web_fetch(domain:<host>)``, ``mcp__server[__tool|__*]`` and a bare tool
name. A relative glob is anchored at the session's workspace; ``*`` is one
path segment, ``**`` any number; a glob without a wildcard also covers
everything below it.

Output: ``load()`` -> ``Policy`` and the questions the gate asks of it --
``friction`` (the deny or ask rule a call meets), ``denies_path`` (a deny
rule on a path a tool resolved), ``allow_may_grant_bash`` (``only_allow``),
``coerce_mode`` / ``launch_decision`` (``disallowed_modes`` and
``require_isolation``), ``domain_within_max``, ``classifier_forced_on``,
``audit_destinations``, ``is_policy_path`` and ``summary_lines``.

Semantics:

* Precedence: managed > CLI flags > user settings > project settings for
  anything that loosens; every layer may tighten. Nothing in the schema
  loosens, so the policy is applied as a floor at the evaluation points
  (``api_client``): a managed deny comes before the session lists, the
  mode, hooks, user allow patterns and every approval, in every mode;
  a managed ask is put to the person in every mode and is a refusal with
  nobody to ask.
* Deny and ask match generously: a bash rule matches its words anywhere in
  any segment, quoted ``sh -c`` body or path-qualified program name; read
  and write rules match every path argument of a tool and every non-option
  word of a shell command. ``only_allow`` matches narrowly: a command is
  auto-allowed only when each of its segments starts with the words of an
  ``only_allow`` bash rule and it carries no substitution or redirection.
* A disallowed mode is replaced by ``default`` (``plan`` when ``default``
  is disallowed too) wherever a mode is set: the permissions object, a
  mode adopted inside a call, ``/permissions``. An interactive session
  started with the mode named explicitly on the command line refuses to
  start instead; an unattended run, and a mode read from a settings file,
  fall back and say so.
* ``require_isolation``: a shell command never runs unsandboxed (the
  ``off``/``auto`` settings are treated as ``bwrap``, whose fallbacks are
  Landlock, Seatbelt or a refusal). On a host without a working layer an
  interactive session refuses to start and an unattended run is capped at
  ``default``. ``container_backend: "bwrap"`` accepts only bubblewrap as
  that layer.
* Trust. A file is TRUSTED when it and every directory above it are owned
  by root or by the account that installed the package (the owner of the
  package directory, when that is not the current user), and none of them
  is group- or world-writable unless it is a sticky directory (a sticky
  directory lets no one replace an entry they do not own). An UNTRUSTED
  file -- user-owned, writable by others, or named by the environment
  variable, which the user controls -- is still read, but only for
  restrictions: it cannot set ``audit_log``, and because several files are
  merged by union of restrictions (and intersection of ``only_allow`` and
  domain ceilings) it can never remove what a trusted file imposes. Each
  untrusted file is announced at startup.
* A trusted file that exists but cannot be read or parsed fails closed:
  ``bypassPermissions`` and ``acceptEdits`` are disallowed and the error
  is shown. A malformed untrusted file is reported and ignored.
* No user setting names a policy path; the set of candidate paths is the
  module constant plus the environment variable.

Rejected: the environment variable as the only source (the user controls
the environment); filtering the advertised tool list (advertising is not
enforcement); trusting a file by mode bits alone (a user-owned 0644 file
is rewritten by its owner at will); treating an untrusted file as absent
(a restriction it adds costs nothing to honour and makes trialling a
policy possible before installing it as root).

Cost: one ``lstat`` per path component of each candidate file per call
that consults the policy (~10 syscalls); the files are parsed again only
when a stamp changes. No file -> ``Policy.active`` is False and every
question answers as if the module did not exist.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import stat
import sys
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

ENV_VAR = "DELFIN_MANAGED_POLICY"
FILE_NAME = "managed_policy.json"
LINUX_PATH = "/etc/delfin/" + FILE_NAME
MACOS_PATH = "/Library/Application Support/DELFIN/" + FILE_NAME
MODES = ("plan", "default", "diff_approval", "acceptEdits",
         "bypassPermissions")
#: Disallowed when a trusted file exists but cannot be used.
FAIL_CLOSED_MODES = frozenset({"bypassPermissions", "acceptEdits"})
#: The modes that run a command or a write nobody approved.
UNATTENDED_MODES = frozenset({"bypassPermissions", "acceptEdits"})
CONTAINER_BACKENDS = ("", "bwrap")
_KEYS = {"version", "permissions", "disallowed_modes", "require_isolation",
         "container_backend", "network", "disable_classifier_off",
         "audit_log", "comment"}

#: Tests set these; production leaves them None.
_PATHS_OVERRIDE: Optional[tuple] = None
_ADMIN_UIDS_OVERRIDE: Optional[frozenset] = None
_ISOLATION_PROBE: Optional[Callable[..., bool]] = None


# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

_KIND_ALIASES = {"bash": "bash", "shell": "bash", "read": "read",
                 "write": "write", "edit": "write", "web_fetch": "web_fetch",
                 "webfetch": "web_fetch"}
_RULE_RE = re.compile(r"^([A-Za-z_]+)\((.*)\)$", re.DOTALL)
_TOOL_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_\-]*$")
_MCP_RE = re.compile(r"^mcp__[A-Za-z0-9_\-]+(?:__(?:[A-Za-z0-9_\-]+|\*))?$")
_HOST_RE = re.compile(r"^[a-z0-9*](?:[a-z0-9.\-]*[a-z0-9])?$")
_NOT_SIMPLE = set("`$(){}<>\\\n\r")
_SEGMENT_SPLIT = re.compile(r"\|\||&&|[;|&\n]")
BASH_TOOLS = frozenset({"bash", "bash_background"})
WRITE_TOOLS = frozenset({
    "write_file", "edit_file", "multi_edit", "apply_patch", "notebook_edit",
    "make_plot", "edit_sheet", "fill_pdf_form", "fill_docx_template",
    "create_docx", "fill_series", "merge_pdfs", "split_pdf", "create_pdf",
    "draft_email",
})
_PATH_KEYS = ("path", "file_path", "filename", "file", "target",
              "notebook_path", "directory", "dir", "root", "source", "src",
              "dest", "destination", "output_path", "input_path", "pdf_path",
              "template_path", "output", "paths", "files", "sources")


@dataclass(frozen=True)
class Rule:
    """One parsed rule; ``value`` is words (bash), a glob (read/write), a
    host (web_fetch) or a tool name (tool)."""
    text: str
    kind: str
    value: Any
    wildcard: bool = False
    source: str = ""


def parse_rule(text: Any) -> "Rule | str":
    """A ``Rule``, or the reason ``text`` is not one."""
    if not isinstance(text, str) or not text.strip():
        return "a rule is a non-empty string"
    raw = text.strip()
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
            try:
                words = tuple(shlex.split(body))
            except ValueError as exc:
                return f"bash rule does not parse: {exc}"
            if not words or any(c in _NOT_SIMPLE or c in ";&|" for c in body):
                return "a bash rule names the words of one command"
            return Rule(raw, "bash", words, wildcard)
        if kind in ("read", "write"):
            return Rule(raw, kind, body)
        if not body.lower().startswith("domain:"):
            return "web_fetch rules are web_fetch(domain:<host>)"
        host = body[len("domain:"):].strip().lower().rstrip(".")
        if not host.isascii() or not _HOST_RE.match(host):
            return f"{host!r} is not a host name"
        return Rule(raw, "web_fetch", host)
    if raw.startswith("mcp__"):
        return (Rule(raw, "tool", raw) if _MCP_RE.match(raw) else
                "an MCP rule is mcp__<server>, mcp__<server>__<tool> or "
                "mcp__<server>__*")
    if _TOOL_RE.match(raw):
        return Rule(raw, "tool", raw)
    return ("not a rule (bash(...), read(...), write(...), "
            "web_fetch(domain:...), mcp__..., or a tool name)")


# ---------------------------------------------------------------------------
# Policy and loading
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Source:
    path: str
    origin: str          # "system", "install" or "environment"
    trusted: bool
    note: str = ""


@dataclass
class Policy:
    sources: tuple = ()
    deny: tuple = ()
    ask: tuple = ()
    #: One tuple of rules per file that sets it; a command must fit each.
    only_allow: tuple = ()
    disallowed_modes: frozenset = frozenset()
    require_isolation: bool = False
    container_backend: str = ""
    #: One tuple of domains per file that sets it; a host must fit each.
    allowed_domains_max: tuple = ()
    disable_classifier_off: bool = False
    audit_log: tuple = ()
    errors: list = field(default_factory=list)

    @property
    def active(self) -> bool:
        return bool(self.sources)


def install_path() -> Optional[Path]:
    """``managed_policy.json`` beside the installed ``delfin`` package."""
    try:
        import delfin
        return Path(delfin.__file__).resolve().parent / FILE_NAME
    except Exception:
        return None


def candidate_paths() -> list:
    """[(Path, origin)] in the order they are read."""
    if _PATHS_OVERRIDE is not None:
        out = [(Path(p), "system") for p in _PATHS_OVERRIDE]
    else:
        system = MACOS_PATH if sys.platform == "darwin" else LINUX_PATH
        out = [(Path(system), "system")]
        inst = install_path()
        if inst is not None:
            out.append((inst, "install"))
    env = os.environ.get(ENV_VAR, "").strip()
    if env:
        out.append((Path(env).expanduser(), "environment"))
    return out


def _admin_uids() -> frozenset:
    if _ADMIN_UIDS_OVERRIDE is not None:
        return frozenset(_ADMIN_UIDS_OVERRIDE)
    uids = {0}
    try:
        import delfin
        owner = os.stat(Path(delfin.__file__).resolve().parent).st_uid
        if owner != os.getuid():
            uids.add(owner)
    except Exception:
        pass
    return frozenset(uids)


def _chain(path: Path) -> list:
    out = [path]
    out.extend(path.parents)
    return out


def trust_of(path: Path) -> tuple:
    """(trusted, reason). Checks the path as written and as resolved."""
    admins = _admin_uids()
    try:
        real = Path(os.path.realpath(path))
    except Exception:
        return False, "cannot be resolved"
    seen = set()
    for p in _chain(Path(os.path.abspath(path))) + _chain(real):
        if str(p) in seen:
            continue
        seen.add(str(p))
        try:
            st = os.lstat(p)
        except OSError as exc:
            return False, f"{p} cannot be inspected ({exc.strerror})"
        if st.st_uid not in admins:
            return False, (f"{p} is owned by uid {st.st_uid}, not by root "
                           "or the account that installed DELFIN")
        if stat.S_ISLNK(st.st_mode):
            continue
        writable = st.st_mode & (stat.S_IWGRP | stat.S_IWOTH)
        if st.st_mode & stat.S_IWGRP and st.st_gid == 0:
            writable &= ~stat.S_IWGRP
        sticky = stat.S_ISDIR(st.st_mode) and st.st_mode & stat.S_ISVTX
        if writable and not sticky:
            return False, f"{p} is group- or world-writable"
    return True, ""


def _stamp(path: Path) -> tuple:
    try:
        st = os.stat(path)
        return (str(path), st.st_mtime_ns, st.st_size, st.st_ino,
                st.st_uid, st.st_mode)
    except OSError as exc:
        return (str(path), exc.errno)


_CACHE: dict = {}
_LOCK = threading.Lock()


def clear_cache() -> None:
    with _LOCK:
        _CACHE.clear()


def _strings(value) -> list:
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value if isinstance(v, str) and v.strip()]
    return []


def _rules(entries, label: str, source: str, errors: list) -> list:
    out = []
    if entries is None:
        return out
    if not isinstance(entries, list):
        errors.append(f"{source}: '{label}' is not a list; none of it is in force")
        return out
    for entry in entries:
        rule = parse_rule(entry)
        if isinstance(rule, str):
            errors.append(f"{source}: {label} rule {entry!r} is not in force: {rule}")
            continue
        out.append(Rule(rule.text, rule.kind, rule.value, rule.wildcard,
                        source))
    return out


def _merge(pol: Policy, data: dict, src: Source) -> None:
    errors = pol.errors
    where = src.path
    for key in data:
        if key not in _KEYS:
            errors.append(f"{where}: unknown key {key!r} is ignored")
    perms = data.get("permissions") or {}
    if not isinstance(perms, dict):
        errors.append(f"{where}: 'permissions' is not an object")
        perms = {}
    for key in perms:
        if key not in ("deny", "ask", "only_allow"):
            errors.append(f"{where}: permissions.{key} is not a managed list "
                          "(deny, ask, only_allow); a managed policy grants nothing")
    pol.deny += tuple(_rules(perms.get("deny"), "deny", where, errors))
    pol.ask += tuple(_rules(perms.get("ask"), "ask", where, errors))
    if perms.get("only_allow") is not None:
        pol.only_allow += (tuple(_rules(perms.get("only_allow"), "only_allow",
                                        where, errors)),)
    modes = set(_strings(data.get("disallowed_modes")))
    for m in sorted(modes - set(MODES)):
        errors.append(f"{where}: disallowed mode {m!r} is not a mode")
    pol.disallowed_modes = pol.disallowed_modes | frozenset(modes & set(MODES))
    if data.get("require_isolation"):
        pol.require_isolation = True
    backend = str(data.get("container_backend") or "").strip().lower()
    if backend and backend not in CONTAINER_BACKENDS:
        errors.append(f"{where}: container_backend {backend!r} is not "
                      "supported (bwrap); isolation is required instead")
        pol.require_isolation = True
    elif backend:
        pol.container_backend = backend
        pol.require_isolation = True
    net = data.get("network") or {}
    if isinstance(net, dict) and net.get("allowed_domains_max") is not None:
        pol.allowed_domains_max += (tuple(
            d.strip().lower().rstrip(".")
            for d in _strings(net.get("allowed_domains_max"))),)
    if data.get("disable_classifier_off"):
        pol.disable_classifier_off = True
    audit = str(data.get("audit_log") or "").strip()
    if audit:
        if not src.trusted:
            errors.append(f"{where}: audit_log is ignored: only a trusted "
                          "file chooses where the audit goes")
        elif not os.path.isabs(audit):
            errors.append(f"{where}: audit_log must be an absolute path")
        else:
            pol.audit_log += (audit,)


def load() -> Policy:
    """The policy in force. Never raises."""
    cands = candidate_paths()
    key = (tuple(_stamp(p) for p, _ in cands), _admin_uids())
    with _LOCK:
        hit = _CACHE.get(key)
    if hit is not None:
        return hit
    pol = Policy(errors=[])
    sources = []
    for path, origin in cands:
        if not os.path.lexists(path):
            if origin == "environment":
                pol.errors.append(f"{path} (from {ENV_VAR}) does not exist")
            continue
        trusted, why = (False, f"named by {ENV_VAR}, which the user sets") \
            if origin == "environment" else trust_of(path)
        if origin != "environment" and not trusted and not why:
            why = "not trusted"
        src = Source(str(path), origin, trusted, why)
        sources.append(src)
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
            if not isinstance(data, dict):
                raise ValueError("the top level is not an object")
        except Exception as exc:
            if trusted:
                pol.errors.append(
                    f"{path} cannot be used ({exc}); failing closed: "
                    + ", ".join(sorted(FAIL_CLOSED_MODES)) + " disallowed")
                pol.disallowed_modes |= FAIL_CLOSED_MODES
            else:
                pol.errors.append(f"{path} cannot be used ({exc}); ignored")
            continue
        _merge(pol, data, src)
    pol.sources = tuple(sources)
    with _LOCK:
        if len(_CACHE) > 16:
            _CACHE.clear()
        _CACHE[key] = pol
    return pol


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def bare_tool_name(name: str) -> str:
    text = str(name or "")
    return text.rsplit("__", 1)[-1] if text.startswith("mcp__") else text


def _tool_matches(value: str, name: str) -> bool:
    if value.startswith("mcp__"):
        if value.endswith("__*"):
            return name.startswith(value[:-1])
        if value.count("__") == 1:
            return name.startswith(value + "__")
        return name == value
    return name == value or bare_tool_name(name) == value


def _words(cmd: str) -> list:
    """Every word of ``cmd`` and of any quoted body inside it."""
    try:
        lex = shlex.shlex(cmd, posix=True, punctuation_chars=True)
        lex.whitespace_split = True
        toks = list(lex)
    except ValueError:
        toks = cmd.split()
    out = []
    for tok in toks:
        out.append(tok)
        if any(c.isspace() for c in tok) and tok != cmd:
            out.extend(_words(tok))
    return out


def _bash_hits(words: tuple, cmd: str) -> bool:
    toks = _words(cmd)
    norm = [os.path.basename(t) if "/" in t else t for t in toks]
    n = len(words)
    for seq in (toks, norm):
        for i in range(len(seq) - n + 1):
            if tuple(seq[i:i + n]) == words:
                return True
    return False


def _glob_regex(glob: str) -> str:
    out, i = [], 0
    while i < len(glob):
        if glob.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
            continue
        if glob.startswith("**", i):
            out.append(".*")
            i += 2
            continue
        c = glob[i]
        out.append("[^/]*" if c == "*" else "[^/]" if c == "?"
                   else re.escape(c))
        i += 1
    return "".join(out)


def _anchor(glob: str, workspace: Any) -> str:
    g = os.path.expanduser(glob.strip())
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
    literal = "/".join(head) or "/"
    try:
        literal = os.path.realpath(literal)
    except Exception:
        pass
    rest = parts[len(head):]
    return literal if not rest else literal.rstrip("/") + "/" + "/".join(rest)


def path_matches(glob: str, path: Any, workspace: Any) -> bool:
    """True when ``path`` matches ``glob`` or lies below a literal glob."""
    try:
        p = os.path.realpath(os.path.expanduser(str(path)))
    except Exception:
        return False
    pattern = _anchor(glob, workspace)
    if not any(ch in pattern for ch in "*?"):
        return p == pattern or p.startswith(pattern.rstrip("/") + "/")
    return re.fullmatch(_glob_regex(pattern), p) is not None


def _call_paths(bare: str, args: dict, workspace: Any) -> list:
    base = Path(workspace or ".")
    cwd = str(args.get("cwd") or "").strip()
    if cwd:
        base = (base / cwd) if not os.path.isabs(cwd) else Path(cwd)
    raw: list = []
    if bare in BASH_TOOLS:
        for tok in _words(str(args.get("command") or "")):
            if not tok or tok.startswith("-") and "=" not in tok:
                continue
            raw.append(tok.split("=", 1)[1] if tok.startswith("-") else tok)
    else:
        for key in _PATH_KEYS:
            val = args.get(key)
            if isinstance(val, str):
                raw.append(val)
            elif isinstance(val, (list, tuple)):
                raw.extend(v for v in val if isinstance(v, str))
    out = []
    for r in raw:
        r = os.path.expanduser(r)
        out.append(r if os.path.isabs(r) else str(base / r))
    return out


def _url_host(url: str) -> str:
    try:
        from urllib.parse import urlsplit
        return (urlsplit(str(url)).hostname or "").lower().rstrip(".")
    except Exception:
        return ""


def _host_in(host: str, entry: str) -> bool:
    entry = entry.strip().lower().rstrip(".")
    if entry.startswith("*."):
        return host.endswith(entry[1:]) and host != entry[2:]
    return host == entry or host.endswith("." + entry)


def _hits(rule: Rule, name: str, args: dict, workspace: Any) -> bool:
    bare = bare_tool_name(name)
    if rule.kind == "tool":
        return _tool_matches(str(rule.value), name)
    if rule.kind == "bash":
        return (bare in BASH_TOOLS
                and _bash_hits(rule.value, str(args.get("command") or "")))
    if rule.kind in ("read", "write"):
        if rule.kind == "write" and not (bare in WRITE_TOOLS
                                         or bare in BASH_TOOLS):
            return False
        return any(path_matches(rule.value, p, workspace)
                   for p in _call_paths(bare, args, workspace))
    if rule.kind == "web_fetch":
        return (bare == "web_fetch"
                and _host_in(_url_host(args.get("url") or ""), rule.value))
    return False


def friction(name: str, args: dict, workspace: Any,
             policy: Optional[Policy] = None) -> Optional[tuple]:
    """("deny"|"ask", Rule) for this call, or None."""
    pol = policy or load()
    if not pol.active:
        return None
    args = args if isinstance(args, dict) else {}
    for rule in pol.deny:
        if _hits(rule, name, args, workspace):
            return ("deny", rule)
    if (pol.allowed_domains_max and bare_tool_name(name) == "web_fetch"
            and not domain_within_max(_url_host(args.get("url") or ""), pol)):
        return ("deny", Rule("network.allowed_domains_max", "web_fetch",
                             "", source=_first_source(pol)))
    for rule in pol.ask:
        if _hits(rule, name, args, workspace):
            return ("ask", rule)
    return None


def denies_path(kind: str, path: Any, workspace: Any) -> Optional[Rule]:
    """The deny rule on ``path`` for a ``kind`` ("read"/"write") access."""
    pol = load()
    for rule in pol.deny:
        if rule.kind in ("read", kind) and rule.kind in ("read", "write") \
                and path_matches(rule.value, path, workspace):
            return rule
    return None


def allow_may_grant_bash(cmd: str) -> bool:
    """False when ``only_allow`` leaves ``cmd`` outside every grant."""
    pol = load()
    if not pol.only_allow:
        return True
    if any(c in _NOT_SIMPLE for c in cmd):
        return False
    segments = [s.strip() for s in _SEGMENT_SPLIT.split(cmd) if s.strip()]
    if not segments:
        return False
    for rules in pol.only_allow:
        for seg in segments:
            try:
                words = tuple(shlex.split(seg))
            except ValueError:
                return False
            if not any(r.kind == "bash" and (
                    words[:len(r.value)] == r.value if r.wildcard
                    else words == r.value) for r in rules):
                return False
    return True


def domain_within_max(host: str, policy: Optional[Policy] = None) -> bool:
    pol = policy or load()
    host = (host or "").strip().lower().rstrip(".")
    for ceiling in pol.allowed_domains_max:
        if not host or not any(_host_in(host, e) for e in ceiling):
            return False
    return True


def network_ceiling_active() -> bool:
    return bool(load().allowed_domains_max)


def classifier_forced_on() -> bool:
    return load().disable_classifier_off


def _first_source(pol: Policy) -> str:
    return pol.sources[0].path if pol.sources else ""


# ---------------------------------------------------------------------------
# Modes and isolation
# ---------------------------------------------------------------------------

def isolation_available(backend: str = "") -> bool:
    """True when a sandbox layer (``backend`` only, when named) works here."""
    if _ISOLATION_PROBE is not None:
        return bool(_ISOLATION_PROBE(backend))
    try:
        from . import api_client as ac
        if backend == "bwrap":
            return bool(ac._bwrap_functional())
        return bool(ac._bwrap_functional() or ac._landlock_functional()
                    or ac._seatbelt_functional())
    except Exception:
        return False


def isolation_unmet(policy: Optional[Policy] = None) -> bool:
    pol = policy or load()
    return pol.require_isolation and not isolation_available(
        pol.container_backend)


def disallowed_modes(policy: Optional[Policy] = None) -> frozenset:
    pol = policy or load()
    out = set(pol.disallowed_modes)
    if pol.require_isolation and isolation_unmet(pol):
        out |= UNATTENDED_MODES
    return frozenset(out)


def coerce_mode(mode: str) -> str:
    """``mode``, or the fallback when the policy disallows it."""
    pol = load()
    if not pol.active:
        return mode
    banned = disallowed_modes(pol)
    if mode not in banned:
        return mode
    return "default" if "default" not in banned else "plan"


def mode_refusal(mode: str) -> Optional[str]:
    """Why ``mode`` cannot be chosen, or None."""
    pol = load()
    if not pol.active or mode not in disallowed_modes(pol):
        return None
    if mode in pol.disallowed_modes:
        why = "the administrator's managed policy disallows it"
    else:
        why = ("the managed policy requires an isolation layer and none "
               "works on this host")
    return (f"{mode} is not available: {why} ({_first_source(pol)}). "
            "No user setting or flag changes that.")


def launch_decision(requested: str, *, explicit: bool,
                    attended: bool) -> tuple:
    """(mode, note, refuse) for a session starting in ``requested``.

    An attended session refuses to start when the person named a
    disallowed mode on the command line, or when isolation is required and
    missing; anything else falls back with a note."""
    pol = load()
    if not pol.active:
        return requested, "", False
    if attended and isolation_unmet(pol):
        note = ("managed policy requires an isolation layer "
                + (f"({pol.container_backend}) " if pol.container_backend
                   else "")
                + "and none works on this host; an interactive session does "
                f"not start ({_first_source(pol)}). Ask the administrator.")
        return requested, note, True
    refusal = mode_refusal(requested)
    if refusal is None:
        return requested, "", False
    if explicit and attended:
        return requested, f"refused --permission-mode {refusal}", True
    fallback = coerce_mode(requested)
    return fallback, (f"managed policy: {requested} -> {fallback}; "
                      + refusal), False


# ---------------------------------------------------------------------------
# The policy file and the audit
# ---------------------------------------------------------------------------

def is_policy_path(path: Any) -> bool:
    """True when ``path`` is one of the candidate policy files."""
    try:
        target = os.path.realpath(os.path.expanduser(str(path)))
    except Exception:
        return False
    for p, _ in candidate_paths():
        try:
            if target in (os.path.realpath(p), os.path.abspath(p)):
                return True
        except Exception:
            continue
    return False


def audit_destinations() -> list:
    """The managed audit files (``{user}``/``{uid}`` expanded)."""
    pol = load()
    try:
        import getpass
        user = getpass.getuser()
    except Exception:
        user = str(os.getuid())
    out = []
    for tpl in pol.audit_log:
        out.append(Path(tpl.replace("{user}", user)
                    .replace("{uid}", str(os.getuid()))))
    return out


def append_audit_line(line: str) -> None:
    """Append one record to every managed destination. Never raises."""
    for dest in audit_destinations():
        try:
            fd = os.open(dest, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o600)
            try:
                os.write(fd, line.encode("utf-8"))
            finally:
                os.close(fd)
        except Exception:
            continue


# ---------------------------------------------------------------------------
# Visibility
# ---------------------------------------------------------------------------

def summary(policy: Optional[Policy] = None) -> dict:
    pol = policy or load()
    return {
        "active": pol.active,
        "sources": [{"path": s.path, "origin": s.origin,
                     "trusted": s.trusted, "note": s.note}
                    for s in pol.sources],
        "deny": [r.text for r in pol.deny],
        "ask": [r.text for r in pol.ask],
        "only_allow": [[r.text for r in rs] for rs in pol.only_allow],
        "disallowed_modes": sorted(disallowed_modes(pol)),
        "require_isolation": pol.require_isolation,
        "isolation_unmet": isolation_unmet(pol),
        "container_backend": pol.container_backend,
        "allowed_domains_max": [list(c) for c in pol.allowed_domains_max],
        "disable_classifier_off": pol.disable_classifier_off,
        "audit_log": [str(p) for p in audit_destinations()],
        "errors": list(pol.errors),
    }


def summary_lines(policy: Optional[Policy] = None) -> list:
    """What the policy enforces, one line each; [] without a policy."""
    s = summary(policy)
    if not s["active"] and not s["errors"]:
        return []
    lines = []
    for src in s["sources"]:
        lines.append(f"managed policy {src['path']} ("
                     + ("trusted" if src["trusted"] else
                        f"NOT TRUSTED, restrictions only: {src['note']}")
                     + ")")
    if s["deny"]:
        lines.append("managed deny   " + ", ".join(s["deny"]))
    if s["ask"]:
        lines.append("managed ask    " + ", ".join(s["ask"]))
    for group in s["only_allow"]:
        lines.append("managed allow  auto-allow at most: "
                     + (", ".join(group) or "nothing"))
    if s["disallowed_modes"]:
        lines.append("managed modes  disallowed: "
                     + ", ".join(s["disallowed_modes"]))
    if s["require_isolation"]:
        lines.append("managed shell  isolation required"
                     + (f" ({s['container_backend']})"
                        if s["container_backend"] else "")
                     + (" -- NOT available on this host" if s["isolation_unmet"]
                        else ""))
    for group in s["allowed_domains_max"]:
        lines.append("managed net    domains at most: " + ", ".join(group))
    if s["disable_classifier_off"]:
        lines.append("managed watcher cannot be switched off")
    for dest in s["audit_log"]:
        lines.append(f"managed audit  {dest}")
    for err in s["errors"]:
        lines.append(f"managed ERROR  {err}")
    return lines


def deny_text(name: str, rule: Rule) -> str:
    return (f"'{bare_tool_name(name)}' is refused by the managed policy rule "
            f"`{rule.text}` ({rule.source or 'managed policy'}), set by the "
            "administrator of this installation. It holds in every mode and "
            "no user setting, flag or approval lifts it. Do not reach the "
            "same effect another way; tell the user what you needed it for.")


__all__ = [
    "ENV_VAR", "LINUX_PATH", "MACOS_PATH", "Policy", "Rule", "Source",
    "allow_may_grant_bash", "append_audit_line", "audit_destinations",
    "candidate_paths", "classifier_forced_on", "clear_cache", "coerce_mode",
    "denies_path", "deny_text", "disallowed_modes", "domain_within_max",
    "friction", "install_path", "is_policy_path", "isolation_available",
    "isolation_unmet", "launch_decision", "load", "mode_refusal",
    "network_ceiling_active", "parse_rule", "path_matches", "summary",
    "summary_lines", "trust_of",
]
