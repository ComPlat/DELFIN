"""What an auto iteration may not do, decided per tool call.

Input: one tool call (bare or ``mcp__server__`` name and its arguments)
made while an auto iteration runs, the workspace, the loop's ``measure``
and ``touch`` lists. Output of :func:`classify_call`: ``None`` -- the call
goes on to the ordinary gate, which decides as it always does -- or a
:class:`Verdict`:

* ``outward``: the call would publish or send something: ``git push``,
  ``gh pr create|merge|comment|...``, issue or release creation, a
  mutating ``gh api``, branch deletion, ``git merge``, leaving the loop's
  branch for the default branch, an upload (curl/wget with a body, scp or
  rsync to a host, ssh), mail, session messages, notifications, webhooks,
  ``worktree_merge``, a future turn outside the loop's budget
  (``cron_create``, ``schedule_wakeup``), an MCP tool whose name sends,
  creates, merges or deletes. Queued for the user (``queue=True``).
* ``widen``: a permission, a writable root or a mode change the user did
  not make (``remember_permission*``, ``exit_plan_mode``,
  ``enter_worktree``). Queued.
* ``measure``: a write to a file the iteration is judged by. Refused.
* ``scope``: a write outside ``touch`` (when given) or into ``.git``.
  Refused.

Semantics. This runs IN ADDITION to the gate, never instead of it: a
``None`` here is not an approval. It over-refuses on purpose -- ``echo
"git push"`` is refused too -- because the cost of a false refusal in an
unattended iteration is one queued item the user reads, and the cost of a
missed one is a publication nobody approved. Shell bodies (``bash -c
'...'``) are unwrapped with the gate's own ``_with_shell_bodies`` when
api_client is importable.

:func:`measure_findings` judges the iteration's file changes at its end:
``hard`` for a declared measure file (the harness reverts the iteration),
soft for what ``edit_signals.classify`` (the experiment-discipline bundle,
when importable) or the built-in rules find in test and grader files (the
iteration is kept and flagged for the user). :func:`mode_rank` orders the
permission postures, which nothing else in the package did.

Never raises: an exception in a classifier reads as "refuse".
"""

from __future__ import annotations

import re
import shlex
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Optional

__all__ = ["Verdict", "Finding", "classify_call", "measure_findings",
           "mode_rank", "widened", "measure_source"]


@dataclass(frozen=True)
class Verdict:
    kind: str           # outward | widen | measure | scope
    reason: str
    queue: bool


@dataclass(frozen=True)
class Finding:
    path: str
    kind: str
    detail: str
    hard: bool


_SHELL_TOOLS = frozenset({"bash", "bash_background"})
_CMD_KEYS = ("command", "cmd", "script", "code")
_OUTWARD_TOOLS = {
    "session_message": "sends a message to another session",
    "push_notification": "sends a notification",
    "remote_trigger": "calls a webhook",
    "cron_create": "schedules a turn outside this loop's budget",
    "schedule_wakeup": "schedules a turn outside this loop's budget",
    "worktree_merge": "merges a worktree and deletes its branch",
    # DELFIN's HPC clusters (cluster_tools.py): each acts on the user's
    # account at a computing centre, which a person approves every time.
    "cluster_submit": "submits jobs to an HPC cluster",
    "cluster_cancel": "cancels jobs on an HPC cluster",
    "cluster_upload": "copies files to an HPC cluster",
    "cluster_run": "runs a command on an HPC cluster's login node",
}
_WIDEN_TOOLS = {
    "remember_permission": "adds a permanent permission",
    "remember_permission_bundle": "adds permanent permissions",
    "exit_plan_mode": "changes the permission mode",
    "enter_worktree": "adds a writable root",
}
_WRITE_TOOLS = frozenset({"write_file", "edit_file", "multi_edit",
                          "notebook_edit", "apply_patch", "edit_sheet",
                          "fill_docx_template", "fill_pdf_form",
                          "create_docx", "create_pdf", "make_plot",
                          "undo_changes", "cluster_fetch"})
_PATH_KEYS = ("path", "file_path", "notebook_path", "output_path", "target")
_MCP_SENDING = re.compile(
    r"(send|post|create|merge|delete|remove|publish|comment|reply|push|"
    r"upload|invite|share|close|reopen|approve|update|write|edit|set)",
    re.IGNORECASE)

_G = r"(?:^|[\s;&|(`'\"])(?:\S*/)?"
_GIT = _G + r"git(?:\s+-[Cc]\s+\S+|\s+--[\w-]+(?:=\S+)?)*\s+"
_RULES: tuple[tuple[re.Pattern, str], ...] = tuple(
    (re.compile(p, re.IGNORECASE if i else 0), why) for i, (p, why) in
    enumerate([
        (_GIT + r"push\b", "pushes to a remote"),
        (_G + r"gh\s+pr\s+(create|merge|comment|review|edit|close|reopen|ready)\b",
         "acts on a pull request"),
        (_G + r"gh\s+(issue|release|gist|repo|secret|variable|label)\s+"
              r"(create|edit|close|delete|comment|upload|fork|rename|set|"
              r"archive|transfer)\b", "creates or changes something on the forge"),
        (_G + r"gh\s+api\b[^;&|\n]*(?:-X\s*(?:POST|PUT|PATCH|DELETE)|"
              r"--method\s+(?:POST|PUT|PATCH|DELETE)|\s-[fF]\s|--field|"
              r"--raw-field|--input)", "sends a forge API request"),
        (_GIT + r"branch\b[^;&|\n]*\s(?:-d|-D|--delete)\b", "deletes a branch"),
        (_GIT + r"merge\b", "merges"),
        (_G + r"curl\b[^;&|\n]*\s(?:-d|--data[\w-]*|-F|--form|-T|"
              r"--upload-file|-X\s*(?:POST|PUT|PATCH|DELETE)|"
              r"--request\s+(?:POST|PUT|PATCH|DELETE))", "uploads data"),
        (_G + r"wget\b[^;&|\n]*--(?:post|body)-(?:data|file)", "uploads data"),
        (_G + r"(?:mail|mailx|sendmail|mutt|ssmtp|msmtp)\b", "sends mail"),
    ]))
_REMOTE_PROGRAMS = ("ssh", "scp", "sftp", "rsync")
_SSH_ARG_OPTS = set("bcDEeFIiJLlmOoPpQRSWw")
_SCP_ARG_OPTS = set("cFiJloPSXD")
_SEGMENT_SPLIT = re.compile(r"\|\||&&|[;|&\n()`]")


def ssh_host_allowed(host: str, entries: Iterable[str]) -> bool:
    """*host* fits *entries*: a plain entry exactly (case folded), an
    explicit ``*.domain`` entry any subdomain of domain (not domain
    itself). No other wildcard matches anything."""
    h = str(host or "").strip().lower().rstrip(".")
    if not h:
        return False
    for raw in entries or ():
        e = str(raw or "").strip().lower().rstrip(".")
        if e.startswith("*.") and len(e) > 2 and "*" not in e[2:]:
            if h.endswith(e[1:]) and h != e[2:]:
                return True
        elif e and "*" not in e and h == e:
            return True
    return False


def _entry_within(entry: str, ceiling: Iterable[str]) -> bool:
    """Every host *entry* admits is admitted by *ceiling*."""
    e = str(entry or "").strip().lower().rstrip(".")
    if e.startswith("*."):
        dom = e[2:]
        for c in ceiling:
            c = str(c).strip().lower().rstrip(".")
            if c.startswith("*.") and (dom == c[2:] or dom.endswith("." + c[2:])):
                return True
        return False
    return ssh_host_allowed(e, ceiling)


def allowed_ssh_hosts(settings: Optional[dict] = None) -> list[str]:
    """The user's ``agent.auto_ssh_hosts`` (a list of strings; anything else
    is no list), narrowed by every ``auto_ssh_hosts_max`` ceiling of the
    managed policy. Never raises; an error allows no host."""
    try:
        if settings is None:
            from delfin.user_settings import load_settings
            settings = load_settings() or {}
        raw = ((settings or {}).get("agent") or {}).get("auto_ssh_hosts")
        if not isinstance(raw, list):
            return []
        entries = [str(x).strip() for x in raw
                   if isinstance(x, str) and x.strip()]
        from . import managed_policy
        for ceiling in managed_policy.load().auto_ssh_hosts_max:
            entries = [e for e in entries if _entry_within(e, ceiling)]
        return entries
    except Exception:
        return []


def _host_of(spec: str) -> str:
    """``[user@]host[:...]`` -> host; bracketed IPv6 kept without brackets."""
    s = spec.split("@", 1)[-1]
    if s.startswith("["):
        return s[1:s.find("]")] if "]" in s else ""
    return s.split(":", 1)[0]


def _remote_hosts(text: str) -> tuple[list[str], str]:
    """(destination hosts, why it cannot be told) of every ssh, scp, sftp
    and rsync in *text*; ([], "") when there is none."""
    hosts: list[str] = []
    for segment in _SEGMENT_SPLIT.split(text):
        try:
            words = shlex.split(segment, posix=True)
        except ValueError:
            words = segment.split()
        idx = next((i for i, w in enumerate(words)
                    if w.rsplit("/", 1)[-1] in _REMOTE_PROGRAMS), None)
        if idx is None:
            continue
        prog = words[idx].rsplit("/", 1)[-1]
        args = words[idx + 1:]
        if prog == "rsync":
            i = 0
            while i < len(args):
                a = args[i]
                rsh = None
                if a == "-e" or a == "--rsh":
                    rsh = args[i + 1] if i + 1 < len(args) else ""
                    i += 1
                elif a.startswith("--rsh="):
                    rsh = a.split("=", 1)[1]
                elif a.startswith("-e") and len(a) > 2:
                    rsh = a[2:]
                if rsh is not None:
                    if any(t.startswith("-") for t in rsh.split()[1:]):
                        return hosts, "an rsync remote shell with options"
                elif a.startswith("rsync://"):
                    hosts.append(_host_of(a[len("rsync://"):].split("/", 1)[0]))
                elif not a.startswith("-") and ":" in a and not a.startswith(("/", "./", "../")):
                    hosts.append(_host_of(a))
                i += 1
            continue
        opts = _SSH_ARG_OPTS if prog == "ssh" else _SCP_ARG_OPTS
        i = 0
        dest_seen = False
        while i < len(args):
            a = args[i]
            if a.startswith("-") and len(a) > 1 and not dest_seen:
                j = 1
                while j < len(a):
                    ch = a[j]
                    if ch in opts:
                        val = a[j + 1:] or (args[i + 1] if i + 1 < len(args) else "")
                        if not a[j + 1:]:
                            i += 1
                        if ch == "J":
                            hosts.extend(_host_of(h) for h in val.split(",") if h)
                        elif ch == "o":
                            key = val.split("=", 1)[0].strip().lower()
                            if key in ("proxycommand", "localcommand"):
                                return hosts, f"an ssh {key}"
                            if key == "proxyjump":
                                hosts.extend(_host_of(h) for h in
                                             val.split("=", 1)[-1].split(",") if h)
                        break
                    j += 1
            elif prog == "ssh":
                if not dest_seen:
                    hosts.append(_host_of(a.replace("ssh://", "", 1)))
                    dest_seen = True       # the rest is the remote command
            else:
                if prog == "sftp" and not dest_seen:
                    hosts.append(_host_of(a.replace("sftp://", "", 1)))
                    dest_seen = True
                elif ":" in a and not a.startswith(("/", "./", "../")):
                    hosts.append(_host_of(a))
            i += 1
        if prog in ("ssh", "sftp") and not dest_seen:
            return hosts, f"{prog} without a destination it can read"
    return [h for h in hosts if h] if hosts else [], ""


_LEAVE_BRANCH = re.compile(
    _GIT + r"(?:switch|checkout)\s+(?!-[bBc]\b)(?:--\S+\s+)*(\S+)")
_DEFAULT_NAMES = frozenset({"main", "master", "trunk"})
_SHELL_WRITE = re.compile(
    r"(\btee\b|\bsed\s+(?:-\S*\s+)*-i|\bperl\s+-\S*i|\brm\b|\bmv\b|"
    r"\bcp\b|\btruncate\b|\bchmod\b|\bln\b|\bdd\b|\bpatch\b|\bgit\s+"
    r"(?:checkout|restore|rm|mv)\b|\bpython\S*\s+-c\b|\bawk\b[^|]*-i\b)")

_REDIRECT = re.compile(r"(?<![0-9&<])>>?\|?\s*([^\s;&|<>]+)|[0-9]>>?\s*([^\s&;|<>]+)")

_last_source = "builtin"

# The protected zone: what keeps an unattended iteration honest and
# contained. An iteration never edits it; a refused edit becomes a change
# request (change_requests) for the supervising person instead.
_GUARD_MODULES = (
    "sandbox", "approval_guard", "approval_code", "approval_answers",
    "cli_gate", "permission_rules", "managed_policy", "principles_guard",
    "landlock_exec", "seatbelt", "socket_guard", "process_guard",
    "egress_proxy", "mcp_isolation", "isolation_posture", "workspace_trust",
    "launch_guard", "verify_guard", "edit_signals", "experiment",
    "experiment_hints", "paired_compare", "comparison_claims", "auto_guard",
    "auto_mode", "loops", "change_requests", "boundaries", "untrusted",
    "output_guard", "contained_run", "terminal_confirm", "file_confirm",
    "kit_confirm", "security_events", "audit_log", "gate_coverage",
    "verification_claims", "turn_verification", "verify_commands",
    "caged_test_runners")
_ZONE_RULES: tuple[tuple[re.Pattern, str], ...] = (
    (re.compile(r"(^|/)delfin/agent/(" + "|".join(_GUARD_MODULES) + r")\.py$"),
     "a sandbox, approval, gate or honesty guard of DELFIN"),
    (re.compile(r"(^|/)delfin/agent/pack/benchmark/"),
     "benchmark tasks and acceptance scripts"),
    (re.compile(r"(^|/)[^/]*principles[^/]*$"), "the principles"),
    (re.compile(r"(^|/)managed_policy\.json$"), "the managed policy"),
    (re.compile(r"(^|/)\.delfin/(settings[^/]*\.json|hooks[^/]*|"
                r"permissions[^/]*|mcp[^/]*\.json|managed[^/]*)$"),
     "DELFIN's settings, hooks or permissions for this workspace"),
    (re.compile(r"(^|/)(\.github/workflows/|\.circleci/|\.gitlab-ci\.ya?ml$|"
                r"Jenkinsfile$|azure-pipelines\.ya?ml$|\.pre-commit-config\.ya?ml$)"),
     "CI configuration"),
)


def protected_reason(path: str) -> str:
    """Why *path* (workspace-relative) is in the protected zone, or ""."""
    rel = str(path or "").replace("\\", "/")
    while rel.startswith("./"):
        rel = rel[2:]
    if not rel:
        return ""
    for rx, why in _ZONE_RULES:
        if rx.search(rel):
            return why
    try:
        from . import edit_signals as _es
        if _es.role(rel) == "grader":
            return ("a grader file: test-runner configuration, a benchmark, "
                    "an acceptance script, a threshold or a snapshot")
        return ""
    except Exception:
        return ("a grader file (test-runner, benchmark, acceptance or "
                "threshold)") if _GRADER_PATH.search(rel) else ""


def _bare(name: str) -> str:
    name = str(name or "")
    return name.rsplit("__", 1)[-1] if name.startswith("mcp__") else name


def _unwrapped(cmd: str) -> str:
    try:
        from .api_client import _with_shell_bodies
        return _with_shell_bodies(cmd)
    except Exception:
        return cmd


def _norm(path: Any, workspace: str) -> str:
    """*path* relative to the workspace, posix, without ./; '' if outside."""
    s = str(path or "").strip()
    if not s:
        return ""
    p = PurePosixPath(s.replace("\\", "/"))
    if p.is_absolute():
        try:
            p = p.relative_to(PurePosixPath(str(workspace).rstrip("/") or "/"))
        except ValueError:
            return s           # outside: compare as given
    parts = [x for x in p.parts if x not in (".", "")]
    return "/".join(parts)


def _covers(prefix: str, rel: str) -> bool:
    prefix = prefix.rstrip("/")
    return rel == prefix or rel.startswith(prefix + "/")


def _write_paths(name: str, args: dict, workspace: str) -> list[str]:
    if name == "cluster_fetch":
        # Its path is on the cluster; it writes where dest says (default
        # cluster/<name>/...).
        return [p for p in [_norm(args.get("dest") or "cluster", workspace)] if p]
    out = []
    for k in _PATH_KEYS:
        if isinstance(args.get(k), str):
            out.append(_norm(args[k], workspace))
    if name == "apply_patch":
        text = str(args.get("patch") or args.get("diff") or "")
        for m in re.finditer(r"^(?:\+\+\+|---)\s+(?:[ab]/)?(\S+)", text,
                             re.MULTILINE):
            if m.group(1) != "/dev/null":
                out.append(_norm(m.group(1), workspace))
    return [p for p in dict.fromkeys(out) if p]


def _command(args: dict) -> str:
    for k in _CMD_KEYS:
        if isinstance(args.get(k), str):
            return args[k]
    return ""


def _shell_verdict(cmd: str, workspace: str, measure: Iterable[str],
                   default_branches: Iterable[str],
                   ssh_hosts: Iterable[str] = ()) -> Optional[Verdict]:
    text = _unwrapped(cmd)
    for rx, why in _RULES:
        if rx.search(text):
            return Verdict("outward", f"the command {why}", True)
    if re.search(_G + r"(?:ssh|scp|sftp|rsync)\b", text):
        hosts, unclear = _remote_hosts(text)
        if unclear:
            return Verdict("outward", f"the command reaches a host through "
                           f"{unclear}, which cannot be checked", True)
        allowed = list(ssh_hosts or ())
        refused = [h for h in hosts if not ssh_host_allowed(h, allowed)]
        if refused or not hosts:
            return Verdict("outward", "the command reaches "
                           + (", ".join(sorted(set(refused))) or "a host")
                           + ", which is not in agent.auto_ssh_hosts", True)
    defaults = set(default_branches) | _DEFAULT_NAMES
    for m in _LEAVE_BRANCH.finditer(text):
        target = m.group(1).strip("'\"")
        if target in defaults:
            return Verdict("outward", f"the command switches to {target}, "
                           "away from the loop's branch", True)
    named: set[str] = set()
    for m in _REDIRECT.finditer(text):
        target = (m.group(1) or m.group(2) or "").strip("'\"")
        if target and not target.startswith("&") and target != "/dev/null":
            named.add(_norm(target, workspace))
    if _SHELL_WRITE.search(text):
        try:
            words = shlex.split(text, posix=True)
        except ValueError:
            words = text.split()
        named |= {_norm(w.strip("'\"<>|;&"), workspace) for w in words}
    if named:
        for m in measure:
            mm = _norm(m, workspace)
            if mm and any(n == mm or _covers(mm, n) for n in named if n):
                return Verdict("measure", f"the command may change {mm}, "
                               "which this loop is judged by", False)
    for n in sorted(x for x in named if x):
        why = protected_reason(n)
        if why:
            return Verdict("protected", f"the command may change {n}, which "
                           f"is in the protected zone ({why})", False)
    return None


def classify_call(name: str, args: Optional[dict], *, workspace: str = "",
                  measure: Iterable[str] = (), touch: Iterable[str] = (),
                  default_branches: Iterable[str] = (),
                  ssh_hosts: Iterable[str] = ()) -> Optional[Verdict]:
    """See the module docstring. ``ssh_hosts``: the hosts ssh, scp, sftp and
    rsync may reach (allowed_ssh_hosts); none by default. Never raises (an
    error refuses)."""
    try:
        return _classify(name, dict(args or {}), str(workspace or ""),
                         list(measure or ()), list(touch or ()),
                         list(default_branches or ()),
                         list(ssh_hosts or ()))
    except Exception as exc:
        return Verdict("outward", f"the call could not be classified ({exc})",
                       True)


def _classify(name, args, workspace, measure, touch, defaults, ssh_hosts=()):
    full = str(name or "")
    bare = _bare(full)
    if bare in _SHELL_TOOLS or (full.startswith("mcp__") and _command(args)):
        return _shell_verdict(_command(args), workspace, measure, defaults,
                              ssh_hosts)
    if bare in _OUTWARD_TOOLS:
        return Verdict("outward", _OUTWARD_TOOLS[bare], True)
    if bare in _WIDEN_TOOLS:
        return Verdict("widen", _WIDEN_TOOLS[bare], True)
    if full.startswith("mcp__") and _MCP_SENDING.search(bare):
        return Verdict("outward", f"the MCP tool {bare} may send or change "
                       "something outside this machine", True)
    if bare in _WRITE_TOOLS:
        for rel in _write_paths(bare, args, workspace):
            for m in measure:
                if _covers(_norm(m, workspace), rel):
                    return Verdict("measure", f"{rel} is part of the measure "
                                   "this loop is judged by", False)
            if rel == ".git" or rel.startswith(".git/"):
                return Verdict("scope", f"{rel} is git's own state", False)
            why = protected_reason(rel)
            if why:
                return Verdict("protected", f"{rel} is in the protected zone "
                               f"({why})", False)
            if touch and not any(_covers(_norm(t, workspace), rel)
                                 for t in touch):
                return Verdict("scope", f"{rel} is outside what this loop may "
                               "change (touch: " + ", ".join(touch) + ")",
                               False)
    return None


# ---------------------------------------------------------------------------
# Postures
# ---------------------------------------------------------------------------

_RANK = {"plan": 0, "default": 1, "diff_approval": 1, "acceptEdits": 2,
         "bypassPermissions": 3}


def mode_rank(mode: str) -> int:
    """plan < default = diff_approval < acceptEdits < bypassPermissions; an
    unknown mode ranks as the widest, so it can never pass for narrower."""
    return _RANK.get(str(mode or ""), 3)


def widened(before: str, after: str) -> bool:
    return mode_rank(after) > mode_rank(before)


# ---------------------------------------------------------------------------
# The iteration's changes, judged at its end
# ---------------------------------------------------------------------------

_TEST_PATH = re.compile(r"(^|/)(tests?/|test_[^/]*$|[^/]*_test\.[a-z]+$|"
                        r"conftest\.py$|[^/]*\.(spec|test)\.[jt]sx?$)")
_GRADER_PATH = re.compile(r"(^|/)(\.github/workflows/|pytest\.ini$|tox\.ini$|"
                          r"setup\.cfg$|\.coveragerc$|[^/]*thresholds?[^/]*$|"
                          r"[^/]*acceptance[^/]*$|__snapshots__/)")
_ASSERT = re.compile(r"^\s*(assert\b|self\.assert|expect\(|assert_)", re.M)
_SKIP = re.compile(r"(pytest\.mark\.(skip|xfail)|@unittest\.skip|\.skip\(|"
                   r"xit\(|it\.skip|# *noqa|pytest\.skip\()")


def measure_source() -> str:
    """Which detector the last :func:`measure_findings` call used:
    ``edit_signals`` or ``builtin``."""
    return _last_source


def _builtin(path: str, old: str, new: str) -> list[tuple[str, str]]:
    out = []
    if _TEST_PATH.search(path):
        if new == "" and old:
            out.append(("test_deleted", "the test file was emptied or deleted"))
        if len(_ASSERT.findall(new)) < len(_ASSERT.findall(old)):
            out.append(("assertion_removed", "an assertion was removed"))
        if len(_SKIP.findall(new)) > len(_SKIP.findall(old)):
            out.append(("skip_added", "a skip or xfail was added"))
    elif _GRADER_PATH.search(path) and old != new:
        out.append(("grader_changed", "a grader or threshold file changed"))
    return out


def measure_findings(changes: dict, *, measure: Iterable[str] = (),
                     session_key: str = "",
                     workspace: str = "") -> list[Finding]:
    """Findings over ``{relpath: (old_text, new_text)}`` (see the module
    docstring); ``old`` None means the iteration created the file, ``new``
    None that it deleted it. A declared measure file is a hard finding
    whatever changed in it.

    With edit_signals importable (#269): ``classify`` per change, with its
    created/deleted flags, and ``measure_edits(session_key)`` -- the
    measure edits the write path journalled this turn, which a diff taken
    at the end cannot attribute (a grader changed and changed back, a write
    outside the diff's view) -- made workspace-relative and de-duplicated
    by (path, kind). Without it, the built-in rules."""
    global _last_source
    declared = [str(m).rstrip("/") for m in measure or ()]
    try:
        from . import edit_signals as _es
        classify = getattr(_es, "classify", None)
    except Exception:
        _es, classify = None, None
    _last_source = "edit_signals" if callable(classify) else "builtin"
    out: list[Finding] = []
    seen: set = set()

    def _add(f: Finding) -> None:
        if (f.path, f.kind) not in seen:
            seen.add((f.path, f.kind))
            out.append(f)

    for path, (old, new) in sorted((changes or {}).items()):
        if old == new:
            continue
        if any(_covers(m, path) for m in declared):
            _add(Finding(path, "measure_file",
                         "a file this loop is judged by changed", True))
            continue
        pairs: list[tuple[str, str]] = []
        if callable(classify):
            try:
                res = classify(path, old, new, created=old is None,
                               deleted=new is None) or {}
                pairs = [(str(k), str(d)) for k, d in res.get("measure", [])]
            except Exception:
                pairs = _builtin(path, old or "", new or "")
        elif old is not None:
            pairs = _builtin(path, old, new or "")
        for kind, detail in pairs:
            _add(Finding(path, kind, detail, False))
    if _es is not None and session_key:
        try:
            for sig in _es.measure_edits(session_key) or []:
                rel = _norm(getattr(sig, "path", ""), workspace)
                hard = any(_covers(m, rel) for m in declared)
                _add(Finding(rel, "measure_file" if hard else
                             str(getattr(sig, "kind", "")),
                             str(getattr(sig, "detail", "") or ""), hard))
        except Exception:
            pass
    return out


def relative(path: str, workspace: str) -> str:
    """Public alias of the normaliser, for the controller."""
    return _norm(path, workspace) or str(Path(path))


# ---------------------------------------------------------------------------
# Verification claims against the harness's own check
# ---------------------------------------------------------------------------

_CLAIM_RE = re.compile(
    r"[^.\n]*\b(tests? (?:now )?pass(?:es|ed)?|all (?:tests )?(?:green|pass)|"
    r"\d+ passed|verified|works now|it works|funktioniert|alle tests "
    r"bestanden|verifiziert)\b[^.\n]*", re.IGNORECASE)
_TESTISH = re.compile(r"\b(pytest|unittest|tox|nox|test|tests|ctest|"
                      r"cargo\s+test|go\s+test|npm\s+test|bun\s+test)\b")


def claim_findings(answer: str, *, check_command: str,
                   check_exit: Optional[int], edited: bool
                   ) -> tuple[list[str], str]:
    """The answer's verification claims the harness's check run does not
    back, and which detector judged them (``verification_claims`` from the
    verify bundle when importable, else ``builtin``).

    Only the harness's own run is known here (the iteration's runs are the
    turn loop's), so only a contradiction is reported: a claim while the
    check failed. No check, or a passing one: nothing. Never raises."""
    if not check_command or check_exit is None or check_exit == 0:
        return [], "builtin"
    try:
        from . import verification_claims as _vc
        find, judge, run_cls = _vc.find_claims, _vc.judge, _vc.Run
    except Exception:
        find = judge = run_cls = None
    if find is not None:
        try:
            kind = "test" if _TESTISH.search(check_command) else "check"
            verdicts = judge(find(answer or ""),
                             [run_cls(kind=kind, outcome="failed")],
                             edited=bool(edited), tried_tests=True)
            return ([f"{v.claim.sentence} — {v.reason}" for v in verdicts],
                    "verification_claims")
        except Exception:
            pass
    found = [" ".join(m.group(0).split())[:200]
             for m in _CLAIM_RE.finditer(answer or "")]
    return found, "builtin"
