"""What the header's branch dropdown lists, and the user's own changes to it.

Input: a checkout (the shared one the dashboard runs from), an optional
``gh`` executable, and the per-user preferences file
``~/.delfin/branch_list.json``. Output: :class:`BranchEntry` rows, in this
order:

1. ``main``, always;
2. the user's pinned entries (any ref: a branch the rule below hides, a
   remote branch, a tag, a commit);
3. branches with an open pull request (``gh pr list``), with its number and
   CI state; an open pull request whose head has no local branch is offered
   as ``origin/<head>``;
4. local branches with commits not in ``origin/main``, with the count.

Not listed: a local branch fully contained in ``origin/main``; a branch whose
pull request is merged or closed (a squash merge leaves the commits "ahead"
in git terms, so git alone cannot tell); anything the user hid. The branch
the checkout is on is listed whatever the rule says, so the label of the
dropdown never names something absent from it.

Semantics of the preferences: per user (under ``$HOME``, never in the
repository) and per repository (keyed by the git common directory, so all
worktrees of one repository share one list). ``hidden`` only stops an entry
from being shown; the ref stays in git. ``main`` cannot be hidden.

Without ``gh`` -- not installed, not logged in, no network, a timeout -- the
list is built from git alone and every entry has ``pr=None``.

Rejected: listing every local branch (352 on the machine this was written
on, 340 of them merged), listing every remote branch (the shared feature
branches of every contributor), and storing the pins in the repository
(each person's list is their own and must not show up in ``git status``).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

BASE = "origin/main"

#: A ref a person may type: no leading dash (an option), no whitespace, no
#: shell or revision-range syntax beyond what a tag or a sha needs.
_REF_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._/@+-]{0,200}$")

_GH_TTL_S = 120.0
_gh_cache: dict[str, tuple[float, Optional[dict]]] = {}
_gh_lock = threading.Lock()


@dataclass(frozen=True)
class PullRequest:
    number: int
    head: str
    state: str          # OPEN | MERGED | CLOSED
    ci: str             # passing | failing | pending | none
    url: str = ""


@dataclass
class BranchEntry:
    ref: str                     # the dropdown value: what git is asked for
    kind: str                    # main | local | remote | tag | commit
    ahead: Optional[int] = None  # commits not in origin/main
    pr: Optional[int] = None
    ci: str = ""
    pinned: bool = False
    current: bool = False
    notes: list = field(default_factory=list)

    @property
    def label(self) -> str:
        parts = [self.ref]
        if self.ahead:
            parts.append(f"+{self.ahead}")
        if self.pr is not None:
            parts.append(f"PR #{self.pr}" + (f" {self.ci}" if self.ci else ""))
        if self.pinned:
            parts.append("pinned")
        if self.current:
            parts.append("checked out here")
        return parts[0] + ("  (" + ", ".join(parts[1:]) + ")"
                           if len(parts) > 1 else "")


# -- git ---------------------------------------------------------------------

def _git(repo, *args: str, timeout: float = 30.0) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env.setdefault("GIT_OPTIONAL_LOCKS", "0")
    env["GIT_TERMINAL_PROMPT"] = "0"
    try:
        return subprocess.run(
            ["git", "-C", str(repo), *args], capture_output=True, text=True,
            timeout=timeout, env=env, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return subprocess.CompletedProcess(args, 1, "", str(exc))


def _out(repo, *args: str, timeout: float = 30.0) -> str:
    r = _git(repo, *args, timeout=timeout)
    return r.stdout.strip() if r.returncode == 0 else ""


def valid_ref(text: str) -> bool:
    """True for something that can only be read as a ref name or a sha."""
    text = text or ""
    return bool(_REF_RE.match(text)) and ".." not in text and not text.endswith(".lock")


def base_ref(repo) -> str:
    """``origin/main`` when it exists, else ``main``."""
    if _git(repo, "rev-parse", "--verify", "--quiet",
            f"refs/remotes/{BASE}").returncode == 0:
        return BASE
    return "main"


def repo_key(repo) -> str:
    """The git common directory, resolved: one key for all worktrees."""
    common = _out(repo, "rev-parse", "--path-format=absolute",
                  "--git-common-dir")
    if not common:
        try:
            return str(Path(repo).resolve())
        except OSError:
            return str(repo)
    try:
        return str(Path(common).resolve())
    except OSError:
        return common


def ref_kind(repo, ref: str) -> str:
    """local | remote | tag | commit, or "" when *ref* does not resolve."""
    for prefix, kind in (("refs/heads/", "local"), ("refs/remotes/", "remote"),
                         ("refs/tags/", "tag")):
        if _git(repo, "show-ref", "--verify", "--quiet",
                prefix + ref).returncode == 0:
            return kind
    if re.fullmatch(r"[0-9a-fA-F]{4,40}", ref or "") and _git(
            repo, "rev-parse", "--verify", "--quiet",
            f"{ref}^{{commit}}").returncode == 0:
        return "commit"
    return ""


def ahead_of_base(repo, ref: str, base: Optional[str] = None) -> Optional[int]:
    base = base or base_ref(repo)
    raw = _out(repo, "rev-list", "--count", f"{base}..{ref}")
    try:
        return int(raw)
    except ValueError:
        return None


def current_branch(repo) -> str:
    return _out(repo, "branch", "--show-current")


def unmerged_local_branches(repo, base: Optional[str] = None) -> list[str]:
    """Local branches with commits not in *base*, newest commit first."""
    base = base or base_ref(repo)
    raw = _out(repo, "for-each-ref", "--sort=-committerdate",
               f"--no-merged={base}", "--format=%(refname:short)",
               "refs/heads")
    return [line.strip() for line in raw.splitlines() if line.strip()]


def merged_local_branches(repo, base: Optional[str] = None) -> list[str]:
    base = base or base_ref(repo)
    raw = _out(repo, "for-each-ref", f"--merged={base}",
               "--format=%(refname:short)", "refs/heads")
    return [line.strip() for line in raw.splitlines() if line.strip()]


# -- gh ----------------------------------------------------------------------

def _ci_state(rollup) -> str:
    if not rollup:
        return "none"
    states = []
    for item in rollup:
        if not isinstance(item, dict):
            continue
        value = str(item.get("conclusion") or item.get("state")
                    or item.get("status") or "").upper()
        states.append(value)
    if any(s in ("FAILURE", "ERROR", "CANCELLED", "TIMED_OUT",
                 "ACTION_REQUIRED", "STARTUP_FAILURE") for s in states):
        return "failing"
    if any(s in ("", "PENDING", "QUEUED", "IN_PROGRESS", "WAITING",
                 "EXPECTED", "REQUESTED") for s in states):
        return "pending"
    return "passing"


def default_gh() -> Optional[str]:
    return shutil.which("gh")


def pull_requests(repo, gh: Optional[str] = "auto", *, refresh: bool = False,
                  timeout: float = 20.0) -> Optional[dict]:
    """``{head branch: PullRequest}`` (newest pull request per head wins), or
    None when gh could not answer. Cached for two minutes per repository."""
    if gh == "auto":
        gh = default_gh()
    if not gh:
        return None
    key = f"{repo_key(repo)}|{gh}"
    with _gh_lock:
        hit = _gh_cache.get(key)
        if hit and not refresh and time.monotonic() - hit[0] < _GH_TTL_S:
            return hit[1]
    env = dict(os.environ)
    env["GH_PROMPT_DISABLED"] = "1"
    env["GIT_TERMINAL_PROMPT"] = "0"
    env.setdefault("GH_NO_UPDATE_NOTIFIER", "1")
    result: Optional[dict] = None
    try:
        proc = subprocess.run(
            [gh, "pr", "list", "--state", "all", "--limit", "300", "--json",
             "number,headRefName,state,statusCheckRollup,url"],
            cwd=str(repo), capture_output=True, text=True, timeout=timeout,
            env=env, check=False)
        if proc.returncode == 0:
            rows = json.loads(proc.stdout or "[]")
            result = {}
            for row in rows if isinstance(rows, list) else []:
                try:
                    pr = PullRequest(
                        number=int(row.get("number")),
                        head=str(row.get("headRefName") or ""),
                        state=str(row.get("state") or "").upper(),
                        ci=_ci_state(row.get("statusCheckRollup")),
                        url=str(row.get("url") or ""))
                except (TypeError, ValueError):
                    continue
                if not pr.head:
                    continue
                old = result.get(pr.head)
                # An open one wins; otherwise the newest number.
                if (old is None or (pr.state == "OPEN") > (old.state == "OPEN")
                        or ((pr.state == "OPEN") == (old.state == "OPEN")
                            and pr.number > old.number)):
                    result[pr.head] = pr
    except (OSError, subprocess.TimeoutExpired, ValueError):
        result = None
    with _gh_lock:
        _gh_cache[key] = (time.monotonic(), result)
    return result


# -- the user's own list -----------------------------------------------------

def _prefs_path() -> Path:
    """Resolved per call, so a redirected home moves it."""
    return Path.home() / ".delfin" / "branch_list.json"


def _load_all() -> dict:
    try:
        data = json.loads(_prefs_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _save_all(data: dict) -> None:
    from delfin.agent.state_paths import write_text_atomic
    write_text_atomic(_prefs_path(), json.dumps(data, indent=1, sort_keys=True))


def load_prefs(repo) -> dict:
    """``{"pinned": [...], "hidden": [...]}`` for *repo*."""
    row = (_load_all().get("repos") or {}).get(repo_key(repo)) or {}
    return {"pinned": [str(x) for x in row.get("pinned") or []],
            "hidden": [str(x) for x in row.get("hidden") or []]}


def _update(repo, *, add: tuple = (), remove: tuple = ()) -> dict:
    data = _load_all()
    repos = data.setdefault("repos", {})
    row = repos.setdefault(repo_key(repo), {"pinned": [], "hidden": []})
    for kind, ref in remove:
        row[kind] = [x for x in row.get(kind) or [] if x != ref]
    for kind, ref in add:
        lst = row.setdefault(kind, [])
        if ref not in lst:
            lst.append(ref)
    _save_all(data)
    return row


def pin(repo, ref: str) -> None:
    _update(repo, add=(("pinned", ref),), remove=(("hidden", ref),))


def unpin(repo, ref: str) -> None:
    _update(repo, remove=(("pinned", ref),))


def hide(repo, ref: str) -> tuple[bool, str]:
    """Stop showing *ref*. Reversible with :func:`unhide`; git is untouched."""
    if ref == "main":
        return False, "main is always listed."
    _update(repo, add=(("hidden", ref),), remove=(("pinned", ref),))
    return True, f"{ref} is hidden from your list (still in git; 'Show again' brings it back)."


def unhide(repo, ref: str) -> None:
    _update(repo, remove=(("hidden", ref),))


def forget(repo, ref: str) -> None:
    """Drop *ref* from both lists (after it was deleted)."""
    _update(repo, remove=(("pinned", ref), ("hidden", ref)))


def fetch_remote_branch(repo, name: str, remote: str = "origin") -> bool:
    """Fetch ``<remote>/<name>`` into its remote-tracking ref. No local
    branch is created and the checkout is not touched."""
    if not valid_ref(name):
        return False
    r = _git(repo, "fetch", "--quiet", "--no-tags", remote,
             f"+refs/heads/{name}:refs/remotes/{remote}/{name}", timeout=90.0)
    return r.returncode == 0


def add_to_list(repo, text: str) -> tuple[bool, str, str]:
    """Pin a ref the automatic rule hides. ``(ok, message, ref)``.

    A local branch, remote-tracking ref, tag or commit is pinned as typed.
    A name only the remote has is fetched first (``origin/<name>``)."""
    name = (text or "").strip()
    if name.startswith("refs/heads/"):
        name = name[len("refs/heads/"):]
    if not valid_ref(name):
        return False, f"Not a branch, tag or commit name: {name!r}", ""
    kind = ref_kind(repo, name)
    ref = name
    if not kind:
        remote_name = name[len("origin/"):] if name.startswith("origin/") else name
        if fetch_remote_branch(repo, remote_name):
            ref = f"origin/{remote_name}"
            kind = "remote"
    if not kind:
        return False, f"{name}: not found locally or on origin.", ""
    pin(repo, ref)
    return True, f"{ref} ({kind}) added to your list.", ref


def create_branch_from_main(repo, name: str) -> tuple[bool, str, str]:
    """Create *name* at ``origin/main`` without switching the checkout.

    ``origin/main`` is fetched first when the network allows; the branch is
    created with ``--no-track`` so a later push goes to its own name."""
    name = (name or "").strip()
    check = _git(repo, "check-ref-format", "--branch", name)
    if not valid_ref(name) or check.returncode != 0:
        return False, f"Not a valid branch name: {name!r}", ""
    if _git(repo, "show-ref", "--verify", "--quiet",
            f"refs/heads/{name}").returncode == 0:
        return False, f"Branch {name} already exists.", ""
    _git(repo, "fetch", "--quiet", "origin", "main", timeout=90.0)
    base = base_ref(repo)
    r = _git(repo, "branch", "--no-track", name, base)
    if r.returncode != 0:
        return False, (r.stderr or r.stdout).strip() or "git branch failed", ""
    pin(repo, name)
    sha = _out(repo, "rev-parse", "--short", name)
    return True, (f"Created {name} at {base} ({sha}); the shared checkout "
                  "stays where it is."), name


# -- the entries -------------------------------------------------------------

def branch_entries(repo, gh: Optional[str] = "auto", *,
                   refresh: bool = False) -> list[BranchEntry]:
    """The dropdown rows for *repo* (see the module docstring)."""
    main = BranchEntry(ref="main", kind="main")
    try:
        if not repo or not Path(repo).exists():
            return [main]
    except OSError:
        return [main]
    if _git(repo, "rev-parse", "--git-dir").returncode != 0:
        return [main]
    base = base_ref(repo)
    prefs = load_prefs(repo)
    hidden = set(prefs["hidden"])
    prs = pull_requests(repo, gh, refresh=refresh) or {}
    current = current_branch(repo)
    entries: list[BranchEntry] = [main]
    seen = {"main"}
    main.current = current == "main"

    def _add(entry: BranchEntry) -> None:
        if entry.ref in seen or (entry.ref in hidden and not entry.current):
            return
        seen.add(entry.ref)
        pr = prs.get(entry.ref) or (prs.get(entry.ref[len("origin/"):])
                                    if entry.ref.startswith("origin/") else None)
        if pr is not None and entry.pr is None:
            entry.pr, entry.ci = pr.number, pr.ci
        entries.append(entry)

    for ref in prefs["pinned"]:
        kind = ref_kind(repo, ref)
        if not kind:
            continue
        _add(BranchEntry(ref=ref, kind=kind, pinned=True,
                         ahead=ahead_of_base(repo, ref, base) or None,
                         current=ref == current))
    local = set(_out(repo, "for-each-ref", "--format=%(refname:short)",
                     "refs/heads").splitlines())
    for head, pr in sorted(prs.items(), key=lambda kv: -kv[1].number):
        if pr.state != "OPEN":
            continue
        if head in local:
            ref, kind = head, "local"
        elif ref_kind(repo, f"origin/{head}") == "remote":
            ref, kind = f"origin/{head}", "remote"
        else:
            ref, kind = f"origin/{head}", "remote"
        _add(BranchEntry(ref=ref, kind=kind, pr=pr.number, ci=pr.ci,
                         ahead=ahead_of_base(repo, ref, base) or None,
                         current=ref == current))
    for name in unmerged_local_branches(repo, base):
        if name == "main":
            continue
        pr = prs.get(name)
        if pr is not None and pr.state in ("MERGED", "CLOSED") and name != current:
            continue
        _add(BranchEntry(ref=name, kind="local",
                         ahead=ahead_of_base(repo, name, base),
                         current=name == current))
    if current and current not in seen:
        _add(BranchEntry(ref=current, kind="local", current=True,
                         ahead=ahead_of_base(repo, current, base) or None))
    return entries


def filter_entries(entries: list[BranchEntry], text: str) -> list[BranchEntry]:
    """Entries whose ref contains *text* (case-insensitive); main stays."""
    needle = (text or "").strip().lower()
    if not needle:
        return list(entries)
    return [e for e in entries if e.ref == "main" or needle in e.ref.lower()]


def hidden_entries(repo) -> list[str]:
    return list(load_prefs(repo)["hidden"])
