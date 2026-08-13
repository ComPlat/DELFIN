"""Per-directory trust for executable definitions a WORKSPACE supplies.

A repository the user merely CHECKS OUT must not be able to run code.
Two kinds of file inside a workspace are executable configuration:

``<ws>/.delfin/settings.json`` and ``<ws>/.delfin/settings.local.json``
    hook commands. Run through a shell with the parent environment, on
    every prompt and every tool call, outside the permission gate and
    outside filesystem isolation.

``<ws>/.delfin/mcp_servers.json``
    MCP server definitions. Spawned with the parent environment while
    the tool surface is being ASSEMBLED -- before the model emits a
    token and before any consent -- and then answer every call routed
    to them.

Before this module the only question asked of such a file was "is this
a REGISTERED OFFICE folder?". Every other directory was honoured, which
includes one cloned five seconds ago. Three earlier fixes closed single
instances of that class (a project MCP entry hijacking a builtin, a
locked scope contributing hooks, a workspace status line running a
command). The class stays open for as long as the question is "what
kind of folder is this?" instead of "did the user say they trust this
one?".

This module asks the second question, once, in one place:

* **Untrusted is the default, and it is not silent.** ``gate()`` reports
  how many definitions of which kind were withheld, from which file,
  and how to trust them. A guard the user cannot see is one they will
  work around.
* **Only the user grants trust.** ``trust_workspace`` refuses any actor
  but ``"user"``. Nothing in the tool path may call it -- a model that
  can grant itself trust has not been gated, only delayed.
* **Per directory, no inheritance.** The key is the RESOLVED absolute
  path and nothing else. A trusted parent does not vouch for a child, a
  child does not vouch for a parent, a sibling vouches for nobody. A
  monorepo of vendored dependencies is the case that makes inheritance
  wrong: trusting the checkout would trust every subdirectory in it.
* **Content-bound.** What was trusted is a digest of the definitions,
  not merely the path. A directory trusted last week whose hook command
  changed today is a NEW decision -- ``git pull`` is precisely how the
  payload arrives.
* **Recorded outside the workspace.** A workspace cannot vouch for
  itself, so the record lives in ``~/.delfin/trusted_workspaces.json``.

Symlinks resolve, deliberately: trust is a statement about a directory
on disk, and a link is a second name for one. Both the grant and the
check resolve, so they agree; a moved directory therefore LOSES its
trust, because the only handle we have on identity is the path, and a
directory that moved is one the user should look at again.

A settings file inside the workspace that is itself a link OUT of the
workspace is refused outright, trusted or not: what the user reviewed
would not be what the loader reads.

Adding a third kind
-------------------
Register it here. Loaders must ask ``gate()`` for the paths they may
read rather than building their own -- that is what stops a fourth kind
from forgetting to ask, and ``tests/test_a_checked_out_repository_
cannot_run_commands.py`` fails when a loader reads a path this module
does not know about.
"""

from __future__ import annotations

import hashlib
import json
import os
import threading
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable


KIND_HOOKS = "hooks"
KIND_MCP_SERVERS = "mcp_servers"

_STORE_VERSION = 1

# Who may grant trust. A single value, checked rather than documented:
# every other caller -- a tool, a hook, a model-emitted command -- is a
# path by which the thing being gated grants itself the gate.
ACTOR_USER = "user"


# ---------------------------------------------------------------------------
# Kinds
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Kind:
    """One family of executable definitions a workspace can supply.

    ``relative_paths`` is the whole answer to "which files inside a
    workspace does this kind live in": the loader asks for them instead
    of spelling them out, so a new file cannot be added at the loader
    and forgotten at the guard -- which is how ``settings.local.json``
    came to be read by the hook loader and missing from the protected
    globs.
    """

    name: str
    relative_paths: tuple[str, ...]
    # Singular noun for the notice, e.g. "hook command".
    noun: str
    # What the user types to grant trust for this kind.
    grant_command: str
    # The executable part of one parsed settings file. Everything else
    # (formatting, unrelated keys, a ``_meta`` timestamp) is excluded so
    # trust survives edits that cannot execute anything.
    extract: Callable[[dict], Any]
    # How many definitions the extracted part contains.
    count: Callable[[Any], int]


def _extract_hooks(raw: dict) -> Any:
    obj = raw.get("hooks")
    return obj if isinstance(obj, dict) else None


def _count_hooks(extracted: Any) -> int:
    if not isinstance(extracted, dict):
        return 0
    total = 0
    for entries in extracted.values():
        if not isinstance(entries, list):
            continue
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            cmds = entry.get("hooks")
            if isinstance(cmds, list):
                total += sum(1 for c in cmds if isinstance(c, dict))
    return total


def _extract_servers(raw: dict) -> Any:
    obj = raw.get("servers")
    return obj if isinstance(obj, dict) else None


def _count_servers(extracted: Any) -> int:
    if not isinstance(extracted, dict):
        return 0
    return sum(1 for cfg in extracted.values()
               if isinstance(cfg, dict) and cfg.get("enabled", True))


_KINDS: dict[str, Kind] = {}
_KINDS_LOCK = threading.Lock()


def register_kind(kind: Kind) -> None:
    with _KINDS_LOCK:
        _KINDS[kind.name] = kind


def kinds() -> tuple[Kind, ...]:
    with _KINDS_LOCK:
        return tuple(_KINDS[name] for name in sorted(_KINDS))


def get_kind(name: str) -> Kind:
    with _KINDS_LOCK:
        try:
            return _KINDS[name]
        except KeyError:
            raise KeyError(
                f"unknown workspace-trust kind {name!r}; "
                f"register it in workspace_trust"
            ) from None


register_kind(Kind(
    name=KIND_HOOKS,
    relative_paths=(".delfin/settings.json", ".delfin/settings.local.json"),
    noun="hook command",
    grant_command="/hooks trust",
    extract=_extract_hooks,
    count=_count_hooks,
))
register_kind(Kind(
    name=KIND_MCP_SERVERS,
    relative_paths=(".delfin/mcp_servers.json",),
    noun="MCP server",
    grant_command="/mcp trust",
    extract=_extract_servers,
    count=_count_servers,
))


def all_relative_paths() -> tuple[str, ...]:
    """Every workspace-relative file any registered kind reads.

    The protected-globs test asserts containment against this, so a kind
    added without protecting its file fails there rather than in the
    field.
    """
    out: list[str] = []
    for kind in kinds():
        for rel in kind.relative_paths:
            if rel not in out:
                out.append(rel)
    return tuple(out)


# ---------------------------------------------------------------------------
# The store
# ---------------------------------------------------------------------------

def _trust_store_path() -> Path:
    """Where the decisions live. OUTSIDE any workspace, on purpose."""
    return Path.home() / ".delfin" / "trusted_workspaces.json"


_STORE_LOCK = threading.Lock()
# (path, mtime_ns, size) -> parsed store. load_hooks runs on every tool
# call, so the store is read on every tool call; this keeps that to a
# stat once the file is warm.
_STORE_CACHE: dict[tuple, dict] = {}


def _read_store() -> dict:
    path = _trust_store_path()
    try:
        st = path.stat()
        key = (str(path), st.st_mtime_ns, st.st_size)
    except OSError:
        return {"version": _STORE_VERSION, "workspaces": {}}
    with _STORE_LOCK:
        cached = _STORE_CACHE.get(key)
    if cached is not None:
        return cached
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        data = {}
    if not isinstance(data, dict):
        data = {}
    workspaces = data.get("workspaces")
    store = {
        "version": data.get("version", _STORE_VERSION),
        "workspaces": workspaces if isinstance(workspaces, dict) else {},
    }
    with _STORE_LOCK:
        _STORE_CACHE.clear()
        _STORE_CACHE[key] = store
    return store


def _write_store(store: dict) -> None:
    path = _trust_store_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(store, indent=2, ensure_ascii=False,
                              sort_keys=True), encoding="utf-8")
    tmp.replace(path)
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass
    with _STORE_LOCK:
        _STORE_CACHE.clear()


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _canonical(workspace: Path | str) -> Path | None:
    """The resolved absolute path a decision is keyed by, or None.

    ``resolve()`` follows symlinks: a link is a second name for a
    directory, not a second directory, and trust is a statement about
    the one on disk. Grant and check resolve identically, so they agree.
    """
    try:
        return Path(workspace).expanduser().resolve()
    except (OSError, RuntimeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# What a workspace offers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Offer:
    """One settings file a workspace ships that would execute if honoured."""

    kind: str
    path: Path              # absolute
    relative: str           # workspace-relative, as the user sees it
    count: int              # definitions inside it
    escapes: bool = False   # a link pointing out of the workspace


def _read_json(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _escapes_workspace(path: Path, root: Path) -> bool:
    """True when *path* resolves outside *root*.

    A workspace may ship ``.delfin`` as a link to somewhere else. What
    the user reviewed would then not be what the loader reads, so such a
    file is refused whether or not the directory is trusted.
    """
    try:
        real = path.resolve()
    except (OSError, RuntimeError, ValueError):
        return True
    try:
        real.relative_to(root)
        return False
    except ValueError:
        return True


def offers(kind_name: str, workspace: Path | str | None) -> list[Offer]:
    """Every file of *kind_name* that exists in *workspace*.

    Reads content to count, never to honour. Safe to call on a workspace
    nobody trusts -- that is the point.
    """
    root = _canonical(workspace) if workspace is not None else None
    if root is None:
        return []
    kind = get_kind(kind_name)
    out: list[Offer] = []
    for rel in kind.relative_paths:
        path = root / rel
        try:
            if not path.is_file():
                continue
        except OSError:
            continue
        if _escapes_workspace(path, root):
            out.append(Offer(kind=kind_name, path=path, relative=rel,
                             count=0, escapes=True))
            continue
        count = kind.count(kind.extract(_read_json(path)))
        out.append(Offer(kind=kind_name, path=path, relative=rel,
                         count=count))
    return out


def digest(kind_name: str, workspace: Path | str | None) -> str:
    """A stable fingerprint of what *workspace* offers for *kind_name*.

    Over the executable part only, canonicalised: reformatting the file
    or touching an unrelated key does not re-ask, changing a command
    does. Empty string when nothing is offered.
    """
    items = [o for o in offers(kind_name, workspace) if not o.escapes]
    if not items:
        return ""
    kind = get_kind(kind_name)
    h = hashlib.sha256()
    for offer in items:
        extracted = kind.extract(_read_json(offer.path))
        try:
            canon = json.dumps(extracted, sort_keys=True,
                               ensure_ascii=False, separators=(",", ":"))
        except (TypeError, ValueError):
            canon = repr(extracted)
        h.update(offer.relative.encode("utf-8"))
        h.update(b"\x00")
        h.update(canon.encode("utf-8"))
        h.update(b"\x00")
    return h.hexdigest()


# ---------------------------------------------------------------------------
# The decision
# ---------------------------------------------------------------------------

# Why a gate answered the way it did.
STATE_NO_WORKSPACE = "no-workspace"
STATE_NOTHING_OFFERED = "nothing-offered"
STATE_TRUSTED = "trusted"
STATE_UNTRUSTED = "untrusted"
STATE_CHANGED = "changed"


@dataclass
class Decision:
    """What a loader may read, what it may not, and what to tell the user."""

    kind: str
    workspace: Path | None
    state: str
    paths: list[Path] = field(default_factory=list)      # honoured
    withheld: list[Offer] = field(default_factory=list)  # refused
    notice: str = ""
    # The same fact in one clause. ``make_record`` truncates a reason at
    # 300 characters, and the actionable part of the full notice -- what
    # to type -- is at the END of it, so the log kept the explanation and
    # cut the instruction.
    short_notice: str = ""

    @property
    def trusted(self) -> bool:
        return self.state == STATE_TRUSTED

    @property
    def withheld_count(self) -> int:
        return sum(o.count for o in self.withheld)


def _plural(n: int, noun: str) -> str:
    return f"{n} {noun}" + ("" if n == 1 else "s")


def _notice_for(kind: Kind, workspace: Path, state: str,
                withheld: list[Offer], trusted_at: str) -> str:
    escaping = [o for o in withheld if o.escapes]
    real = [o for o in withheld if not o.escapes]
    parts: list[str] = []
    if real:
        where = ", ".join(o.relative for o in real)
        total = sum(o.count for o in real)
        was = "was" if total == 1 else "were"
        if state == STATE_CHANGED:
            parts.append(
                f"{_plural(total, kind.noun)} in {where} under {workspace} "
                f"{was} NOT loaded: you trusted this directory on "
                f"{trusted_at or 'an earlier date'} and the definitions have "
                f"changed since. Review the file(s) and re-trust with "
                f"`{kind.grant_command}`."
            )
        else:
            parts.append(
                f"{_plural(total, kind.noun)} in {where} under {workspace} "
                f"{was} NOT loaded: this directory is not trusted. Untrusted "
                f"is the default -- a checked-out repository must not be able "
                f"to run commands. Review the file(s) and run "
                f"`{kind.grant_command}` to load them."
            )
    for o in escaping:
        parts.append(
            f"{o.relative} under {workspace} is a link out of the workspace "
            f"and was NOT read: what you would review is not what would run."
        )
    return " ".join(parts)


def _short_notice_for(kind: Kind, state: str, withheld: list[Offer]) -> str:
    total = sum(o.count for o in withheld if not o.escapes)
    links = sum(1 for o in withheld if o.escapes)
    bits: list[str] = []
    if total:
        why = ("the definitions changed since you trusted them"
               if state == STATE_CHANGED else "the directory is not trusted")
        bits.append(f"{_plural(total, kind.noun)} withheld — {why}; "
                    f"run `{kind.grant_command}`")
    if links:
        bits.append(f"{links} settings file(s) link out of the workspace")
    return "; ".join(bits)


# One announcement per (kind, workspace, digest). load_hooks runs on
# every tool call, so recording each refusal would bury the log it is
# meant to inform -- and a line the user scrolls past is as good as no
# line at all.
_ANNOUNCED: set[tuple[str, str, str]] = set()
_ANNOUNCE_LOCK = threading.Lock()


def _announce(kind: Kind, workspace: Path, decision: "Decision",
              fingerprint: str) -> None:
    key = (kind.name, str(workspace), fingerprint)
    with _ANNOUNCE_LOCK:
        if key in _ANNOUNCED:
            return
        _ANNOUNCED.add(key)
    try:
        from . import audit_log as _audit
        _audit.append(_audit.make_record(
            tool="workspace_trust",
            decision="denied",
            mode="trust",
            path=str(workspace),
            reason=decision.short_notice or decision.notice,
            extra={
                "trust_kind": kind.name,
                "withheld": [o.relative for o in decision.withheld],
                "definitions": decision.withheld_count,
                "notice": decision.notice,
            },
        ))
    except Exception:
        pass


def reset_announcements() -> None:
    """Forget what has already been reported. Used by tests."""
    with _ANNOUNCE_LOCK:
        _ANNOUNCED.clear()
    with _STORE_LOCK:
        _STORE_CACHE.clear()


def gate(kind_name: str, workspace: Path | str | None) -> Decision:
    """Which of *workspace*'s *kind_name* files a loader may read.

    The single question every loader asks. It never raises: a broken
    store, an unreadable directory or a malformed settings file all
    resolve to "not trusted", because the failure mode of a trust check
    has to be refusal.
    """
    kind = get_kind(kind_name)
    root = _canonical(workspace) if workspace is not None else None
    if root is None:
        return Decision(kind=kind_name, workspace=None,
                        state=STATE_NO_WORKSPACE)
    try:
        found = offers(kind_name, root)
    except Exception:
        found = []
    if not found:
        return Decision(kind=kind_name, workspace=root,
                        state=STATE_NOTHING_OFFERED)
    readable = [o for o in found if not o.escapes]
    escaping = [o for o in found if o.escapes]

    fingerprint = digest(kind_name, root)
    record = (_read_store().get("workspaces") or {}).get(str(root)) or {}
    granted = (record.get("kinds") or {}).get(kind_name) or {}
    trusted_at = str(granted.get("trusted_at", ""))

    if granted.get("digest") == fingerprint and fingerprint:
        decision = Decision(
            kind=kind_name, workspace=root, state=STATE_TRUSTED,
            paths=[o.path for o in readable], withheld=escaping,
        )
    else:
        state = STATE_CHANGED if granted else STATE_UNTRUSTED
        decision = Decision(kind=kind_name, workspace=root, state=state,
                            paths=[], withheld=found)
    if decision.withheld:
        decision.notice = _notice_for(kind, root, decision.state,
                                      decision.withheld, trusted_at)
        decision.short_notice = _short_notice_for(
            kind, decision.state, decision.withheld)
        _announce(kind, root, decision, fingerprint)
    return decision


def workspace_paths(kind_name: str,
                    workspace: Path | str | None) -> list[Path]:
    """The absolute paths a loader may read. Convenience over ``gate``."""
    return gate(kind_name, workspace).paths


def pending_notices(workspace: Path | str | None) -> list[str]:
    """Every kind's notice for *workspace*, for a caller that shows all."""
    out: list[str] = []
    for kind in kinds():
        note = gate(kind.name, workspace).notice
        if note:
            out.append(note)
    return out


# ---------------------------------------------------------------------------
# Granting and revoking -- the user, and nobody else
# ---------------------------------------------------------------------------

class TrustActorError(PermissionError):
    """Raised when something other than the user tries to grant trust."""


def _require_user(actor: str, verb: str) -> None:
    if actor != ACTOR_USER:
        raise TrustActorError(
            f"only the user may {verb} a workspace (actor={actor!r}). "
            "Trust is not something the agent, a tool or a hook can grant "
            "itself -- ask the user to run the trust command."
        )


def trust_workspace(
    workspace: Path | str,
    kind_names: list[str] | tuple[str, ...] | None = None,
    *,
    actor: str,
) -> dict:
    """Record that the USER trusts *workspace* for the given kinds.

    ``actor`` is mandatory and must be ``"user"``. It is a keyword with
    no default so that granting trust cannot be an incidental call: the
    only correct value has to be typed at the call site, and every call
    site is checked by
    ``tests/test_a_checked_out_repository_cannot_run_commands.py``.

    Records the resolved absolute path, which kinds were trusted, a
    digest of exactly what was trusted, and when. Returns the stored
    record.
    """
    _require_user(actor, "trust")
    root = _canonical(workspace)
    if root is None:
        raise ValueError(f"cannot resolve workspace {workspace!r}")
    if not root.is_dir():
        raise ValueError(f"not a directory: {root}")
    names = tuple(kind_names) if kind_names else tuple(
        k.name for k in kinds())
    for name in names:
        get_kind(name)                      # unknown kind -> KeyError, loudly

    store = _read_store()
    workspaces = dict(store.get("workspaces") or {})
    record = dict(workspaces.get(str(root)) or {})
    record.setdefault("first_trusted_at", _now_iso())
    record["granted_by"] = ACTOR_USER
    named = str(Path(workspace).expanduser())
    if named != str(root):
        # Kept for the listing only. The DECISION is keyed by the
        # resolved path; this just lets the user recognise what they
        # typed when the two differ.
        record["as_named"] = named
    granted = dict(record.get("kinds") or {})
    for name in names:
        found = offers(name, root)
        granted[name] = {
            "digest": digest(name, root),
            "files": [o.relative for o in found if not o.escapes],
            "definitions": sum(o.count for o in found if not o.escapes),
            "trusted_at": _now_iso(),
        }
    record["kinds"] = granted
    workspaces[str(root)] = record
    _write_store({"version": _STORE_VERSION, "workspaces": workspaces})
    reset_announcements()
    return {"workspace": str(root), **record}


def revoke_workspace(
    workspace: Path | str,
    kind_names: list[str] | tuple[str, ...] | None = None,
    *,
    actor: str,
) -> bool:
    """Withdraw trust. True when something was actually removed."""
    _require_user(actor, "untrust")
    root = _canonical(workspace)
    if root is None:
        return False
    store = _read_store()
    workspaces = dict(store.get("workspaces") or {})
    record = dict(workspaces.get(str(root)) or {})
    if not record:
        return False
    if not kind_names:
        workspaces.pop(str(root), None)
    else:
        granted = dict(record.get("kinds") or {})
        removed = [n for n in kind_names if granted.pop(n, None) is not None]
        if not removed:
            return False
        if granted:
            record["kinds"] = granted
            workspaces[str(root)] = record
        else:
            workspaces.pop(str(root), None)
    _write_store({"version": _STORE_VERSION, "workspaces": workspaces})
    reset_announcements()
    return True


def trusted_record(workspace: Path | str) -> dict | None:
    """What is on file for *workspace*, or None. Read-only."""
    root = _canonical(workspace)
    if root is None:
        return None
    record = (_read_store().get("workspaces") or {}).get(str(root))
    if not isinstance(record, dict):
        return None
    return {"workspace": str(root), **record}


def describe(workspace: Path | str | None) -> str:
    """One human-readable paragraph about *workspace*'s trust state."""
    root = _canonical(workspace) if workspace is not None else None
    if root is None:
        return "No workspace — only your own ~/.delfin configuration is read."
    lines: list[str] = [f"Workspace: {root}"]
    record = trusted_record(root) or {}
    granted = record.get("kinds") or {}
    for kind in kinds():
        found = offers(kind.name, root)
        total = sum(o.count for o in found if not o.escapes)
        state = gate(kind.name, root).state
        if not found:
            lines.append(f"  {kind.name}: none in this workspace")
            continue
        when = str((granted.get(kind.name) or {}).get("trusted_at", ""))
        suffix = f" (trusted {when})" if state == STATE_TRUSTED and when else ""
        lines.append(
            f"  {kind.name}: {_plural(total, kind.noun)} in "
            f"{', '.join(o.relative for o in found)} — {state}{suffix}"
        )
    return "\n".join(lines)


__all__ = [
    "ACTOR_USER",
    "KIND_HOOKS",
    "KIND_MCP_SERVERS",
    "Decision",
    "Kind",
    "Offer",
    "TrustActorError",
    "all_relative_paths",
    "describe",
    "digest",
    "gate",
    "get_kind",
    "kinds",
    "offers",
    "pending_notices",
    "register_kind",
    "reset_announcements",
    "revoke_workspace",
    "trust_workspace",
    "trusted_record",
    "workspace_paths",
]
