"""Session-persistent memory of operator refusals.

``terminal_confirm.py`` keeps ``last_refusal_reason`` -- ONE reason,
cleared by the next dialog ("a reason never outlives its refusal").
That is right for display, but leaves the session without a memory:
two turns after a refusal it no longer knows that, and why, something
was refused -- least of all after compaction. The night of 26 Sep 2026
a session asked three times for the same file outside its workspace,
each time worded differently, and each question woke the operator.

This module is that memory. It stores, per session, what a refusal was
ABOUT -- the path or directory for a file tool, the program plus target
for bash -- never the wording of the command. ``matches`` recognizes
the same request in a different spelling (``cat X``, ``sed -n 1,50p X``,
``read_file X``, ``head X`` -> the same target X), conservatively: only
unambiguous hits.

It never decides "allow". Its only effect is that a repeated question
can be AVOIDED by handing the session the old refusal with its reason;
the human may change their mind at any time and grant something
themselves. The refusal memory is a reminder, not a gate.
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

# Refusals are small; a handful covers a session. Bounded so a long
# session cannot grow this without limit (and the block it feeds stays
# small after compaction).
MAX_ENTRIES = 16

# Programs that read one file given as an argument. Used to pull the
# target out of a bash command; anything not listed yields no match --
# that is the conservative direction.
_READER_PROGRAMS = (
    "cat", "head", "tail", "sed", "grep", "less", "more", "nl",
    "awk", "sort", "wc", "md5sum", "sha256sum", "cp", "diff",
    "python", "python3", "jq",
)


class Refusal:
    """One recorded refusal: what it was about, why, when."""

    __slots__ = ("tool", "target", "reason", "time", "is_dir")

    def __init__(self, tool: str, target: str, reason: str, time: str,
                 is_dir: bool = False) -> None:
        self.tool = str(tool or "")
        self.target = str(target or "")
        self.reason = str(reason or "")
        self.time = str(time or "")
        self.is_dir = bool(is_dir)

    def to_dict(self) -> dict:
        return {"tool": self.tool, "target": self.target,
                "reason": self.reason, "time": self.time,
                "is_dir": self.is_dir}

    @classmethod
    def from_dict(cls, d: dict) -> "Refusal":
        return cls(
            tool=d.get("tool", ""),
            target=d.get("target", ""),
            reason=d.get("reason", ""),
            time=d.get("time", ""),
            is_dir=bool(d.get("is_dir", False)),
        )


def _norm(p: str) -> str:
    """Canonical form of a path for comparison: absolute, no redundant
    separators or dot segments. ``.``-relative stays distinguishable."""
    if not p:
        return ""
    if p.startswith("./"):
        p = p[2:]
    if p == ".":
        return "."
    p = re.sub(r"/+", "/", p)
    if not p.startswith(("/", "~", ".")):
        p = "./" + p
    parts = []
    for seg in p.split("/"):
        if seg in ("", "."):
            continue
        if seg == ".." and parts and parts[-1] not in ("", ".."):
            parts.pop()
        else:
            parts.append(seg)
    lead = "/" if p.startswith("/") else ""
    return lead + "/".join(parts) or ("." if not p.startswith("/") else "/")


def _contains(directory: str, path: str) -> bool:
    """True when ``path`` lies inside ``directory``."""
    d = _norm(directory).rstrip("/") + "/"
    return _norm(path).startswith(d)


def _extract_read_target(command: str) -> str | None:
    """The file a reader command is about, or None when unclear.

    Splits off a leading env/pipeline prefix, takes the program name,
    skips the program's flag arguments, and returns the first
    non-flag argument only when the program is a known reader AND has
    exactly one such argument left -- two or more is ambiguous (a
    compound ``cat a > b`` is not the same ask) and matches nothing.
    """
    text = command.strip()
    if not text:
        return None
    # Drop a leading env assignment prefix (FOO=1 cmd ...) and quotes
    # around the whole command.
    while re.match(r"^\w+=\S*\s+", text):
        text = text.split(None, 1)[1]
    tokens = text.split()
    if not tokens:
        return None
    prog = tokens[0]
    prog = prog.rsplit("/", 1)[-1]
    if prog not in _READER_PROGRAMS:
        return None
    # A sed range/address argument ("1,50p", "$a", "2q") is not a file;
    # tokens starting with a digit or $ are skipped. A "/" would clash
    # with absolute paths, so slash-addresses are not recognized.
    addr = re.compile(r"^[\d$]")
    rest = [t for t in tokens[1:]
            if not t.startswith("-") and not addr.match(t)]
    if len(rest) != 1:
        return None
    return rest[0]


def _targets_of(tool: str, args: dict) -> list[str]:
    """All paths the call is about, or [] when unclear."""
    if not isinstance(args, dict):
        return []
    if tool in ("read_file", "write_file", "edit_file", "read_document"):
        for key in ("path", "file"):
            v = args.get(key)
            if isinstance(v, str) and v.strip():
                return [v]
    if tool == "bash":
        t = _extract_read_target(str(args.get("command", "")))
        return [t] if t else []
    return []


class RefusalMemory:
    """Bounded, serializable per-session record of refusals."""

    MAX_ENTRIES = MAX_ENTRIES

    def __init__(self, store: Path | str | None = None,
                 entries: list[Refusal] | None = None) -> None:
        self.store = Path(store) if store else None
        self.entries: list[Refusal] = list(entries or [])

    # -- recording --------------------------------------------------
    def record(self, refusal: Refusal) -> None:
        """Add a refusal; a same-tool-same-target entry is replaced, the
        list stays newest-last and bounded to MAX_ENTRIES."""
        keep = [r for r in self.entries
                if not (r.tool == refusal.tool
                        and _norm(r.target) == _norm(refusal.target))]
        keep.append(refusal)
        self.entries = keep[-self.MAX_ENTRIES:]
        self._save()

    def _save(self) -> None:
        if self.store is None:
            return
        try:
            self.store.parent.mkdir(parents=True, exist_ok=True)
            tmp = self.store.with_suffix(".tmp")
            tmp.write_text(json.dumps(self.to_dict()), encoding="utf-8")
            os.replace(tmp, self.store)
        except Exception:
            # Memory is a convenience: a store that cannot be written
            # must never fail the refusal it records.
            pass

    # -- lookup -----------------------------------------------------
    def matches(self, tool: str, args: dict) -> Refusal | None:
        """The refusal this call repeats, or None.

        Conservative: only calls with exactly one unambiguous target
        match. A refused directory covers paths beneath it; a refused
        file covers only itself. The tool does not have to match --
        asking for the same file with ``cat`` after a refused
        ``read_file`` is the same ask, just spelled differently.
        """
        targets = _targets_of(tool, args)
        if not targets:
            return None
        for target in targets:
            norm = _norm(target)
            if not norm:
                continue
            for r in reversed(self.entries):
                if r.is_dir:
                    if _contains(r.target, target):
                        return r
                elif _norm(r.target) == norm:
                    return r
        return None

    def note(self, refusal: Refusal) -> str:
        """The English reminder handed back with a repeated ask."""
        time = refusal.time or "(time unknown)"
        reason = refusal.reason or "(no reason given)"
        hint = (""
                if refusal.is_dir
                else " The refusal suggests asking for what you need "
                     "inside your own workspace, or asking the "
                     "operator to grant that specific path.")
        return (f"The operator already refused this at {time}: "
                f"'{reason}'. Do not ask again; the human may still "
                f"grant it themselves if they choose.{hint}")

    # -- serialization ----------------------------------------------
    def to_dict(self) -> dict:
        return {"entries": [r.to_dict() for r in self.entries]}

    @classmethod
    def from_dict(cls, d: dict) -> "RefusalMemory":
        entries = [Refusal.from_dict(e)
                   for e in (d or {}).get("entries", [])
                   if isinstance(e, dict)]
        return cls(entries=entries[-cls.MAX_ENTRIES:])

    @classmethod
    def load(cls, store: Path | str) -> "RefusalMemory":
        try:
            return cls.from_dict(
                json.loads(Path(store).read_text(encoding="utf-8")))
        except Exception:
            return cls()
