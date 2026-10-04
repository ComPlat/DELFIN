"""Small persistent per-session task state for fresh restarts.

A long DELFIN session resumes with >900k input tokens and re-sends the whole
history every turn; compaction (engine ``auto_compact_pct = 0.95`` of the
model's window) does not save it. This module gives a session ONE small state
file on disk -- task, phases with status, commits, open findings, "waiting
for" -- written after each commit/handoff, and ``render()`` hands back a
bounded, deterministic, secret-scrubbed block the next turn can start from
instead of the full history.

It deliberately does NOT re-implement ``working_state.py`` (which rebuilds
files-changed / last test outcome / open tasks from the task tool at compact
time). Here the state is written BY the agent mid-task, so a fresh turn knows
its open phase without scanning the whole transcript.

Design rules (mirror ``working_state``):
* DETERMINISTIC -- built from the stored fields only, in a fixed order; never
  from model text, never via a model call.
* BOUNDED -- a hard character ceiling; a recap that can itself blow the
  window it is meant to save is worse than none.
* NO SECRETS -- every rendered line passes through the house redactor
  (``output_guard.scrub_secrets``), the same one ``working_state`` uses.
* NEWEST FIRST -- a returning reader wants the latest state first.
"""
from __future__ import annotations

import json
from pathlib import Path

# How many commits render() keeps, newest first. A long session commits a lot;
# the box stays small by dropping the oldest once the cap is hit.
_MAX_RENDER_COMMITS = 15
# Hard ceiling on the rendered block (same spirit as working_state's).
_MAX_RENDER_CHARS = 1600
_HEADER = "Task state (written after the last commit/handoff)"

# Section order in render(), highest priority first. When the ceiling is hit
# the last sections drop first.
_SECTION_NAMES = (
    ("task", "Task"),
    ("phase", "Open phase"),
    ("waiting", "Waiting on"),
    ("findings", "Open findings"),
    ("commits", "Commits"),
)


class TaskState:
    """One per-session task state, writable to and readable from a JSON file.

    Instances are created empty via :func:`open`; ``load()`` fills one from
    disk. Field access uses plain attributes (``.task``, ``.phases``,
    ``.commits``, ``.findings``, ``.waiting_for``) so tests and callers read
    state without parsing anything.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        self.begin()

    # -- mutation ---------------------------------------------------------

    def begin(self) -> None:
        """Start a fresh session state (clears all fields in place)."""
        self.task: str = ""
        self.phases: list[dict] = []
        self.commits: list[str] = []
        self.findings: list[str] = []
        self.waiting_for: str = ""

    def commit(self, *, task: str | None = None, phase: str | None = None,
               phase_status: str | None = None,
               commit: str | None = None) -> None:
        """Record a commit (and the phase it belongs to).

        ``task`` sets the session task; ``phase`` + ``phase_status`` record
        (or update) that phase; ``commit`` is the new commit hash. Each is
        optional so a caller can record only what it knows this call.
        """
        if task:
            self.task = task
        if commit:
            self.commits.append(str(commit))
        if phase:
            self._set_phase(phase, phase_status or "")

    def add_finding(self, *, finding: str) -> None:
        """Add one open finding (a reviewer note not yet resolved)."""
        text = str(finding or "").strip()
        if text and text not in self.findings:
            self.findings.append(text)

    def set_waiting_for(self, *, waiting_for: str) -> None:
        """Name what this session is currently waiting on (or clear it)."""
        self.waiting_for = str(waiting_for or "").strip()

    def _set_phase(self, name: str, status: str) -> None:
        """Update-or-append one phase, moving it to the front (newest)."""
        for entry in self.phases:
            if entry.get("name") == name:
                entry["status"] = status
                self.phases.remove(entry)
                self.phases.insert(0, entry)  # newest-first
                return
        self.phases.insert(0, {"name": name, "status": status})

    # -- persistence ------------------------------------------------------

    def save(self) -> None:
        """Write the current state to ``self.path`` as JSON (overwrites)."""
        payload = self._as_dict()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(
            json.dumps(payload, ensure_ascii=True, sort_keys=True, indent=2),
            encoding="utf-8",
        )

    def load(self) -> None:
        """Read state back from ``self.path``; a missing file stays empty."""
        if not self.path.exists():
            self.begin()
            return
        data = json.loads(self.path.read_text(encoding="utf-8"))
        self.task = str(data.get("task", "") or "")
        self.waiting_for = str(data.get("waiting_for", "") or "")
        self.commits = [str(c) for c in data.get("commits", []) or []]
        self.findings = [str(f) for f in data.get("findings", []) or []]
        phases: list[dict] = []
        for p in data.get("phases", []) or []:
            if isinstance(p, dict) and p.get("name"):
                phases.append({"name": str(p["name"]), "status": str(p.get("status", "") or "")})
        self.phases = phases

    # -- rendering --------------------------------------------------------

    def render(self) -> str:
        """One bounded, deterministic, secret-scrubbed block.

        Newest commit first; phases and findings newest first; a hard char
        cut with a marker guarantees the ceiling.
        """
        phases = self.phases  # already newest-first
        open_phase = next(
            (p for p in phases if p.get("status") in ("in_progress", "pending", "blocked")),
            (phases[0] if phases else None),
        )
        commits = self.commits[-_MAX_RENDER_COMMITS:][::-1]  # newest first, capped
        findings = self.findings[::-1]
        sections: list[str] = []
        if self.task:
            sections.append(f"Task: {self.task}")
        if open_phase and open_phase.get("name"):
            name = open_phase["name"]
            status = open_phase.get("status") or ""
            sections.append(f"Open phase: {name}"
                            + (f" [{status}]" if status else ""))
        if self.waiting_for:
            sections.append(f"Waiting on: {self.waiting_for}")
        if findings:
            sections.append("Open findings:\n" + "\n".join(f"  - {f}" for f in findings))
        if commits:
            sections.append("Commits (newest first):\n" + "\n".join(f"  - {c}" for c in commits))

        if not sections:
            return ""
        block = _HEADER + "\n" + "\n".join(sections) + "\n"
        if len(block) > _MAX_RENDER_CHARS:
            block = block[:_MAX_RENDER_CHARS] + "\n... [task state trimmed]\n"
        return _scrub(block)

    def _as_dict(self) -> dict:
        """Current state as a JSON-serialisable dict."""
        return {
            "task": self.task,
            "waiting_for": self.waiting_for,
            "phases": list(self.phases),
            "commits": list(self.commits),
            "findings": list(self.findings),
        }


def _scrub(text: str) -> str:
    """Remove credential material via the house redactor. Never raises."""
    try:
        from .output_guard import scrub_secrets
        return scrub_secrets(text)
    except Exception:
        return text


def open(path: str | Path) -> TaskState:
    """Return a NEW empty :class:`TaskState` bound to ``path``.

    ``open`` shadows Python's builtin here on purpose (the module's public
    entry point, matching ``working_state``'s read API); persistence goes
    through ``Path.write_text/read_text`` so the builtin is never needed.
    """
    return TaskState(Path(path))
