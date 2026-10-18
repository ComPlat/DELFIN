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
# Per-value clip for the new free-text fields, so a hostile or runaway
# value burns the ceiling in a bounded number of lines.
_MAX_FIELD_CHARS = 300
# Newest decisions/done entries kept.
_MAX_DECISIONS = 10
_MAX_DONE = 15
_HEADER = "Task state (written after the last commit/handoff)"
# Default store filename under a workspace. working_state reads the
# same store at compaction (the second writer of the same state).
DEFAULT_FILENAME = "task_state.json"

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
        self.success_criterion: str = ""
        self.decisions: list[dict] = []
        self.done: list[dict] = []
        self.person_instruction: str = ""
        self.person_instruction_is_fresh: bool = False

    def commit(self, *, task: str | None = None, phase: str | None = None,
               phase_status: str | None = None,
               commit: str | None = None) -> None:
        """Record a commit (and the phase it belongs to).

        ``task`` sets the session task; ``phase`` + ``phase_status`` record
        (or update) that phase; ``commit`` is the new commit hash. Each is
        optional so a caller can record only what it knows this call.

        Recording a commit is the agent ACTING on the state; the person's
        instruction (if any) is therefore no longer fresh afterwards --
        a stored plan never resumes against a newer instruction, and an
        executed instruction is not presented as live steering.
        """
        if task:
            self.task = task
        if commit:
            self.commits.append(str(commit))
            self.person_instruction_is_fresh = False
        if phase:
            self._set_phase(phase, phase_status or "")

    def set_person_instruction(self, text: str) -> None:
        """The person's newest instruction, verbatim.

        It REPLACES the stored plan as the goal (the person's newest word
        outranks every stored plan) and is marked fresh until the agent
        next acts on the state. Only the person's own channel sets this;
        a session/operator message must never reach this setter.
        """
        self.person_instruction = str(text or "").strip()
        self.person_instruction_is_fresh = bool(self.person_instruction)
        if self.person_instruction:
            self.task = self.person_instruction

    def set_success_criterion(self, text: str) -> None:
        """What 'done' means for the session task."""
        self.success_criterion = str(text or "").strip()

    def add_decision(self, *, decision: str, reason: str = "") -> None:
        """Record one decision with why (newest first, capped)."""
        entry = {"decision": str(decision or "").strip(),
                 "reason": str(reason or "").strip()}
        if not entry["decision"]:
            return
        if entry not in self.decisions:
            self.decisions.insert(0, entry)
            del self.decisions[_MAX_DECISIONS:]

    def add_done(self, *, what: str, evidence: str = "") -> None:
        """Record completed work WITH its evidence (newest first, capped)."""
        entry = {"what": str(what or "").strip(),
                 "evidence": str(evidence or "").strip()}
        if not entry["what"]:
            return
        if entry not in self.done:
            self.done.insert(0, entry)
            del self.done[_MAX_DONE:]

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
        """Write the current state to ``self.path`` as JSON (atomic).

        Write to a temp file in the same directory, then rename over the
        target, so a crash or a concurrent reader never sees a torn JSON
        file (the main way a store corrupts).
        """
        payload = self._as_dict()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_name(self.path.name + ".tmp")
        tmp.write_text(
            json.dumps(payload, ensure_ascii=True, sort_keys=True, indent=2),
            encoding="utf-8",
        )
        tmp.replace(self.path)

    def load(self) -> None:
        """Read state back from ``self.path``; a missing or corrupt file
        degrades to an empty state instead of raising (best-effort, like
        working_state — a broken store must never break the caller)."""
        if not self.path.exists():
            self.begin()
            return
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            self.begin()
            return
        self.task = str(data.get("task", "") or "")
        self.waiting_for = str(data.get("waiting_for", "") or "")
        self.commits = [str(c) for c in data.get("commits", []) or []]
        self.findings = [str(f) for f in data.get("findings", []) or []]
        self.success_criterion = str(
            data.get("success_criterion", "") or "")
        self.person_instruction = str(
            data.get("person_instruction", "") or "")
        # Freshness is a session-live property: a state loaded from disk
        # (restart, other model) is by definition not mid-act, so a
        # restored instruction is presented as the standing goal, not as
        # live steering. An old file without the field stays False.
        self.person_instruction_is_fresh = bool(
            data.get("person_instruction_is_fresh", False))
        # A state loaded from disk (restart, other model) is by
        # definition not mid-act: the instruction is presented as the
        # standing goal, never as live steering, however it was saved.
        self.person_instruction_is_fresh = False
        decisions: list[dict] = []
        for d in data.get("decisions", []) or []:
            if isinstance(d, dict) and d.get("decision"):
                decisions.append({"decision": str(d["decision"]),
                                  "reason": str(d.get("reason", "") or "")})
        self.decisions = decisions[:_MAX_DECISIONS]
        done: list[dict] = []
        for d in data.get("done", []) or []:
            if isinstance(d, dict) and d.get("what"):
                done.append({"what": str(d["what"]),
                             "evidence": str(d.get("evidence", "") or "")})
        self.done = done[:_MAX_DONE]
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
            sections.append(f"Task: {_flat(self.task)}")
        if open_phase and open_phase.get("name"):
            name = _flat(open_phase["name"])
            status = _flat(open_phase.get("status") or "")
            sections.append(f"Open phase: {name}"
                            + (f" [{status}]" if status else ""))
        if self.success_criterion:
            sections.append(
                f"Success criterion: {_flat(self.success_criterion)}")
        if self.waiting_for:
            sections.append(f"Waiting on: {_flat(self.waiting_for)}")
        if self.person_instruction:
            marker = " (fresh)" if self.person_instruction_is_fresh else ""
            sections.append(
                f"Person's latest instruction{marker}: "
                f"{_flat(self.person_instruction)}")
        if self.decisions:
            sections.append(
                "Decisions (newest first):\n"
                + "\n".join(
                    f"  - {_flat(d.get('decision', ''))}"
                    + (f" — {_flat(d.get('reason', ''))}" if d.get("reason") else "")
                    for d in self.decisions))
        if self.done:
            sections.append(
                "Done (newest first):\n"
                + "\n".join(
                    f"  - {_flat(d.get('what', ''))}"
                    + (f" [{_flat(d.get('evidence', ''))}]"
                       if d.get("evidence") else "")
                    for d in self.done))
        if findings:
            sections.append("Open findings:\n" + "\n".join(f"  - {_flat(f)}" for f in findings))
        if commits:
            sections.append("Commits (newest first):\n" + "\n".join(f"  - {_flat(c)}" for c in commits))

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
            "success_criterion": self.success_criterion,
            "decisions": list(self.decisions),
            "done": list(self.done),
            "person_instruction": self.person_instruction,
            "person_instruction_is_fresh": self.person_instruction_is_fresh,
        }


def _scrub(text: str) -> str:
    """Remove credential material via the house redactor. Never raises."""
    try:
        from .output_guard import scrub_secrets
        return scrub_secrets(text)
    except Exception:
        return text


def _flat(text: str) -> str:
    """Flatten a stored value for rendering: newlines become " / ", so a
    hostile value can never open a new line that could pass for a section
    header the code did not write (render's lines are the only ones the
    reader may trust)."""
    return " / ".join(str(text or "").splitlines())


def open(path: str | Path) -> TaskState:
    """Return a NEW empty :class:`TaskState` bound to ``path``.

    ``open`` shadows Python's builtin here on purpose (the module's public
    entry point, matching ``working_state``'s read API); persistence goes
    through ``Path.write_text/read_text`` so the builtin is never needed.
    """
    return TaskState(Path(path))
