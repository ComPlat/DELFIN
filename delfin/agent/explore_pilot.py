"""SA1 pilot: a measured, narrow read-only explore helper.

This module is a MEASUREMENT HARNESS, not a new capability. DELFIN already
has sub-agents (``delfin/agent/subagents.py``, preset ``explore`` at
:1260-1276: mode "plan", read-only tool set). The wave task forbids them by
default; this pilot measures whether ONE narrow use pays off, and the
helper becomes a recommendation (a skill) only if the measurement says so.

Phase 1 — the guard. The pilot setting allows ONLY:
- the ``explore`` sub-agent type (read-only filesystem + grep + web search),
- at most one sub-agent running per session,
- no nesting (a delegate may not start another delegate),
- no background runs (the harness measures one blocking run per arm).

Everything else is refused with the reason stated — a refusal names a
ROUTE, so every refusal here says what is not part of the pilot setting.

The real enforcement point for sub-agent calls is the engine's
``subagent_runner`` binding (engine.py:961-992, a protected file). This
guard lives at the wrapper the pilot harness calls; the wiring question
(engine hook vs. wrapper-only) is tracked in the package's PR notes.

Phase 2 — the measurement harness (``run_pilot``) runs the same small set
of read-heavy tasks twice: once plain, once allowed to call the explore
sub-agent, and compares tokens, wall time, tool calls and correctness.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

from . import subagents as _sa

# The one type the pilot setting allows.
PILOT_TYPE = "explore"

#: Background runs (reserve_running) are long-lived: a reserved slot counts
#: as running until the delegate is reaped, and is exactly where two
#: delegates could overlap.
BOTH_CHANNELS = "pilot setting allows at most one explore sub-agent running per session"

#: The real runner, bound at import so tests can monkeypatch it and so the
#: refusal paths below demonstrably never reach it.
_real_run_subagent = _sa.run_subagent


def check_delegate(subagent_type: str, *, slot: Any = None,
                   depth: int = 0, isolation: str = "") -> tuple[bool, str]:
    """Should this delegation be allowed under the pilot setting?

    Returns ``(ok, reason)``. ``reason`` is "" when allowed; otherwise a
    sentence stating what the pilot setting does not allow. Pure: no I/O,
    no state of its own — the one-running slot is passed in.
    """
    if subagent_type != PILOT_TYPE:
        return False, (
            f"pilot setting allows only the {PILOT_TYPE!r} sub-agent type; "
            f"{subagent_type!r} is refused")
    if depth > 0:
        return False, (
            "pilot setting allows no nesting: a delegate may not start "
            "another delegate")
    if isolation:
        return False, (
            f"pilot setting allows no isolation={isolation!r}: the pilot "
            "runs read-only in the harness's own workspace")
    if slot is not None and not _slot_free(slot):
        return False, (
            f"{BOTH_CHANNELS}; one is already running")
    return True, ""


# ---------------------------------------------------------------------------
# The one-running-per-session slot
# ---------------------------------------------------------------------------

def pilot_slot() -> Any:
    """A fresh slot object for one pilot run (context manager)."""
    return _PilotSlot()


def _slot_free(slot: Any) -> bool:
    return not getattr(slot, "claimed", False)


def claim_slot(slot: Any) -> None:
    """Mark the slot as taken (one sub-agent running)."""
    slot.claimed = True


def reap_slot(slot: Any) -> None:
    """Release the slot (the delegate finished or was reaped)."""
    slot.claimed = False


class _PilotSlot:
    """One delegation channel of the pilot setting.

    ``depth`` is a property of the slot, not of the call: the slot says
    who is allowed to delegate through it (the session, depth 0) and the
    guard is installed with the SAME slot for anything it lets run, so a
    delegate that re-enters the pilot entry is seen at depth 1 and
    refused — not handed a fresh unconstrained slot.
    """

    def __init__(self, depth: int = 0) -> None:
        self.claimed = False
        self.depth = depth

    def __enter__(self) -> "_PilotSlot":
        return self

    def __exit__(self, *exc) -> None:
        self.claimed = False


# ---------------------------------------------------------------------------
# The guarded entry the harness calls
# ---------------------------------------------------------------------------

def run_pilot_subagent(*, subagent_type: str, description: str, prompt: str,
                       parent_client: Any, parent_perms: Any,
                       model: str = "", inherit_context: bool = True,
                       depth: int = 0, _pilot_slot: Any = None,
                       **kwargs: Any) -> Any:
    """``run_subagent`` under the pilot setting.

    Every refusal returns a ``SubagentResult`` with the reason in
    ``error`` — the same shape the real runner returns for its own
    refusals — and the real runner is never called. The slot and its
    depth travel with the call: the delegate's own runner is built with
    ``depth=slot.depth + 1`` baked in, so nesting is refused on the real
    path, not only when a caller passes depth by hand.

    Like the real ``run_subagent``, this never raises: a delegate that
    crashes returns an error result, so the harness never loses an arm.
    """
    slot = _pilot_slot if _pilot_slot is not None else _PilotSlot()
    ok, why = check_delegate(subagent_type, slot=slot, depth=depth)
    if not ok:
        return SimpleNamespace(subagent_type=subagent_type,
                               description=description,
                               final_text="", error=why)
    claim_slot(slot)
    try:
        return _real_run_subagent(
            subagent_type=subagent_type, description=description,
            prompt=prompt, parent_client=parent_client,
            parent_perms=parent_perms, model=model,
            inherit_context=inherit_context, **kwargs)
    except Exception as exc:  # the real runner's contract: never raises
        return SimpleNamespace(subagent_type=subagent_type,
                               description=description,
                               final_text="", error=f"delegate failed: {exc}")
    finally:
        slot.claimed = False


def child_runner(parent_perms: Any, slot: Any, runner: Any = None) -> Any:
    """The runner installed for a delegate the pilot let run.

    Depth is a property of whose runner it is: anything a delegate tries
    to start goes through a slot one level deeper, where the guard
    refuses it. The returned callable binds (slot, depth) so the guard
    cannot be forgotten on the child path.
    """
    child_slot = _PilotSlot(depth=slot.depth + 1)
    return lambda **kw: run_pilot_subagent(
        _pilot_slot=child_slot, parent_client=None,
        parent_perms=parent_perms, **kw) if runner is None else runner(
        _pilot_slot=child_slot, **kw)


# ---------------------------------------------------------------------------
# Phase 2 — the measurement harness
# ---------------------------------------------------------------------------

import re as _re
from dataclasses import dataclass, field as _field


_CITATION = _re.compile(
    r"\b([A-Za-z0-9_./\-]+\.[a-z]{1,4}):(\d+)\b")


@dataclass
class PilotTask:
    task_id: str
    prompt: str


@dataclass
class CitationReport:
    claims: int = 0
    bad: list[str] = _field(default_factory=list)


@dataclass
class ArmResult:
    arm: str
    success: bool = False
    input_tokens: int = 0
    output_tokens: int = 0
    tool_calls: int = 0
    duration_s: float = 0.0
    final_text: str = ""
    citations_ok: bool = False
    citation_report: CitationReport = _field(default_factory=CitationReport)


@dataclass
class TaskRow:
    task_id: str
    arms: dict[str, ArmResult] = _field(default_factory=dict)


def verify_citations(text: str, *, root: Any = None,
                     expect: str = "") -> CitationReport:
    """Every ``path:line`` claim in ``text`` must exist under ``root`` and
    the line must contain ``expect`` when that is given (the line is read
    from the real file — no mock stands in for the filesystem).

    A claim whose file is missing, whose line number is out of range or
    whose content contradicts ``expect`` is listed in ``bad`` with the
    reason; the report counts how many claims were checked.
    """
    report = CitationReport()
    for m in _CITATION.finditer(str(text or "")):
        rel, line_no = m.group(1), int(m.group(2))
        report.claims += 1
        path = (Path(root) / rel) if root else Path(rel)
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            report.bad.append(f"{rel}:{line_no}: file missing or unreadable")
            continue
        if not (1 <= line_no <= len(lines)):
            report.bad.append(f"{rel}:{line_no}: out of range "
                              f"({len(lines)} lines)")
            continue
        if expect and expect not in lines[line_no - 1]:
            report.bad.append(f"{rel}:{line_no}: line does not mention "
                              f"{expect!r}")
    return report


def run_pilot(tasks: list[PilotTask], *, runner: Any,
              root: Any = ".") -> list[TaskRow]:
    """Run every task twice — arm ``plain`` and arm ``explore`` — through
    ``runner(task, arm=...)``, and check each arm's final text for
    correct ``path:line`` citations against ``root``.

    ``runner`` returns a dict with success, input_tokens, output_tokens,
    tool_calls, duration_s and final_text; in production it is
    ``benchmark_runner.run_task`` (protected, read-only for this package)
    called once per arm. Here the shape is the contract, not the runner.
    """
    rows: list[TaskRow] = []
    for task in tasks:
        row = TaskRow(task_id=task.task_id)
        for arm in ("plain", "explore"):
            raw = runner(task, arm=arm) or {}
            report = verify_citations(raw.get("final_text", ""), root=root)
            row.arms[arm] = ArmResult(
                arm=arm, success=bool(raw.get("success")),
                input_tokens=int(raw.get("input_tokens") or 0),
                output_tokens=int(raw.get("output_tokens") or 0),
                tool_calls=int(raw.get("tool_calls") or 0),
                duration_s=float(raw.get("duration_s") or 0.0),
                final_text=str(raw.get("final_text") or ""),
                citations_ok=(report.claims == 0 or not report.bad),
                citation_report=report)
        rows.append(row)
    return rows


@dataclass
class PilotVerdict:
    explore_wins: bool | None = None
    equal_correctness: bool | None = None
    detail: str = ""


def pilot_verdict(rows: list[TaskRow]) -> PilotVerdict:
    """Does the explore arm win? Wins = strictly fewer tool calls and
    less wall time at EQUAL correctness (every arm either fully cited or
    both arms equally wrong). One task's noise is not a win: the verdict
    is None while the evidence is mixed."""
    v = PilotVerdict()
    wins = ties = losses = 0
    correct_all = True
    for row in rows:
        plain, pilot = row.arms.get("plain"), row.arms.get("explore")
        if plain is None or pilot is None:
            continue
        if plain.citations_ok != pilot.citations_ok:
            correct_all = False
        if not (plain.success and pilot.success):
            continue
        if pilot.tool_calls < plain.tool_calls and pilot.duration_s <= plain.duration_s:
            wins += 1
        elif pilot.tool_calls > plain.tool_calls or pilot.duration_s > plain.duration_s:
            losses += 1
        else:
            ties += 1
    v.equal_correctness = correct_all
    if not rows:
        v.detail = "no tasks measured"
        return v
    if wins and not losses:
        v.explore_wins = True
    elif losses and not wins:
        v.explore_wins = False
    v.detail = (f"{len(rows)} tasks: {wins} explore-wins, {ties} ties, "
                f"{losses} losses; correctness equal: {correct_all}")
    return v


# ---------------------------------------------------------------------------
# Phase 3 — the production runner adapter (SP1)
# ---------------------------------------------------------------------------

# The keys run_pilot reads from every arm result (explore_pilot.py:254-257).
# The production runner is benchmark_runner.run_task (protected, called
# never edited), which returns a BenchmarkResult that carries NO
# final_text — only ``text_excerpt``, truncated to 400 chars on success
# (benchmark.py:2012-2013). A citation answer read through the excerpt is
# graded on 400 characters, so the adapter captures the FULL text the
# real ``run_once`` produced and forwards it as ``final_text``.

_ARMS = ("plain", "explore")


def make_pilot_runner(*, model: str, run_task: Any = None,
                      full_text: str = "",
                      engine_factory_builder: Any = None,
                      default_task: Any = None,
                      **run_task_kwargs: Any) -> Any:
    """Build the runner ``run_pilot`` expects around
    ``benchmark_runner.run_task``.

    The returned callable is ``runner(task, arm=...) -> dict`` with the
    six keys run_pilot reads (success, input_tokens, output_tokens,
    tool_calls, duration_s, final_text).

    - ``model``: the model both arms run under (one measurement = one
      model; the arms differ ONLY in the delegation permission).
    - ``run_task``: the production runner; default
      ``benchmark_runner.run_task`` (imported lazily so a test can hand a
      fake without importing the protected module first).
    - ``full_text``: the full final text of the LAST run's model output,
      captured by the run_once wrapper below. A caller that already has
      it (tests) may pass it directly.
    - ``engine_factory_builder(arm) -> factory``: optional hook so the
      explore arm can hand run_task an engine whose permissions bind
      ``run_pilot_subagent`` (the pilot guard) while the plain arm's
      engine is built untouched. ``None``: no engine_factory is passed
      and both arms use run_task's default factory — the wiring used by
      the CLI entry point, which builds its own engines.
    - ``default_task``: benchmark.Task defaults (task_class, mode,
      permission_mode) applied to tasks that do not carry them.
    - extra kwargs are forwarded to run_task verbatim (backend,
      provider, max_tokens, ...).
    """
    if run_task is None:
        from .benchmark_runner import run_task as _production_run_task
        run_task = _production_run_task
    task_defaults = dict(default_task) if isinstance(default_task, dict) \
        else (vars(default_task).copy() if default_task is not None else {})

    def _to_benchmark_task(task: Any) -> Any:
        """Map a PilotTask (or a Task-like dict) onto the shape run_task
        needs. A real benchmark.Task passes through unchanged; a
        PilotTask gets its defaults from ``default_task``."""
        from .benchmark import Task as _Task
        if isinstance(task, _Task):
            return task
        tid = getattr(task, "task_id", None)
        if tid is None:
            tid = getattr(task, "id", "")
        prompt = getattr(task, "prompt", "")
        fields = dict(task_defaults)
        fields.setdefault("task_class", "pilot")
        # solo: the pilot measures a read-heavy answering task; the
        # benchmark's modes are solo|dashboard|office|research
        # (benchmark.py:301-303) and solo is the one with no fixture
        # workspace dependency.
        fields.setdefault("mode", "solo")
        fields["id"] = tid
        fields["prompt"] = prompt
        # Only valid Task keys survive: an unknown key in the YAML is a
        # spelling error the caller must see, not a silent default.
        import dataclasses as _dc
        known = {f.name for f in _dc.fields(_Task)}
        dropped = sorted(set(fields) - known - {"id", "prompt"})
        if dropped:
            raise ValueError(
                f"task {tid}: unknown benchmark.Task field(s) {dropped}")
        fields.pop("id", None)
        fields.pop("prompt", None)
        return _Task(id=tid, prompt=prompt, **fields)

    def _arm_run_once_factory(captured: dict) -> Any:
        """A run_once that delegates to the real one and keeps the FULL
        text — the excerpt in BenchmarkResult is not enough for citation
        grading."""
        from .cli import _run_once as _real_run_once

        def run_once(engine, prompt, **kw):
            raw = _real_run_once(engine, prompt, **kw)
            captured["full_text"] = str(raw.get("text") or "")
            return raw
        return run_once

    def runner(task: Any, arm: str = "plain") -> dict:
        if arm not in _ARMS:
            raise ValueError(
                f"unknown arm {arm!r}; the pilot measures {_ARMS}")
        captured: dict = {}
        kwargs = dict(run_task_kwargs)
        kwargs["model"] = model
        kwargs["run_once"] = _arm_run_once_factory(captured)
        if engine_factory_builder is not None:
            if arm == "explore":
                # The wiring the package exists to measure: the engine
                # built for this arm binds the pilot guard, so the task
                # MAY delegate (explore only, one at a time, depth 0).
                kwargs["engine_factory"] = engine_factory_builder("explore")
            # The plain arm passes nothing: run_task's default factory
            # builds an engine with no subagent_runner binding, so a
            # delegate attempt is refused upstream (engine.py:972-974).
        result = run_task(_to_benchmark_task(task), **kwargs)
        final_text = captured.get("full_text") or ""
        if not final_text:
            # A fake caller without a real run_once falls back to the
            # excerpt (or an explicit full_text) rather than hiding the
            # answer entirely — but the fallback is visible.
            final_text = full_text or str(
                getattr(result, "text_excerpt", "") or "")
        return {
            "success": bool(getattr(result, "success", False)),
            "input_tokens": int(getattr(result, "input_tokens", 0) or 0),
            "output_tokens": int(getattr(result, "output_tokens", 0) or 0),
            "tool_calls": int(getattr(result, "tool_calls", 0) or 0),
            "duration_s": float(getattr(result, "duration_s", 0.0) or 0.0),
            "final_text": final_text,
        }

    return runner
