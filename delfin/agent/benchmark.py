"""Canned-task benchmark suite for iterative agent optimisation.

Workflow:

1. Define tasks in ``delfin/agent/pack/benchmark/tasks.yaml`` (prompt +
   expected/forbidden signals + budget caps).
2. Run the suite against a model+mode combo — the runner (separate
   module) sends each prompt and collects the trajectory.
3. ``score_outcome(task, trajectory)`` produces a deterministic
   ``BenchmarkResult`` (success bool + 0-100 quality + components).
4. ``write_run`` appends each result to a JSONL file under
   ``~/.delfin/benchmark_runs/<ts>_<model>.jsonl``.
5. ``compare_runs(baseline_path, candidate_path)`` produces a per-task
   delta table — the "did this profile knob change help?" verdict.

The data plane is deliberately model-free: it knows nothing about the
engine or providers, only about (prompt, trajectory, outcome).  The
runner that drives the engine and supplies the trajectory lives in a
separate module so we can unit-test the scoring rubric without a live
provider.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Optional

from . import pricing as _pricing

try:
    import yaml as _yaml
except ImportError:                                            # pragma: no cover
    _yaml = None  # type: ignore[assignment]


_DEFAULT_TASKS_PATH = (
    Path(__file__).resolve().parent / "pack" / "benchmark" / "tasks.yaml"
)
_DEFAULT_RUNS_DIR = Path.home() / ".delfin" / "benchmark_runs"


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Signal:
    """One regex signal against the trajectory."""

    pattern: str
    against: str = "any"        # text | action | tool_name | any
    optional: bool = False


@dataclass(frozen=True)
class Task:
    """One canned benchmark task."""

    id: str
    task_class: str
    mode: str                   # solo | dashboard | plan | quick | …
    prompt: str
    expected_signals: tuple[Signal, ...] = ()
    forbidden_signals: tuple[Signal, ...] = ()
    max_duration_s: float = 60.0
    max_cost_usd: float = 0.10
    max_tool_calls: int = 5
    # Optional behaviour tag: "plan" | "scout" | "verify" | "ask".
    # When set, ``score_outcome`` computes one trace-derived behaviour
    # flag for this task (did it plan/scout/verify/ask?) into
    # ``BenchmarkResult.behavior`` — orthogonal to the quality rubric.
    behavior: str = ""


@dataclass
class Trajectory:
    """What the runner observed during a single task execution.

    The scorer reads only this shape — the runner can build it from the
    engine event stream however it likes."""

    text: str = ""                  # concatenated assistant text
    actions: list[str] = field(default_factory=list)
    tool_calls: list[dict] = field(default_factory=list)
    duration_s: float = 0.0
    cost_usd: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0
    error: str = ""
    # Gate denials observed during the run; ``None`` when the runner had
    # no way to look. Unobserved is not zero.
    denials: Optional[int] = None
    # Where the checkout under test lives. An absolute prefix in a tool
    # input is a routing detail, exactly like the transport namespace on a
    # tool NAME, and it is dropped for the same reason -- see
    # ``_strip_checkout_prefix``. Empty means "not known", and then nothing
    # is rewritten.
    checkout_root: str = ""

    def as_string(self) -> str:
        parts = [self.text]
        for a in self.actions:
            parts.append(f"\nACTION: {a}")
        root = str(self.checkout_root or "")
        for c in self.tool_calls:
            name = _tool_semantic_name(c.get("name", ""))
            rendered = _strip_checkout_prefix(c.get("input", ""), root)
            parts.append(f"\nTOOL: {name}({rendered})")
        return "".join(parts)


@dataclass
class BenchmarkResult:
    """Score-card for one (task, model, profile) execution.

    With ``n_samples > 1`` this is an aggregated record over N
    repeated executions of the same task; ``quality_stdev`` and
    ``success_rate`` then expose run-to-run noise.  ``components``,
    ``duration_s`` etc are medians across the N replicates and
    ``cost_usd`` is the SUM of the N runs (the real spend), not a
    median — useful for budget tracking.
    """

    task_id: str
    task_class: str
    model: str
    profile_name: str = ""
    mode: str = ""
    ts: float = 0.0
    success: bool = False
    quality_0_100: int = 0
    components: dict[str, int] = field(default_factory=dict)
    duration_s: float = 0.0
    cost_usd: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0
    tool_calls: int = 0
    matched_signals: list[str] = field(default_factory=list)
    violated_signals: list[str] = field(default_factory=list)
    missing_signals: list[str] = field(default_factory=list)
    #: Matched in some replicates and missed in others. Only ever set by
    #: ``aggregate_replicates``; a single run has nothing to be flaky
    #: about. Kept apart from ``missing_signals`` because "never matched"
    #: and "matched twice of three" call for different next steps.
    flaky_signals: list[str] = field(default_factory=list)
    budget_violations: list[str] = field(default_factory=list)
    error: str = ""
    # A run whose request never reached the model. Not a score: the
    # endpoint had no capacity, or the connection failed, so nothing
    # about the model was observed. Booking that as a failing model
    # lowers a baseline permanently and silently, and the next real
    # regression then passes under it.
    unmeasured: bool = False
    # --- replicate-aware fields (N=1 → trivially the only sample) ---
    n_samples: int = 1
    quality_stdev: float = 0.0
    success_rate: float = 0.0           # fraction of N samples with success=True
    per_run_quality: list[int] = field(default_factory=list)
    per_run_success: list[bool] = field(default_factory=list)
    # --- forensic fields for pattern-bug-vs-real-fail diagnosis ---
    text_excerpt: str = ""              # first ≤400 chars of model output
    tool_names: list[str] = field(default_factory=list)  # tools the model actually called
    # --- trace-derived behaviour flags (behavioural-parity eval) ---
    # Only populated for tasks carrying a ``behavior:`` tag.  Per-run this
    # is {flag: 0|1}; aggregated across replicates it is {flag: rate}.
    # Kept OUT of ``components`` on purpose — ``quality = sum(components)``,
    # so mixing behaviour flags in would corrupt the quality score.
    behavior: dict[str, float] = field(default_factory=dict)
    # --- cost side: what the guards charge for what they catch ---
    # See caveat_count. Aggregated across replicates these are medians,
    # like duration and quality.
    caveats: int = 0
    answer_chars: int = 0
    # Gate denials during the run. ``None`` means NOT OBSERVED -- the
    # runner could not reach the audit log -- and must never be reported
    # as zero, which is the same sentence as "nothing was refused".
    denials: Optional[int] = None


# ---------------------------------------------------------------------------
# Task loading
# ---------------------------------------------------------------------------


def _coerce_signal(raw: Any) -> Signal:
    if isinstance(raw, Signal):
        return raw
    if isinstance(raw, str):
        return Signal(pattern=raw)
    if isinstance(raw, dict):
        return Signal(
            pattern=str(raw.get("pattern", "")),
            against=str(raw.get("against", "any")),
            optional=bool(raw.get("optional", False)),
        )
    raise TypeError(f"Cannot coerce signal: {raw!r}")


def _coerce_task(raw: dict) -> Task:
    expected = tuple(_coerce_signal(s) for s in raw.get("expected_signals") or [])
    forbidden = tuple(_coerce_signal(s) for s in raw.get("forbidden_signals") or [])
    return Task(
        id=str(raw["id"]),
        task_class=str(raw.get("task_class", "")),
        mode=str(raw.get("mode", "solo")),
        prompt=str(raw.get("prompt", "")),
        expected_signals=expected,
        forbidden_signals=forbidden,
        max_duration_s=float(raw.get("max_duration_s", 60.0)),
        max_cost_usd=float(raw.get("max_cost_usd", 0.10)),
        max_tool_calls=int(raw.get("max_tool_calls", 5)),
        behavior=str(raw.get("behavior", "") or "").strip().lower(),
    )


def load_tasks(path: Path | None = None) -> list[Task]:
    """Load tasks from a YAML file. Defaults to the packaged suite PLUS
    every ``tasks_*.yaml`` sibling.

    The default-path case scans the parent directory of the main
    ``tasks.yaml`` and concatenates what it finds, so a suite is picked up
    by existing — neither the CLI nor this function needs to know about
    each one. The glob covers hand-written suites as well as the
    auto-generated ``tasks_auto_*`` ones: a suite that has to be
    registered somewhere is a suite that gets written and then silently
    never runs.

    Explicit ``path`` only loads that single file (caller controls).
    """
    if path is not None:
        return _load_one(path)
    out = _load_one(_DEFAULT_TASKS_PATH)
    parent = _DEFAULT_TASKS_PATH.parent
    for sibling in sorted(parent.glob("tasks_*.yaml")):
        if sibling != _DEFAULT_TASKS_PATH:
            out.extend(_load_one(sibling))
    return out


def _load_one(p: Path) -> list[Task]:
    if not p.exists():
        return []
    if _yaml is None:                                           # pragma: no cover
        raise RuntimeError("PyYAML is required to load benchmark tasks")
    raw = _yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    items = raw.get("tasks") if isinstance(raw, dict) else raw
    if not isinstance(items, list):
        return []
    return [_coerce_task(t) for t in items if isinstance(t, dict)]


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


_SIGNAL_AGAINST_VALUES = {"text", "action", "tool_name", "any"}


# Negation markers that waive a forbidden-pattern hit when they appear
# close to the match: the answer is REJECTING the term, not recommending
# it. Window is deliberately tight so a negation elsewhere in a long
# paragraph cannot launder a genuine recommendation.
# Word-boundary matching: "cannot" must NOT count as a negation of a
# nearby term (refusal detectors rely on it), while "NOT", "nicht",
# "never" etc. as standalone words do.
# Tool names reach the trajectory namespaced by the transport
# (``mcp__kit-coding__write_file``). The namespace is a routing detail,
# not semantics, and a pattern written against the tool's real name must
# not miss because of it.
_TOOL_NAMESPACE_RE = re.compile(r"^mcp__[^_]+(?:_[^_]+)*__")


def _tool_semantic_name(name: str) -> str:
    """Tool name without its transport namespace prefix."""
    return _TOOL_NAMESPACE_RE.sub("", str(name or ""))


def _strip_checkout_prefix(rendered: Any, root: str) -> str:
    """Tool input with the checkout's absolute prefix removed.

    Same argument as ``_tool_semantic_name`` one field over: where the
    checkout happens to live is a routing detail, and a pattern written
    against a repo-relative path must not miss because of it.

    Measured 2026-08-14. ``gen_report_unverified`` expects
    ``TOOL: write_file([^\\n]{0,120}?tests/fixtures/user_project_workspace/
    export.py)``. The gap exists so a match cannot spill across a long
    input into an unrelated path -- but it is counted from the START of
    the rendered input, which begins with the absolute file path. From the
    ordinary checkout that prefix is 61 characters and the signal fires;
    running the identical code from a deeper directory made it 135 and the
    task FAILED. The suite reported a product defect that did not exist,
    and the only thing that had changed was the working directory.

    A measurement whose answer depends on where it is run is not a
    measurement. The prefix is dropped so the window measures what it was
    written to measure: the distance from the tool call to the path.
    """
    s = str(rendered or "")
    root = str(root or "").rstrip("/")
    if not root:
        return s
    return s.replace(root + "/", "").replace(root, "")


# Markdown emphasis is presentation, not content: an answer that writes
# "werde ich **nicht** ausführen" must match a pattern for "nicht
# ausführen". Emphasis markers are removed before matching so a model's
# formatting choice cannot decide whether a signal fires.
# Only ``*`` and backticks: underscores carry meaning here — stripping
# them would break identifiers like user_project_workspace or write_file.
_EMPHASIS_RE = re.compile(r"(\*{1,3}|`{1,3})(?=\S)|(?<=\S)(\*{1,3}|`{1,3})")


def _strip_emphasis(text: str) -> str:
    """Drop markdown emphasis/code markers, preserving every other char."""
    try:
        return _EMPHASIS_RE.sub("", text or "")
    except Exception:
        return text or ""


_NEGATION_RE = re.compile(
    r"(?i)\b(?:not|nicht|never|niemals|kein|keine|falsch|wrong|"
    r"invalid|statt|avoid|vermeide|no such)\b|instead of")
_NEGATION_WINDOW = 160


def _is_word_char(ch: str) -> bool:
    return ch.isalnum() or ch == "_"


def _match_is_negated(haystack: str, start: int, end: int) -> bool:
    """True when a forbidden-pattern match sits in an explicit negation
    context within the surrounding window. The matched span itself is
    excluded so a pattern can never self-negate."""
    lo = max(0, start - _NEGATION_WINDOW)
    hi = min(len(haystack), end + _NEGATION_WINDOW)
    # A window edge inside a word invents a word boundary: 'notes' cut to
    # 'not' matches the negation marker and silently waives a real
    # violation. Drop the partial words at both edges.
    while lo < start and lo > 0 and _is_word_char(
            haystack[lo - 1]) and _is_word_char(haystack[lo]):
        lo += 1
    while hi > end and hi < len(haystack) and _is_word_char(
            haystack[hi]) and _is_word_char(haystack[hi - 1]):
        hi -= 1
    ctx = haystack[lo:start] + " " + haystack[end:hi]
    return bool(_NEGATION_RE.search(ctx))


def _signal_matches(
    signal: Signal, traj: Trajectory, *, waive_negated: bool = False,
) -> bool:
    pat = signal.pattern
    if not pat:
        return False
    try:
        rx = re.compile(pat)
    except re.error:
        return False
    against = signal.against if signal.against in _SIGNAL_AGAINST_VALUES else "any"
    haystacks: list[str]
    if against == "text":
        haystacks = [traj.text]
    elif against == "action":
        haystacks = list(traj.actions)
    elif against == "tool_name":
        haystacks = [_tool_semantic_name(c.get("name", ""))
                     for c in traj.tool_calls]
    else:                                                       # any
        haystacks = [traj.as_string()]
    # Emphasis is formatting, not content — match the words, not the
    # markdown around them. The tool_name channel carries bare names.
    if against != "tool_name":
        haystacks = [_strip_emphasis(h) for h in haystacks]
    if not waive_negated:
        return any(rx.search(h or "") for h in haystacks)
    for h in haystacks:
        for m in rx.finditer(h or ""):
            if not _match_is_negated(h or "", m.start(), m.end()):
                return True
    return False


# ---------------------------------------------------------------------------
# Behaviour classifier  (behavioural-parity eval)
#
# Three of the four target behaviours are ORDER-derived ("plan BEFORE act",
# "read BEFORE edit", "run AFTER edit"), which the presence-only signal
# rubric above cannot express.  This block derives four booleans from the
# already-ordered ``Trajectory.tool_calls`` using the agent's REAL tool
# vocabulary (snake_case: ``read_file``/``edit_file``/``bash`` — the tool
# names this harness actually emits).  Pure functions, unit-tested against
# synthetic trajectories before any live/token spend.
# ---------------------------------------------------------------------------

# Read/scout tools (knowing state before acting).
_READ_TOOLS = frozenset({
    "read_file", "grep_file", "list_files", "search_docs", "read_section",
    "list_docs", "list_sections", "search_calcs", "get_calc_info",
    "calc_summary", "find_definition", "find_references", "notebook_read",
    "project_introspect", "view_image", "glob_files",
})
# File-mutating tools.
_MUTATE_TOOLS = frozenset({
    "write_file", "edit_file", "multi_edit", "apply_patch", "notebook_edit",
})
# Execution tools (running code = verifying).  ``bash`` is dual-use and
# split below into read-only vs acting by inspecting its command.
_EXEC_TOOLS = frozenset({"bash", "bash_background", "run_tests"})
_PLAN_TOOLS = frozenset({"task_create", "exit_plan_mode"})
_ASK_TOOLS = frozenset({"ask_user_question"})

# First token of a read-only shell command (scouting via bash, not mutating).
_RO_BASH_CMDS = frozenset({
    "ls", "cat", "head", "tail", "grep", "rg", "find", "pwd", "wc",
    "which", "file", "stat", "tree", "diff", "less", "more", "column",
    "sort", "uniq", "cut", "nl", "basename", "dirname", "realpath",
})
_RO_GIT_SUB = frozenset({
    "status", "log", "diff", "show", "branch", "ls-files", "rev-parse",
})
# Anything that mutates the filesystem disqualifies a bash cmd from "read-only".
_BASH_WRITE_RE = re.compile(
    r"(^|\s)(rm|mv|cp|mkdir|touch|tee|dd|chmod|chown|ln)\s"
    r"|>>?|;\s*>|sed\s+-i|\bpip\s+install|\bgit\s+(add|commit|push|checkout|reset)"
)


def _normalize_tool_name(name: str) -> str:
    """Strip an MCP server namespace so classification matches the bare tool.

    Live backends (KIT-Toolbox, MCP servers) advertise tools as
    ``mcp__<server>__<tool>`` (e.g. ``mcp__delfin-docs__read_file``).  The
    classifier reasons about the bare tool (``read_file``), so normalise
    first — otherwise a genuine read/edit/run is invisible to the flags.
    """
    if name.startswith("mcp__"):
        parts = name.split("__")
        if len(parts) >= 3 and parts[-1]:
            return parts[-1]
    return name


def _input_str(inp: Any) -> str:
    """Best-effort string view of a tool_call input (dict or str)."""
    if inp is None:
        return ""
    if isinstance(inp, str):
        return inp
    if isinstance(inp, dict):
        try:
            return json.dumps(inp, ensure_ascii=False)
        except (TypeError, ValueError):
            return str(inp)
    return str(inp)


def _is_readonly_bash(cmd: str) -> bool:
    """True if a bash command only reads/inspects, never mutates."""
    c = (cmd or "").strip()
    if not c:
        return False
    if _BASH_WRITE_RE.search(c):
        return False
    # First meaningful token (skip a leading VAR=val env assignment).
    tokens = c.split()
    head = tokens[0]
    if "=" in head and not head.startswith("/"):
        tokens = tokens[1:]
        head = tokens[0] if tokens else ""
    head = head.rsplit("/", 1)[-1]                # strip any path prefix
    if head == "git":
        sub = tokens[1] if len(tokens) > 1 else ""
        return sub in _RO_GIT_SUB
    return head in _RO_BASH_CMDS


def _basename(path: str) -> str:
    return path.replace("\\", "/").rsplit("/", 1)[-1]


def _edited_paths(calls: list[dict]) -> list[str]:
    """Extract edited-file basenames from MUTATE calls' inputs."""
    out: list[str] = []
    for c in calls:
        if c["name"] not in _MUTATE_TOOLS:
            continue
        m = re.search(r'"(?:file_path|path|notebook_path)"\s*:\s*"([^"]+)"',
                      c["input"])
        if m:
            base = _basename(m.group(1))
            if base and base not in out:
                out.append(base)
    return out


def _classify_calls(tool_calls: list[dict]) -> list[dict]:
    """Annotate each ordered tool call with behaviour ``kinds``.

    Kinds: ``mutate`` | ``plan`` | ``ask`` | ``read`` (scout) |
    ``exec_act`` (acting execution: run_tests / non-read-only bash).
    A read-only bash is ``read`` (scouting), never ``exec_act``.
    """
    out: list[dict] = []
    for i, c in enumerate(tool_calls or []):
        name = _normalize_tool_name(str((c or {}).get("name") or ""))
        inp = _input_str((c or {}).get("input"))
        kinds: set[str] = set()
        if name in _MUTATE_TOOLS:
            kinds.add("mutate")
        if name in _PLAN_TOOLS:
            kinds.add("plan")
        if name in _ASK_TOOLS:
            kinds.add("ask")
        if name in _READ_TOOLS:
            kinds.add("read")
        if name == "subagent" and "explore" in inp.lower():
            kinds.add("read")                     # explore subagent = scouting
        if name in _EXEC_TOOLS:
            if name in ("bash", "bash_background") and _is_readonly_bash(inp):
                kinds.add("read")                 # read-only shell = scouting
            else:
                kinds.add("exec_act")
        out.append({"idx": i, "name": name, "input": inp, "kinds": kinds})
    return out


def _first_idx(calls: list[dict], kind: str) -> int | None:
    for c in calls:
        if kind in c["kinds"]:
            return c["idx"]
    return None


def _last_idx(calls: list[dict], kind: str) -> int | None:
    found: int | None = None
    for c in calls:
        if kind in c["kinds"]:
            found = c["idx"]
    return found


def _behavior_planned(calls: list[dict]) -> bool:
    """A PLAN tool call precedes the first mutate AND the first acting exec."""
    plan_idx = _first_idx(calls, "plan")
    if plan_idx is None:
        return False
    mut_idx = _first_idx(calls, "mutate")
    act_idx = _first_idx(calls, "exec_act")
    if mut_idx is not None and plan_idx > mut_idx:
        return False
    if act_idx is not None and plan_idx > act_idx:
        return False
    return True


def _behavior_scouted(calls: list[dict]) -> bool:
    """A read/scout precedes the first mutate (read before edit)."""
    read_idx = _first_idx(calls, "read")
    mut_idx = _first_idx(calls, "mutate")
    if read_idx is None:
        return False
    return mut_idx is None or read_idx < mut_idx


def _behavior_verified(calls: list[dict]) -> bool:
    """After the LAST mutate, the agent runs something OR reads an edited
    file back.  No mutate at all → nothing to verify → False."""
    last_mut = _last_idx(calls, "mutate")
    if last_mut is None:
        return False
    edited = _edited_paths(calls)
    for c in calls:
        if c["idx"] <= last_mut:
            continue
        if "exec_act" in c["kinds"]:
            return True                           # ran it (pytest/script/…)
        # read-back: read_file or read-only bash referencing an edited file
        if "read" in c["kinds"] and edited:
            if any(b in c["input"] for b in edited):
                return True
    return False


# Measured 2026-08-14 against kit.qwen3.5-397b-A17b. Given "Behebe den
# fehlschlagenden Test." the agent investigated, found nothing to fix, and
# wrote "**Frage:** Wo befindet sich ..." — an answer that LABELS its own
# question — and this pattern scored it as not having asked. It knew
# `welche[rs]?`, `spezifizier` and `meinst du`, and none of `wo`, `wie`,
# `was für`, `soll ich` or the label `Frage:`. `was für` is listed in the
# ask tasks' own expected_signals, so it was judged necessary there and
# forgotten here.
#
# German only on the recognising side, as everywhere in this framework:
# these match what a USER or the model writes, so they stay German while
# the code around them stays English. Word boundaries on every new term —
# a bare `wo` would fire inside "wollen", `wie` inside "wieder".
_ASK_TEXT_RE = re.compile(
    r"(?is)("
    r"which|what|welche[rsnm]?|which file|which molecule|specify|"
    r"do you mean|which one|should i|please (?:confirm|clarify|specify)|"
    r"could you (?:clarify|specify)|need(?: to know| more)|"
    r"can you (?:clarify|specify)"
    # --- German, added after the live measurement above ---
    r"|\bwas für\b|\bworum\b|\bwomit\b|\bworauf\b"
    r"|\bwo\b|\bwie\b|\bwoher\b|\bwohin\b"
    r"|\bsoll ich\b|\bsollen wir\b|\bmöchtest du\b|\bwillst du\b"
    r"|\bfrage\s*:|\brückfrage\b"
    r"|\b(?:benötige|brauche)\s+ich\b|\bich\s+(?:benötige|brauche)\b"
    r"|\bbitte\s+(?:gib|nenne|sag|teile)\b"
    r"|\bunklar\b|\bpräzisier|\bkonkretisier"
    r").*\?"
)


def _behavior_asked(calls: list[dict], traj: Trajectory) -> bool:
    """Asked-correctly: emitted a clarifying question (tool or prose) AND
    did NOT guess (no mutate, no acting exec).  Guessing overrides asking."""
    guessed = (_first_idx(calls, "mutate") is not None
               or _first_idx(calls, "exec_act") is not None)
    if guessed:
        return False
    asked_tool = _first_idx(calls, "ask") is not None
    asked_prose = bool(_ASK_TEXT_RE.search(traj.text or ""))
    return asked_tool or asked_prose


def behavior_flags(task: Task, traj: Trajectory) -> dict[str, int]:
    """Compute the ONE trace-derived behaviour flag for ``task``.

    Returns ``{}`` for untagged tasks, else a single-key dict such as
    ``{"planned": 1}`` / ``{"scouted": 0}``.  Only the flag matching the
    task's ``behavior`` tag is computed (a scout task says nothing about
    planning), so downstream rates are per-behaviour and unpolluted.
    """
    beh = (getattr(task, "behavior", "") or "").strip().lower()
    if not beh:
        return {}
    calls = _classify_calls(traj.tool_calls)
    if beh == "plan":
        return {"planned": int(_behavior_planned(calls))}
    if beh == "scout":
        return {"scouted": int(_behavior_scouted(calls))}
    if beh == "verify":
        return {"verified": int(_behavior_verified(calls))}
    if beh == "ask":
        return {"asked": int(_behavior_asked(calls, traj))}
    return {}


def behavior_rates(
    records: list[dict] | list[BenchmarkResult],
) -> dict[str, dict[str, float]]:
    """Aggregate behaviour flags across a run into per-flag rates.

    Returns ``{flag: {"rate": mean, "n": count}}`` over every record that
    carries that flag (a task's untagged behaviours contribute nothing).
    ``rate`` is the fraction of tasks that exhibited the behaviour.
    """
    sums: dict[str, float] = {}
    counts: dict[str, int] = {}
    for r in records:
        beh = r.behavior if isinstance(r, BenchmarkResult) else (r.get("behavior") or {})
        if not isinstance(beh, dict):
            continue
        for k, v in beh.items():
            if v is None:
                continue
            try:
                fv = float(v)
            except (TypeError, ValueError):
                continue
            sums[k] = sums.get(k, 0.0) + fv
            counts[k] = counts.get(k, 0) + 1
    return {
        k: {"rate": sums[k] / counts[k], "n": counts[k]}
        for k in sorted(sums) if counts[k]
    }


# Errors that mean the request never reached the model. Each was seen in
# a real run: the KIT gateway reporting no deployment for the requested
# model, a proxy connection error, an auth rejection, a gateway 5xx. The
# distinguishing property is that no token was ever generated, so there
# is nothing to score -- as opposed to a model that answered badly, or a
# tool that failed inside an otherwise live turn.
_TRANSPORT_FAILURE_MARKERS: tuple[str, ...] = (
    "no deployments available",
    "server connection error",
    "connection error",
    "connection refused",
    "connection reset",
    "service unavailable",
    "bad gateway",
    "gateway timeout",
    "temporarily unavailable",
    "rate limit",
    "429",
    "502", "503", "504",
    "authentication", "unauthorized", "invalid api key",
    "engine init failed",
)


def is_unmeasured(traj: "Trajectory") -> bool:
    """True when the run says nothing about the model.

    Observed 2026-08-12: two office tasks came back at quality 35 with a
    zero pass rate, and their entire recorded output was the harness's
    own three retry banners. The endpoint had no capacity. Scored as
    model failures, they took a suite from 9/11 to 7/11 and wrote that
    into the file baselines are compared against.

    The text check is deliberately paired with a second condition -- no
    tool call and no output beyond the retry notices -- so a turn that
    ran, did work and then hit a transport error near the end is still
    scored on what it did.
    """
    try:
        err = str(getattr(traj, "error", "") or "").lower()
        if not err:
            return False
        if not any(m in err for m in _TRANSPORT_FAILURE_MARKERS):
            return False
        if getattr(traj, "tool_calls", None):
            return False
        text = str(getattr(traj, "text", "") or "")
        stripped = _RETRY_NOTICE_RE.sub("", text).strip()
        return len(stripped) < 40
    except Exception:
        return False


_RETRY_NOTICE_RE = re.compile(
    r"(?m)^\s*[\u23f3\u26a0\u23f9\U0001f6d1][^\n]*$")


# ---------------------------------------------------------------------------
# The cost side of a guard
# ---------------------------------------------------------------------------
#
# Pass rate and quality are HIT-side numbers. Almost every fix in this
# framework has the same shape -- a path that used to be silent now says
# something -- and every one of those costs text, a check or a refusal.
# None of that can move a pass/fail number, so a suite can go on reading
# 11/11 while the answers grow a hedge apron nobody finishes reading. An
# answer carrying three caveats is as silent as one carrying none.
#
# These are the markers the framework itself emits, so this counts what
# was actually appended rather than guessing at wording:
_CAVEAT_MARKERS = (
    "> ⚠️ ",        # figure ledger, count-vs-enumeration,
                              # truncated output, ambiguous column
    "[verify] Caveat:",       # grounding, functional claims
)


def caveat_count(text: str) -> int:
    """How many caveats the framework appended to one answer.

    Counted, not judged: a caveat is worth its space exactly once. The
    number belongs beside the pass rate because it is the one quantity a
    wave of "say it out loud" fixes can push in the wrong direction
    without a single test going red.
    """
    body = str(text or "")
    return sum(body.count(marker) for marker in _CAVEAT_MARKERS)


def score_outcome(
    task: Task,
    traj: Trajectory,
    *,
    model: str = "",
    profile_name: str = "",
    ts: float | None = None,
) -> BenchmarkResult:
    """Deterministic rubric: success + quality 0-100 + component breakdown."""

    matched: list[str] = []
    missing: list[str] = []
    violated: list[str] = []
    budget_violations: list[str] = []

    # 1. Expected signals — each non-optional MUST match for success.
    success_required_ok = True
    optional_hits = 0
    optional_total = 0
    for idx, sig in enumerate(task.expected_signals):
        label = f"{task.id}.expected[{idx}]"
        hit = _signal_matches(sig, traj)
        if hit:
            matched.append(label)
            if sig.optional:
                optional_hits += 1
        else:
            if sig.optional:
                missing.append(label + ":optional")
            else:
                missing.append(label)
                success_required_ok = False
        if sig.optional:
            optional_total += 1

    # 2. Forbidden signals — any match flips success to False. A match
    # inside an explicit NEGATION context is waived: an answer that names
    # a fake keyword in order to warn against it ("the keywords are NOT
    # Nactel/Nactorb") shows exactly the grounded behavior the suite
    # rewards, and must not score as if it recommended the fake.
    for idx, sig in enumerate(task.forbidden_signals):
        label = f"{task.id}.forbidden[{idx}]"
        if _signal_matches(sig, traj, waive_negated=True):
            violated.append(label)

    unmeasured = is_unmeasured(traj)
    # An unmeasured run is not a pass and not a fail. `success` stays
    # False so nothing counts it as working; the flag is what keeps the
    # aggregates from counting it as broken.
    success = bool(success_required_ok and not violated and not traj.error)

    # 3. Budget checks (don't flip success on their own — visible in score).
    if traj.duration_s > task.max_duration_s > 0:
        budget_violations.append(
            f"duration_s {traj.duration_s:.1f}>{task.max_duration_s:.1f}"
        )
    if traj.cost_usd > task.max_cost_usd > 0:
        budget_violations.append(
            f"cost_usd {traj.cost_usd:.4f}>{task.max_cost_usd:.4f}"
        )
    if len(traj.tool_calls) > task.max_tool_calls > 0:
        budget_violations.append(
            f"tool_calls {len(traj.tool_calls)}>{task.max_tool_calls}"
        )

    # 4. Component scoring (0-100):
    #    success_pts  (40)   binary did-it-work
    #    routing_pts  (20)   expected signals + optional bonus
    #    speed_pts    (15)   within duration budget, scaled
    #    cost_pts     (10)   within cost budget, scaled
    #    clean_pts    (15)   no forbidden, no error, within tool-call budget
    success_pts = 40 if success else 0
    n_required = sum(1 for s in task.expected_signals if not s.optional)
    n_required_hit = n_required - sum(
        1 for m in missing if not m.endswith(":optional")
    )
    if n_required > 0:
        routing_base = int(round(15 * n_required_hit / n_required))
    else:
        routing_base = 15
    routing_bonus = (
        int(round(5 * optional_hits / max(1, optional_total)))
        if optional_total > 0
        else 0
    )
    routing_pts = routing_base + routing_bonus

    if task.max_duration_s > 0 and traj.duration_s >= 0:
        speed_ratio = min(1.0, traj.duration_s / task.max_duration_s)
        speed_pts = int(round(15 * (1.0 - speed_ratio)))
    else:
        speed_pts = 15

    # Cost points require a cost that was actually measured. The guard
    # used to be ``traj.cost_usd >= 0``, which no cost can fail, so a run
    # that reported 0.0 because nobody knew its model's price collected
    # the full 10 points for thrift it never demonstrated. A run IS
    # measured when the runner observed real spend (cost > 0), or when the
    # provider genuinely bills no USD — a local or quota-funded model
    # earns its zero. Anything else is unmeasured and scores nothing, with
    # the reason recorded so the missing points are explainable.
    if task.max_cost_usd > 0:
        price = _pricing.resolve(model)
        if traj.cost_usd > 0 or price.state == _pricing.NON_BILLING:
            cost_ratio = min(1.0, max(0.0, traj.cost_usd) / task.max_cost_usd)
            cost_pts = int(round(10 * (1.0 - cost_ratio)))
        else:
            cost_pts = 0
            budget_violations.append(f"cost_usd unmeasured ({price.reason})")
    else:
        cost_pts = 10

    clean_pts = 15
    if violated:
        clean_pts -= 8
    if traj.error:
        clean_pts -= 5
    if len(traj.tool_calls) > task.max_tool_calls > 0:
        clean_pts -= 2
    clean_pts = max(0, clean_pts)

    components = {
        "success_pts": success_pts,
        "routing_pts": routing_pts,
        "speed_pts": speed_pts,
        "cost_pts": cost_pts,
        "clean_pts": clean_pts,
    }
    quality = max(0, min(100, sum(components.values())))

    # A passing task needs only a glimpse; a FAILING one has to be
    # diagnosable from the record alone — 400 chars regularly cut off the
    # sentence that tripped a signal, leaving the failure unexplainable
    # without re-running.
    _ex_cap = 400 if (success and not violated) else 4000
    excerpt = (traj.text or "")[:_ex_cap].strip()
    tool_names = [str(c.get("name", "")) for c in traj.tool_calls
                  if c.get("name")]

    return BenchmarkResult(
        task_id=task.id,
        task_class=task.task_class,
        model=model,
        profile_name=profile_name,
        mode=task.mode,
        ts=ts if ts is not None else time.time(),
        success=success,
        quality_0_100=quality,
        components=components,
        duration_s=traj.duration_s,
        cost_usd=traj.cost_usd,
        input_tokens=traj.input_tokens,
        output_tokens=traj.output_tokens,
        tool_calls=len(traj.tool_calls),
        matched_signals=matched,
        violated_signals=violated,
        missing_signals=missing,
        budget_violations=budget_violations,
        error=traj.error,
        unmeasured=unmeasured,
        text_excerpt=excerpt,
        tool_names=tool_names,
        behavior=behavior_flags(task, traj),
        caveats=caveat_count(traj.text),
        answer_chars=len(str(traj.text or "")),
        denials=traj.denials,
    )


# ---------------------------------------------------------------------------
# Replicate aggregation
# ---------------------------------------------------------------------------


def _median(values: list[float]) -> float:
    """Numerically-stable median for small lists.  Empty → 0.0."""
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2 == 0:
        return (s[mid - 1] + s[mid]) / 2.0
    return float(s[mid])


def _stdev(values: list[float]) -> float:
    """Sample standard deviation (Bessel's correction).  Empty/single → 0."""
    n = len(values)
    if n < 2:
        return 0.0
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    return var ** 0.5


def aggregate_replicates(
    results: list[BenchmarkResult],
) -> BenchmarkResult:
    """Collapse N per-run results for the SAME task into one aggregated
    BenchmarkResult.

    Rules:
      - ``success`` = majority vote (ceil(N/2) successes required to pass)
      - ``quality_0_100`` = median across N runs
      - ``duration_s`` / ``tool_calls`` / ``input_tokens`` / ``output_tokens``
        = median (typical observation, not skewed by one outlier)
      - ``cost_usd`` = SUM (real total spend across replicates)
      - ``components`` = median per component
      - ``matched_signals`` / ``missing_signals`` / ``violated_signals`` /
        ``budget_violations`` = union across replicates (so a flaky
        violation is still visible)
      - ``quality_stdev`` exposes run-to-run noise; ``success_rate``
        the fraction of passes; ``per_run_*`` keep the raw samples for
        deeper analysis.

    UNMEASURED replicates are excluded from every one of those, and the
    aggregate is unmeasured only when none of them measured anything.
    ``is_unmeasured`` guarded the single-run path and this one dropped
    the flag, so with ``--repeats`` the whole mechanism was bypassed:
    observed live, an engine that never started -- a missing provider
    argument -- came back 0/11 at quality 35 with no NOT MEASURED notice
    at all, and would have been written to the file baselines compare
    against. One hop in the middle is all a guard needs to lose.

    A partial outage is the more common shape: two replicates ran, the
    third hit an endpoint with no capacity. Medianing that third one's
    zero in punishes the model for the network.

    Raises ``ValueError`` on empty input or task_id mismatch.
    """
    if not results:
        raise ValueError("aggregate_replicates needs at least one result")
    first = results[0]
    if not all(r.task_id == first.task_id for r in results):
        raise ValueError("aggregate_replicates: all results must share task_id")

    all_unmeasured = all(r.unmeasured for r in results)
    if not all_unmeasured:
        results = [r for r in results if not r.unmeasured]
    n = len(results)
    qualities = [int(r.quality_0_100) for r in results]
    durations = [float(r.duration_s) for r in results]
    tool_counts = [int(r.tool_calls) for r in results]
    in_toks = [int(r.input_tokens) for r in results]
    out_toks = [int(r.output_tokens) for r in results]
    success_flags = [bool(r.success) for r in results]
    n_pass = sum(1 for f in success_flags if f)

    # Components: median per field across replicates
    comp_keys: set[str] = set()
    for r in results:
        comp_keys.update((r.components or {}).keys())
    components = {
        k: int(round(_median([float((r.components or {}).get(k, 0)) for r in results])))
        for k in comp_keys
    }

    # Union of signal labels (preserves order of first occurrence)
    def _union(field_name: str) -> list[str]:
        out: list[str] = []
        for r in results:
            for s in getattr(r, field_name, []) or []:
                if s not in out:
                    out.append(s)
        return out

    def _matched_every_time(missing: list[str]) -> list[str]:
        """Signals that matched in EVERY scored replicate.

        The union was applied to matched and missing alike, so a signal
        that matched twice out of three came back in BOTH lists. Observed
        live on workflow_plan_before_act: four signals listed as matched
        and as missing at once, and the audit printed five missing
        patterns where exactly one was genuinely absent — which is how a
        reader gets sent to the wrong defect. (It sent me there.)

        The union stays right for ``missing`` and ``violated``: a signal
        the model dropped even once is not something the run can vouch
        for. So ``matched`` becomes the complement, and the two lists are
        disjoint by construction.
        """
        gone = set(missing)
        return [s for s in _union("matched_signals") if s not in gone]

    def _flaky(missing: list[str]) -> list[str]:
        """Matched sometimes and missed sometimes — the run-to-run noise
        that "missing" alone flattens into a verdict."""
        gone = set(missing)
        return [s for s in _union("matched_signals") if s in gone]

    # Pick the first non-empty excerpt — gives forensic value without
    # bloating storage with N copies; tool_names unioned the same way
    # as signals (a flaky tool-call still surfaces).
    excerpt = ""
    for r in results:
        if r.text_excerpt:
            excerpt = r.text_excerpt
            break
    tool_names_union: list[str] = []
    for r in results:
        for n_name in (r.tool_names or []):
            if n_name and n_name not in tool_names_union:
                tool_names_union.append(n_name)

    # Behaviour flags: mean per flag across replicates → a 0..1 rate.
    beh_sums: dict[str, float] = {}
    beh_counts: dict[str, int] = {}
    for r in results:
        for k, v in (r.behavior or {}).items():
            if v is None:
                continue
            try:
                beh_sums[k] = beh_sums.get(k, 0.0) + float(v)
                beh_counts[k] = beh_counts.get(k, 0) + 1
            except (TypeError, ValueError):
                continue
    behavior_agg = {
        k: beh_sums[k] / beh_counts[k] for k in beh_sums if beh_counts[k]
    }

    return BenchmarkResult(
        task_id=first.task_id,
        task_class=first.task_class,
        model=first.model,
        profile_name=first.profile_name,
        mode=first.mode,
        ts=first.ts,
        success=(n_pass * 2 >= n),                   # majority (ties = pass)
        quality_0_100=int(round(_median([float(q) for q in qualities]))),
        components=components,
        duration_s=_median(durations),
        cost_usd=sum(float(r.cost_usd) for r in results),
        input_tokens=int(_median([float(x) for x in in_toks])),
        output_tokens=int(_median([float(x) for x in out_toks])),
        tool_calls=int(_median([float(x) for x in tool_counts])),
        matched_signals=_matched_every_time(_union("missing_signals")),
        violated_signals=_union("violated_signals"),
        missing_signals=_union("missing_signals"),
        flaky_signals=_flaky(_union("missing_signals")),
        budget_violations=_union("budget_violations"),
        error="; ".join(sorted({r.error for r in results if r.error}))[:500],
        n_samples=n,
        quality_stdev=_stdev([float(q) for q in qualities]),
        success_rate=n_pass / n,
        per_run_quality=list(qualities),
        per_run_success=list(success_flags),
        text_excerpt=excerpt,
        tool_names=tool_names_union,
        behavior=behavior_agg,
        unmeasured=all_unmeasured,
        caveats=int(_median([float(r.caveats) for r in results])),
        answer_chars=int(_median([float(r.answer_chars) for r in results])),
        # Unobserved stays unobserved: a median over a list with holes in
        # it would report a number for samples that never looked.
        denials=(None if any(r.denials is None for r in results)
                 else int(_median([float(r.denials or 0) for r in results]))),
    )


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------


def _slug(s: str) -> str:
    out = "".join(c if c.isalnum() or c in "-._" else "_" for c in s)
    return out[:64] or "unknown"


def write_run(
    results: Iterable[BenchmarkResult],
    *,
    model: str,
    runs_dir: Path | None = None,
    run_id: str | None = None,
) -> Path:
    """Persist all results from one benchmark run as JSONL."""

    d = runs_dir or _DEFAULT_RUNS_DIR
    d.mkdir(parents=True, exist_ok=True)
    ts = int(time.time())
    rid = run_id or f"{ts}_{_slug(model)}"
    path = d / f"{rid}.jsonl"
    with path.open("w", encoding="utf-8") as f:
        for r in results:
            f.write(json.dumps(asdict(r), ensure_ascii=False) + "\n")
    return path


def read_run(path: Path) -> list[dict]:
    if not path.exists():
        return []
    out: list[dict] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return out


# ---------------------------------------------------------------------------
# Aggregation + comparison
# ---------------------------------------------------------------------------


def summarise_run(results: list[dict] | list[BenchmarkResult]) -> dict[str, Any]:
    """Roll-up across a single run: pass-rate + avg quality + totals."""

    rows = [
        asdict(r) if isinstance(r, BenchmarkResult) else r for r in results
    ]
    if not rows:
        return {
            "n_tasks": 0, "n_pass": 0, "n_fail": 0, "pass_rate": 0.0,
            "n_unmeasured": 0, "unmeasured_tasks": [],
            "avg_quality": 0.0, "total_cost_usd": 0.0,
            "total_duration_s": 0.0, "total_tool_calls": 0,
            "avg_caveats": 0.0, "max_caveats": 0,
            "avg_output_tokens": 0.0, "avg_answer_chars": 0.0,
            "total_denials": None,
        }
    # A task whose request never reached the model is excluded from both
    # the rate and the average, and counted separately. Including it
    # writes an endpoint outage into the number a baseline is compared
    # against -- which lowers the bar permanently, and silently, so the
    # next real regression passes under it. Observed 2026-08-12: two
    # tasks whose entire output was the harness's own retry banners took
    # a suite from 9/11 to 7/11.
    scored = [r for r in rows if not r.get("unmeasured")]
    n_unmeasured = len(rows) - len(scored)
    n = len(scored)
    n_pass = sum(1 for r in scored if r.get("success"))
    return {
        "n_tasks": n,
        "n_pass": n_pass,
        "n_unmeasured": n_unmeasured,
        "unmeasured_tasks": sorted(
            str(r.get("task_id") or "?") for r in rows if r.get("unmeasured")),
        "pass_rate": (n_pass / n) if n else 0.0,
        "avg_quality": (sum(int(r.get("quality_0_100") or 0)
                            for r in scored) / n) if n else 0.0,
        "total_cost_usd": sum(float(r.get("cost_usd") or 0) for r in rows),
        "total_duration_s": sum(float(r.get("duration_s") or 0) for r in rows),
        "total_tool_calls": sum(int(r.get("tool_calls") or 0) for r in rows),
        # The cost side. These do not judge the run, they are what a
        # later run is compared against: a suite can hold 11/11 while
        # every answer grows another hedge, and only these move.
        "avg_caveats": (sum(int(r.get("caveats") or 0)
                            for r in scored) / n) if n else 0.0,
        "max_caveats": max((int(r.get("caveats") or 0) for r in scored),
                           default=0),
        "avg_output_tokens": (sum(int(r.get("output_tokens") or 0)
                                  for r in scored) / n) if n else 0.0,
        "avg_answer_chars": (sum(int(r.get("answer_chars") or 0)
                                 for r in scored) / n) if n else 0.0,
        # Kept apart from the rate: a task that fails IS a cost signal,
        # but a task nobody could measure is not.
        "n_fail": n - n_pass,
        # None when any scored task did not observe denials at all. A
        # zero here would read as "nothing was refused", which is a
        # different sentence from "nobody looked".
        "total_denials": (
            None if any(r.get("denials") is None for r in scored)
            else sum(int(r.get("denials") or 0) for r in scored)),
    }


def compare_runs(
    baseline: list[dict] | Path,
    candidate: list[dict] | Path,
) -> dict[str, Any]:
    """Per-task delta + roll-up.  Used to answer "did the knob help?"

    Direction conventions:
      - quality_0_100 / success → higher is better
      - cost_usd / duration_s / tool_calls → lower is better

    Returns ``{"per_task": [...], "summary": {...}, "verdict": "..."}``.
    Verdict is one of {"better", "worse", "mixed", "neutral", "thin"} where
    ``thin`` means fewer than 3 overlapping tasks (too little signal).
    """

    if isinstance(baseline, Path):
        baseline = read_run(baseline)
    if isinstance(candidate, Path):
        candidate = read_run(candidate)

    by_id_old = {r.get("task_id"): r for r in baseline}
    by_id_new = {r.get("task_id"): r for r in candidate}
    overlap = sorted(set(by_id_old) & set(by_id_new))

    per_task: list[dict] = []
    n_better = n_worse = n_neutral = 0
    for tid in overlap:
        o = by_id_old[tid]
        n = by_id_new[tid]
        d_quality = int(n.get("quality_0_100") or 0) - int(o.get("quality_0_100") or 0)
        d_cost = float(n.get("cost_usd") or 0) - float(o.get("cost_usd") or 0)
        d_dur = float(n.get("duration_s") or 0) - float(o.get("duration_s") or 0)
        d_tools = int(n.get("tool_calls") or 0) - int(o.get("tool_calls") or 0)
        # Classification per task:
        # better if quality strictly up OR (quality flat AND cost or duration
        # meaningfully down ≥10% or absolute 1s/1¢).
        improved = False
        regressed = False
        if d_quality > 0:
            improved = True
        elif d_quality < 0:
            regressed = True
        else:
            base_cost = float(o.get("cost_usd") or 0)
            base_dur = float(o.get("duration_s") or 0)
            cost_drop = d_cost < -0.01 or (
                base_cost > 0 and d_cost / base_cost <= -0.10
            )
            dur_drop = d_dur < -1.0 or (
                base_dur > 0 and d_dur / base_dur <= -0.10
            )
            cost_rise = d_cost > 0.01 or (
                base_cost > 0 and d_cost / base_cost >= 0.10
            )
            dur_rise = d_dur > 1.0 or (
                base_dur > 0 and d_dur / base_dur >= 0.10
            )
            if cost_drop or dur_drop:
                improved = True
            elif cost_rise or dur_rise:
                regressed = True
        if improved and not regressed:
            cls = "better"
            n_better += 1
        elif regressed and not improved:
            cls = "worse"
            n_worse += 1
        else:
            cls = "neutral"
            n_neutral += 1
        per_task.append({
            "task_id": tid,
            "class": cls,
            "old_quality": int(o.get("quality_0_100") or 0),
            "new_quality": int(n.get("quality_0_100") or 0),
            "d_quality": d_quality,
            "d_cost_usd": d_cost,
            "d_duration_s": d_dur,
            "d_tool_calls": d_tools,
            "old_success": bool(o.get("success")),
            "new_success": bool(n.get("success")),
        })

    summary_old = summarise_run(baseline)
    summary_new = summarise_run(candidate)
    summary = {
        "n_overlap": len(overlap),
        "old": summary_old,
        "new": summary_new,
        "n_better": n_better,
        "n_worse": n_worse,
        "n_neutral": n_neutral,
    }

    if len(overlap) < 3:
        verdict = "thin"
    elif n_better >= max(3, int(0.8 * len(overlap))) and n_worse == 0:
        verdict = "better"
    elif n_worse >= max(3, int(0.5 * len(overlap))):
        verdict = "worse"
    elif n_better > n_worse + 1:
        verdict = "better"
    elif n_worse > n_better + 1:
        verdict = "worse"
    elif n_better > 0 and n_worse > 0:
        verdict = "mixed"
    else:
        verdict = "neutral"

    return {"per_task": per_task, "summary": summary, "verdict": verdict}


# ---------------------------------------------------------------------------
# Markdown export + profile-commit linking
# ---------------------------------------------------------------------------


def runs_dir() -> Path:
    """Public accessor for the default runs directory."""
    return _DEFAULT_RUNS_DIR


def run_timestamp(path: Path) -> float:
    """Return the ts of the first record (≈ run start) or 0 if unreadable.

    Used to bracket ``git log`` so we can auto-find which profile commit
    sits between a baseline and a candidate run.
    """
    records = read_run(path)
    if records:
        try:
            return float(records[0].get("ts") or 0)
        except (TypeError, ValueError):
            return 0.0
    return 0.0


def find_profile_commits_between(
    old_ts: float,
    new_ts: float,
    *,
    profile_file: str = "delfin/agent/model_profiles.py",
    repo_root: str | Path | None = None,
    timeout_s: float = 5.0,
) -> list[str]:
    """Return ``<short-hash> <subject>`` lines for commits that touched
    the profile file in the (old_ts, new_ts] window.

    Used to auto-annotate compare reports so we know which knob change
    drove the observed Δ.  Returns ``[]`` on any git failure (missing
    repo, no commits, timeout) — never raises.
    """

    if old_ts <= 0 or new_ts <= 0 or new_ts <= old_ts:
        return []
    cwd = str(repo_root) if repo_root else os.getcwd()
    # Pad the window by 1 s on each side: git stores commit timestamps at
    # second precision, so a commit at exactly `new_ts` would otherwise
    # be excluded by `--until` (which is strict "before").
    try:
        old_iso = datetime.fromtimestamp(max(0.0, old_ts - 1.0)).isoformat()
        new_iso = datetime.fromtimestamp(new_ts + 1.0).isoformat()
    except (OSError, ValueError, OverflowError):
        return []
    try:
        out = subprocess.run(
            [
                "git", "log", "--oneline",
                f"--since={old_iso}", f"--until={new_iso}",
                "--", profile_file,
            ],
            capture_output=True, text=True, timeout=timeout_s, cwd=cwd,
        )
    except (subprocess.SubprocessError, OSError):
        return []
    if out.returncode != 0:
        return []
    return [line.strip() for line in out.stdout.splitlines() if line.strip()]


_VERDICT_EMOJI = {
    "better": "[+]", "worse": "[-]", "mixed": "[~]",
    "neutral": "[=]", "thin": "[?]",
}
_CLASS_MARK = {"better": "[+]", "worse": "[-]", "neutral": "[=]"}


def format_compare_markdown(
    cmp_result: dict,
    *,
    baseline_path: Path | None = None,
    candidate_path: Path | None = None,
    include_profile_commits: bool = True,
    repo_root: str | Path | None = None,
) -> str:
    """Render a ``compare_runs`` result as a markdown report.

    Suitable for pasting into a PR body, memory entry, or chat message.
    If both baseline_path and candidate_path are supplied AND
    ``include_profile_commits=True``, the report also lists any commits
    that touched ``delfin/agent/model_profiles.py`` between the two run
    timestamps — the "which knob change drove this delta?" annotation.
    """

    summary = cmp_result.get("summary") or {}
    verdict = str(cmp_result.get("verdict") or "neutral")
    old = summary.get("old") or {}
    new = summary.get("new") or {}

    lines: list[str] = ["# Benchmark Comparison Report", ""]
    lines.append(
        f"**Verdict: {verdict.upper()}** "
        f"{_VERDICT_EMOJI.get(verdict, '[?]')}"
    )
    lines.append(
        f"  {summary.get('n_better', 0)} better, "
        f"{summary.get('n_worse', 0)} worse, "
        f"{summary.get('n_neutral', 0)} neutral  "
        f"(n={summary.get('n_overlap', 0)})"
    )
    if baseline_path is not None:
        lines.append(f"  Baseline:  `{baseline_path}`")
    if candidate_path is not None:
        lines.append(f"  Candidate: `{candidate_path}`")
    lines.append("")

    # Summary table
    lines.append("## Summary")
    lines.append("")
    lines.append("|              | Baseline   | Candidate  | Δ          |")
    lines.append("|--------------|------------|------------|------------|")
    op = float(old.get("pass_rate") or 0)
    np_ = float(new.get("pass_rate") or 0)
    lines.append(
        f"| Pass rate    | {op:>9.0%}  | {np_:>9.0%}  | "
        f"{(np_ - op) * 100:+8.0f}pp |"
    )
    oq = float(old.get("avg_quality") or 0)
    nq = float(new.get("avg_quality") or 0)
    lines.append(
        f"| Avg quality  | {oq:>10.1f} | {nq:>10.1f} | {nq - oq:+10.1f} |"
    )
    oc = float(old.get("total_cost_usd") or 0)
    nc = float(new.get("total_cost_usd") or 0)
    lines.append(
        f"| Total cost   | ${oc:>9.4f} | ${nc:>9.4f} | "
        f"${nc - oc:+9.4f} |"
    )
    od = float(old.get("total_duration_s") or 0)
    nd = float(new.get("total_duration_s") or 0)
    lines.append(
        f"| Total time   | {od:>9.1f}s | {nd:>9.1f}s | {nd - od:+9.1f}s |"
    )
    lines.append("")

    # Per-task table
    per_task = cmp_result.get("per_task") or []
    if per_task:
        lines.append("## Per-task")
        lines.append("")
        lines.append(
            "| Task                          | Class   | Quality | Δq   | "
            "Δcost     | Δdur      |"
        )
        lines.append(
            "|-------------------------------|---------|---------|------|"
            "-----------|-----------|"
        )
        for row in per_task:
            mark = _CLASS_MARK.get(row.get("class") or "", "[?]")
            tid = str(row.get("task_id") or "")[:30]
            oq2 = int(row.get("old_quality") or 0)
            nq2 = int(row.get("new_quality") or 0)
            dq = int(row.get("d_quality") or 0)
            dc = float(row.get("d_cost_usd") or 0)
            dd = float(row.get("d_duration_s") or 0)
            lines.append(
                f"| `{tid:<29}` | {mark:<7} | {oq2:>3}→{nq2:<3} | "
                f"{dq:+5d} | {dc:+9.4f} | {dd:+8.2f}s |"
            )
        lines.append("")

    # Profile-commit annotation
    if (include_profile_commits and baseline_path is not None
            and candidate_path is not None):
        try:
            old_ts = run_timestamp(baseline_path)
            new_ts = run_timestamp(candidate_path)
            commits = find_profile_commits_between(
                old_ts, new_ts, repo_root=repo_root,
            )
        except Exception:
            commits = []
        if commits:
            lines.append("## Profile changes between runs")
            lines.append("")
            lines.append(
                "Commits that modified "
                "`delfin/agent/model_profiles.py`:"
            )
            lines.append("")
            for c in commits:
                lines.append(f"- `{c}`")
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


# ---------------------------------------------------------------------------
# A/B convenience — compact drift view over two runs
# ---------------------------------------------------------------------------


def list_runs(runs_dir: Path | None = None) -> list[Path]:
    """All persisted benchmark run files, oldest first (mtime, then name).

    Convenience for "compare the two most recent runs" flows.  Never
    raises — a missing/unreadable directory yields ``[]``.
    """
    base = runs_dir or _DEFAULT_RUNS_DIR
    try:
        files = [p for p in base.glob("*.jsonl") if p.is_file()]
    except Exception:
        return []

    def _key(p: Path) -> tuple[float, str]:
        try:
            return (p.stat().st_mtime, p.name)
        except OSError:
            return (0.0, p.name)

    return sorted(files, key=_key)


def _resolve_run_path(run: str, base: Path) -> Path:
    """Accept a full path, a filename in ``base``, or a bare run id."""
    p = Path(run)
    if p.exists():
        return p
    cand = base / run
    if cand.exists():
        return cand
    return base / f"{run}.jsonl"


def ab_compare(
    run_a: str,
    run_b: str,
    *,
    runs_dir: Path | None = None,
) -> dict[str, Any]:
    """Compact A/B view over two persisted runs (A = baseline, B = candidate).

    ``run_a`` / ``run_b`` may be full paths, filenames, or bare run ids
    resolved against the runs directory.  Built on :func:`compare_runs`
    plus :func:`behavior_rates`; missing files behave like empty runs
    (verdict ``thin``).  Returns::

        {run_a, run_b, verdict, n_overlap, n_better, n_worse, n_neutral,
         per_task_delta:  {task_id: {d_quality, d_cost_usd, class}},
         score_delta:     candidate avg_quality − baseline avg_quality,
         cost_delta:      candidate total cost  − baseline total cost,
         behaviour_flag_changes: {flag: {old, new, delta}}}   # changed only
    """
    base = runs_dir or _DEFAULT_RUNS_DIR
    path_a = _resolve_run_path(run_a, base)
    path_b = _resolve_run_path(run_b, base)
    rec_a = read_run(path_a)
    rec_b = read_run(path_b)
    cmp_result = compare_runs(rec_a, rec_b)
    summary = cmp_result["summary"]

    per_task_delta = {
        row["task_id"]: {
            "d_quality": row["d_quality"],
            "d_cost_usd": row["d_cost_usd"],
            "class": row["class"],
        }
        for row in cmp_result["per_task"]
    }

    rates_a = behavior_rates(rec_a)
    rates_b = behavior_rates(rec_b)
    flag_changes: dict[str, dict[str, float]] = {}
    for flag in sorted(set(rates_a) | set(rates_b)):
        old = float(rates_a.get(flag, {}).get("rate", 0.0))
        new = float(rates_b.get(flag, {}).get("rate", 0.0))
        if abs(new - old) > 1e-9:
            flag_changes[flag] = {"old": old, "new": new, "delta": new - old}

    return {
        "run_a": str(path_a),
        "run_b": str(path_b),
        "verdict": cmp_result["verdict"],
        "n_overlap": summary["n_overlap"],
        "n_better": summary["n_better"],
        "n_worse": summary["n_worse"],
        "n_neutral": summary["n_neutral"],
        "per_task_delta": per_task_delta,
        "score_delta": (float(summary["new"]["avg_quality"])
                        - float(summary["old"]["avg_quality"])),
        "cost_delta": (float(summary["new"]["total_cost_usd"])
                       - float(summary["old"]["total_cost_usd"])),
        "behaviour_flag_changes": flag_changes,
    }


def format_ab_note(ab: dict) -> str:
    """Short plain-text rendering of an :func:`ab_compare` result — the
    "Benchmark drift" note for the eval report (and chat)."""
    verdict = str(ab.get("verdict") or "neutral")
    lines = [
        f"verdict: {verdict.upper()} "
        f"({int(ab.get('n_better') or 0)} better / "
        f"{int(ab.get('n_worse') or 0)} worse / "
        f"{int(ab.get('n_neutral') or 0)} neutral, "
        f"n={int(ab.get('n_overlap') or 0)})",
        f"score delta: {float(ab.get('score_delta') or 0.0):+.1f} avg quality; "
        f"cost delta: ${float(ab.get('cost_delta') or 0.0):+.4f}",
    ]
    flags = ab.get("behaviour_flag_changes") or {}
    if flags:
        parts = [
            f"{flag} {float(v.get('old') or 0.0):.0%}→"
            f"{float(v.get('new') or 0.0):.0%}"
            for flag, v in sorted(flags.items())
        ]
        lines.append("behaviour: " + ", ".join(parts))
    regressed = sorted(
        tid for tid, d in (ab.get("per_task_delta") or {}).items()
        if d.get("class") == "worse"
    )
    if regressed:
        head = ", ".join(regressed[:5])
        more = f" (+{len(regressed) - 5} more)" if len(regressed) > 5 else ""
        lines.append(f"regressed: {head}{more}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Audit — pattern-bug-vs-real-fail diagnosis
# ---------------------------------------------------------------------------


def audit_run(
    records: list[dict] | Path,
    *,
    tasks: list[Task] | None = None,
) -> list[dict]:
    """Return one audit entry per FAILED task in a run.

    A failure is anything with ``success == False``.  For each, the
    entry surfaces enough state to decide pattern-bug-vs-real-fail in
    one glance:

      - task_id / model
      - quality / success_rate / quality_stdev (low σ + low rate
        = deterministic; high σ = flaky)
      - missing / violated signal labels (which check failed)
      - text_excerpt (what the model actually wrote)
      - tool_names (which tools the model actually called)
      - hint_pattern_bug: True if the model output looks reasonable
        (excerpt is non-empty AND no error AND quality_stdev < 5),
        which means the signals are likely the problem, not the model.
    """
    if isinstance(records, Path):
        records = read_run(records)
    task_by_id = {t.id: t for t in (tasks or [])}
    if not task_by_id:
        try:
            task_by_id = {t.id: t for t in load_tasks()}
        except Exception:
            task_by_id = {}

    out: list[dict] = []
    for r in records:
        if r.get("success"):
            continue
        rate = float(r.get("success_rate") or 0)
        stdev = float(r.get("quality_stdev") or 0)
        excerpt = str(r.get("text_excerpt") or "")
        error = str(r.get("error") or "")
        violated = list(r.get("violated_signals") or [])
        # Pattern-bug heuristic: deterministic fail + non-empty reasonable-
        # looking output + no engine error + NO forbidden violation.  A
        # violated forbidden signal means the model did something explicitly
        # disallowed → that's real misbehaviour, never a pattern bug.
        hint_pattern_bug = (
            rate <= 0.34
            and stdev < 5.0
            and len(excerpt) >= 30
            and not error
            and not violated
        )
        task = task_by_id.get(r.get("task_id") or "")
        signal_defs: dict[str, str] = {}
        if task is not None:
            for idx, sig in enumerate(task.expected_signals):
                signal_defs[f"{task.id}.expected[{idx}]"] = (
                    f"{sig.pattern}   (against={sig.against})"
                )
            for idx, sig in enumerate(task.forbidden_signals):
                signal_defs[f"{task.id}.forbidden[{idx}]"] = (
                    f"{sig.pattern}   (against={sig.against})"
                )
        out.append({
            "task_id": r.get("task_id") or "",
            "model": r.get("model") or "",
            "quality": int(r.get("quality_0_100") or 0),
            "success_rate": rate,
            "quality_stdev": stdev,
            "missing_signals": list(r.get("missing_signals") or []),
            "violated_signals": violated,
            "tool_names": list(r.get("tool_names") or []),
            "text_excerpt": excerpt,
            "error": error,
            "signal_defs": signal_defs,
            "hint_pattern_bug": hint_pattern_bug,
        })
    return out


def format_audit_report(entries: list[dict]) -> str:
    """Render audit entries as a developer-friendly text report."""
    if not entries:
        return "No failed tasks — nothing to audit.\n"
    lines: list[str] = []
    bug_entries = [e for e in entries if e["hint_pattern_bug"]]
    real_entries = [e for e in entries if not e["hint_pattern_bug"]]
    lines.append(f"=== AUDIT: {len(entries)} failed task(s) "
                 f"({len(bug_entries)} likely pattern-bug, "
                 f"{len(real_entries)} likely real-fail) ===\n")
    for group_name, group in (
        ("SUSPECTED PATTERN BUG", bug_entries),
        ("LIKELY REAL FAIL OR FLAKY", real_entries),
    ):
        if not group:
            continue
        lines.append(f"--- {group_name} ---")
        for e in group:
            lines.append(f"\n  task   : {e['task_id']}    (model={e['model']})")
            lines.append(f"  quality: q={e['quality']}  "
                         f"rate={e['success_rate']:.0%}  "
                         f"σ={e['quality_stdev']:.1f}")
            if e["tool_names"]:
                lines.append(f"  tools  : {', '.join(e['tool_names'])}")
            if e["missing_signals"]:
                lines.append("  missing signals:")
                for label in e["missing_signals"]:
                    if label.endswith(":optional"):
                        continue
                    defn = e["signal_defs"].get(label, "(definition not found)")
                    lines.append(f"    {label}")
                    lines.append(f"      pattern: {defn}")
            if e["violated_signals"]:
                lines.append("  VIOLATED signals:")
                for label in e["violated_signals"]:
                    defn = e["signal_defs"].get(label, "(definition not found)")
                    lines.append(f"    {label}")
                    lines.append(f"      pattern: {defn}")
            if e["error"]:
                lines.append(f"  error  : {e['error'][:200]}")
            if e["text_excerpt"]:
                lines.append("  model output (≤400 chars):")
                for line in e["text_excerpt"].splitlines()[:8]:
                    lines.append(f"    │ {line[:120]}")
                if len(e["text_excerpt"].splitlines()) > 8:
                    lines.append("    │ …")
        lines.append("")
    return "\n".join(lines) + "\n"


__all__ = [
    "Signal",
    "Task",
    "Trajectory",
    "BenchmarkResult",
    "load_tasks",
    "score_outcome",
    "behavior_flags",
    "behavior_rates",
    "write_run",
    "read_run",
    "summarise_run",
    "compare_runs",
    "aggregate_replicates",
    "list_runs",
    "ab_compare",
    "format_ab_note",
    "runs_dir",
    "run_timestamp",
    "find_profile_commits_between",
    "format_compare_markdown",
    "audit_run",
    "format_audit_report",
]
