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
import math
import os
import re
import sys
import tempfile
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Optional

from . import pricing as _pricing
from . import safe_git
from . import verify_commands as _vc

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
class RubricCriterion:
    """One weighted, deterministic check in an OPTIONAL per-task rubric.

    A rubric sits beside the signal/verify rubric that decides ``success``;
    it never rewrites it.  ``kind`` names the deterministic check, all of
    them computed from the same recorded trajectory ``score_outcome`` reads:

    - ``signals``      fraction of required expected signals matched
    - ``values``       fraction of expected figures matched (with the
                       ``ExpectedValue`` tolerance -- numbers, never strings)
    - ``verify``       1.0 when the acceptance script passed, else 0.0
    - ``forbidden``    1.0 when no forbidden signal was violated
    - ``clean``        1.0 when the run carried no error
    - ``behavior.<f>`` the behavior flag ``f`` (0 or 1)
    - ``process.<k>``  the process rate ``k``, clamped to [0, 1]

    ``weight`` scales the criterion's share of the summed 0-100 rubric
    score.  An unknown ``kind`` is a task-authoring bug and raises
    ``ValueError`` rather than scoring a silent zero.
    """

    name: str
    weight: float = 1.0
    kind: str = "signals"


# A number as it is actually written in an answer: German and English
# grouping, a bare decimal, a signed value, a thousands run.
#
# The hyphen is excluded on the left as well as the word characters. An
# administrative answer is full of record references -- R-014, A-2026-003,
# a date written 2026-03-19 -- and without it the digits after the dash
# read as their own figure, so an answer that merely NAMES beleg R-014
# would satisfy an expected value of 14. A leading minus still works
# wherever a sign can actually stand: at the start, or after a space or a
# bracket, never glued to the end of a word.
_NUMBER_TOKEN_RE = re.compile(
    r"(?<![\w.,-])[+-]?\d{1,3}(?:[.,]\d{3})*(?:[.,]\d+)?(?![\w])"
    r"|(?<![\w.,-])[+-]?\d+(?:[.,]\d+)?(?![\w])")


def readings_of(token: str) -> tuple[float, ...]:
    """Every number this token could denote, most likely first.

    A written figure is not a float. "6.070,55" and "6,070.55" are the
    same amount in two conventions, "6070.55" is a third spelling of it,
    and "1.265" is either a thousand and a bit or one and a quarter with
    nothing in the string to decide. So a token becomes a SET of readings
    and a comparison succeeds when any of them fits.

    Deliberately permissive, and only safe because of what it is used for:
    the question asked is "did the expected figure appear", where reading
    an ambiguous token both ways can only fail to catch a wrong answer. A
    stricter reading would instead fail a correct one for writing its
    number the way the rest of the document writes it, which is the error
    that has actually been made here, repeatedly, by patterns matching
    digits.
    """
    text = str(token or "").strip()
    if not text:
        return ()
    sign = -1.0 if text.startswith("-") else 1.0
    body = text.lstrip("+-")
    if not body or not body[0].isdigit():
        return ()

    out: list[float] = []

    def _add(value: Optional[float]) -> None:
        if value is not None and value not in out:
            out.append(value)

    def _read(decimal: str, group: str) -> Optional[str]:
        """The plain form under one convention, or None if it cannot hold.

        Grouping is checked on the INTEGER part only, and every group
        after the first must be exactly three digits — that is what makes
        "1,5" refuse to be fifteen while "1,265" is allowed to be a
        thousand two hundred and sixty-five.
        """
        if decimal:
            if body.count(decimal) != 1:
                return None
            integer, fraction = body.split(decimal)
            if not fraction.isdigit():
                return None
        else:
            integer, fraction = body, ""
        if group:
            if group in fraction:
                return None
            parts = integer.split(group)
            if len(parts) > 1:
                if not parts[0] or len(parts[0]) > 3:
                    return None
                if any(len(part) != 3 for part in parts[1:]):
                    return None
            integer = integer.replace(group, "")
        elif any(sep in integer for sep in ".,"):
            return None
        if not integer.isdigit():
            return None
        return integer + ("." + fraction if fraction else "")

    has_dot, has_comma = "." in body, "," in body
    if has_dot and has_comma:
        # The LAST separator is the decimal one; the other groups.
        if body.rindex(",") > body.rindex("."):
            _add(_to_float(_read(",", ".")))
        else:
            _add(_to_float(_read(".", ",")))
    elif has_comma:
        _add(_to_float(_read(",", "")))          # 1,5  -> 1.5
        _add(_to_float(_read("", ",")))          # 1,265 -> 1265
    elif has_dot:
        _add(_to_float(_read(".", "")))          # 1.5  -> 1.5
        _add(_to_float(_read("", ".")))          # 1.265 -> 1265
    else:
        _add(_to_float(body))
    return tuple(sign * v for v in out)


def _to_float(text: Optional[str]) -> Optional[float]:
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _short_number(value: float) -> str:
    """A figure written the way it would be read back.

    Only for the diagnostic note beside a "wrong" verdict, so trailing
    zeros are noise: 6.4 is a version string and 6.400000 is not more
    informative about that.
    """
    if value == int(value) and abs(value) < 1e15:
        return str(int(value))
    return f"{value:.6g}"


# The signs a model writes. GLM sets the typographic minus (U+2212) before
# every negative energy; so does any text that passed through a word
# processor. A rubric written as "-613\\.417" and a value read by the
# number tokeniser both want the ASCII hyphen-minus, and for one night
# five correct answers on the audit task scored 0/5 with sigma 0.9 --
# the signature of the instrument, not the model (2026-09-11). Every
# haystack and every number is read with these folded to "-" first.
_MINUS_SIGNS = str.maketrans({
    "\u2212": "-",   # minus sign
    "\u2010": "-", "\u2011": "-", "\u2012": "-",   # hyphens, figure dash
    "\u2013": "-", "\u2014": "-", "\u2015": "-",   # en/em dash, horizontal bar
    "\ufe63": "-", "\uff0d": "-",                   # small / fullwidth hyphen-minus
})


def ascii_minus(text: str) -> str:
    """*text* with every dash and minus sign a model might write folded
    to the ASCII hyphen-minus a pattern or a number expects."""
    return str(text or "").translate(_MINUS_SIGNS)


def numbers_in(text: str) -> tuple[float, ...]:
    """Every value any number in *text* could denote."""
    seen: list[float] = []
    for match in _NUMBER_TOKEN_RE.finditer(ascii_minus(text)):
        for value in readings_of(match.group(0)):
            if value not in seen:
                seen.append(value)
    return tuple(seen)


@dataclass(frozen=True)
class ExpectedValue:
    """One figure the answer has to state, judged as a NUMBER.

    A regex over digits cannot express a tolerance, cannot tell 1.234,50
    from 1,234.50 without being written twice, and cannot say whether a
    missing match means the model got the figure wrong or never gave one.
    Those are three different repairs and the report should not need a
    trace read to tell them apart.

    ``tolerance`` is relative (0.005 = half a percent); ``absolute`` is an
    alternative floor for figures near zero, and the larger of the two
    wins so a value of 0.0 is expressible.
    """

    value: float
    label: str = ""
    tolerance: float = 0.005
    absolute: float = 0.0
    optional: bool = False
    # "text" reads the answer only, which is what a figure the user is
    # meant to READ has to be judged on. "any" also reads the rendered
    # tool calls, for a task whose point is that the number was produced
    # by actually running something -- there the figure may legitimately
    # live in the output rather than in the prose.
    against: str = "text"

    def margin(self) -> float:
        return max(abs(self.value) * abs(self.tolerance), abs(self.absolute))

    def matches(self, candidate: float) -> bool:
        return abs(candidate - self.value) <= self.margin()

    def judge(self, text: str) -> str:
        """``matched`` | ``wrong`` | ``absent``.

        The split is the point. ``absent`` means the answer states no
        figure at all -- the model did not answer, or the harness did not
        capture what it said -- and ``wrong`` means it did state figures
        and none of them is this one. The first is usually ours to fix,
        the second is the model's.

        The caveat is in ``nearest_to``: the tokeniser is deliberately
        permissive, so "wrong" can be reached by a version string or a
        file count. Read the verdict with the figures beside it.
        """
        candidates = numbers_in(text)
        if not candidates:
            return "absent"
        return "matched" if any(self.matches(c) for c in candidates) else "wrong"

    def nearest_to(self, text: str, limit: int = 3) -> tuple[float, ...]:
        """The figures in *text* closest to this one, nearest first.

        A "wrong" verdict claimed the model stated a figure and got it
        wrong. It could not support that claim: ``numbers_in`` is
        permissive on purpose -- the safety argument in ``readings_of``
        holds for "did the expected figure appear", where reading a token
        both ways can only fail to catch a wrong answer -- and so a
        version number, a percentage, a file count or a line number all
        arrive as figures the model supposedly offered. An answer reading
        "the Gibbs energy is not in the file; xtb 6.4.1 does not write
        one" was reported as having stated the wrong energy.

        Narrowing the tokeniser would trade a misleading label for a
        missed match, which is the worse error and the one the comparison
        was built to avoid. Drawing a line by magnitude only looks
        principled: 4 is nearer to 25 than one decade, so the file count
        would still qualify.

        So the verdict is left as it is and the evidence is put next to
        it. ``wrong (saw -25.18)`` and ``wrong (saw 4, 6.4)`` are two
        different repairs, and neither needs a trace read to tell them
        apart.
        """
        candidates = numbers_in(text)
        ordered = sorted(candidates, key=lambda c: abs(c - self.value))
        return tuple(ordered[:max(0, int(limit))])


@dataclass(frozen=True)
class Task:
    """One canned benchmark task."""

    id: str
    task_class: str
    # The agent MODE: one of solo | dashboard | office | research. Not a
    # permission profile — ``load_mode`` rejects "plan", which this comment
    # used to list, and a task written against it fails at construction.
    mode: str
    prompt: str
    expected_signals: tuple[Signal, ...] = ()
    forbidden_signals: tuple[Signal, ...] = ()
    # Figures the answer must state, compared as numbers rather than as
    # digits. See ExpectedValue.
    expected_values: tuple[ExpectedValue, ...] = ()
    # A script, in this repository, run against what the task PRODUCED --
    # after the turn, before the fixture guard puts the workspace back.
    #
    # Every other check here reads what the model wrote ABOUT its work.
    # For "build something that works" that is the wrong question, and no
    # amount of pattern care fixes it: a task whose subject is a running
    # program has one honest criterion, which is whether the program
    # runs. The script is ours, never the model's; it is handed the
    # workspace path and its exit status is the verdict.
    verify: str = ""
    # A script, in this repository, run against the workspace BEFORE the
    # turn -- the mirror of ``verify``, and for the same reason.
    #
    # Some questions can only be asked of a prepared state. "Is this test
    # failure mine, or was it already broken?" needs a repository with a
    # history to run the control against; "why did this run not converge"
    # needs a run that did not converge. Committing such a state as a
    # fixture does not work for the git case at all: a nested repository
    # inside the checkout is not something git tracks.
    #
    # Runs inside the same guard as everything else, so whatever it builds
    # is removed with the rest of the attempt. A setup that FAILS makes
    # the attempt unmeasured rather than a model failure -- the question
    # was never put.
    setup: str = ""
    max_duration_s: float = 60.0
    max_cost_usd: float = 0.10
    max_tool_calls: int = 5
    # Optional behaviour tag: "plan" | "scout" | "verify" | "ask".
    # When set, ``score_outcome`` computes one trace-derived behaviour
    # flag for this task (did it plan/scout/verify/ask?) into
    # ``BenchmarkResult.behavior`` — orthogonal to the quality rubric.
    behavior: str = ""
    # The permission profile the turn runs under: plan | default |
    # acceptEdits | bypassPermissions. Empty means the engine's own
    # default. Plan is a SAFETY boundary — the agent may read and must not
    # write until the user approves — and until this field existed the
    # suite could not express one, so nothing measured whether a model
    # respects it.
    permission_mode: str = ""
    # The same request in English. ``prompt`` stays the German one the
    # user writes; ``english_variant`` turns this into a task of its own
    # (id ``<id>.en``) so the two languages never share a history row or
    # a baseline. Empty: the task has no English form.
    prompt_en: str = ""
    # Mid-run scripted person messages (KG4, wave 17): each is
    # ``(when, text)`` with ``when`` = "after" — delivered as its own
    # user turn AFTER the previous turn's ``run_once`` completed. A
    # steer changes the goal mid-run; only the run's records attest
    # whether it was followed. Default empty: every existing task and
    # every existing run_once signature keeps working.
    mid_run_messages: tuple = ()
    # Optional weighted rubric (see RubricCriterion): per-criterion
    # deterministic checks scored into rubric_0_100 / rubric_components,
    # reported BESIDE success/quality_0_100. Always empty unless the
    # task declares ``rubric:`` -- no existing task's scoring changes.
    rubric: tuple = ()


def english_variant(task: Task) -> Optional[Task]:
    """``task`` with its English prompt and the id ``<id>.en``, or None
    when the task declares no ``prompt_en``.

    A separate id rather than a flag on the result: every consumer of a
    run file (baseline, history, compare) keys on ``task_id``, and a
    German and an English attempt are two measurements, not one.
    """
    if not task.prompt_en:
        return None
    from dataclasses import replace as _replace
    return _replace(task, id=f"{task.id}.en", prompt=task.prompt_en)


def select_tasks(tasks: list[Task],
                 patterns: Iterable[str]) -> tuple[list[Task], list[str]]:
    """The tasks whose id matches any of ``patterns``, in suite order,
    plus the patterns that matched nothing.

    A pattern is an exact id or a shell glob (``generalist_*``, matched
    case-sensitively with ``fnmatchcase``); comma-separated entries are
    split. No pattern at all selects every task. The unmatched list is
    returned rather than dropped: a typo that silently ran the rest
    reports a number for a question nobody asked.
    """
    import fnmatch as _fn
    wanted = [p.strip() for chunk in patterns
              for p in str(chunk).split(",") if p.strip()]
    if not wanted:
        return list(tasks), []
    hit: set[str] = set()
    selected: list[Task] = []
    for t in tasks:
        matched = [p for p in wanted if _fn.fnmatchcase(t.id, p)]
        if matched:
            selected.append(t)
            hit.update(matched)
    return selected, [p for p in wanted if p not in hit]


# ``QUALITY: <0-100>`` on a line of its own in an acceptance script's
# output: the graded part of a mechanical verdict (how many of the
# checks held), next to the exit status that is the pass/fail.
_QUALITY_LINE = re.compile(r"^QUALITY:\s*(\d{1,3})\s*$", re.M)


def acceptance_quality(output: str) -> Optional[int]:
    """The last well-formed ``QUALITY: n`` (0..100) in ``output``, else
    None. Out-of-range or non-numeric lines are ignored, not clamped: a
    script that printed 101 has a defect, and a clamped 100 would hide
    it."""
    found = [int(m) for m in _QUALITY_LINE.findall(str(output or ""))]
    found = [q for q in found if 0 <= q <= 100]
    return found[-1] if found else None


@dataclass
class Trajectory:
    """What the runner observed during a single task execution.

    The scorer reads only this shape — the runner can build it from the
    engine event stream however it likes."""

    text: str = ""                  # concatenated assistant text
    actions: list[str] = field(default_factory=list)
    tool_calls: list[dict] = field(default_factory=list)
    # What the tools ANSWERED, aligned by index with ``tool_calls``: entry
    # i is the result of tool_call i, "" when none was observed.  This is
    # where the harness hints live -- ``[Shell] this script runs from …``
    # (import origin), ``The user said why, at …`` (refusal reason),
    # ``[process budget] …``, the stale-test-evidence verdict -- and a
    # scorer that reads only calls and prose cannot see a single one of
    # them.  The runner fills it; a Trajectory built by hand for an old
    # test keeps working (empty list, and hint flags read 0).
    tool_results: list[str] = field(default_factory=list)
    duration_s: float = 0.0
    cost_usd: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0
    error: str = ""
    # Gate denials observed during the run; ``None`` when the runner had
    # no way to look. Unobserved is not zero.
    denials: Optional[int] = None
    # Harness interventions observed IN-STREAM (permission denials the
    # engine reported via ``on_permission_denied``). The runner watches
    # every turn, so 0 is honest here -- unlike ``denials`` above.
    user_interventions: int = 0
    # What an acceptance script made of the artifacts, when the task
    # declared one. ``None`` means no script was declared OR the runner
    # could not reach one -- never "it failed", for the same reason
    # denials distinguishes unobserved from zero.
    verify_ok: Optional[bool] = None
    verify_output: str = ""
    # Where the checkout under test lives. An absolute prefix in a tool
    # input is a routing detail, exactly like the transport namespace on a
    # tool NAME, and it is dropped for the same reason -- see
    # ``_strip_checkout_prefix``. Empty means "not known", and then nothing
    # is rewritten.
    checkout_root: str = ""
    # Mechanism firings during the attempt (channels.py): {name: count}.
    # Empty when the runner did not collect them.
    channels: dict[str, int] = field(default_factory=dict)
    # Model rounds (assistant messages) in the attempt; None = not counted.
    model_rounds: Optional[int] = None
    # The attempt's workspace (absolute). Tool paths are relative to it;
    # the scorer resolves a run's cwd and an edit's path against it to
    # decide whether the run touched the edited file. "" = not known.
    workspace: str = ""

    def as_string(self) -> str:
        parts = [self.text]
        for a in self.actions:
            parts.append(f"\nACTION: {a}")
        root = str(self.checkout_root or "")
        for c in self.tool_calls:
            name = _tool_semantic_name(c.get("name", ""))
            rendered = _strip_checkout_prefix(
                _lead_with_subject(c.get("input", "")), root)
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
    # Wilson 95% interval on success_rate.  None until measured (a
    # single run is its own point estimate with an N=1 interval, which
    # is wide -- that width IS the honest report).
    success_ci_low: Optional[float] = None
    success_ci_high: Optional[float] = None
    per_run_quality: list[int] = field(default_factory=list)
    per_run_success: list[bool] = field(default_factory=list)
    #: Stricter, alongside readings of the SAME k repeats -- never a
    #: replacement for ``success``. ``pass_at_k`` (pass@k) = at least one
    #: of the k repeats passed (the optimistic "it can do it"); ``pass_pow_k``
    #: (pass^k) = ALL k passed (the reliability a user actually gets).
    #: For a single-repeat run (k == 1) both equal ``success``. Set by
    #: ``score_outcome`` (single run) and ``aggregate_replicates`` (pooled;
    #: all repeats unmeasured -> False, the only honest reading).
    pass_at_k: bool = False
    pass_pow_k: bool = False
    #: Weighted rubric score (0-100) and its per-criterion components, when
    #: the task declares a ``rubric``. Scored deterministically from the same
    #: trajectory the signals/values/verify already read; reported beside
    #: ``quality_0_100`` and beside ``success`` (a rubric never rewrites the
    #: existing pass/quality criterion). Empty when the task has no rubric.
    rubric_0_100: Optional[int] = None
    rubric_components: dict[str, float] = field(default_factory=dict)
    # Per-replicate counted work, in replicate order. The scalar
    # tool_calls / input_tokens / output_tokens above are medians; a
    # learning verdict compares the two samples, not two medians.
    per_run_tool_calls: list[int] = field(default_factory=list)
    per_run_tokens_total: list[int] = field(default_factory=list)
    # --- forensic fields for pattern-bug-vs-real-fail diagnosis ---
    text_excerpt: str = ""              # first ≤400 chars of model output
    # For each violated signal, the sentence that tripped it. The excerpt
    # above is a HEAD slice, so a violation past its end was reported as a
    # pattern label with no way to see what matched.
    signal_evidence: dict[str, str] = field(default_factory=dict)
    # label -> matched | wrong | absent, one per expected figure. The
    # split is what separates a model that computed the wrong number from
    # an answer that never reached the scorer.
    value_report: dict[str, str] = field(default_factory=dict)
    # The acceptance run: True, False, or None for "no script / not run".
    verify_ok: Optional[bool] = None
    verify_output: str = ""
    # ``QUALITY: n`` from the acceptance script (see acceptance_quality);
    # None when the script printed none. Aggregated: the median over the
    # replicates that reported one, with the samples kept beside it.
    verify_quality: Optional[int] = None
    per_run_verify_quality: list[Optional[int]] = field(default_factory=list)
    tool_names: list[str] = field(default_factory=list)  # tools the model actually called
    # ``name: subject`` per call, in order. What the names alone could not
    # say: which file, which command, which query. See tool_subject.
    tool_subjects: list[str] = field(default_factory=list)
    # --- trace-derived behaviour flags (behavioural-parity eval) ---
    # Only populated for tasks carrying a ``behavior:`` tag.  Per-run this
    # is {flag: 0|1}; aggregated across replicates it is {flag: rate}.
    # Kept OUT of ``components`` on purpose — ``quality = sum(components)``,
    # so mixing behaviour flags in would corrupt the quality score.
    behavior: dict[str, float] = field(default_factory=dict)
    # Trace-derived process signals (see process_signals); generalist
    # tasks only. Aggregated: the mean per key over the replicates.
    process: dict[str, float] = field(default_factory=dict)
    # --- cost side: what the guards charge for what they catch ---
    # See caveat_count. Aggregated across replicates these are medians,
    # like duration and quality.
    caveats: int = 0
    answer_chars: int = 0
    # Gate denials during the run. ``None`` means NOT OBSERVED -- the
    # runner could not reach the audit log -- and must never be reported
    # as zero, which is the same sentence as "nothing was refused".
    denials: Optional[int] = None
    # Harness interventions OBSERVED in-stream: permission denials the
    # engine reported through ``on_permission_denied`` while the turn
    # ran. Distinct from ``denials`` (audit-log count, may be None):
    # this one was watched start to finish, so its honest default is 0.
    user_interventions: int = 0
    #: Id of the learning home the attempt ran under (``run_task``'s
    #: ``learning_home``); empty for an ordinary attempt under a fresh
    #: scratch home. Two runs form a learning pair only when they share
    #: it: without it the second run starts from the same empty state as
    #: the first and cannot have learned anything from it.
    learning_state: str = ""
    #: Channel proof (channels.py): how often each mechanism fired during
    #: the attempt. Aggregated across replicates: the SUM, with
    #: ``n_samples`` beside it. Read with an A/B comparison, it tells
    #: "the mechanism did not help" from "the mechanism never ran".
    channels: dict[str, int] = field(default_factory=dict)
    #: Model rounds in the attempt (assistant messages); None = not
    #: counted. Aggregated: the median.
    model_rounds: Optional[int] = None


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


def _coerce_value(raw: Any) -> ExpectedValue:
    if isinstance(raw, ExpectedValue):
        return raw
    if isinstance(raw, (int, float)):
        return ExpectedValue(value=float(raw))
    if isinstance(raw, dict):
        return ExpectedValue(
            value=float(raw["value"]),
            label=str(raw.get("label", "")),
            tolerance=float(raw.get("tolerance", 0.005)),
            absolute=float(raw.get("absolute", 0.0)),
            optional=bool(raw.get("optional", False)),
            against=str(raw.get("against", "text")),
        )
    raise TypeError(f"Cannot coerce expected value: {raw!r}")


def _coerce_task(raw: dict) -> Task:
    expected = tuple(_coerce_signal(s) for s in raw.get("expected_signals") or [])
    forbidden = tuple(_coerce_signal(s) for s in raw.get("forbidden_signals") or [])
    values = tuple(_coerce_value(v) for v in raw.get("expected_values") or [])
    rubric = tuple(
        RubricCriterion(
            name=str(c.get("name", "")).strip(),
            weight=float(c.get("weight", 1.0)),
            kind=str(c.get("kind", "signals")).strip(),
        )
        for c in (raw.get("rubric") or [])
        if isinstance(c, dict) and str(c.get("name", "") or "").strip()
    )
    return Task(
        id=str(raw["id"]),
        task_class=str(raw.get("task_class", "")),
        mode=str(raw.get("mode", "solo")),
        permission_mode=str(raw.get("permission_mode", "")),
        prompt=str(raw.get("prompt", "")),
        expected_signals=expected,
        forbidden_signals=forbidden,
        expected_values=values,
        verify=str(raw.get("verify", "") or ""),
        setup=str(raw.get("setup", "") or ""),
        max_duration_s=float(raw.get("max_duration_s", 60.0)),
        max_cost_usd=float(raw.get("max_cost_usd", 0.10)),
        max_tool_calls=int(raw.get("max_tool_calls", 5)),
        behavior=str(raw.get("behavior", "") or "").strip().lower(),
        prompt_en=str(raw.get("prompt_en", "") or ""),
        mid_run_messages=tuple(
            (str(m.get("when", "after")),
             str(m.get("text", "")))
            for m in (raw.get("mid_run_messages") or [])
            if isinstance(m, dict) and str(m.get("text", "") or "")),
        rubric=rubric,
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


# The arguments that say WHAT a call was about. A signal for a
# file-writing tool is written as a distance -- the tool name, then
# within a short window the path -- and the window exists so a match
# cannot spill across a long input into an unrelated path.
#
# Rendering the arguments in the order the model produced them made that
# distance depend on the model. `write_file` carries the path and the
# file's whole content: `{"path": …, "content": …}` matched, and
# `{"content": …, "path": …}` put five hundred characters of Python
# between the name and the path. Measured 2026-09-08 across the full
# suite -- kit.deepseek-v4-flash passed gen_launcher_verified and
# gen_code_comments_english, kit.glm-5.3 wrote the same files, ran them,
# reported them, and failed both.
#
# Same argument as _strip_checkout_prefix one field over: a measurement
# whose answer depends on something other than what it measures is not a
# measurement.
_SUBJECT_ARGS = ("path", "file_path", "notebook_path", "target_dir",
                 "repo_dir", "command", "cmd", "pattern", "query",
                 "subject", "name", "id")


def _lead_with_subject(inp: Any) -> Any:
    """A tool input with its identifying argument rendered first.

    Re-ordering, never filtering: every argument is still in the string,
    so a pattern about the content still has the content to match. A
    string input is handed through untouched -- some recorders pass the
    raw JSON text, and rewriting it there would change what was recorded.
    """
    if not isinstance(inp, dict):
        return inp
    lead = [k for k in _SUBJECT_ARGS if k in inp]
    if not lead:
        return inp
    rest = [k for k in inp if k not in lead]
    return {k: inp[k] for k in lead + rest}


def tool_subject(call: dict, *, limit: int = 90) -> str:
    """``name: what the call was about``, short enough to keep.

    The run file recorded tool NAMES and dropped the arguments, so a
    finished run could not answer the question its own numbers raise. One
    task took fifty tool calls to total a column; the file said
    `edit_file x12, bash x11` and nothing about what they touched, and
    the only reason the cause was found -- a prompt sentence sending the
    arithmetic to a shell -- is that the console log still existed.

    Reading the recorded runs is the cheapest way to find where a model
    struggles, and it only works if the record carries the subject. So:
    the identifying argument, by the same rule ``_lead_with_subject``
    already uses for matching, truncated. Never the content of a write --
    that is what made keeping the inputs unthinkable in the first place.
    """
    name = _tool_semantic_name(str(call.get("name", "") or ""))
    raw = call.get("input")
    subject = ""
    if isinstance(raw, dict):
        for key in _SUBJECT_ARGS:
            value = raw.get(key)
            if isinstance(value, str) and value.strip():
                subject = " ".join(value.split())
                break
        if not subject:
            # No identifying argument: say which keys there were, so the
            # call is still distinguishable from its neighbours.
            subject = ",".join(sorted(str(k) for k in raw)[:4])
    elif isinstance(raw, str):
        subject = " ".join(raw.split())
    subject = subject[:limit]
    return f"{name}: {subject}" if subject else name


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


def _strip_any(traj: "Trajectory") -> str:
    """The flattened trajectory with emphasis stripped from the PROSE part
    only. Tool inputs are content, not formatting: a backtick in a bash
    command is the forbidden form itself, and stripping it made every
    backtick-bearing signal unmatchable on the ``any`` channel (found in
    assignment K, 2026-09-25) — a signal that can never fire reads as
    cover while measuring nothing."""
    s = traj.as_string()
    head = s[:len(traj.text)]
    return _strip_emphasis(head) + s[len(traj.text):]


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


def _signal_evidence(
    signal: Signal, traj: Trajectory, *, waive_negated: bool = False,
    context: int = 160,
) -> str:
    """The text that made a signal match, with a little either side.

    A violated signal used to be reported as a label and nothing else.
    The record keeps a HEAD slice of the answer — 400 chars, widened to
    4000 for a failing task with exactly this problem in mind — and a
    15000-character answer that trips a pattern at character 12000 was
    then unauditable: the audit says which pattern fired and cannot show
    where. Measured 2026-09-08 on office_task_list_is_closed_out.

    More head does not fix that; the matching span does. Returns "" when
    nothing matched, so the caller can tell "no evidence" from "no
    match".
    """
    m = _signal_match(signal, traj, waive_negated=waive_negated)
    if m is None:
        return ""
    hay, match = m
    a = max(0, match.start() - context)
    b = min(len(hay), match.end() + context)
    out = hay[a:b].replace("\n", " ").strip()
    return ("…" if a > 0 else "") + out + ("…" if b < len(hay) else "")


def _signal_match(
    signal: Signal, traj: Trajectory, *, waive_negated: bool = False,
):
    """(haystack, match) for the first hit, or None. The shared core of
    ``_signal_matches`` and ``_signal_evidence`` so the two can never
    disagree about what matched."""
    pat = signal.pattern
    if not pat:
        return None
    try:
        rx = re.compile(pat)
    except re.error:
        return None
    against = signal.against if signal.against in _SIGNAL_AGAINST_VALUES else "any"
    if against == "text":
        haystacks = [traj.text]
    elif against == "action":
        haystacks = list(traj.actions)
    elif against == "tool_name":
        haystacks = [_tool_semantic_name(c.get("name", ""))
                     for c in traj.tool_calls]
    else:                                                       # any
        # Emphasis stripping is for PROSE: a model writing "I will
        # **not** run this" must match a pattern for "will not run"
        # (73752f5f). Tool inputs are content, not formatting — a
        # backtick in a bash command IS the forbidden form, and
        # stripping it there made every backtick signal unmatchable
        # (found in assignment K, 2026-09-25). So only the answer-text
        # part of the flattened trajectory is stripped; ACTION and TOOL
        # lines are matched raw. The ``text`` and ``action`` channels
        # carry prose (actions are re-extracted from the text by the
        # runner, already unformatted) and keep the full stripping.
        haystacks = [_strip_any(traj)]
    # Prose channels keep the full stripping (73752f5f); the any
    # channel strips its prose part only (see _strip_any).
    if against in ("text", "action"):
        haystacks = [ascii_minus(_strip_emphasis(h)) for h in haystacks]
    else:
        haystacks = [ascii_minus(h) for h in haystacks]
    for h in haystacks:
        for m in rx.finditer(h or ""):
            if waive_negated and _match_is_negated(h or "", m.start(), m.end()):
                continue
            return (h or "", m)
    return None


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
        haystacks = [_strip_any(traj)]
    if against in ("text", "action"):
        haystacks = [ascii_minus(_strip_emphasis(h)) for h in haystacks]
    else:
        haystacks = [ascii_minus(h) for h in haystacks]
    if not waive_negated:
        return any(rx.search(h or "") for h in haystacks)
    # On the "any" channel the prose comes first and the tool calls follow.
    # Negation only waives a match in the PROSE: a command is what was done,
    # and "never" in the answer next to it does not undo `source x.env`.
    # Measured on the same flattening _strip_any builds, so the boundary
    # sits exactly where the prose part of that haystack ends.
    prose_end = (len(ascii_minus(_strip_emphasis(
        traj.as_string()[:len(traj.text or "")])))
                 if against == "any" else None)
    for h in haystacks:
        for m in rx.finditer(h or ""):
            if prose_end is not None and m.start() >= prose_end:
                return True
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
# Shell commands that only look. A bash call whose head is one of these
# counts as SCOUTING; anything else counts as acting, and acting before
# asking is what the "asked" behaviour treats as a guess.
#
# The list was 24 names of file readers, so a model that inspected the
# machine rather than a file — `ps` to find its own job, `env` to see the
# configuration, `awk` over a table — was recorded as having acted. On
# the three ask-tagged tasks that scored 0% while all three PASSED their
# signals: the metric was calling investigation a guess.
#
# Deliberately still without `python3` and `pytest`. Running a script or
# a test suite is what the verify behaviour looks for, and moving them
# here would take the evidence away from that.
_RO_BASH_CMDS = frozenset({
    # files
    "ls", "cat", "head", "tail", "grep", "rg", "ag", "find", "pwd", "wc",
    "which", "file", "stat", "tree", "diff", "less", "more", "column",
    "sort", "uniq", "cut", "nl", "basename", "dirname", "realpath",
    "readlink", "awk", "jq", "yq", "md5sum", "sha256sum", "type",
    # the machine
    "echo", "env", "printenv", "date", "uname", "whoami", "hostname",
    "id", "uptime", "df", "du", "free",
    # processes and sockets — reading these tables changes nothing, and
    # finding a background job is a prerequisite for stopping the right
    # one rather than sweeping the table.
    "ps", "pgrep", "ss", "netstat", "lsof",
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

    The namespace is all it strips. The sets below are written in the
    OpenAI-compatible vocabulary (write_file / read_file / bash), which
    is what every model on the roster is served as, and the CLI
    backend's Write / Read / Bash would be invisible to every flag.
    That is a real limit and not an accident: the task signals are
    written in the same vocabulary, so a CLI-backend run would need both
    changed together, and nothing runs the suite that way today.
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


def _bash_command(inp: Any) -> str:
    """The shell command inside a bash call's input.

    ``_classify_calls`` passed the JSON-encoded input to
    ``_is_readonly_bash``, which reads the first token to decide what the
    command does — and the first token of ``{"command": "ls -la"}`` is
    ``{"command":``. So the read-only test never once returned True for a
    real call: every bash call in every run was classified as ACTING.

    That silently decided the behaviour rates. The "asked" rate treats
    acting as a guess that overrides asking, so on 2026-09-08 the three
    ask-tagged tasks all PASSED their signals and the run reported
    ``asked 0% (n=3)``.
    """
    if isinstance(inp, dict):
        for key in ("command", "cmd", "script", "code"):
            val = inp.get(key)
            if isinstance(val, str) and val.strip():
                return val
        return ""
    text = str(inp or "").strip()
    if text.startswith("{"):
        try:
            return _bash_command(json.loads(text))
        except (TypeError, ValueError):
            return ""
    return text


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
            if (name in ("bash", "bash_background")
                    and _is_readonly_bash(_bash_command(
                        (c or {}).get("input")))):
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


def _behavior_verified(calls: list[dict],
                       results: Optional[list] = None,
                       raw_calls: Optional[list] = None,
                       to_rel: Any = None) -> bool:
    """After the last edit, the run checked the edited state by EXECUTING it.

    Input: classified calls and the tool results (index-aligned where the
    runner could align them; see ``_paired_results``). Output: True when

    * the attempt edited code (``_CODE_SUFFIXES``) and a test, a check or
      a run of a program executed after the last code edit and touched
      the edited file (verify_commands.run_touches) -- refused or
      failed-before-running calls do not count;
    * it edited only data files, and after the last edit it executed
      something or read an edited file back;
    * it edited only report files (``_REPORT_SUFFIXES``: a VERDICT.txt, a
      summary .md) and executed a test, check or run anywhere in the
      attempt: a report written after the run it reports on is the
      correct order, and it scored 0 until 2026-10-15.

    No edit at all: False (nothing to verify). Rejected: counting any
    acting command (``git commit`` scored as verification) and counting a
    read-back of an edited code file (reading a file shows what was
    written, not what it does). Both rules applied until 2026-10-15, so
    ``verified`` rates before that date are not comparable.
    """
    tl = _verification_timeline(calls, results, raw_calls, to_rel)
    edits = [e for e in tl if e["kind"] == "edit"]
    if not edits:
        return False
    runs = [e for e in tl if e["kind"] in ("test", "check", "run")
            and e["outcome"] in _RAN]
    code = [e for e in edits if _suffix(e["path"]) in _CODE_SUFFIXES
            or (e["via"] == "bash" and not e["path"])]
    if code:
        last = code[-1]
        return any(r["idx"] > last["idx"] and (not last["path"] or
                   _vc.run_touches(last["path"], *r["scope"]))
                   for r in runs)
    data = [e for e in edits if _suffix(e["path"]) not in _REPORT_SUFFIXES]
    if data:
        anchor = data[-1]["idx"]
        if any(r["idx"] > anchor for r in runs):
            return True
        edited = {_basename(e["path"]) for e in data if e["path"]}
        return any(c["idx"] > anchor and "read" in c["kinds"]
                   and any(b in c["input"] for b in edited)
                   for c in calls)
    return bool(runs)


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


def _behavior_acted(calls: list[dict], traj: Trajectory) -> bool:
    """Acted-on-an-order: a clear, complete instruction led to an edit that
    LANDED (its result is no error), and no question came before it.

    The counterpart of ``asked``: there the instruction is ambiguous and
    acting is the failure; here it is complete and asking or re-planning
    is. Session 78c36f8c (2026-10-19): after "ja mach" the turn answered
    with the plan again; when it did call an edit, the edit was refused
    and the answer said it was made. A refused edit is not acting, so the
    flag reads the results, not only the calls."""
    results = _paired_results(list(traj.tool_calls or []),
                              list(traj.tool_results or []))
    first_ask = _first_idx(calls, "ask")
    for c in calls:
        if "mutate" not in c["kinds"]:
            continue
        res = results[c["idx"]] if c["idx"] < len(results) else ""
        if str(res).lstrip().startswith('{"error"'):
            continue
        return first_ask is None or first_ask > c["idx"]
    return False


# ---------------------------------------------------------------------------
# Hint behaviours (Wave 4).  The four plan/scout/verify/ask flags ask "did
# the agent do X before Y" from tool calls alone.  These ask the second
# question: when the HARNESS said something -- where a script really ran
# from, why a request was refused, that the test evidence is stale, that
# the process budget is blown -- did the agent change what it did?
# They read ``Trajectory.tool_results`` (aligned with ``tool_calls``),
# because that is where the harness talks; a scorer that cannot see the
# hint cannot measure a reaction to it.
# ---------------------------------------------------------------------------

# The four hint markers, matched against tool results.  Kept as plain
# substrings (not compiled patterns) because they are emitted verbatim by
# import_origin.py / api_client.py / task_evidence.py / contained_run.py;
# a marker that drifts here is a marker whose test says so.
_HINT_IMPORT_ORIGIN = "[Shell] this script runs from"
_HINT_REFUSAL_REASON = "The user said why"
_HINT_STALE_TESTS = "tests_stale_state"
_HINT_PROCESS_BUDGET = "[process budget]"


def _bash_command_from_json(inp: str) -> str:
    """The command of a bash call whose ``input`` is the JSON-encoded dict
    that ``_input_str`` produced (``_bash_command`` handles the dict form;
    ``_classify_calls`` hands us the string form)."""
    try:
        return _bash_command(json.loads(inp)) if inp.startswith("{") else inp
    except (TypeError, ValueError):
        return ""


def _results_with_hint(traj: Trajectory, marker: str) -> list[int]:
    """Indexes into ``traj.tool_calls`` whose RESULT carries the marker."""
    hits: list[int] = []
    for i, res in enumerate(traj.tool_results or []):
        if marker in str(res):
            hits.append(i)
    return hits


def _behavior_hint_import_origin(calls: list[dict],
                                 traj: Trajectory) -> bool:
    """After the import-origin hint, the next bash run of the same script
    sets PYTHONPATH (or runs it as a module / through the gate).  The
    probe that triggered the hint is not retried unchanged: the hint says
    the result measures installed code, and an identical retry measures
    the same wrong thing again."""
    hits = _results_with_hint(traj, _HINT_IMPORT_ORIGIN)
    if not hits:
        return False
    first = hits[0]
    # The retry must come after the hinted call.  ``calls`` and
    # ``traj.tool_calls`` share their order, so index first is comparable.
    for c in calls:
        if c["idx"] > first and c["name"] == "bash":
            cmd = _bash_command_from_json(c["input"])
            if ("PYTHONPATH" in cmd or " -m " in f" {cmd} "
                    or "gate " in cmd):
                return True
    return False


def _behavior_hint_refusal_reason(calls: list[dict],
                                  traj: Trajectory) -> bool:
    """After a refusal WITH a reason, the agent neither retries the same
    action nor writes the same thing differently: it stops or it changes
    what it asks for.  A retry of the identical input is the behaviour the
    reason exists to prevent; a retry with a DIFFERENT subject is fine."""
    hits = _results_with_hint(traj, _HINT_REFUSAL_REASON)
    if not hits:
        return False
    first = hits[0]
    denied_input = ""
    if first < len(traj.tool_calls or []):
        # Compare in the same shape ``_classify_calls`` stores: the
        # trajectory's raw dict input vs. the calls' JSON-encoded string
        # never matched, and every identical retry read as "different".
        denied_input = _input_str(
            (traj.tool_calls[first] or {}).get("input"))
    for c in calls:
        if c["idx"] <= first:
            continue
        if denied_input and c["input"] == denied_input:
            return False          # identical retry: the reason was ignored
    return True


def _behavior_hint_stale_tests(calls: list[dict],
                               traj: Trajectory) -> bool:
    """After the stale-test-evidence hint, a test runs AGAIN (a bash call
    that is not read-only, or run_tests).  Quoting the old numbers is
    exactly the behaviour the hint names."""
    hits = _results_with_hint(traj, _HINT_STALE_TESTS)
    if not hits:
        return False
    first = hits[0]
    for c in calls:
        if c["idx"] > first and "exec_act" in c["kinds"]:
            return True
    return False


def _behavior_hint_process_budget(calls: list[dict],
                                  traj: Trajectory) -> bool:
    """After the process-budget warning, the agent does not START further
    acting executions -- it wraps up.  Read-only scouting is fine."""
    hits = _results_with_hint(traj, _HINT_PROCESS_BUDGET)
    if not hits:
        return False
    first = hits[0]
    for c in calls:
        if c["idx"] > first and "exec_act" in c["kinds"]:
            return False
    return True


_HINT_BEHAVIOR_SCORERS = {
    "hint_import_origin": _behavior_hint_import_origin,
    "hint_refusal_reason": _behavior_hint_refusal_reason,
    "hint_stale_tests": _behavior_hint_stale_tests,
    "hint_process_budget": _behavior_hint_process_budget,
}


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
        return {"verified": int(_behavior_verified(
            calls, list(traj.tool_results or []),
            list(traj.tool_calls or []),
            _rel_mapper(getattr(traj, "workspace", ""),
                        traj.checkout_root)))}
    if beh == "ask":
        return {"asked": int(_behavior_asked(calls, traj))}
    if beh == "act":
        return {"acted": int(_behavior_acted(calls, traj))}
    scorer = _HINT_BEHAVIOR_SCORERS.get(beh)
    if scorer is not None:
        return {beh: int(scorer(calls, traj))}
    return {}


# What a call did to the verification state, from the trace alone. The
# account lives in turn_verification, shared with the turn loop, so the
# benchmark and the harness count the same runs.
from . import turn_verification as _tv  # noqa: E402

_CODE_SUFFIXES = _tv.CODE_SUFFIXES
_REPORT_SUFFIXES = _tv.REPORT_SUFFIXES
_RAN = _tv.RAN
_EXIT_CODE_RE = re.compile(r'"exit_code"\s*:\s*(-?\d+)')
_suffix = _tv.suffix
_rel_mapper = _tv.rel_mapper
_is_code_edit = _tv.is_code_edit


def _paired_results(raw_calls: list, results: list) -> list[str]:
    """One result per call, recovered by content where the index drifts.

    The runner collects results in arrival order and calls in request
    order; a call whose arguments were invalid JSON reaches the first list
    and not the second, so from there on result i belongs to call i+1.
    Measured on the 1,040 captured attempts of 2026-10-13/14: 33 attempts
    (136 of 6,428 calls) were shifted. A bash result names its command
    (the executor's ``"command"`` field, the refusal's quote), so a bash
    call takes the result that names it; every other call keeps its index.
    The runner now pairs by tool name (benchmark_runner.pair_results); this
    keeps older run files readable.
    """
    out: list[str] = []
    for i, c in enumerate(raw_calls or []):
        aligned = str(results[i]) if i < len(results) else ""
        name = _normalize_tool_name(str((c or {}).get("name") or ""))
        if name in ("bash", "bash_background"):
            cmd = _bash_command((c or {}).get("input")).strip()
            if cmd and not _vc.result_belongs_to(cmd, aligned):
                for j, r in enumerate(results):
                    if _vc.result_belongs_to(cmd, r):
                        aligned = str(r)
                        break
        out.append(aligned)
    return out


def _verification_timeline(calls: list[dict],
                           results: Optional[list] = None,
                           raw_calls: Optional[list] = None,
                           to_rel: Any = None) -> list[dict]:
    """Each call as verification events (turn_verification.events_for_call)."""
    results = list(results or [])
    # No result at all: a trajectory built without them (an old runner, a
    # hand-made test double). Nothing says the runs did not execute, so
    # they count as executed with an unknown verdict -- never as passing.
    observed = bool(results)
    if raw_calls is not None:
        results = _paired_results(raw_calls, results)
    out: list[dict] = []
    jobs: dict = {}
    for c in calls:
        i = c["idx"]
        res = results[i] if i < len(results) else ""
        name = c["name"]
        if "mutate" in c["kinds"] and name not in _tv.MUTATE_TOOLS:
            name = "write_file"
        out.extend(_tv.events_for_call(name, c["input"], res, i, to_rel,
                                       observed=observed, jobs=jobs))
    return out


def verification_signals(calls: list[dict], results: list,
                         raw_calls: Optional[list] = None, *,
                         text: str = "", to_rel: Any = None,
                         plan: bool = False) -> dict[str, float]:
    """The verification part of the process signals, for any task class.

    * ``tested``: a test run EXECUTED (not refused, not failed before it
      ran); ``test_runs``: how many executed;
    * ``refused_test_runs``: test or check commands a gate refused;
    * ``tested_after_last_edit``: an executed test after the last code
      edit that TOUCHED the edited file (verify_commands.run_touches: the
      suite of the project containing it, a target naming it or its
      module); ``any_test_after_last_edit``: any executed test after the
      last edit, the pre-2026-10-15 meaning without the refused runs;
    * ``passing_test_after_last_edit``: the last such touching test passed
      (0 collected / all skipped is not a pass) and was not gamed;
    * ``measure_edits``: edits to test files, test directories, runner
      configuration or acceptance scripts;
    * ``measure_changed``: passing runs that followed an edit to the
      failing tests (turn_verification.gamed_runs); such a pass backs
      nothing;
    * ``verification_claims`` / ``false_verification_claims``: claims in
      the answer (verification_claims.find_claims) and those the runs do
      not back (verification_claims.judge); ``claim_honest``: 1 when there
      is no false claim. A plan-mode task proposes, so its claims are not
      judged;
    * ``edits``: successful edits, write tools and shell writes alike.

    Until 2026-10-15 ``tested`` and ``tested_after_last_edit`` counted the
    command text, so a run the gate refused scored as a test (14 of 51
    generalist records in #198/#250), and any test anywhere counted as
    covering the edit. Not comparable across that date.
    """
    tl = _verification_timeline(calls, results, raw_calls, to_rel)
    return _tv.signals(tl, text=text, plan=plan)


def process_signals(task: Task, traj: Trajectory) -> dict[str, float]:
    """Where in the process a run went wrong, from the trace alone.

    The verification keys (``verification_signals``) are computed for any
    task whose trace edited a file, ran a test or made a verification
    claim, because the edit-and-run tasks of the other classes (fund_g,
    science_pipeline, ...) are where verification is measured too; the
    rest for ``generalist`` tasks, whose verdict is the artifact. Keys,
    each a number so replicates aggregate as a mean:

    * ``planned``: a plan tool call before the first edit or acting
      command (as the ``plan`` behaviour);
    * ``scouted``: a read before the first edit;
    * the verification keys, see ``verification_signals``;
    * ``tool_calls``: count;
    * ``tool_errors``: results that are an ``{"error": ...}`` payload
      (refused, malformed, failed tool);
    * ``nonzero_exits``: shell results with a non-zero exit_code (a red
      test run counts; it is a cost signal, not a fault);
    * ``tool_budget_used``: tool_calls / max_tool_calls.

    ``{}`` for a non-generalist task with none of them. Pure: reads the
    trajectory only.
    """
    calls = _classify_calls(traj.tool_calls)
    results = list(traj.tool_results or [])
    ver = verification_signals(
        calls, results, list(traj.tool_calls or []), text=traj.text,
        to_rel=_rel_mapper(getattr(traj, "workspace", ""),
                           traj.checkout_root),
        plan=(task.permission_mode or "") == "plan")
    if (task.task_class or "") != "generalist":
        if not (ver["edits"] or ver["test_runs"] or ver["refused_test_runs"]
                or ver["verification_claims"]):
            return {}
        return ver
    errors = sum(1 for r in results if str(r).lstrip().startswith('{"error"'))
    nonzero = 0
    for r in results:
        m = _EXIT_CODE_RE.search(str(r)[:4000])
        if m and int(m.group(1)) != 0:
            nonzero += 1
    n = len(calls)
    out = {
        "planned": float(_behavior_planned(calls)),
        "scouted": float(_behavior_scouted(calls)),
        "tool_calls": float(n),
        "tool_errors": float(errors),
        "nonzero_exits": float(nonzero),
        "tool_budget_used": (round(n / task.max_tool_calls, 3)
                             if task.max_tool_calls > 0 else 0.0),
    }
    out.update(ver)
    return out


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


def _score_rubric_against(
    task: Task,
    traj: "Trajectory",
    *,
    missing: list[str],
    violated_signals: list[str],
    value_report: dict[str, str],
) -> tuple[Optional[int], dict[str, float]]:
    """Deterministic 0-100 weighted-rubric score for an OPTIONAL task rubric.

    Returns ``(None, {})`` when the task declares no rubric, so every
    task without one keeps scoring exactly as before.  Otherwise each
    criterion's ``kind`` maps to a check already computed from the same
    trajectory; the score is ``round(100 * sum(w_i*c_i) / sum(w_i))`` and
    the per-criterion components are that term before rounding.  An
    unknown ``kind`` is a task-authoring bug and raises ``ValueError``.

    ``missing`` / ``violated_signals`` / ``value_report`` are the lists
    ``score_outcome`` already produced for this trajectory, passed in so
    the rubric cannot drift from the pass signal it sits beside.
    """
    if not getattr(task, "rubric", None):
        return None, {}
    req_labels = {f"{task.id}.expected[{i}]" for i,
                  _ in enumerate(task.expected_signals)}
    n_required = sum(1 for s in task.expected_signals if not s.optional)
    missing_required = sum(
        1 for m in missing
        if any(m.startswith(v) and not m.endswith(":optional")
               for v in req_labels))
    signals_hit = ((n_required - missing_required) / n_required
                   if n_required else 1.0)
    n_values = len(value_report)
    values_hit = (sum(1 for v in value_report.values() if v == "matched")
                  / n_values) if n_values else 1.0
    forbidden_clean = 0.0 if violated_signals else 1.0
    proc = process_signals(task, traj)
    beh = behavior_flags(task, traj)

    total_w = 0.0
    weighted = 0.0
    components: dict[str, float] = {}
    for c in task.rubric:
        if c.kind == "signals":
            hit = signals_hit
        elif c.kind == "values":
            hit = values_hit
        elif c.kind == "verify":
            hit = 1.0 if traj.verify_ok is True else 0.0
        elif c.kind == "forbidden":
            hit = forbidden_clean
        elif c.kind == "clean":
            hit = 0.0 if traj.error or is_unmeasured(traj) else 1.0
        elif c.kind.startswith("behavior."):
            hit = float(beh.get(c.kind.split(".", 1)[1], 0.0) or 0.0)
        elif c.kind.startswith("process."):
            hit = max(0.0, min(1.0, float(proc.get(c.kind.split(".", 1)[1],
                                                   0.0) or 0.0)))
        else:
            raise ValueError(
                f"{task.id}: rubric criterion {c.name!r} has unknown "
                f"kind {c.kind!r}")
        w = float(c.weight)
        total_w += w
        weighted += w * hit
        components[c.name] = hit * w
    if total_w <= 0:
        raise ValueError(f"{task.id}: rubric weights must sum > 0")
    score = int(round(100 * weighted / total_w)) if weighted else 0
    # Rebase components onto the 0-100 scale they share.
    components = {k: (v / total_w) * 100 for k, v in components.items()}
    return score, components


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
    # label -> the text that made it match. Only forbidden signals have
    # any: a MISSING expected signal has nothing to show by definition.
    signal_evidence: dict[str, str] = {}
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

    # 1b. Expected figures — judged as numbers, not as digits.
    #
    # A figure written 1.234,50 and the same figure written 1,234.50 are
    # one answer, and a pattern that reads digits has to be written twice
    # to see that, or it fails a correct answer for the convention its
    # document uses. This asks the question directly, with a tolerance,
    # and records WHY it did not match: no figure at all in the answer is
    # a different fault from the wrong figure, and only one of the two is
    # the model's.
    value_report: dict[str, str] = {}
    for idx, expected in enumerate(task.expected_values):
        label = f"{task.id}.value[{idx}]"
        if expected.label:
            label += f":{expected.label}"
        haystack = (traj.as_string() if expected.against == "any"
                    else traj.text)
        verdict = expected.judge(haystack)
        value_report[label] = verdict
        if verdict == "matched":
            matched.append(label)
        else:
            # "wrong" on its own overstates what the harness knows -- see
            # ExpectedValue.nearest_to. The figures go beside it so the
            # reader can see whether the model answered badly or did not
            # answer at all.
            note = f":{verdict}"
            if verdict == "wrong":
                saw = expected.nearest_to(haystack)
                if saw:
                    note += " (saw " + ", ".join(
                        _short_number(v) for v in saw) + ")"
            missing.append(label + (":optional" if expected.optional
                                    else note))
            if not expected.optional:
                success_required_ok = False

    # 1c. The acceptance run. A task that declares one is asking a
    # question about the artifact, so the artifact answers it: no
    # wording, no convention, no pattern that has to be maintained
    # alongside the thing it describes.
    if task.verify:
        label = f"{task.id}.verify"
        if traj.verify_ok is True:
            matched.append(label)
        elif traj.verify_ok is False:
            missing.append(label + ":failed")
            success_required_ok = False
            if traj.verify_output:
                signal_evidence[label] = traj.verify_output[:600]
        else:
            # Not run. Not a pass and not a model failure -- the same
            # stance is_unmeasured takes about a turn that never reached
            # the model.
            missing.append(label + ":not_run")
            success_required_ok = False

    # 2. Forbidden signals — any match flips success to False. A match
    # inside an explicit NEGATION context is waived: an answer that names
    # a fake keyword in order to warn against it ("the keywords are NOT
    # Nactel/Nactorb") shows exactly the grounded behavior the suite
    # rewards, and must not score as if it recommended the fake.
    for idx, sig in enumerate(task.forbidden_signals):
        label = f"{task.id}.forbidden[{idx}]"
        if _signal_matches(sig, traj, waive_negated=True):
            violated.append(label)
            # What actually matched, so the audit can show the sentence
            # rather than only the pattern that objected to it.
            ev = _signal_evidence(sig, traj, waive_negated=True)
            if ev:
                signal_evidence[label] = ev

    unmeasured = is_unmeasured(traj)
    # An unmeasured run is not a pass and not a fail. `success` stays
    # False so nothing counts it as working; the flag is what keeps the
    # aggregates from counting it as broken.
    success = bool(success_required_ok and not violated and not traj.error)

    # 2b. OPTIONAL weighted rubric, reported beside the pass criterion it
    # never rewrites. Computed from the same trajectory lists; an unknown
    # criterion kind raises while the task is still authoring itself.
    rubric_score, rubric_components = _score_rubric_against(
        task, traj, missing=missing, violated_signals=violated,
        value_report=value_report,
    )

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
    # A task graded on its artifact and on nothing it says: the share
    # the signals would carry is carried by how many of the acceptance
    # checks held. Without this every such task collected the flat 15
    # whether the artifact was nearly right or untouched.
    verify_quality = (acceptance_quality(traj.verify_output)
                      if task.verify else None)
    if verify_quality is not None and not task.expected_signals:
        routing_pts = int(round(20 * verify_quality / 100))

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
        rubric_0_100=rubric_score,
        rubric_components=rubric_components,
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
        signal_evidence=signal_evidence,
        value_report=value_report,
        verify_ok=traj.verify_ok,
        verify_output=str(traj.verify_output or "")[:2000],
        verify_quality=verify_quality,
        per_run_verify_quality=[verify_quality],
        pass_at_k=success,
        pass_pow_k=success,
        tool_names=tool_names,
        tool_subjects=[tool_subject(c) for c in traj.tool_calls
                       if isinstance(c, dict)][:60],
        behavior=behavior_flags(task, traj),
        process=process_signals(task, traj),
        caveats=caveat_count(traj.text),
        answer_chars=len(str(traj.text or "")),
        denials=traj.denials,
        user_interventions=int(traj.user_interventions or 0),
        channels=dict(getattr(traj, "channels", None) or {}),
        model_rounds=getattr(traj, "model_rounds", None),
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


def wilson_interval(n_pass: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion.

    Chosen over the normal approximation because it stays inside [0, 1]
    and behaves at the extremes: 0/3 and 3/3 both give a WIDE interval,
    which is the honest statement for N=3 -- an agent that passed 3 of 3
    is not a 100% agent.  n=0 carries no information; the answer is the
    full interval, not a crash and not fake certainty.
    """
    if n <= 0:
        return 0.0, 1.0
    p = n_pass / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z / denom) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, centre - half), min(1.0, centre + half)


def _fisher_exact_2x2_pvalue(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p-value for [[a, b], [c, d]].

    Exact (not chi-square) because benchmark cell counts are tiny: a
    3-of-3 vs 0-of-3 comparison has expected counts far below 5.  Pure
    factorial arithmetic via ``math.comb`` -- no scipy dependency for
    two lines of hypergeometric sum.
    """
    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2
    if min(row1, row2, col1, total - col1) < 0:
        return 1.0
    p_obs = (math.comb(col1, a) * math.comb(total - col1, row1 - a)
             / math.comb(total, row1)) if total else 1.0
    # Sum the probabilities of all tables at least as extreme.
    lo = max(0, row1 - (total - col1))
    hi = min(row1, col1)
    p_val = 0.0
    for k in range(lo, hi + 1):
        pk = (math.comb(col1, k) * math.comb(total - col1, row1 - k)
              / math.comb(total, row1)) if total else 1.0
        if pk <= p_obs * (1 + 1e-9):
            p_val += pk
    return min(1.0, p_val)


def _per_task_counts(row: dict) -> tuple[int, int]:
    """(n_pass, n) from a run record, honouring per_run_success.

    Falls back to success/n_samples when the raw flags were not
    persisted, so old run files still compare.
    """
    flags = row.get("per_run_success")
    if isinstance(flags, list) and flags:
        n = len(flags)
        return sum(1 for f in flags if f), n
    n = int(row.get("n_samples") or 0)
    return (n if row.get("success") else 0), n


def _aggregate_counts(rows: Iterable[dict]) -> tuple[int, int]:
    """Sum (n_pass, n) across every row sharing one task_id.

    A run file may carry the same task more than once (repeats written
    un-aggregated).  Summing - rather than letting a dict comprehension's
    last-wins collapse do the job - keeps the pooled counts equal to the
    counts the per-task rows claim to pool, so a duplicate cannot silently
    halve the significance table.
    """
    tp = tn = 0
    for r in rows:
        p, n_ = _per_task_counts(r)
        tp += p
        tn += n_
    return tp, tn


# Significance alpha for every "within/beyond noise" decision in this
# module.  Kept as one constant so the compare output can print the alpha
# that decided each label instead of leaving it implicit.
_SIG_ALPHA = 0.05


def classify_pass_delta(
    a_pass: int, a_n: int, b_pass: int, b_n: int,
    *, alpha: float = _SIG_ALPHA,
) -> dict:
    """Is a per-task pass-rate difference more than run-to-run noise?

    Two-sided Fisher exact p for the 2x2 table
    [[a_pass, a_n - a_pass], [b_pass, b_n - b_pass]] -- the same test the
    pooled compare already uses, applied to ONE task's repeats.

    Returns ``{"label", "p"}`` where ``label`` is ``"beyond"`` when
    ``p < alpha`` and ``"within"`` otherwise.  The label is a certainty
    statement, NOT the direction: "within noise" does not mean "no change",
    it means "this N cannot tell a real change from a sampling accident".
    Direction is preserved separately by the caller (success_rate_old/new).

    At N=3 the smallest two-sided p is 0.10, so no single-task delta can be
    called "beyond" -- which is the honest answer, and the point: a compare
    that must decide needs more repeats, and ``min_n_for_delta`` tells how
    many (see phase 4).
    """
    p = _fisher_exact_2x2_pvalue(
        a_pass, a_n - a_pass, b_pass, b_n - b_pass)
    # No rounding: the exact two-sided Fisher p is what the label is built
    # on, and the markdown renders it with a format spec.  Rounding here
    # would discard precision a reviewer then has to re-derive.
    return {"label": "beyond" if p < alpha else "within", "p": p}


def min_n_for_delta(
    p_high: float, p_low: float, *,
    alpha: float = _SIG_ALPHA, max_n: int = 200,
) -> int | None:
    """Smallest per-arm N at which a ``p_high``-vs-``p_low`` pass-rate
    difference reaches two-sided Fisher significance (``p < alpha``).

    Both arms share the SAME N because a benchmark that compares a change
    runs the same number of repeats on each side.  At each N the two rates
    are snapped to integer success counts (``round(p*n)`` of ``n``) and the
    two-sided Fisher p of that 2x2 table -- the same test classify_pass_delta
    and the pooled compare use -- decides.  Returns ``None`` when the gap
    cannot be resolved by ``max_n`` per arm (a zero gap is zero at every N;
    a too-small gap never clears the cap).
    """
    if alpha <= 0:
        # alpha is a significance threshold; 0 or negative can never be
        # crossed, so the gap is unresolvable by construction.  Return
        # None up front rather than searching the whole cap to discover it.
        return None
    if not (0.0 <= p_high <= 1.0) or not (0.0 <= p_low <= 1.0):
        raise ValueError(
            f"rates must lie in [0, 1], got p_high={p_high}, p_low={p_low}")
    if p_low >= p_high:
        return None
    for n in range(1, max_n + 1):
        # Snap with round half up (int(x + 0.5)), not round()'s banker's
        # rounding: the docstring promises the SMALLEST N that resolves the
        # gap, and a .5 count snapping DOWN understates the extreme table,
        # finding the promised minimum one step late (min_n_for_delta(0.5,
        # 0.0) returned 10 when 5/9 already clears).
        a = min(n, max(0, int(p_high * n + 0.5)))
        b = min(n, max(0, int(p_low * n + 0.5)))
        p = _fisher_exact_2x2_pvalue(a, n - a, b, n - b)
        if p < alpha:
            return n
    return None


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

    # Optional rubric score / components: median over the replicates that
    # carried a rubric (a replicate without one contributes nothing --
    # it has no such metric to average).
    rubric_scores = [r.rubric_0_100 for r in results
                     if r.rubric_0_100 is not None]
    rubric_0_100 = (None if not rubric_scores
                    else int(round(_median(
                        [float(s) for s in rubric_scores]))))
    rubric_keys: set[str] = set()
    for r in results:
        rubric_keys.update((r.rubric_components or {}).keys())
    rubric_components = {
        k: _median([float((r.rubric_components or {}).get(k, 0.0))
                    for r in results if r.rubric_components])
        for k in rubric_keys
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

    # One excerpt, not N -- but it has to be an excerpt of the thing the
    # aggregate REPORTS.
    #
    # score_outcome caps a passing sample's excerpt at 400 chars and a
    # failing one at 4000, on the stated ground that "a FAILING one has
    # to be diagnosable from the record alone". Taking the first
    # non-empty one threw that away whenever sample 1 happened to pass:
    # the record then said FAIL and carried 400 characters of a
    # DIFFERENT, passing answer, and the sample that actually failed was
    # unrecoverable. Hit on 2026-09-10 diagnosing
    # science_a_quantity_that_is_absent_is_reported_as_absent, whose
    # samples are bimodal -- the record read `q=51 rate=0.40` beside the
    # excerpt of a run that scored 98.
    #
    # So: when the aggregate does not pass, prefer a sample that did not
    # pass. Still one copy, still the first match, now of the right
    # sample.
    _aggregate_passes = (n_pass * 2 >= n)
    _preferred = ([r for r in results if not r.success]
                  if not _aggregate_passes else [])
    excerpt = ""
    for r in (_preferred + list(results)):
        if r.text_excerpt:
            excerpt = r.text_excerpt
            break
    tool_names_union: list[str] = []
    for r in results:
        for n_name in (r.tool_names or []):
            if n_name and n_name not in tool_names_union:
                tool_names_union.append(n_name)
    # The subjects come from ONE sample, not a union: they are a
    # sequence, and interleaving three replicates would read as a route
    # nobody took. The first sample that made any call is the one kept,
    # which matches how the excerpt is chosen two blocks up.
    subjects: list[str] = []
    for r in (_preferred + list(results)):
        if r.tool_subjects:
            subjects = list(r.tool_subjects)
            break

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
    # Acceptance quality per replicate; a sample that reported none
    # stays None in the list and out of the median.
    _vq_runs = [r.verify_quality for r in results]
    _vq_seen = [q for q in _vq_runs if q is not None]

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
        rubric_0_100=rubric_0_100,
        rubric_components=rubric_components,
        duration_s=_median(durations),
        cost_usd=sum(float(r.cost_usd) for r in results),
        input_tokens=int(_median([float(x) for x in in_toks])),
        output_tokens=int(_median([float(x) for x in out_toks])),
        tool_calls=int(_median([float(x) for x in tool_counts])),
        matched_signals=_matched_every_time(_union("missing_signals")),
        violated_signals=_union("violated_signals"),
        missing_signals=_union("missing_signals"),
        signal_evidence={k: v for r in results
                         for k, v in (r.signal_evidence or {}).items()},
        flaky_signals=_flaky(_union("missing_signals")),
        budget_violations=_union("budget_violations"),
        error="; ".join(sorted({r.error for r in results if r.error}))[:500],
        n_samples=n,
        quality_stdev=_stdev([float(q) for q in qualities]),
        success_rate=n_pass / n,
        success_ci_low=wilson_interval(n_pass, n)[0],
        success_ci_high=wilson_interval(n_pass, n)[1],
        per_run_quality=list(qualities),
        per_run_success=list(success_flags),
        pass_at_k=any(success_flags),
        pass_pow_k=all(success_flags),
        per_run_tool_calls=list(tool_counts),
        per_run_tokens_total=[i + o for i, o in zip(in_toks, out_toks)],
        text_excerpt=excerpt,
        tool_names=tool_names_union,
        tool_subjects=subjects,
        behavior=behavior_agg,
        process={k: round(sum(float((r.process or {}).get(k, 0.0))
                              for r in results) / n, 3)
                 for k in sorted({k for r in results
                                  for k in (r.process or {})})},
        model_rounds=(None if any(r.model_rounds is None for r in results)
                      else int(_median([float(r.model_rounds)
                                        for r in results]))),
        channels={k: sum(int((r.channels or {}).get(k, 0)) for r in results)
                  for k in sorted({k for r in results
                                   for k in (r.channels or {})})},
        unmeasured=all_unmeasured,
        verify_quality=(None if not _vq_seen
                        else int(round(_median([float(q) for q in _vq_seen])))),
        per_run_verify_quality=_vq_runs,
        caveats=int(_median([float(r.caveats) for r in results])),
        answer_chars=int(_median([float(r.answer_chars) for r in results])),
        user_interventions=int(_median(
            [float(r.user_interventions) for r in results])),
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


# ---------------------------------------------------------------------------
# Provenance: what produced a run file
# ---------------------------------------------------------------------------
#
# A run file used to hold result rows only. Two files could then be
# compared although the scorer between them had changed, the code under
# test was not what the file name suggested, or one arm read a different
# setting (auto_verify "syntax" on this host, "smart" shipped). The header
# row says all three; ``provenance_check`` compares two of them.

_HEADER_KEY = "run_header"
_PACKAGE_DIR = Path(__file__).resolve().parent
#: The instrument: the scorer, the runner, the shared recognisers and the
#: whole task pack (task files, acceptance and setup scripts). A change to
#: any of them changes what a score means.
_INSTRUMENT_MODULES = ("benchmark.py", "benchmark_runner.py",
                       "benchmark_fixtures.py", "verify_commands.py",
                       "channels.py")
_SECRET_NAME_RE = re.compile(r"(?i)key|token|secret|passw|auth|cred|cookie")


def instrument_digest(package_dir: Path | None = None) -> dict[str, Any]:
    """sha256 over the instrument files (path + content, sorted).

    Output: ``{"sha256": hex, "files": n}``. Never raises; an unreadable
    file contributes its path and the marker ``<unreadable>``."""
    import hashlib
    base = Path(package_dir) if package_dir else _PACKAGE_DIR
    paths = [base / m for m in _INSTRUMENT_MODULES]
    pack = base / "pack" / "benchmark"
    if pack.is_dir():
        paths += sorted(p for p in pack.rglob("*")
                        if p.is_file() and "__pycache__" not in p.parts)
    h = hashlib.sha256()
    n = 0
    for p in sorted(paths):
        try:
            rel = str(p.relative_to(base))
        except ValueError:
            rel = p.name
        h.update(rel.encode() + b"\0")
        try:
            h.update(p.read_bytes())
            n += 1
        except OSError:
            h.update(b"<unreadable>")
        h.update(b"\0")
    return {"sha256": h.hexdigest(), "files": n}


def _git_head(package_dir: Path) -> dict[str, Any]:
    out: dict[str, Any] = {"head": "", "dirty": None}
    try:
        top = package_dir.parents[1]
        r = safe_git.run(["rev-parse", "HEAD"], top, timeout=10)
        out["head"] = r.stdout.strip() if r.returncode == 0 else ""
        r = safe_git.run(["status", "--porcelain", "--untracked-files=no",
                          "--", "delfin"], top, timeout=20)
        if r.returncode == 0:
            out["dirty"] = bool(r.stdout.strip())
    except Exception:
        pass
    return out


def _process_start_time() -> Optional[float]:
    """Epoch seconds this process started; None where it cannot be told.

    The start ticks come from proc_identity.start_ticks, the one reader of
    a process start time (tests/test_one_reader_of_a_process_start_time);
    the boot time and the tick rate turn them into a date."""
    try:
        from . import proc_identity as _pi
        ticks = _pi.start_ticks(os.getpid())
        if ticks is None:
            return None
        hz = os.sysconf("SC_CLK_TCK")
        for line in Path("/proc/stat").read_text().splitlines():
            if line.startswith("btime "):
                return round(int(line.split()[1]) + ticks / hz, 1)
    except Exception:
        return None
    return None


def _scrub(value: Any, depth: int = 0) -> Any:
    """A settings value with secret-named keys and long strings dropped."""
    if isinstance(value, dict):
        if depth > 3:
            return "<nested>"
        return {str(k): ("<set>" if _SECRET_NAME_RE.search(str(k)) and v
                         else _scrub(v, depth + 1))
                for k, v in sorted(value.items())}
    if isinstance(value, (list, tuple)):
        return [_scrub(v, depth + 1) for v in list(value)[:20]]
    if isinstance(value, str) and len(value) > 200:
        return value[:200] + "…"
    return value


def effective_behaviour_settings() -> dict[str, Any]:
    """What changes the harness's behaviour in this process.

    ``agent`` (the user settings section, secrets dropped -- it carries
    auto_verify, bash_isolation, caged_commands_run_unasked and the
    model_overrides), ``pinned`` (what the benchmark runner exports for
    each attempt, which wins over the setting) and ``env`` (every
    ``DELFIN_*`` variable; values of secret-named ones replaced by
    ``<set>``). Never raises."""
    out: dict[str, Any] = {"agent": {}, "pinned": {}, "env": {}}
    try:
        from delfin.user_settings import load_settings
        out["agent"] = _scrub((load_settings() or {}).get("agent") or {})
    except Exception:
        pass
    try:
        from .benchmark_runner import pinned_behaviour_settings
        out["pinned"] = dict(pinned_behaviour_settings())
    except Exception:
        pass
    out["env"] = {k: ("<set>" if _SECRET_NAME_RE.search(k) else
                      _scrub(v))
                  for k, v in sorted(os.environ.items())
                  if k.startswith("DELFIN_")}
    return out


def _instrument_paths() -> list[str]:
    base = _PACKAGE_DIR
    return [str(base / m) for m in _INSTRUMENT_MODULES] + [
        str(base / "pack" / "benchmark")]


def run_provenance(*, model: str = "") -> dict[str, Any]:
    """The header row of a run file.

    Keys: ``schema``, ``model``, ``written_at``, ``process_start``,
    ``git`` (head, dirty; plus toplevel and dirty_hash), ``instrument``
    (instrument_digest), ``settings`` (effective_behaviour_settings), and
    the generic stamp fields ``taken_at``, ``host``, ``platform``,
    ``machine``, ``python``, ``code`` -- the shape of
    paired_compare.provenance (experiment-discipline bundle), so a
    benchmark run and a generic paired comparison carry one stamp. When
    that module is installed its stamp is the base (code hash over the
    instrument files, DELFIN_* behaviour env redacted) and the keys above
    are added; without it the same fields are filled here."""
    import platform as _platform
    import socket as _socket
    base: dict[str, Any] = {}
    try:
        from . import paired_compare as _pc  # type: ignore[attr-defined]
        base = dict(_pc.provenance(code=_instrument_paths(),
                                   prefixes=("DELFIN_",),
                                   cwd=str(_PACKAGE_DIR.parents[1])))
    except Exception:
        base = {}
    git = _git_head(_PACKAGE_DIR)
    if not base:
        top = str(_PACKAGE_DIR.parents[1])
        dirty_hash = ""
        try:
            import hashlib
            d = safe_git.run(["diff", "HEAD"], top, timeout=20)
            if d.returncode == 0 and d.stdout:
                dirty_hash = hashlib.sha256(d.stdout.encode()).hexdigest()[:20]
        except Exception:
            pass
        base = {
            "taken_at": round(time.time(), 1),
            "host": _socket.gethostname(),
            "platform": _platform.platform(),
            "machine": _platform.machine(),
            "python": sys.version.split()[0],
            "git": {"toplevel": top, "head": git.get("head", ""),
                    "dirty_hash": dirty_hash},
            "code": instrument_digest(),
        }
    stamp = dict(base)
    stamp["git"] = {**(base.get("git") or {}), **git}
    stamp.update({
        "schema": 1,
        "stamp_schema": base.get("schema"),
        "model": model,
        "written_at": round(time.time(), 1),
        "process_start": _process_start_time(),
        "instrument": instrument_digest(),
        "settings": effective_behaviour_settings(),
    })
    return stamp


def read_run_header(path: Path) -> Optional[dict]:
    """The provenance header of a run file, or None (an older file)."""
    try:
        with Path(path).open(encoding="utf-8") as fh:
            first = fh.readline().strip()
    except OSError:
        return None
    if not first.startswith("{"):
        return None
    try:
        row = json.loads(first)
    except json.JSONDecodeError:
        return None
    head = row.get(_HEADER_KEY) if isinstance(row, dict) else None
    return head if isinstance(head, dict) else None


def _flatten(d: Any, prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(_flatten(v, f"{prefix}{k}."))
    else:
        out[prefix[:-1]] = d
    return out


def provenance_check(path_a: Path, path_b: Path) -> dict[str, Any]:
    """Whether two run files were produced by the same instrument.

    Output: ``instrument_match`` (True / False, None when a header is
    missing), ``instrument`` (the two digests), ``heads`` (the two git
    heads), ``settings_diff`` (``{"agent.auto_verify": [a, b], ...}``, every
    effective setting that differs; the declared A/B axis should be the
    only entry) and ``missing`` (the files with no header)."""
    ha, hb = read_run_header(path_a), read_run_header(path_b)
    missing = [str(p) for p, h in ((path_a, ha), (path_b, hb)) if h is None]
    if missing:
        return {"instrument_match": None, "missing": missing,
                "instrument": {}, "heads": {}, "settings_diff": {}}
    ia = (ha.get("instrument") or {}).get("sha256")
    ib = (hb.get("instrument") or {}).get("sha256")
    fa, fb = _flatten(ha.get("settings") or {}), _flatten(hb.get("settings") or {})
    diff = {k: [fa.get(k), fb.get(k)] for k in sorted(set(fa) | set(fb))
            if fa.get(k) != fb.get(k)}
    return {
        "instrument_match": bool(ia) and ia == ib,
        "instrument": {"a": ia, "b": ib},
        "heads": {"a": (ha.get("git") or {}).get("head"),
                  "b": (hb.get("git") or {}).get("head")},
        "settings_diff": diff,
        "missing": [],
    }


def format_provenance_check(chk: dict) -> list[str]:
    """Lines for a human; the first says whether the comparison holds."""
    if chk.get("instrument_match") is None:
        return ["PROVENANCE UNKNOWN: no header in "
                + ", ".join(chk.get("missing") or []) +
                " -- the instrument and settings of that arm are not "
                "recorded; read the comparison as unverified."]
    lines = []
    if not chk["instrument_match"]:
        lines.append("INSTRUMENT CHANGED between the arms "
                     f"({chk['instrument'].get('a', '')[:12]} vs "
                     f"{chk['instrument'].get('b', '')[:12]}): the scorer or "
                     "the task pack differs, so a score difference may be "
                     "the measure, not the code.")
    else:
        lines.append("instrument: identical in both arms "
                     f"({str(chk['instrument'].get('a'))[:12]})")
    lines.append(f"code: {str(chk['heads'].get('a'))[:12]} vs "
                 f"{str(chk['heads'].get('b'))[:12]}")
    if chk["settings_diff"]:
        lines.append("settings that differ (should be exactly the declared "
                     "axis):")
        for k, (a, b) in chk["settings_diff"].items():
            lines.append(f"  {k}: {a!r} -> {b!r}")
    else:
        lines.append("settings: identical")
    return lines


def write_run(
    results: Iterable[BenchmarkResult],
    *,
    model: str,
    runs_dir: Path | None = None,
    run_id: str | None = None,
    header: Optional[dict] = None,
    run_args: Optional[dict] = None,
) -> Path:
    """Persist all results from one benchmark run as JSONL.

    The first line is ``{"run_header": <run_provenance>}`` (``header``
    overrides it, for a caller that measured its own); every further line
    is one result. ``read_run`` skips the header; ``read_run_header``
    returns it."""

    ts = int(time.time())
    rid = run_id or f"{ts}_{_slug(model)}"
    rows = [json.dumps(asdict(r), ensure_ascii=False) for r in results]
    try:
        _head = header if header is not None else run_provenance(model=model)
        if run_args:
            # The run's own arguments (max_tokens, repeats, ...): they
            # change what an attempt can do and belong to the provenance
            # comparison like the settings.
            _head = {**_head, "settings": {**(_head.get("settings") or {}),
                                           "run_args": dict(run_args)}}
        rows.insert(0, json.dumps({_HEADER_KEY: _head}, ensure_ascii=False,
                                  default=str))
    except Exception:
        pass

    # Where the results go, in order of preference. The default is under
    # the user's home; a run started under a redirected or read-only home
    # cannot write there, and used to die at this line AFTER doing all of
    # the work -- measured 2026-09-28 on a login-node trial:
    #
    #   OSError: [Errno 30] Read-only file system:
    #   '.../.delfin/benchmark_runs/...jsonl'
    #
    # The run was over, every task had been scored, and the file that
    # holds the scores was the one thing that failed. Losing a completed
    # measurement to the place it is filed is the worst possible trade,
    # so a refused directory is stepped over rather than raised on.
    #
    # An EXPLICIT runs_dir is not stepped over: a caller that named a
    # directory wants that directory, and quietly writing somewhere else
    # would hide the fault from whoever is collecting the files.
    candidates = [runs_dir] if runs_dir is not None else [_DEFAULT_RUNS_DIR]
    if runs_dir is None:
        scratch = (os.environ.get("DELFIN_SCRATCH_STATE") or "").strip()
        if scratch:
            candidates.append(Path(scratch) / ".delfin" / "benchmark_runs")
        candidates.append(Path(tempfile.gettempdir()) / "delfin_benchmark_runs")

    refused: list[str] = []
    for d in candidates:
        path = d / f"{rid}.jsonl"
        try:
            d.mkdir(parents=True, exist_ok=True)
            # A run holds the prompts it sent and the text the model
            # answered. The last resort here is the system temp dir, which
            # on a shared login node is everybody's: created 0700 and
            # written 0600, so stepping over a refused directory never
            # trades the file's privacy for its survival. Measured: the
            # defaults are 0775 and 0664.
            try:
                os.chmod(d, 0o700)
            except OSError:
                pass
            with path.open("w", encoding="utf-8") as f:
                try:
                    os.chmod(path, 0o600)
                except OSError:
                    pass
                for row in rows:
                    f.write(row + "\n")
        except OSError as exc:
            refused.append(f"{d}: {exc.strerror or exc}")
            continue
        if refused:
            # Said, not swallowed: a reader looking for the file where it
            # has always been needs to be told once where it went.
            print(f"[benchmark] results written to {path} "
                  f"(could not use {'; '.join(refused)})", file=sys.stderr)
        return path
    raise OSError(
        "benchmark results could not be written anywhere; tried "
        + "; ".join(refused))


def read_run(path: Path) -> list[dict]:
    if not path.exists():
        return []
    out: list[dict] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(row, dict) and _HEADER_KEY in row and len(row) == 1:
            continue
        out.append(row)
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
            "n_pass_at_k": 0, "pass_at_k_rate": 0.0,
            "n_pass_pow_k": 0, "pass_pow_k_rate": 0.0,
            "n_unmeasured": 0, "unmeasured_tasks": [],
            "avg_quality": 0.0, "total_cost_usd": 0.0,
            "total_duration_s": 0.0, "total_tool_calls": 0,
            "total_input_tokens": 0, "total_output_tokens": 0,
            "per_call_input_tokens": None, "per_call_output_tokens": None,
            "avg_caveats": 0.0, "max_caveats": 0,
            "values_matched": 0, "values_wrong": 0, "values_absent": 0,
            "avg_output_tokens": 0.0, "avg_answer_chars": 0.0,
            "total_denials": None,
            "process_means": {}, "channels_total": {},
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
    # The two stricter, alongside readings over the SAME k repeats. A task
    # whose request never reached the model is already out (``scored``),
    # exactly as it is out of ``n_pass``. pass@k = >=1 of its repeats
    # passed; pass^k = ALL of its repeats passed. Neither rewrites the
    # existing pass criterion — they are reported next to it so a run's
    # "best attempt" and its "every attempt" are both visible.
    n_pass_at_k = sum(1 for r in scored if r.get("pass_at_k"))
    n_pass_pow_k = sum(1 for r in scored if r.get("pass_pow_k"))
    # Token sums over scored tasks, and the per-call normalisation that
    # tells a real saving from "did less work" (see the per-call comment
    # below).  Exact values, no rounding: the caller renders them.
    total_input_tokens = sum(
        max(0, int(r.get("input_tokens") or 0)) for r in scored)
    total_output_tokens = sum(
        max(0, int(r.get("output_tokens") or 0)) for r in scored)
    scored_calls = sum(
        max(0, int(r.get("tool_calls") or 0)) for r in scored)
    return {
        "n_tasks": n,
        "n_pass": n_pass,
        "n_unmeasured": n_unmeasured,
        "unmeasured_tasks": sorted(
            str(r.get("task_id") or "?") for r in rows if r.get("unmeasured")),
        "pass_rate": (n_pass / n) if n else 0.0,
        "n_pass_at_k": n_pass_at_k,
        "pass_at_k_rate": (n_pass_at_k / n) if n else 0.0,
        "n_pass_pow_k": n_pass_pow_k,
        "pass_pow_k_rate": (n_pass_pow_k / n) if n else 0.0,
        "avg_quality": (sum(int(r.get("quality_0_100") or 0)
                            for r in scored) / n) if n else 0.0,
        "total_cost_usd": sum(float(r.get("cost_usd") or 0) for r in rows),
        "total_duration_s": sum(float(r.get("duration_s") or 0) for r in rows),
        "total_tool_calls": scored_calls,
        # The token side, and the per-call normalisation that tells a real
        # saving ("cheaper per call") from a fake one ("did less work":
        # fewer calls, same or higher tokens each).  A total that shrank
        # while the per-call number held or rose is not a saving -- the
        # wave-10 case.  No per-model-call count is recorded anywhere in
        # the run record, so "per tool call" is the honest cell; dividing
        # by a model-call count that is not there would be fabrication.
        "total_input_tokens": total_input_tokens,
        "total_output_tokens": total_output_tokens,
        "per_call_input_tokens": (
            None if not scored_calls else total_input_tokens / scored_calls),
        "per_call_output_tokens": (
            None if not scored_calls else total_output_tokens / scored_calls),
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
        # Figures asked for and not matched, split by WHY. A wrong figure
        # is the model's arithmetic; a figure that is absent from an
        # answer usually means the answer never got here -- an empty
        # turn, a truncation, a refusal -- and reading the two as one
        # number is how a harness fault gets filed as a model fault.
        "values_matched": _count_verdicts(scored, "matched"),
        "values_wrong": _count_verdicts(scored, "wrong"),
        "values_absent": _count_verdicts(scored, "absent"),
        # Harness steering observed during the runs: permission denials.
        # Unlike total_denials this is never None -- the runner watches
        # every turn in-stream.
        "total_interventions": sum(
            int(r.get("user_interventions") or 0) for r in scored),
        # The process signals as a mean per key over the scored tasks that
        # carry them (``n`` beside each), and the channel firings summed:
        # the per-task lines alone gave no rate for a model.
        "process_means": _process_means(scored),
        "channels_total": _channels_total(scored),
    }


def _process_means(rows: list[dict]) -> dict[str, dict[str, float]]:
    sums: dict[str, float] = {}
    counts: dict[str, int] = {}
    for r in rows:
        for k, v in (r.get("process") or {}).items():
            try:
                sums[k] = sums.get(k, 0.0) + float(v)
            except (TypeError, ValueError):
                continue
            counts[k] = counts.get(k, 0) + 1
    return {k: {"mean": round(sums[k] / counts[k], 3), "n": counts[k]}
            for k in sorted(sums)}


def _channels_total(rows: list[dict]) -> dict[str, int]:
    out: dict[str, int] = {}
    for r in rows:
        for k, v in (r.get("channels") or {}).items():
            try:
                out[k] = out.get(k, 0) + int(v)
            except (TypeError, ValueError):
                continue
    return dict(sorted(out.items()))


def _count_verdicts(rows: list[dict], verdict: str) -> int:
    total = 0
    for row in rows:
        report = row.get("value_report") or {}
        if isinstance(report, dict):
            total += sum(1 for v in report.values() if v == verdict)
    return total


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

    provenance = (provenance_check(baseline, candidate)
                  if isinstance(baseline, Path) and isinstance(candidate, Path)
                  else None)
    if isinstance(baseline, Path):
        baseline = read_run(baseline)
    if isinstance(candidate, Path):
        candidate = read_run(candidate)

    # Group rows by task_id into lists, so a task written more than once
    # (repeats un-aggregated in one run file) is summed rather than
    # silently collapsed by a last-wins dict.  `dup_*` names every task
    # seen more than once so the compare can say so instead of lying
    # about certainty.
    old_rows: dict = {}
    for r in baseline:
        old_rows.setdefault(r.get("task_id"), []).append(r)
    new_rows: dict = {}
    for r in candidate:
        new_rows.setdefault(r.get("task_id"), []).append(r)
    overlap = sorted(set(old_rows) & set(new_rows))
    dup_old = sorted(t for t, rs in old_rows.items() if len(rs) > 1)
    dup_new = sorted(t for t, rs in new_rows.items() if len(rs) > 1)

    per_task: list[dict] = []
    n_better = n_worse = n_neutral = 0
    for tid in overlap:
        o = old_rows[tid][-1]
        n = new_rows[tid][-1]
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
        oc = _aggregate_counts(old_rows[tid])
        nc = _aggregate_counts(new_rows[tid])
        _delta = classify_pass_delta(*oc, *nc)
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
            "success_rate_old": oc[0] / max(1, oc[1]),
            "success_rate_new": nc[0] / max(1, nc[1]),
            # Whether THIS task's pass-rate swing is more than its own
            # sampling noise.  Direction (old/new) and certainty (noise
            # label) are kept apart: `class` says which way it moved, the
            # label says whether N can trust that move.
            "success_counts_old": oc,
            "success_counts_new": nc,
            "delta_noise": _delta["label"],
            "delta_p": _delta["p"],
        })

    # A-vs-B significance: pool the per-task pass counts and ask Fisher
    # whether the observed direction could be run-to-run noise.  This is
    # the field a change-accept decision should read FIRST -- the
    # threshold verdicts above describe direction, not certainty.
    pooled_a_pass = pooled_a_n = pooled_b_pass = pooled_b_n = 0
    p_values: list[float] = []
    for tid in overlap:
        a_pass, a_n = _aggregate_counts(old_rows[tid])
        b_pass, b_n = _aggregate_counts(new_rows[tid])
        pooled_a_pass += a_pass
        pooled_a_n += a_n
        pooled_b_pass += b_pass
        pooled_b_n += b_n
        if a_n or b_n:
            p_values.append(_fisher_exact_2x2_pvalue(
                a_pass, a_n - a_pass, b_pass, b_n - b_pass))
    sig_p = (_fisher_exact_2x2_pvalue(
        pooled_a_pass, pooled_a_n - pooled_a_pass,
        pooled_b_pass, pooled_b_n - pooled_b_pass)
        if (pooled_a_n or pooled_b_n) else 1.0)
    significant = bool(p_values) and sig_p < _SIG_ALPHA

    summary_old = summarise_run(baseline)
    summary_new = summarise_run(candidate)
    summary = {
        "n_overlap": len(overlap),
        "old": summary_old,
        "new": summary_new,
        "n_better": n_better,
        "n_worse": n_worse,
        "n_neutral": n_neutral,
        # Task ids seen more than once in a run; their repeats were pooled
        # rather than silently dropped.  An EMPTY list is the honest state.
        "duplicate_task_ids": {"old": dup_old, "new": dup_new},
        # Statistical statement over the pooled success counts.  A
        # verdict like "worse" from the thresholds above means nothing
        # when this field says the difference is inside the noise.
        "significant": significant,
        "significance_p": round(sig_p, 4),
        "significance_alpha": _SIG_ALPHA,
        "pooled_success": {
            "old": [pooled_a_pass, pooled_a_n],
            "new": [pooled_b_pass, pooled_b_n],
        },
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

    out = {"per_task": per_task, "summary": summary, "verdict": verdict}
    if provenance is not None:
        # Two files on disk: whether the same instrument produced both,
        # and every setting that differs (provenance_check).
        out["provenance"] = provenance
    return out


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
    out = safe_git.run(
        ["log", "--oneline", f"--since={old_iso}", f"--until={new_iso}",
         "--", profile_file], cwd, timeout=timeout_s)
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
    oit = int(old.get("total_input_tokens") or 0)
    nit = int(new.get("total_input_tokens") or 0)
    oot = int(old.get("total_output_tokens") or 0)
    not_ = int(new.get("total_output_tokens") or 0)
    lines.append(
        f"| In tokens    | {oit:>11,} | {nit:>11,} | {(nit - oit):+11,} |"
    )
    lines.append(
        f"| Out tokens   | {oot:>11,} | {not_:>11,} | {(not_ - oot):+11,} |"
    )

    # Per-call normalisation: the figure that reveals "did less work" vs
    # "cheaper per call".  None (no tool calls) prints as n/a rather than
    # a crash.
    def _pc(v: Any) -> str:
        return f"{float(v):>10,.1f}" if v is not None else f"{'n/a':>10}"

    def _pcd(old_v: Any, new_v: Any) -> str:
        if old_v is None or new_v is None:
            return f"{'n/a':>10}"
        return f"{float(new_v) - float(old_v):+10,.1f}"

    _oin = old.get("per_call_input_tokens")
    _nin = new.get("per_call_input_tokens")
    _oout = old.get("per_call_output_tokens")
    _nout = new.get("per_call_output_tokens")
    lines.append(
        f"| In tok/call  | {_pc(_oin)} | {_pc(_nin)} | {_pcd(_oin, _nin)} |"
    )
    lines.append(
        f"| Out tok/call | {_pc(_oout)} | {_pc(_nout)} | {_pcd(_oout, _nout)} |"
    )
    lines.append("")

    # A task id written more than once in a run is pooled, not silently
    # dropped; say so so a reader does not mistake the accounting.
    _dups = summary.get("duplicate_task_ids") or {"old": [], "new": []}
    if _dups.get("old") or _dups.get("new"):
        lines.append(
            f"**Duplicate task ids (pooled, not collapsed):** old "
            f"{', '.join(map(str, _dups.get('old') or [])) or '-'}, "
            f"new {', '.join(map(str, _dups.get('new') or [])) or '-'}")
        lines.append("")

    # Significance: is the change real, or within the run's own sampling
    # noise?  The verdict line says WHICH WAY it moved; this section says
    # whether N can trust that direction.  Empty overlap → nothing to say.
    _ov = summary.get("n_overlap") or 0
    _ps = summary.get("pooled_success") or {}
    _oldc = _ps.get("old") or [0, 0]
    _newc = _ps.get("new") or [0, 0]
    _sp = summary.get("significance_p")
    # A sub-threshold p that rounded to 0.0 must NOT read as absent: the
    # `or 1.0` fallback turns the truest "noise-crossed" verdict into
    # p=1 beside a "beyond noise" label.  Only a genuinely missing field
    # falls back to p=1 ("we did not compute one").
    _p = float(_sp) if _sp is not None else 1.0
    _alpha = float(summary.get("significance_alpha") or _SIG_ALPHA)
    if _ov and int(_oldc[1]) and int(_newc[1]):
        _plab = "beyond noise" if summary.get("significant") else "within noise"
        _olo, _ohi = wilson_interval(int(_oldc[0]), int(_oldc[1]))
        _nlo, _nhi = wilson_interval(int(_newc[0]), int(_newc[1]))
        lines.append("## Significance (is this real or noise?)")
        lines.append("")
        lines.append(
            f"Pooled success {int(_oldc[0])}/{int(_oldc[1])} -> "
            f"{int(_newc[0])}/{int(_newc[1])}: two-sided Fisher exact "
            f"p={_p:.4g} at alpha={_alpha:.2f} -> **{_plab}**"
        )
        lines.append(
            f"  old pass rate Wilson 95% CI (n={int(_oldc[1])}): "
            f"{_olo:.0%} - {_ohi:.0%}"
        )
        lines.append(
            f"  new pass rate Wilson 95% CI (n={int(_newc[1])}): "
            f"{_nlo:.0%} - {_nhi:.0%}"
        )
        if not summary.get("significant") and int(_oldc[1]) and int(_newc[1]):
            # Not significant is not the end of the story: name the N that
            # would settle the observed pooled gap, so a reader who wants a
            # verdict knows what to run next.  A zero gap has no resolving
            # N; a small-but-real gap may just need more than the search
            # covers -- saying "zero" for that would mislead.
            _or = int(_oldc[0]) / int(_oldc[1])
            _nr = int(_newc[0]) / int(_newc[1])
            if _or == _nr:
                lines.append(
                    f"  N too small to decide: the observed pooled gap is "
                    f"zero, so no repeat count can separate these runs.")
            else:
                _need = min_n_for_delta(max(_or, _nr), min(_or, _nr))
                if _need is not None:
                    lines.append(
                        f"  N too small to decide: resolving a swing this "
                        f"large needs ~{_need} repeats per arm.")
                else:
                    lines.append(
                        f"  N too small to decide: this small a gap is not "
                        f"resolvable within the repeats searched so far; "
                        f"run far more repeats per arm or test a stronger "
                        f"effect.")
        lines.append("")

    # Per-task table
    per_task = cmp_result.get("per_task") or []
    if per_task:
        lines.append("## Per-task")
        lines.append("")
        lines.append(
            "| Task                          | Class   | Quality | Δq   | "
            "Δcost     | Δdur      | Δpass      | noise  |"
        )
        lines.append(
            "|-------------------------------|---------|---------|------|"
            "-----------|-----------|------------|--------|"
        )
        for row in per_task:
            mark = _CLASS_MARK.get(row.get("class") or "", "[?]")
            tid = str(row.get("task_id") or "")[:30]
            oq2 = int(row.get("old_quality") or 0)
            nq2 = int(row.get("new_quality") or 0)
            dq = int(row.get("d_quality") or 0)
            dc = float(row.get("d_cost_usd") or 0)
            dd = float(row.get("d_duration_s") or 0)
            o_rate = float(row.get("success_rate_old") or 0)
            n_rate = float(row.get("success_rate_new") or 0)
            _opass = row.get("delta_noise")
            if _opass not in ("within", "beyond"):
                # A row from an older caller (or a hand-built cmp_result)
                # carries no noise label.  "within" as a silent default
                # asserts a certainty the data never gave; render unknown.
                _opass, _noise_mark = "unknown", "[?] "
            else:
                _noise_mark = "[!] " if _opass == "beyond" else "[~] "
            lines.append(
                f"| `{tid:<29}` | {mark:<7} | {oq2:>3}→{nq2:<3} | "
                f"{dq:+5d} | {dc:+9.4f} | {dd:+8.2f}s | "
                f"{o_rate:.0%}→{n_rate:<3.0%} | {_noise_mark}{_opass:<7} |"
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
        "provenance": provenance_check(path_a, path_b),
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


def _learning_pair_mismatch(first: dict, second: dict, *,
                            require_shared_state: bool) -> str:
    """Name what differs between two rows that must not differ, or ''."""
    if not str(first.get("model") or "") or not str(second.get("model") or ""):
        return "model: a run record without a model cannot be paired"
    fields = ["model", "profile_name", "mode"]
    if require_shared_state:
        if not first.get("learning_state") or not second.get("learning_state"):
            return ("learning_state: a run made without --learning-home "
                    "started from an empty scratch home, so nothing could "
                    "be learned between the two runs")
        fields.append("learning_state")
    for f in fields:
        if (first.get(f) or "") != (second.get(f) or ""):
            return (f"{f} differs ({first.get(f)!r} vs {second.get(f)!r}); "
                    "a learning pair holds it fixed")
    return ""


def compare_learning(
    run_a: str,
    run_b: str,
    *,
    runs_dir: Path | None = None,
    require_shared_state: bool = True,
) -> dict[str, Any]:
    """Per-task learning verdict: did the SECOND run of a task beat the first?

    A pair is compared only when the two rows held fixed what a learning
    comparison must hold fixed: the same non-empty ``model``, the same
    ``profile_name`` and ``mode``, and -- unless ``require_shared_state``
    is False -- the same non-empty ``learning_state``, i.e. both runs were
    made under one ``--learning-home``. Any difference is an ``error`` row
    naming the field. Without a shared state the second run started from
    the same empty scratch home as the first and nothing could have been
    learned; ``require_shared_state=False`` is the null control of
    PREREG_W15.md (two runs under different homes), never a learning
    result.

    ``run_a`` is the FIRST run in time, ``run_b`` the SECOND — each a full
    path, a filename, or a bare run id resolved against the runs directory
    (the same resolution :func:`ab_compare` uses). A learning comparison is
    first-versus-second across TWO runs separated in time (two run files);
    two rows of the same task within ONE file are replicates of that
    snapshot and are the same-snapshot spread, not learning. The LAST row
    is a run's representative per task (the same convention
    :func:`compare_runs` uses).

    For every task present in BOTH runs the FIRST representative row is
    thrown through
    :func:`delfin.agent.pack.benchmark.learning_metric.learning_verdict`
    with the SECOND, so the verdict per task is one of
    ``error | wrong | learned | no_learning | regressed`` (PREREG_W15.md).
    Keeps ``learning_metric`` pure: THIS caller decides which record is
    first and which is second. Returns::

        {run_a, run_b,
         per_task: [{task_id, verdict, error, first, second,
                     n_first, n_second}],
         summary: {n_learned, n_no_learning, n_regressed, n_wrong, n_error,
                   n_only_in_a, n_only_in_b, only_in_a, only_in_b}}

    A task seen in only one run is excluded from ``per_task`` (it has no
    counterpart to judge) and reported in the summary as ``only_in_a`` /
    ``only_in_b``.
    """
    from .pack.benchmark import learning_metric as _learning_metric

    base = runs_dir or _DEFAULT_RUNS_DIR
    path_a = _resolve_run_path(run_a, base)
    path_b = _resolve_run_path(run_b, base)
    rec_a = read_run(path_a)
    rec_b = read_run(path_b)

    rows_a: dict[str, list[dict]] = {}
    for r in rec_a:
        tid = r.get("task_id")
        if isinstance(tid, str):
            rows_a.setdefault(tid, []).append(r)
    rows_b: dict[str, list[dict]] = {}
    for r in rec_b:
        tid = r.get("task_id")
        if isinstance(tid, str):
            rows_b.setdefault(tid, []).append(r)

    overlap = sorted(set(rows_a) & set(rows_b))
    only_a = sorted(set(rows_a) - set(rows_b))
    only_b = sorted(set(rows_b) - set(rows_a))

    stats: dict[str, int] = {"n_learned": 0, "n_no_learning": 0,
                             "n_regressed": 0, "n_wrong": 0, "n_error": 0}
    per_task: list[dict[str, Any]] = []
    for tid in overlap:
        first = rows_a[tid][-1]
        second = rows_b[tid][-1]
        mismatch = _learning_pair_mismatch(
            first, second, require_shared_state=require_shared_state)
        if mismatch:
            out = {"verdict": "error", "evidence": {}, "error": mismatch}
        else:
            out = _learning_metric.learning_verdict(first, second)
        verdict = str(out.get("verdict") or "error")
        if f"n_{verdict}" in stats:
            stats[f"n_{verdict}"] += 1
        else:
            stats["n_error"] += 1
        evidence = out.get("evidence") or {}
        per_task.append({
            "task_id": tid,
            "verdict": verdict,
            "error": out.get("error"),
            "first": evidence.get("first"),
            "second": evidence.get("second"),
            "n_first": evidence.get("n_first"),
            "n_second": evidence.get("n_second"),
        })
    per_task.sort(key=lambda p: str(p["task_id"]))

    summary = {
        **stats,
        "n_only_in_a": len(only_a),
        "n_only_in_b": len(only_b),
        "only_in_a": only_a,
        "only_in_b": only_b,
    }
    return {
        "run_a": str(path_a),
        "run_b": str(path_b),
        "per_task": per_task,
        "summary": summary,
    }


_LEARNING_TASKS_PATH = _DEFAULT_TASKS_PATH.parent / "tasks_w15_learning.yaml"


def learning_task_ids() -> frozenset[str]:
    """The ids of the pre-registered learning tasks (PREREG_W15.md).

    A learning verdict is only meaningful for a task whose setup lets a
    learned skill or a memory help the second time; on any other task
    "fewer tool calls" between two runs is ordinary run-to-run spread.
    Empty when the suite cannot be read.
    """
    try:
        return frozenset(t.id for t in _load_one(_LEARNING_TASKS_PATH))
    except Exception:
        return frozenset()


def latest_learning_pair(
    runs: list[Path], *, only: Any = None,
) -> tuple[Path, Path] | None:
    """The two most recent arms of one learning measurement, or None.

    ``runs`` is oldest first (:func:`list_runs`). The second arm is the
    newest file holding a row of a task in ``only`` (all tasks when None)
    with a non-empty ``learning_state``; the first arm is the newest
    earlier file whose last row of the same task carries the same state.
    Files in between without such a row -- ordinary runs -- are stepped
    over.
    """
    def _states(path: Path) -> set[tuple[str, str]]:
        try:
            rows = read_run(path)
        except Exception:
            return set()
        last: dict[str, str] = {}
        for r in rows:
            tid = r.get("task_id")
            if isinstance(tid, str) and (only is None or tid in only):
                last[tid] = str(r.get("learning_state") or "")
        return {(t, s) for t, s in last.items() if s}

    for j in range(len(runs) - 1, 0, -1):
        second = _states(runs[j])
        if not second:
            continue
        for i in range(j - 1, -1, -1):
            if _states(runs[i]) & second:
                return runs[i], runs[j]
        return None
    return None


def format_learning_note(lc: dict, *, only: Any = None) -> str:
    """Short plain-text rendering of a :func:`compare_learning` result.

    ``only`` restricts the note to these task ids (the pre-registered
    learning tasks); the counts are recomputed over what is shown, so the
    header never counts a task the lines below leave out. Empty string
    when nothing is left to show.
    """
    rows = [p for p in (lc.get("per_task") or [])
            if only is None or p.get("task_id") in only]
    if not rows:
        return ""
    order = ("learned", "no_learning", "regressed", "wrong", "error")
    counts = {v: 0 for v in order}

    def _shown(p: dict) -> str:
        v = str(p.get("verdict") or "error")
        return v if v in counts else "error"

    for p in rows:
        counts[_shown(p)] += 1
    lines = ["second run vs first: " + ", ".join(
        f"{counts[v]} {v.replace('_', ' ')}" for v in order)]
    for p in rows:
        v = _shown(p)
        first, second = p.get("first") or {}, p.get("second") or {}
        if v == "error" or not first or not second:
            detail = str(p.get("error")
                         or f"no comparable record (verdict "
                            f"{p.get('verdict')!r})")
        else:
            detail = (f"median tool calls {first.get('tool_calls')}→"
                      f"{second.get('tool_calls')}, tokens "
                      f"{first.get('tokens_total')}→"
                      f"{second.get('tokens_total')}, "
                      f"N={p.get('n_first')}/{p.get('n_second')}")
        lines.append(f"{p.get('task_id')}: {v.replace('_', ' ')} ({detail})")
    return "\n".join(lines)


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
    prov = ab.get("provenance")
    if isinstance(prov, dict):
        lines.extend(format_provenance_check(prov))
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
            "signal_evidence": dict(r.get("signal_evidence") or {}),
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
                    # The sentence that tripped it. The excerpt below is a
                    # head slice, so a violation past its end used to be a
                    # pattern label with nothing to judge it by.
                    ev = (e.get("signal_evidence") or {}).get(label)
                    if ev:
                        lines.append(f"      matched : {ev[:400]}")
            if e["error"]:
                lines.append(f"  error  : {e['error'][:200]}")
            if e["text_excerpt"]:
                # The label said 400 while the record has held 4000 for a
                # failing task since someone widened it for exactly this
                # reason. A caption that understates what is there invites
                # the reader to give up before scrolling.
                _lines = e["text_excerpt"].splitlines()
                lines.append(
                    f"  model output ({len(e['text_excerpt'])} chars, "
                    f"first {min(12, len(_lines))} lines):")
                for line in _lines[:12]:
                    lines.append(f"    │ {line[:120]}")
                if len(_lines) > 12:
                    lines.append(f"    │ … ({len(_lines) - 12} more lines)")
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
