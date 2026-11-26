"""A before/after benchmark measurement that states its own noise.

Input: two git refs, a task selection, a repeat count. Output: per task, the
pass rate on each side, the difference between them, and an INTERVAL on that
difference -- plus the interval's containment of zero, read out as
``no measurable difference`` when it holds.

WHY THIS EXISTS. The package apparatus could compare two persisted runs
(``benchmark.ab_compare``) and schedule a nightly round (``bench_watch``), but
nothing ran the SAME task N times on EACH side of a change. A comparison of one
run against one run cannot separate an improvement from noise: an agent that
passes 1 of 1 and an agent that passes 3 of 3 are the same measurement, and the
first would have been reported as the better one. Three times in one wave the
baseline was already right and a single probe said otherwise.

THE FOUR THINGS THAT MAKE A PAIR VALID, and are enforced here rather than
assumed:

* one change            -- the sides are the base and head commits, each its
                           own worktree, so the code differs by exactly the
                           branch and by nothing else;
* one instrument        -- same model, same sampling, same task set, same
                           repeats. A differing model or a differing task TEXT
                           is REFUSED, not averaged (see ``AbRefused``);
* one population        -- the same task ids on both sides, and a task that was
                           not measured on a side is reported as not measured,
                           never as a failure;
* a stated error bar    -- the interval below, printed with its method and the
                           reason it was chosen, because an interval whose
                           method is unstated is not evidence.

WHAT IS NOT CLAIMED. The sampling parameters are ASKED of the endpoint, and an
endpoint may refuse them (``sampling.REFUSED``); the report says so rather than
calling the pair deterministic. Repeats reduce noise, they do not remove it,
and at ``--repeat 3`` the interval is wide -- that width IS the finding.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Optional, Sequence

#: Two-sided interval, the same z ``benchmark.wilson_interval`` uses for
#: ``success_rate``. Kept here as a named constant so the report can state it.
Z = 1.96

#: How a difference of proportions gets its interval. Stated, not implied: the
#: whole point of the verb is that a reader can check the arithmetic.
INTERVAL_METHOD = (
    "Newcombe hybrid score interval (method 10) on the difference of two "
    "Wilson score intervals, z=1.96"
)
INTERVAL_WHY = (
    "Wilson per side because at --repeat 3 a normal-approximation interval "
    "runs outside [0,1] and collapses at 0/3 and 3/3; Newcombe on the "
    "difference because it is closed-form (no RNG, so the same counts give the "
    "same interval -- reproducible without a seed) and is built from the same "
    "Wilson intervals this repository already uses for success_rate. A "
    "bootstrap over three binary outcomes has four distinct values and would "
    "resample the grid, not the data. The price is slight conservatism at tiny "
    "N, which is the direction an under-powered comparison should err."
)

#: The words a covered zero is reported with. A constant because it is a
#: contract: a summary that says "improved" beside an interval containing 0 is
#: the defect this module was written against.
NO_DIFFERENCE = "no measurable difference"

# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def newcombe_interval(pass_a: int, n_a: int, pass_b: int, n_b: int,
                      z: float = Z) -> tuple[float, float]:
    """Interval on ``p_b - p_a`` from the two Wilson intervals (Newcombe 1998).

    ``pass_a``/``n_a`` is the base side, ``pass_b``/``n_b`` the head side, so a
    positive difference means the head passed more. Newcombe's construction
    (method 10 in the paper) is

        l = d - sqrt((p_a - l_a)^2 + (u_b - p_b)^2)
        u = d + sqrt((u_a - p_a)^2 + (p_b - l_b)^2)

    with ``d = p_b - p_a`` and ``[l_x, u_x]`` the Wilson interval of side x.

    An empty side carries no information: the answer is the whole possible
    range of a difference of proportions, ``(-1, 1)``, not a narrower one and
    not a crash. That interval contains 0, so an empty side can never be read
    as a change -- which is the correct statement about a side nobody measured.
    """
    if n_a <= 0 or n_b <= 0:
        return -1.0, 1.0
    from .benchmark import wilson_interval

    la, ua = wilson_interval(int(pass_a), int(n_a), z)
    lb, ub = wilson_interval(int(pass_b), int(n_b), z)
    pa = pass_a / n_a
    pb = pass_b / n_b
    d = pb - pa
    lo = d - math.sqrt((pb - lb) ** 2 + (ua - pa) ** 2)
    hi = d + math.sqrt((ub - pb) ** 2 + (pa - la) ** 2)
    return max(-1.0, lo), min(1.0, hi)


#: A verdict. ``no difference`` is its own word on purpose -- "the difference is
#: 0" and "the difference is not resolvable" are different claims and a run at
#: N=3 almost never supports the first.
VERDICT_IMPROVED = "improved"
VERDICT_REGRESSED = "regressed"
VERDICT_NO_DIFFERENCE = "no difference"
VERDICT_UNMEASURED = "not measured"


@dataclass
class SideResult:
    """What one side's attempts on one task came back with."""

    passes: list[bool] = field(default_factory=list)
    unmeasured: int = 0
    errors: list[str] = field(default_factory=list)
    model: str = ""
    task_digest: str = ""
    module_file: str = ""
    worktree: str = ""

    @property
    def measured(self) -> int:
        return len(self.passes)

    @property
    def n_pass(self) -> int:
        return sum(1 for p in self.passes if p)

    @property
    def rate(self) -> Optional[float]:
        return (self.n_pass / self.measured) if self.measured else None


@dataclass
class TaskAB:
    """One task's paired comparison, with everything the reader needs."""

    task_id: str
    base: SideResult
    head: SideResult
    selection_rule: str = ""

    @property
    def diff(self) -> Optional[float]:
        if self.base.rate is None or self.head.rate is None:
            return None
        return self.head.rate - self.base.rate

    @property
    def interval(self) -> tuple[float, float]:
        return newcombe_interval(self.base.n_pass, self.base.measured,
                                 self.head.n_pass, self.head.measured)

    @property
    def verdict(self) -> str:
        if self.base.rate is None or self.head.rate is None:
            return VERDICT_UNMEASURED
        lo, hi = self.interval
        if lo <= 0.0 <= hi:
            return VERDICT_NO_DIFFERENCE
        return VERDICT_IMPROVED if lo > 0.0 else VERDICT_REGRESSED

    def as_dict(self) -> dict[str, Any]:
        lo, hi = self.interval
        return {
            "task_id": self.task_id,
            "n_base": self.base.measured,
            "n_head": self.head.measured,
            "pass_base": self.base.n_pass,
            "pass_head": self.head.n_pass,
            "rate_base": self.base.rate,
            "rate_head": self.head.rate,
            "diff": self.diff,
            "interval_low": lo,
            "interval_high": hi,
            "contains_zero": lo <= 0.0 <= hi,
            "verdict": self.verdict,
            "unmeasured_base": self.base.unmeasured,
            "unmeasured_head": self.head.unmeasured,
            "task_digest": self.base.task_digest,
            "selection_rule": self.selection_rule,
            "errors": (self.base.errors + self.head.errors)[:6],
        }


class AbRefused(RuntimeError):
    """The two sides are not a pair, so no difference is reported.

    Raised for a differing MODEL and for a differing TASK (same id, a different
    question). Each of those makes the difference between the sides the wrong
    difference: it would be read as the effect of the change when it is the
    effect of the instrument.

    NOT raised for a differing sampling seed. The seed is GIVEN to both sides
    through ``DELFIN_SAMPLING_SEED``, so the two are asked for the same sample;
    whether an endpoint honoured it is reported, not compared. The docstring
    used to list a seed mismatch as a refusal while no code checked one, which
    is a promise a reader would rely on.
    """


# ---------------------------------------------------------------------------
# Task selection: two routes, each one stated
# ---------------------------------------------------------------------------

#: Changed file -> the capability it touches. Ordered: the first matching rule
#: wins, so the specific patterns sit above the general ones. ``rule`` is the
#: capability name that appears in the report, and ``select`` is the strategy:
#: ``"suite"`` takes everything (the turn path is every task's path),
#: ``"skill:<name>"`` and ``"specialist"`` and ``"file:<path>"`` narrow.
_CHANGE_RULES: tuple[tuple[str, str, str], ...] = (
    # A packaged skill: the tasks that declare it are the ones that exercise it.
    (r"^delfin/agent/pack/skills/(?P<skill>[A-Za-z0-9_\-]+)\.md$",
     "skill", "skill"),
    # A specialist preset (or the registry that builds the presets): the tasks
    # whose prompt asks for a sub-agent are the ones that reach it.
    (r"^delfin/agent/pack/agents/", "specialist", "specialist"),
    (r"^delfin/agent/subagents\.py$", "specialist", "specialist"),
    # One benchmark suite file: exactly the tasks it defines.
    (r"^delfin/agent/pack/benchmark/(?P<file>[A-Za-z0-9_\-]+\.ya?ml)$",
     "task-pack", "file"),
    # The measurement code itself: every task is graded by it, so every task is
    # evidence about it.
    (r"^delfin/agent/bench", "bench", "suite"),
    (r"^delfin/agent/benchmark", "bench", "suite"),
    # The turn path: every task goes through it.
    (r"^delfin/agent/(api_client|engine|tool_catalog|prompt_loader)\.py$",
     "turn-path", "suite"),
    (r"^delfin/agent/cli\.py$", "turn-path", "suite"),
)


def _pack_dir() -> Path:
    """The directory the packaged benchmark suites live in.

    Read from the loader's own constant so it moves when the loader moves; the
    derivation is the fallback for a build that renames it. Two rules for one
    location would be a defect, so the loader's answer is used when present.
    """
    from . import benchmark as bm
    default = getattr(bm, "_DEFAULT_TASKS_PATH", None)
    if default is not None:
        return Path(default).parent
    return Path(bm.__file__).resolve().parent / "pack" / "benchmark"


def task_metadata(pack_dir: Path | str | None = None) -> dict[str, dict[str, str]]:
    """``{task_id: {"file": <yaml name>, "skill": <skill name>}}``.

    Read from the suite files, not from :class:`~.benchmark.Task`: ``skill`` is
    a key the task file carries for a human and a selection rule, and
    ``_coerce_task`` deliberately drops it (it changes nothing about how the
    task runs). Reading it here keeps that drop and still lets ``--changes``
    find the tasks a skill change affects.

    ``pack_dir`` defaults to the INSTALLED pack. A caller comparing two refs of
    a checkout can pass that checkout's own ``pack/benchmark`` so the selection
    describes the trees being measured rather than the one running the command.

    Unreadable or malformed files are skipped: a selection route that crashes
    on one bad suite file would take the whole verb down.
    """
    out: dict[str, dict[str, str]] = {}
    try:
        import yaml as _yaml
    except ImportError:                                        # pragma: no cover
        return out
    base = Path(pack_dir) if pack_dir is not None else _pack_dir()
    for path in sorted(base.glob("tasks*.yaml")):
        try:
            raw = _yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except (OSError, ValueError):
            continue
        items = raw.get("tasks") if isinstance(raw, dict) else raw
        if not isinstance(items, list):
            continue
        for item in items:
            if not isinstance(item, dict):
                continue
            tid = str(item.get("id") or "").strip()
            if not tid:
                continue
            out[tid] = {"file": path.name,
                        "skill": str(item.get("skill") or "").strip()}
    return out


def select_by_tags(tasks: Sequence[Any],
                   tags: Iterable[str]) -> tuple[list[Any], list[str]]:
    """Tasks whose ``behavior`` tag or ``task_class`` matches any of ``tags``.

    Both are capability words the suite already carries: ``behavior`` names the
    agent BEHAVIOUR a task measures (verify, plan, scout, ask, act -- see the
    ``tasks_auto_behavior.yaml`` family) and ``task_class`` names the work it
    belongs to (chemistry, office, generalist, ...). Matching is on the whole
    lowercased word, so ``plan`` does not select ``planning_notes``.

    Returns ``(selected, unmatched)``. The unmatched list is returned rather
    than dropped, for the reason ``benchmark.select_tasks`` gives: a typo that
    silently ran the rest reports a number for a question nobody asked.
    """
    wanted = {str(t).strip().lower() for t in tags if str(t).strip()}
    if not wanted:
        return [], []
    hit: set[str] = set()
    for task in tasks:
        words = {str(getattr(task, "behavior", "") or "").strip().lower(),
                 str(getattr(task, "task_class", "") or "").strip().lower()}
        matched = wanted & {w for w in words if w}
        if matched:
            hit.update(matched)
    selected = [
        t for t in tasks
        if ({str(getattr(t, "behavior", "") or "").strip().lower(),
             str(getattr(t, "task_class", "") or "").strip().lower()} & wanted)
    ]
    return selected, sorted(wanted - hit)


def select_by_changes(tasks: Sequence[Any],
                      changed: Iterable[str],
                      metadata: Optional[dict[str, dict[str, str]]] = None,
                      *,
                      pack_dir: Path | str | None = None,
                      ) -> tuple[list[tuple[Any, str]], list[str]]:
    """Tasks the changed files can be evidence about, each with its rule.

    ``changed`` is the list of files a branch changed, relative to the
    repository root. Each file is mapped through :data:`_CHANGE_RULES` to a
    capability and a selection strategy; the task then carries the capability
    and the file that selected it, so the report can say WHY a task ran.

    A changed file no rule covers selects the whole suite and says so: the
    conservative reading of "this change is somewhere under the agent" is that
    any task may be affected, and silently selecting nothing would report "no
    difference" over an empty population.

    Returns ``(selected, unmatched_files)``. ``unmatched_files`` is the
    fallback-triggered list -- named, because a rule table that quietly grows a
    catch-all is a rule table nobody has to maintain.
    """
    import re

    meta = (metadata if metadata is not None
            else task_metadata(pack_dir))
    by_id = {str(getattr(t, "id", "")): t for t in tasks}
    chosen: dict[str, str] = {}
    unmatched: list[str] = []
    for rel in changed:
        rel = str(rel or "").strip()
        if not rel:
            continue
        rule_name = ""
        for pattern, capability, strategy in _CHANGE_RULES:
            m = re.match(pattern, rel)
            if not m:
                continue
            rule_name = f"{capability} ({rel})"
            if strategy == "suite":
                for tid in by_id:
                    chosen.setdefault(tid, rule_name)
            elif strategy == "skill":
                skill = m.group("skill")
                for tid, info in meta.items():
                    if info.get("skill") == skill and tid in by_id:
                        chosen.setdefault(
                            tid, f"skill:{skill} ({rel})")
            elif strategy == "specialist":
                # A specialist change is reached by the tasks that DELEGATE:
                # the ones whose prompt names a sub-agent.
                for task in tasks:
                    prompt = str(getattr(task, "prompt", "") or "").lower()
                    if "subagent" in prompt or "sub-agent" in prompt:
                        chosen.setdefault(
                            str(task.id), f"specialist ({rel})")
            elif strategy == "file":
                name = m.group("file")
                for tid, info in meta.items():
                    if info.get("file") == name and tid in by_id:
                        chosen.setdefault(tid, f"task-pack:{name} ({rel})")
            break
        else:
            unmatched.append(rel)
            for tid in by_id:
                chosen.setdefault(tid, f"unmapped ({rel})")
    return [(by_id[tid], chosen[tid]) for tid in by_id if tid in chosen], unmatched


def changed_files(repo: Path | str, base: str, head: str) -> list[str]:
    """Files that differ between ``base`` and ``head``, via git.

    ``...`` (three dots) so the comparison is against the MERGE BASE: what the
    branch changed, not what the branch is behind by. A rebase onto a newer
    main otherwise reports every file main moved as if the branch had touched
    it.
    """
    from . import safe_git
    r = safe_git.run(["diff", "--name-only", f"{base}...{head}"], str(repo),
                     timeout=60)
    if r.returncode != 0:
        raise AbRefused(
            f"could not list the files between {base} and {head}: "
            f"{(r.stderr or '').strip()[:300]}")
    return [ln.strip() for ln in (r.stdout or "").splitlines() if ln.strip()]


# ---------------------------------------------------------------------------
# One side: its worktree and its attempts
# ---------------------------------------------------------------------------


@dataclass
class Side:
    """One side of the comparison: a ref, its commit, its worktree."""

    name: str
    ref: str
    commit: str
    path: Optional[Path] = None
    worktree: Any = None


def resolve_commit(repo: Path | str, ref: str) -> str:
    """The commit ``ref`` names, or :class:`AbRefused` naming the ref.

    Refused rather than defaulted: measuring a silent fallback (``HEAD`` when
    ``origin/main`` does not resolve) reports a difference against the wrong
    baseline, which is the failure this whole module exists to prevent.
    """
    from . import safe_git
    if not safe_git.is_repository(repo):
        # A copied tree has no refs and no remote, so "fetch it" is advice the
        # user cannot follow. Said separately, because the cause is different:
        # nothing is missing here, the comparison is simply not possible.
        raise AbRefused(
            f"{str(repo)!r} is not a git repository, so {ref!r} names no "
            f"commit. A comparison needs a checkout with both sides in it; a "
            f"copy of a tree has neither.")
    r = safe_git.run(["rev-parse", "--verify", "--quiet", f"{ref}^{{commit}}"],
                     str(repo), timeout=30)
    commit = (r.stdout or "").strip()
    if r.returncode != 0 or not commit:
        raise AbRefused(
            f"ref {ref!r} does not name a commit here. Fetch it, or name one "
            f"that exists; a comparison needs both sides.")
    return commit


def make_worktree(repo: Path | str, commit: str, *, parent: Path | str | None = None,
                  prefix: str = "ab"):
    """A fresh detached worktree at ``commit``. Returns the ``WorktreeInfo``.

    ``base_ref`` is the COMMIT, not the ref name: two sides of one change must
    be pinned to what was resolved, not to a name that could move under them
    between the two creations.
    """
    from . import worktree as _wt
    return _wt.enter_worktree(repo, branch_prefix=prefix, parent=parent,
                              base_ref=commit)


#: ``Task`` fields that annotate a task instead of asking it, and are
#: therefore outside its digest. Each side hashes its own task, so a field one
#: side's dataclass does not have would read as a changed question: every
#: comparison against a commit older than the field would be refused, and the
#: refusal would blame the task's text. Only fields that change nothing the
#: model is asked or judged on may be listed -- ``capability`` says which
#: product capability a task is meant to measure, which is bookkeeping about
#: the suite, not part of the question it puts to the model.
#:
#: The list is a deny-list on purpose. An allow-list would silently ignore a
#: field added later that DOES shape the question, and the refusal it guards
#: would stop firing exactly when it should.
#:
#: Defined above the worker source because the worker is standalone text that
#: cannot import this module; it is substituted into that text below, so there
#: is one list and no third place for the two copies of the rule to drift.
NOT_THE_QUESTION: tuple[str, ...] = ("capability",)

#: The program run in a side's worktree. Kept as source, not as a file inside
#: the tree: the BASE side predates this module, so a helper file committed here
#: would not exist there -- and a helper that exists on one side only is an
#: instrument that differs between the arms.
#:
#: ``-P`` and ``PYTHONPATH`` are what make it measure the right tree; the
#: comment on the runner below says why, and the runner refuses the result when
#: the import came from somewhere else.
_WORKER_SOURCE = r'''
import json, sys
NOT_THE_QUESTION = __NOT_THE_QUESTION__
task_id, model, provider, backend, repeats, out_path = sys.argv[1:7]
payload = {"error": "", "delfin_file": "", "passes": [], "unmeasured": 0,
           "n_samples": 0, "task_digest": "", "success": None, "model": ""}


def digest(task):
    """sha256 over the task as THIS tree's code defines it.

    Duplicated from bench_ab.task_digest on purpose. The worker may be
    running against the BASE side, where bench_ab.py does not exist yet;
    importing the parent's copy would digest the parent's dataclass and
    compare a question against itself. What must be compared is what each
    side's own code would ask, so each side hashes its own.

    Fields named in ``NOT_THE_QUESTION`` are left out for the same reason in
    both copies: they annotate a task instead of asking it, and a side that
    predates one would otherwise report every task as a changed question.
    """
    import dataclasses, hashlib
    try:
        body = dataclasses.asdict(task)
    except TypeError:
        body = {"id": str(getattr(task, "id", "")),
                "prompt": str(getattr(task, "prompt", "") or "")}
    for field_name in NOT_THE_QUESTION:
        body.pop(field_name, None)
    return hashlib.sha256(
        json.dumps(body, sort_keys=True, default=str).encode("utf-8")
    ).hexdigest()


try:
    import delfin
    payload["delfin_file"] = getattr(delfin, "__file__", "") or ""
    from delfin.agent import benchmark as bm
    from delfin.agent import benchmark_runner as br
    task = next((t for t in bm.load_tasks() if t.id == task_id), None)
    if task is None:
        raise RuntimeError(
            "task %r is not in the suite this tree carries" % task_id)
    payload["task_digest"] = digest(task)
    result = br.run_task(task, model=model, backend=backend or "api",
                         provider=provider, repeats=int(repeats))
    flags = list(result.per_run_success or [])
    payload["n_samples"] = int(result.n_samples or 0)
    payload["model"] = str(result.model or "")
    payload["success"] = bool(result.success)
    # `unmeasured` is the field that answers "did anything measure?", and it is
    # asked FIRST, because `per_run_success` cannot answer it on its own:
    # aggregate_replicates marks a result unmeasured only when EVERY replicate
    # was, and in that case it KEEPS them carrying success=False -- so
    # per_run_success is [False, False], not empty. Reading the flags alone
    # therefore scored an endpoint outage as two failures of the model: measured
    # 2026-11-25 with no provider key present, the verb printed "0/2 vs 0/2,
    # no measurable difference" and exited 0 for a comparison that called
    # nothing. An unmeasured attempt is not a failure of the model; it is an
    # attempt that never put the question.
    if getattr(result, "unmeasured", False) or not flags:
        payload["unmeasured"] = max(1, int(result.n_samples or 1))
        payload["passes"] = []
    else:
        payload["passes"] = [bool(x) for x in flags]
    payload["error"] = str(result.error or "")[:400]
except BaseException as exc:
    payload["error"] = ("%s: %s" % (type(exc).__name__, exc))[:400]
finally:
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh)
'''
# The worker is text, so the field list is carried into it literally: one
# definition, substituted where the worker is defined rather than kept in two.
_WORKER_SOURCE = _WORKER_SOURCE.replace(
    "__NOT_THE_QUESTION__", repr(tuple(NOT_THE_QUESTION)))


@dataclass
class RunSpec:
    """Everything one attempt needs, resolved once."""

    repo: Path
    side: Side
    repeat: int
    model: str
    provider: str = ""
    backend: str = "api"
    max_tokens: int = 0
    seed: str = ""
    temperature: str = ""
    timeout_s: float = 3600.0


def task_digest(task: Any) -> str:
    """sha256 over the task AS DELFIN EXECUTES IT.

    The fields that define the question -- prompt, mode, class, signals,
    expected figures, verify/setup scripts, budgets, permission mode, rubric --
    in a stable serialization. Two sides offering the same task id with
    different text are two questions, and a difference between their pass rates
    is not the effect of the change under test.

    ``skill:`` is deliberately outside the digest: it changes no field the
    runner reads, so a skill rename is not a new question. ``NOT_THE_QUESTION``
    names the fields left out for the same kind of reason.
    """
    try:
        payload = dataclasses.asdict(task)
    except TypeError:                                          # pragma: no cover
        payload = {"id": str(getattr(task, "id", "")),
                   "prompt": str(getattr(task, "prompt", "") or "")}
    for field_name in NOT_THE_QUESTION:
        payload.pop(field_name, None)
    blob = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


#: Provider -> the environment variable the credential loader reads for it.
#: Mirrors ``api_client.build_client``'s per-provider ``key_env_var`` and
#: ``credentials._WELL_KNOWN_KEYS``, and a test keeps the three in step. Only
#: the provider under test crosses into a side's worker: see ``_child_env``.
PROVIDER_KEY_NAMES: dict[str, str] = {
    "kit": "KIT_TOOLBOX_API_KEY",
    "blablador": "BLABLADOR_API_KEY",
    "openai": "OPENAI_API_KEY",
    "claude": "ANTHROPIC_API_KEY",
    "anthropic": "ANTHROPIC_API_KEY",
}

#: The provider whose key crosses when the caller names none. ``build_client``
#: with ``backend="api"`` and no provider resolves ``APIClient``, which reads
#: ``ANTHROPIC_API_KEY`` -- so this is that default, named rather than assumed.
DEFAULT_KEY_PROVIDER = "claude"


def _child_env(worktree: Path, home: Path, seed: str,
               temperature: str, *, parent: Optional[dict] = None,
               provider: str = "") -> dict[str, str]:
    """The environment a side's worker runs in.

    ``PYTHONPATH`` first on the import path and the SIDE'S OWN tree opening it:
    that is how the attempt measures the ref it belongs to instead of whatever
    ``delfin`` the parent process imported. ``HOME`` is a fresh directory per
    attempt, so nothing the agent would read about a previous attempt -- a
    memory store, an outcome history -- exists to be read.

    THAT FRESH HOME ALSO HIDES THE CREDENTIAL FILE, and a credential is what an
    attempt cannot do without. ``credentials.load_credential`` resolves env >
    keyring > file; with the file behind the redirected home the environment
    variable is the only route left across the process boundary, so the one key
    of the provider under test is carried over explicitly. Without it every
    attempt died in engine init and the verb reported a pass rate for a run
    that called nothing -- measured on this branch before the fix.

    ONLY THAT KEY. ``parent`` is read for the named provider's variable and no
    other: a benchmark measures an agent that can run shell commands, and a
    GitHub token in its environment is a token it can spend. The parent
    environment is passed in rather than read here so a test can state what it
    hands over.
    """
    parent = os.environ if parent is None else parent
    env = {
        "PATH": parent.get("PATH", "/usr/bin:/bin"),
        "HOME": str(home),
        "LANG": parent.get("LANG", "C.UTF-8"),
        "LC_ALL": parent.get("LC_ALL", "C.UTF-8"),
        "PYTHONPATH": str(worktree),
        "PYTHONDONTWRITEBYTECODE": "1",
        "DELFIN_NO_BROWSER": "1",
        # Keep the measurement's own state inside the fresh home, and out of
        # any directory that belongs to a person.
        "DELFIN_SCRATCH_STATE": str(home / ".delfin"),
    }
    key_name = PROVIDER_KEY_NAMES.get(
        (provider or DEFAULT_KEY_PROVIDER).strip().lower(), "")
    key_value = parent.get(key_name, "") if key_name else ""
    if key_value:
        env[key_name] = key_value
    if seed:
        env["DELFIN_SAMPLING_SEED"] = str(seed)
    if temperature:
        env["DELFIN_SAMPLING_TEMPERATURE"] = str(temperature)
    return env


def subprocess_runner(task: Any, spec: RunSpec) -> SideResult:
    """Run one side's attempts on one task, in that side's worktree.

    A SUBPROCESS, not an in-process call, and that is the whole design:

    * the code under test must be the SIDE'S code -- a CWD change does not
      rebind an already-imported ``delfin``, so an in-process run would measure
      the parent's tree twice and report zero difference by construction;
    * the fresh ``HOME`` is then a real process boundary rather than a patched
      module attribute;
    * a side that hangs or aborts costs one attempt, not the comparison.

    ``-P`` matters: a benchmark workspace is itself a copy of the repository, so
    plain ``python -m`` would put the WORKSPACE first on ``sys.path`` and import
    the agent under test as the agent doing the testing (observed: six runs of
    one member refused by exactly that check).
    """
    if spec.side.path is None:
        raise AbRefused(f"side {spec.side.name} has no worktree")
    wt = Path(spec.side.path)
    import tempfile
    with tempfile.TemporaryDirectory(prefix=f"ab-{spec.side.name}-home-") as h:
        home = Path(h)
        out = home / "attempt.json"
        cmd = [sys.executable, "-P", "-c", _WORKER_SOURCE,
               str(getattr(task, "id", "")), spec.model, spec.provider,
               spec.backend, str(max(1, int(spec.repeat))), str(out)]
        r = subprocess.run(cmd, cwd=str(wt), env=_child_env(
            wt, home, spec.seed, spec.temperature, provider=spec.provider),
            capture_output=True, text=True, timeout=spec.timeout_s)
        payload: dict[str, Any] = {}
        try:
            payload = json.loads(out.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            payload = {"error": (
                "the attempt worker wrote no result (exit %s): %s"
                % (r.returncode, ((r.stderr or r.stdout or "")
                                  .strip()[-400:]))[:400])}
    imported = str(payload.get("delfin_file") or "")
    if imported and not Path(imported).resolve().is_relative_to(wt.resolve()):
        # The measurement-of-the-wrong-tree defect, caught where it happens
        # rather than after the numbers look plausible.
        raise AbRefused(
            f"side {spec.side.name} imported delfin from {imported}, which is "
            f"not in its worktree {wt} -- the two sides would not be two.")
    passes = [bool(x) for x in (payload.get("passes") or [])]
    return SideResult(
        passes=passes,
        unmeasured=int(payload.get("unmeasured") or 0),
        errors=[str(payload.get("error") or "")][:1] if payload.get("error")
        else [],
        model=str(payload.get("model") or spec.model),
        task_digest=str(payload.get("task_digest") or task_digest(task)),
        module_file=imported,
        worktree=str(wt),
    )


# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------


@dataclass
class ABReport:
    """The measurement, in the shape a caller prints and a test asserts."""

    base: Side
    head: Side
    tasks: list[TaskAB] = field(default_factory=list)
    repeat: int = 1
    model: str = ""
    seed: str = ""
    temperature: str = ""
    budget_limit: float = 0.0
    budget_spent: float = 0.0
    budget_refused: list[str] = field(default_factory=list)
    unmatched_tags: list[str] = field(default_factory=list)
    unmatched_files: list[str] = field(default_factory=list)
    selection: str = ""

    @property
    def n_no_difference(self) -> int:
        return sum(1 for t in self.tasks
                   if t.verdict == VERDICT_NO_DIFFERENCE)

    @property
    def n_improved(self) -> int:
        return sum(1 for t in self.tasks if t.verdict == VERDICT_IMPROVED)

    @property
    def n_regressed(self) -> int:
        return sum(1 for t in self.tasks if t.verdict == VERDICT_REGRESSED)

    @property
    def n_unmeasured(self) -> int:
        return sum(1 for t in self.tasks if t.verdict == VERDICT_UNMEASURED)

    def as_dict(self) -> dict[str, Any]:
        return {
            "base": {"ref": self.base.ref, "commit": self.base.commit},
            "head": {"ref": self.head.ref, "commit": self.head.commit},
            "repeat": self.repeat,
            "model": self.model,
            "seed": self.seed,
            "temperature": self.temperature,
            "same_commit": self.base.commit == self.head.commit,
            "interval_method": INTERVAL_METHOD,
            "interval_why": INTERVAL_WHY,
            "budget": {"limit": self.budget_limit, "spent": self.budget_spent,
                       "refused": list(self.budget_refused)},
            "selection": self.selection,
            "unmatched_tags": list(self.unmatched_tags),
            "unmatched_files": list(self.unmatched_files),
            "tasks": [t.as_dict() for t in self.tasks],
            "n_no_difference": self.n_no_difference,
            "n_improved": self.n_improved,
            "n_regressed": self.n_regressed,
            "n_unmeasured": self.n_unmeasured,
        }


def _fmt_rate(side: SideResult) -> str:
    if side.rate is None:
        return "—"
    return f"{side.n_pass}/{side.measured}"


def format_report(rep: ABReport) -> str:
    """The report as lines. States the method before it states any number.

    A reader who sees an interval before seeing how it was built has no way to
    check it; the header carries the method, the reason, and the repeat count,
    and only then the table.
    """
    d = rep.as_dict()
    lines: list[str] = []
    lines.append(f"bench ab: {rep.base.ref} ({rep.base.commit[:12]}) "
                 f"→ {rep.head.ref} ({rep.head.commit[:12]})")
    if d["same_commit"]:
        lines.append("  NOTE: both sides are the SAME COMMIT — this measures "
                     "sampling noise alone, and any difference reported here "
                     "is a false positive by construction.")
    lines.append(f"  model {rep.model or '<unset>'} · repeat {rep.repeat} per "
                 f"task per side · seed {rep.seed or '<provider default>'}"
                 f"{' · temperature ' + rep.temperature if rep.temperature else ''}")
    lines.append(f"  interval: {INTERVAL_METHOD}")
    lines.append(f"  why: {INTERVAL_WHY}")
    lines.append(f"  selection: {rep.selection or '<all tasks>'}")
    if rep.unmatched_tags:
        lines.append(f"  ⚠️  tags that matched no task: "
                     f"{', '.join(rep.unmatched_tags)}")
    if rep.unmatched_files:
        lines.append(f"  ⚠️  changed files no capability rule covers "
                     f"(whole suite selected): "
                     f"{', '.join(rep.unmatched_files[:8])}"
                     f"{' …' if len(rep.unmatched_files) > 8 else ''}")
    if rep.budget_limit:
        lines.append(f"  budget: {rep.budget_spent:.2f} of "
                     f"{rep.budget_limit:.2f} core-hours"
                     + (f", {len(rep.budget_refused)} attempt(s) refused"
                        if rep.budget_refused else ""))
    lines.append("")
    lines.append(f"  {'task':<38} {'base':>7} {'head':>7} {'Δ rate':>7}  "
                 f"{'interval':>18}  verdict")
    for t in rep.tasks:
        row = t.as_dict()
        lo, hi = row["interval_low"], row["interval_high"]
        delta = "—" if row["diff"] is None else f"{row['diff']:+.2f}"
        verdict = row["verdict"]
        if verdict == VERDICT_NO_DIFFERENCE:
            verdict = NO_DIFFERENCE
        elif verdict == VERDICT_UNMEASURED:
            verdict = (f"not measured (base {row['unmeasured_base']} / head "
                       f"{row['unmeasured_head']} attempt(s) reached nothing)")
        lines.append(f"  {t.task_id:<38} {_fmt_rate(t.base):>7} "
                     f"{_fmt_rate(t.head):>7} {delta:>7}  "
                     f"{f'[{lo:+.2f}, {hi:+.2f}]':>18}  {verdict}")
        if t.selection_rule:
            lines.append(f"       selected by: {t.selection_rule}")
        for err in row["errors"][:2]:
            lines.append(f"       error: {err[:160]}")
    for reason in rep.budget_refused[:6]:
        lines.append(f"  refused by budget: {reason}")
    lines.append("")
    lines.append(f"  {len(rep.tasks)} task(s): {rep.n_improved} improved, "
                 f"{rep.n_regressed} regressed, {rep.n_no_difference} "
                 f"{NO_DIFFERENCE}, {rep.n_unmeasured} not measured")
    lines.append(f"  A task is only called improved or regressed when the "
                 f"interval on its difference EXCLUDES zero; at repeat "
                 f"{rep.repeat} most real differences will not clear that bar, "
                 f"and that is the honest answer rather than a miss.")
    return "\n".join(lines)


def _check_pair(rep_tasks: Sequence[TaskAB], model: str,
                seed: str) -> None:
    """Refuse a comparison whose sides are not one instrument.

    Read from what the ATTEMPTS reported, not from the arguments the caller
    passed: a caller can pass one model and a side can still resolve another,
    and the run is the record of what happened.
    """
    models = {m for t in rep_tasks for m in
              (t.base.model, t.head.model) if m}
    if len(models) > 1:
        raise AbRefused(
            "the two sides did not use the same model ("
            + ", ".join(sorted(models)) + "). A difference between them is a "
            "difference between models, not between the code under test.")
    if model and models and model not in models:
        raise AbRefused(
            f"the sides used {', '.join(sorted(models))}, not the requested "
            f"{model!r}; the comparison would be reported under the wrong "
            f"instrument.")
    differing = [t.task_id for t in rep_tasks
                 if t.base.task_digest and t.head.task_digest
                 and t.base.task_digest != t.head.task_digest]
    if differing:
        raise AbRefused(
            "the same task id is a DIFFERENT question on the two sides: "
            + ", ".join(differing[:6])
            + (f" (+{len(differing) - 6} more)" if len(differing) > 6 else "")
            + ". A task whose text changed between the refs is a new "
            "measurement, not a comparable one.")
    del seed  # carried in the report; the sides share the parent's value


def run(*,
        base: str = "origin/main",
        head: str = "HEAD",
        model: str = "",
        provider: str = "",
        backend: str = "api",
        repeat: int = 3,
        task_ids: Sequence[str] = (),
        tags: Sequence[str] = (),
        changes: str = "",
        changed: Optional[Sequence[str]] = None,
        tasks: Optional[Sequence[Any]] = None,
        repo: Path | str | None = None,
        budget: float | None = None,
        accounting: Any = None,
        runner: Optional[Callable[[Any, RunSpec], SideResult]] = None,
        worktree_factory: Optional[Callable[..., Any]] = None,
        worktree_parent: Path | str | None = None,
        seed: str = "",
        temperature: str = "",
        log: Callable[[str], None] = print,
        ) -> ABReport:
    """Run the paired comparison and return the measurement.

    ``runner`` and ``worktree_factory`` are injectable so a test measures the
    ORCHESTRATION (the gate, the refusal, the arithmetic, the cleanup) without
    a model and without two real checkouts. The defaults are the real ones, and
    a test that only ever passes fakes proves nothing about them -- which is
    why ``tests/test_bench_ab.py`` also drives the real worktree factory once.

    ``budget`` defaults to the nightly round's own daily limit
    (``bench_nightly.DEFAULT_DAILY_BUDGET_CORE_HOURS``), and the gate is
    applied per ATTEMPT with the same rule the nightly round uses: the estimate
    joins the day's spend before the attempt starts, and an attempt that would
    push the total past the limit is refused and named.
    """
    repo_path = Path(repo or os.getcwd())
    if repeat < 1:
        raise AbRefused(f"--repeat must be at least 1 (got {repeat})")
    if not model:
        raise AbRefused(
            "a comparison needs --model: without one each side may resolve its "
            "own default, and two models measured against each other is not an "
            "A/B of the code.")

    # Seed: ask both sides for the same sampling. Read from the parent's own
    # configuration when the caller names none, so the number in the report is
    # the one the attempts were given rather than a guess.
    from . import sampling as _sampling
    cfg = _sampling.configured()
    seed = str(seed or cfg.get("seed") or "")
    temperature = str(temperature or cfg.get("temperature") or "")

    if accounting is None:
        from .bench_holdout_glue import GlueAccounting
        accounting = GlueAccounting()
    if budget is None:
        from .bench_nightly import DEFAULT_DAILY_BUDGET_CORE_HOURS
        budget = DEFAULT_DAILY_BUDGET_CORE_HOURS
    budget = float(budget)
    runner = runner or subprocess_runner
    factory = worktree_factory or make_worktree

    # --- task selection ---------------------------------------------------
    from . import benchmark as bm
    # ``tasks`` is the caller's population, when it has one. It exists so a
    # test measures the ORCHESTRATION against a fixed set instead of the
    # installed pack (whose contents change between waves), and so a caller can
    # pin a population explicitly. Default: the packaged suite.
    all_tasks = list(tasks) if tasks is not None else bm.load_tasks()
    if not all_tasks:
        raise AbRefused("no benchmark tasks were found; nothing to compare")
    by_id = {str(t.id): t for t in all_tasks}

    selection_parts: list[str] = []
    unmatched_tags: list[str] = []
    unmatched_files: list[str] = []
    rules: dict[str, str] = {}
    chosen: list[Any] = []

    if task_ids:
        picked, missing = bm.select_tasks(all_tasks, list(task_ids))
        chosen = picked
        rules = {str(t.id): "task id / glob" for t in picked}
        selection_parts.append(f"tasks {', '.join(task_ids)}")
        if missing:
            unmatched_tags.extend(f"task:{m}" for m in missing)
    if tags:
        picked, missing = select_by_tags(all_tasks, tags)
        for t in picked:
            if t not in chosen:
                chosen.append(t)
        for t in picked:
            rules.setdefault(str(t.id), f"tag ({', '.join(tags)})")
        selection_parts.append(f"tags {', '.join(tags)}")
        unmatched_tags.extend(missing)
    if changes or changed is not None:
        files = (list(changed) if changed is not None
                 else changed_files(repo_path, base, changes or head))
        picked, unmatched = select_by_changes(
            all_tasks, files,
            pack_dir=(repo_path / "delfin" / "agent" / "pack" / "benchmark"))
        for task, rule in picked:
            if task not in chosen:
                chosen.append(task)
            rules.setdefault(str(task.id), rule)
        selection_parts.append(
            f"{len(files)} changed file(s) between {base} and "
            f"{changes or head}")
        unmatched_files.extend(unmatched)
    if not (task_ids or tags or changes or changed is not None):
        chosen = list(all_tasks)
        selection_parts.append("all tasks")
        rules = {str(t.id): "all tasks" for t in chosen}

    if not chosen:
        raise AbRefused(
            "the selection matched no task"
            + (f" (unmatched: {', '.join(unmatched_tags)})"
               if unmatched_tags else "")
            + ". A comparison over an empty population reports a difference "
              "of zero and means nothing.")

    # --- the two sides ----------------------------------------------------
    base_commit = resolve_commit(repo_path, base)
    head_commit = resolve_commit(repo_path, head)
    sides = [Side("base", base, base_commit), Side("head", head, head_commit)]
    log(f"bench ab: {base} ({base_commit[:12]}) → {head} ({head_commit[:12]})"
        f" · {len(chosen)} task(s) × {repeat} repeat(s) × 2 sides")

    report = ABReport(base=sides[0], head=sides[1], repeat=repeat, model=model,
                      seed=seed, temperature=temperature,
                      budget_limit=budget,
                      unmatched_tags=sorted(set(unmatched_tags)),
                      unmatched_files=sorted(set(unmatched_files)),
                      selection="; ".join(selection_parts))

    opened: list[Any] = []
    try:
        for side in sides:
            info = factory(repo_path, side.commit, parent=worktree_parent)
            side.worktree = info
            side.path = Path(getattr(info, "path", info))
            opened.append(info)

        # --- the attempts, under the nightly round's budget rule ----------
        spent = 0.0
        pairs = {str(t.id): {} for t in chosen}
        for side in sides:
            for task in chosen:
                tid = str(task.id)
                est = float(accounting.estimate(task))
                if est < 0:
                    # A negative estimate would reduce the day's spend and let
                    # later work fit for free; the nightly round refuses it for
                    # the same reason.
                    report.budget_refused.append(
                        f"{side.name}/{tid}: negative estimate {est:.3f}")
                    pairs[tid][side.name] = SideResult(
                        errors=[f"negative resource estimate {est:.3f}"])
                    continue
                if spent + est > budget:
                    report.budget_refused.append(
                        f"{side.name}/{tid}: estimate {est:.3f} + spent "
                        f"{spent:.3f} exceeds the {budget:.3f} core-hour limit")
                    pairs[tid][side.name] = SideResult(
                        errors=[f"over the daily budget "
                                f"({budget:.2f} core-hours)"])
                    continue
                accounting.charge(est)
                spent += est
                spec = RunSpec(repo=repo_path, side=side, repeat=repeat,
                               model=model, provider=provider, backend=backend,
                               seed=seed, temperature=temperature)
                try:
                    pairs[tid][side.name] = runner(task, spec)
                except AbRefused:
                    raise
                except Exception as exc:                       # noqa: BLE001
                    pairs[tid][side.name] = SideResult(
                        errors=[f"{type(exc).__name__}: {exc}"[:300]])
        report.budget_spent = spent

        for task in chosen:
            tid = str(task.id)
            got = pairs.get(tid) or {}
            report.tasks.append(TaskAB(
                task_id=tid,
                base=got.get("base") or SideResult(),
                head=got.get("head") or SideResult(),
                selection_rule=rules.get(tid, "")))

        _check_pair(report.tasks, model, seed)
    finally:
        from . import worktree as _wt
        # Never remove the repository itself or a directory that CONTAINS it.
        # `run(worktree_factory=...)` is injectable, and a factory answering
        # with the user's checkout would otherwise have its answer deleted --
        # the failure mode that cost this package's own working tree once (see
        # tests/test_bench_ab.py::_factory). The identity of a real side is
        # checked separately, below, by the attempt that imported from it.
        repo_resolved = repo_path.resolve()
        for info in opened:
            path = Path(getattr(info, "path", info)).resolve()
            if path == repo_resolved or repo_resolved.is_relative_to(path):
                continue
            try:
                _wt.exit_worktree(info, keep_if_changed=False)
            except Exception:                                  # noqa: BLE001
                pass
    return report


def _cmd_options_from_args(args: Any) -> dict[str, Any]:
    """The CLI's ``args`` as :func:`run`'s keywords. One place, so the parser and
    the call cannot drift apart."""
    def _list(name: str) -> list[str]:
        raw = getattr(args, name, None) or []
        if isinstance(raw, str):
            raw = [raw]
        out: list[str] = []
        for item in raw:
            out.extend(p.strip() for p in str(item).split(",") if p.strip())
        return out

    return {
        "base": getattr(args, "base", "") or "origin/main",
        "head": getattr(args, "head", "") or "HEAD",
        "model": (getattr(args, "model", "") or "").strip(),
        "provider": getattr(args, "provider", "") or "",
        "backend": getattr(args, "backend", "") or "api",
        "repeat": int(getattr(args, "repeat", 3) or 3),
        "task_ids": _list("ab_task"),
        "tags": _list("tags"),
        "changes": getattr(args, "changes", "") or "",
        "budget": getattr(args, "budget", None),
        "seed": getattr(args, "seed", "") or "",
        "temperature": getattr(args, "temperature", "") or "",
    }
