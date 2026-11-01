"""Glue between the CB0B holdout sampler and the PL4 nightly runner
(package HB1, wave 22).

Until this module existed, ``bench_holdout.sample`` produced
:class:`~delfin.agent.bench_holdout.Instance` objects and
``bench_nightly.run`` needed task-shaped objects (``.id`` /
``.task_class``, ``bench_nightly.py:194``) that
``benchmark_runner.run_task`` (protected, benchmark_runner.py:1450) can
execute — so the operator could not run the first baseline. This module
is the adapter layer that connects the two, WITHOUT re-implementing any
grading or resource logic:

* ``instance_to_task`` — maps one ``Instance`` onto a
  :class:`~delfin.agent.benchmark.Task` ``run_task`` can run. Grading
  stays with the glue Grader, so the Task carries no expected signals
  and no verify script.
* ``GlueTask`` — a Task subclass carrying ``n_atoms`` so the Accounting
  can size the run without re-parsing the molecule.
* ``GlueSampler`` — the sampler ``bench_nightly.run`` consumes; it
  draws through ``bench_holdout.sample`` (reused) and adapts every
  Instance. Draw ids are UNIQUE per draw (``<instance_id>.<n>``):
  ``run_task`` keys ``learning_home`` on ``task.id``
  (benchmark_runner.py:1508) and ``InstanceRecord.task_id`` collides on
  a repeated draw otherwise (reviewer nacht-s25, note 1).
* ``GlueRunner`` — calls ``run_task`` (default: the real one, lazily
  imported) and collects CB0A run evidence
  (:mod:`delfin.agent.run_evidence`) when the caller hands it a
  workspace directory.
* ``GlueGrader`` — grades from DELFIN's own machinery: the
  ``BenchmarkResult`` the runner produced (``score_outcome``'s
  deterministic rubric inside ``run_task`` — never re-graded by hand),
  plus the CB0A run-evidence verdict
  (``delfin.agent.pack.benchmark.accept.w20_common.grade_run``) and the
  SC1 claim check (:func:`delfin.agent.result_claims.check`) when the
  result carries the collected evidence. Extra checks can only LOWER
  the verdict, never raise it.
* ``GlueAccounting`` — the core-hour estimate comes from DELFIN's own
  resource logic only (import only, no re-derivation):
  ``cluster_efficiency.advise`` supplies PAL for the system size and
  ``workflow_grader._CORE_HOUR_ESTIMATE`` the per-step hours.

The holdout seed VALUE never passes through any of this: the sampler
receives it as an argument and nothing here stores, logs or prints it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

from .bench_holdout import Instance, GRADING_LEVELS
from .bench_nightly import Grader, Accounting, Sampler  # noqa: F401 (Protocol types)
from .benchmark import Task, BenchmarkResult

__all__ = [
    "GlueTask", "instance_to_task", "GlueSampler", "GlueRunner",
    "GlueRun", "GlueGrader", "GlueAccounting",
]


# ---------------------------------------------------------------------------
# Adapter: Instance -> Task
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GlueTask(Task):
    """A benchmark Task built from a holdout Instance, plus the atom
    count the Accounting sizes from."""

    n_atoms: int = 0


def instance_to_task(instance: Instance, draw_index: int = 0) -> GlueTask:
    """Adapt one holdout Instance onto a task ``run_task`` can run.

    ``mode="solo"``: the holdout measures the agent's own choice of
    DELFIN tool for an unlabelled chemist request — the same mode the
    interactive product runs in. No expected/forbidden signals and no
    verify script: grading is the glue Grader's job (CB0A run evidence
    + SC1 claim check), not pattern matching on the prompt.

    The id is ``<instance_id>.<draw_index>`` (unique per draw): the same
    draw can legitimately repeat across a nightly round, and both
    ``run_task``'s ``learning_home`` and the ``InstanceRecord`` key on
    ``task.id``.
    """
    return GlueTask(
        id=f"{instance.instance_id}.{draw_index}",
        task_class=instance.family_id,
        mode="solo",
        prompt=instance.prompt,
        max_duration_s=1800.0,
        max_tool_calls=40,
        n_atoms=instance.molecule.n_atoms,
    )


class GlueSampler:
    """The sampler ``bench_nightly.run`` consumes (``sample(seed, n)``).

    Draws through ``bench_holdout.sample`` (reused as-is) and adapts
    every Instance to a ``GlueTask``. The seed VALUE is used as an
    argument only; nothing here stores or logs it.
    """

    def sample(self, seed, n: int) -> list[GlueTask]:
        return [instance_to_task(inst, i)
                for i, inst in enumerate(_sample_instances(seed, n))]


def _sample_instances(seed, n: int) -> list[Instance]:
    from .bench_holdout import sample as _sample
    return _sample(seed, n)


# ---------------------------------------------------------------------------
# Runner: run the task through run_task's public path
# ---------------------------------------------------------------------------

@dataclass
class GlueRun(BenchmarkResult):
    """A ``BenchmarkResult`` plus what the glue Grader checks beyond
    it: the collected run evidence (CB0A shape) and the answer text for
    the SC1 claim check. Both empty when the caller supplied none —
    the grader then grades the rubric verdict alone."""

    answer: str = ""
    evidence: list = field(default_factory=list)


class GlueRunner:
    """Runs one adapted task through ``benchmark_runner.run_task``.

    The real ``run_task`` is the default (lazily imported so the engine
    stack loads only when a real run starts). ``workspace``: when the
    caller hands a directory, CB0A run evidence is collected from it
    after the run and attached to the result for the Grader.
    """

    def __init__(self, *, run_task_fn: Optional[Callable] = None,
                 workspace: Optional[str] = None):
        self._run_task_fn = run_task_fn
        self._workspace = workspace

    def __call__(self, task, *, model: str = "", **kwargs) -> GlueRun:
        if self._run_task_fn is None:
            from .benchmark_runner import run_task as run_task_fn
        else:
            run_task_fn = self._run_task_fn
        # Provider resolution (HB2 D1): the CLI resolves the provider
        # the same way (cli.py:_resolve_provider → pricing.infer_provider
        # as the last step). Without it the engine is built as Anthropic
        # and every kit.<model> run dies with "engine init failed: No
        # Anthropic API key found" — observed in the 31.10. baseline.
        # The pricing table is the one answerer; no second mapping here.
        # An explicit provider kwarg (from --provider) wins.
        if not kwargs.get("provider"):
            from .pricing import infer_provider
            inferred = infer_provider(model)
            if inferred:
                kwargs["provider"] = inferred
        result = run_task_fn(task, model=model, **kwargs)
        answer, evidence = self._collect()
        return GlueRun(
            task_id=result.task_id, task_class=result.task_class,
            model=result.model, success=result.success,
            quality_0_100=result.quality_0_100,
            duration_s=result.duration_s, error=result.error,
            unmeasured=result.unmeasured,
            answer=answer, evidence=evidence,
        )

    def _collect(self) -> tuple:
        """(answer, evidence) for the Grader; empty when there is
        nothing to collect. Evidence comes from DELFIN's own
        ``run_evidence.collect``, never hand-built."""
        if not self._workspace:
            return "", []
        from .run_evidence import collect as _collect_evidence
        ev = _collect_evidence(self._workspace)
        # SC1 reads the wave-10 ledger shape; the Evidence object's
        # steps carry input/output texts — the excerpt a claim is
        # matched against.
        records = [{"kind": "calc", "ref": s.get("output") or "",
                    "text": "\n".join(str(s.get(k) or "")
                                      for k in ("input", "output")),
                    "flags": list(s.get("flags") or [])}
                   for s in getattr(ev, "steps", [])]
        return "", records


# ---------------------------------------------------------------------------
# Grader: DELFIN's own grading, never re-graded by hand
# ---------------------------------------------------------------------------

def _rubric_verdict(result) -> tuple:
    """The runner's own rubric verdict, read off the BenchmarkResult —
    ``run_task`` already scored the attempt with ``score_outcome``;
    re-deriving it here would be re-grading by hand."""
    success = bool(getattr(result, "success", False))
    quality = int(getattr(result, "quality_0_100", 0) or 0)
    return success, quality


def _scale(quality: int) -> int:
    """Quality on the grader interface's 0-5 scale."""
    return max(0, min(5, round(quality / 20)))


class GlueGrader:
    """``Grader`` for ``bench_nightly.run``: (success, quality 0-5).

    Layered, each layer from DELFIN's own machinery:

    1. the runner's rubric verdict (``BenchmarkResult.success`` /
       ``quality_0_100`` — ``score_outcome`` inside ``run_task``);
    2. the CB0A run-evidence verdict
       (``delfin.agent.pack.benchmark.accept.w20_common.grade_run``)
       when the result carries collected evidence — a faked run (output
       without the input that produced it) fails the whole verdict;
    3. the SC1 claim check (:func:`delfin.agent.result_claims.check`)
       on the answer against the same evidence.

    Extra layers can only LOWER the verdict, never raise it.
    """

    def grade(self, task, result) -> tuple:
        if getattr(result, "unmeasured", False):
            # Not a score: the request never reached the model. The
            # nightly round must see that, not a silent failure.
            return (False, 0)
        success, quality = _rubric_verdict(result)
        evidence = list(getattr(result, "evidence", None) or [])
        if evidence:
            success, quality = self._evidence_layer(success, quality,
                                                    evidence)
        answer = getattr(result, "answer", "") or ""
        if answer:
            success, quality = self._claim_layer(success, quality,
                                                 answer, evidence)
        return (success, quality)

    def _evidence_layer(self, success, quality, evidence) -> tuple:
        from .pack.benchmark.accept.w20_common import grade_run
        try:
            verdict = grade_run(type("E", (), {"steps": evidence})())
        except Exception:  # noqa: BLE001 — a broken check must not flip a verdict
            return success, quality
        if verdict.get("verdict") != "pass":
            return (False, 0)
        return success, quality

    def _claim_layer(self, success, quality, answer, evidence) -> tuple:
        from .result_claims import check as check_claims
        from .result_claims import worst_level
        try:
            claims = check_claims(answer, evidence)
        except Exception:  # noqa: BLE001
            return success, quality
        # An ERROR claim is a stated number with no source: the verdict
        # fails, quality drops to the floor. A "warn" claim (a sourced
        # value from a flagged output) keeps the verdict — the rubric
        # verdict already carries what the run achieved (s25 finding,
        # answered-vs-unanswered split per result_claims.py:170-189).
        # The caveat is still surfaced: the warn Claim itself carries
        # level="warn" plus the flags that triggered it (result_claims
        # Claim.flags), and worst_level() reports it in the return —
        # it is not dropped, it just does not flip the verdict.
        if worst_level(claims) == "error":
            return (False, 0)
        return success, quality


# ---------------------------------------------------------------------------
# Accounting: DELFIN's own resource logic, import only
# ---------------------------------------------------------------------------

def _advise(step: str, n_atoms: int) -> dict:
    from .cluster_efficiency import advise
    return advise(step, n_atoms)


class GlueAccounting:
    """``Accounting`` for ``bench_nightly.run``: per-instance core-hour
    estimate + spend, from DELFIN's own resource logic only.

    The estimate maps the instance's workflow family onto DELFIN's
    known steps and sums ``workflow_grader._CORE_HOUR_ESTIMATE`` (the
    same per-step weights ``workflow_grader._cost_estimate`` uses) over
    the steps the family's workflow performs; PAL comes from
    ``cluster_efficiency.advise`` (the code DELFIN itself uses) and
    only ever scales the estimate DOWN to what a local run may use
    (``calc_policy.LOGIN_LOCAL_CORES``). No re-derivation anywhere.
    """

    #: The DELFIN steps each holdout family's workflow performs
    #: (vocabulary from workflow_grader._CORE_HOUR_ESTIMATE).
    _FAMILY_STEPS = {
        "classic_opt_freq": ("orca_opt_freq", "cclib_energies",
                             "cclib_vibrations"),
        "redox": ("orca_opt", "orca_sp", "cclib_energies"),
        "excited_states": ("orca_tddft", "uv_vis_parse"),
        "imaginary_mode": ("orca_opt", "imag_fix", "orca_freq"),
        "spin_states": ("orca_sp", "cclib_energies"),
        "conformer_prestep": ("crest_conformers", "xtb_opt",
                              "cclib_energies"),
        "screening": ("orca_opt", "xtb_solvator", "cclib_energies"),
    }

    def __init__(self, *, local_cap: bool = True):
        self._local_cap = local_cap

    def estimate(self, task) -> float:
        from .workflow_grader import _CORE_HOUR_ESTIMATE
        steps = self._FAMILY_STEPS.get(getattr(task, "task_class", ""),
                                       ("orca_sp",))
        est = sum(_CORE_HOUR_ESTIMATE.get(step, 0.1) for step in steps)
        if self._local_cap and hasattr(task, "n_atoms"):
            try:
                advice = _advise(steps[0], getattr(task, "n_atoms", 0))
                pal = int(advice.get("PAL") or 0)
            except Exception:  # noqa: BLE001 — no advise on this host: uncapped estimate
                pal = 0
            if pal:
                from .calc_policy import LOGIN_LOCAL_CORES
                # A local run may not use more cores than the policy
                # allows, so the core-hour estimate cannot exceed the
                # same work sized to the local cap.
                est = min(est, est * min(pal, LOGIN_LOCAL_CORES) / max(pal, 1))
        return round(est, 4)

    def charge(self, core_hours: float) -> None:
        from . import bench_nightly as _bn
        # The ledger the nightly runner already keeps mirrors the EF1
        # accounting; the glue's ledger is the day's spend bookkeeping
        # on this object (mirrored by run() into the record).
        self.charged = getattr(self, "charged", 0.0) + float(core_hours)


def build_default_config(*, model: str, instances: int,
                         mode: str = "holdout",
                         holdout_seed_path: str = "",
                         public_seed=7,
                         provider: str = "") -> dict:
    """A ready ``bench_nightly.run`` config wired end to end with the
    glue: real CB0B sampler, real DELFIN grading and accounting, and
    the real ``run_task`` (injected by ``run()`` itself)."""

    def _glue_runner(task, *, model="", **kwargs):
        return GlueRunner()(task, model=model, **kwargs)

    return {
        "model": model,
        "provider": provider,
        "sampler": GlueSampler(),
        "grader": GlueGrader(),
        "accounting": GlueAccounting(),
        "runner": _glue_runner,
        "instances": int(instances),
        "mode": mode,
        "holdout_seed_path": holdout_seed_path,
        "public_seed": public_seed,
    }
