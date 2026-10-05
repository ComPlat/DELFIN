"""Model-free optimization campaigns.

A :class:`Campaign` is a closed optimization loop over a bounded
candidate space (a parameter grid, e.g. substituent choices on a
molecular anchor).  It makes every decision itself — which point to
measure next, whether to stop, what to record — and wakes the agent
only on failure, stagnation or budget exhaustion.  All chemistry runs
through DELFIN's existing pieces:

* :func:`delfin.tools._runner.run_step` runs the ``xtb_sp`` adapter
  (native xTB, parses the HOMO-LUMO gap itself);
* the scheduler integration wraps candidates as
  :class:`~delfin.workflows.engine.classic.WorkflowJob` via
  :func:`delfin.tools._runner.step_as_workflow_job`;
* job status uses ``delfin.agent.job_monitor.register_agent_job`` /
  ``delfin.agent.cli_jobs.collect_job_rows`` — no private squeue calls;
* the Bayes-vs-random comparison reuses
  ``delfin.agent.benchmark.wilson_interval`` /
  ``_fisher_exact_2x2_pvalue``.

This module is intentionally library-minimal: acquisition is Expected
Improvement on a ``sklearn`` Gaussian process, seeding is
``scipy.stats.qmc.LatinHypercube``.  No new dependency.

Observations are a plain ``{candidate_index: gap_eV}`` mapping so a
campaign can be resumed from its files alone.
"""

from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

import numpy as np

__all__ = [
    "CampaignBudget",
    "CampaignDecision",
    "DecisionLog",
    "DecisionReason",
    "TargetGap",
    "campaign_should_wake",
    "expected_improvement",
    "score_candidate",
    "select_next",
    "Campaign",
]


# ---------------------------------------------------------------------------
# decision log
# ---------------------------------------------------------------------------

#: Decision kinds — every decision the loop makes is one of these.
ACQUIRE = "acquire"
EVALUATE = "evaluate"
STOP = "stop"
WAKE = "wake"


class DecisionReason:
    """Why a decision was taken. ``campaign_should_wake`` keys off it."""

    ACQUIRE = "acquire"
    EVALUATE = "evaluate"
    BUDGET = "budget"
    STAGNATION = "stagnation"
    FAILURE = "failure"


@dataclass
class CampaignDecision:
    """One logged decision: what was decided, why, and the data."""

    kind: str                    # ACQUIRE / EVALUATE / STOP / WAKE
    reason: str                  # a DecisionReason
    payload: Optional[dict] = None
    ts: float = field(default_factory=time.time)


class DecisionLog:
    """Append-only JSONL log of every campaign decision.

    One file per campaign; readable back without the campaign object,
    so the loop is auditable from its files alone.
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def record(self, decision: CampaignDecision) -> None:
        row = {
            "kind": decision.kind,
            "reason": decision.reason,
            "payload": decision.payload,
            "ts": decision.ts,
        }
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(row, sort_keys=True) + "\n")

    @staticmethod
    def entries(path: str | Path) -> list[dict]:
        out: list[dict] = []
        p = Path(path)
        if not p.exists():
            return out
        for line in p.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            out.append({
                "kind": row.get("kind"),
                "reason": row.get("reason"),
                "payload": row.get("payload"),
            })
        return out


def campaign_should_wake(*, reason: str) -> bool:
    """Should this decision wake the language model?

    The loop runs without one; the model is woken ONLY on failure,
    stagnation or budget exhaustion.  Everything else is the loop
    deciding for itself.
    """
    return reason in (
        DecisionReason.BUDGET,
        DecisionReason.STAGNATION,
        DecisionReason.FAILURE,
    )


# ---------------------------------------------------------------------------
# budget
# ---------------------------------------------------------------------------


@dataclass
class CampaignBudget:
    """Hard evaluation cap. ``spend()`` returns False once exhausted."""

    max_evaluations: int
    _used: int = 0

    @property
    def exhausted(self) -> bool:
        return self._used >= self.max_evaluations

    def spend(self) -> bool:
        """Consume one evaluation slot; False once the cap is reached."""
        if self.exhausted:
            return False
        self._used += 1
        return True


# ---------------------------------------------------------------------------
# objective
# ---------------------------------------------------------------------------


@dataclass
class TargetGap:
    """Objective on the HOMO-LUMO gap.

    With a target value, closeness to the target is what "better"
    means (a plain ``maximize`` degenerates toward overexposed
    derivatives and makes the Bayes-vs-random comparison unreadable).
    Without one, the goal is minimizing the gap.
    """

    target_eV: Optional[float]

    def score(self, gap_eV: float) -> float:
        """Lower score is better. 0 is the perfect score."""
        if self.target_eV is None:
            return abs(float(gap_eV))
        return abs(float(gap_eV) - self.target_eV)


def score_candidate(*, params: dict, value: float, target: TargetGap) -> dict:
    """Record one observation: parameters, value, score, together."""
    return {
        "params": dict(params),
        "value": float(value),
        "score": target.score(value),
    }


# ---------------------------------------------------------------------------
# acquisition
# ---------------------------------------------------------------------------


#: Fixed seed for the seed design AND the random baseline: a campaign
#: must be reproducible run to run — the Bayes-vs-random comparison is
#: only a statement about the acquisition if both sides sample the same
#: deterministic design.
_CAMPAIGN_SEED = 0


def _latin_hypercube_points(space_size: int, n: int,
                            rng: Optional[np.random.Generator] = None,
                            ) -> list[int]:
    """Seed design: ``n`` distinct indices via scipy qmc Latin hypercube.

    Falls back to even spreading if qmc is unavailable, so the module
    never hard-depends on a scipy minor feature.  Deterministic by
    default (fixed module seed) — pass ``rng`` for a different stream.
    """
    n = min(n, space_size)
    if n <= 0:
        return []
    try:
        from scipy.stats.qmc import LatinHypercube
    except ImportError:
        stride = max(1, space_size // max(n, 1))
        return sorted({min(i * stride, space_size - 1)
                       for i in range(n)})[:n]
    if rng is None:
        rng = np.random.default_rng(_CAMPAIGN_SEED)
    sampler = LatinHypercube(d=1, seed=rng)
    sample = sampler.random(n).ravel()
    return sorted({int(s * space_size) for s in sample})[:n]


def _seed_indices(remaining: Sequence, n: int) -> list:
    """Deterministic seed design over the REMAINING candidates.

    ``remaining`` is the candidate list after blocking observations and
    failures.  The fixed-seed LHS design is generated ONCE over that
    list — so with ``n=1`` per round, the next round gets the NEXT
    design point instead of re-drawing the same first point (a fixed
    seed with n=1 would otherwise re-propose the same candidate until
    the budget burns on one index).  Already-visited positions are
    consumed: the design state lives in the caller, not here.
    """
    design = _latin_hypercube_points(len(remaining),
                                     min(n, len(remaining)))
    return [remaining[p] for p in design]


def _gp_surrogate(x: np.ndarray, y: np.ndarray):
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, WhiteKernel
    kernel = Matern(length_scale=0.3, nu=2.5) + WhiteKernel(noise_level=0.05)
    # Training points are 1-D indices normalized to [0, 1]; normalize_y
    # keeps predictions on the observed score scale.
    return GaussianProcessRegressor(kernel=kernel, normalize_y=True,
                                    random_state=0).fit(x, y)


def expected_improvement(gap_eV: float, y_best: float,
                         y_mean: float, y_std: float) -> float:
    """EI for minimization of ``score`` = |gap − target|.

    ``gap_eV`` is not used in the formula itself — it documents the
    caller's contract: the surrogate models SCORES, not raw gaps, so
    the target mode is respected inside the surrogate too.
    """
    if y_std <= 0.0:
        return 0.0
    xi = 0.01
    z = (y_best - y_mean - xi) / y_std
    # Normal CDF/PDF via erf — no scipy.optimize dependency here.
    from math import erf, sqrt, pi, exp
    cdf = 0.5 * (1.0 + erf(z / sqrt(2.0)))
    pdf = exp(-0.5 * z * z) / sqrt(2.0 * pi)
    return (y_best - y_mean - xi) * cdf + y_std * pdf


def select_next(space: Sequence, observations: dict,
                n: int = 1, target: Optional[TargetGap] = None,
                seed_used: int = 0, excluded: Optional[set] = None,
                seed_offset: int = 0) -> list:
    """Pick the next ``n`` candidates.

    With no (or too few) observations, Latin-hypercube seeding over the
    REMAINING candidates; ``seed_offset`` advances the deterministic
    design so consecutive seeding rounds do not re-draw the same first
    point (a fixed seed with n=1 would otherwise re-propose the same
    candidate until the budget burns on one index).  With data, EI on a
    GP surrogate over the SCORES — never proposing a point already
    observed.  ``excluded`` marks candidates the loop already burned
    budget on WITHOUT getting a value (failed conversion, failed xtb):
    re-proposing them would spend the budget on a candidate that
    provably yields nothing.
    """
    blocked = set(observations) | set(excluded or set())
    untried = [i for i in range(len(space)) if i not in blocked]
    if not untried:
        return []
    if len(observations) < max(2, seed_used):
        # Advance the fixed design by the offset, then take the next n.
        design = _latin_hypercube_points(len(untried),
                                         len(untried))
        design = [untried[p] for p in design[seed_offset % len(design):]]
        return design[:n]
    target = target or TargetGap(None)
    x = np.array([i / max(len(space) - 1, 1)
                  for i in sorted(observations)]).reshape(-1, 1)
    y = np.array([target.score(observations[i])
                  for i in sorted(observations)])
    gp = _gp_surrogate(x, y)
    y_best = float(np.min(y))
    x_new = np.array([i / max(len(space) - 1, 1)
                      for i in untried]).reshape(-1, 1)
    mean, std = gp.predict(x_new, return_std=True)
    eis = [expected_improvement(space[i], y_best, float(m), float(s))
           for i, m, s in zip(untried, mean, std)]
    order = sorted(range(len(untried)), key=lambda j: -eis[j])
    return [untried[j] for j in order[:n]]


# ---------------------------------------------------------------------------
# the loop itself
# ---------------------------------------------------------------------------


def _xtb_binary_env() -> Optional[dict]:
    """PATH prefix carrying DELFIN's resolved xtb binary, if any.

    DELFIN resolves xtb through ``qm_runtime.find_tool_executable`` (see
    ``delfin/agent/doctor.py:_check_binaries``), but the native adapter
    ``delfin/tools/adapters/xtb_native.py:_run_xtb`` gates on
    ``shutil.which("xtb")`` only.  On nodes where xtb lives under
    ``~/.delfin/qm_tools`` and is NOT on PATH, the campaign therefore
    adds the resolved binary's directory to the child env's PATH —
    without touching the adapter.
    """
    try:
        from delfin import qm_runtime
    except Exception:
        return None
    try:
        path = qm_runtime.find_tool_executable("xtb")
    except Exception:
        return None
    if not path:
        return None
    parent = str(Path(path).parent)
    env = dict(os.environ)
    existing = env.get("PATH", "")
    if parent in existing.split(os.pathsep):
        return None
    env["PATH"] = parent + os.pathsep + existing
    return env


def _evaluate_with_xtb(geom_xyz: str, work_dir: str | Path,
                       *, charge: int = 0, mult: int = 1,
                       cores: int = 1) -> dict:
    """One chemistry evaluation: xTB single point, gap parsed by DELFIN.

    Returns ``{"value": gap_eV, "error": None}`` or
    ``{"value": None, "error": "..."}`` on failure.
    """
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    geometry = work / "geometry.xyz"
    geometry.write_text(geom_xyz, encoding="utf-8")

    from delfin.tools._registry import get as get_adapter
    adapter = get_adapter("xtb_sp")
    if adapter is None:
        return {"value": None, "error": "xtb_sp adapter not registered"}
    result = adapter.execute(
        work, geometry=geometry, cores=cores,
        charge=charge, mult=mult, method="gfn2",
    )
    status = getattr(result, "status", None)
    if status is not None and getattr(status, "name", "") == "FAILED":
        return {"value": None,
                "error": str(getattr(result, "error", "xtb failed"))}
    data = getattr(result, "data", None) or {}
    gap = data.get("homo_lumo_gap_eV")
    if gap is None:
        return {"value": None, "error": "no HOMO-LUMO gap in xtb output"}
    return {"value": float(gap), "error": None}


def _smiles_to_geometry(smiles: str, work_dir: str | Path) -> dict:
    """SMILES -> XYZ via DELFIN's existing converter (no rebuild).

    Uses ``delfin.smiles_converter.smiles_to_xyz`` (RDKit ETKDG, see
    delfin/manta/single_structure.py:349).  Returns
    ``{"xyz": str, "error": None}`` or ``{"xyz": None, "error": ...}``.
    """
    try:
        from delfin.smiles_converter import smiles_to_xyz
    except Exception as exc:                       # pragma: no cover
        return {"xyz": None, "error": f"converter import failed: {exc}"}
    try:
        xyz_content, error = smiles_to_xyz(smiles)
    except Exception as exc:
        return {"xyz": None, "error": f"SMILES conversion failed: {exc}"}
    if not xyz_content:
        return {"xyz": None,
                "error": str(error or "SMILES conversion returned nothing")}
    return {"xyz": xyz_content, "error": None}


def _evaluate_candidate(smiles: str, work_dir: str | Path,
                        log: Optional[DecisionLog] = None,
                        index: int = -1,
                        *, charge: int = 0, mult: int = 1,
                        cores: int = 1) -> dict:
    """Full candidate chain: SMILES -> XYZ -> xtb_sp -> gap.

    Every sub-step failure is returned as ``{"value": None, "error"}``
    so the loop logs it as a FAILURE and wakes the model — the only
    failure path that is allowed to.
    """
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    geom = _smiles_to_geometry(smiles, work)
    if log is not None:
        log.record(CampaignDecision(
            kind="convert",
            reason="ok" if geom["xyz"] else "failed",
            payload={"index": index, "smiles": smiles,
                     "error": geom["error"]}))
    if geom["xyz"] is None:
        return {"value": None, "error": geom["error"]}
    return _evaluate_with_xtb(geom["xyz"], work, charge=charge,
                              mult=mult, cores=cores)


@dataclass
class Campaign:
    """One campaign: space, objective, budget, log, loop."""

    space: Sequence[str]                  # e.g. substituent SMILES
    target: TargetGap
    budget: CampaignBudget
    log: DecisionLog
    work_dir: str | Path
    seed_used: int = 2
    stagnation_rounds: int = 3
    stagnation_tolerance: float = 0.01
    observations: dict = field(default_factory=dict)   # index -> gap_eV
    # Candidates whose evaluation produced NO value (failed conversion
    # or failed xtb).  Budget was burned; the acquisition must never
    # re-propose them.
    failed: set = field(default_factory=set)
    # How many points of the deterministic seed design are already
    # consumed — what makes consecutive one-point seeding rounds walk
    # the design instead of re-drawing its first point.
    _seed_offset: int = 0
    _best_score: float = math.inf
    # The scheduler is set by the caller (a workflow step or the agent
    # runtime).  Anything duck-typing Scheduler.schedule_once works --
    # tests pass a recorder; the real thing is
    # delfin.agent.scheduler.get_scheduler().
    scheduler: Any = None
    workspace: str = ""

    # --- decisions ---------------------------------------------------

    def _decide(self, kind: str, reason: str, payload: dict) -> None:
        self.log.record(CampaignDecision(kind=kind, reason=reason,
                                         payload=payload))

    def _wake_reason(self) -> Optional[str]:
        """The condition (if any) that justifies waking the model."""
        if self.budget.exhausted:
            return DecisionReason.BUDGET
        # Stagnation needs data: a campaign that has measured nothing
        # is not stagnated (best_score() is inf, and inf <= inf + tol
        # would otherwise fire "stagnation" on an empty loop).
        if self.observations and \
                self.best_score() <= self._best_score + self.stagnation_tolerance:
            return DecisionReason.STAGNATION
        return None

    def report_to_scheduler(self, *, delay_seconds: int = 60) -> Optional[dict]:
        """Set the LLM-free wake-up if (and only if) one is justified.

        The model is woken ONLY on failure, stagnation or budget
        exhaustion — the same rule ``campaign_should_wake`` pins.  A
        healthy loop never sets a wake-up and returns ``None``.  The
        wake-up itself goes through DELFIN's own
        ``Scheduler.schedule_once`` (delfin/agent/scheduler.py:360),
        with the campaign folder as workspace so the woken turn finds
        its own decision log and observations.
        """
        reason = self._wake_reason()
        if reason is None or not campaign_should_wake(reason=reason):
            return None
        self._decide(WAKE, reason, {
            "best_score": self.best_score(),
            "n_observations": len(self.observations),
        })
        if self.scheduler is None:
            return None
        summary = (
            f"campaign wake-up ({reason}): "
            f"{len(self.observations)}/{self.budget.max_evaluations} "
            f"evaluations done, best score {self.best_score():.3f} "
            f"(target {self.target.target_eV} eV). "
            f"Decision log: {self.log.path}"
        )
        return self.scheduler.schedule_once(
            delay_seconds=delay_seconds, prompt=summary,
            reason=reason, workspace=str(self.workspace or self.work_dir),
        )

    def best_score(self) -> float:
        if not self.observations:
            return math.inf
        return min(self.target.score(v) for v in self.observations.values())

    # --- the loop ----------------------------------------------------

    def step(self) -> list[int]:
        """One acquisition decision; returns the picked indices."""
        if self.budget.exhausted:
            self._decide(STOP, DecisionReason.BUDGET, {})
            return []
        picks = select_next(self.space, self.observations, n=1,
                            target=self.target, seed_used=self.seed_used)
        self._decide(ACQUIRE, DecisionReason.ACQUIRE, {"picked": picks})
        for p in picks:
            if not self.budget.spend():
                break
            self._decide(EVALUATE, DecisionReason.EVALUATE, {"index": p})
        return picks

    # --- chemistry evaluation ------------------------------------------

    def evaluate_candidates(self, *, strategy: str = "ei") -> list[dict]:
        """Run candidate evaluations until budget or stagnation stops.

        One cycle per call: acquire, evaluate, log.  Repeated calls
        continue an existing campaign (resume) and are bounded by the
        same budget.  ``strategy`` selects the acquisition family —
        ``"ei"`` (default) is the GP/EI loop, ``"random"`` is the
        fixed-seed Latin-hypercube baseline for the demo comparison.
        """
        if self.budget.exhausted:
            self._decide(STOP, DecisionReason.BUDGET, {})
            return []
        evaluations: list[dict] = []
        while not self.budget.exhausted:
            blocked = set(self.observations) | self.failed
            if strategy == "random":
                # Baseline: the fixed-seed LHS design over the remaining
                # space, advanced by the seed offset — no surrogate, no
                # EI.  The offset is what makes consecutive one-point
                # rounds walk the design instead of re-drawing point 0.
                untried = [i for i in range(len(self.space))
                           if i not in blocked]
                design = _latin_hypercube_points(len(untried), len(untried))
                if not design:
                    break
                picks = [untried[design[self._seed_offset % len(design)]]]
            else:
                picks = select_next(self.space, self.observations, n=1,
                                    target=self.target,
                                    seed_used=self.seed_used,
                                    excluded=self.failed,
                                    seed_offset=self._seed_offset)
            if not picks:
                break
            for p in picks:
                if self.budget.exhausted or not self.budget.spend():
                    break
                self._seed_offset += 1
                evaluations.append(self._evaluate_one(p))
            self._decide(ACQUIRE, DecisionReason.ACQUIRE,
                         {"strategy": strategy, "picked": picks})
            if self._stagnated():
                self._decide(STOP, DecisionReason.STAGNATION, {
                    "best_score": self.best_score()})
                break
        self._best_score = self.best_score() if self.observations \
            else math.inf
        return evaluations

    def _stagnated(self) -> bool:
        """Best score failed to improve over the last window."""
        if not self.observations:
            return False
        scores = [self.target.score(v)
                  for v in self.observations.values()]
        if len(scores) < self.stagnation_rounds:
            return False
        window = scores[-self.stagnation_rounds:]
        return max(window) - min(window) < self.stagnation_tolerance

    def _evaluate_one(self, index: int) -> dict:
        """Evaluate one candidate and log the observation."""
        work = Path(self.work_dir) / f"candidate_{index}"
        smiles = str(self.space[index])
        # The real chemistry: DELFIN's existing converter for the
        # geometry, then the registered xtb_sp adapter for the gap.
        result = _evaluate_candidate(smiles, work, self.log, index)
        self._decide(EVALUATE, DecisionReason.EVALUATE, {
            "index": index,
            "params": {"smiles": smiles},
            "value": result.get("value"),
            "error": result.get("error"),
        })
        if result.get("value") is None:
            # A burned candidate: budget was spent and no value came
            # back.  Record it as failed so the acquisition never
            # re-proposes it — re-running a provably useless candidate
            # is the fastest way to waste a budget cap.
            self.failed.add(index)
            self._decide(WAKE, DecisionReason.FAILURE, {
                "index": index, "error": result.get("error")})
            return {"index": index, "error": result.get("error")}
        self.observations[index] = float(result["value"])
        return {"index": index, "value": float(result["value"]),
                "score": self.target.score(result["value"])}

    def best_observations(self) -> dict:
        """index -> {params, value, score} for every observed candidate."""
        out: dict = {}
        for idx, value in self.observations.items():
            out[idx] = score_candidate(
                params={"smiles": str(self.space[idx])},
                value=value, target=self.target)
        return out
