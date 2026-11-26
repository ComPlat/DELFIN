"""A closed loop on the simulator: the measurement decides the next conditions.

Input: a kernel, an actor, the reactor device's id, the two condition channels
and a round budget. Output: a run record — every round's proposed conditions,
every write and what the kernel answered, what was measured — and one
``stopped`` reason.

One round is:

1. propose the next conditions *inside the envelope the manifest declares*
   (each channel's ``min``/``max``);
2. write them through :meth:`Kernel.set`, the one path to a driver's write.
   A setpoint further than the channel's ``max_step`` allows is reached by a
   RAMP of writes, because that is what the limit means: the reactor must
   not be shocked. A ramp is not a way around a refusal — each step is
   judged on its own, and a step the kernel refuses stops the loop;
3. write ``run=start``, read ``yield``, record the round.

What the loop does NOT do, and why each matters:

* it never writes a driver directly, so the manifest pin, the interlocks, the
  rate limit, ``max_step``, the lock and the e-stop apply to it exactly as
  they apply to a person's call;
* it never grants an approval. A ``hazardous`` channel stops the loop with
  ``approval_required`` and the request id, so a person can approve *that
  action* and the loop continue by hand. A loop that could answer its own
  approval request would make the kernel's approval step decorative;
* it never proposes outside the manifest. Clamping is the quiet form of
  retrying around the kernel: it turns "this step is forbidden" into "this
  other step ran instead", and the record would then show a measurement of
  conditions nobody chose. The proposal goes to the kernel as proposed, and a
  refusal is recorded with the proposed value still in it;
* it never retries a refused write. One refusal is one stop.
* it does not sleep on its own. The pause between rounds is injected, so a
  test runs 15 rounds instantly and the same loop can pace a real bench.

Stopping, in the order checked: ``estop`` (before each round, and a write
that lands during a stop is undone by the kernel and reported as ``estop``),
``refused`` (any other kernel refusal: interlock, out of range, rate, lock,
manifest), ``approval_required``, ``sink`` (the measurement could not be
read), ``budget`` (the round budget ran out). The record keeps the kernel's
``code`` and ``reason`` verbatim, so a reader never has to infer why a loop
halted.

The optimiser. Expected improvement on a Gaussian process (``GPBayes``): the
surfaces a self-driving lab meets are smooth and expensive, and EI is the
standard answer for that. Written on numpy, which DELFIN already requires
(``numpy>=1.26,<2.3``, BSD-3-Clause), rather than adding scikit-optimize or
botorch — an extra dependency for sixty lines of textbook GP, and a
dependency is a licence and a version to carry forever. Rejected
alternative: EI with a fitted mean function. On this two-dimensional ridge
with Halton starts it made no difference to the round count a zero mean did
not already reach, and a fitted mean is one more thing to be wrong about.

``GridSearch`` is the control arm: walk the lattice in a fixed order. Any
claimed benefit of the GP is measured against it, not asserted.
"""

from __future__ import annotations

import math
from typing import Any, Callable, Optional

import numpy as np

from . import surface as _surface

#: ``stopped`` values. A reader of a record uses these to tell "the search
#: finished" from "the kernel said no" without parsing reasons.
STOP_ESTOP = "estop"
STOP_REFUSED = "refused"
STOP_APPROVAL = "approval_required"
STOP_SINK = "sink"
STOP_BUDGET = "budget"


def _halton(n: int, base: int) -> np.ndarray:
    """The first ``n`` points of the van der Corput sequence in ``base``."""
    out = np.empty(n, dtype=float)
    for i in range(n):
        f, r, idx = 1.0, 0.0, i + 1
        while idx > 0:
            f /= base
            r += f * (idx % base)
            idx //= base
        out[i] = r
    return out


def _norm_cdf(z: np.ndarray) -> np.ndarray:
    """Standard normal CDF, from the stdlib's ``math.erfc``.

    numpy has no ``erfc`` and scipy is not a dependency the loop may assume
    (``import scipy`` fails on a bare installation), so the special function
    comes from ``math`` and is applied elementwise. The pool is a few
    thousand points per round, which is why a vectorised library is not
    needed for two cheap scalar calls.
    """
    flat = np.asarray(z, dtype=float).reshape(-1)
    out = np.fromiter((0.5 * math.erfc(-v / math.sqrt(2.0)) for v in flat),
                      dtype=float, count=flat.size)
    return out.reshape(np.shape(z))


def _norm_pdf(z: np.ndarray) -> np.ndarray:
    return np.exp(-0.5 * np.asarray(z, dtype=float) ** 2) / math.sqrt(2.0 * math.pi)


class GridSearch:
    """The control arm: propose lattice points in a fixed order."""

    name = "grid"

    def __init__(self, points):
        self._points = np.asarray(points, dtype=float).reshape(-1, 2)
        self._next = 0

    def propose(self, history: list[dict]) -> Optional[dict]:
        if self._next >= len(self._points):
            return None
        x = self._points[self._next]
        self._next += 1
        return {"setpoint": float(x[0]), "duration": float(x[1])}


class GPBayes:
    """Expected improvement on a Gaussian process over the condition space.

    Input: the envelope's ``lo``/``hi`` and the rounds so far. Output: the
    point with the largest expected improvement, or a spread-out Halton point
    while there is no data.

    ``length_scale`` is expressed as a FRACTION of each axis' range, so one
    number covers a degC axis and a minutes axis. It is a stated constant,
    not a silently tuned one; ``n_starts`` Halton points do the initial
    spreading, so the search does not depend on where it happens to start.
    """

    name = "gp-ei"

    def __init__(self, lo, hi, *, length_scale: float = 0.25,
                 n_starts: int = 5, jitter: float = 1e-6, density: int = 41):
        self.lo = np.asarray(lo, dtype=float).reshape(2)
        self.hi = np.asarray(hi, dtype=float).reshape(2)
        self.span = np.maximum(self.hi - self.lo, 1e-12)
        self.length_scale = float(length_scale)
        self.n_starts = int(n_starts)
        self.jitter = float(jitter)
        self.density = int(density)

    def _unit(self, points) -> np.ndarray:
        arr = np.asarray(points, dtype=float).reshape(-1, 2)
        return (arr - self.lo) / self.span

    def _pool(self) -> np.ndarray:
        steps = tuple(self.span / max(self.density - 1, 1))
        return _surface.lattice(tuple(self.lo), tuple(self.hi), steps)

    def propose(self, history: list[dict]) -> Optional[dict]:
        pool = self._pool()
        tried = {(round(float(r["setpoint"]), 9), round(float(r["duration"]), 9))
                 for r in history}
        # A point already run is not worth a run: EI at a sampled point has
        # ~zero variance, but the pool's lattice can still land on it, and a
        # repeated experiment spends a round and teaches nothing.
        keep = np.array([(round(float(a), 9), round(float(b), 9)) not in tried
                         for a, b in pool], dtype=bool)
        free = pool[keep]
        if not len(free):
            return None
        if not history:
            u = np.column_stack([_halton(self.n_starts, 2),
                                 _halton(self.n_starts, 3)])
            start = self.lo + u[0] * self.span
            if (round(float(start[0]), 9), round(float(start[1]), 9)) in tried:
                free_view = self._unit(free)
                start = free[int(np.argmin((free_view ** 2).sum(1)))]
            return self._conditions(start)
        xs = self._unit([[r["setpoint"], r["duration"]] for r in history])
        ys = np.array([float(r["yield"]) for r in history], dtype=float)
        # Standardise by the observed spread: the kernel amplitude is then 1
        # and the jitter has one meaning whatever the surface's units are.
        spread = float(ys.std()) or 1.0
        ys = (ys - float(ys.mean())) / spread
        diff = (xs[:, None, :] - xs[None, :, :]) / self.length_scale
        k = np.exp(-0.5 * (diff ** 2).sum(-1))
        k = k + (self.jitter + 1e-10) * np.eye(len(xs))
        try:
            chol = np.linalg.cholesky(k)
        except np.linalg.LinAlgError:      # pragma: no cover - guarded above
            chol = np.linalg.cholesky(k + 1e-6 * np.eye(len(xs)))
        alpha = np.linalg.solve(chol.T, np.linalg.solve(chol, ys))
        pool_u = self._unit(free)
        cross = (xs[None, :, :] - pool_u[:, None, :]) / self.length_scale
        ks = np.exp(-0.5 * (cross ** 2).sum(-1))
        mu = ks @ alpha
        v = np.linalg.solve(chol, ks.T)
        var = np.maximum(1.0 - (v ** 2).sum(0), 1e-12)
        sd = np.sqrt(var)
        z = (mu - float(ys.max())) / sd
        ei = sd * (z * _norm_cdf(z) + _norm_pdf(z))
        return self._conditions(free[int(np.argmax(ei))])

    def _conditions(self, x) -> dict:
        return {"setpoint": float(x[0]), "duration": float(x[1])}


class Loop:
    """The closed loop. ``run()`` performs it and returns the record."""

    def __init__(self, kernel, actor, *, device: str,
                 setpoint_channel: str = "setpoint",
                 duration_channel: str = "duration",
                 run_channel: str = "run",
                 result_channel: str = "yield",
                 optimiser=None, rounds: int = 15,
                 sleep: Callable[[float], None] = lambda _s: None,
                 settle_s: float = 0.0):
        self.kernel = kernel
        self.actor = actor
        self.device = device
        self.setpoint_channel = setpoint_channel
        self.duration_channel = duration_channel
        self.run_channel = run_channel
        self.result_channel = result_channel
        self.rounds = int(rounds)
        self.sleep = sleep
        self.settle_s = float(settle_s)
        if optimiser is None:
            lo, hi = self.envelope()
            optimiser = GPBayes(lo, hi)
        self.optimiser = optimiser

    def _device_spec(self):
        manifest = self.kernel.manifest
        if manifest is None:
            raise ValueError(f"no manifest: {self.kernel.load_error or ''}")
        dev = manifest.devices.get(self.device)
        if dev is None:
            raise ValueError(f"no device {self.device!r} in the manifest")
        return dev

    def envelope(self) -> tuple[tuple[float, float], tuple[float, float]]:
        """``(lo, hi)`` per AXIS: ``((sp_min, min_min), (sp_max, max_min))``."""
        dev = self._device_spec()
        lo: list[float] = []
        hi: list[float] = []
        for name in (self.setpoint_channel, self.duration_channel):
            ch = dev.writes.get(name)
            if ch is None or ch.min is None or ch.max is None:
                raise ValueError(
                    f"{self.device}.{name} needs a declared numeric min and "
                    "max to be searched")
            lo.append(float(ch.min))
            hi.append(float(ch.max))
        return (lo[0], lo[1]), (hi[0], hi[1])

    def candidates(self, steps: tuple[float, float] = (10.0, 10.0)) -> np.ndarray:
        """The lattice of the envelope, at the given steps."""
        lo, hi = self.envelope()
        return _surface.lattice(lo, hi, steps)

    def run(self) -> dict:
        dev = self._device_spec()
        record: dict[str, Any] = {
            "device": self.device, "optimiser": self.optimiser.name,
            "rounds_requested": self.rounds, "rounds": [], "stopped": None,
        }
        for n in range(self.rounds):
            if self.kernel.estopped(dev):
                return self._stop(record, STOP_ESTOP, code="estop",
                                  detail="an e-stop is active; no write was "
                                         "attempted")
            proposal = self.optimiser.propose(record["rounds"])
            if proposal is None:
                return self._stop(record, STOP_BUDGET,
                                  detail="the optimiser had no further point "
                                         "to propose")
            entry: dict[str, Any] = {"n": n + 1, "proposed": dict(proposal),
                                     "writes": []}
            record["rounds"].append(entry)
            stop = self._round(entry, proposal)
            if stop:
                return self._finish(record, stop)
            if self.settle_s:
                self.sleep(self.settle_s)
        return self._stop(record, STOP_BUDGET)

    def _finish(self, record: dict, stop: dict) -> dict:
        """End the record with the stop a round reported."""
        return self._stop(record, stop["stopped"],
                          code=stop.get("code"), detail=stop.get("detail"),
                          approval_request=stop.get("approval_request"))

    def _stop(self, record: dict, reason: str, **extra) -> dict:
        record["stopped"] = reason
        record.update({k: v for k, v in extra.items() if v is not None})
        # Give the device back. A loop that ended -- finished or refused --
        # must not leave the reactor locked to its own session: the next
        # person or job would be refused with `locked` by a session that no
        # longer exists. A hazardous device goes to its safe state first
        # (``Kernel.release`` does that), which is also the right end state.
        try:
            record["released"] = self.kernel.release(self.actor)
        except Exception as exc:            # never mask the stop reason
            record["released"] = None
            record["release_error"] = str(exc)
        best = best_measured(record)
        record["best"] = best
        return record

    def _round(self, entry: dict, proposal: dict) -> Optional[dict]:
        """Perform one round; return a stop dict or None when it completed."""
        stop = self._ramp(entry, self.setpoint_channel, proposal["setpoint"])
        if stop:
            return stop
        stop = self._write(entry, self.duration_channel, proposal["duration"])
        if stop:
            return stop
        stop = self._write(entry, self.run_channel, "start")
        if stop:
            return stop
        got = self.kernel.read(self.device, self.actor)
        if not got.get("ok"):
            return {"stopped": STOP_SINK,
                    "code": got.get("code") or "read_failed",
                    "detail": got.get("reason") or ""}
        rd = (got.get("readings") or {}).get(self.result_channel) or {}
        value = rd.get("value")
        if value is None or not rd.get("fresh"):
            return {"stopped": STOP_SINK, "code": "stale_or_missing",
                    "detail": f"{self.result_channel} is not a fresh reading "
                              "after run=start"}
        entry["setpoint"] = float(proposal["setpoint"])
        entry["duration"] = float(proposal["duration"])
        entry["yield"] = float(value)
        entry["mean_temperature"] = (
            (got.get("readings") or {}).get("mean_temperature") or {}).get("value")
        return None

    def _ramp(self, entry: dict, channel: str, target: float):
        """Reach ``target`` in steps the channel's ``max_step`` allows.

        A channel without ``max_step`` takes one write. The current value is
        read first; an unknown current value is written as a single step and
        the kernel decides (it refuses a step it cannot check).

        A target OUTSIDE the channel's declared range is NOT ramped: it is
        written as proposed, once, so the kernel refuses it as proposed. A
        ramp towards an unreachable target would produce a series of legal
        writes and then a refusal at the end, and the record would show a
        search that ran experiments it never chose.
        """
        spec = self._device_spec().writes[channel]
        if spec.min is not None and spec.max is not None and (
                not spec.min <= float(target) <= spec.max):
            return self._write(entry, channel, target)
        step = spec.max_step
        if step is None:
            return self._write(entry, channel, target)
        got = self.kernel.read(self.device, self.actor)
        cur = ((got.get("readings") or {}).get(channel) or {})
        value = cur.get("value")
        if not got.get("ok") or value is None:
            return self._write(entry, channel, target)
        value = float(value)
        if abs(target - value) <= step:
            return self._write(entry, channel, target)
        direction = 1.0 if target > value else -1.0
        while abs(target - value) > step:
            value += direction * step
            stop = self._write(entry, channel, round(value, 6))
            if stop:
                return stop
        return self._write(entry, channel, target)

    def _write(self, entry: dict, channel: str, value) -> Optional[dict]:
        """One kernel write, with the channel's OWN declared unit.

        The unit is read from the manifest rather than guessed: the kernel
        refuses a numeric setpoint without one, and guessing would be the
        loop supplying a fact it does not have.
        """
        spec = self._device_spec().writes[channel]
        out = self.kernel.set(self.actor, self.device, channel, value,
                              None if spec.is_state else spec.unit)
        entry["writes"].append({"channel": channel, "value": value,
                                "unit": None if spec.is_state else spec.unit,
                                "ok": bool(out.get("ok")),
                                "code": out.get("code")})
        if out.get("ok"):
            return None
        code = str(out.get("code") or "refused")
        stop = {"stopped": STOP_APPROVAL if code == "approval_required"
                else STOP_REFUSED, "code": code,
                "detail": out.get("reason") or ""}
        if out.get("approval_request"):
            stop["approval_request"] = out["approval_request"]
        return stop


def best_measured(record: dict) -> Optional[dict]:
    """The best completed round of a loop record, or None."""
    done = [r for r in record.get("rounds") or [] if "yield" in r]
    return max(done, key=lambda r: r["yield"]) if done else None


def rounds_to(record: dict, target: float) -> Optional[int]:
    """The round at which the best-so-far FIRST reached ``target``.

    This, not the final best, is what a comparison between two searches is
    about: how many runs it took before the search could have stopped with
    that answer.
    """
    best = -math.inf
    for r in record.get("rounds") or []:
        if "yield" not in r:
            continue
        best = max(best, float(r["yield"]))
        if best >= target:
            return int(r["n"])
    return None
