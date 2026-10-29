# PR_NOTES — Endpoint-aware concurrency for model requests (task e01, wave3c)

## What was delivered in this worktree

**`delfin/agent/endpoint_gate.py` (new, complete, 22 tests pass).** A
host-wide, per-endpoint limit on concurrent model requests, shared by every
DELFIN process of the user:

- **Semaphore over file locks.** `<runtime_dir>/endpoint_gates/<key>/slot-N`,
  one slot file per slot, each taken with a non-blocking `flock`. The runtime
  directory is the service's own (`delfin.agent.service.rundir.ensure_runtime_dir`:
  `/tmp/delfin-<uid>`, 0700, refused unless safely the user's own), so every
  process of the user shares it and no other user can touch it. A slot held by
  a dead process frees itself — the kernel closes the flock'd descriptor when
  the process dies; there is no stale-holder sweep.
- **Keyed by endpoint base URL** (lower-cased, scheme-preserving via a
  `https_...` key, trailing `/` normalised). `DELFIN_ENDPOINT_PER_MODEL_LIMIT=1`
  adds the model to the key for endpoints that limit per model. Keys longer
  than ~160 chars get a sha256 suffix.
- **Adaptive limit (AIMD), shared.** Starts at `agent.endpoint_max_concurrency`
  (default 16; `0` disables the gate entirely). Halved on HTTP 429, on 503,
  or on a stall (`TimeoutError("no byte from the endpoint in Ns")` — exactly
  what the stream loop raises at api_client.py:32395), never below 1. Raised
  by 1 after 8 consecutive quick successes (first byte ≤ 10 s), never above
  the configured maximum; any congestion signal resets the climb. The state
  lives in `state.json` next to the slots, written atomically
  (`state_paths.write_text_atomic`) under a `state.lock` flock (bounded 5 s —
  a stuck state lock skips the update, it never blocks a request). The
  configured value is a CEILING: lowering it clamps at once, raising it lets
  the limit climb back.
- **Never deadlocks, never fails a request.** `acquire()` polls slot files,
  emits the wait notice `waiting for a slot on <base_url>: N in use` (surfaced
  by the stream loop's existing `waiting`/notice mechanism), and raises
  `SlotWaitTimeout` after a bounded wait (900 s default) instead of hanging.
  Release is `finally`-guaranteed via `with gate.slot()`. Every flock/OSError,
  a missing runtime dir, a corrupt state file → the request goes out ungated.
  A re-entrant `acquire` on one gate instance raises instead of leaking the
  outer slot's flock; a second request builds its own gate.
- **Config**: `agent.endpoint_max_concurrency` in the user settings (read via
  `delfin.user_settings.load_settings`, same pattern as other `agent.*` keys)
  plus env override `DELFIN_ENDPOINT_MAX_CONCURRENCY` (for tests/ops).

## Measured behaviour (all from runs in this worktree, venv python 3.12)

- `tests/test_the_endpoint_gate_limits_every_process_of_the_user.py`:
  **22 passed** (ruff F821/F811/F601/F402/F823: clean).
- Cross-process limit: 3 real child processes, limit 1 → measured maximum
  concurrent holders **1** (timeline sweep over actual hold spans).
- Contention + wait notice: parent holds the slot, 2 workers wait, parent
  releases → both workers acquired after >0.5 s wait (they really queued).
- Killed holder: child holds slot 0 of 1 → parent's acquire times out
  (SlotWaitTimeout); `kill` the child → parent acquires slot 0 within 5 s
  (the kernel freed the flock).
- AIMD shared state: 429/503/stall each halve 16→8 across processes through
  `state.json`; 8 quick successes climb 8→9; deterministic 400 leaves the
  limit untouched; climb capped at the configured max; halving resets
  `quick_run`.
- Limit 0 / no base URL: gate disabled, `slot()` yields None, request ungated.
- Regression: `test_a_stalled_endpoint_does_not_end_the_turn.py` (5 passed),
  `test_the_watcher_waits_as_long_as_the_endpoint_takes.py` (5 passed) — the
  suites that own the retry/stall behaviour the gate feeds back into.

## BLOCKED: the one remaining step — wiring the gate into `OpenAIClient.stream_response`

`delfin/agent/api_client.py` is on the agent's own protected-path list
(`matches_path_protected`, api_client.py:9857). Both `edit_file` and
`apply_patch` were refused: *"requires explicit user confirmation but no
confirm_callback is configured"* — the same block wave3/t02 recorded. I did
not route around it (no python-c writer, no bash rewrite). The orchestrator
can apply this exact diff; it touches three places, all in
`OpenAIClient.stream_response`:

```diff
--- a/delfin/agent/api_client.py
+++ b/delfin/agent/api_client.py
@@ (before `for _round in range(_MAX_TOOL_ROUNDS + 1):` at line 32085)
         # Host-wide, per-endpoint concurrency gate (delfin.agent.
         # endpoint_gate): every DELFIN process of the user shares one
         # limit per endpoint, so 60 sessions on one host no longer queue
         # their overload INSIDE the endpoint (measured 2026-10: requests
         # answered nothing for 600 s, every process retried on its own).
         try:
             from .endpoint_gate import (
                 for_client as _gate_for,
                 report_result as _gate_report,
                 SlotWaitTimeout as _gate_timeout,
             )
         except Exception:
             _gate_for = None
             _gate_report = None
             _gate_timeout = Exception

@@ (in the round loop, immediately before `_round_t0 = time.monotonic()`
    at line 32356 — i.e. after kwargs are built, before the request thread
    `_ct.start()` at 32374 spawns the request)
+            _gate = _gate_for(
+                str(getattr(self, "_base_url", "") or ""), self.model)
+            if _gate is not None:
+                try:
+                    _gate.acquire()
+                except _gate_timeout as _gw_exc:
+                    # Bounded wait expired with the endpoint full: same
+                    # path as a stalled request — transient, retried by
+                    # the existing ladder (523-2510 area).
+                    raise TimeoutError(str(_gw_exc)) from None
+                if _gate.wait_notice:
+                    yield StreamEvent(type="notice", text=(
+                        "\n⏳ " + _gate.wait_notice + "\n"))
             _round_t0 = time.monotonic()

@@ (release + feedback: wrap the existing round try/except — simplest
    correct form is two hooks)
   (a) in the success path, right before `_stream_attempt = 0` at line
       32639 (a round that produced a response):
+            if _gate is not None:
+                _gate.release()
+                _gate_report(_gate, exception=None,
+                             first_byte_s=(
+                                 _round_ttft - _round_t0
+                                 if _round_ttft is not None else None))
       (then delete the now-duplicated release below if using try/finally)
   (b) in the transient-retry branch (after `_stream_attempt += 1` /
       `continue` decision at 32503-32510) and on any re-raise out of the
       round (line 32557):
+            if _gate is not None:
+                _gate_report(_gate, exception=_stream_exc)
+                _gate.release()
```

Shape notes for whoever applies it:
- `stream_message` (the public entry, line 31429) funnels into this same
  round loop, so the wrap covers streamed and non-streamed fallbacks
  (32563, 32616) as long as those inner `.create` calls happen while the
  slot is still held — they are inside the same round.
- `_gate_for` returns None when disabled (limit 0, no base URL, no runtime
  dir): every hook above is then a no-op — the loop behaves exactly as
  today, which the existing api_client/engine retry+stall suites confirm.
- The stall→halve link works without further code: the stream loop raises
  `TimeoutError("no byte from the endpoint in Ns")` (32395), which
  `endpoint_gate.halve_for_exception` recognises verbatim.

## Test plan after applying the wiring (recommended order)

1. `pytest tests/test_the_endpoint_gate_limits_every_process_of_the_user.py -q`
2. `pytest tests/test_a_stalled_endpoint_does_not_end_the_turn.py tests/test_the_watcher_waits_as_long_as_the_endpoint_takes.py tests/test_a_turn_that_died_on_the_endpoint_continues_itself.py -q`
3. One end-to-end check with a fake OpenAI-compatible endpoint
   (`http.server` returning SSE) + 3 multiprocessing children driving
   `OpenAIClient.stream_message`, asserting in-flight ≤ limit and one
   `⏳ waiting for a slot` notice observed.

## Refusals the sandbox issued (exact)

1. `edit_file delfin/agent/api_client.py` → *"targets the agent's own safety
   layer … requires explicit user confirmation but no confirm_callback is
   configured"*.
2. `apply_patch` on the same file → identical refusal.
3. `bash` invoking the venv python with a pipe: blocked once with
   *"read denied: …/bin/python3.12 is outside the allowed workspace roots"*
   when the command would print the interpreter path — worked around by
   using `run_tests` instead; `venv/bin/python script.py` (without printing
   the path) is allowed and was used for the manual probe.

Scratch probes live in `./scratch/` (kh.py, killed-holder reproduction);
nothing was written outside the worktree.
