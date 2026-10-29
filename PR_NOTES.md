# PR — A model KIT stopped serving now switches to the next KIT model, bounded

Branch: `swarm/wave6/m01`. Two commits: RED test, then the fix.

## Problem

Measured 2026-10-29 for about an hour: KIT answered every request to one
model with HTTP 400 and the text

    Model 'glm-5.3-flash' is temporarily unavailable. Please try again later.

DELFIN classified that text as **gateway exhaustion** because
`"temporarily unavailable"` sits in `_INFRA_EXHAUSTION_MARKERS`
(`api_client.py`). `_is_transient_api_error` returned True, so the round
retried the *same* model up to `_STREAM_RETRY_MAX = 3` times and then
`raise`d — the model-switch path (`_kit_fallback_model`) is only reached
for "Model not found", and only once per turn. Ten unattended agents stood
still for an hour.

## Fix

Distinguish a **MODEL** being unavailable from the **gateway** being
exhausted. The two have different causes and different remedies.

- **Model unavailable** (`Model 'X' is temporarily unavailable`,
  `X is not available`): a property of the model, fixed by switching to the
  next family KIT serves.
- **Gateway exhaustion** (`too many connections`, `service temporarily
  unavailable`, `connection pool exhausted`, …): a property of the gateway;
  switching models does not help, so it keeps retrying the same model.

`api_client.py`:
- New `_is_model_unavailable_error(exc)` — 400 whose body pairs an
  availability phrase (`is temporarily unavailable`, `is not available`,
  …) with a **model scope** (the word "model", or a quoted model id). The
  gateway's own `service temporarily unavailable` names no model and stays
  a gateway retry.
- `_is_transient_api_error` now returns False for a model-unavailable
  message, so it does **not** path into the same-model transient retry.
- The stream handler (`stream_message`) switches via the existing
  `_kit_fallback_model` on both "Model not found" and "temporarily
  unavailable". An unavailable model gets **at most one short retry**
  (1.5 s, `_MODEL_UNAVAILABLE_RETRY_S`) on itself first — a blip — then
  switches, said in a stream notice exactly as "Model not found" was.
- `_model_fallback_done` (once per turn) became `_model_fallback_count`
  capped at **`_MODEL_FALLBACK_MAX = 3`** per turn. A headless run is one
  long turn; the cap lets it switch onward when the replacement is down
  too, without looping the roster forever; when the budget is spent the
  branch ends the round like any other unfixable model error.
- Fallback reuses `_kit_fallback_model`, which already limits itself to
  `model_capabilities._KIT_PREFERENCE` — only `kit.*` chat models
  (`kit.glm-5.3-flash`, `kit.glm-5.3`, `kit.deepseek-v4-flash`,
  `kit.mistral-small-4-119b-a8b`), never `azure.*`, `google.*`, `alias.*`,
  `standard-*`, and never embedding/reranker/speech/image. Verified:
  `_kit_fallback_model` returns `""` for any non-`kit` provider, and
  `kit_models_served` is built from `_KIT_PREFERENCE`. No change needed to
  `model_capabilities.kit_best_model` or the _KIT_PREFERENCE roster.

`turn_metrics.py`:
- `record_request` now carries a `model_switches` array
  `[{before, after, reason}]`, written to the `.requests.jsonl` so the
  switch is visible in the turn's metrics. The stream notice (a `notice`
  event) surfaces it live in `chat` text and `stream-json` output.

## Recording decision (rejected alternative)

The task asked to also show the switch in the JSON result of
`chat -p --output-format json`. An initial version added `model_switches`
to `repl.TURN_KEYS` and `_run_once`'s return dict. That broke the pinned
contract `tests/test_a_headless_run_reports_a_state_a_program_can_read.py`
(`test_the_result_keys_are_pinned`, `test_under_supervised_the_request…`,
`test_the_interactive_turn_carries_the_same_keys`): the interactive path
(`repl.TurnResult().to_dict()`) has no such key, and `TURN_KEYS` is pinned
to `headless == interactive == TURN_KEYS`. Extending the schema to three
places plus the documented headless schema in `docs/USER_MANUAL.md` was
out of scope for a bounded fix. Instead the switch is surfaced (a) in the
turn metrics `.requests.jsonl` (machine-readable) and (b) as the stream
`notice` event, which the `stream-json` format emits as a JSON line. The
reverted cli.py/repl.py edits are not in this branch.

## Tests

New `tests/test_a_model_kit_stops_serving_switches_models.py` (8 checks,
fake `create()` like `test_kit_runs_on_its_best_model…`):
- unavailable-model 400 → the round continues on the best other KIT model,
  notice says it, client ends on the new model, `_turn_model_switches`
  recorded;
- gateway-exhaustion 400 (`too many connections`) → retries the **same**
  model (no switch, no `_turn_model_switches`);
- the switch is bounded (`_MODEL_FALLBACK_MAX`, no infinite loop);
- classifier unit checks incl. `is not available` / `is currently
  unavailable`, and that gateway messages are not model-unavailable.

Green: `tests/test_a_model_kit_stops_serving_switches_models.py`,
`tests/test_transient_retry.py`, `tests/test_kit_runs_on_its_best_model…`,
`tests/test_turn_metrics.py` (33 pass together).

Run the fast gate: `pytest -q -m "not slow"` on this branch before merge.
The heavily-sandboxed `test_the_flags_a_run_can_be_bounded_by.py` (part of
the headless/JSON run) passed; the 3 headless-contract tests fail ONLY
while cli.py/repl.py carried the reverted `model_switches` key, and those
two files were `git checkout`'d back to HEAD, so the contract is intact
again.

## What was not changed

- `model_capabilities.py`: fallback scope already correct, no edit.
- The gateway-marker list: `too many connections` etc. still retry the same
  model.
- `docs/USER_MANUAL.md` headless schema (intentionally; see above).
