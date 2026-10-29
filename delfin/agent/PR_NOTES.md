# PR_NOTES — Headless run continues an announced action

## What is being fixed

A headless DELFIN run (`delfin-agent chat -p … --output-format json`, and
the `run` one-shot door that shares the same path) used to stop for good
after exactly one `_run_once` turn. When the model's final answer ended by
announcing a step it did not take ("I'll continue the review now", "Next I
will …") and then produced nothing more, the run stopped — "done" — with
that intention as its entire output. In a 60-agent headless run several
agents stopped this way and produced no result at all.

An attended session has a guard for this: `_announces_an_action` +
the in-turn `_announced` auto-continue (`delfin/agent/api_client.py`) send
such a turn back once. But that path is part of `stream_response`'s
attended loop and never fires from a one-shot headless run, which calls
`_run_once` a single time and returns.

## Change

In `delfin/agent/cli.py`:

1. **`_run_once` tracks `tools_after_last_text`.** A per-turn flag set True
   on each `on_tool_use` and reset False on each text chunk. True at
   return means the turn's final text was followed by a tool call — the
   turn ended ON a tool round, not on its last message. Read back off the
   engine (a private attribute), not the returned dict, whose keys are
   pinned to `repl.TURN_KEYS`.

2. **`_continue_announced`** is a new loop invoked from `cmd_run` right
   after the first `_run_once`. While the turn's final answer announces an
   action and that final message carried no tool call, it runs one more
   `_run_once` with a nudge prompt, bounded by
   `agent.announce_continuations` (env `DELFIN_ANNOUNCE_CONTINUATIONS`
   over settings; default 5; 0 disables).

   Each continuation is a full turn: its tokens are billed, its rounds
   draw the per-turn tool-round budget, and `run_limits.ceiling_reason`
   counts it against the run's ceilings.

   The loop stops on:
   - a real final answer (the continuation's text does not announce an
     action — with or without tool calls),
   - a non-`done` stop reason or a turn error,
   - a reached run ceiling (checked before each continuation, so no nudge
     is spent against an already-spent ceiling),
   - a no-progress continuation that calls no tool and still announces —
     nudging again would only repeat the same sentence. A real-final
     no-tool continuation is the task finishing, not a stall: it is
     caught by the "real answer" stop, not by this one.

3. **`_resolve_announce_continuations`** returns the per-run continuation
   bound: env var first, then `agent.announce_continuations` from the
   user settings, else 5. Invalid values keep the default; 0 (via `"0"`)
   disables. Read per run; never written into anyone's settings file.

4. **Emit.** `cmd_run` passes `**_emit` to the first `_run_once` and
   `emit=_emit.get("emit")` to `_continue_announced`, once each; a
   continuation emits the nudge as a `{"type": "notice", …}` line so a
   stream-json reader can tell "the model's answer" from "the CLI nudged".
   No duplicated emit wiring around `_run_once`.

5. **Semantics of the decision — final message only.** The decision to
   continue is about the turn's FINAL assistant message, never the
   turn-wide `tool_calls` count. A run that made many tool calls and then
   closed on an announcement is exactly the measured failure; it continues
   because its last text was the announcement and no tool followed it. A
   turn that ended ON a tool round is already mid-work and gets no nudge.

**Attended sessions are unchanged.** The terminal and the dashboard drive
`repl.run_turn` / `stream_response`'s own guard; the headless loop lives on
`cmd_run` only, which attended sessions never pass through.

## Tests

`tests/test_a_headless_run_continues_an_announced_action.py` (new):

- `test_many_tools_then_announcing_close_gets_one_continuation` — the
  measured failure: one turn with tool calls, closing on an announcing
  text (no tool after it) → exactly one continuation, then the nudge's
  real answer ends the run.
- `test_continuations_are_bounded_by_the_configurable_cap` — progress
  turns that keep announcing stop at the configured cap (2), not at the
  10 scripted.
- `test_no_progress_stops_the_continuations` — a continuation that calls
  no tool and still announces is the stall: one nudge, then stop.
- `test_a_ceiling_ends_the_continuations` — a reached run ceiling stops
  before any continuation runs.
- `test_zero_disables_the_continuations` — `DELFIN_ANNOUNCE_CONTINUATIONS=0`
  disables.
- `test_a_real_final_answer_never_gets_a_continuation` — a final non-
  announcing answer is returned untouched.
- `test_an_error_or_engine_stop_is_left_alone` — an errored/stopped turn
  is not nudged.
- `test_resolution_default_env_and_settings` — resolution order: env over
  settings over default; invalid env keeps the default.

Fixtures use a fake engine (`_AnnouncingEngine`) and a scripted
`_run_once`, so no real ORCA/engine is involved. The announcement phrase
was chosen to be one the real `_announces_an_action` guard recognises
("I'll continue the review now"); the earlier draft phrase "Continuing the
review now…" ends on a plain sentence that the guard's last-sentence rule
does not classify as announcing, so it could not drive the loop — a
fixture bug, fixed.

Two pre-existing failures observed on this tree before/independent of this
change (background run of `test_orca_knowledge_c…`, same for
`test_knowledge_curation` / `test_knowledge_agent_surface`) are unrelated
to the headless loop; they are reported and not touched here.

## Rejected alternatives

- **Keying the decision on `out["tool_calls"]` being empty.** That is the
  turn-wide count and would have left the measured failure (tools, then an
  announcing close) unfixed.
- **Extending api_client's in-turn `_announced` guard into the headless
  path.** That guard serves attended sessions; the headless continuation
  loop belongs to the CLI (`cmd_run`), and touching api_client would reach
  attended behaviour that must stay unchanged.
