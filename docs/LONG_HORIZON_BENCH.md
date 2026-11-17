# Long-Horizon Benchmark (F2)

A yardstick for compaction and long-session changes in DELFIN. Every
compaction change (F3, operator-B) must pass this benchmark before it is
adopted — adopt only if it measures better, never on argument.

This package only MEASURES compaction and long sessions. It changes
nothing in `in_turn_compaction.py`, `compaction_log.py`, `api_client.py`
or `engine.py`.

## The tasks

Six tasks, each needing a long session (100+ turns, several compactions),
built from DELFIN's real work. Each task's `lh_` block in
`delfin/agent/pack/benchmark/tasks_w28_lh.yaml` names its **planted facts**
(the values it must quote exactly, late), **decisions** and **refusals**
(early constraints that must still bind at the end), **expected_turns**,
**min_compactions** and **thresholds** (what "good enough to adopt" means).

1. `lh-code-change-fact` — a multi-file refactor audit of
   `delfin/orca_recovery.py`; three phrase→strategy-code mappings are given
   at start (facts.txt) and must be reprinted exactly from memory near the
   end (late probe).
2. `lh-calc-campaign-numbers` — a calculation campaign whose numbers must be
   quoted exactly late.
3. `lh-user-decision-binds` — a user decision early that must still govern
   the end.
4. `lh-refusal-holds` — a refusal early that must still hold late.
5. `lh-session-restart-continuity` — state that must survive a restart.
6. `lh-goal-focus` — work must stay on the task (goal drift measured).

## The metrics

Computed from the session files only — no model call to grade. Machine-
readable on disk as `.longhorizon_metrics.json` plus a short `report.md`:

- late recall of planted facts (exact numbers count only if exact),
- decisions / refusals still honoured,
- goal drift (work outside the task),
- redone work (the same tool call with the same arguments repeated after a
  compaction),
- tokens, cost, cache hit rate,
- number and size of compactions,
- wall time.

The producer is `delfin/agent/longhorizon/metrics.py`; the runner
`delfin/agent/longhorizon/runner.py` ties a task's `lh_` metadata, a
replayed (or live) transcript and the accept rule together and writes the
artifacts. The grade rule lives in
`delfin/agent/pack/benchmark/accept/lh_metrics_accept.py` (the harness's
`verify`) and is mirrored in `runner._grade` so the CLI and a test agree.

## Two modes

### (a) Replay mode — deterministic, CI

`python -m delfin.agent.longhorizon.runner <task_id> <transcript.json>
[workspace]` replays a scripted transcript through the REAL record writers
into a redirected store (never the real `~/.delfin`) and reduces it.
Deterministic; the CI control is a pair of transcripts of the same task that
differ only in whether a planted fact survives to the late records — the
recalling one must score higher (`tests/test_lh_replay.py`).

### (b) Live mode — operator-run against a real model

The OPERATOR runs the same `runner` CLI against a real model under its own
account and drops the resulting run directories into `.gate/lh_live/`. Live
runs are long (several hours); builders never start them themselves.

## Baseline (today's main)

Metrics of current `main` for at least two live tasks, so every future
compaction change has a number to beat. Filled from operator-run results in
`.gate/lh_live/` once they land:

| task | model | main commit | recall_late | redone_work | turns | compactions | wall_s |
|---|---|---|---|---|---|---|---|
| _pending_ | _pending_ | _pending_ | | | | | |

Method and commit are named with every number, and a metric is compared only
within the same method/commit.
