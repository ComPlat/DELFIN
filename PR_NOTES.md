# PR_NOTES — Reviewed lessons carried across tasks

## What this adds

`delfin/agent/lessons.py` turns a maintainer's Markdown lesson catalog into
untrusted data that reaches a new task's context automatically. A mistake
found once is not made again: the lessons whose scope matches the task are
injected, fenced as DATA, at prompt-build time.

**Input** (configured, off by default): `agent.lessons_dir` — a directory of
Markdown catalog files, one lesson per bullet, each bullet carrying a short
scope tag (inline `[tag1,tag2]` and/or a `## <tag>` section heading):

```
## orca-cards
- [orca-cards,dispersion] Match the dispersion correction to the functional.
## sandbox
- [sandbox] The scratch TMPDIR must resolve inside the bwrap layer.
```

**Output**: for each new task, at most 8 lessons (default) whose scope
matches the task, within 1600 characters, added to the volatile prompt tail
as a `--- Lessons (scope-matched, data) ---` block wrapped in an
`[UNTRUSTED EXTERNAL CONTENT — notes …]` fence (same helper as recalled
memory / briefing). Empty when unsupported, off, or nothing matches.

## Design decisions (and anything they ruled out)

- **Reuses the existing memory/briefing path** — `PromptLoader._load_lessons_context`
  is the sibling of `_load_briefing_context`, wrapped with the same
  `_fence_notes` untrusted fence and injected beside the briefing in
  `compose_sections`. No second recall mechanism; no new config parser
  (reads `agent.lessons_dir` from the existing `user_settings.load_settings`).
- **Selection is deterministic and model-free** — scope tags (+ hyphen
  fragments of them, so a task about "orca" surfaces an `orca-cards`
  catalog) matched against the lowercased task text; candidates ranked by
  matched-tag count, then file order, then line. No BM25/embedding, no LLM
  call, so earlier lessons can never steer the ranking.
- **Lessons are untrusted data, not instructions** — each body is stripped of
  control characters (`repl_render.strip_control`), wrapped in a nonce'd
  notes fence, and gated through the existing `instruction_shapes` injection
  guard (the same one the memory store uses): a lesson shaped like an
  authority ("ignore all previous instructions", a role marker, a standing
  consent) is refused and never reaches the prompt. Only files inside the
  configured directory are read; a lesson source carrying a path separator
  is dropped. No network.
- **The prompt's stable head is unchanged** — lessons are injected only in
  `LAYER_VOLATILE` (the per-turn tail), never in the cached stable prefix,
  so the prefix-cache key and the layer-order contract are untouched. Both
  injection sites in `compose_sections` (solo and subagent/other-role paths)
  follow the briefing.

## Measured prompt cost (with / without lessons, real pack)

`PromptLoader.prompt_size_report()` on the real solo pack (`role_id=solo_agent`,
`mode_id=solo`), same `task_text` with `agent.lessons_dir` off vs on. Size is
house-estimated at `chars // 4` tokens.

| task | lessons | total_chars | stable_chars | volatile_chars | lessons_chars |
|---|---|---|---|---|---|
| dispersion correction in ORCA geometry | off | 57 652 | 52 647 | 5 005 | 0 |
| dispersion correction in ORCA geometry | on | 58 644 | 52 647 | 5 997 | 990 |
| refactor pipeline, branch first, run tests | off | 57 274 | 48 778 | 8 496 | 0 |
| refactor pipeline, branch first, run tests | on | 57 274 | 48 778 | 8 496 | 0 |
| sandbox rejected scratch dir after git rebase | off | 60 793 | 54 371 | 6 422 | 0 |
| sandbox rejected scratch dir after git rebase | on | 61 533 | 54 371 | 7 162 | 738 |

Readings:

- **stable_chars is byte-identical on/off for every task** (52 647 / 48 778 /
  54 371): the lessons block never enters the cached stable head — exactly
  the requirement that the prefix cache keeps hitting.
- **The cost appears only when a lesson actually matches**: the pipeline
  task matched nothing in the sample catalog, and `on` adds 0 chars. The
  ORCA task matched the dispersion lesson (+990 chars ≈ 248 tokens); the
  sandbox task matched sandbox+git lessons (+738 chars ≈ 185 tokens).
  Both stay well under the 1600-char budget.
- Off-by-default: with no `agent.lessons_dir`, `_load_lessons_context`
  returns `""` and the prompt is unchanged.

## Tests

- `tests/test_lessons.py` (13): bullet/section-scope parsing, control-char
  stripping, selection by scope tag and tag-fragment, deterministic order,
  max-lessons and character-budget caps, injection-shaped lesson refused,
  plain imperative kept, path-traversal source dropped, off-by-default.
- `test_recalled_notes_are_fenced_as_data.py`, `test_briefing.py`,
  `test_agent_prompt_loader.py`, `test_lessons.py`: 102 passed.
- `test_prompt_cache_alignment.py`, `test_prompt_identity.py`,
  `test_a_prompt_head_that_does_not_move.py`, session-prompt and
  cold-request tests: 142 passed — the stable/volatile layer order and the
  unchanged prompt head are pinned.
- `test_a_memory_cannot_carry_an_instruction.py` (injection guard) and the
  wider memory suite remain green.

## Files

- `delfin/agent/lessons.py` (new) — parse + selection + guard.
- `delfin/agent/prompt_loader.py` — `_load_lessons_context`, injected in the
  volatile tail of `compose_sections` at both sites.
- `tests/test_lessons.py` (new).

Not changed: `delfin/manta/`, `delfin/fffree/`, `delfin/smiles_converter.py`,
`delfin/converter_backend.py`, the SMILES->3D builders,
`delfin/review*.py`, `delfin/dashboard/tab_review.py`, review sets, ratings.
