# DELFIN Agent Memory Layers

Reference: what memory layers the DELFIN agent has, what each is for,
when it is retrieved, and where its limits are. Use this to decide where
a piece of information belongs and which retrieval path serves a query.

## Overview

| Layer | Module(s) | Written when | Retrieved when |
|---|---|---|---|
| Working memory | engine compaction + live-state block | every turn (transient) | every turn, automatic |
| Shared state | `project_memory` (DELFIN.MD / AGENTS.md) | by the user, by hand | every turn + every subagent, automatic |
| Facts (semantic) | `memory_store` (typed notes, BM25) | `remember` tool, `/memorize`, auto-distill | every turn via BM25 (score > 0), automatic |
| Procedure (skills) | skills + `skill_proposals` | `skill_propose_patch`, session-end learning | on `/skill` invocation, on demand |
| Episodes (past sessions) | `episodes` (markdown) + `session_index` (SQLite FTS5) | at session end, automatic | `episodes`: at session load, automatic; `session_index`: on `session_search` call, on demand |
| Curation / hygiene | `memory_tidy`, `memory_nudge` | proposals only (never deletes, retires by moving) | nudge: mid-session by work; tidy: manual CLI |

## Working memory — what currently applies

- **Task:** hold the state of *this* conversation: recent messages,
  live task list, current plan, uncommitted facts of the current job.
- **Retrieval:** always. Compaction keeps the tail verbatim and
  summarises the rest when the context window fills (95 % trigger).
- **Limits:** transient — nothing here survives the session. Anything
  durable must be promoted to facts (`remember`) or episodes (session
  end) explicitly. The `memory_nudge` prompt is the mid-session reminder
  for that promotion; it triggers on *work done* (tool calls + text),
  never on message count.

## Shared state — what every part must know

- **Task:** project-level rules that apply to every turn and every
  subagent: conventions, ownership, "how we work here".
- **Carrier:** `DELFIN.MD` / `DELFIN.md` / `AGENTS.md`, discovered by
  walking up from the working directory (deepest first, 8000-char
  budget, 6000 per file). Loaded by `project_memory.load_project_memory`.
- **Retrieval:** every turn and in every subagent prompt
  (`subagents.py`), automatic.
- **Limits:** small, hand-written, no ranking. This is not storage for
  facts or history — use facts for that. Other frameworks' instruction
  files are deliberately not loaded.

## Facts — semantic knowledge

- **Task:** durable, self-contained statements (user preferences,
  project constraints, failure→fix pairs), one per note, typed
  (user / feedback / project / reference).
- **Carrier:** markdown files with frontmatter in the per-project
  memory directory; `memory_store` reads/writes them atomically.
- **Retrieval:** every turn, automatically: `format_memory_context`
  ranks notes with BM25 over the current conversation and injects the
  relevant ones (score > 0). On demand via `session`-local recall in
  `history_search` (BM25 over this session's transcript).
- **Limits:** recall is keyword-driven (BM25) — a fact phrased with
  none of the current conversation's terms is not recalled. The store
  is curated, not accumulated: `memory_distill` extracts few facts per
  session, `memory_tidy` proposes merges/retirement (retiring moves a
  note aside, nothing is deleted). Without periodic tidy, near-duplicate
  notes can crowd the recall budget.

## Procedure — skills

- **Task:** validated, reusable playbooks ("how to do X here").
- **Carrier:** SKILL.md files, user-global (`~/.delfin/skills/`) or
  built-in (`delfin/agent/pack/skills/`); proposals pending review via
  `skill_proposals` (with evidence — no evidence, no proposal; never
  auto-activated).
- **Retrieval:** explicit — `/skill <name>` or a matching invocation;
  the agent lists them via `/skills`. Skills are NOT auto-recalled.
- **Limits:** a skill that is never invoked never influences behavior;
  the skill listing in the prompt is the only passive surface. Proposals
  need human acceptance before they change anything.

## Episodes — past sessions

Two carriers with a deliberate split:

1. **`episodes` (markdown, per project):** one compact record per
   finished session (`<memory dir>/episodes/<date>_<sid8>.md`, newest
   100 kept). Retrieved automatically at session load
   (`recall_episodes`, gated on relevance), so a new session starts
   with what recent sessions in this project concluded.
2. **`session_index` (SQLite FTS5, `~/.delfin/session_index.sqlite`):**
   full-text index over all archived sessions, queried through the
   `session_search` tool when the agent explicitly looks for how an
   *earlier* session solved something.

- **Split rule:** episodes = automatic short-term priming of *this*
  project; session_index = on-demand deep search across *all* sessions.
  When both match, episodes give the summary, session_index the detail.
- **Limits:** the index is built at session end (`index_at_session_end`)
  and swallows failures by design — a failed index run means those
  sessions stay unsearchable until a later index pass catches up
  (see the backfill on session start). `session_search` only sees
  sessions that were indexed.

## Curation — keeping the layers honest

- **`memory_tidy`:** proposes merging near-duplicates and retiring
  stale notes (moved aside, never deleted; un-landed suggestions are
  detected against the repo). Manual: `delfin agent memory-tidy`;
  a hint surfaces in the CLI when the budget is tight.
- **`memory_nudge`:** the mid-session, work-triggered reminder to save
  what was just learned before it is compacted away. Pure decision
  function, writes nothing.
- **Limits:** neither ever applies changes automatically — every merge
  or retirement needs acceptance.

## Cross-layer rules

- Promotion is explicit: working memory → facts via `remember` /
  `/memorize` (or the session-end distill), working memory → episodes
  automatically at session end.
- Nothing is ever deleted: retire = move aside (`memory_tidy`), skills
  change only through accepted proposals, episodes are capped by moving
  the oldest out.
- Automatic layers (shared state, facts recall, episodes priming,
  compaction) run every turn; on-demand layers (skills, `session_search`,
  `history_search`) run only when invoked — the agent must know to ask.
