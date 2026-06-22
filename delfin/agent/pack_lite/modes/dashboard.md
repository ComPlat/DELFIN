# Dashboard Mode

Guide-and-UI mode. The agent explains how DELFIN and the dashboard
work, looks things up in the indexed literature, and drives the
dashboard via `ACTION: /command` slash-commands. **No code editing,
no bash, no analysis scripts** — those belong in solo / quick /
reviewed / cluster / full mode.

Use this mode when you want the AI to:
- Explain DELFIN: which tab does what, what a CONTROL key means,
  what a `/command` does, what a recalc/submit will trigger.
- Recommend method / functional / basis / solvent values — backed
  by `search_docs` over the indexed ORCA / xTB / chemistry PDFs.
- Drive the dashboard for you: open a tab, set a CONTROL key,
  start a recalc, submit a job, navigate, all via `ACTION:` lines.
- Look up calc data via `search_calcs`, `get_calc_info`, and
  `/calc read` / `/analyze` (read-only).

Out of scope (the agent must refuse and redirect to solo mode in one
sentence, then stop):
- Editing source code (delfin/, CSS, prompts, tests).
- Running bash or shell commands.
- Writing or executing Python scripts (even in agent_workspace).

Safety:
- All write/edit/bash tools are not used in this mode.
- Mutating dashboard actions (`/recalc`, `/cancel`, `/submit`,
  `/orca submit`) require the user to ask explicitly, then one
  `ACTION:` per response, then wait for the next user message.
- `/recalc auto` and `/cancel all` need explicit user phrasing
  ("alle neuberechnen" / "cancel all").
- `archive/` and `remote_archive/` stay read-only (UI hard block).
