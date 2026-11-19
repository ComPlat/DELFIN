# TUI parity with `delfin-agent chat`

Phase 4's "done when": every command of the REPL's `/help` and every key of
its `/keys` works in `delfin-agent tui`, or the difference is written here.
The session behind the window is the REPL's own loop (`TerminalAgent`, in
the worker), so a typed line means in the TUI what it means in chat:
messages, `/commands`, `!shell`, `# notes`, and the numbered next steps.

Checked with the fake engine in tmux (`tui/test/e2e/`), in the headless
tests (`tui/test/`), and with a real KIT model (see the end).

## Commands (`/help`)

All commands run in the session, in Python, exactly as in chat. Their
output arrives as one notice, with its lines kept. The window adds pickers
for some, and runs three itself.

| Command | In the TUI |
|---|---|
| `/help` | as in chat |
| `/keys` | **the window's own**: its keys as bound (`tui.json` rebindings show) |
| `/editor [text]` | **the window's own**: it hands the terminal to `$VISUAL`/`$EDITOR` and takes it back. The chat rules apply: a private scratch draft pre-filled with `text`; a failed, empty or unchanged draft sends nothing |
| `/status`, `/cost`, `/usage`, `/context` | as in chat; the footer also shows this turn's tokens and cost, and the session's cost |
| `/compact`, `/export`, `/undo`, `/trace`, `/clear` | as in chat |
| `/exit`, `/quit` | as in chat: the session ends and saves, and the window leaves with it |
| `/model [name]` | with a name, as in chat. Bare, it opens the **model picker**: the provider's models, with where each runs and what it costs |
| `/mode [name]` | with a name, as in chat. Bare, it opens the **mode picker** |
| `/effort [level]` | with a level, as in chat. Bare, it opens the **effort picker** |
| `/permissions [...]` | as in chat (the posture and the rules). The posture is also in the palette's picker and on shift+tab |
| `/doctor`, `/fix`, `/init`, `/mcp`, `/trust`, `/tools`, `/hooks`, `/agents`, `/skills`, `/commands` | as in chat |
| `/session`, `/rewind`, `/tasks`, `/memories`, `/forget`, `/plans` | as in chat. To open a past session in the window: `/resume` or the palette's **resume picker** (the sessions running in the service first, then the workspace's saved ones) |
| `/undo-file`, `/pending`, `/approve`, `/reject` | as in chat |
| `/git`, `/bash`, `/attention` | as in chat |
| `/jobs` | the window's own screen: the agents' monitors, DELFIN's job monitor and each logged-in cluster's jobs, grouped (`jobs.list`); `/clusters` is its screen for the clusters and their login |
| user commands, skills, `/explore` `/review` `/plan` `/delegate` | as in chat (the order a typed line tries them in); all are completed as you type and listed in the palette |
| `!command`, `# note` | as in chat; a `!` at the start of the prompt shows shell mode (the marker `!` in its own colour, the box's border says so, no completion; backspace at the line's start or esc leaves it), and the line goes out as typed. A shell line is never steered into a turn: it runs when the turn ends. In the rows it reads `$ command` |

## Keys (`/keys`)

| Key | Chat | TUI |
|---|---|---|
| esc | end the turn | the same; with a question open it refuses (as `n`, as in chat) |
| enter | queue during a turn | the same; at the prompt it sends |
| ctrl+g | steer into the running turn | the same |
| shift+tab | next approval posture: default → acceptEdits → plan → bypassPermissions (bypass only where the policy allows it and isolation holds each command; its first time in a session asks for y, or a click on `[y] switch`, and switches nothing until then) | the same (from a posture outside the cycle: default) |
| ctrl+o | the last tool result in full | **Changed on purpose** (the user's ask, 2026-10-14): a tool call is its headline and what came back; ctrl+o shows or hides what every call printed (the common agent-CLI key), a click one call's. The last result in full: the palette. With a question open: its whole preview, as before |
| ctrl+t | the task list on/off | the task list stands above the prompt while a task is open (`tasks` on session.info/session.updated); ctrl+t folds it to one line |
| ctrl+l | redraw | writes the session again at the terminal's width (the last `replayTurns` turns) |
| backspace / ctrl+u | edit / clear the line | the same, plus ctrl+w, ctrl+k, ctrl+a/e, alt+b/f, undo ctrl+- |
| up / down | history; on an empty line with jobs listed, mark a job | history, from the prompt's first/last line. **Not here: marking a job** -- `/jobs` and `/bash` list them |
| tab | complete a /command or @path; on an empty line take the next suggestion | completes as you type, Tab takes the selected one. **Not here: Tab cycling the next-step suggestions** -- type the suggestion's number, as in chat |
| ctrl+r | history search | the same; enter takes the match into the prompt |
| ctrl+d | leave (empty line) | leave the window; the session keeps running in the agent service |
| ctrl+c | interrupt ladder | end the turn; else clear what is typed; else, twice, leave the window |

New in the TUI: shift/alt+enter and ctrl+j (and a trailing `\`) for a new
line; ctrl+p for the command palette; `N` to refuse a question with a reason
for the model; `?` lists every option of a question with what it does;
large pastes stand in as `[Pasted text #1 +N lines]`.

## Questions

| | Chat | TUI |
|---|---|---|
| confirm | `y n A k R e a ?` | the same keys, exactly the ones the session offers; `A`/`k` show the rule and its reach and need a second `y`; `N` = no + reason |
| plan | `d e n` | the same; `N` = reject + reason (it goes into the turn) |
| ask | `1`-`9`, esc dismisses | the same; multi-select toggles, enter sends |
| typeahead | -- | keys typed in a burst as a question opens are not taken as answers |

## Not the same, on purpose

- **Leaving.** In chat, leaving ends the session. In the TUI it never
  does: the session runs on in the agent service, and the next
  `delfin-agent tui` (or another window) finds it. `/exit` still ends it.
- **The banner** is chat's, minus its wordmark box (the window draws its
  own head).
- **Resizing** re-writes the session at the new width (chat leaves the
  old rows as they were). Under GNU screen, which cannot clear its
  scrollback, a divider is printed instead.
