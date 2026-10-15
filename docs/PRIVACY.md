# Privacy notice

This notice describes which personal data the DELFIN software processes,
where it is stored, who can read it, how long it is kept, and how to
withdraw it. It describes the behaviour of the code on `main`; where the
code does not yet do what a privacy notice should be able to promise, the
gap is stated and marked *planned*.

DELFIN is software that you or your institution run. Whoever operates a
deployment (a research group, a computing centre) decides which model
provider, archive and shared folders are configured, and is the controller
for the data processed there.

**Controller and contact:** `<contact to be filled by the maintainers>`

## Summary

| Data | Where it goes | Who can read it | Kept |
|------|---------------|-----------------|------|
| Prompts, conversation, tool results, memory text, file contents the agent reads | the configured model provider | the provider, under its own terms | per the provider; locally see below |
| Session transcripts, tool traces, timing records | `~/.delfin/` on your machine | your account | 90 days (sessions), 30 days (derived state) |
| Memories, private lessons, outcome history | `~/.delfin/` on your machine | your account | until you delete them; outcome history: newest 1000 records |
| Bug report (account name, host name, workspace path, full session, touched files) | local archive, copied to the shared archive when a transfer target is configured | the archive directory's group | not deleted, except the per-user fallback folder (30 days) |
| Shared knowledge entry (author handle, date, evidence pointers, text) | review queue or a git branch, then the repository | the queue's readers; after merge, everyone who can read the repository | git history keeps it |
| Git identity (`user.name`, `user.email`) of a share commit | the branch DELFIN prepares; the remote once you push it | everyone who can read the repository | git history keeps it |
| Admin's OS account name (`decided_by`) | the queue record of a decided share | the queue's readers | not deleted |

DELFIN has no built-in channel that reports usage to its developers. The
files named "telemetry" and "metrics" in the code are local files under
`~/.delfin/`. Data reaches the maintainers only through a bug report copied
to an archive that your settings name (section 4) or a share you confirm
(section 5).

## 1. Data kept on your machine

Stored under `~/.delfin/` in your home directory:

- `agent_sessions/` — saved sessions: the conversation and the engine's
  message list, so that a session can be resumed.
- tool traces and turn and request timings, per session.
- `agent_memory.json` and `projects/<slug>/memory/` — memories the agent or
  you saved; `knowledge/private/lessons.json` — your private lessons.
- `outcome_history.jsonl` — one record per agent cycle: the task text
  (first 200 characters), provider, model, verdict, cost, duration and
  refused commands. The session briefing quotes recent failures from it
  into later prompts of the same account.
- `audit.log` and `failure_log.jsonl` — permission and failure records.

The memory store and the outcome history are written with mode 0600. The
session, trace and audit directories are reset to owner-only modes
(0700/0600) by the start-up maintenance (`delfin/agent/state_paths.py`).
The dashboard's CLI-backend runner appends every approved command to
`$XDG_CACHE_HOME/delfin/agent-audit.jsonl` (mode 0600).

## 2. What is sent to the model provider

Each model request carries the system prompt (role instructions, the
memories and knowledge entries selected for the task, the session briefing),
the conversation so far, and the results of the tools the agent called:
the contents of files it read, the output of commands it ran, search
results, and images you attach. Whatever the agent reads can therefore
reach the provider. Credentials DELFIN holds and credential-shaped strings
are redacted before sending; this redaction is pattern-based and
best-effort (`docs/SECURITY_MODEL.md`, residual risks).

With a local model server (Ollama, vLLM, LM Studio) these requests stay on
the machine or network that runs it. With a hosted provider (Anthropic,
OpenAI, KIT Toolbox) they are processed by that provider under your
account and its terms. Other outbound requests — web search and fetch,
browser libraries, one-time downloads — are listed in the README under
"Network access and offline use".

## 3. Shared folders and archives configured by a deployment

Two features write into locations shared with other people. Both are
configured in your settings; neither is active without that configuration.

- the bug archive (`$DELFIN_BUG_ARCHIVE`, `agent.bug_archive_dir`, and the
  remote `<transfer.remote_path>/AGENT_BUGS`), section 4;
- the knowledge review queue (`$DELFIN_KNOWLEDGE_QUEUE`,
  `agent.knowledge.queue_dir`), section 5.

## 4. Bug reports

A bug report is filed by you: the dashboard's bug-report button, `/bug` in
the terminal, or `delfin-agent bug report`. The terminal and CLI routes show
a summary (counts of messages, refused commands and files, the size and
the archive) and ask for confirmation; `--yes` files without asking. The
dashboard button files directly and shows where the report went.

A report contains:

- your OS account name, the host name, the workspace path, the session id,
  the provider, model, mode and permission settings, versions and the git
  commit;
- your description, the full conversation, the engine's message list, the
  system prompt (including the memories it carried), the tool trace with
  inputs and outputs, the refused commands, turn and request timings;
- an allow-listed excerpt of your settings, including the transfer target
  (host, user, remote path, port; no password);
- copies of up to 40 files the session touched (at most 2 MB each, 25 MB
  in total); files that look like secrets are skipped.

Before the report leaves your account, the credentials DELFIN holds and
credential-shaped values are redacted from the text files, and the count is
written to `REDACTIONS.txt`. Known limits: binary files are copied without
redaction, and if the redactor itself fails the text is written unredacted
rather than losing the report (`delfin/agent/bug_report.py`,
`_guarded_write`).

Where it is stored and who can read it: the report is written to
`$DELFIN_BUG_ARCHIVE`, else `agent.bug_archive_dir`, else
`~/.delfin/agent_bugs`. It is made readable to the group that owns the
archive directory. When the settings name a transfer target
(`transfer.host` and `transfer.remote_path`), the report is copied with
rsync over SSH into `<remote_path>/AGENT_BUGS` on that host, readable by
that directory's group. That copy is made whenever the target is set; it
is not asked separately (*planned*: a separate consent setting for the
upload, with a content preview in the dashboard).

Maintainer side: the optional bug watcher (`agent.bug_watcher.enabled`, off
by default) sends the text of new reports in the archive (the rendered
report, up to 8000 characters, and the tail of the tool trace) to the model
provider the maintainer configured, in a read-only session, to draft a
diagnosis. This is a transfer of the reporter's data to that provider.

## 5. Shared knowledge

Sharing a lesson is a separate act that only you can confirm, with a
one-time code shown on your own screen; no agent tool can start or answer
it. You choose an author handle (empty by default, `agent.knowledge.handle`),
and a check refuses credentials, email addresses, your home path, account
and machine names, and further identifiers for public targets
(`docs/USER_MANUAL.md`, "Shared memory").

What a shared entry carries: the lesson text, the author handle, the date,
and up to three evidence pointers. Where it goes:

- **Global layer**: a record in the review queue — by default
  `~/.delfin/knowledge/queue` (owner only), or a shared folder set by
  `$DELFIN_KNOWLEDGE_QUEUE` or `agent.knowledge.queue_dir`, where records
  are created group-readable (0640) so a knowledge admin can review them.
  When an admin decides, the record gains the status, the reason and the
  admin's OS account name (`decided_by`).
- **Project layer**: a commit on a new local branch
  `delfin-knowledge/<entry>` in that repository. The commit is made with
  your git identity (`user.name`, `user.email` from your git
  configuration). Nothing is pushed; once you push it, the commit and its
  author fields are visible to everyone who can read the remote.

## 6. Retention

| Location | Current behaviour |
|----------|-------------------|
| Saved sessions (`~/.delfin/agent_sessions`) | deleted after 90 days (`agent.session_retention_days`) |
| Tool traces, timings, subagent records, rotated audit logs | deleted after 30 days (`agent.state_retention_days`) |
| `~/.delfin/agent_bugs` (per-user fallback archive) | deleted after 30 days |
| Configured and shared bug archives, remote `AGENT_BUGS` | not deleted by DELFIN (*planned*: a retention period) |
| Knowledge review queue records | not deleted by DELFIN (*planned*: a retention period) |
| Memories, private lessons | kept until you delete them |
| Outcome history | the newest 1000 records |
| `audit.log`, `failure_log.jsonl` | trimmed by size and rotation |
| `$XDG_CACHE_HOME/delfin/agent-audit.jsonl` | not trimmed |
| Shared knowledge in a repository | permanent in git history |

## 7. Withdrawal and deletion

- **Memories**: `/forget <name>` deletes one memory; `/memory` lists them.
- **Local state**: the files under `~/.delfin/` are yours; deleting them
  removes the data from your machine. Data already sent to a hosted model
  provider is governed by that provider.
- **Bug reports**: a report in a shared or remote archive can be removed by
  someone with write access to that directory. Ask the deployment's contact
  (above) to delete a report you filed.
- **Shared knowledge**: `/unshare <entry>` prepares a branch that removes
  your entry, or only your provenance from an entry several people wrote.
  A removal does not erase the text from the git history.

Gaps (*planned*):

- an entry still waiting in the review queue cannot be withdrawn by its
  author; ask the knowledge admin;
- `/unshare ... as <handle>` checks the handle you type, not an
  authenticated identity;
- there is no single command that removes every entry, provenance line and
  queue record of one person; an admin removes them by hand.

## 8. Questions

Privacy questions and deletion requests: `<contact to be filled by the
maintainers>`. Security problems are reported as described in
[SECURITY.md](../SECURITY.md).
