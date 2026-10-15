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
| Bug report (account name, host name, workspace path, full session, touched text files; binary files only by size and sha256) | local archive; copied to the shared remote archive only with consent (`agent.bug_upload`) | the archive directory's group | not deleted, except the per-user fallback folder (30 days) |
| Shared knowledge entry (author handle, date, evidence pointers, text) | review queue or a git branch, then the repository | the queue's readers; after merge, everyone who can read the repository | git history keeps it |
| Git identity (`user.name`, `user.email`) of a share commit | the branch DELFIN prepares; the remote once you push it | everyone who can read the repository | git history keeps it |
| Status and reason of a decided share; the record file is then owned by the deciding admin's OS account | the queue record | the queue's readers | not deleted |

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
a summary (counts of messages, refused commands and files, the size, the
archive and the remote copy) and ask for confirmation; `--yes` files
without asking. The dashboard's Bug Report button opens a row that names
the local archive and, when an upload would follow, the remote host; Send
files the report.

A report contains:

- your OS account name, the host name, the workspace path, the session id,
  the provider, model, mode and permission settings, versions and the git
  commit;
- your description, the full conversation, the engine's message list, the
  system prompt (including the memories it carried), the tool trace with
  inputs and outputs, the refused commands, turn and request timings;
- an allow-listed excerpt of your settings, including the transfer target
  (host, user, remote path, port; no password);
- copies of up to 40 text files the session touched (at most 2 MB each,
  25 MB in total); files that look like secrets are skipped, and a binary
  file (a NUL byte in its first 8 KiB) is not copied but listed with its
  size and sha256.

Before the report leaves your account, the credentials DELFIN holds and
credential-shaped values are redacted from every text file of the report,
and the count is written to `REDACTIONS.txt`. This redaction does not
depend on the output-guard settings. If the redactor fails on a file, that
file is replaced by a placeholder that carries none of its text, and
`REDACTION_FAILED.txt` names each replaced file and the reason; the rest of
the report is kept (`delfin/agent/bug_report.py`, `redact_for_report`,
`_guarded_write`). The redaction is pattern-based: a credential DELFIN does
not hold and that has no known shape is not recognised.

Where it is stored and who can read it: the report is written to
`$DELFIN_BUG_ARCHIVE`, else `agent.bug_archive_dir`, else
`~/.delfin/agent_bugs`. It is made readable to the group that owns the
archive directory. When the settings name a transfer target
(`transfer.host` and `transfer.remote_path`), the report can be copied with
rsync over SSH into `<remote_path>/AGENT_BUGS` on that host, readable by
that directory's group. That copy has its own consent, `agent.bug_upload`
(also in the Settings tab):

- `ask` (default; a missing or unknown value counts as `ask`): copied only
  after a yes for this report. `/bug` and `delfin-agent bug report` ask a
  separate `[y/N]` naming `user@host:path/AGENT_BUGS`; in the dashboard,
  Send is pressed under the line naming the host. With nobody to ask
  (`--yes`, an unattended run) the report stays local, and the printed line
  says where it is and which setting would upload it.
- `never`: every report stays local, and no screen names a remote host.
- `always`: copied without asking, also in unattended runs.

Maintainer side: the optional bug watcher (`agent.bug_watcher.enabled`, off
by default) can send the text of new reports in the archive (the rendered
report, up to 8000 characters, and the tail of the tool trace) to the model
provider the maintainer configured, in a read-only session, to draft a
diagnosis. This is a transfer of the reporter's data to that provider, and
it has its own setting, `agent.bug_watcher.send_reports`:

- `ask` (default; a missing or unknown value counts as `ask`): `delfin-agent
  bug triage` and `bug watch --once` ask one `[y/N]` per run, naming the
  number of reports and the provider and model; the unattended watcher
  sends nothing.
- `never`: no report content is sent.
- `always`: sent without asking, also by the unattended watcher.

A report that is withheld stays pending and gets no diagnosis. Every start
of the watcher and of a triage prints where report content goes and under
which setting.

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
  When an admin decides, the record gains the status and the reason; the
  admin's account rewrites the file and then owns it. No admin name is
  written into the record.
  A knowledge admin is recognised by the OS account of the running process
  (the password-database entry of its uid, not the `$USER` variable) or by
  the GitHub login, never by an author handle.
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
| Knowledge review queue records, decided records and withdrawal tombstones | not deleted by DELFIN (*planned*: a retention period) |
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
- **Shared knowledge waiting for review**: `/unshare <queue id>` or
  `/unshare <lesson id>` (for a lesson whose status is queued) withdraws it,
  after you confirm with a one-time code. The record is rewritten as a
  tombstone that keeps only its id, its creation date, the status
  `withdrawn` and the withdrawal date; the text, the author handle, the
  evidence and the topic are removed. Only the OS account that queued the
  record can withdraw it: the record file's owner must be the account of
  the running process, and the handle you type is not consulted. A record
  that has already been decided cannot be withdrawn this way.
- **Shared knowledge in a repository**: `/unshare <entry>` prepares a
  branch that removes your entry, or only your provenance from an entry
  several people wrote. A removal does not erase the text from the git
  history.

Gaps (*planned*):

- an accepted or rejected queue record cannot be withdrawn or purged by its
  author (after the decision the file belongs to the admin's account); ask
  the knowledge admin;
- `/unshare <entry> as <handle>` for an entry in a repository checks the
  handle you type, not an authenticated identity; what stops a wrong
  removal is the review of the pull request;
- there is no single command that removes every entry, provenance line and
  queue record of one person; an admin removes them by hand.

## 8. Questions

Privacy questions and deletion requests: `<contact to be filled by the
maintainers>`. Security problems are reported as described in
[SECURITY.md](../SECURITY.md).
