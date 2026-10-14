# DELFIN agent security model

This document describes the threat model of the DELFIN agent (`delfin-agent`,
the dashboard's agent tab and the scheduled runs) as the code on `main`
implements it: what is protected, against whom, which layer provides which
guarantee on which host, and what is not covered. Reporting a vulnerability
is described in [SECURITY.md](../SECURITY.md).

Every flag, setting, environment variable, module, symbol, command and
security-event kind named in backticks below is checked against the code by
`tests/test_the_security_model_names_what_exists.py`, so a rename on either
side fails the suite instead of leaving this page wrong.

To see what applies on your own host before a session starts, run
`delfin-agent capabilities` (the boundary, including the isolation mechanism
this host would actually use) and `delfin-agent gate explain -- <command>`
(the gate's verdict for one command, without running it). Both are described
in the [User Manual](USER_MANUAL.md#seeing-the-boundary-before-a-session-starts).

## What DELFIN protects

| Asset | The promise |
|-------|-------------|
| The operator's files outside the workspace roots | A command or file tool writes only inside the workspace and the directories granted with `--add-dir`; it reads outside only where `--read-dir` or an approval allows it. |
| Secrets | The model provider's key, stored credentials (`~/.delfin/credentials.json`), SSH, cloud and package-index credentials are not readable by the agent's commands and do not reach a transcript. |
| Other users | Every store the agent writes for approvals, the audit log and the egress proxy is created owner-only (0600 files, 0700 directories); an approval answer is accepted only from a file owned by the same account. |
| The approval channel | A pending approval is answered by the person, never by the session that is waiting or by any process it started. |
| Git configuration | A command cannot plant a program that git runs later (hooks, filter and diff drivers, `core.hooksPath`) without the person being asked. |
| Network egress | An isolated command reaches the network only through the egress proxy, which lets through an allow-list of domains; the cloud metadata address is refused; uploading workspace data is asked or, unattended, refused. |
| The agent's own safety layer | The files that define the permission system and the settings that decide what runs unattended are protected: writing them is asked in every mode, including `bypassPermissions`. |

## Against whom

**In scope: the model.** The agent executes text a language model wrote. The
model can be mistaken (it misreads a path, deletes the wrong directory,
"fixes" a test by editing it), and it can be prompt-injected: a file, a job
log, a web page, an MCP tool result or another session's message can contain
instructions that the model follows. Every layer below is built on the
assumption that a command may be adversarial even though the person never
asked for it, and that the person approves a command by its *text*, which an
interpreter, a symlink or an encoding can make differ from its *act*.

**Out of scope:**

- *A malicious co-tenant running as the same user.* Every process of the
  account — a terminal, an editor plugin, another agent — can read and write
  what the account can. DELFIN draws its boundary around the commands *it*
  starts; it cannot keep the account from itself. The approval channel is the
  clearest case: the answer checks establish that an answer came from this
  account, not that it came from a person (see residual risks).
- *Root, or a compromised account.* Anyone with root, or with the user's
  credentials, is outside every boundary described here.
- *The model provider.* Everything the model is shown is sent to the
  provider you configured; DELFIN redacts credentials from tool results
  (`output_guard.scrub_secrets` via `api_client._redact_tool_result`) but does
  not hide workspace content from the provider.

## The layers

A tool call passes these layers in order. The first two are text-based and
apply on every host; the third depends on what the host offers.

### 1. Permission mode and launch posture

`plan` (read-only), `default` (destructive actions ask), `diff_approval`,
`acceptEdits` and `bypassPermissions` (nothing inside the roots asks). An
interactive session starts in `plan` unless `--permission-mode` or the
`kit.default_mode` setting says otherwise. `launch_guard.resolve_posture`
refuses `--permission-mode bypassPermissions` without `--unattended`, and
never inherits a persisted `bypassPermissions`: a terminal reaches bypass
only by typing both flags. Settings, hooks and MCP servers that come from a
checked-out repository are ignored until the directory is trusted (`/trust`).

### 2. The permission gate

`api_client._run_permission_gate` decides allow / ask / refuse for every tool
call, before it runs. It applies, among others:

- the shell deny patterns (`api_client._DEFAULT_BASH_DENY_PATTERNS` plus the
  user's `kit.deny_patterns`, matched as case-insensitive regular
  expressions) — these refuse in every mode;
- the path deny globs for writes (`api_client._DEFAULT_PATH_DENY_GLOBS`:
  `.git/**`, `.env`, keys, `.ssh`, credential files) and the secret-path
  checks for reads;
- the protected paths (`api_client._DEFAULT_PATH_PROTECTED_GLOBS`: the
  permission layer, the confirmation UI, the settings modules): a write is
  asked in every mode, and refused where nobody can be asked;
- the git configuration files (`api_client._GIT_CONFIG_PROTECTED_GLOBS`:
  `.gitattributes`, `.gitmodules`, a `.git` gitfile, `HEAD`) and `git config`
  writes of keys that name a program: asked in every mode;
- the approval routes (`approval_guard`): a shell command that answers an
  approval, or a write into an answer location, is refused in every mode;
- the egress scan: a shell command that sends data out is asked, and refused
  under `bypassPermissions`; network tools are asked outside bypass;
- the session tool lists (`--allowed-tools`, `--disallowed-tools`), which can
  only narrow.

Behind the rules sits the **action classifier** (`action_classifier.judge`):
a second model that is asked "did the user ask for this, and is it safe?"
about actions the rules allowed without a human seeing them. It can only add
friction — ask or refuse — and never grants what a rule refused. Its input
is the user's messages and the agent's tool calls; tool outputs, files and
web pages are never shown to it. It is on by default
(`agent.classifier.enabled`) with `agent.classifier.scope` = `unattended`
(bypass actions, messages to other agents, network calls); `all` also judges
every auto-allowed command. `agent.classifier.user_rules` takes rules that
start with "never allow" or "always allow" and are applied without the model;
never-allow wins. `DELFIN_CLASSIFIER` overrides the switch for one process.

### 3. Isolation of the shell command

`api_client._bash_isolation_argv` decides how a shell command is run.
`agent.bash_isolation` is `auto` by default: the command is isolated
**wherever this host can hold it**, in attended and unattended sessions
alike; `--isolate` forces isolation and refuses commands where the host
cannot provide it; `--no-isolate` (or `agent.bash_isolation` = `off`) runs
them unisolated. A session locked to one folder is always isolated, and its
commands are refused on a host where nothing can hold them.

The mechanisms, in order of preference:

- **bubblewrap** — `/` read-only, the workspace roots and the repository's git
  metadata writable, a fresh `/tmp`, the credential locations of
  `sandbox._HOME_SECRET_DIRS` masked, `.git/hooks` and the approval answer
  locations bound read-only. Around every command (filesystem isolation or
  not) the **process cage** (`api_client._process_cage_options`): its own PID
  namespace, so nothing it starts outlives it; `--die-with-parent`;
  `--new-session` (no controlling terminal to type into); the systemd user
  manager, D-Bus, tmux and screen servers, `~/.ssh` and the SSH agent masked.
  `DELFIN_PROCESS_CAGE` = `off` in the environment of the process that runs
  the agent turns the cage off.
- **Landlock** (`landlock_exec`, Linux 5.13+, no user namespaces needed) —
  writes only beneath the roots, the git metadata, a private temp directory
  and `/dev`; the credential locations unreadable. The helper runs with
  `python -I`, imports nothing from DELFIN, and refuses the command (exit 126)
  if the policy cannot be applied.
- **Seatbelt** (`seatbelt`, macOS `sandbox-exec`) — writes only in the roots,
  the git metadata (hooks write-denied) and a private temp directory; the
  credential locations unreadable; Unix sockets only in the roots.

In `default` and `acceptEdits`, a command off the auto-allow list runs without
a question when bubblewrap with the process cage, Landlock read confinement
and a network mode other than open all hold it
(`api_client._caged_command_may_run_unasked`). Its word list
(`api_client._CAGED_UNASKED_EXCLUDED_WORDS`) only adds questions and carries
no security weight: the cage, the Landlock read confinement and the network
mode are what bound such a command. A bubblewrap that fails at run time (a
nested sandbox) runs nothing in its place, and withdraws this admission for
the rest of the process (`api_client._note_wrap_failure`).

Two guards complete the Linux mechanisms:

- The **socket guard** (`socket_guard.Supervisor`) receives every `connect()`
  of the command through seccomp user notification and performs the connect
  itself on a copy of the socket, so the address cannot change after the
  check. For an isolated command, Unix sockets are allowed only beneath the
  workspace roots, the temp directory and (under Landlock)
  `socket_guard.DEFAULT_ALLOWED`, the cluster's authentication and name
  services; for an unisolated command on a host without the process cage,
  every socket except the session doors (`api_client._guard_without_cage`).
  Abstract Unix sockets, datagram Unix sockets and io_uring are refused; the
  cloud metadata addresses and link-local ranges are refused in every
  network mode. If the supervisor dies, every later connect fails.
- The **process guard** (`process_guard.protect`) marks the terminal agent
  and the dashboard launcher non-dumpable, so a command of the same user
  cannot read their `/proc/<pid>/environ` or attach a debugger to read the
  key from memory. `DELFIN_PROCESS_GUARD` = `off` leaves them debuggable.
  At start-up a provider key exported in the shell environment is moved into
  the 0600 credential store (`process_guard.put_exported_keys_away`), and
  every agent child gets an environment without provider keys and secrets
  (`api_client._scrubbed_bash_env`).

### 4. The egress proxy

`egress_proxy` runs inside the agent process, outside every sandbox, on a
Unix socket in a 0700 directory with a per-process token. An isolated
command reaches the network only through it (`api_client._sandbox_network_args`):
`CONNECT host:port` and plain `http://` to the domains in
`egress_proxy.DEFAULT_ALLOWED_DOMAINS` (PyPI, conda, GitHub) plus
`agent.sandbox_network.allowed_domains`, on ports 80 and 443 or
`agent.sandbox_network.allowed_ports`, plus the SLURM controller so job
submission works. Everything else answers 403 and is recorded.
`agent.sandbox_network.mode` is `proxy` (default), `none`, or `open`
(unrestricted; the metadata addresses are still refused by the socket guard).
The proxy holds only for isolated commands: an unisolated command has the
host's network.

### 5. Other tools

- **Web fetch** (`web_tools`) refuses private, loopback and metadata
  addresses and re-checks every redirect target
  (`web_tools._GuardedRedirectHandler`).
- **MCP servers** run as processes DELFIN did not write a command line for;
  the shell isolation does not reach them. A stdio server declared with
  `roots`/`read_roots` starts in a namespace holding only those paths
  (`mcp_isolation`; `agent.mcp_isolation` = `builtin` does this for DELFIN's
  own servers). A server that declares neither runs uncontained, and the
  banner, `/mcp` and `delfin-agent doctor` say so.
- The dashboard's approval of a CLI-backend `Bash` call runs through
  `sandbox.run_agent_command`: an allow-list, then bubblewrap or firejail,
  selected with `DELFIN_AGENT_SANDBOX`.

## What each layer guarantees on which host

"Isolated command" below means a shell command under `agent.bash_isolation`
= `auto` on that host. Layers 1, 2 and 5 hold on every host.

| Host | Writes | Credential reads | Sessions (tmux, SSH agent, systemd, D-Bus) | Network of an isolated command | Detached processes | Locked session / `--isolate` |
|------|--------|------------------|------|---------|-----------|------|
| **Linux, bubblewrap** (with seccomp user notification) | roots + git metadata only; `.git/hooks` and answer locations read-only | masked | masked; with seccomp user notification, connects also judged by the socket guard | egress proxy only; without seccomp user notification, no network unless `mode` = `open` | killed with the command (PID namespace) | runs isolated |
| **Linux, Landlock** (bubblewrap refused; with seccomp user notification) | roots + git metadata + private temp + `/dev`; hooks writable (only the text check on hook paths stops a planted hook) | unreadable | refused by the socket guard | egress proxy only | can outlive the command; the lifeline and `delfin-agent stop-all` end what they find | runs isolated |
| **Linux, partial** (Landlock without seccomp user notification) | as above | unreadable | every Unix socket refused by a plain seccomp filter | **not restricted** | can outlive the command | **refused** while the network mode is `proxy` or `none` |
| **macOS, Seatbelt** | roots + git metadata + private temp; hooks write-denied | unreadable | Unix sockets only in the roots and the private temp | TCP only to the proxy's loopback port; in `open` mode the metadata address is not refused by address | can outlive the command; no process guard (Linux only) | runs isolated |
| **None** (no bubblewrap, Landlock or Seatbelt) | **not restricted** — only the permission gate's reading of the command text | **not restricted** for the shell; the file tools keep their secret-path checks | where seccomp user notification exists, the session doors are refused by the socket guard; otherwise not restricted | **not restricted** | can outlive the command | **refused** |

On a host in the last row an unattended (`bypassPermissions`) run still
proceeds and records `isolation_missing`; a locked session refuses its shell
commands and records the same kind. `delfin-agent capabilities` prints which
row applies.

## Where a decision is recorded

Every refusal and every routed approval is recorded as a security event
(`security_events`): an in-memory ring buffer for the dashboard's containment
panel, and for blocked events also the durable attention inbox. The kinds
most relevant to this model:

- isolation: `isolation`, `isolation_missing`, `no_sandbox`, `no_sandbox_mcp`
- the gate: `deny_pattern`, `secret_path`, `self_mod`, `outside_ws`,
  `denied_path_via_bash`, `git_config_program`
- egress: `egress`, `egress_blocked`
- approvals: `approval_self_answer`, `approval_from_outside`,
  `approval_timeout`
- the classifier: `classifier_ask`, `classifier_refuse`, `classifier_error`
- credentials: `key_exported`

Code-modifying tool calls are also appended to the audit log (`audit_log`,
`~/.delfin/audit.log`).

## Residual risks

- **Text-based checks without isolation.** On a host in the "None" row, and
  for any command run with `--no-isolate`, the only thing between a command
  and the filesystem is the gate's reading of its text. A command that builds
  a path or a program at run time (a variable, `base64 -d | sh`, an
  interpreter) can do what its text does not show. The same holds for
  approvals: the person approves the text.
- **The network of a partially isolated command.** On the "Linux, partial"
  row an attended or unattended command has the host's network; only a
  locked session or `--isolate` refuses instead.
- **Same-uid limits.** The approval answer checks (owner, mode, age, id)
  prove that an answer came from this account. A process of the same account
  that DELFIN did not start — and, without isolation, a command that strips
  `DELFIN_AGENT_SESSION` from its environment and builds the answer path at
  run time — can write an answer. Answers from outside a terminal session
  are recorded as `approval_from_outside`, not refused. Without the process
  guard (non-Linux) a command of the same user can read the agent's
  environment.
- **Classifier failures under bypass.** When the classifier cannot give a
  verdict under `bypassPermissions`, the rule result stands except for a
  narrow high-risk class (`api_client._bypass_failure_asks`: deletes, history
  rewrites, publishing, uploads, remote triggers, writes to agent-instruction
  or git configuration files), which asks or is refused.
- **Git configuration outside the workspace.** `~/.gitconfig`, a file named
  by `GIT_CONFIG_GLOBAL` and `git config --edit` are not covered by the git
  configuration checks; on the Landlock row, hooks are guarded by command
  text only.
- **Files reached through a symlink inside an allowed folder.** Any executor
  that walks a folder the gate admitted and opens what it finds can follow a
  link placed there to a file outside the roots, unless it checks each entry
  itself. The executors not yet shown to do so are listed in `OPEN` in
  tests/test_the_gate_scan_watches_the_real_executors.py.
- **MCP servers** without declared roots are outside every containment.
- **Kernel and sandbox bugs.** Bubblewrap, Landlock, seccomp and Seatbelt
  are only as strong as their implementations. DELFIN does not defend against
  a kernel exploit.
- **What the model is shown is sent to the provider.** Credential redaction
  is pattern-based and best-effort.

## How to run safest

1. Keep isolation on (the default). Run `delfin-agent capabilities` once per
   host and check the isolation line names `bwrap`, `Landlock` or `Seatbelt`;
   if it says `none`, use `--isolate` so commands are refused rather than run
   unisolated, or move to a host that offers one.
2. Run unattended only with `--permission-mode bypassPermissions --unattended`
   on a host with isolation, and give it a budget (`--max-budget-usd`,
   `--max-run-seconds`). For a run with no terminal that should still be able
   to ask, use `--supervised` and answer from `delfin-agent approvals`.
3. Grant the least: `--add-dir` for each directory it must write,
   `--read-dir` for what it only reads, `--allowed-tools` /
   `--disallowed-tools` to narrow the tool surface.
4. Write deny rules for what must never run: a regular expression in
   `kit.deny_patterns` (`~/.delfin/settings.json`) refuses a shell command in
   every mode; a "never allow …" line in `agent.classifier.user_rules`
   refuses a matching action without asking the model.
5. Leave `agent.sandbox_network.mode` at `proxy` and add only the domains a
   project needs to `agent.sandbox_network.allowed_domains`.
6. Do not export provider keys in the shell; store them with
   `delfin-agent credentials`. Keep the process guard and the process cage on.
7. Declare `roots` for every MCP server, and trust (`/trust`) only
   repositories whose `.delfin/` settings and hooks you have read.
8. Read the containment panel (or the attention inbox) after an unattended
   run.

## Changes in review

These open pull requests change parts of this model. They are not on `main`
and nothing above relies on them; this section is updated when they land.

- #201 — headless runs take the interactive checks; a session starts in
  `default`.
- #204 — structured permission rules (deny, ask, allow) in the gate.
- #206 — files git reads as configuration asked on every write route.
- #207 — an out-of-process approval answer carries a code shown to the
  person.
- #208 — an unisolated session says so and asks before bypass.
- #209 — the script content scan follows `cd`.
- #210 — the action classifier fails closed.
- #211 — child processes reach the network only through the egress proxy.
- #212 — commands off the auto-allow list run unasked when the cage holds
  them.
