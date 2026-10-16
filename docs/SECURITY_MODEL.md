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
| Secrets | The model provider's key, stored credentials (the OS keyring, or `~/.delfin/credentials.json`), SSH, cloud and package-index credentials are not readable by the agent's commands and do not reach a transcript. |
| Other users | Every store the agent writes for approvals, the audit log and the egress proxy is created owner-only (0600 files, 0700 directories); an approval answer is accepted only from a file owned by the same account. |
| The approval channel | A pending approval is answered by the person, never by the session that is waiting or by any process it started. |
| Git configuration | A command cannot plant a program that git runs later (hooks, filter and diff drivers, `core.hooksPath`) without the person being asked; a program entry a command adds to the repository config anyway is undone after the command; DELFIN's own git reads of the workspace run no program the repository names (`safe_git`). |
| Network egress | An isolated command, and every program it starts, reaches the network only through the egress proxy, which lets through an allow-list of domains and asks about the rest where someone can be asked; the cloud metadata address is refused; uploading workspace data is asked or, unattended, refused. |
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
apply on every host; the third depends on what the host offers. An
administrator's managed policy (section 6) sets a floor at each of them.

### 1. Permission mode and launch posture

`plan` (read-only), `default` (destructive actions ask), `diff_approval`,
`acceptEdits` and `bypassPermissions` (nothing inside the roots asks). A
session starts in `default` (`launch_guard.DEFAULT_MODE`) unless
`--permission-mode` or the `kit.default_mode` setting says otherwise.
`launch_guard.resolve_posture` refuses `--permission-mode bypassPermissions`
without `--unattended`, and never inherits a persisted `bypassPermissions`:
a terminal reaches bypass only by typing both flags. A one-shot run
(`chat -p`, a piped prompt, `delfin-agent run`) takes the same set-up as a
terminal session (grants, isolation flags, posture, session restore,
`--supervised`) before the engine runs; where a terminal falls back to
`default`, a one-shot run asking for bypass without `--unattended` exits
with status 2 instead of doing less than it was asked. Settings, hooks and
MCP servers that come from a checked-out repository are ignored until the
directory is trusted (`/trust`).

When no filesystem layer holds a shell command (none works on the host, or
`--no-isolate` / `agent.bash_isolation` = `off`),
`isolation_posture.apply_at_start` settles the start posture: an attended
`bypassPermissions` session asks once, naming the consequence (no keeps
`default`); with nobody to ask it starts in `default`, says so on stderr and
records `isolation_missing`. In that lowered session a command off the
auto-allow list that passed every rule is judged by the action classifier
and refused on any answer but allow
(`api_client._unisolated_fallback_verdict`). The banner, `/status` and the
turn status line show `unisolated` (`isolation_posture.isolation_status`).

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
- the git configuration files: in the work tree
  `api_client._GIT_CONFIG_PROTECTED_GLOBS` (`.gitattributes`, `.gitmodules`,
  a `.git` gitfile, `HEAD`); outside it every file git reads as
  configuration — the global and system files (`GIT_CONFIG_GLOBAL`,
  `~/.gitconfig`, the XDG location), the repository's config,
  config.worktree and info/attributes files, and every file an include,
  conditional include or `core.attributesFile` names, followed through
  includes. A
  write on any route (the file tools, a redirect, `tee`, `cp`, `sed -i`) is
  asked in every mode and refused where nobody can be asked; so are
  `git config` writes of keys that name a program and `git config --edit`.
  Both are recorded as `git_config_program`;
- the permission rules (`permission_rules`): a `permissions` block with
  `deny`, `ask` and `allow` lists in the user's settings file, in the forms
  `bash(<words>)`, `bash(<words>:*)`, `read(<glob>)`, `write(<glob>)`,
  `web_fetch(domain:<host>)`, an MCP tool or a bare tool name. deny refuses
  in every mode, `bypassPermissions` included; ask asks in every mode and
  refuses with nobody to ask; deny and ask match a command's words in every
  segment, substitution and `sh -c` body. allow only answers a question the
  gate would otherwise ask, after every hard check, and never covers a
  compound command, a wrapper (`env`, `xargs`, `sudo`) or inline code
  (`python -c`). A repository's `.delfin/settings.json` may add deny and ask
  rules; its allow list is ignored. `/permissions` lists the rules in force;
- the approval routes (`approval_guard`): a shell command that answers an
  approval, or a write into an answer location, is refused in every mode.
  An answer written outside the waiting process is accepted only with the
  one-time code (`approval_code`) that process showed on its human channel
  (the terminal frame, the dashboard row, or the controlling terminal of a
  `--supervised` run); the answer carries an HMAC of the code over the
  request id and the decision, never the code. A run with no terminal at
  all writes the code to an owner-only file under `~/.delfin/approval_codes`,
  which the path deny globs, `approval_guard` and every cage hide from the
  agent. A wrong answer is moved aside and recorded
  (`approval_code_refused`); the third expires the request;
- the egress scan: a shell command that sends data out is asked, and refused
  under `bypassPermissions`; network tools are asked outside bypass;
- the script content scan (`api_client._scan_bash_script_payloads`): a
  script, a `python -m` module or a makefile the command runs is read where
  the shell will be when it runs it, following `cd`, `pushd`/`popd`,
  subshells and `-C` (`api_client._ShellDirWalk`). Where that directory
  cannot be known (`cd $VAR` with the variable unset, `cd -` with no earlier
  `cd`) the command is asked in every mode and refused with nobody to ask
  (`script_dir_unknown`);
- the session tool lists (`--allowed-tools`, `--disallowed-tools`), which can
  only narrow.

A command off the auto-allow list is asked in `default` and `acceptEdits`,
except where the filesystem cage holds it
(`api_client._caged_command_may_run_unasked`): bubblewrap is the isolation in
use, the process cage is on, Landlock works, the network mode is not `open`,
and the text names no network client, git, package manager, job scheduler,
signal, privilege change, deletion or file authoring. Such a command runs
unasked (`caged_unasked`) with its reads confined by Landlock to the roots,
the system directories and the toolchain prefixes, the agent's own files and
the git metadata read-only, and no direct TCP route. `plan`,
`diff_approval` and `bypassPermissions` are unchanged by it;
`agent.caged_commands_run_unasked` = `false` turns it off.

Behind the rules sits the **action classifier** (`action_classifier.judge`):
a second model that is asked "did the user ask for this, and is it safe?"
about actions the rules allowed without a human seeing them. It can only add
friction — ask or refuse — and never grants what a rule refused. Its input
is the user's messages and the agent's tool calls; tool outputs, files and
web pages are never shown to it. It is on by default
(`agent.classifier.enabled`) with `agent.classifier.scope` = `unattended`
(bypass actions, messages to other agents, network calls); `all` also judges
every auto-allowed command; a session nobody can be asked in also judges
auto-allowed actions that leave the workspace (`git fetch`, a write to a
granted directory). `agent.classifier.user_rules` takes rules that start
with "never allow" or "always allow" and are applied without the model;
never-allow wins.

The classifier fails closed where nobody looks. When it gives no verdict (a
model error after one retry, a timeout, a garbled answer), an action that
leaves the workspace (`api_client._leaves_workspace`: network and remote
tools, messages to other agents, writes outside the workspace or into
agent-instruction and git configuration files, commands that reach the
network, and under bypass every command off the auto-allow list) is asked
when someone can be asked and refused when nobody can; work inside the
workspace keeps the rule result. After `api_client._CLASSIFIER_BREAKER_FAILURES`
failed calls in a row it is marked down (`classifier_down`), probed at most
once a minute, and `api_client.watcher_status_line` shows DOWN in the banner
and `/status`. `DELFIN_CLASSIFIER` overrides the switch for one process,
except in a process an agent session started (`DELFIN_AGENT_SESSION` is set
in every agent child), so the agent cannot switch off its own watcher.

### 3. Isolation of the shell command

`api_client._bash_isolation_argv` decides how a shell command is run.
`agent.bash_isolation` is `auto` by default: the command is isolated
**wherever this host can hold it**, in attended and unattended sessions
alike; `--isolate` forces isolation and refuses commands where the host
cannot provide it; `--no-isolate` (or `agent.bash_isolation` = `off`) runs
them unisolated. A session locked to one folder is always isolated, and its
commands are refused on a host where nothing can hold them.

`agent.bash_isolation` = `container` runs `bash`, `bash_background` and
`run_tests` in a container instead (`container_exec`, Apptainer/Singularity,
Podman or Docker, chosen by `agent.container.runtime`, image
`agent.container.image`). The container mounts the writable roots and the
git metadata read-write, the read-only roots, the git hooks and the
approval answer locations read-only, covers the credential locations of
`sandbox._HOME_SECRET_DIRS` that lie inside a mounted root, and mounts
nothing else of `$HOME`. It runs as the user with no capabilities,
no-new-privileges and a read-only root filesystem, never pulls an image at
run time, and has no network unless `agent.sandbox_network.mode` is `open`
(the egress proxy is not reachable from inside it, so `proxy` means none).
Where no runtime or image is available, the mode falls back to `auto`
below, with a recorded `isolation` event, never to a plain shell. A command
admitted unasked inside the cage (layer 2) is never sent to a container.

The mechanisms of `auto`, in order of preference:

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
  the credential store (`process_guard.put_exported_keys_away`), and
  every agent child gets an environment without provider keys and secrets
  (`api_client._scrubbed_bash_env`).

The credential store is the operating system's keyring (macOS Keychain,
Secret Service, Windows Credential Manager, through the optional `keyring`
package) where a secure backend answers in time, and the 0600 file
`~/.delfin/credentials.json` otherwise (`agent.credentials.store` = `auto`,
`keyring` or `file`, read only from the user's settings file, never from a
project's). The agent's tools reach neither store: the file and its name
index are deny paths and masked in the cages; the shell deny patterns
refuse the keyring clients (`secret-tool lookup`, `security
find-generic-password`, `keyring get`, the Secret Service bus name); the
session bus is unset and its socket refused by the socket guard; and every
value the keyring holds is redacted from tool results like any other
credential (`credentials.held_secret_values`). Moving keys between the
stores happens only on request (`delfin-agent credentials`).

### 4. The egress proxy

`egress_proxy` runs inside the agent process, outside every sandbox, on a
Unix socket in a 0700 directory with a per-process token. An isolated
command reaches the network only through it (`api_client._sandbox_network_args`):
`CONNECT host:port` and plain `http://` to the domains in
`egress_proxy.DEFAULT_ALLOWED_DOMAINS` (PyPI, conda, GitHub) plus
`agent.sandbox_network.allowed_domains` and the hosts the configuration
already names (the model endpoint, the package index hosts, the workspace's
HTTP(S) git remotes), on ports 80 and 443 or
`agent.sandbox_network.allowed_ports`. Direct TCP is opened only for
protocols that are not HTTP: the SLURM controller, SSH git remotes and
license servers. A command the gate let run unasked inside the cage gets no
direct TCP route at all.

The proxy (`egress_proxy.EgressProxy`) decides per session, through a
per-session token: a port outside the allowed set is refused; a
`web_fetch(domain:<host>)` permission rule denies, asks or allows next; a
known host passes, and an
upload to it above `agent.sandbox_network.upload_ask_bytes` (1 MiB by
default) is asked or, with nobody to ask, passes and is recorded
(`egress_upload`); an unknown host is asked when the session has an
approval channel and is not in `plan` or bypass, and refused otherwise with
a message naming the rule and setting that would let it through. Refusals
answer 403 and are recorded. `agent.sandbox_network.mode` is `proxy`
(default), `none`, or `open` (unrestricted; the metadata addresses are still
refused by the socket guard). The first value a process reads is a ceiling:
a stricter one applies at once, a looser one only at the next start, so a
command that edits the settings cannot open the network for its own session.

The proxy holds for programs a command starts, not only for its text. What
holds them depends on the host (`api_client.child_network_enforcement`): the
socket guard; a network namespace with a forwarder to the proxy inside
(bubblewrap without the socket guard); Landlock's TCP restriction (ABI 4:
TCP connects only to the forwarder's port, UDP not covered — partial);
Seatbelt. With `agent.bash_isolation` = `off`, `mode` = `open`, or on a host
with none of these, the command has the host's network, and the banner, the
security panel and `delfin-agent doctor` say so.

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
  else (mode `allowlist`) Landlock with the network denied as far as the
  host holds it, selected with `DELFIN_AGENT_SANDBOX`; `sandbox.confinement`
  names the state applied (`net=off`, `off-tcp` or `open`).

### 6. The managed policy

On a shared installation an administrator can set restrictions that no
user setting, project setting, flag or agent action loosens
(`managed_policy`). The file is `/etc/delfin/managed_policy.json` (Linux),
`/Library/Application Support/DELFIN/managed_policy.json` (macOS), or
a file of that name beside the installed package; `DELFIN_MANAGED_POLICY`
names one more file, which can only add restrictions. It applies at the
points of layers 1 to 4:

- `permissions` with `deny`, `ask` and `only_allow` lists, in the syntax
  of the user's permission rules and matched by the same implementation
  (`permission_rules.load` merges the managed deny and ask ahead of the
  user's and the project's). A managed deny refuses before the session
  lists, the mode, hooks, every allow and every approval, in every mode; a
  managed ask asks in every mode and refuses with nobody to ask.
  `only_allow` matches narrowly (`permission_rules.within_only_allow`) and
  bounds the auto-allow list, the user's allow rules and the caged
  admission alike;
- `disallowed_modes`: a disallowed mode becomes `default` (or `plan`)
  wherever a mode is set; an attended session that names one on the
  command line refuses to start (`managed_policy.launch_decision`);
- `require_isolation`: no shell command runs unsandboxed, `--no-isolate`
  and `agent.bash_isolation` = `off` are ignored, and on a host with no
  working layer an attended session refuses to start and an unattended one
  is capped at `default`; `container_backend` = `bwrap` accepts only
  bubblewrap as that layer;
- the domain ceiling (`allowed_domains_max` under `network`): the egress
  proxy and `web_fetch` refuse every host outside it before any user rule
  or known host, and `agent.sandbox_network.mode` = `open` becomes `proxy`;
- `disable_classifier_off`: the action classifier stays on against the
  settings and `DELFIN_CLASSIFIER`;
- `audit_log`: a copy of every audit line to the administrator's
  destination.

A file is trusted when it and every directory above it are owned by root
or by the account that installed the package and none is group- or
world-writable except a sticky directory. An untrusted file (user-owned,
writable by others, or named by the environment variable) is still read,
for restrictions only: several files merge by the union of restrictions
and the intersection of `only_allow` and domain ceilings, so it never
removes what a trusted file imposes, and it cannot set `audit_log`. A
trusted file that cannot be read or parsed fails closed (`bypassPermissions`
and `acceptEdits` disallowed). The policy file is unwritable through every
write route (`managed_policy.is_policy_path`); `delfin-agent capabilities`,
the start banner and `/permissions` name each file, its trust and what it
enforces. Its refusals are recorded as `managed_policy_deny` and
`managed_policy_ask`.

## What each layer guarantees on which host

"Isolated command" below means a shell command under `agent.bash_isolation`
= `auto` on that host. Layers 1, 2 and 5 hold on every host.

| Host | Writes | Credential reads | Sessions (tmux, SSH agent, systemd, D-Bus) | Network of an isolated command | Detached processes | Locked session / `--isolate` |
|------|--------|------------------|------|---------|-----------|------|
| **Container** (`agent.bash_isolation` = `container`, a runtime and image present) | roots + git metadata only; hooks and answer locations read-only; nothing else of `$HOME` mounted | covered where inside a mounted root, otherwise not mounted | not mounted | none; the runtime's own network only in `open` mode (the proxy is not reachable inside) | the container runs on if its Podman or Docker client is killed with SIGKILL; not in the process cage | runs in the container |
| **Linux, bubblewrap** (with seccomp user notification) | roots + git metadata only; `.git/hooks` and answer locations read-only | masked | masked; with seccomp user notification, connects also judged by the socket guard | egress proxy only; without seccomp user notification, a network namespace whose only way out is a forwarder to the proxy | killed with the command (PID namespace) | runs isolated |
| **Linux, Landlock** (bubblewrap refused; with seccomp user notification) | roots + git metadata + private temp + `/dev`; hooks writable (only the text check on hook paths stops a planted hook) | unreadable | refused by the socket guard | egress proxy only | can outlive the command; the lifeline and `delfin-agent stop-all` end what they find | runs isolated |
| **Linux, partial** (Landlock without seccomp user notification) | as above | unreadable | every Unix socket refused by a plain seccomp filter | Landlock ABI 4: TCP only to the proxy's forwarder port and the direct endpoints (by port, not address; UDP not covered); older ABI: **not restricted** | can outlive the command | **refused** while the network mode is `proxy` or `none` |
| **macOS, Seatbelt** | roots + git metadata + private temp; hooks write-denied | unreadable | Unix sockets only in the roots and the private temp | TCP only to the proxy's loopback port; in `open` mode the metadata address is not refused by address | can outlive the command; no process guard (Linux only) | runs isolated |
| **None** (no bubblewrap, Landlock or Seatbelt — for example Windows, or Linux with neither layer working) | **not restricted** — only the permission gate's reading of the command text | **not restricted** for the shell; the file tools keep their secret-path checks | where seccomp user notification exists, the session doors are refused by the socket guard; otherwise not restricted | egress proxy only, held by the socket guard, where seccomp user notification exists; otherwise **not restricted** | can outlive the command | **refused** |

On a host in the last row a `bypassPermissions` start is asked once when a
person is present; unattended it starts in `default`
(`isolation_posture.apply_at_start`), where a command off the auto-allow
list is asked or, with nobody to ask, judged by the action classifier. Both
record `isolation_missing`. A locked session refuses its shell commands
unless the person approves one command unconfined, and records the same
kind. `delfin-agent capabilities` prints which
row applies.

## Where a decision is recorded

Every refusal and every routed approval is recorded as a security event
(`security_events`): an in-memory ring buffer for the dashboard's containment
panel, and for blocked events also the durable attention inbox. The kinds
most relevant to this model:

- isolation: `isolation`, `isolation_missing`, `no_sandbox`, `no_sandbox_mcp`
- the gate: `deny_pattern`, `secret_path`, `self_mod`, `outside_ws`,
  `denied_path_via_bash`, `git_config_program`, `script_dir_unknown`,
  `caged_unasked`
- permission rules: `permission_rule_deny`, `permission_rule_ask`,
  `permission_rule_allow`
- the managed policy: `managed_policy_deny`, `managed_policy_ask`
- egress: `egress`, `egress_blocked`, `egress_asked`, `egress_upload`,
  `network`, `network_missing`
- approvals: `approval_self_answer`, `approval_from_outside`,
  `approval_code_refused`, `approval_timeout`
- the classifier: `classifier_ask`, `classifier_refuse`, `classifier_error`,
  `classifier_held`, `classifier_down`, `classifier_up`
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
  row Landlock ABI 4 restricts TCP by port, not by address, and does not
  cover UDP; with an older ABI the command has the host's network. Only a
  locked session or `--isolate` refuses instead.
- **Truncation on older kernels.** Landlock handles `truncate(2)` from ABI 3
  (Linux 6.2) on (`landlock_exec._write_rights`). On the Landlock rows with
  an older kernel a command can truncate a file outside the roots that the
  account may write, although it cannot otherwise write to it.
- **Same-uid limits.** The approval answer checks (owner, mode, age, id)
  prove that an answer came from this account; the one-time code proves
  that its writer saw the human channel. A process of the same account that
  DELFIN did not start and that can read that channel — the terminal's
  screen, the dashboard, or, for a run with no terminal, the code file —
  can still write a valid answer. Answers from outside a terminal session
  are recorded as `approval_from_outside`, not refused. Without the process
  guard (non-Linux) a command of the same user can read the agent's
  environment.
- **Classifier failures inside the workspace.** When the classifier cannot
  give a verdict, an action that stays inside the workspace keeps the rule
  result in every mode; only actions that leave it are asked or refused.
  Under bypass a narrow high-risk class (`api_client._bypass_failure_asks`:
  deletes, history rewrites, publishing, uploads, remote triggers, writes to
  agent-instruction or git configuration files) asks or is refused as well.
- **Commands run unasked inside the cage.** Under `default` and
  `acceptEdits` a command off the auto-allow list can run without a
  question where bubblewrap and Landlock hold it; what it can do is bounded
  by the cage, not by anyone reading it.
- **Git configuration.** The git configuration checks read the paths a
  write names; a program the command starts that writes a git configuration
  file is not seen by them. Three further layers cover that case:
  bubblewrap binds `hooks/`, `config`, the per-worktree config and `info/` of
  every git directory read-only, and Seatbelt denies writes below them
  (`_GIT_READ_ONLY_IN_THE_SANDBOX`; cost: `git config` from inside the
  sandbox fails); on every host, `git_config_watch` compares `config` and
  the per-worktree config before and after each shell command and background job
  and undoes a new or changed program-naming entry, keeping other changes;
  and DELFIN's own git reads of the workspace go through `safe_git`, which
  overrides every program-naming key the repository defines and runs git
  with no fsmonitor, hooks, pager, textconv or external diff. Remaining:
  on the Landlock row a hook or an `info/attributes` file can be written
  (Landlock cannot make a path read-only inside a writable root); the
  watch acts after the command, so a background job's entry exists until
  the job ends; and git run by the user in a terminal reads whatever the
  config holds at that moment. Modules not yet on `safe_git` are listed in
  `tests/test_delfins_own_git_reads_go_through_safe_git.py`.
- **Files reached through a symlink inside an allowed folder.** Any executor
  that walks a folder the gate admitted and opens what it finds can follow a
  link placed there to a file outside the roots, unless it checks each entry
  itself. The executors not yet shown to do so are listed in `OPEN` in
  tests/test_the_gate_scan_watches_the_real_executors.py.
- **MCP servers** without declared roots are outside every containment.
- **Containers.** Apptainer is covered only by argument tests, not by a
  real run; Podman and Docker run for real in CI. Rootless Docker's
  `--user` can map to a subordinate uid. A container outlives a client
  killed with SIGKILL.
- **The managed policy binds DELFIN's own gate**: native tools, MCP tools
  and the CLI backend's hook. A program the user runs outside DELFIN is
  governed by the operating system's permissions only. An untrusted policy
  file can only add restrictions, and the user can delete it.
- **The keyring.** On macOS a keychain item DELFIN created is readable by
  the same interpreter without a prompt; a command that hides the call
  from the shell deny patterns is stopped only where the cage blocks the
  keychain, and its output is still redacted. A keyring call that times
  out is abandoned, not cancelled.
- **Kernel and sandbox bugs.** Bubblewrap, Landlock, seccomp and Seatbelt
  are only as strong as their implementations. DELFIN does not defend against
  a kernel exploit.
- **What the model is shown is sent to the provider.** Credential redaction
  is pattern-based and best-effort.

## How to run safest

1. Keep isolation on (the default). Run `delfin-agent capabilities` once per
   host and check the isolation line names `bwrap`, `Landlock`, `Seatbelt`
   or a container;
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
   `delfin-agent credentials` (in the OS keyring where one is usable). Keep
   the process guard and the process cage on.
7. Declare `roots` for every MCP server, and trust (`/trust`) only
   repositories whose `.delfin/` settings and hooks you have read.
8. Read the containment panel (or the attention inbox) after an unattended
   run.
9. On a shared installation, put the restrictions every user must keep
   into a root-owned managed policy rather than into each user's settings.

## Changes in review

Open pull requests that change parts of this model are listed here. They
are not on `main` and nothing above relies on them; this section is updated
when they land.

None at present.
