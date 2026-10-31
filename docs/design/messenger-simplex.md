# SimpleX: private session chats and optional forum

Pair the phone once. Then type `/sync` in each desired TUI or dashboard session.
New sessions default to **off**. `/sync off` immediately stops phone routing for
that session; `/sync status` reports its setting. The setting is kept by session
ID, not by title. A resumed dashboard needs `/sync` to connect its live tab again.
An activated service session remains addressable while its worker is available.
The background gateway discovers enabled sessions and invites the paired user to
one private group per session, plus DELFIN Forum. Accept the SimpleX invitations.
Phone `sync` refreshes invitations; it cannot activate another session.
`neu NAME` on the phone starts a new session instead (see below).

## Everyday use

Ordinary text in a session chat addresses exactly that session. Several sessions
can run concurrently. Every notification includes the name and stable session ID.
`name Research` renames the session and its SimpleX group without changing routing.
A name in the forum changes the forum title. Renaming is not a recipient switch.

Commands below are plain text; a leading slash is optional:

| Text | Effect |
| --- | --- |
| `status` | Current short state, without asking the model |
| `update` | Ask the agent for a brief progress report |
| `pending` | Show currently open questions and approvals |
| `stop` | Cancel the current turn |
| `sessions` | List locally enabled sessions |
| `neu NAME` / `neu NAME: task` | Start a session in its own folder, invite to its chat, optionally with a first task |
| `name NAME` | Name this session/chat |
| `forum on` / `forum off` | Join/leave forum routing; in the forum itself, enable/disable it globally |
| `meldungen kurz` | Questions, approvals and brief completion/error messages (default) |
| `meldungen fragen` | Only questions and approvals automatically |
| `meldungen aus` | No automatic notifications; `status` and `pending` still work |
| `dateien` | List files explicitly made available by this session |
| `datei NAME` | Download one of those files |
| `austausch on` / `austausch off` | In the forum: enable/disable explicit agent announcements (default off) |

No complete transcript, reasoning, tool log or automatic file dump is sent.
An individual room can remain enabled while its forum participation is disabled.
Forum text addresses enabled participants only, with a delivery result for each.
`stop` in the forum cancels their current turns, not the service itself.

## Questions and approvals

Use SimpleX's **Reply** gesture on the concrete request, then type `yes`/`ja`,
`no`/`nein`, an offered option key or its full label. For a question, reply with
your selection or free text. No slash command or code is required. A bare `yes`
without a reply reference never approves anything. Images are not approvals; a
spoken decision word in reply counts like the typed one (see below).

The message shows the session, action/plan/question and the actual local options.
A yes selects the normal approval (for TUI plans, default mode); alternative modes
must be selected explicitly. Dashboard plans offer their existing execute option.
A shell command can be changed before it runs: reply `ändern: NEW COMMAND`
(also `change:` or `c:`, optionally in backticks), typed only. The message shows
the command back exactly as it will run; it runs once, as typed, and still passes
every check a command passes. A command with typographic quotes or dashes
(autocorrect) is refused with a hint.
Permanent rules require another reply to a separate confirmation message. The
complete concrete preview must fit within 1400 characters for phone approval;
longer previews remain available locally and can still be denied remotely.

Internally each outgoing request message maps to an exact session/request ID,
complete request fingerprint and activation generation. Codes expire after ten
minutes and are consumed before dispatch. Changed, resolved, expired or replaced
requests cannot apply to a later request. Disabling/re-enabling a session invalidates
previous requests. After a gateway restart use `pending` to get current messages.
There is no forum-wide approval. Compatibility `/approve CODE`, `/deny CODE` and
`/answer CODE text` paths retain the same binding checks.

Answers go through the existing service `permission.reply`/ProtocolBroker, or the
live dashboard's existing confirmation, question and plan handlers. That clears
the local pending request and wakes its waiting thread. A concurrent later answer
is rejected. Timeout/continuation behavior is the existing broker's behavior;
phone expiry does not grant the blocked action or resolve a different request.

## Images, files and voice

Attach images/documents/files in the intended session chat. The paired sender,
room, session activation and transfer generation are checked. Files are limited
to 8 MiB, four simultaneous transfers and 100 stored incoming files per session.
The gateway cache is capped at 100 entries and expires files after 24 hours.
Transfers cross the gateway/service user boundary in bounded chunks with a digest.

Incoming files are stored under the session workspace in
`.delfin-messenger/<session hash>/inbox`. Files are data, never automatically
executed. The agent receives the path and can use existing `view_image`,
`read_document` or `read_file` tools under its normal permissions. Image analysis
requires a vision-capable model. Downloads come only from the corresponding
`outbox`: ask the agent to prepare the desired artifact, then use `dateien` and
`datei NAME`. Images can also be downloaded as file attachments. Symlink,
hardlink, device, directory and traversal downloads are rejected. Files are sent
as binary data on explicit request, not automatically redacted like text summaries.

Voice transcription is optional and local. `sudo bash scripts/pi/pi_messenger.sh
--whisper` installs it: faster-whisper 1.2.1 in a venv of its own
(`/opt/delfin-messenger/whisper-venv`; it needs numpy 2, DELFIN keeps numpy below 2)
with `--require-hashes` from `whisper-requirements-<arch>.txt`, and the
faster-whisper-small model at a pinned revision, every file checked against its
SHA-256, under `/opt/delfin-messenger/models`. The gateway unit then carries
`DELFIN_WHISPER_MODEL` and `DELFIN_WHISPER_PYTHON`. No model download or external
speech API is invoked at run time; a separate isolated process (`python -I`) has a
120-second deadline.

A voice message does what the same words typed would do, and no more:
- without a reply reference it is the session's next prompt, and the phone is shown
  the transcript;
- as a reply to a question, the transcript is the answer;
- as a reply to an approval, only a transcript that is exactly a decision word
  (`ja`, `nein`, `yes`, `no`) or an offered option's label decides it, through the
  same path as the typed word; anything else (a sentence, "ja aber nicht jetzt")
  frees nothing and says what was understood;
- a lasting rule is never confirmed by voice.
Without working transcription, the recording is delivered as an audio attachment
with an explanation; the agent is not assumed to understand audio automatically.

## Explicit agent announcements

Normal agent output is never fed into another agent. With `austausch on` in the
forum, a participating agent can explicitly publish a short UTF-8 `forum.txt` in
its own outbox (maximum 1400 bytes). The gateway shares it with forum participants
as **another agent's information**, never a user command or permission decision.
The source session is named. Content fingerprints prevent repeat delivery, and
limits of one announcement per source per minute and ten per gateway per hour
bound automated chatter. The gateway does not request replies or relay transcripts.
Turn this off with `austausch off`, or remove a session with `forum off`.

## Sessions started from the phone

`neu NAME` starts a session in `~/DELFIN-Handy/NAME` of the agent-service user
(`DELFIN_MESSENGER_WORKSPACES` in the bridge unit). The phone supplies the name
only: it is reduced to letters, digits, `.`, `_` and `-` (no leading dot, 40
characters), an existing entry is reused only when it is a real directory, never
a symlink, and no folder, model or permission field is accepted from the gateway.
The session starts with the permission profile `default` -- it asks before
commands and edits, whatever the configured profile is -- and is enabled for the
messenger at once, so its chat follows without `/sync` at the computer. At most
three starts per ten minutes. `neu NAME: task` sends the task as the first prompt.
Provisioning with `--no-phone-sessions` leaves the folder root out; the bridge
then refuses every start.

## Setup

One run on the machine that hosts the agent service (a Pi, a server, a
workstation; any Linux with systemd and glibc >= 2.34):

```sh
sudo bash scripts/pi/pi_messenger.sh --whisper --pair
```

It installs the pinned simplex-chat build and lets its daemon create the bot
profile, makes the gateway venv (Python 3.11 from `MESSENGER_PYTHON`, uv or
`python3.11`; websockets with `--require-hashes`, DELFIN copied root-owned from
the checkout the script lies in), writes the units and, with `--pair`, prints the
invitation, as a QR code too when `qrencode` is installed, the link last (a
phone's terminal shows the end of the output). It waits up to 60 s for a daemon
still starting. Scan the code or paste the link in SimpleX; once the phone
connects, the script prints the security code (the chat core's `/_get code`,
best effort) to compare on the phone, and the bot sends a first message.
`--lang=en` sets every phone text to English (`DELFIN_MESSENGER_LANG`);
commands are understood in both languages (`notify short/questions/off`,
`exchange on/off`, `new`, `files`, `change:` beside the German words). `--name=LABEL` names the installation
(default: the short host name): the bot is `DELFIN-LABEL`, every session chat and
the forum carry the label, so one phone paired with a Pi, a server and a login
node tells them apart. The label is display only; rooms stay bound by chat and
session ID. The bot profile takes its name on first start only. Only that contact can complete pairing.
Without `--whisper` voice messages arrive as audio files; `--pair` can be given
on a later run. After a DELFIN update, run the script again: it copies the new
code and restarts daemon, bridge and gateway, never the agent service.

The agent service runs as user `delfin` (`MESSENGER_SERVICE_USER` names
another), as its own unit `delfin-agent-service` when that user has the DELFIN
venv `~/delfin/.venv` (the Pi layout; `MESSENGER_BRIDGE_PYTHON` names another
interpreter). The bridge uses the same interpreter. A service a terminal started
first is used as it is; the unit takes over when it ends. The daemon and
gateway use port 5227. The bridge runs as the agent-service user with that user's own
group: the agent service lets a client act only when it can read the client's
`/proc/<pid>/environ`, which Linux allows only with uid and gid equal (with group
`delfin-gw` every phone command was refused, -32003). Its 0660 socket is under
`/run/delfin-messenger-bridge`, set to group `delfin-gw` and setgid (2750) by
the unit as root, so the gateway reaches the socket and nothing else; the service
user's home is the bridge's one writable path. The service runtime directory
stays private. Activate desired sessions locally with `/sync`.
The direct contact remains a setup channel. `revoke` with the same environment
removes routing authorization immediately; already dispatched commands are not
undone. `status` lists stored chat bindings.

## Security and validation

Only the paired human's received messages are accepted, including inside groups.
Other members, sent-message echoes and unbound rooms cannot command sessions.
UID and agent-process checks protect local sockets. The daemon's API port is
kept for `delfin-gw` by an nftables rule stored in
`/etc/delfin-messenger/messenger.nft`: the daemon unit loads it before every
start (as root, and does not start if it fails), and the host ruleset includes
it, so neither a reboot nor a reload of `/etc/nftables.conf` (which begins with
`flush ruleset`) leaves the port open to the agent's user. Group names stay
within the daemon's 50-character limit; a name it refuses is retried once with
the `validName` it offers. A failed command is logged with its reason
(`journalctl -u delfin-messenger-gateway`) and told on the phone in words. The activation store is
excluded from agent file tools and shell sandboxes. Room IDs cannot be rebound to
another session. Rate/size limits and persistent incoming message IDs prevent
replay across gateway restarts. Dispatch is at most once; after a connection
failure, check `status` before retrying. Partial broadcasts are reported.

Wire commands and envelopes follow pinned SimpleX v7.0.3 primary references:
[commands](https://github.com/simplex-chat/simplex-chat/blob/v7.0.3/bots/api/COMMANDS.md),
[types](https://github.com/simplex-chat/simplex-chat/blob/v7.0.3/bots/api/TYPES.md),
[events](https://github.com/simplex-chat/simplex-chat/blob/v7.0.3/bots/api/EVENTS.md).

Automated tests cover simulated SimpleX traffic, actual local sockets, the real
DELFIN service/approval broker with a fake model, dashboard gate resolution,
activation, stale replies and file isolation. A live SimpleX daemon/phone test,
real media transfers and an installed speech model remain deployment acceptance
checks. These tests do not provision or deploy services or certify compliance
with a particular external security standard.
