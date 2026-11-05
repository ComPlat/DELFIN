# SimpleX: private session chats and optional forum

Connect your phone in three steps, no admin needed:

1. In your session type `/pair` (terminal: `delfin-agent messenger pair`).
2. Scan the QR code in SimpleX (or paste the link: + > paste link).
3. Type the 6-digit code DELFIN sends to the phone: `/pair 123456`.

Then `hol NAME` / `get NAME` on the phone brings a running session to it,
`neu NAME` starts one, and `/sync` (or `/sync NAME`) in a session connects
it. `hilfe` / `help` on the phone lists the commands.

## Several people, one account (people mode)

On a box in the launcher's people mode every person has `~/people/<name>` as
home, with their own agent service (`~/people/<name>/.run`), sessions,
activation records and phone folder (`~/people/<name>/DELFIN-Handy`). Each
person pairs their own phones; each phone is bound to exactly one person and
reaches only that person's sessions, chats, approvals, files and forum.

**What this boundary is.** Everyone in people mode is the same Unix user.
The separation is enforced by the gateway and the bridge (application level),
not by the operating system: a process of the account can read every person's
files. It keeps phones from reaching another person's sessions by routing or
mistake; it does not protect people of one account from each other.

How a pairing works:

- `/pair` writes a request into `~/.delfin/messenger_pairing` of the person's
  home (0700 directory, 0600 files, masked from agent file tools and shells
  like the activation records). It carries `sha256(nonce)`; the nonce stays in
  the memory of the process that asked. `delfin-agent messenger pair` refuses
  inside an agent's process tree; `/pair` in a service session arrives through
  the agent service socket, which refuses agent processes.
- The bridge (the agent-service user) collects every person's requests; the
  person is the directory a request was found in. The gateway asks the daemon
  for a one-time invitation (`/_connect`, random connection keys), records it
  bound to the person and the request (10 minutes, single use), and the bridge
  writes the link into the person's 0600 answer file. The session shows it as
  a live notice that is not written to the session's event log; the audit log
  gets a 12-character hash of the link only.
- A phone that connects through it is bound to that person, **unverified**,
  and receives a 6-digit code. Until the code is typed in the session that
  started `/pair` (the verify request must carry the nonce), the phone hears
  only "verify first"; nothing it sends is routed. Three wrong codes or 10
  minutes without one remove the contact. The SimpleX security code is shown
  in the session too, for comparison on the phone.
- On verification the person's other phones are told ("new phone connected"),
  the person's pairing notes (`/pair phones`) get a line, and the phone is
  invited to the person's forum and enabled sessions.
- A used, expired or superseded invitation binds nothing: a second phone
  through the same link, a late scan and a replayed request id are refused and
  audited. A new `/pair` supersedes the person's open invitation.
- Limits: 3 phones per person, 3 invitations per person per 10 minutes.

Revocation: `/pair revoke LABEL|all` in a session or `delfin-agent messenger
unpair LABEL|all`; `entkoppeln` / `unpair` on the phone itself; the admin
`sudo bash scripts/pi/pi_messenger.sh --people` and `--revoke=NAME` (`owner`
for the account owner). A revoked phone is removed from its groups (or the bot
leaves a group it cannot remove it from; the next sync makes a new one), its
contact is deleted, the person's other phones are told. Lost phone: revoke,
then pair the new one. A binding never changes in place: a phone moves to
another person only through revocation and a new pairing.

The admin path `sudo bash scripts/pi/pi_messenger.sh --pair` pairs the account
owner (the sessions outside `~/people`) and needs the code as well: it is asked
at that terminal, or later with `--verify` (asked for, never an argument, so
it stays out of the shell history). An existing single-person installation is
migrated on the first write: its paired contact becomes the owner's verified
phone, its chats the owner's; nothing is re-paired.

Every invitation, connection, verification, failed code, revocation, refused
message and `hol`/`weg` is recorded with a `messenger_*` kind
(`delfin/messenger/audit.py`) in the containment panel, the attention inbox
(refusals) and the gateway's hash-chained audit log. No link, code or security
code is written to a log; log lines are masked once more before they are
written. The gateway's store carries an HMAC under `store.key` (0600): a file
whose mapping was edited without the key, a corrupt file or a format-1 file
after the key exists refuses the whole store (nobody is served) and is
audited. The key lives beside the store, so this detects edits by anything
but the gateway's own user; it is not a boundary against that user.

**Forum: one per person.** Each person gets their own forum with their own
enabled sessions; phones of different people never share a group. Text in a
forum reaches only that person's enabled sessions (`@NAME text` one of them);
explicit agent announcements (`austausch on`) are relayed only between that
person's sessions. A shared forum would show one person's session results to
everybody and feed one person's agent output into another person's agents as
context; neither is wanted. Posts are labelled `<session name> (<person>)`.

**Sessions to the phone by name.** `sessions` lists the person's own running
sessions, numbered (`📱`: connected). `hol N` / `hol NAME` (`get`) switches one
on for the messenger, invites the person's phones to its chat and shows the
session a notice (`/sync status`: "enabled (from phone LABEL)"); it is
audited. A session switched off in the session (`/sync off`) cannot be brought
from the phone. `weg` / `off` in a session's chat takes it off again; `hol`
may bring it back. The display name is the session's title in the agent
service (`/sync NAME`, `name NAME` on the phone, the TUI's rename): one field,
at most 50 characters; the chat follows it. Every message in a session's chat
starts with `[name]`.

What a person of people mode cannot do from the phone: schedules and
workspace freedom profiles stay the account owner's (one store per account).

## Sessions and the phone

Type `/sync` in each desired TUI or dashboard session.
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
| `sessions` | Your running sessions, numbered (📱: on the phone) |
| `hol N` / `hol NAME` (`get`) | Bring a running session of yours to the phone |
| `weg` (`off`) | In a session's chat: take it off the phone |
| `@NAME text` | In the forum: to that one session |
| `handys` (`phones`) | Your paired phones |
| `entkoppeln` (`unpair`) | In the direct chat: this phone unpairs itself |
| `hilfe` (`help`) | The commands, briefly |
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
(`DELFIN_MESSENGER_WORKSPACES` in the bridge unit); for a person of people mode
in `~/people/<name>/DELFIN-Handy/NAME`, in that person's own agent service. The phone supplies the name
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

Only verified phones' received messages are accepted, including inside groups,
and only in chats of the phone's own person. Other members, unverified or
revoked phones, sent-message echoes and unbound rooms cannot command sessions.
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
