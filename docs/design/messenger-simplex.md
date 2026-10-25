# SimpleX rooms for DELFIN

The phone has a **DELFIN Forum** and a **separate private group per session**.
Each group contains the bot and the paired human. DELFIN labels every outgoing
message with the session title and stable session ID; the agents do not need
separate SimpleX accounts. A group ID is permanently bound to a session ID.
Renaming a session never changes its recipient, and a missing session never
falls back to another agent. The original direct contact is for setup only.

## What reaches the phone

Short results, failures, questions and pending approval previews. No transcript,
model reasoning, tool log or file attachment is forwarded. `status` reads a short
service summary without asking the model. `update` sends a request for a brief
progress report into the current turn. Reports are redacted and length-limited.

## Phone commands

- `/sessions`: session IDs, titles and state.
- `/sync`: create the forum and individual chats for sessions on this service.
  Accept the group invitations in SimpleX. Repeat to include new sessions.
- `/chat SESSION_ID`: create the chat for one existing local session.
- `/forum`: create just the shared forum.
- In an individual chat, ordinary text goes only to that session.
- In the forum, ordinary text goes to every connected session, once each.
  The confirmation lists each recipient, including failures. Agent output is
  never automatically sent as a prompt to other agents.
- `/status`, `/update`, `/pending`: one session in its chat; all connected
  sessions in the forum, with labels.
- `/stop`: cancel the current turn in the chat's session; in the forum, cancel
  the current turns of connected sessions. This does not shut down the service.
- `/answer CODE text`: answer the exact question in its individual chat.
- `/approve CODE`, `/deny CODE`: decide the exact request in its individual chat.
  There is no collective approval command in the forum.

Free text never implicitly answers a question or grants an approval. `@name`
redirection is rejected: use the intended session's chat.

## Approvals

The service's `session.pending` endpoint reads live requests without attaching,
waking a worker or fetching its history. The bridge issues a random code bound
to the session, request ID and complete request description. It expires after
10 minutes and is consumed before a decision is forwarded. Replaced, resolved,
wrong-session and previously used requests are refused. Bridge restarts
invalidate codes: `/pending` obtains fresh ones.

Decisions go through the existing `permission.reply` and ProtocolBroker. A
confirmation may be allowed once or denied; a plan may be approved into default
mode or rejected. Permanent allow rules, accept-edits mode and blanket approvals
are not exposed. A preview too long for the phone cannot be approved there;
use the dashboard for the complete action. Questions use the broker's text
answer path, never its approval choice path.

## Setup

Install DELFIN with `pip install ".[agent,messenger]"` into
`/opt/delfin-messenger/bin/venv` and provision with `scripts/pi/pi_messenger.sh`.
The agent service must run on the same host as user `delfin`, or set
`MESSENGER_SERVICE_USER` and optionally `MESSENGER_BRIDGE_PYTHON` when provisioning.
Restart the agent service after upgrading: the bridge needs `session.pending`.

The provisioner gives daemon and gateway port 5227. The bridge unit runs as the
agent-service user with group `delfin-gw`; its 0660 socket lives in the separate
0750 directory `/run/delfin-messenger-bridge`. The service's own runtime directory
stays private. The bridge uses the host's `/tmp` to reach the service socket.

Create a SimpleX bot profile in the daemon database before pairing. Then:

```sh
sudo -u delfin-gw env DELFIN_MESSENGER_DIR=/var/lib/delfin-messenger/pairing /opt/delfin-messenger/bin/venv/bin/python -m delfin.messenger pair
```

Open the one-time link in the phone app and compare the connection security code.
Only the connection belonging to that invitation may complete pairing. Send
`/sync` to the paired bot to get the forum and session invitations. Alternatively,
use `python -m delfin.messenger forum` or `connect SESSION_ID` with the same user
and environment. There is no dashboard/TUI connection button yet; these commands
also work from a terminal.

`revoke` with the same environment immediately removes routing authorization;
commands already dispatched are not undone. `status` lists the stored chat
bindings. Old single-contact pairing files remain valid, but free text in the
old direct chat no longer dispatches tasks: create session rooms first.

## Delivery and validation

Only messages from the paired contact are accepted, including inside groups.
Sent-message echoes, other group members and unbound groups cannot command
agents. Rate/size limits and persistent message IDs prevent replay after a
restart. Dispatch is at most once: after a connection failure, check `/status`
before repeating a command; partially delivered broadcasts are not retried
silently. Revocation is checked before each broadcast recipient and notification.

SimpleX v7.0.3 uses textual commands inside the WebSocket `cmd` field. Group,
message and invitation formats follow the pinned upstream references:
[commands](https://github.com/simplex-chat/simplex-chat/blob/v7.0.3/bots/api/COMMANDS.md),
[types](https://github.com/simplex-chat/simplex-chat/blob/v7.0.3/bots/api/TYPES.md).

Automated tests use simulated SimpleX traffic, actual local sockets and a real
DELFIN service/approval broker with a fake model. A live SimpleX daemon and phone
acceptance test remain necessary; provisioning and deployment are not performed
by these tests.
