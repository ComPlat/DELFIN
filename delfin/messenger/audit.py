"""The messenger's security events: invitations, pairings, revocations, refusals.

Each event goes to the two existing records: the containment panel and
attention inbox (``delfin.agent.security_events``, blocked events) and the
0600 hash-chained audit log of the process's HOME (``delfin.agent.audit_log``;
the gateway's HOME is its state directory). A detail never carries an
invitation link, a confirmation code or a security code: callers pass the
12-character link tag (pairing.link_tag), and :func:`event` masks links and
long digit runs once more before anything is written.

The gateway runs from a venv with websockets only; an import that fails
there costs the record, never the message being handled.
"""
from __future__ import annotations

import sys

from .pairing import mask

#: Every kind the messenger records, with the panel's label. Registered in
#: delfin.agent.security_events._KINDS (a test pins the two equal).
KINDS = {
    "messenger_invitation_issued": "Messenger: one-time phone invitation issued",
    "messenger_invitation_refused": "Messenger: phone invitation refused (limit or rate)",
    "messenger_invitation_expired": "Messenger: phone invitation expired unused",
    "messenger_invitation_replayed": "Messenger: used or expired invitation presented again",
    "messenger_phone_connected": "Messenger: phone connected, not yet verified",
    "messenger_phone_verified": "Messenger: phone verified and bound to a person",
    "messenger_verification_failed": "Messenger: wrong phone confirmation code",
    "messenger_phone_revoked": "Messenger: phone unpaired",
    "messenger_unverified_message": "Messenger: message from a phone not yet verified",
    "messenger_cross_person_refused": "Messenger: phone addressed another person's chat",
    "messenger_store_refused": "Messenger: pairing store failed its integrity check",
    "messenger_session_attached": "Messenger: session brought to the phone",
    "messenger_session_detached": "Messenger: session taken off the phone",
}


def event(kind: str, detail: str, *, blocked: bool = False) -> None:
    """Record one messenger event. Never raises."""
    text = mask(detail)[:300]
    try:
        from delfin.agent import security_events
        security_events.record(kind, "messenger", text, blocked=blocked)
    except Exception as exc:  # noqa: BLE001 - a missing module costs the record only
        print(f"[messenger] security event not recorded: {type(exc).__name__}",
              file=sys.stderr, flush=True)
    try:
        from delfin.agent import audit_log
        audit_log.append({"tool": "messenger", "event": kind,
                          "decision": "denied" if blocked else "auto", "detail": text})
    except Exception as exc:  # noqa: BLE001
        print(f"[messenger] audit line not written: {type(exc).__name__}",
              file=sys.stderr, flush=True)
