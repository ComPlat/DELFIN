"""`python -m delfin.messenger pair|verify|run|revoke|people|status`.

Run as the gateway's user (``delfin-gw``; pi_messenger.sh does it):

    pair         the admin path: pair the account OWNER's phone. Shows the
                 invitation, waits for the phone, sends it a 6-digit code
                 and asks for it here; until it is typed the phone is
                 bound but may do nothing (``verify CODE`` later)
    verify [CODE]
                 confirm the owner's phone that waits for its code (asked
                 for when left out; pi_messenger.sh --verify)
    run          serve: daemon events in, gateway decisions, answers out
    revoke [P]   unpair person P's phones; ``revoke --owner`` the account
                 owner's; neither: everybody
    people       every person's phones: label, state, since (no codes)
    status       pairing, rooms, bridge reachability

People pair their own phones from their sessions (``/pair``,
delfin.messenger.pairing); this module is the admin's side.
"""

from __future__ import annotations

import argparse
import asyncio
import sys
import time

from . import audit, pairing, simplex, store
from .bridge_client import BridgeClient, BridgeDown, default_socket_path
from .gateway import Rooms, bot_name, person_name, _PERSON
from .lang import t
from .people import OWNER, valid as valid_person

#: How long a just started daemon is given to listen. Pairing right after
#: the unit's start was refused on a Pi 5 (ConnectionRefusedError).
DAEMON_WAIT_S = 60.0


def _show_qr(invitation: str) -> None:
    """Draw the invitation for the phone's scanner; the link stays printed."""
    import shutil
    import subprocess
    tool = shutil.which("qrencode")
    if tool is None:
        print("(a QR code needs qrencode: apt install qrencode)")
        return
    sys.stdout.flush()
    subprocess.run([tool, "-t", "ansiutf8", invitation], check=False, timeout=10)


def _ask_code() -> str | None:
    """The code typed at this terminal, or None without one (no terminal)."""
    try:
        return input("Type the 6-digit code DELFIN sent to the phone: ")
    except (EOFError, OSError, KeyboardInterrupt):
        return None


def _pair() -> int:
    if store.refused():
        print("the pairing store failed its integrity check; see the audit log", file=sys.stderr)
        return 1
    contact_id = store.paired_contact()
    if contact_id:
        print(f"already paired (contact {contact_id}); revoke first to re-pair")
        return 1

    async def run() -> int:
        async with simplex.Daemon(wait=DAEMON_WAIT_S) as daemon:
            await daemon.subscribe()
            invitation = await daemon.pending_invitation()
            conn = getattr(daemon, "invitation_id", None)
            if not invitation or conn is None:
                print("No valid pairing invitation received.", file=sys.stderr)
                return 1
            try:
                store.add_invitation(str(conn), person=OWNER, request="admin", by="admin",
                                     tag=pairing.link_tag(invitation))
            except store.PairingRefused as exc:
                print(f"pairing refused: {exc.reason}", file=sys.stderr)
                return 1
            audit.event("messenger_invitation_issued",
                        f"owner: admin invitation {pairing.link_tag(invitation)}")
            # The link last: a phone's terminal shows the end of the output,
            # and a QR code is taller than its screen.
            _show_qr(invitation)
            print("Scan the code with SimpleX, or copy this link into SimpleX "
                  "(+ > paste link):")
            print(invitation)
            print("waiting for the phone to connect ...")
            async for event in daemon.events():
                new_id = simplex.new_contact_id(event)
                if not new_id:
                    continue
                connection = (event.get("contact") or {}).get("activeConn") or {}
                if str(connection.get("connId")) != str(conn):
                    continue
                inv, why = store.consume_invitation(str(conn), by="admin")
                if inv is None:
                    print(f"invitation {why}; not paired", file=sys.stderr)
                    return 1
                code = pairing.new_code()
                from .activation import clean_name
                name = clean_name(simplex.contact_name_of(event.get("contact") or {}))
                store.bind(new_id, person=OWNER, label=name or f"phone-{new_id}",
                           request="admin", code=code)
                audit.event("messenger_phone_connected",
                            f"contact {new_id} connected for the owner; waits for its code")
                print(f"connected: contact {new_id} (not yet allowed to do anything)")
                security = await daemon.contact_code(new_id)
                if security:
                    print("security code (on the phone: the contact > "
                          f"Verify security code; both must match):\n  {security}")
                else:
                    print("security code: compare it on the phone (the contact > "
                          "Verify security code); this daemon did not report it")
                await daemon.send_text(new_id, t(
                    'DELFIN-Kopplung: Code {code}. Am Rechner eingeben, wo die Kopplung gestartet wurde '
                    '(in der Sitzung: /pair {code}). Gilt 10 Minuten; bis dahin wird hier nichts ausgeführt.').format(code=code))
                for _ in range(pairing.MAX_ATTEMPTS):
                    typed = _ask_code()
                    if typed is None:
                        print("not verified yet: within 10 minutes run "
                              "sudo bash pi_messenger.sh --verify and type the code")
                        return 0
                    result, _rec = store.check_code(new_id, typed, request="admin")
                    if result == "verified":
                        return await _welcome_owner(daemon, new_id)
                    print({"wrong": "wrong code", "failed": "too many wrong codes: removed",
                           "expired": "too late: removed", "none": "nothing waits"}.get(result, result))
                    audit.event("messenger_verification_failed",
                                f"owner phone {new_id}: {result}", blocked=True)
                    if result != "wrong":
                        await _quiet(daemon.delete_contact(new_id))
                        return 1
                return 1
        return 1

    return asyncio.run(run())


async def _quiet(awaitable):
    try:
        await awaitable
    except Exception:  # noqa: BLE001 - best effort
        pass


async def _welcome_owner(daemon, cid: str) -> int:
    audit.event("messenger_phone_verified", f"contact {cid} verified for the owner (admin path)")
    print(f"paired: contact {cid}")
    await daemon.send_text(cid, t("Verbunden mit {bot} ✓").format(bot=bot_name()))
    # The forum for all agents from the start, through the gateway's own
    # path: only the owner's verified phones are invited.
    token = _PERSON.set(OWNER)
    try:
        await Rooms(daemon).ensure_room()
        print("forum created: accept its invitation on the phone")
    except Exception as exc:
        print(f"forum not created now ({type(exc).__name__}); "
              "the gateway creates it with the first /sync")
    finally:
        _PERSON.reset(token)
    return 0


def _verify(code: str | None) -> int:
    if code is None:
        code = _ask_code()
        if code is None:
            print("no code typed", file=sys.stderr)
            return 1
    cid = next((c for c, r in store.phones(OWNER) if r.get("state") == "unverified"
                and r.get("request") == "admin"), None)
    if cid is None:
        print("no owner phone waits for its code", file=sys.stderr)
        return 1
    result, _rec = store.check_code(cid, code, request="admin")
    if result != "verified":
        audit.event("messenger_verification_failed", f"owner phone {cid}: {result}", blocked=True)
        print(result, file=sys.stderr)
        return 1
    audit.event("messenger_phone_verified", f"contact {cid} verified for the owner (admin path)")
    print(f"paired: contact {cid}; the gateway greets it and opens the forum")
    return 0


def _run() -> int:
    async def run() -> int:
        async with simplex.Daemon(wait=DAEMON_WAIT_S) as daemon:
            print(f"serving {len(store.persons())} person(s) with a verified phone")
            await Rooms(daemon).run()
        return 0

    try:
        return asyncio.run(run())
    except KeyboardInterrupt:
        return 0


def _revoke(person: str | None = None, *, owner: bool = False) -> int:
    if owner:
        person = OWNER
    elif person is None:
        if store.revoke():
            audit.event("messenger_phone_revoked", "every phone unpaired by the admin")
            print("unpaired")
            return 0
        print("was not paired")
        return 1
    who = person
    if not valid_person(who):
        print(f"not a person name: {person!r}", file=sys.stderr)
        return 2
    gone = store.revoke_person(who, by="the admin")
    for entry in gone:
        audit.event("messenger_phone_revoked",
                    f"contact {entry['contact']} of {person_name(who)} unpaired by the admin")
    print(f"unpaired {len(gone)} phone(s) of {who or 'the owner'}; "
          "the gateway ends their group memberships")
    return 0 if gone else 1


def _people() -> int:
    if store.refused():
        print("the pairing store failed its integrity check", file=sys.stderr)
        return 1
    rows = sorted(store.contacts().items(), key=lambda kv: (kv[1].get("person") != OWNER,
                                                            kv[1].get("person"), kv[0]))
    if not rows:
        print("no phone paired")
    for cid, rec in rows:
        since = time.strftime("%Y-%m-%d %H:%M", time.localtime(float(rec.get("since") or 0)))
        print(f"{rec.get('person') or 'owner':<16} {rec.get('label') or '-':<24} "
              f"{rec.get('state'):<10} since {since}  contact {cid}")
    return 0


def _status() -> int:
    if store.refused():
        print("pairing store: REFUSED (integrity check failed)")
        return 1
    print(f"paired: {store.paired_contact() or 'no'} (owner); "
          f"persons with a verified phone: {len(store.persons())}")
    for chat, room in store.rooms().items():
        print(f"{chat}: {room['title']} -> {room['session_id'] or 'ALL CONNECTED SESSIONS'}"
              f" ({room.get('person') or 'owner'})")
    path = default_socket_path()
    try:
        BridgeClient(path).ask({"cmd": "status"})
        print(f"bridge: reachable at {path}")
    except (BridgeDown, OSError) as exc:
        print(f"bridge: not reachable ({exc})")
        return 1
    return 0


def _setup(session_id: str | None = None) -> int:
    if not store.paired_contact():
        print("Pair the phone first.", file=sys.stderr)
        return 1
    async def run():
        async with simplex.Daemon(wait=DAEMON_WAIT_S) as daemon:
            token = _PERSON.set(OWNER)
            try:
                rooms = Rooms(daemon)
                if session_id:
                    await rooms.ensure_room(session_id)
                else:
                    await rooms.sync()
            finally:
                _PERSON.reset(token)
        print("Accept the SimpleX group invitations: forum + one chat per session.")
        return 0
    return asyncio.run(run())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="delfin.messenger",
        description="the SimpleX gateway for the DELFIN agents")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("pair", help="pair the account owner's phone (admin)")
    verify = sub.add_parser("verify", help="confirm the owner's phone with its code")
    # Asked at the terminal when left out: a code on the command line
    # would stay in the shell's history.
    verify.add_argument("code", nargs="?")
    sub.add_parser("run", help="serve the paired phones")
    revoke = sub.add_parser("revoke", help="unpair a person's phones (--owner: the account "
                                           "owner's; no name: everybody)")
    revoke.add_argument("person", nargs="?")
    revoke.add_argument("--owner", action="store_true",
                        help="the account owner's phones (a person may be called owner)")
    sub.add_parser("people", help="every person's phones")
    sub.add_parser("status", help="pairing, rooms and bridge state")
    sub.add_parser("forum", help="create forum and chats for current sessions")
    connect = sub.add_parser("connect", help="create a fixed chat for one session")
    connect.add_argument("session_id")
    args = parser.parse_args(argv)
    if args.command in {"forum", "connect"}:
        return _setup(getattr(args, "session_id", None))
    if args.command == "verify":
        return _verify(args.code)
    if args.command == "revoke":
        if args.owner and args.person:
            parser.error("revoke: a person name or --owner, not both")
        return _revoke(args.person, owner=args.owner)
    return {"pair": _pair, "run": _run, "people": _people,
            "status": _status}[args.command]()


if __name__ == "__main__":
    sys.exit(main())
