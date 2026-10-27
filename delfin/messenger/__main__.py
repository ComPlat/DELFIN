"""`python -m delfin.messenger pair|run|revoke|status`.

    pair     show the connection invitation, wait for the phone to scan
             it, save the new contact as the one allowed sender
    run      serve: daemon events in, gateway decisions, answers out
    revoke   unpair (the phone may no longer command the agent)
    status   whether paired, whether the bridge is reachable
"""

from __future__ import annotations

import argparse
import asyncio
import sys

from . import simplex, store
from .bridge_client import BridgeClient, BridgeDown, default_socket_path
from .gateway import Rooms, bot_name
from .lang import t

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


def _pair() -> int:
    contact_id = store.paired_contact()
    if contact_id:
        print(f"already paired (contact {contact_id}); revoke first to re-pair")
        return 1

    async def run() -> int:
        async with simplex.Daemon(wait=DAEMON_WAIT_S) as daemon:
            await daemon.subscribe()
            invitation = await daemon.pending_invitation()
            if not invitation or getattr(daemon, "invitation_id", None) is None:
                print("No valid pairing invitation received.", file=sys.stderr)
                return 1
            # The link last: a phone's terminal shows the end of the output,
            # and a QR code is taller than its screen.
            _show_qr(invitation)
            print("Scan the code with SimpleX, or copy this link into SimpleX "
                  "(+ > paste link):")
            print(invitation)
            print("waiting for the phone to connect ...")
            async for event in daemon.events():
                new_id = simplex.new_contact_id(event)
                if new_id:
                    connection = (event.get("contact") or {}).get("activeConn") or {}
                    if connection.get("connId") != getattr(daemon, "invitation_id", None):
                        continue
                    store.save(new_id)
                    print(f"paired: contact {new_id}")
                    code = await daemon.contact_code(new_id)
                    if code:
                        print("security code (on the phone: the contact > "
                              f"Verify security code; both must match):\n  {code}")
                    else:
                        print("security code: compare it on the phone (the contact > "
                              "Verify security code); this daemon did not report it")
                    await daemon.send_text(new_id, t("Verbunden mit {bot} \u2713").format(bot=bot_name()))
                    return 0
        return 1

    return asyncio.run(run())


def _run() -> int:
    contact_id = store.paired_contact()
    if not contact_id:
        print("not paired; run `python -m delfin.messenger pair` first",
              file=sys.stderr)
        return 1

    async def run() -> int:
        async with simplex.Daemon(wait=DAEMON_WAIT_S) as daemon:
            print(f"serving as the phone's agent (contact {contact_id})")
            await Rooms(daemon).run()
        return 0

    try:
        return asyncio.run(run())
    except KeyboardInterrupt:
        return 0


def _revoke() -> int:
    if store.revoke():
        print("unpaired")
        return 0
    print("was not paired")
    return 1


def _status() -> int:
    contact_id = store.paired_contact()
    print(f"paired: {contact_id or 'no'}")
    for chat, room in store.rooms().items():
        print(f"{chat}: {room['title']} -> {room['session_id'] or 'ALL CONNECTED SESSIONS'}")
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
            rooms = Rooms(daemon)
            if session_id:
                await rooms.ensure_room(session_id)
            else:
                await rooms.sync()
        print("Accept the SimpleX group invitations: forum + one chat per session.")
        return 0
    return asyncio.run(run())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="delfin.messenger",
        description="the SimpleX gateway for the DELFIN agents")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("pair", help="show the invitation, pair the phone")
    sub.add_parser("run", help="serve the paired phone")
    sub.add_parser("revoke", help="unpair")
    sub.add_parser("status", help="pairing, rooms and bridge state")
    sub.add_parser("forum", help="create forum and chats for current sessions")
    connect = sub.add_parser("connect", help="create a fixed chat for one session")
    connect.add_argument("session_id")
    args = parser.parse_args(argv)
    if args.command in {"forum", "connect"}:
        return _setup(getattr(args, "session_id", None))
    return {"pair": _pair, "run": _run,
            "revoke": _revoke, "status": _status}[args.command]()


if __name__ == "__main__":
    sys.exit(main())
