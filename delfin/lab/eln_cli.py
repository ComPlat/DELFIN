"""``delfin-agent eln ...``: the person's side of an ELN write.

Reading (allowed everywhere): ``pending``, ``show ID``, ``collections``.

The person's action is ``approve ID``. What guards it, measured on this tree:
the process marker -- ``_person_only`` below and
:func:`delfin.lab.eln_write.grant` both refuse when
``DELFIN_AGENT_SESSION`` is set, and the tool side
(``delfin/lab/eln_tools.py``) files a request and then waits, so no tool
call sends a write. A second layer is in place: the command-text gate in
``delfin/agent/approval_guard.py`` refuses ``eln approve|deny`` and a Python
``eln_write.grant(...)`` call written as a command, as it refuses
``delfin-agent approvals approve``. The durable fix -- a grant no descendant of an agent process can
mint, rather than one that merely kept an environment variable -- belongs
with ``delfin/lab/approvals.py``, which carries the same single-layer shape
for lab actions (identical on ``origin/main``).

Why a command and not a tool. An agent that could approve its own ELN
write has no approval -- the two tools that matter here (``eln_draft``,
``eln_write``) file a request and then wait, so the only route from a
filed request to a sent one is this command, run by a person in their own
terminal. The signature underneath (:func:`delfin.lab.eln_write.grant`)
refuses the agent-session marker as well, so a shell that somehow reached
it still produces no valid grant.

The shape is ``delfin/lab/cli.py``'s, one file over, for the same reasons:
``--help`` lists what exists, every action prints JSON a person can read,
and a refusal exits non-zero with its reason on stderr.
"""

from __future__ import annotations

import argparse
import json
import sys


def _print(obj) -> None:
    print(json.dumps(obj, indent=1, default=str, sort_keys=True))


def _person_only(what: str) -> None:
    from delfin.agent.approval_guard import inside_agent_session
    if inside_agent_session():
        print(f"blocked: `eln {what}` is a person's action; this process was "
              "started by an agent session", file=sys.stderr)
        raise SystemExit(3)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="delfin-agent eln",
        description="Chemotion ELN: the person's approvals and readings")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("pending", help="Writes waiting for your approval")
    s = sub.add_parser("show", help="One request in full")
    s.add_argument("request_id")
    a = sub.add_parser("approve", help="Approve one write (person only)")
    a.add_argument("request_id")
    sub.add_parser("collections",
                   help="The collections the ELN lets this account read")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd == "pending":
        from . import eln_write
        rows = []
        for rec in eln_write.pending():
            rows.append({
                "id": rec.get("id"),
                "what": eln_write.describe(
                    rec.get("element", ""),
                    {"name": rec.get("name"),
                     "collection_id": rec.get("collection_id")}),
                "server": rec.get("base_url"),
                "expires_in_s": max(0, int(rec.get("expires", 0)
                                           - eln_write._now())),
                "digest": str(rec.get("digest", ""))[:12],
            })
        _print(rows)
        return 0
    if args.cmd == "show":
        from . import store
        rec = store.snapshot().get("eln_writes", {}).get(str(args.request_id))
        if rec is None:
            print(f"no ELN write request {args.request_id!r}", file=sys.stderr)
            return 1
        _print(rec)
        return 0
    if args.cmd == "approve":
        _person_only("approve")
        from . import eln_write
        try:
            rec = eln_write.grant(args.request_id)
        except eln_write.WriteRefused as exc:
            print(f"refused: {exc}", file=sys.stderr)
            return 1
        print(f"approved {rec['id']}: {rec['element']} "
              f"{rec.get('name')!r} (single use, digest "
              f"{str(rec.get('digest', ''))[:12]})")
        return 0
    if args.cmd == "collections":
        from . import eln_api, eln_tools
        try:
            client = eln_tools._client()
            _print(client.collections())
        except eln_api.ClientError as exc:
            print(f"refused: {exc}", file=sys.stderr)
            return 1
        return 0
    return 2
