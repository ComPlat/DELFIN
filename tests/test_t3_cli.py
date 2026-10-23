"""T3 phase 3 — the delivery-receipt surfaces reach the real CLI entry point.

The storage core (session_messages.status()/ls(), scoped to the calling
sender) is committed and green in test_t3_phase3.py. A sender must also be
able to ask what happened to its mail from the terminal: that is what
``delfin-agent messages ls [--from KEY]`` and ``delfin-agent messages
status <message_id> --from KEY`` are for.

These controls drive delfin.agent.cli.main(argv) — parser, dispatch, exit
code, stdout — exactly as the ``delfin-agent`` command runs it, against the
message core. The security boundary the operator set is a hard contract
here: both commands are scoped to the caller's own ``from`` key. A caller
that names no key, or a key that is not the sender, sees nothing — the CLI
must never list another session's or the operator's mailbox, and ``unknown``
must not reveal that a message under a foreign id exists.

They are RED on the current tree, where ``messages`` is not a registered
subcommand at all (build_parser() has no ``messages`` parser and
``_route_argv`` would fall back to a chat session). They go green once
.gate/t3_cli.patch is built by the operator.
"""

from __future__ import annotations

import pytest

from delfin.agent import cli
from delfin.agent import session_messages as M
from delfin.agent import session_presence as P


@pytest.fixture(autouse=True)
def _dirs(tmp_path, monkeypatch):
    monkeypatch.setattr(M, "_DIR", tmp_path / "inbox")
    monkeypatch.setattr(P, "_DIR", tmp_path / "presence")
    P._last_written.clear()
    P._git_cache.clear()


def _main(argv, capsys):
    rc = cli.main(argv)
    out = capsys.readouterr()
    return rc, out.out, out.err


def test_messages_ls_from_the_terminal_lists_the_senders_mail(capsys):
    """The real CLI lists the mail the calling sender knows about, with its
    receipt, through delfin.agent.cli.main()."""
    sent = M.send("nacht-s17", "do you have the hash?", from_key="nacht-s16")
    M.send("nacht-s18", "second", from_key="nacht-s16")
    rc, out, err = _main(["messages", "ls", "--from", "nacht-s16"], capsys)
    assert rc == 0, err
    assert sent["id"][:12] in out   # a message id the sender got is listed
    assert "nacht-s17" in out
    assert "nacht-s18" in out


def test_messages_status_from_the_terminal_reports_the_receipt(capsys):
    """status(message_id) through the CLI reports queued, then delivered once
    the recipient takes the mailbox."""
    sent = M.send("nacht-s17", "heads-up", from_key="nacht-s16")
    rc, out, err = _main(["messages", "status", sent["id"],
                          "--from", "nacht-s16"], capsys)
    assert rc == 0, err
    assert "queued" in out
    M.take("nacht-s17")
    rc, out, err = _main(["messages", "status", sent["id"],
                          "--from", "nacht-s16"], capsys)
    assert rc == 0, err
    assert "delivered" in out


def test_messages_status_is_unknown_for_a_foreign_or_unknown_id(capsys):
    """A caller never learns that a message under a foreign id exists: status
    reports unknown, and the CLI never prints the body or a hint."""
    b = M.send("nacht-s17", "B's secret", from_key="nacht-s18")
    rc, out, err = _main(["messages", "status", b["id"],
                          "--from", "nacht-s16"], capsys)
    assert rc == 0, err
    assert "unknown" in out
    assert "B's secret" not in out


def test_messages_ls_without_a_sender_sees_nothing(capsys):
    """The security boundary the operator set: a caller that names no key
    sees nothing — the CLI does not fall back to 'everyone's mail'."""
    M.send("nacht-s17", "secret", from_key="nacht-s16")
    rc, out, err = _main(["messages", "ls"], capsys)
    assert rc == 0, err
    assert "secret" not in out
    assert "nacht-s17" not in out


def test_messages_is_a_registered_subcommand():
    """The parser knows `messages` — it must not be routed as a chat prompt."""
    from delfin.agent import cli as agent_cli
    parser = agent_cli.build_parser()
    assert "messages" in agent_cli._subcommand_names(parser)
