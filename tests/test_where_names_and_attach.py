"""/`delfin-agent where` names every running session and says how to rejoin it.

The person's rule (2026-11-12) and wave55/e2: ``delfin-agent where`` must, for
every running session, print the session's name (the title the dashboard and
the TUI both show — the same name everywhere), the exact one-line terminal
command to attach (``delfin-agent tui --resume <id>``), and the dashboard link
when the dashboard has a kept page for it.

Before this change ``where`` printed one line per open session that showed only
a host-mangled ``key``, its branch and "seen Nm ago" — no name, no attach
command, no dashboard link. Finding a session again took several commands.
"""

from __future__ import annotations

import time

from delfin.agent import where


def _row(**over):
    row = {
        "host": "login02",
        "key": "login02-79b8abcd",
        "session_id": "0f3a9c21deadbeef",
        "title": "Butan Konformere",
        "workspace": "/home/me/proj",
        "branch": "swarm/main",
        "updated_at": time.time() - 30,
    }
    row.update(over)
    return row


def test_a_running_session_is_named_by_its_title():
    out = where.format_text({}, [_row()])
    assert "Butan Konformere" in out
    # the name shown is the title, not the host-mangled key
    assert "'login02-79b8abcd'" not in out


def test_where_prints_the_one_line_attach_command():
    out = where.format_text({}, [_row()])
    assert "delfin-agent tui --resume 0f3a9c21deadbeef" in out


def test_no_attach_command_when_the_session_has_no_id():
    out = where.format_text({},
                            [_row(session_id="", title="", key="login02-raw")])
    assert "delfin-agent tui --resume" not in out
    # the key is still the fallback name when there is no title
    assert "login02-raw" in out


def test_where_prints_the_dashboard_link_for_a_kept_session():
    record = {"host": "login02", "port": 8867, "resume_path": "/delfin",
              "running": True}
    from unittest import mock
    with mock.patch.object(where, "kept_session_names",
                           return_value=["Butan Konformere"]):
        out = where.format_text(record, [_row()])
    # the running session's own line carries the dashboard link
    assert "http://localhost:8867/delfin?session=Butan%20Konformere" in out


def test_no_session_dashboard_link_when_the_session_is_not_kept():
    record = {"host": "login02", "port": 8867, "resume_path": "/delfin",
              "running": True}
    from unittest import mock
    with mock.patch.object(where, "kept_session_names", return_value=["Other"]):
        out = where.format_text(record, [_row()])
    # the running session still lists, named and with its attach command,
    # but carries no dashboard link of its own
    assert "Butan Konformere" in out
    assert "delfin-agent tui --resume 0f3a9c21deadbeef" in out
    assert "session=Butan" not in out
