"""Drafting an email, deliberately without sending it.

Sending is irreversible, happens under the user's own identity, and has
no read-back — and reading back what was written is the discipline that
makes the rest of the office layer trustworthy. It also runs on a shared
academic machine, where the credential that would authorise it
authenticates the human to everything, not to one API.

So this writes a .eml the user opens and sends themselves. Their mail
client is a better final confirmation than any preview we could build: it
shows the real rendering, the real sender, and the send button is theirs.
No credential, no network, no new trust — and the file goes through the
same write gate as every other document.

The recipient is validated here rather than left to the client. An
address the client silently drops is worse than one it rejects, and the
whole point of a draft is that it opens already correct.
"""

from __future__ import annotations

import email
import json
import pathlib
import tempfile

import pytest

from delfin.agent import api_client as A
from delfin.agent.office import OfficeError, draft_email


@pytest.fixture
def folder():
    d = pathlib.Path(tempfile.mkdtemp(prefix="office_"))
    (d / "anhang.pdf").write_bytes(b"%PDF-1.4 stub")
    return d


def _perms(folder, mode="acceptEdits"):
    perms = A.KitToolPermissions(workspace=folder, agent_role="office_agent")
    perms.mode = mode
    perms.confirm_callback = lambda n, a, prev: True
    return perms


def _run(folder, args, mode="acceptEdits"):
    executor = A._DocToolExecutor.__new__(A._DocToolExecutor)
    return json.loads(executor.execute("draft_email", args,
                                       _perms(folder, mode)))


# ---------------------------------------------------------------------------
# The file
# ---------------------------------------------------------------------------

def test_a_draft_is_a_file_the_user_opens(folder):
    out = _run(folder, {"path": "anfrage.eml", "to": "verwaltung@kit.edu",
                        "subject": "Reisekostenantrag Juni",
                        "body": "Guten Tag,\n\nanbei der Antrag."})
    assert (folder / "anfrage.eml").exists()
    assert out["to"] == ["verwaltung@kit.edu"]
    assert "NOT sent" in out["note"], (
        "the result must say the send is still the user's to make")


def test_the_file_opens_as_a_real_message(folder):
    _run(folder, {"path": "a.eml", "to": "a@kit.edu", "cc": "b@kit.edu",
                  "subject": "Betreff", "body": "Viele Grüße"})
    msg = email.message_from_bytes((folder / "a.eml").read_bytes())
    assert msg["Subject"] == "Betreff"
    assert msg["To"] == "a@kit.edu"
    assert msg["Cc"] == "b@kit.edu"
    assert msg["Message-ID"] and msg["Date"]


def test_umlauts_survive_the_round_trip(folder):
    """German administrative text is the normal case here, not an edge one."""
    _run(folder, {"path": "u.eml", "to": "a@kit.edu", "subject": "Grüße",
                  "body": "Müller, Straße, Öl"})
    # message_from_bytes yields a legacy Message; the modern policy is what
    # gives get_content(). Parse with it rather than reaching for an
    # attribute the default parser does not have.
    from email import policy

    msg = email.message_from_bytes((folder / "u.eml").read_bytes(),
                                   policy=policy.default)
    assert "Müller" in msg.get_content()
    assert "Grüße" in msg["Subject"]


def test_several_recipients_are_split(folder):
    out = _run(folder, {"path": "m.eml", "to": "a@kit.edu, b@kit.edu",
                        "subject": "s", "body": "b"})
    assert out["to"] == ["a@kit.edu", "b@kit.edu"]


# ---------------------------------------------------------------------------
# Attachments come from the folder
# ---------------------------------------------------------------------------

def test_a_relative_attachment_is_found(folder):
    """The agent names files the way it sees the folder. Resolved against
    the process directory instead, they are not there -- and the writer
    skipped them silently, so the draft went out without the document it
    was about."""
    out = _run(folder, {"path": "a.eml", "to": "a@kit.edu", "subject": "s",
                        "body": "b", "attachments": ["anhang.pdf"]})
    assert out["attachments"] == ["anhang.pdf"]
    assert out["attachments_skipped"] == []


def test_an_attachment_outside_the_folder_is_refused(folder):
    out = _run(folder, {"path": "a.eml", "to": "a@kit.edu", "subject": "s",
                        "body": "b", "attachments": ["/etc/hostname"]})
    assert "error" in out


def test_a_missing_attachment_is_reported_not_hidden(folder):
    out = _run(folder, {"path": "a.eml", "to": "a@kit.edu", "subject": "s",
                        "body": "b", "attachments": ["gibtsnicht.pdf"]})
    assert out.get("attachments_skipped"), (
        "a dropped attachment must be named; silence reads as attached")


# ---------------------------------------------------------------------------
# The same gates as every other document write
# ---------------------------------------------------------------------------

def test_writing_outside_the_folder_is_refused(folder):
    assert "error" in _run(folder, {"path": "/etc/leak.eml", "to": "a@b.de",
                                    "subject": "s", "body": "b"})


def test_plan_mode_refuses_it(folder):
    assert "error" in _run(folder, {"path": "a.eml", "to": "a@b.de",
                                    "subject": "s", "body": "b"}, mode="plan")


def test_an_existing_draft_is_not_overwritten(folder):
    args = {"path": "a.eml", "to": "a@b.de", "subject": "s", "body": "b"}
    _run(folder, args)
    assert "error" in _run(folder, args)


# ---------------------------------------------------------------------------
# The recipient is checked here, not left to the mail client
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", ["kein-at-zeichen", "a@b", "@kit.edu", "a b@k.de"])
def test_a_malformed_address_is_refused(folder, bad):
    with pytest.raises(OfficeError):
        draft_email(folder / "x.eml", to=bad, subject="s", body="b")


def test_a_draft_needs_a_recipient_and_a_subject(folder):
    with pytest.raises(OfficeError):
        draft_email(folder / "x.eml", to=None, subject="s", body="b")
    with pytest.raises(OfficeError):
        draft_email(folder / "y.eml", to="a@b.de", subject="  ", body="b")


def test_only_eml_is_accepted(folder):
    with pytest.raises(OfficeError):
        draft_email(folder / "x.txt", to="a@b.de", subject="s", body="b")


# ---------------------------------------------------------------------------
# Nothing here sends
# ---------------------------------------------------------------------------

def test_the_office_layer_has_no_send_path():
    """The capability is deliberately absent, not merely unused: an
    irreversible outward action with no read-back does not belong behind a
    tool the model can reach."""
    import delfin.agent.office as office

    source = pathlib.Path(office.__file__).read_text(encoding="utf-8")
    for forbidden in ("smtplib", "SMTP(", "sendmail", "starttls"):
        assert forbidden not in source, forbidden


def test_the_tool_description_says_it_does_not_send():
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI

    tool = next(t["function"] for t in _DOC_TOOLS_OPENAI
                if t["function"]["name"] == "draft_email")
    assert "Never sends" in tool["description"]


def test_office_may_use_it_and_the_surface_still_fits():
    from delfin.agent.api_client import (_tool_denied_for_role,
                                         tool_schema_token_report)

    assert not _tool_denied_for_role("office_agent", "draft_email")
    # 9_125 -> 9_133. Column paging (`start_col` on read_document) is the
    # capability that raised it, and it is paid for as far as it can be:
    # read_document's own description and two of its parameter texts were
    # tightened, returning 12 of the 20 tokens the new parameter costs.
    # The remaining 8 are the measured price of the change, so the ceiling
    # moves by 8 and not by a round number.
    #
    # It buys the remedy for a limit that previously had none: the reader
    # said "showing 40 of 87 columns" and there was no way to reach the
    # other 47 — the slice always began at column 1.
    #
    # 9_133 -> 9_287: sum_column, measured at 154 tokens after its texts
    # were cut to the shortest form that still says what it refuses and
    # how to get past the refusal. The reason it is worth a raise is
    # written next to the canonical ceiling in test_tool_schema_budget.py;
    # this assertion tracks it so the office surface is checked from the
    # email side too.
    assert tool_schema_token_report()["total_tokens"] <= 9_287
