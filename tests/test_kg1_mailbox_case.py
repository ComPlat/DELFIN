"""A message reaches the session no matter how its address was cased.

Wave 16: an idle ``delfin-agent chat`` session announced itself as
``welle-16`` (the name it was started with), but a sender handed the
address ``Welle-16`` wrote ``Welle-16.jsonl`` while the idle prompt read
``welle-16.jsonl`` — the message sat there and was never taken, silently.

``_reserved_canonical`` (session_messages.py) folds only the RESERVED
addresses ("operator") to one mailbox; every other key keeps the case it
was handed, and ``_known_key`` compares the presence record's key
case-sensitively. A case variant of a real session's name must map to the
SAME mailbox file a read by that session opens, and must be deliverable —
a true typo (a name nobody announced) stays refused.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_messages as msgs


class _Theme:
    def dim(self, t): return t
    def bold(self, t): return t
    def red(self, t): return t


class _Transcript:
    def __init__(self): self.lines = []
    theme = _Theme()
    def chrome(self, line): self.lines.append(line)


class _Opts:
    def __init__(self, name=""): self.session_name = name


def _agent(name="welle-16", session_id="abcdef0123456789"):
    from delfin.agent.repl import TerminalAgent
    a = object.__new__(TerminalAgent)
    a.opts = _Opts(name)
    a.engine = type("E", (), {"session_id": session_id})()
    a.transcript = _Transcript()
    return a


@pytest.fixture
def inbox(tmp_path, monkeypatch):
    monkeypatch.setattr(msgs, "_DIR", tmp_path / "session_inbox")
    return tmp_path


class TestOneMailboxPerName:
    def test_case_variants_share_one_file(self, inbox):
        assert msgs._inbox("Welle-16") == msgs._inbox("welle-16")

    def test_the_exact_name_still_writes_its_own_file(self, inbox):
        assert msgs._inbox("welle-16") != msgs._inbox("welle-17")

    def test_case_variant_send_is_taken_by_the_idle_prompt(self, inbox):
        # The session announces as "welle-16"; the sender writes "Welle-16".
        msgs.send("Welle-16", "wake: your partner pushed")
        got = _agent(name="welle-16")._operator_messages("")
        assert "wake: your partner pushed" in got

    def test_the_exact_name_is_taken_as_before(self, inbox):
        msgs.send("welle-16", "still delivered")
        got = _agent(name="welle-16")._operator_messages("")
        assert "still delivered" in got

    def test_a_true_typo_is_not_a_case_variant(self, inbox):
        assert msgs._inbox("welle-16") != msgs._inbox("welle-16x")


class TestDeliverability:
    def test_a_case_variant_of_a_known_session_is_deliverable(self, inbox):
        msgs.send("welle-16", "first contact")   # makes the inbox known
        assert msgs.deliverable("Welle-16") is True

    def test_an_unknown_name_stays_refused(self, inbox):
        assert msgs.deliverable("no-such-session") is False

    def test_a_presence_record_is_identity_for_its_EXACT_name(self, inbox,
                                                              monkeypatch):
        # Old decision (phase 2): the record's key was compared through the
        # case-folding _canonical_key, so the upper-case variant of an
        # announced name was deliverable. That broke the earlier t3
        # adversarial pin (tests/test_t3_adversarial_phase2.py:
        # test_known_session_case_variant_not_deliverable): a sender must
        # not reach nacht-s12 by spelling NACHT-S12. New decision: identity
        # is the record's EXACT key; only the mailbox (send/take) folds.
        import delfin.agent.session_presence as P
        monkeypatch.setattr(P, "_DIR", inbox / "session_presence")
        P.announce("welle-16", session_id="abcdef01", title="W16",
                   workspace=".")
        assert msgs.deliverable("welle-16") is True
        assert msgs.deliverable("Welle-16") is False
        # The fold still serves delivery: a send to the variant lands in
        # the ONE mailbox the exact-name session reads.
        msgs.send("Welle-16", "wake: folded mailbox")
        got = _agent(name="welle-16")._operator_messages("")
        assert "wake: folded mailbox" in got

    def test_reserved_addresses_stay_case_insensitive(self, inbox):
        assert msgs.deliverable("OPERATOR") is True

    def test_a_case_variant_of_an_ANNOUNCED_name_is_not_deliverable(self,
                                                                    inbox,
                                                                    monkeypatch):
        # The t3 adversarial pin (tests/test_t3_adversarial_phase2.py,
        # test_known_session_case_variant_not_deliverable): a case variant
        # of a known session NAME must not count as deliverable. A sender
        # should not reach nacht-s12 by spelling it NACHT-S12, and a
        # mailbox keyed by the variant must never be taken by the real
        # session. The case fold serves only the MAILBOX (send/take land
        # in one file), never identity: without an inbox or presence
        # record under the variant, deliverable("NACHT-S12") stays False.
        import delfin.agent.session_presence as P
        monkeypatch.setattr(P, "_DIR", inbox / "session_presence")
        P.announce("nacht-s12", session_id="abcdef02", title="S12",
                   workspace=".")
        assert msgs.deliverable("nacht-s12") is True
        assert msgs.deliverable("NACHT-S12") is False
