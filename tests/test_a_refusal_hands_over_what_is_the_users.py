"""A refusal about consent offers the command; one about content never does.

Input: a refused command and the kind of refusal. Output: a line the user
can run with ``!``, or "".

Measured 2026-09-28: two agent sessions spent eleven and thirteen
attempts each on a push that could not work from their host, trying other
transports, other proxies, other temp directories. Nothing in a refusal
said "this one is yours", so they kept looking for a way it could be
theirs -- which is the behaviour a refusal must never provoke. The user
pushed by hand in the end.

The line this file defends:

  consent   the user may do it and the agent may not without being asked.
            A push names a branch on the project's own remote, and that
            boundedness is why it qualifies.

  content   the command itself is the problem. Offering it to a person
            who pastes it is the same harm by another hand, and a rule
            that can be satisfied by asking somebody else is not a rule.

Outbound transfers and writes to other machines sit in CONTENT although
both refusals mention a human, because neither is bounded by anything a
reader could check before running it.
"""

from __future__ import annotations

import pytest

from delfin.agent import handover


@pytest.mark.parametrize("kind", sorted(handover._CONSENT_KINDS))
def test_a_consent_refusal_offers_the_command(kind):
    line = handover.for_user("git push origin HEAD:work/x", kind=kind)
    assert "! git push origin HEAD:work/x" in line
    assert "run it yourself" in line


@pytest.mark.parametrize("kind", sorted(handover._CONTENT_KINDS))
def test_a_content_refusal_offers_nothing(kind):
    assert handover.for_user("rm -rf /tmp/anything", kind=kind) == ""


@pytest.mark.parametrize("kind", [
    "outbound_unattended", "remote_write",
])
def test_the_two_that_read_as_consent_are_still_refused(kind):
    """Both refusals say a human would make them fine. Neither is bounded
    by anything the person pasting it could check first."""
    assert kind in handover._CONTENT_KINDS
    assert handover.for_user("nc -X connect host 22", kind=kind) == ""


def test_an_unknown_kind_is_treated_as_content():
    """The safe default for a rule nobody has classified is the strict
    one -- a refusal reaches this before anybody decides what it is."""
    assert handover.for_user("anything at all", kind="a_new_rule") == ""
    assert handover.for_user("anything at all", kind="") == ""


def test_the_two_sets_do_not_overlap():
    assert not (handover._CONSENT_KINDS & handover._CONTENT_KINDS)


def test_only_bounded_kinds_are_handed_over():
    """A guard on the list itself: adding a kind here is a security
    decision, and this names the two that were argued for."""
    assert handover._CONSENT_KINDS == {"push_unrequested", "git_role"}


def test_an_empty_command_offers_nothing():
    assert handover.for_user("", kind="push_unrequested") == ""
    assert handover.for_user("   ", kind="push_unrequested") == ""


def test_the_command_is_offered_unchanged():
    """A handover that edits the command is a different command."""
    cmd = "git push origin HEAD:refs/heads/work/j1-grounding"
    assert cmd in handover.for_user(cmd, kind="push_unrequested")


def test_nothing_here_runs_anything():
    import inspect

    src = inspect.getsource(handover)
    for forbidden in ("subprocess", "os.system", "Popen", "socket",
                      "urllib", "requests", "eval("):
        assert forbidden not in src, forbidden


# -- the gate still decides, and still refuses ----------------------------

def test_the_push_gate_still_returns_a_refusal():
    """The handover is appended to a refusal. If it ever replaced one,
    the gate would have stopped refusing."""
    import pathlib

    src = (pathlib.Path(handover.__file__).resolve().parent
           / "api_client.py").read_text(encoding="utf-8")
    at = src.index("blocked: `git push` publishes to a shared remote")
    window = src[at:at + 900]
    assert 'kind="push_unrequested"' in window, (
        "the push refusal does not hand the command over")
    assert "+ _handover(" in window, (
        "the handover replaces the refusal instead of following it")


def test_no_content_refusal_calls_the_handover():
    """The grep that would catch a future wiring mistake."""
    import pathlib
    import re

    src = (pathlib.Path(handover.__file__).resolve().parent
           / "api_client.py").read_text(encoding="utf-8")
    used = set(re.findall(r'_handover\([^)]*kind="([a-z_]+)"', src))
    assert used <= handover._CONSENT_KINDS, (
        f"a content refusal hands its command over: "
        f"{sorted(used - handover._CONSENT_KINDS)}")
