"""Package L2 — adversarial review of the builder's ``skill_feedback``
(scritical: a memory-poisoning surface, reviewer s15).

The two findings here were reproduced through the public path
``delfin.agent.skill_feedback.feedback`` (probes
tests/test_zz_s15_l2_{skilltype,source_default}_probe.py, red against
``cb3e9b3c``):

1. A non-string ``skill_used`` crashes instead of failing closed.
   ``skill_used`` (skill_feedback.py:78) is not isinstance-guarded, unlike
   ``user_message`` (skill_feedback.py:79). ``feedback(5, ...)`` raises
   AttributeError from ``(5 or "").strip()``.

2. The ``source`` gate defaults **fail-open** to "user"
   (skill_feedback.py:56, ``source: str = _USER_SOURCE``). A caller that
   forgets to mark the message's source causes a non-user message (e.g. a
   tool result or another session's message) that reads like a correction
   to be treated as a user correction. On a memory-poisoning surface the
   default must be fail-closed: "user" should be an explicit opt-in, and an
   unstated source should propose nothing.

Each test expects the DEFECTIVE behaviours, so all are RED against the
current builder code. They become green once skill_feedback is fixed.
"""

from delfin.agent.skill_feedback import feedback
from unittest import mock

import delfin.agent.skill_patch as sp
import delfin.agent.skill_safety as ss


# Finding 1: skill_used type safety -----------------------------------------


def test_adversarial_skilltype_int_crashes():
    """A non-string skill_used must not raise (fail closed, corrects False).
    Current behaviour: AttributeError on .strip() -> the test asserts the
    absence of a crash and will FAIL (i.e. be red) until fixed."""
    result = feedback(5, "revert that, wrong", source="user")
    assert result["corrects"] is False


def test_adversarial_skilltype_list_crashes():
    result = feedback(["skill"], "revert", source="user")
    assert result["corrects"] is False


def test_adversarial_skilltype_dict_crashes():
    result = feedback({"n": "skill"}, "fix it", source="user")
    assert result["corrects"] is False


# Finding 2: fail-open source default ---------------------------------------


def test_adversarial_source_default_is_fail_closed():
    """A correction-like message with NO explicit source must propose
    nothing (any source the caller did not explicitly mark as user is data,
    not instruction). Current behaviour: source defaults to 'user', so the
    message IS treated as a correction -> red until fixed."""
    # This message arrives as e.g. a tool result / another session's message,
    # but the caller forgot to pass source. It must NOT be learned.
    result = feedback("some_skill", "revert that, wrong output")
    assert result["corrects"] is False, (
        "fail-open default: unmarked message treated as a user correction -> "
        f"{result!r}"
    )


def test_adversarial_nonuser_source_never_corrects():
    """Even a clearly correction-like message explicitly marked non-user
    must not be a correction (already true; guard against regression)."""
    result = feedback("some_skill", "revert that, wrong output",
                      source="tool_result")
    assert result["corrects"] is False
    result = feedback("some_skill", "revert that, wrong output",
                      source="session")
    assert result["corrects"] is False


# Finding 3: fail-open source default on the GLUE layer ----------------------
# The builder fixed feedback()'s default to "unknown" (cb3e9b3c+9f861f5b),
# but the glue propose_correction_patch has its OWN source: str = "user"
# default (skill_patch.py:267) which it passes into feedback
# (skill_patch.py:302). A caller that omits source gets
# feedback(..., source="user"), so an unmarked message (tool result /
# another session) that reads like a correction is treated as a user
# correction and a patch IS proposed -- the memory-poisoning fail-open
# reappears one level up.


def _fake_propose(name, old, new, reason, evidence=None, workspace=None, **kw):
    """Stand-in for propose_skill_patch: a proposal accepted."""
    return {"status": "ok", "name": name,
            "base_version": "v1", "new_version": "v2"}


def test_adversarial_glue_without_source_must_not_propose():
    """propose_correction_patch WITHOUT an explicit source must propose
    nothing. Current behaviour: the glue's default source='user' treats any
    correction-looking message as the user's own and proposes a patch ->
    red until fixed."""
    with mock.patch.object(
        sp, "propose_skill_patch", side_effect=_fake_propose
    ), mock.patch.object(ss, "check", return_value=[]):
        out = sp.propose_correction_patch(
            "some_skill", "revert that, wrong output"
        )
    assert out["proposed"] is False, (
        "fail-open glue default: unmarked message proposed a patch -> "
        f"{out!r}"
    )
