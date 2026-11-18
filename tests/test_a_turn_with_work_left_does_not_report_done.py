"""A turn that stops with work left must not report itself as finished.

The engine stops a turn and says why: the task ledger was still open
(``end_turn_open_tasks``), the ledger could not be read
(``end_turn_tasks_unknown``), the model announced a tool call and made
none (``promised_no_tool_call``), or a quota wait was longer than the
run's remaining budget (``quota_wait_exceeds_budget``). The comment at
the stop site states the intent: "the supervisor sees WHY it ended and
resumes it".

None of the four was in ``_ENGINE_STOP``, and ``classify`` returns
``"done"`` for anything unmapped. So all four reached the CLI's JSON, the
service supervisor and every measurement as ``stop_reason: "done"``: a
run that stopped with its task list open was indistinguishable from one
that finished. For a benchmark harness or an unattended supervisor that
difference is the whole question.

The second half of this file is the trap that comes with the fix: the
announcement continuation asked ``stop_reason != "done"`` and would have
switched itself off for exactly the turns that most need it. It asks what
it means now (``ended_by_itself``), so reporting the truth did not cost a
behaviour.
"""

from __future__ import annotations

import inspect

from delfin.agent import turn_outcome as TO


def _reason(engine_stop: str) -> str:
    return TO.classify(error="", engine_stop_reason=engine_stop,
                       denials=[], ended_on_denial=False)


def test_an_open_task_ledger_is_not_done():
    for engine_stop in ("end_turn_open_tasks", "end_turn_tasks_unknown",
                        "promised_no_tool_call"):
        assert _reason(engine_stop) == "incomplete", engine_stop


def test_a_quota_wait_longer_than_the_budget_is_a_ceiling():
    assert _reason("quota_wait_exceeds_budget") == "budget"


def test_a_plain_end_is_still_done_and_a_ceiling_still_its_own_word():
    assert _reason("end_turn") == "done"
    assert _reason("") == "done"
    assert _reason("max_tool_rounds") == "max_turns"
    assert _reason("run_budget_s") == "timeout"


def test_every_stop_reason_has_an_exit_code():
    """The contract of this module: each reason a program can read carries
    an exit code. A new reason without one is a reason nobody can act on."""
    assert set(TO.STOP_REASONS) == set(TO.EXIT_CODES)
    assert TO.EXIT_CODES["incomplete"] not in {
        v for k, v in TO.EXIT_CODES.items() if k != "incomplete"}
    assert TO.EXIT_CODES["done"] == 0


def test_a_turn_that_ended_by_itself_may_still_be_continued():
    assert TO.ended_by_itself("done")
    assert TO.ended_by_itself("incomplete")
    assert TO.ended_by_itself("")
    for ceiling in ("max_turns", "budget", "timeout", "error", "refused",
                    "needs_approval"):
        assert not TO.ended_by_itself(ceiling), ceiling


def test_the_continuation_asks_what_it_means_not_whether_it_is_done():
    """Pinned on the source: a bare ``!= "done"`` here is the regression
    this file exists to prevent, and it would be invisible in behaviour
    until a turn happened to end with its task list open."""
    from delfin.agent import cli

    for fn in (cli._continue_announced,):
        src = inspect.getsource(fn)
        assert "ended_by_itself" in src, (
            "the continuation must ask whether the turn ended by itself")
        assert '!= "done"' not in src, (
            "comparing against done alone switches the continuation off for "
            "a turn that ended with work left")
