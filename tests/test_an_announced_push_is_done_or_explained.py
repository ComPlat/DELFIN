"""An announced action is carried out, or the user is told why not.

Field reports 20261005-135825 and 20261005-140408, one session six
minutes apart: "kann nicht pushen. agent auch bei erlaubnis warum?" and
"warum macht er es nicht und hoehrt einfach auf?".

What the transcript shows, at engine-message level:

    [108] tool      blocked: `git push` ... the user has not asked for a
                    push since their last message
    [112] user      push den branch
    [113] assistant The user explicitly asked to push the branch. Let me
                    push `jn/agent-chat-resizer` to the remote.
    [114] user      push
    [115] assistant Explicit push confirmed twice. Let me push the branch now.
    [116] user      push den brnach
    [117] assistant User confirmed the push explicitly. Pushing now.
    [118] user      push
    [119] assistant User has explicitly asked for the push multiple times.
                    Let me push the branch now.

Four turns ended on an announcement with no tool call in them. The
consent matcher was not at fault -- `_asks_for_push` returns True for
all four spellings, measured. Two mechanisms were:

1. The continuation detector requires OPEN tasks. The task list at
   [110] was `{"completed": 5}`, so a turn that ended announcing a push
   with every task done was never nudged. The condition asked "is there
   work left on the list" when the question is "did this turn announce
   an action it did not take".

2. The session was in PLAN MODE ([122] refused `bash` read-only, [124]
   the system note). The plan-mode refusal is addressed to the model --
   "call exit_plan_mode with the plan" -- and says nothing to the
   person, who saw "Let me push the branch now" four times and no push.
   The model worked it out at [125] and still did not say it.
"""

from __future__ import annotations

from delfin.agent import turn_continuation as TC


_ANNOUNCED = "Let me push the branch now."
_REAL_ANSWERS = (
    "The user explicitly asked to push the branch. Let me push "
    "`jn/agent-chat-resizer` to the remote.",
    "Explicit push confirmed twice. Let me push the branch now.",
    "User confirmed the push explicitly. Pushing `jn/agent-chat-resizer` "
    "to the remote now.",
    "User has explicitly asked for the push multiple times. Let me push "
    "the branch now.",
)


def test_an_announcement_is_followed_up_with_no_open_tasks():
    """The defect, as the report produced it: every task completed and
    the turn still ended on "Let me push the branch now"."""
    cont, note = TC.should_continue(_ANNOUNCED, 0)
    assert cont, (
        "a turn that announces an action it did not take is the thing to "
        "follow up on; whether a task list still has rows is a different "
        "question")
    assert note


def test_every_answer_from_the_report_is_followed_up():
    for answer in _REAL_ANSWERS:
        cont, _ = TC.should_continue(answer, 0)
        assert cont, answer


def test_open_tasks_still_follow_up():
    """Unchanged: the case the detector was built for."""
    assert TC.should_continue(_ANNOUNCED, 3)[0]


def test_a_hand_back_is_never_followed_up():
    """The guard that keeps this from nagging. A turn that gives the work
    back to the user is finished, with or without open tasks."""
    for answer in ("Let me know what you would like next.",
                   "I'll keep you posted when CI finishes.",
                   "What should I do about the failing test?",
                   "I am done.",
                   ""):
        for tasks in (0, 4):
            assert not TC.should_continue(answer, tasks)[0], (answer, tasks)


def test_a_plan_mode_refusal_is_said_to_the_user():
    """Report 2. In plan mode nothing runs, so nudging the model to
    continue cannot help -- the person has to be told, in the answer they
    read, that plan mode is the reason and how to leave it."""
    said = TC.blocked_by_plan_mode(_ANNOUNCED, refusals=2)
    assert said, "the user is told nothing about why the action did not run"
    low = said.lower()
    assert "plan" in low
    assert "exit_plan_mode" in low or "/mode" in low, (
        "a reason that names no way out is what provoked four repeats")


def test_nothing_is_said_when_plan_mode_refused_nothing():
    assert not TC.blocked_by_plan_mode(_ANNOUNCED, refusals=0)


def test_nothing_is_said_when_the_turn_announced_nothing():
    """A read-only investigation in plan mode is the normal case and must
    not carry a warning."""
    assert not TC.blocked_by_plan_mode(
        "The splitter is bound in tab_agent.py:14920.", refusals=3)
