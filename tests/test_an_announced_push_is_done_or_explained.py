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


def test_a_participle_closing_the_turn_is_an_announcement():
    """One of the four turns ended "Pushing `jn/agent-chat-resizer` to the
    remote now." -- no first person, no future auxiliary, the same
    meaning, and the detector missed exactly that turn."""
    assert TC.announced_action("Pushing `jn/agent-chat-resizer` to the remote now.")
    assert TC.announced_action("Running the covering tests now.")


def test_the_participle_form_stays_narrow():
    """It must end the text and carry "now", or a report that merely
    contains a participle would be read as an announcement."""
    for text in (
            "Pushing it would need a grant the user has not given, so I stopped.",
            "Pushed the branch; CI is green now.",
            "The splitter is bound in tab_agent.py:14920.",
            "Task 3: pushing the branch",
            ""):
        assert not TC.announced_action(text), text


def test_the_plan_mode_redirect_also_speaks_to_the_user():
    """The redirect that exists today is addressed to the MODEL: it tells
    it to call exit_plan_mode. In report 20261005-140408 the model did
    not, and the person was left with "Let me push the branch now" and no
    sign of why nothing ran. The same place must also say it to them.

    Asserted on the source, because the branch sits inside the streaming
    loop: reaching it needs a provider round whose results carry the
    plan-mode refusal, which is the loop's own integration test, not a
    unit of this module.
    """
    import pathlib

    from delfin.agent import api_client

    src = pathlib.Path(api_client.__file__).read_text(encoding="utf-8")
    i = src.index("_plan_redirect_sent = True")
    window = src[i:i + 2600]
    assert 'yield StreamEvent(type="notice"' in window, (
        "the redirect steers the model and tells the person nothing")
    assert "Plan mode is on" in window
    assert "/mode solo" in window or "exit_plan_mode" in window, (
        "a reason that names no way out is what provoked the four repeats")
    assert "_plan_mode_refused_this_turn" in window


def test_the_sentence_reaches_the_answer_and_not_only_a_notice():
    """A notice scrolls past; the answer is what the user scrolls back
    to. The engine appends the sentence when the client reports that
    plan mode refused something this turn, and clears the flag so the
    next turn does not inherit it."""
    import pathlib

    from delfin.agent import engine as E

    src = pathlib.Path(E.__file__).read_text(encoding="utf-8")
    i = src.index("return full_response + _guard_note")
    window = src[max(0, i - 1400):i]
    assert "blocked_by_plan_mode" in window, (
        "the sentence is written and nobody calls it")
    assert "_plan_mode_refused_this_turn = False" in window, (
        "the flag must be cleared, or every later turn carries the note")
