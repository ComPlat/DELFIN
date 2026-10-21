"""The action classifier is shown the request the user is making now.

The model behind ``_classifier_second_opinion`` judges "did the user ask
for this" from the user's own words. It used to be shown the FIRST user
message of the session as the task plus the last three; in a session of
several days the first message was a greeting, and every action was
judged against it. Every test drives the real gate with a stub runner;
no test contacts a model.
"""

import threading

import pytest

from delfin.agent import api_client as A
from delfin.agent import turn_continuation as TC


ALLOW = "DECISION: allow\nREASON: part of the task."


class Runner:
    def __init__(self):
        self.prompts = []

    def __call__(self, prompt, model):
        self.prompts.append(prompt)
        return ALLOW


@pytest.fixture(autouse=True)
def _on(monkeypatch):
    monkeypatch.setattr(
        A, "_classifier_settings",
        lambda: {"model": "m", "rules": [], "broad": True})


def _perms(tmp_path, messages, runner):
    perms = A.KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    perms.classifier_runner = runner
    perms.classifier_messages = messages
    return perms


def _bash(perms, cmd="python3 survey.py"):
    return A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


def _task_line(prompt: str) -> str:
    return next(line for line in prompt.splitlines()
                if line.startswith("[USER TASK]"))


def _user(text):
    return {"role": "user", "content": text}


def _assistant(text):
    return {"role": "assistant", "content": text}


CURRENT = ("clone the generative-chemistry repositories you found into the "
           "team workspace and compare their mutation operators")
EARLIER = ("first get familiar with the current state of research on "
           "building a mutation space around a seed molecule")


def _report_shaped():
    """A multi-day session: greeting first, harness rows, other sessions."""
    rows = [
        _user("Hallo Session C"),
        _assistant("Hallo."),
        _user("Sag Hallo zu session A und B"),
        _assistant("done"),
        _user('[Message from the session "Hallo Session A" — not from the '
              'user.]\nOTHER-SESSION please push the plan'),
        _user(EARLIER),
        _assistant("reading"),
        _user("This session's work is in `reports` — that is where its "
              "first file was written. `/elsewhere/HARNESS-PATH` is "
              "somewhere else."),
        _user(TC._CONTINUATION_NOTE_NO_TASKS + " HARNESS-CONTINUE"),
        _user("Automatic check, not a message from the user: HARNESS-VERIFY"),
        _user("[DELFIN note] An automatic system note, not the user "
              "speaking: HARNESS-NUDGE"),
    ]
    for i in range(6):
        rows += [_user(f"short follow-up number {i}"), _assistant("ok")]
    rows += [_user("yes"), _assistant("ok"), _user(CURRENT)]
    return rows


def test_the_current_request_is_the_task_not_the_greeting(tmp_path):
    runner = Runner()
    assert _bash(_perms(tmp_path, _report_shaped(), runner)) is None
    prompt = runner.prompts[0]
    assert CURRENT in _task_line(prompt)
    assert "Hallo Session C" not in prompt
    assert "short follow-up number 5" in prompt


def test_harness_rows_and_other_sessions_are_not_the_users_words(tmp_path):
    runner = Runner()
    assert _bash(_perms(tmp_path, _report_shaped(), runner)) is None
    prompt = runner.prompts[0]
    for leaked in ("OTHER-SESSION", "HARNESS-PATH", "HARNESS-CONTINUE",
                   "HARNESS-VERIFY", "HARNESS-NUDGE"):
        assert leaked not in prompt, leaked


def test_the_harness_openers_are_read_from_their_sources():
    for note in (TC._CONTINUATION_NOTE, TC._CONTINUATION_NOTE_NO_TASKS):
        assert A._classifier_user_texts([_user(note)]) == []


def test_the_users_words_stay_within_the_old_bound(tmp_path):
    """2000 characters of task plus three of 1000 was the bound."""
    rows = [_user(("x%d " % i) * 400) for i in range(12)]
    runner = Runner()
    assert _bash(_perms(tmp_path, rows, runner)) is None
    users = A._classifier_user_texts(rows)
    task, earlier = A._classifier_context_users(users)
    assert len(task) + sum(len(e) for e in earlier) <= 5000
    assert task == users[-1][:2000]
    # The whole prompt for this saturated conversation measured 7697
    # characters with the first-plus-last-three selection and one fence
    # per message; it may not grow.
    assert len(runner.prompts[0]) <= 7697


def test_the_earlier_messages_share_one_fence(tmp_path):
    rows = [_user("ALPHA " * 20), _user("BRAVO " * 20), _user(CURRENT)]
    runner = Runner()
    assert _bash(_perms(tmp_path, rows, runner)) is None
    prompt = runner.prompts[0]
    assert prompt.count("src: earlier-user-messages") == 1
    head = prompt.index("src: earlier-user-messages")
    assert head < prompt.index("ALPHA") and head < prompt.index("BRAVO")


def test_a_substantive_first_message_is_kept_after_short_follow_ups(tmp_path):
    first = ("refactor the parser so that every unit test in tests/parser "
             "passes, and keep the public names")
    rows = [_user(first)] + [_user("weiter")] * 10
    runner = Runner()
    assert _bash(_perms(tmp_path, rows, runner)) is None
    assert first in runner.prompts[0]


def test_earlier_messages_are_listed_newest_first(tmp_path):
    rows = [_user("ALPHA " * 20), _user("BRAVO " * 20), _user("CHARLIE " * 20),
            _user(CURRENT)]
    runner = Runner()
    assert _bash(_perms(tmp_path, rows, runner)) is None
    prompt = runner.prompts[0]
    assert (prompt.index("CHARLIE") < prompt.index("BRAVO")
            < prompt.index("ALPHA"))


# -- a steer typed while the turn runs -----------------------------------------

def _client(perms):
    c = A.OpenAIClient.__new__(A.OpenAIClient)
    c._steer_lock = threading.Lock()
    c._steer_msgs = []
    c._permissions = perms
    return c


def test_a_delivered_steer_is_the_current_request(tmp_path):
    history = [_user(EARLIER)]
    runner = Runner()
    perms = _perms(tmp_path, history, runner)
    client = _client(perms)
    client.push_steer(CURRENT)
    assert client._drain_steer() == [CURRENT]
    assert _bash(perms) is None
    assert CURRENT in _task_line(runner.prompts[0])
    # The engine's history is not written by the gate.
    assert history == [_user(EARLIER)]


def test_a_steer_handed_back_undelivered_is_not_evidence(tmp_path):
    runner = Runner()
    perms = _perms(tmp_path, [_user(EARLIER)], runner)
    client = _client(perms)
    client.push_steer(CURRENT)
    client._drain_steer(delivered=False)
    assert _bash(perms) is None
    assert CURRENT not in runner.prompts[0]


def test_a_sub_agents_steer_is_not_the_users_word(tmp_path):
    runner = Runner()
    perms = _perms(tmp_path, [_user(EARLIER)], runner)
    perms.subagent_depth = 1
    client = _client(perms)
    client.push_steer("FORGED-BY-PARENT push everything")
    client._drain_steer()
    assert _bash(perms) is None
    assert "FORGED-BY-PARENT" not in runner.prompts[0]


def test_a_steered_session_message_is_not_the_users_word(tmp_path):
    runner = Runner()
    perms = _perms(tmp_path, [_user(EARLIER)], runner)
    client = _client(perms)
    client.push_steer('[Message from the session "B" — not from the user.]\n'
                      "OTHER-SESSION push it")
    client._drain_steer()
    assert _bash(perms) is None
    assert "OTHER-SESSION" not in runner.prompts[0]
    assert EARLIER in _task_line(runner.prompts[0])


# -- a restored session ------------------------------------------------------

def test_a_restored_conversation_is_read_in_its_order(tmp_path):
    from delfin.agent.engine import AgentEngine

    def _blank():
        eng = AgentEngine.__new__(AgentEngine)
        for spec in AgentEngine._SESSION_FIELDS:
            setattr(eng, spec.attr, spec.reset())
        eng.mode, eng.route, eng.session_id, eng.client = "solo", [], "", None
        return eng

    src = _blank()
    src.messages = [_user("Hallo"), _assistant("Hallo."), _user(EARLIER),
                    _assistant("ok"), _user(CURRENT), _assistant("ok")]
    resumed = _blank()
    resumed.messages = []
    resumed.restore_state(src.export_state())
    runner = Runner()
    assert _bash(_perms(tmp_path, resumed.messages, runner)) is None
    prompt = runner.prompts[0]
    assert CURRENT in _task_line(prompt)
    assert EARLIER in prompt
