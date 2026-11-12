"""A question names exactly the answers it takes, on every surface.

Field report: the agent asked a question with 3 options; the hint said
"1-9", a generic range. Every surface that renders an ask_user_question --
the phone message, the TUI, the REPL dialog, the dashboard card -- now
takes the hint from one function (approval_answers.answer_hint):

* "Reply 1, 2 or 3", "1-N" for N > 3, free text allowed and said once;
* multi-select says how to pick several ("1,3");
* a number outside the range is refused with the real range
  (approval_answers.range_problem), never accepted or taken as free text.

Tested for N = 1, 2, 3 and 9, and multi-select.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import approval_answers as AA
from delfin.messenger import render

ROOT = Path(__file__).resolve().parents[1]


def _payload(n, multi=False):
    return {"question": "Which?", "options": [{"label": f"opt{i}"} for i in range(1, n + 1)],
            "multiSelect": multi}


@pytest.mark.parametrize("n,want", [
    (1, "Reply 1, or type your own answer."),
    (2, "Reply 1 or 2, or type your own answer."),
    (3, "Reply 1, 2 or 3, or type your own answer."),
    (9, "Reply 1-9, or type your own answer."),
])
def test_the_hint_names_the_real_range(n, want):
    assert AA.answer_hint(n, False) == want


def test_multi_select_says_how_to_pick_several():
    assert AA.answer_hint(3, True) == ("Reply one or more of 1, 2 or 3, separated by commas "
                                       "(e.g. 1,3), or type your own answer.")


@pytest.mark.parametrize("n,text,problem", [
    (3, "4", "Not an option: reply 1, 2 or 3."),
    (3, "0", "Not an option: reply 1, 2 or 3."),
    (9, "12", "Not an option: reply 1-9."),
    (1, "2", "Not an option: reply 1."),
    (3, "2", None), (3, "opt2", None), (3, "lieber Toluol", None),
])
def test_a_number_outside_the_range_is_refused(n, text, problem):
    assert AA.range_problem(text, _payload(n)) == problem


def test_multi_select_numbers_are_checked_each():
    assert AA.range_problem("1,4", _payload(3, True)) == "Not an option: reply 1, 2 or 3."
    assert AA.range_problem("1,3", _payload(3, True)) is None


@pytest.mark.parametrize("n", [1, 2, 3, 9])
def test_the_phone_question_card_numbers_the_options_and_names_the_range(n, monkeypatch):
    monkeypatch.setenv("DELFIN_MESSENGER_LANG", "en")
    text = render.question_card({"text": "Which?", "choices": [f"opt{i}" for i in range(1, n + 1)]})
    for i in range(1, n + 1):
        assert f"{i}. opt{i}" in text
    assert AA.answer_hint(n, False).rstrip(".") in text
    assert "1-9" not in text or n == 9


def test_the_tui_shows_the_real_range():
    src = (ROOT / "tui" / "src" / "footer" / "approval.tsx").read_text()
    assert '"[1-9] answer  "' not in src and '"[1-9] [space] pick  "' not in src
    assert "answerRange(" in src


def test_the_worker_refuses_an_out_of_range_number(monkeypatch):
    from delfin.agent import terminal_confirm as tc
    from delfin.agent.service import approvals
    monkeypatch.setattr(tc, "_publish_pending", lambda *a, **k: None)
    broker = approvals.ProtocolBroker(lambda *a, **k: None)
    req = broker._enqueue(tc.ConfirmRequest(kind=tc.ASK, tool="ask_user_question",
                                            payload=_payload(3)))
    rid = broker.request_id(req)
    with pytest.raises(approvals.InvalidAnswer) as refused:
        broker.answer(rid, text="7")
    assert "reply 1, 2 or 3" in str(refused.value)
    assert not req.resolved
    assert broker.answer(rid, text="2")[0] and req.decision["answers"] == ["opt2"]
