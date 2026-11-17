"""A question to the person is shown once, and so is its hint.

Seen on a real session: under "Code Agent is asking you: ..." the card
read "Reply to the Code Agent." twice -- the gate put its hint into the
card's body and the card drew the hint again -- and the question itself
stood twice in the chat: at the end of the agent's message ("QUESTION:
...") and again in the card.

Now the hint is drawn once (the card's hint line), and a question the
agent's own message already ends with (compared as normalised text) is not
repeated by the card. Checked in the in-kernel tab's chat as it renders,
and in the live view of a service session (its chat plus its question
card). Control: a question that is NOT at the end of the message is still
shown in the card.
"""

from __future__ import annotations

import html
import re
import time

import pytest

pytest.importorskip("ipywidgets")

from delfin.agent import scheduler as S                          # noqa: E402
from delfin.dashboard import live_sessions as LS                 # noqa: E402
from delfin.dashboard import tab_agent                           # noqa: E402
from delfin.dashboard.context import DashboardContext            # noqa: E402

QUESTION = "Which solvent model should the scan use?"


def _visible_text(value: str) -> str:
    """The chat's text as a reader sees it: the inert template's markup
    without tags, entities decoded, whitespace collapsed."""
    text = re.sub(r"<[^>]+>", " ", value)
    return " ".join(html.unescape(text).split())


@pytest.fixture(autouse=True)
def _scheduler(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))


def _tab(tmp_path):
    for name in ("calc", "archive", "office", "work"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc", archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "work")
    return tab_agent.create_tab(ctx)


def _ask(refs, reply: str) -> str:
    state = refs["state"]
    state["chat_messages"].append({"role": "assistant", "content": reply,
                                   "role_label": "Code Agent"})
    assert refs["question_gate"]("solo_agent", reply) is True
    refs["refresh_chat"]()
    return _visible_text(refs["chat_html"].value)


def test_the_hint_and_a_closing_question_appear_once(tmp_path):
    _w, refs = _tab(tmp_path)
    try:
        shown = _ask(refs, f"I looked at the inputs.\n\nQUESTION: {QUESTION}")
        assert shown.count("Code Agent is asking you") == 1
        assert shown.count("Reply to the Code Agent.") == 1, shown
        assert shown.count(QUESTION) == 1, shown
    finally:
        refs["shutdown"]()


def test_a_question_inside_the_message_is_still_in_the_card(tmp_path):
    _w, refs = _tab(tmp_path)
    try:
        # Control: the question is followed by more text, so the card is
        # where it stands on its own.
        shown = _ask(refs, f"QUESTION: {QUESTION}\n\nMeanwhile I keep the old settings.")
        assert shown.count(QUESTION) == 2, shown          # the message, and the card
        assert shown.count("Reply to the Code Agent.") == 1, shown
    finally:
        refs["shutdown"]()


def test_the_live_view_does_not_repeat_a_question_its_message_ends_with(tmp_path):
    from delfin.agent.service.client import ServiceClient
    from delfin.dashboard import live_tab
    from tui_service import ServiceProcess

    (tmp_path / "svc").mkdir()
    svc = ServiceProcess(tmp_path / "svc")
    try:
        tui = ServiceClient(svc.socket, client="delfin-tui")
        svc.clients.append(tui)
        dash = LS.LiveClient(svc.socket)
        sid = svc.session(tui)
        _w, refs = live_tab.build_live_tab(dash, sid)
        dash.prompt(sid, "ask aloud")
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline and not refs["approval_box"].children:
            time.sleep(0.05)
        refs["render_now"]()
        card = refs["approval_box"].children[0]
        card_text = _visible_text(card.children[0].value)
        chat = _visible_text(refs["chat"].value)
        question = refs["transcript"].pending and next(iter(refs["transcript"].pending.values()))["question"]
        assert question and question in chat                      # said in the message
        assert (chat + " " + card_text).count(question) == 1, (chat, card_text)
        refs["close"]()
        dash.close()
    finally:
        svc.stop()
