"""A go-ahead is acted on, and a choice is decided or asked as a choice.

Session 78c36f8c (2026-10-19, gpt-5-nano building a browser game): after
4 of 7 go-aheads ("Ja", "ja mach", "Patch a + B") the turn called no tool
that acts and answered with the plan again; the user ended up writing
"warum machst du nciht einfach??". Choices (width 10, 12 or 16) came as
numbered prose although ask_user_question was bound -- "sowas kannst du
doch als rückfrage machen". Both answers end in a question, which is
exactly where the announced-action trigger stands down.

Driven through the real tool loop with a fake streaming endpoint.
"""

from __future__ import annotations

import pytest

from delfin.agent import go_ahead as G
from delfin.agent import mcp_client as M


# -- the detectors ----------------------------------------------------------

_OFFER = "Ich kann die Breite per Overlay einstellbar machen. Soll ich?"
_DONE = "Fertig: Breite 10, Overlay eingebaut, Syntax-Check grün."


@pytest.mark.parametrize("text", [
    "Ja", "ja bitte den muss man ja nciht sehen", "ok", "los", "go",
])
def test_a_bare_yes_answers_an_offer(text):
    assert G.is_go_ahead(text, _OFFER)


@pytest.mark.parametrize("text", ["Ja", "ok", "passt", "go", "gut"])
def test_a_bare_yes_after_a_finished_report_is_an_acknowledgement(text):
    assert not G.is_go_ahead(text, _DONE)
    assert not G.is_go_ahead(text)


@pytest.mark.parametrize("text", [
    "ja mach", "Patch a + B", "du sollst alles angehen umd umsetzten",
    "warum machst du nciht einfach??", "go ahead", "Yes, do it",
    "du musst das fixen und selbst testen dann siehst du das doch oder?",
])
def test_an_order_counts_whatever_came_before(text):
    assert G.is_go_ahead(text, _DONE)


@pytest.mark.parametrize("text", [
    "was ist der stand?", "funktinert es jetzt?", "ja?", "nein, warte",
    "Wie funktioniert apply_patch?", "der kasten ist zu klein oder", "",
    "ja " + "sehr lange Beschreibung " * 10,
])
def test_a_question_or_a_no_is_not(text):
    assert not G.is_go_ahead(text)


def test_a_turn_that_acted_did_not_ignore_it():
    assert G.ignored_go_ahead("ja mach", ["read_file"])
    assert G.ignored_go_ahead("ja", ["read_file"], _OFFER)
    assert not G.ignored_go_ahead("ja mach", ["read_file", "edit_file"])
    assert not G.ignored_go_ahead("ja mach", ["mcp__x__bash"])


_PROSE_CHOICE = (
    "Gute Frage. Vorschläge (bitte eine Zahl wählen):\n\n"
    "- Breite 10 (outer=5) – Standard\n- Breite 12 (outer=6)\n"
    "- Breite 16 (outer=8)\n\nWelche Breite soll es sein?")


@pytest.mark.parametrize("text", [
    _PROSE_CHOICE,
    "Zwei Wege:\n\n1. Patch A: URL-Parameter\n2. Patch B: Overlay\n\n"
    "Sag mir, welchen du willst.",
    "Soll ich Option A oder Option B nehmen?",
])
def test_a_choice_written_as_prose_is_recognised(text):
    assert G.is_choice_question(text)


@pytest.mark.parametrize("text", [
    "Erledigt:\n\n- Breite ist 10\n- Mittelbalken entfernt\n\nTeste es unter "
    "http://127.0.0.1:8038/.",
    "Soll ich den Server neu starten?",
    "",
    "Was soll das Spiel können?",
])
def test_a_report_or_an_open_question_is_not(text):
    assert not G.is_choice_question(text)


# -- through the tool loop --------------------------------------------------

class _Delta:
    def __init__(self, content=None):
        self.content = content
        self.tool_calls = None


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


def _say(text):
    return _Stream([_Chunk([_Choice(_Delta(content=text), finish="stop")],
                           usage=_Usage())])


@pytest.fixture
def client(monkeypatch, tmp_path):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    from delfin.agent.api_client import KitToolPermissions, create_client
    cl = create_client(backend="api", provider="ollama",
                       model="qwen2.5-coder:7b", cwd=str(tmp_path))
    cl.set_permissions(KitToolPermissions(
        workspace=tmp_path, mode="default", task_session_id="sess-ga"))
    return cl


def _run(client, user, answers, before=None):
    sent = []

    def _create(**kwargs):
        sent.append(kwargs.get("messages") or [])
        return _say(answers[min(len(sent) - 1, len(answers) - 1)])

    client.client.chat.completions.create = _create
    history = list(before or []) + [{"role": "user", "content": user}]
    events = list(client.stream_message("sys", history, max_tokens=100))
    notices = "".join(getattr(e, "text", "") or "" for e in events
                      if getattr(e, "type", "") == "notice")
    return sent, notices


def _rows(sent):
    return " ".join(str(m.get("content") or "") for m in sent[-1])


_PLAN_AGAIN = ("Was ich jetzt mache:\n- BX/BZ aus dem URL-Parameter lesen\n"
               "- Overlay vor dem Start\n\nSoll ich das so umsetzen?")


def test_a_plan_after_a_go_ahead_is_sent_back_to_act(client):
    sent, notices = _run(client, "ja mach", [_PLAN_AGAIN, "Erledigt."])
    assert len(sent) == 2
    assert "go-ahead given" in notices
    assert "The user already told you to go ahead" in _rows(sent)


def test_a_bare_yes_to_an_offer_is_sent_back_to_act(client):
    sent, notices = _run(client, "ja", [_PLAN_AGAIN, "Erledigt."],
                         before=[{"role": "user", "content": "mach es breiter"},
                                 {"role": "assistant", "content": _OFFER}])
    assert len(sent) == 2 and "go-ahead given" in notices


def test_an_ok_after_a_finished_report_ends_the_turn(client):
    sent, notices = _run(client, "ok", ["Gern."],
                         before=[{"role": "user", "content": "mach es breiter"},
                                 {"role": "assistant", "content": _DONE}])
    assert len(sent) == 1 and "go-ahead" not in notices


def test_it_fires_once_per_turn(client):
    sent, _ = _run(client, "ja mach",
                   ["Soll ich?", "Soll ich wirklich?", "Soll ich jetzt?"])
    assert len(sent) == 2


def test_a_question_about_state_is_answered_not_pushed(client):
    sent, notices = _run(client, "was ist der stand?",
                         ["Stand: Breite 6, Overlay fehlt noch."])
    assert len(sent) == 1 and "go-ahead" not in notices


def test_a_prose_choice_is_asked_with_the_dialog_when_there_is_one(client):
    client._permissions.ask_user_callback = lambda q: {"answers": ["10"]}
    sent, notices = _run(client, "Wie breit soll es sein? Mach Vorschläge",
                         [_PROSE_CHOICE, "Ich nehme 10."])
    assert len(sent) == 2
    assert "choice written as text" in notices
    assert "ask_user_question" in _rows(sent)


def test_without_a_dialog_the_note_does_not_name_it(client):
    sent, _ = _run(client, "Wie breit soll es sein? Mach Vorschläge",
                   [_PROSE_CHOICE, "Ich nehme 10."])
    assert len(sent) == 2
    assert "sensible default" in _rows(sent)
    assert "with ask_user_question" not in _rows(sent)


def test_plan_mode_is_left_alone(client):
    client._permissions.mode = "plan"
    sent, _ = _run(client, "ja mach", [_PLAN_AGAIN])
    assert len(sent) == 1
