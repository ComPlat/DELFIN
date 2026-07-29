"""An explicit delegation request stays visible until it is honoured.

Field case 20260729-155339: the user named five areas to delegate
("Architektur, Tetris, Snake, Benutzeroberfläche und Tests"), the agent
used no sub-agent and never said why — although the prompt rule
requiring exactly that was already in place. A static rule was too weak,
so the open request becomes a per-turn block that disappears the moment
a sub-agent runs.
"""

from __future__ import annotations

from delfin.agent.engine import AgentEngine


class _Stub(AgentEngine):
    """Bare instance: the block builders only read message/tool state."""

    def __init__(self):            # noqa: D107 - deliberately no super()
        self.messages: list[dict] = []
        self._session_tool_names: set[str] = set()
        self._delegation_satisfied = False


def _eng(*user_texts: str) -> _Stub:
    e = _Stub()
    e.messages = [{"role": "user", "content": t} for t in user_texts]
    return e


def test_request_detected_in_both_languages():
    for text in ("nutze Sub-Agenten für die Aufgabe",
                 "please use subagents for this",
                 "Nutze nach Möglichkeit Sub-Agenten für klar getrennte "
                 "Aufgabenbereiche",
                 "delegiere die Teilaufgaben",
                 "delegate the independent parts"):
        assert _eng(text)._delegation_was_requested(), text


def test_ordinary_requests_do_not_trigger():
    for text in ("baue ein Voila Dashboard mit zwei Spielen",
                 "split the work into steps and work through them",
                 "arbeite sorgfältig und teste alles"):
        assert not _eng(text)._delegation_was_requested(), text


def test_block_is_shown_while_the_request_is_open():
    e = _eng("nutze Subagents für Tetris und Snake")
    block = e._build_unmet_delegation_block()
    assert "Open request: delegation" in block
    assert "subagent" in block


def test_block_disappears_once_a_subagent_ran():
    e = _eng("nutze Subagents für Tetris und Snake")
    e._session_tool_names = {"write_file", "mcp__kit-coding__subagent"}
    assert e._build_unmet_delegation_block() == ""
    # And stays gone for the rest of the session.
    e._session_tool_names = {"write_file"}
    assert e._build_unmet_delegation_block() == ""


def test_no_block_without_a_request():
    e = _eng("baue mir eine App")
    assert e._build_unmet_delegation_block() == ""


def test_builder_never_raises_on_broken_state():
    e = _Stub()
    e.messages = [{"role": "user"}, {"bad": "shape"}, None]  # type: ignore
    assert e._build_unmet_delegation_block() == ""
