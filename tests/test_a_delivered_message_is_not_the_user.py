"""A message another session delivered is not the user's message.

Driven 2026-09-16 with two sessions side by side: the delivered text went
through the input box like a typed message, so (1) it became the receiving
session's title, and the title then wrapped itself into every later message
("[Message from the session "[Message from the session ..."); (2) its
English wrapper pinned the session's language to English; (3) "answer if it
needs one" read as an invitation, and the two sessions greeted each other in
a loop.
"""
import inspect

from delfin.agent import session_messages as SM
from delfin.dashboard import agent_sessions as AS


def test_a_delivered_message_does_not_name_the_session():
    delivered = SM.render({"from": "abc", "from_title": "other", "text": "Hallo"})
    state = {"chat_messages": [{"role": "user", "content": delivered},
                               {"role": "user", "content": "Lies bitte notes.txt"}]}
    assert AS._session_title(state) == "Lies bitte notes.txt"
    assert AS._session_title({"chat_messages": [{"role": "user", "content": delivered}]}) == "New session"


def test_a_delivered_message_does_not_pin_the_language():
    from delfin.agent.engine import AgentEngine
    src = inspect.getsource(AgentEngine._note_session_language)
    assert '"[Message from the session"' in src


def test_the_wrapper_asks_for_a_reply_only_to_a_request():
    text = SM.render({"from": "abc", "from_title": "other", "text": "Danke!"})
    assert "only if it asks you for something" in text
    assert "acknowledgement gets no reply" in text
    assert "if it needs one" not in text
