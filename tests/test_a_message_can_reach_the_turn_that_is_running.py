"""Typing during a turn: queue it, or put it inside the turn.

Queue-not-inject was the original decision and it stays the default. A
queued message is never lost and always lands in a context the user could
see; an injected one arrives in the middle of a tool loop, in a context
they cannot, and the model may act on it under whatever permissions are
in force at that moment. Those are different promises.

So they are different keys rather than a session flag. Enter queues,
Ctrl+G steers, both always available — a mode you can forget you are in
is the wrong place to keep a distinction that matters per message.

The engine end of this already existed and is not touched here: the
client's steer inbox is drained between tool rounds AND at a turn that
would otherwise end, so a message sent late continues the turn instead of
being answered next time. What was missing was a way to reach it from a
keyboard.
"""

from __future__ import annotations

import io

from delfin.agent import repl_keys as rk
from delfin.agent import repl_commands as rc
from delfin.agent.repl import ReplOptions, TerminalAgent


class _Tty(io.StringIO):
    def isatty(self):
        return True


class _Engine:
    """Enough engine for the controller, with a recording steer inbox."""

    def __init__(self, *, accepts: bool = True):
        self.accepts = accepts
        self.steered: list[str] = []
        self.session_id = "sid"
        self.messages: list[dict] = []
        self.kit_permissions = None
        self.mode = "solo"
        self.client = type("C", (), {"model": "m"})()
        self.provider = "kit"

    def steer(self, text):
        if not self.accepts:
            return False
        self.steered.append(text)
        return True

    def get_status(self):
        return {}


def _agent(**kw):
    err = _Tty()
    agent = TerminalAgent(_Engine(**kw), out=_Tty(), err=err,
                          opts=ReplOptions(color="never"))
    return agent, agent.engine, err


# ---------------------------------------------------------------------------
# The decoder: two keys, one buffer
# ---------------------------------------------------------------------------

def test_ctrl_g_produces_a_steer_carrying_the_line():
    d = rk.KeyDecoder()
    d.feed("stop and read the log")
    events = d.feed("\x07")
    steers = [e for e in events if e.kind == rk.STEER]
    assert len(steers) == 1
    assert steers[0].text == "stop and read the log"
    assert d.buffer == "", "the line is consumed, exactly as Enter consumes it"


def test_enter_is_still_a_queue_and_not_a_steer():
    """The control. Without it, this could pass by turning Enter into steer."""
    d = rk.KeyDecoder()
    d.feed("later please")
    events = d.feed("\r")
    assert [e.kind for e in events] == [rk.SUBMIT]
    assert rk.STEER not in [e.kind for e in events]


def test_ctrl_g_inside_a_paste_is_a_character():
    d = rk.KeyDecoder()
    events = d.feed("\x1b[200~a\x07b\x1b[201~")
    assert rk.STEER not in [e.kind for e in events]
    assert "\x07" in d.buffer


def test_an_empty_line_steers_nothing():
    d = rk.KeyDecoder()
    events = d.feed("\x07")
    steers = [e for e in events if e.kind == rk.STEER]
    assert len(steers) == 1 and steers[0].text == ""


# ---------------------------------------------------------------------------
# The controller: where it goes
# ---------------------------------------------------------------------------

def test_a_steered_line_reaches_the_engine_and_not_the_queue():
    agent, engine, err = _agent()
    agent._steer("check the failing test first")
    assert engine.steered == ["check the failing test first"]
    assert agent.queued == []
    assert "running turn" in err.getvalue()


def test_a_queued_line_does_not_reach_the_engine():
    agent, engine, _err = _agent()
    agent.queued.append("afterwards, tidy up")
    assert engine.steered == []


def test_a_backend_without_the_inbox_says_so_and_keeps_the_message():
    """The message is the user's. Dropping it loses work; pretending it
    landed is worse, because the next thing they do assumes it did."""
    agent, engine, err = _agent(accepts=False)
    agent._steer("do it now")
    assert engine.steered == []
    assert agent.queued == ["do it now"], "not lost"
    out = err.getvalue()
    assert "no message mid-turn" in out
    assert "queued" in out


def test_a_raising_backend_is_the_same_as_a_refusing_one():
    agent, engine, err = _agent()

    def _boom(_text):
        raise RuntimeError("client is gone")

    engine.steer = _boom
    agent._steer("do it now")
    assert agent.queued == ["do it now"]


def test_an_empty_steer_is_not_sent():
    agent, engine, err = _agent()
    agent._steer("   ")
    assert engine.steered == []
    assert agent.queued == []
    assert err.getvalue() == ""


# ---------------------------------------------------------------------------
# Discoverability — a key nobody can find is a key nobody uses
# ---------------------------------------------------------------------------

def test_the_queue_notice_names_the_other_key():
    """The moment someone learns Ctrl+G exists is the moment they have
    just queued something and might have wanted it sooner."""
    import inspect
    from delfin.agent import repl
    src = inspect.getsource(repl.TerminalAgent._on_key)
    assert "ctrl+g" in src


def test_the_keys_command_exists_and_names_them():
    out = rc.BUILTINS["/keys"].handler(None, "").output
    assert "ctrl+g" in out
    assert "esc" in out
    assert "shift+tab" in out


def test_every_key_the_decoder_produces_is_documented():
    """The anti-drift check.

    A new key added to the decoder and not to this table is a feature
    only its author knows about.
    """
    documented = {event for event, _key, _desc in rc._key_rows()}
    produced = {
        getattr(rk, name) for name in rk.__all__
        if name.isupper() and isinstance(getattr(rk, name), str)
    }
    missing = produced - documented
    assert not missing, f"keys nobody is told about: {sorted(missing)}"


def test_help_lists_the_keys_command():
    from delfin.agent import help_gen
    text = help_gen.generate_help(rc.palette_rows())
    assert "/keys" in text


# ---------------------------------------------------------------------------
# What the engine end already guarantees, pinned so it stays true
# ---------------------------------------------------------------------------

def test_the_inbox_is_drained_after_the_last_tool_call_too():
    """A message sent while the model was writing its final answer must
    continue the turn, not be answered in the next one."""
    import inspect
    from delfin.agent import api_client
    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    assert src.count("_drain_steer()") >= 2, (
        "one drain between rounds is not enough: a turn that ends without "
        "another tool call would swallow the message")


def test_an_undelivered_steer_is_named_rather_than_dropped():
    import inspect
    from delfin.agent import api_client
    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    assert "Not delivered" in src
