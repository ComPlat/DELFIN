"""The worker's session (ProtocolAgent) and approval broker (ProtocolBroker)
subclass the REPL's TerminalAgent and TerminalConfirmBroker and replace
some of their private methods. A rename in repl.py or terminal_confirm.py
would turn such a replacement into a method nobody calls, with no failure
anywhere: the base's own version would then run in the worker (reading a
terminal it does not have, drawing a screen nobody sees).

Pinned in both directions: every name listed as an override exists on the
base class, and every method the subclass defines is listed either as an
override or as its own -- a new override has to be added here.
"""

from __future__ import annotations

import inspect

import pytest

from delfin.agent.repl import TerminalAgent
from delfin.agent.service import approvals, worker
from delfin.agent.terminal_confirm import TerminalConfirmBroker

AGENT_OVERRIDES = {
    "_announce_presence", "_begin_turn_resets_the_abort", "_checkpoint_session",
    "_continue_after_announcement", "_continue_after_endpoint_failure",
    "_continue_after_length", "_endpoint_read_key", "_handle_line",
    "_install_completer", "_install_sigwinch", "_line_ready", "_load_history",
    "_recall_later", "_restore_sigwinch", "_save_history", "_show_user_input",
    "_teardown_screen", "_wake_text", "turn",
}
AGENT_OWN = {
    "_back_at_the_prompt", "_follow_ups_off", "_leave", "_own_stream", "_read_next",
    "_release_stream", "_status", "_take_back_steers", "_took_line", "leave_now",
    "post", "refresh_state", "state", "steer_text",
}
BROKER_OVERRIDES = {"_audit_dialog", "_enqueue", "_wait", "abort_all", "reset_abort",
                    "resolve", "take"}
BROKER_OWN = {"_announce_resolved", "_apply", "_suggestion", "_waiting_changed", "answer",
              "close", "describe", "pending", "request_id"}

PAIRS = [
    pytest.param(worker._agent_class, TerminalAgent, AGENT_OVERRIDES, AGENT_OWN, id="agent"),
    pytest.param(lambda: approvals.ProtocolBroker, TerminalConfirmBroker, BROKER_OVERRIDES,
                 BROKER_OWN, id="broker"),
]


def _methods(cls) -> set[str]:
    return {n for n, v in vars(cls).items() if not n.startswith("__")
            and (inspect.isfunction(v) or isinstance(v, (property, staticmethod, classmethod)))}


@pytest.mark.parametrize("sub, base, overrides, own", PAIRS)
def test_every_override_names_a_method_the_base_has(sub, base, overrides, own):
    cls = sub()
    assert base in cls.__mro__
    missing = sorted(n for n in overrides if not hasattr(base, n))
    assert not missing, f"{base.__name__} no longer has {missing}: the override is dead"
    assert overrides <= _methods(cls)


@pytest.mark.parametrize("sub, base, overrides, own", PAIRS)
def test_every_method_of_the_subclass_is_pinned(sub, base, overrides, own):
    cls = sub()
    defined = _methods(cls)
    assert defined - overrides - own == set(), "list it as an override or as the subclass's own"
    assert not {n for n in own if hasattr(base, n)}, "an own method that overrides the base"
