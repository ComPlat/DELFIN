"""One flag for every thread turned an expiry into a refusal.

An expired confirmation must not be recorded as a denial: absence is not
a decision, and a refusal closes the path for the rest of the session.
The brokers say so, and mark it with ``last_timed_out`` — which the gate
reads off the bound callback's ``__self__``.

That flag was ONE value on the broker, and the broker's own docstring
says requests arrive "on whatever thread is executing tools: the turn
worker, a subagent, a background job". So they overlap:

    subagent   asks about /etc/hosts     → nobody answers → expires,
                                           flag = True
    main turn  asks about repl.py        → the user approves,
                                           flag = False
    subagent   reads the flag            → False → "user denied"

and a path the user never saw is closed for the rest of the session.
That is the exact outcome the flag exists to prevent.

The flag is per thread now. It keeps its name and stays an attribute,
because the gate reaches it through ``__self__``; only its storage moved.
"""

from __future__ import annotations

import threading

import pytest


@pytest.mark.parametrize("module_name", ["kit_confirm", "terminal_confirm"])
def test_two_threads_do_not_share_one_expiry(module_name):
    import importlib
    mod = importlib.import_module(f"delfin.agent.{module_name}")
    broker = _make(mod)

    seen: dict = {}
    ready = threading.Event()
    released = threading.Event()

    def _expired_here():
        broker.last_timed_out = True
        ready.set()
        released.wait(5)
        seen["expired_thread"] = broker.last_timed_out

    t = threading.Thread(target=_expired_here, daemon=True)
    t.start()
    ready.wait(5)

    # Another thread answers its own request in the meantime.
    broker.last_timed_out = False
    released.set()
    t.join(5)

    assert seen["expired_thread"] is True, (
        "an answer on one thread cleared another thread's expiry — the "
        "expired request is then recorded as a refusal")
    assert broker.last_timed_out is False, "and this thread keeps its own"


@pytest.mark.parametrize("module_name", ["kit_confirm", "terminal_confirm"])
def test_the_gate_still_reads_it_off_the_bound_method(module_name):
    """The gate does `perms.confirm_callback.__self__.last_timed_out`.
    Wrapping the callback in a lambda breaks that silently, and so would
    moving the flag off the instance."""
    import importlib
    mod = importlib.import_module(f"delfin.agent.{module_name}")
    broker = _make(mod)
    cb = broker.callback
    assert getattr(cb, "__self__", None) is broker
    broker.last_timed_out = True
    assert cb.__self__.last_timed_out is True
    broker.last_timed_out = False
    assert cb.__self__.last_timed_out is False


@pytest.mark.parametrize("module_name", ["kit_confirm", "terminal_confirm"])
def test_a_fresh_thread_starts_undecided(module_name):
    import importlib
    mod = importlib.import_module(f"delfin.agent.{module_name}")
    broker = _make(mod)
    broker.last_timed_out = True
    out: dict = {}

    def _fresh():
        out["value"] = broker.last_timed_out

    t = threading.Thread(target=_fresh, daemon=True)
    t.start()
    t.join(5)
    assert out["value"] is False, (
        "a thread that has asked nothing must not inherit somebody "
        "else's expiry")


def _make(mod):
    for name in ("KitConfirmBroker", "TerminalConfirmBroker"):
        cls = getattr(mod, name, None)
        if cls is not None:
            return cls()
    raise AssertionError(f"no broker in {mod.__name__}")
