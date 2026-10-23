"""When a turn ends, the next prompt is already waiting, in grey.

The terminal has had this since 2026-09-19: the open tasks the agent
itself keeps are offered under the answer, and Tab on an empty line takes
the top one. It fills the line; sending stays the user's.

The dashboard had the same suggestions as BUTTONS beside the box. Good,
and not the same thing — a button is a thing to find, and the box is
where the eye already is. The browser's own word for "grey text that
disappears the moment you type" is the placeholder, which this field was
already using for a fixed hint.

So after a turn the placeholder IS the suggestion, verbatim. Verbatim
matters: Tab copies the placeholder into the value, and a prefix like
"Next: " or a trailing "(Tab)" would have to be cut back off in
JavaScript. String surgery in two languages over one value is how the
two halves drift apart. When there is nothing open, the old hint comes
back — an empty grey line would say the agent had nothing to propose,
which is not the same as having nothing to say.

Same collector as the buttons and as the terminal, for the same reason
the buttons already share it: a second reading of the task list would
drift, and this codebase has paid for that twice.
"""

from __future__ import annotations

import pytest

from delfin.dashboard.tab_agent import input_placeholder


HINT = ("Message the agent... (Enter = send, Shift+Enter = newline)\n"
        "KIT tip: say 'also work in /path' to grant write access.")


def test_the_top_suggestion_becomes_the_grey_line():
    got = input_placeholder(["run the full suite on the rebased chain",
                             "push the branch and dispatch CI"], HINT)
    assert got == "run the full suite on the rebased chain"


def test_it_is_verbatim_so_tab_can_copy_it():
    """No prefix, no suffix: what Tab pastes is what the agent meant."""
    step = "measure the refusal corpus again"
    assert input_placeholder([step], HINT) == step


def test_nothing_open_brings_the_hint_back():
    assert input_placeholder([], HINT) == HINT
    assert input_placeholder(None, HINT) == HINT


def test_a_blank_suggestion_is_not_a_suggestion():
    assert input_placeholder(["", "   "], HINT) == HINT


def test_a_very_long_step_is_cut_to_one_readable_line():
    long = "x" * 400
    got = input_placeholder([long], HINT)
    assert len(got) < len(long)
    assert got.endswith("…")


def test_a_newline_would_break_the_grey_line():
    """A placeholder is one line in every browser that matters; a step
    carrying a newline would show half of itself."""
    got = input_placeholder(["do this\nthen that"], HINT)
    assert "\n" not in got
    assert "do this" in got and "then that" in got


def test_it_never_raises():
    for junk in (["ok"], [], None, ("tuple",), [None], [123]):
        input_placeholder(junk, HINT)


def test_the_tab_key_takes_it_and_does_not_send():
    """Driven at the source the browser actually runs: the handler has to
    copy the placeholder through the native setter, or ipywidgets never
    sees the value and Send would post an empty box."""
    import inspect

    from delfin.dashboard import tab_agent as T

    src = inspect.getsource(T)
    # Up to the NEXT key branch, not a fixed window: the Enter branch
    # right below does click Send, and a window that runs into it would
    # fail this for the wrong reason.
    branch = src.split("e.key === 'Tab'")[1].split("if (e.key ===")[0]
    assert "placeholder" in branch, branch[:300]
    assert "HTMLTextAreaElement.prototype" in branch, (
        "setting .value directly leaves ipywidgets unaware of it")
    assert "new Event('input'" in branch
    assert "send-row" not in branch, (
        "Tab fills the box; sending stays the user's")


def test_the_buttons_and_the_grey_line_read_the_same_list():
    import inspect

    from delfin.dashboard import tab_agent as T

    src = inspect.getsource(T)
    assert src.count("from delfin.agent.task_ticker import next_steps") == 1
    assert "input_placeholder(" in src, (
        "the grey line must come from the same refresh as the buttons")
