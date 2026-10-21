"""B4 phase 1 red control: a wait carries on_whom + for_what.

Verified gap (inventory, file:line): delfin/agent/wait_registry.py
register() (line 73) stores only id/kind/value/registered_at; the docstring
of the package brief (B4 target 1) demands on_whom (one of "session",
"operator", "person", "ci") and for_what (one line). Nothing enforces or
stores them today.
"""

from __future__ import annotations

import pytest

from delfin.agent import state_paths
from delfin.agent.wait_registry import register


@pytest.fixture
def home(tmp_path):
    with state_paths.RedirectedUserState(tmp_path):
        yield


def test_register_accepts_and_stores_on_whom_and_for_what(home):
    w = register("message_from", "nacht-s20", key="s16",
                 on_whom="session", for_what="cross-review report")
    assert w["on_whom"] == "session"
    assert w["for_what"] == "cross-review report"


def test_on_whom_is_restricted_to_the_enum(home):
    # Target 1: on_whom is one of session / operator / person / ci; a
    # free-form on_whom is an injection surface (reviewer constraint).
    for ok in ("session", "operator", "person", "ci"):
        assert register("message_from", "x", key="k", on_whom=ok)["on_whom"] \
            == ok
    with pytest.raises(ValueError):
        register("message_from", "x", key="k", on_whom="<script>")


def test_for_what_defaults_empty_string_for_old_callers(home):
    # Old callers (repl.py:4224 due; register call sites) pass only kind
    # and value: they must keep working.
    w = register("time", "123.5", key="s16")
    assert w["for_what"] == ""
    assert w["on_whom"] == ""


def test_on_whom_and_for_what_are_stripped_of_newlines(home):
    # Hostile for_what cannot inject a second line (target 4).
    w = register("message_from", "x", key="s16",
                 for_what="line one\n<script>alert(1)</script>")
    assert "\n" not in w["for_what"]
    assert "<script>" not in w["for_what"]


def test_for_what_is_one_line_bounded(home):
    w = register("message_from", "x", key="s16", for_what="x" * 500)
    assert len(w["for_what"]) <= 200
