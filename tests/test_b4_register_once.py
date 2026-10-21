"""Red control for the register_once helper (B4 open item: no caller).

A turn that ends blocked IS a wait — but register() would append a new
row on every blocked turn, stacking duplicates. register_once must not.
"""
import pytest

from delfin.agent import state_paths
from delfin.agent import wait_registry


@pytest.fixture
def home(tmp_path):
    with state_paths.RedirectedUserState(tmp_path):
        yield tmp_path / ".delfin"


def test_register_once_exists_and_registers():
    w = wait_registry.register_once(
        "message_from", "operator",
        on_whom="operator", for_what="waiting for the operator",
        key="sess-a")
    assert w["kind"] == "message_from"
    assert w["on_whom"] == "operator"


def test_register_once_does_not_stack_duplicates():
    for _ in range(3):
        wait_registry.register_once(
            "message_from", "operator",
            on_whom="operator", for_what="same wait",
            key="sess-b")
    waits = wait_registry.SessionWaits("sess-b").list()
    assert len(waits) == 1, waits


def test_register_once_allows_a_second_distinct_wait():
    wait_registry.register_once(
        "message_from", "operator", on_whom="operator",
        for_what="one", key="sess-c")
    wait_registry.register_once(
        "time", "9999999999", on_whom="person",
        for_what="deadline", key="sess-c")
    waits = wait_registry.SessionWaits("sess-c").list()
    assert len(waits) == 2, waits


def test_register_once_differs_only_in_for_what_still_dedupes():
    # The identity is (kind, value, on_whom): a changed for_what is the
    # same wait reworded, not a second one.
    wait_registry.register_once(
        "message_from", "operator", on_whom="operator",
        for_what="a", key="sess-d")
    wait_registry.register_once(
        "message_from", "operator", on_whom="operator",
        for_what="b", key="sess-d")
    waits = wait_registry.SessionWaits("sess-d").list()
    assert len(waits) == 1, waits
