"""Tests for the ask_user_question tool dispatch + executor."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from delfin.agent.api_client import (
    KitToolPermissions, _DocToolExecutor,
)


def _make_perms(**kwargs) -> KitToolPermissions:
    tmp = tempfile.mkdtemp(prefix="askuser_")
    return KitToolPermissions(workspace=Path(tmp), mode="default", **kwargs)


def _exec(arguments, perms):
    return _DocToolExecutor().execute(
        "ask_user_question", arguments, permissions=perms,
    )


def test_returns_error_without_callback():
    perms = _make_perms()  # no ask_user_callback
    out = _exec({
        "question": "x?",
        "options": [{"label": "a"}, {"label": "b"}],
    }, perms)
    payload = json.loads(out)
    assert "error" in payload
    assert "not available" in payload["error"]


def test_returns_error_with_too_few_options():
    perms = _make_perms(ask_user_callback=lambda _: {"answers": ["a"]})
    out = _exec({
        "question": "x?",
        "options": [{"label": "only-one"}],
    }, perms)
    payload = json.loads(out)
    assert "error" in payload


def test_returns_error_with_empty_question():
    perms = _make_perms(ask_user_callback=lambda _: {"answers": ["a"]})
    out = _exec({
        "question": "  ",
        "options": [{"label": "a"}, {"label": "b"}],
    }, perms)
    payload = json.loads(out)
    assert "error" in payload


def test_callback_invoked_with_normalised_args():
    captured: dict = {}

    def cb(args):
        captured.update(args)
        return {"answers": ["yes"]}

    perms = _make_perms(ask_user_callback=cb)
    out = _exec({
        "question": "delete?",
        "header": "danger",
        "options": [
            {"label": "yes", "description": "do it"},
            {"label": "no"},
        ],
        "multiSelect": False,
    }, perms)
    payload = json.loads(out)
    assert payload["answers"] == ["yes"]
    assert payload["multiSelect"] is False
    assert captured["question"] == "delete?"
    assert captured["header"] == "danger"
    assert len(captured["options"]) == 2
    assert captured["options"][0] == {"label": "yes", "description": "do it"}


def test_multi_select_passed_through():
    perms = _make_perms(
        ask_user_callback=lambda _: {"answers": ["a", "b"]},
    )
    out = _exec({
        "question": "pick many?",
        "options": [
            {"label": "a"}, {"label": "b"}, {"label": "c"},
        ],
        "multiSelect": True,
    }, perms)
    payload = json.loads(out)
    assert payload["answers"] == ["a", "b"]
    assert payload["multiSelect"] is True


def test_single_select_truncates_multi_answer():
    perms = _make_perms(
        ask_user_callback=lambda _: {"answers": ["a", "b"]},
    )
    out = _exec({
        "question": "single only?",
        "options": [
            {"label": "a"}, {"label": "b"},
        ],
        "multiSelect": False,
    }, perms)
    payload = json.loads(out)
    assert payload["answers"] == ["a"]


def test_callback_exception_returns_error():
    def bad(_args):
        raise RuntimeError("boom")
    perms = _make_perms(ask_user_callback=bad)
    out = _exec({
        "question": "x?",
        "options": [{"label": "a"}, {"label": "b"}],
    }, perms)
    payload = json.loads(out)
    assert "error" in payload
    assert "boom" in payload["error"]


def test_invalid_callback_return_caught():
    perms = _make_perms(ask_user_callback=lambda _: "not-a-dict")
    out = _exec({
        "question": "x?",
        "options": [{"label": "a"}, {"label": "b"}],
    }, perms)
    payload = json.loads(out)
    assert "error" in payload


def test_option_without_label_rejected():
    perms = _make_perms(ask_user_callback=lambda _: {"answers": []})
    out = _exec({
        "question": "x?",
        "options": [{"label": "a"}, {"description": "no label here"}],
    }, perms)
    payload = json.loads(out)
    assert "error" in payload
