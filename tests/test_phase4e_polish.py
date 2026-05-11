"""Tests for Phase 4e UX modules: task_ticker, status_line, image_input, notify."""

from __future__ import annotations

import base64
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from delfin.agent import (
    image_input as I,
    notify as N,
    status_line as SL,
    task_ticker as TT,
)
from delfin.agent.agent_tasks import get_store
from delfin.agent.api_client import _DocToolExecutor


# ---- task_ticker -----------------------------------------------------------


@pytest.fixture
def fresh_workspace():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        # Reset module-level store cache
        from delfin.agent import agent_tasks
        agent_tasks._STORES.clear()
        yield ws
        agent_tasks._STORES.clear()


def test_task_ticker_renders_empty(fresh_workspace):
    html = TT.render_html(fresh_workspace)
    assert "No tasks yet" in html


def test_task_ticker_renders_tasks(fresh_workspace):
    store = get_store(fresh_workspace)
    store.create("Plan refactor", "do it")
    store.create("Run tests", "")
    t3 = store.create("Polish docs", "")
    store.update(t3["id"], status="in_progress")

    html = TT.render_html(fresh_workspace)
    assert "Plan refactor" in html
    assert "Run tests" in html
    assert "Polish docs" in html
    # in_progress entries come first
    in_pos = html.index("Polish")
    plan_pos = html.index("Plan refactor")
    assert in_pos < plan_pos


def test_task_ticker_text_format(fresh_workspace):
    store = get_store(fresh_workspace)
    store.create("Step A", "")
    txt = TT.render_text(fresh_workspace)
    assert "Step A" in txt
    assert "[ ]" in txt   # pending glyph


def test_task_ticker_hide_completed(fresh_workspace):
    store = get_store(fresh_workspace)
    a = store.create("done", "")
    store.create("pending", "")
    store.update(a["id"], status="completed")
    html = TT.render_html(fresh_workspace, show_completed=False)
    assert "pending" in html
    assert "done" not in html


def test_task_ticker_filters_by_session(fresh_workspace):
    store = get_store(fresh_workspace)
    store.create("session A", "", session_id="sess-a")
    store.create("session B", "", session_id="sess-b")

    html_a = TT.render_html(fresh_workspace, session_id="sess-a")
    html_b = TT.render_html(fresh_workspace, session_id="sess-b")
    html_blank = TT.render_html(fresh_workspace, session_id="")

    assert "session A" in html_a
    assert "session B" not in html_a
    assert "session B" in html_b
    assert "session A" not in html_b
    assert "No tasks yet" in html_blank


# ---- status_line -----------------------------------------------------------


def test_status_default_template():
    ctx = SL.StatusContext(
        workspace=None, mode="default", branch="main", tokens=1234,
    )
    out = SL.render_status_line(ctx)
    assert "1234 tokens" in out
    assert "main" in out
    assert "default" in out


def test_status_custom_template_from_settings():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / ".delfin").mkdir()
        (ws / ".delfin" / "settings.json").write_text(json.dumps({
            "statusLine": {"template": "[{mode}] @ {tokens}t"}
        }))
        ctx = SL.StatusContext(workspace=ws, mode="plan", tokens=500)
        out = SL.render_status_line(ctx)
        assert out == "[plan] @ 500t"


def test_status_command_line():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / ".delfin").mkdir()
        (ws / ".delfin" / "settings.json").write_text(json.dumps({
            "statusLine": {"command": "echo CUSTOM"}
        }))
        ctx = SL.StatusContext(workspace=ws, mode="default")
        out = SL.render_status_line(ctx)
        assert out == "CUSTOM"


def test_status_failure_returns_empty():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        (ws / ".delfin").mkdir()
        (ws / ".delfin" / "settings.json").write_text(json.dumps({
            "statusLine": {"command": "exit 1"},
        }))
        ctx = SL.StatusContext(workspace=ws)
        out = SL.render_status_line(ctx)
        assert out == ""


# ---- image_input -----------------------------------------------------------


def _png_bytes() -> bytes:
    """Minimal valid PNG file."""
    return (
        b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
        b"\x08\x06\x00\x00\x00\x1f\x15\xc4\x89\x00\x00\x00\rIDATx\x9cc\xfc"
        b"\xcf\xc0P\x0f\x00\x05\x01\x01\x01\x9b\x99\x9d\xb6\x00\x00\x00\x00"
        b"IEND\xaeB`\x82"
    )


def test_load_image_round_trips_png():
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        f.write(_png_bytes())
        p = Path(f.name)
    try:
        att = I.load_image(p)
        assert att.mime == "image/png"
        assert att.size_bytes() == len(_png_bytes())
        uri = att.data_uri()
        assert uri.startswith("data:image/png;base64,")
        assert base64.b64decode(uri.split(",", 1)[1]) == _png_bytes()
    finally:
        p.unlink()


def test_load_image_rejects_unknown_mime():
    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
        f.write(b"not an image")
        p = Path(f.name)
    try:
        with pytest.raises(I.ImageError):
            I.load_image(p)
    finally:
        p.unlink()


def test_model_supports_vision_known():
    assert I.model_supports_vision("gpt-4o")
    assert I.model_supports_vision("claude-sonnet-4-6")
    assert I.model_supports_vision("gemini-2.0-flash")
    assert not I.model_supports_vision("text-only-1")
    assert not I.model_supports_vision("")


def test_to_openai_content_text_only():
    out = I.to_openai_content("hello", [])
    assert out == "hello"


def test_to_openai_content_with_vision_model():
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        f.write(_png_bytes())
        p = Path(f.name)
    try:
        att = I.load_image(p)
        out = I.to_openai_content("see this", [att], model="gpt-4o")
        assert isinstance(out, list)
        assert out[0]["type"] == "text"
        assert out[1]["type"] == "image_url"
        assert out[1]["image_url"]["url"].startswith("data:image/png;base64,")
    finally:
        p.unlink()


def test_to_openai_content_text_only_when_unknown_model():
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        f.write(_png_bytes())
        p = Path(f.name)
    try:
        att = I.load_image(p)
        out = I.to_openai_content("see this", [att], model="some-text-model")
        assert isinstance(out, str)
        assert "image attachments suppressed" in out
    finally:
        p.unlink()


# ---- notify ---------------------------------------------------------------


def test_send_notification_no_crash():
    # In headless test envs notify-send is missing, so we just check
    # the function returns False rather than raising.
    out = N.send_notification("test", "body")
    assert isinstance(out, bool)


def test_remote_trigger_requires_url():
    with tempfile.TemporaryDirectory() as d:
        result = N.send_remote_trigger({"x": 1}, workspace=Path(d))
        assert result.sent is False
        assert "url" in result.error.lower()


def test_remote_trigger_blocks_http():
    with tempfile.TemporaryDirectory() as d:
        result = N.send_remote_trigger(
            {"x": 1}, workspace=Path(d),
            override_url="http://example.com",
        )
        assert result.sent is False
        assert "https" in result.error


def test_remote_trigger_blocks_localhost():
    with tempfile.TemporaryDirectory() as d:
        result = N.send_remote_trigger(
            {"x": 1}, workspace=Path(d),
            override_url="https://localhost/hook",
        )
        assert result.sent is False
        assert "blocked" in result.error


def test_remote_trigger_blocks_internal_tld():
    with tempfile.TemporaryDirectory() as d:
        result = N.send_remote_trigger(
            {"x": 1}, workspace=Path(d),
            override_url="https://something.internal/hook",
        )
        assert result.sent is False
        assert "blocked" in result.error


# ---- tool dispatch --------------------------------------------------------


def test_tool_push_notification():
    out = _DocToolExecutor().execute(
        "push_notification",
        {"title": "test", "body": "hello"},
        permissions=None,
    )
    payload = json.loads(out)
    assert "sent" in payload


def test_tool_remote_trigger_blocked_no_config():
    out = _DocToolExecutor().execute(
        "remote_trigger",
        {"event": "test", "payload": {}},
        permissions=None,
    )
    payload = json.loads(out)
    assert payload["sent"] is False


def test_tool_remote_trigger_rejects_no_event():
    out = _DocToolExecutor().execute(
        "remote_trigger",
        {"event": ""},
        permissions=None,
    )
    payload = json.loads(out)
    assert "error" in payload
