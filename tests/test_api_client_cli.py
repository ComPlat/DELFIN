"""Tests for CLIClient empty-text guard.

The Claude CLI tags each user content block with cache_control. An
empty text block triggers an Anthropic API 400:
    "cache_control cannot be set for empty text blocks"
CLIClient.stream_message must skip the send entirely in that case.
"""
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.api_client import CLIClient


@pytest.fixture
def client(tmp_path):
    """Build a CLIClient without spawning a subprocess."""
    c = CLIClient.__new__(CLIClient)
    c.claude_path = "claude"
    c.model = "sonnet"
    c.permission_mode = ""
    c.cwd = str(tmp_path)
    c.mcp_config = ""
    c.allowed_tools = None
    c.extra_dirs = None
    c._proc = None
    c._session_id = ""
    return c


def _drain(gen):
    """Consume a generator and return events as a list."""
    return list(gen)


def test_empty_user_content_skips_send(client):
    """Empty user message must not spawn the CLI subprocess."""
    with patch.object(client, "_ensure_proc") as ensure:
        events = _drain(
            client.stream_message(
                system="sys",
                messages=[{"role": "user", "content": ""}],
            )
        )
    assert events == []
    ensure.assert_not_called()


def test_whitespace_only_user_content_skips_send(client):
    """Whitespace-only user message must not spawn the CLI subprocess."""
    with patch.object(client, "_ensure_proc") as ensure:
        events = _drain(
            client.stream_message(
                system="sys",
                messages=[{"role": "user", "content": "   \n\t  "}],
            )
        )
    assert events == []
    ensure.assert_not_called()


def test_no_user_message_in_history_skips_send(client):
    """A history with only assistant messages must skip the send."""
    with patch.object(client, "_ensure_proc") as ensure:
        events = _drain(
            client.stream_message(
                system="sys",
                messages=[{"role": "assistant", "content": "earlier reply"}],
            )
        )
    assert events == []
    ensure.assert_not_called()


def test_non_string_user_content_treated_as_empty(client):
    """A list-typed content (rare) is treated as empty for safety."""
    with patch.object(client, "_ensure_proc") as ensure:
        events = _drain(
            client.stream_message(
                system="sys",
                messages=[{"role": "user", "content": []}],
            )
        )
    assert events == []
    ensure.assert_not_called()


def test_signal_stop_sends_sigint_to_running_proc(client):
    """signal_stop sends SIGINT to a live subprocess, leaves _proc reference
    in place so the next stream_message can detect the running session."""
    import signal as _signal
    fake_proc = MagicMock()
    fake_proc.poll.return_value = None  # still running
    client._proc = fake_proc
    client.signal_stop()
    fake_proc.send_signal.assert_called_once_with(_signal.SIGINT)
    # Reference preserved (caller owns lifecycle for resume)
    assert client._proc is fake_proc


def test_signal_stop_no_proc_is_safe(client):
    """signal_stop on a client without a subprocess is a no-op."""
    client._proc = None
    client.signal_stop()  # must not raise


def test_signal_stop_on_dead_proc_skips(client):
    """signal_stop ignores already-exited subprocess."""
    fake_proc = MagicMock()
    fake_proc.poll.return_value = 0  # exited
    client._proc = fake_proc
    client.signal_stop()
    fake_proc.send_signal.assert_not_called()


def test_signal_stop_falls_back_to_terminate(client):
    """If SIGINT raises, signal_stop terminates and clears the proc handle."""
    fake_proc = MagicMock()
    fake_proc.poll.return_value = None
    fake_proc.send_signal.side_effect = OSError("EPERM")
    client._proc = fake_proc
    client.signal_stop()
    fake_proc.terminate.assert_called_once()
    assert client._proc is None


def test_base_client_signal_stop_is_noop():
    """APIClient and other backends without a subprocess inherit a no-op."""
    from delfin.agent.api_client import _BaseClient
    base = _BaseClient.__new__(_BaseClient)
    # Must not raise
    base.signal_stop()


def test_non_empty_user_content_writes_stdin(client):
    """A real user message must write a stream-json line to CLI stdin."""
    fake_proc = MagicMock()
    fake_proc.stdin = MagicMock()
    fake_proc.stdout = iter([])  # empty stdout ends _read_turn quickly
    fake_proc.poll.return_value = 0
    fake_proc.stderr.read.return_value = ""
    with patch.object(client, "_ensure_proc", return_value=fake_proc):
        # We don't care about events — just that stdin got the JSON line.
        try:
            list(
                client.stream_message(
                    system="sys",
                    messages=[{"role": "user", "content": "hello"}],
                )
            )
        except RuntimeError:
            pass  # _read_turn may raise on empty stdout, that's fine
    fake_proc.stdin.write.assert_called()
    written = fake_proc.stdin.write.call_args_list[0].args[0]
    assert '"text": "hello"' in written or '"text":"hello"' in written
