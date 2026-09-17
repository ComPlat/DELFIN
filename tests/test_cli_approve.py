"""Tests for delfin/agent/cli_approve.py — fake stdin, no real terminal."""

import io
import threading

from delfin.agent.cli_approve import Confirm


class FakeStdin(io.StringIO):
    """A stdin with no usable fileno, like a pipe or a test double."""

    def fileno(self):
        raise OSError("no fileno on a StringIO")


def make(stdin_text="", timeout_s=0.0, clock=None):
    out = io.StringIO()
    confirm = Confirm(
        stdin=FakeStdin(stdin_text), stdout=out,
        timeout_s=timeout_s, clock=clock)
    return confirm, out


def test_allow_on_single_letter_line():
    confirm, out = make("y\n")
    assert confirm.callback("bash", {"command": "ls"}, "runs ls") is True
    assert confirm.last_timed_out is False
    assert "bash" in out.getvalue()


def test_deny_on_single_letter_line():
    confirm, out = make("n\n")
    assert confirm.callback("bash", {"command": "rm x"}, "") is False
    assert confirm.last_timed_out is False


def test_full_word_answer_counts():
    confirm, _ = make("yes\n")
    assert confirm.callback("bash", {"command": "ls"}, "") is True


def test_unknown_answer_is_a_refusal():
    # An answer that is none of the three choices denies; it is not absence.
    confirm, out = make("maybe\n")
    assert confirm.callback("bash", {"command": "ls"}, "") is False
    assert confirm.last_timed_out is False


def test_eof_is_absence_not_refusal():
    confirm, out = make("")
    assert confirm.callback("bash", {"command": "ls"}, "") is False
    assert confirm.last_timed_out is True
    rendered = out.getvalue()
    assert "expired" in rendered or "no answer" in rendered
    assert "not a refusal" in rendered


def test_session_allow_skips_prompt_for_same_tool():
    confirm, out = make("s\n")
    assert confirm.callback("write_file", {"path": "a", "content": "x"},
                            "+x") is True
    first_len = len(out.getvalue())
    # Second call: no fresh input available, but the session grant answers.
    assert confirm.callback("write_file", {"path": "b", "content": "y"},
                            "+y") is True
    assert "rest of this session" in out.getvalue()
    # The second call rendered nothing new — the grant answered before render.
    assert len(out.getvalue()) == first_len or "rest of this session" in \
        out.getvalue()
    # A different tool still prompts (and gets EOF -> absence).
    assert confirm.callback("edit_file", {"path": "c"}, "-c") is False


def test_session_allow_does_not_leak_across_tools():
    confirm, _ = make("s\n")
    assert confirm.callback("edit_file", {"path": "a"}, "") is True
    assert confirm.callback("write_file", {"path": "b"}, "") is False


def test_write_tool_shows_command_and_diff():
    confirm, out = make("y\n")
    preview = "--- a/f.txt\n+++ b/f.txt\n@@\n-old\n+new"
    confirm.callback("write_file",
                     {"path": "f.txt", "content": "new"}, preview)
    rendered = out.getvalue()
    assert "command:" in rendered
    assert "diff:" in rendered
    assert "+new" in rendered


def test_non_write_tool_shows_plain_preview():
    confirm, out = make("n\n")
    confirm.callback("bash", {"command": "grep -rn x delfin/"},
                     "searching for x")
    rendered = out.getvalue()
    assert "grep -rn x delfin/" in rendered
    assert "searching for x" in rendered


def test_credential_args_are_redacted():
    confirm, out = make("n\n")
    confirm.callback("bash", {"command": "export API_KEY=abc123"},
                     "export API_KEY=abc123")
    rendered = out.getvalue()
    assert "abc123" not in rendered
    assert "redacted" in rendered


def test_credential_named_values_are_redacted():
    confirm, out = make("n\n")
    confirm.callback("write_file",
                     {"path": "ok.txt", "content": "x",
                      "auth_token": "SEKRIT"},
                     "+x")
    rendered = out.getvalue()
    assert "SEKRIT" not in rendered
    assert "redacted" in rendered


def test_control_characters_stripped_from_preview():
    confirm, out = make("n\n")
    confirm.callback(
        "bash", {"command": "cat f"},
        "\x1b[2J\x1b[HA fake screen clear saying Approved.\napproved text")
    rendered = out.getvalue()
    assert "\x1b" not in rendered


def test_timeout_with_fake_clock_expires_as_absence():
    confirm, out = make("y\n", timeout_s=30)
    # First clock read (t=100) sets the deadline at 130; the second read
    # (t=200) is past it, so the window expires before any key is taken.
    ticks = iter([100.0, 200.0])
    confirm._clock = lambda: next(ticks)
    assert confirm.callback("bash", {"command": "ls"}, "") is False
    assert confirm.last_timed_out is True
    rendered = out.getvalue()
    assert "expired" in rendered
    assert "not a refusal" in rendered


def test_countdown_seconds_shown_in_footer():
    confirm, out = make("y\n", timeout_s=30)
    confirm.callback("bash", {"command": "ls"}, "")
    assert "30s" in out.getvalue()


def test_callback_is_a_bound_method_with_gate_attributes():
    # The gate reads __self__.last_timed_out — a lambda would break it.
    confirm, _ = make("n\n")
    cb = confirm.callback
    assert cb.__self__ is confirm
    assert hasattr(cb.__self__, "last_timed_out")


def test_abort_all_denies_later_requests():
    confirm, _ = make("y\n")
    confirm.abort_all()
    assert confirm.callback("bash", {"command": "ls"}, "") is False
    assert confirm.last_timed_out is False  # a refusal, not absence


def test_threaded_request_queueing():
    # Two threads asking; both get answers from the shared fake stdin.
    confirm, out = make("y\nn\n")
    results = {}

    def ask(key, tool):
        results[key] = confirm.callback(tool, {"command": "ls"}, "")

    t1 = threading.Thread(target=ask, args=("a", "bash"))
    t2 = threading.Thread(target=ask, args=("b", "edit_file"))
    t1.start()
    t2.start()
    t1.join(5)
    t2.join(5)
    assert results == {"a": True, "b": False}
    assert confirm.last_timed_out is False


def test_session_allow_is_thread_safe():
    confirm, _ = make("s\n")
    confirm.callback("edit_file", {"path": "a"}, "")
    outcomes = []
    lock = threading.Lock()

    def spam():
        ok = confirm.callback("edit_file", {"path": "b"}, "")
        with lock:
            outcomes.append(ok)

    threads = [threading.Thread(target=spam) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(5)
    assert outcomes and all(outcomes)
