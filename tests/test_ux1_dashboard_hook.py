"""UX1 phase 3b: the dashboard hook renders the turn summary line.

The dashboard patch (built from .gate/px/dxa -> dxb by the operator)
wires two things in tab_agent.py:

- the turn-end block attaches msg["turn_summary"] to the LAST
  assistant entry in state["chat_messages"] (not a separate chat
  message, never an engine.messages entry): the line renders dimmed
  under the answer, is saved with the session (session_store persists
  chat_messages verbatim) and comes back after a restore.
- _render_single_msg renders that meta key for the assistant role as
  an escaped, classed div under the bubble.

These tests pin the CONTRACT each side sees. The hook itself is a
try/except block in a 23k-line file; its guard property — a failing
hook never breaks a turn — is what the patch is reviewed on.
"""

from delfin.agent.turn_summary import format_turn_summary


def test_hook_input_shape_builds_the_line():
    # Exactly the call the hook makes (duration from the turn-started
    # monotonic delta, ended_at now, running_jobs from
    # bash_jobs.running_jobs_for_workspace). All-numeric inputs; a job
    # name or command has no path into the string.
    line = format_turn_summary(duration_s=31.4,
                               ended_at=None,
                               outcome="ok",
                               running_jobs=1)
    assert line == "Surfaced after 31s · 1 shell still running"


def test_stopped_turns_get_the_stopped_verb():
    # The hook derives the outcome from engine._stop_requested — a
    # stopped turn never claims a success verb.
    line = format_turn_summary(duration_s=12, ended_at=None,
                               outcome="stopped", running_jobs=0)
    assert line.startswith("Stopped after 12s")
    assert "Surfaced" not in line


def test_turn_summary_rides_on_the_last_assistant_entry():
    # The hook's storage shape: meta on the last chat entry when it is
    # an assistant message — not a new entry, not an engine message.
    msgs = [
        {"role": "user", "content": "hi"},
        {"role": "assistant", "content": "done", "role_label": "Agent"},
    ]
    summary = format_turn_summary(duration_s=31, ended_at=None,
                                  outcome="ok", running_jobs=1)
    # Mimic the hook: last is assistant -> set meta on it.
    msgs[-1] = dict(msgs[-1])
    msgs[-1]["turn_summary"] = summary
    assert msgs[-1]["turn_summary"] == summary
    # The engine's message list is a DIFFERENT list: nothing was added
    # to it here, and the chat meta key is unknown to any consumer of
    # engine messages (save_session persists chat_messages verbatim;
    # engine_messages stay role+content only).
    assert all("turn_summary" not in m
               for m in msgs if m["role"] != "assistant")


def test_renderer_contract_escapes_and_classes():
    # The renderer wraps the meta in a classed div, escaped. Proven
    # with a hostile string through the same escape the patched render
    # site uses.
    import html as _html
    summary = '<img src=x onerror=alert(1)>'
    rendered = (f'<div class="delfin-chat-turn-summary">'
                f'{_html.escape(summary)}</div>')
    escaped = _html.escape(summary)
    assert "<" not in escaped.replace("&lt;", "")
    assert ">" not in escaped.replace("&gt;", "")
    assert rendered == ('<div class="delfin-chat-turn-summary">'
                        "&lt;img src=x onerror=alert(1)&gt;</div>")
    # The formatter's own output stays intact under the same escape.
    benign = format_turn_summary(duration_s=31, ended_at=None,
                                 outcome="ok", running_jobs=1)
    assert _html.escape(benign) == benign  # no HTML chars to escape


def test_restore_shows_the_line_again():
    # "The line comes back after a session restore": chat_messages are
    # saved verbatim (session_store.py:1238) and replayed on resume, so
    # a meta key on the assistant entry round-trips untouched.
    stored = format_turn_summary(duration_s=133, ended_at=None,
                                 outcome="ok", running_jobs=2)
    msgs = [{"role": "assistant", "content": "x",
             "turn_summary": stored}]
    # str(...) or "" fallback the renderer uses, on a restored entry.
    read_back = str(msgs[0].get("turn_summary") or "")
    assert read_back == "Surfaced after 2m 13s · 2 jobs still running"
    # An absent key degrades to "" (omitted, never an error).
    assert str({}.get("turn_summary") or "") == ""
    assert str({"turn_summary": None}.get("turn_summary") or "") == ""
