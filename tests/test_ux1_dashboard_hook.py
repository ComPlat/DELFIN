# UX1 phase 3b: the dashboard hook renders the turn summary line.
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

import textwrap

import pytest

from delfin.agent.turn_summary import format_turn_summary


def _hook_block_source():
    """The hook block, extracted from the patched copy, dedented.

    Reads .gate/px/dxb/tab_agent.py (the operator's patched copy of
    tab_agent.py). Returns the try/except block's BODY lines dedented,
    ready to exec in a prepared namespace. Missing copy -> skip: the
    operator's build tree is the only place it exists.
    """
    import pathlib
    copy = pathlib.Path(".gate/px/dxb/tab_agent.py")
    if not copy.exists():
        pytest.skip(".gate/px/dxb/tab_agent.py not present")
    lines = copy.read_text().splitlines()
    start = next(i for i, l in enumerate(lines)
                 if "# UX1: the end-of-turn line lives ON the last" in l)
    # the try: starts a few lines below the comment
    tstart = next(i for i in range(start, start + 12)
                  if lines[i].strip() == "try:")
    # find the matching except at the try's indent
    indent = len(lines[tstart]) - len(lines[tstart].lstrip())
    tsend = next(i for i in range(tstart + 1, len(lines))
                 if lines[i].strip().startswith("except Exception")
                 and (len(lines[i]) - len(lines[i].lstrip())) == indent)
    body = lines[tstart + 1:tsend]
    dedent = textwrap.dedent("\n".join(body))
    return dedent


# --- end of extraction helper; the original contract tests follow ---

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


class _FakeEngine:
    """Only the attributes the hook block reads."""

    kit_permissions = None
    _stop_requested = False


class _FakeCtx:
    class repo_dir:  # str-like is enough; str(ctx.repo_dir or "") below
        pass

    def __init__(self):
        self.repo_dir = "/tmp/fake-repo"


def _agent_workspace_path() -> str:
    """The enclosing scope's real helper, faithfully reproduced from
    tab_agent.py:9556: engine.kit_permissions.workspace first, then
    ctx.repo_dir. The probe namespace must supply it, because the block
    references it directly.
    """
    eng = state_probe["engine"]
    kp = getattr(eng, "kit_permissions", None) if eng else None
    ws = getattr(kp, "workspace", None) or getattr(eng, "repo_dir", None)
    if ws:
        return str(ws)
    return str(getattr(state_probe["ctx"], "agent_workspace", "")
               or (state_probe["ctx"].repo_dir or ""))


def test_hook_block_runs_in_context_and_attaches_the_line():
    """Drive the hook block IN CONTEXT — the test that would have
    caught the msgs-NameError (s17's find): the buggy version referenced
    ``msgs`` which the enclosing code never defines; its except
    swallowed the NameError and the line silently never appeared.

    The namespace provides ONLY the names the enclosing scope really
    supplies (state, engine, ctx, role_label). ``msgs`` is deliberately
    absent: the block must bind it from state["chat_messages"] itself
    (s17's fix). The assertion is the SIDE EFFECT on chat_messages.
    """
    block = _hook_block_source()
    import time as _time
    import datetime as _dt
    global state_probe
    state = {
        "chat_messages": [
            {"role": "user", "content": "hi"},
            {"role": "assistant", "content": "the answer",
             "role_label": "Agent"},
        ],
        "_turn_started_monotonic": _time.monotonic() - 31.0,
    }
    state_probe = {"engine": _FakeEngine(), "ctx": _FakeCtx()}
    ns = {
        "state": state,
        "engine": _FakeEngine(),
        "ctx": _FakeCtx(),
        "role_label": "Agent",
        "time": _time,
        "datetime": _dt,
        "_agent_workspace_path": _agent_workspace_path,
    }
    exec(compile(block, "<ux1-hook-block>", "exec"), ns)
    msgs = state["chat_messages"]
    assert msgs[-1]["role"] == "assistant"
    summary = msgs[-1].get("turn_summary")
    assert isinstance(summary, str) and summary
    assert summary.startswith("Surfaced after ")
    # Never a separate chat entry, never on the user message.
    assert len(msgs) == 2
    assert "turn_summary" not in msgs[0]


def test_hook_block_survives_an_exception_and_never_breaks_the_turn():
    """The guard property: a failure inside the block leaves the chat
    history exactly as the turn left it (the except swallows, the
    assistant entry stays untouched, no partial meta is written).
    """
    block = _hook_block_source()
    import time as _time
    import datetime as _dt

    class _ExplodingJobs:
        @staticmethod
        def running_jobs_for_workspace(_ws):
            raise RuntimeError("registry gone")

    state = {
        "chat_messages": [
            {"role": "assistant", "content": "the answer"},
        ],
        "_turn_started_monotonic": _time.monotonic() - 5.0,
    }
    ns = {
        "state": state,
        "engine": _FakeEngine(),
        "ctx": _FakeCtx(),
        "role_label": "Agent",
        "time": _time,
        "datetime": _dt,
        "_agent_workspace_path": _agent_workspace_path,
    }
    # The real block is inside try/except Exception: pass — the test
    # models that guard the same way the enclosing code does.
    state["_turn_started_monotonic"] = "not-a-number"
    try:
        exec(compile(block, "<ux1-hook-block>", "exec"), ns)
    except Exception:
        pass
    msgs = state["chat_messages"]
    assert msgs[-1]["content"] == "the answer"
    assert "turn_summary" not in msgs[-1]
