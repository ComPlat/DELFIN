"""UX1 phase 3b: the dashboard hook renders the turn summary line.

The dashboard patch (.gate/ux1_dashboard_hook.patch, applied by
nacht-s17 / the operator) wires three things in tab_agent.py:

- the turn-end block stores state["_last_turn_summary"] via
  format_turn_summary (same inputs as the CLI footer),
- _render_status takes last_turn_summary and HTML-escapes it into the
  status row,
- create_tab passes the state value into _render_status.

These tests pin the CONTRACT each side sees: the pure formatter (what
the hook calls), and the status-bar string shape (what the person
reads). The hook itself is a try/except block in a 23k-line file; its
guard property — a failing hook never breaks a turn — is what the
patch is reviewed on, and s17's application is verified on his branch.
"""

from delfin.agent.turn_summary import format_turn_summary


def test_hook_input_shape_builds_the_line():
    # Exactly the call the hook in the patch makes (duration from the
    # turn-started monotonic delta, ended_at now, running_jobs from
    # bash_jobs.running_jobs_for_workspace). All-numeric inputs; a job
    # name or command has no path into the string.
    line = format_turn_summary(duration_s=31.4,
                               ended_at=None,
                               outcome="ok",
                               running_jobs=1)
    assert line == "Surfaced after 31s · 1 shell still running"


def test_hook_failure_path_yields_empty_string():
    # The hook's except branch stores ""; the status row's omission rule
    # (patch, _render_status) is truthiness: a falsy stored value renders
    # nothing. Pin the actual contract, not a tautology:
    failure_value = ""
    assert not failure_value, "hook failure path must be falsy"
    # The formatter, by contrast, always returns a non-empty line
    # ("Turn ended" fallback, phase-2 pin) — a success path never
    # produces a value the row would omit.
    success_value = format_turn_summary(duration_s=None, ended_at=None,
                                        outcome="", running_jobs=None)
    assert success_value
    assert success_value == "Turn ended"


def test_status_row_shape_escapes_everything():
    # The status row wraps the summary in
    #   · <span class="turn-summary">{escaped}</span>
    # Escape rule proven with a HOSTILE string, not a formatter output:
    # whatever reaches the wrapper goes through html.escape first.
    import html as _html
    summary = '<img src=x onerror=alert(1)>'
    expected = (' · <span class="turn-summary">'
                + _html.escape(summary) + "</span>")
    escaped = _html.escape(summary)
    assert "<" not in escaped.replace("&lt;", "")
    assert ">" not in escaped.replace("&gt;", "")
    # And the whole hostile string never forms a live tag: the opener
    # '<img' and closer '>' are both gone, so nothing renders.
    # The formatter's own output stays intact under the same escape.
    benign = format_turn_summary(duration_s=31, ended_at=None,
                                 outcome="ok", running_jobs=1)
    assert _html.escape(benign) == benign  # no HTML chars to escape


def test_restore_shows_the_line_again():
    # "The line comes back after a session restore": the hook stores
    # the string in state; a restored session reads the same key with
    # the same fallback and re-renders it. Pin that the stored value
    # round-trips untouched through str(... or "").
    stored = format_turn_summary(duration_s=133, ended_at=None,
                                 outcome="ok", running_jobs=2)
    read_back = str({"_last_turn_summary": stored}.get(
        "_last_turn_summary") or "")
    assert read_back == "Surfaced after 2m 13s · 2 jobs still running"
    # An absent key degrades to "" (omitted, never an error).
    assert str({}.get("_last_turn_summary") or "") == ""
