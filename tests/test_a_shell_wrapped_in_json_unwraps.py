"""A shell is what a shell says, not the box it came in.

Measured on 2026-09-20, in a real session driven through a pty:

    ⏺ kit-coding:bash  find delfin/agent -name '*.py' | wc -l
      ⎿ 1 lines, 190 B
          {"exit_code": 0, "elapsed_s": 0.095, "stdout": …  under delfin/agent", "cwd": "."}

The answer — 136 — is nowhere on the screen. The tool returned its
result as a ONE-LINE JSON envelope whose carrying field is ``stdout``;
``_result_tail`` keeps the last three lines, a one-line result has one,
and ``truncate_middle`` cuts its middle — which is the output. So the
line that is watched during the turn loses exactly the answer, while
``/trace`` shows the same payload and truncates from the back.

This is not data loss — the model reads the full payload — which makes
it sharper, not harmless: whoever wants to know what the tool said must
read the model's prose and trust it.

The fix decides by CONTENT, not by tool name: a one-line JSON object
with exactly one carrying field (``stdout`` — the shape bash tools
return; ``output``/``text``/``result`` — the shape native tools return)
is unwrapped and the CARRYING FIELD is treated as the output. Any other
shape keeps today's behaviour, so a provider whose envelope looks
different is not silently mangled.
"""

from __future__ import annotations

import json

from delfin.agent.repl_render import tool_result_line

#: The shape a bash tool returns: the payload rides in "stdout".
ENVELOPE_STDOUT = json.dumps({
    "exit_code": 0, "elapsed_s": 0.095,
    "stdout": "136\n", "cwd": ".",
})

#: The shape a native (non-shell) tool returns when it wraps at all.
ENVELOPE_OUTPUT = json.dumps({"ok": True, "output": "x y z"})


def _plain(text: str) -> list[str]:
    import re
    return [re.sub(r"\x1b\[[0-9;]*m", "", ln)
            for ln in (text or "").splitlines() if ln.strip()]


def _screen(output: str, meta: dict | None = None, width: int = 80) -> list[str]:
    return _plain(tool_result_line(
        "kit-coding:bash", output, meta=meta if meta is not None
        else {"ok": True, "chars": len(output)}, width=width))


def test_a_wrapped_stdout_is_on_the_screen():
    """The one number the screen was missing: the carrying field's own
    lines, not the envelope's byte count."""
    lines = _screen(ENVELOPE_STDOUT)
    assert "136" in "\n".join(lines), lines
    assert not any("190 B" in ln or "lines," in ln for ln in lines), lines


def test_a_wrapped_result_keeps_multi_line_tail_behaviour():
    """A carrying field with real lines keeps the tail discipline: the
    last three lines, cut to the width — not the envelope's one."""
    payload = "\n".join(f"row {i}" for i in range(1, 8))
    out = json.dumps({"exit_code": 0, "stdout": payload})
    lines = _screen(out, meta={"ok": True, "chars": len(out)})
    joined = "\n".join(lines)
    assert "row 7" in joined and "row 6" in joined
    assert "row 1" not in joined          # earlier rows are still gone
    assert all(len(ln) <= 80 for ln in lines)


def test_a_json_object_that_is_not_a_shell_is_untouched():
    """Decide by shape, not by tool name: multi-field JSON that carries
    no output-like field (a search result, a summary table) keeps the
    count-and-tail behaviour it has today."""
    out = json.dumps({"matches": 3, "files": ["a.py", "b.py"]})
    lines = _screen(out)
    # "untouched" means the count header stays and the object appears as
    # the tail, verbatim -- today's behaviour -- rather than being
    # unwrapped into its fields (which would print no header at all).
    assert any("1 line," in ln for ln in lines), lines
    assert not any(ln.startswith("      matches") for ln in lines), lines


def test_a_json_array_stays_an_envelope():
    arrays = '["a", "b", "c"]'
    lines = _screen(arrays)
    assert any("1 line," in ln for ln in lines), lines


def test_a_single_field_that_is_not_an_output_is_not_unwrapped():
    """``exit_code`` alone is a flag, not a payload. Unwrapping it would
    print a bare ``0`` — less than nothing."""
    lines = _screen('{"exit_code": 0}')
    assert any("1 line," in ln for ln in lines), lines


def test_a_native_style_envelope_unwraps_too():
    lines = _screen(ENVELOPE_OUTPUT)
    assert "x y z" in "\n".join(lines), lines


def test_an_unwrapped_result_is_not_json_so_it_stays_plain():
    """A plain (non-JSON) one-line output must not gain or lose
    anything: still the last line, still the count."""
    lines = _screen("136\n")
    assert "136" in "\n".join(lines), lines
    assert any("1 line" in ln for ln in lines), lines


def test_pretty_printed_json_is_not_a_one_line_envelope():
    """An envelope the provider pretty-printed over several lines is
    today's case and stays today's case — the shape test is
    deliberately about the ONE-line form, the only one whose single
    line is the middle ``truncate_middle`` destroys."""
    pretty = json.dumps({"exit_code": 0, "stdout": "136\n"}, indent=2)
    lines = _screen(pretty)
    # multi-line: the existing tail behaviour shows its last lines
    assert "136" in "\n".join(lines), lines


def test_a_failed_envelope_is_still_the_blocked_line():
    out = json.dumps({"exit_code": 1, "stdout": "boom"})
    line = tool_result_line(
        "kit-coding:bash", out, meta={"ok": False, "error": "no"},
        width=80)
    assert "blocked" in _plain(line)[0]
    assert "stdout" not in _plain(line)[0]


def test_a_wrapped_result_cannot_smuggle_control_sequences():
    dirty = json.dumps({"exit_code": 0,
                        "stdout": "ok\x1b[2J\x1b]0;pwned\x07"})
    rendered = tool_result_line(
        "kit-coding:bash", dirty, meta={"ok": True}, width=80)
    assert "\x1b" not in rendered
