"""The generated TS->GOAT viewer JavaScript must be structurally whole.

Regression for the P0 white-viewer bug of 2026-11-02: a stick-restyle
removal left dangling ``viewer.setStyle({serial: serial},`` calls inside
three JS generators, so the whole ``(function(){...})()`` failed to
parse -- the viewer stayed blank and no panel message appeared, because
the syntax error only reached the browser console.

Method: a lexer in Python walks each generated script, skipping string,
template, regex and comment bodies, and requires every ``(``, ``[`` and
``{`` to close in order and every literal to terminate. That is the
class of defect the bug was (an unclosed call). Rejected alternative:
``node --check`` behind a skip -- no JS runtime is a DELFIN dependency
and a test must run on every machine
(tests/test_a_test_is_bound_to_delfin_not_to_a_machine.py).
"""

import pytest

_OPEN = {"(": ")", "[": "]", "{": "}"}
_CLOSE = {v: k for k, v in _OPEN.items()}
# After these, a '/' starts a regex literal rather than a division.
_REGEX_AFTER = set("(,=:[!&|?{};+-*%<>~^") | {""}


def _structure_errors(js: str) -> list[str]:
    stack: list[tuple[str, int]] = []
    i, n = 0, len(js)
    last = ""  # last significant character outside literals
    while i < n:
        ch = js[i]
        if ch in " \t\r\n":
            i += 1
            continue
        if js.startswith("//", i):
            j = js.find("\n", i)
            i = n if j < 0 else j
            continue
        if js.startswith("/*", i):
            j = js.find("*/", i + 2)
            if j < 0:
                return [f"unterminated comment at {i}"]
            i = j + 2
            continue
        if ch in "'\"`":
            j = i + 1
            while j < n and js[j] != ch:
                if js[j] == "\\":
                    j += 1
                elif js[j] == "\n" and ch != "`":
                    return [f"unterminated string at {i}"]
                j += 1
            if j >= n:
                return [f"unterminated string at {i}"]
            i, last = j + 1, ch
            continue
        if ch == "/" and last in _REGEX_AFTER:
            j, in_class = i + 1, False
            while j < n:
                c = js[j]
                if c == "\\":
                    j += 2
                    continue
                if c == "\n":
                    return [f"unterminated regex at {i}"]
                if c == "[":
                    in_class = True
                elif c == "]":
                    in_class = False
                elif c == "/" and not in_class:
                    break
                j += 1
            i = j + 1
            while i < n and js[i].isalpha():  # flags
                i += 1
            last = "/"
            continue
        if ch in _OPEN:
            stack.append((ch, i))
        elif ch in _CLOSE:
            if not stack or stack[-1][0] != _CLOSE[ch]:
                return [f"unexpected {ch!r} at {i}: {js[max(0, i - 60):i + 1]!r}"]
            stack.pop()
        last = ch
        i += 1
    return [f"unclosed {c!r} at {p}: {js[p:p + 60]!r}" for c, p in stack]


def _check(js: str, name: str) -> None:
    errors = _structure_errors(js)
    assert not errors, f"{name} produced malformed JavaScript: {errors[:3]}"


def test_the_checker_catches_the_original_defect():
    """Control: the shape of the 2026-11-02 bug must be refused."""
    good = "(function(){ viewer.setStyle({serial: s}, {stick: {}}); })();"
    bad = "(function(){ viewer.setStyle({serial: s},\n var bt = 1; })();"
    assert _structure_errors(good) == []
    assert _structure_errors(bad)
    assert _structure_errors("var a = 'x;")
    assert _structure_errors("x.split(/\\s+/); var b = [1, 2;") != []
    assert _structure_errors("x.split(/\\s+/); var b = a / 2;") == []


def test_bootstrap_js_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_bootstrap_js
    _check(
        ts_goat_viewer_bootstrap_js(
            "vid-test", arrows=[0.1, 0.2, 0.3] * 3),
        "ts_goat_viewer_bootstrap_js")


def test_step_js_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_step_js
    _check(ts_goat_viewer_step_js("vid-test", 1.0), "ts_goat_viewer_step_js")


@pytest.mark.parametrize("running", [True, "pause", False])
def test_frames_js_is_valid(running):
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_frames_js
    _check(
        ts_goat_viewer_frames_js("vid-test", running=running),
        "ts_goat_viewer_frames_js")


def test_center_js_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_center_js
    _check(ts_goat_viewer_center_js("vid-test"), "ts_goat_viewer_center_js")


def test_bootstrap_without_arrows_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_bootstrap_js
    _check(
        ts_goat_viewer_bootstrap_js("vid-test", arrows=None),
        "ts_goat_viewer_bootstrap_js(no arrows)")


def test_highlight_js_is_valid():
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_highlights_js
    _check(
        ts_goat_viewer_highlights_js("vid-test", [0, 1, 8]),
        "ts_goat_viewer_highlights_js")
