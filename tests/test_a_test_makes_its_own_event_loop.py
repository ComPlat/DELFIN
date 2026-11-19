"""No test drives a coroutine through ``asyncio.get_event_loop()``.

``asyncio.run`` closes the loop it made AND unsets it. So the first test
in a run that uses it leaves the thread with no current loop, and every
later ``get_event_loop()`` raises "There is no current event loop in
thread 'MainThread'". The test that does it passes alone and fails in the
run that matters.

This has now happened twice, which is why it is a scanner and not a
comment. The first time,
``tests/test_the_server_ends_what_nobody_watches.py`` wrote the lesson
down in its own ``_run`` helper: "Ten of these passed alone and failed
in the run that matters." The second time, 2026-11-18, a new test in
``test_a_dashboard_login_serves_named_accounts.py`` reintroduced it and
reached main -- it was reviewed by running its own file, where it passes,
and CI then failed shard 4 of a later pull request with the same
RuntimeError. Written-down knowledge that nothing enforces does not
bind.

``asyncio.run(coro)`` is the replacement, and it needs no fixture.
"""

from __future__ import annotations

import re
from pathlib import Path

TESTS = Path(__file__).resolve().parent

#: The call, not the words. A line that merely mentions the name -- this
#: file, and the comment in the file that first recorded the lesson -- is
#: not a call, and a scanner that tripped on its own explanation would be
#: the defect it is meant to prevent. (Learned on a sibling guard the day
#: before, which refused the comment that documented it.)
_CALL = re.compile(r"\basyncio\s*\.\s*get_event_loop\s*\(")


def _offending_lines(path: Path) -> list[tuple[int, str]]:
    out: list[tuple[int, str]] = []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return out
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if line.startswith("#"):
            continue                      # an explanation is not a call
        code = line.split("#", 1)[0]
        if _CALL.search(code):
            out.append((number, line[:100]))
    return out


def test_no_test_uses_get_event_loop():
    found: list[str] = []
    for path in sorted(TESTS.rglob("*.py")):
        if path.name == Path(__file__).name:
            continue
        for number, line in _offending_lines(path):
            found.append(f"{path.relative_to(TESTS)}:{number}: {line}")
    assert not found, (
        "asyncio.get_event_loop() in a test raises once any earlier test "
        "has called asyncio.run, which closes and unsets the loop. Use "
        "asyncio.run(coro). Found:\n" + "\n".join(found))


def test_the_scanner_sees_a_real_call(tmp_path):
    """A guard that cannot fail is not a guard."""
    bad = tmp_path / "test_bad.py"
    bad.write_text("import asyncio\n"
                   "def test_x():\n"
                   "    asyncio.get_event_loop().run_until_complete(None)\n")
    assert [n for n, _ in _offending_lines(bad)] == [3]


def test_the_scanner_ignores_an_explanation(tmp_path):
    """The sibling guard of 2026-11-17 refused the comment that documented
    it. This one is pinned not to."""
    fine = tmp_path / "test_fine.py"
    fine.write_text("import asyncio\n"
                    "# asyncio.run rather than get_event_loop(): the loop\n"
                    "# is unset by then.\n"
                    "def test_x():\n"
                    "    return asyncio.run(None)  # not get_event_loop()\n")
    assert _offending_lines(fine) == []
