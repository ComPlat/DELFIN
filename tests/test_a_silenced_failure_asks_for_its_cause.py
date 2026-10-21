"""A failure silenced rather than removed asks for its cause, and a fix
answer names its root cause.

Input: the object-side findings of ``edit_signals`` (a change's added
lines, beyond the removed ones) and its test runs. Output under test:
``edit_signals.turn_caveat``.

Pinned here:

* in a fix turn, an added ``except Exception: pass`` / bare ``except:``
  that does nothing, an empty ``catch {}``, ``|| true`` and ``set +e``,
  ``continue-on-error``, an added sleep or retry loop, a raised timeout,
  a silenced warning, a rerun plugin or flaky marker, and a test re-run
  unchanged until it passed are each named with the file;
* no alarm when the code already had the construct (only the excess of
  added over removed lines counts), for new code written with its own
  handler (only a handler around a statement the old text had counts),
  when the user asked for a workaround,
  when the answer names the cause in the sentence that names the file,
  and outside a fix turn;
* a fix turn that changed the object and whose answer names no cause
  gets the root-cause line; an answer stating the cause, or that it was
  not found, gets none.
"""

from __future__ import annotations

import pytest

from delfin.agent import edit_signals as E

S = "sess-symptom"

OLD = "def load(path):\n    data = read(path)\n    return parse(data)\n"


@pytest.fixture(autouse=True)
def _fresh():
    E._TURNS.clear()
    yield
    E._TURNS.clear()


def _red():
    E.note_command(S, "pytest -q", "1 failed, 9 passed")


@pytest.mark.parametrize("path,old,new,kind", [
    ("src/io.py", OLD, "def load(path):\n    try:\n        data = read(path)\n"
     "    except Exception:\n        pass\n    return parse(data)\n",
     "failure swallowed"),
    ("src/io.py", OLD, "def load(path):\n    try:\n        data = read(path)"
     "\n    except: return None\n    return parse(data)\n",
     "failure swallowed"),
    ("web/app.js", "load();\n", "try { load(); } catch (e) {}\n",
     "failure swallowed"),
    ("scripts/run.sh", "make check\n", "make check || true\n",
     "failure swallowed"),
    ("scripts/run.sh", "#!/bin/sh\nmake\n", "#!/bin/sh\nset +e\nmake\n",
     "failure swallowed"),
    ("src/io.py", OLD, OLD.replace("    data = read(path)\n",
                                   "    time.sleep(2)\n    data = read(path)\n"),
     "sleep or retry added"),
    ("src/io.py", OLD, OLD.replace("    data = read(path)\n",
                                   "    for attempt in range(5):\n"
                                   "        data = read(path)\n"),
     "sleep or retry added"),
    ("src/net.py", "TIMEOUT = 0\nr = get(url, timeout=10)\n",
     "TIMEOUT = 0\nr = get(url, timeout=120)\n", "timeout raised"),
    ("src/io.py", OLD, "import warnings\nwarnings.filterwarnings('ignore')\n"
     + OLD, "warning or error silenced"),
    ("pytest.ini", "[pytest]\naddopts = -q\n",
     "[pytest]\naddopts = -q --reruns 3\n", "rerun of failing tests added"),
])
def test_a_symptom_fix_in_a_fix_turn_is_named(path, old, new, kind):
    _red()
    E.note_change(S, path, old, new)
    cav = E.turn_caveat(S, "the test fails, fix it", "Fixed; it passes now.")
    assert "A failure is silenced rather than removed" in cav, cav
    assert kind in cav and path in cav
    assert "Name the cause it hides" in cav


def test_a_test_rerun_unchanged_until_green_is_named():
    E.note_command(S, "pytest tests/test_a.py", "1 failed")
    E.note_command(S, "pytest tests/test_a.py", "1 failed")
    E.note_command(S, "pytest tests/test_a.py", "3 passed")
    cav = E.turn_caveat(S, "fix the flaky test", "It passes now.")
    assert "flaky test re-run until it passed" in cav


def test_a_rerun_after_a_change_is_no_rerun_until_green():
    E.note_command(S, "pytest tests/test_a.py", "1 failed")
    E.note_change(S, "src/a.py", "x = 1\n", "x = 2\n")
    E.note_command(S, "pytest tests/test_a.py", "3 passed")
    assert E.symptom_fixes(S) == []


def test_code_that_already_had_it_is_not_flagged():
    _red()
    old = ("def f():\n    try:\n        g()\n    except Exception:\n"
           "        pass\n    time.sleep(1)\n")
    new = old.replace("g()", "h()")
    E.note_change(S, "src/f.py", old, new)
    assert E.symptom_fixes(S) == []


def test_new_code_with_its_own_best_effort_handler_is_not_flagged():
    # A new helper written with a handler of its own is ordinary code; the
    # symptom fix is a handler put around a statement that was there.
    _red()
    E.note_change(S, "src/io.py", OLD, OLD + "\n\ndef note(x):\n    try:\n"
                  "        LOG.write(x)\n    except Exception:\n        pass\n")
    assert E.symptom_fixes(S) == []


def test_a_narrow_handler_that_handles_is_not_flagged():
    _red()
    E.note_change(S, "src/io.py", OLD, OLD.replace(
        "    data = read(path)\n",
        "    try:\n        data = read(path)\n    except FileNotFoundError:\n"
        "        raise ConfigMissing(path)\n"))
    assert E.symptom_fixes(S) == []


@pytest.mark.parametrize("request_text", [
    "Add a retry around the download for now, the server is flaky.",
    "Just a quick fix please, raise the timeout.",
    "Bitte eine vorübergehende Notlösung, wir reparieren es später.",
    "Suppress the deprecation warning in the CLI.",
])
def test_no_alarm_when_the_user_asked_for_a_workaround(request_text):
    _red()
    E.note_change(S, "src/io.py", OLD, OLD.replace(
        "    data = read(path)\n", "    time.sleep(2)\n    data = read(path)\n"))
    assert "silenced" not in E.turn_caveat(S, request_text, "Done.")


def test_no_alarm_when_the_answer_names_the_cause_beside_the_file():
    _red()
    E.note_change(S, "src/io.py", OLD, OLD.replace(
        "    data = read(path)\n", "    time.sleep(2)\n    data = read(path)\n"))
    answer = ("The sleep in io.py is a workaround because the root cause is "
              "the NFS attribute cache, which the server team fixes.")
    assert E.turn_caveat(S, "fix the failing load test", answer) == ""


def test_no_alarm_outside_a_fix_turn():
    E.note_change(S, "src/io.py", OLD, OLD.replace(
        "    data = read(path)\n", "    time.sleep(2)\n    data = read(path)\n"))
    assert E.turn_caveat(S, "Add a polling loop to the watcher.", "Done.") == ""


def test_a_fix_answer_without_a_cause_gets_the_root_cause_line():
    E.note_change(S, "src/calc.py", "return x * x + 1\n", "return x * x\n")
    cav = E.turn_caveat(S, "Fix the bug in area().", "Fixed area(). Tests pass.")
    assert "names no root cause" in cav
    assert "state it in one sentence, or say it was not found" in cav


@pytest.mark.parametrize("answer", [
    "Fixed. The bug was an extra +1 in area().",
    "Root cause: area() added 1 to every result.",
    "area() returned x*x+1 because of a copy-paste from perimeter().",
    "The root cause was not found; the change works around it.",
    "Behoben: Die Ursache war ein zusätzliches +1 in area().",
    "It failed because the index was off by one.",
])
def test_a_fix_answer_that_names_its_cause_gets_nothing(answer):
    E.note_change(S, "src/calc.py", "return x * x + 1\n", "return x * x\n")
    assert E.turn_caveat(S, "Fix the bug in area().", answer) == ""


def test_no_root_cause_line_when_only_tests_or_docs_changed():
    E.note_change(S, "tests/test_calc.py", None, "def test_a():\n    pass\n")
    assert E.turn_caveat(S, "Fix the docs typo.", "Done.") == ""
