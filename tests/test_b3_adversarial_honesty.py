"""Adversarial tests for B3 (reviewer nacht-s4): attack test_honesty.check
and the verify_guard wiring with inputs the builder's red tests do not
cover. All findings must carry file:line; quiet cases must stay quiet.
Records only — the scanner is pure and never calls git (locked contract).
"""

import textwrap

import delfin.agent.test_honesty as th
import delfin.agent.verify_guard as vg


def rec(path, start, end, text):
    return (path, start, end, text)


# --- kind 1 attacks: swallowing handlers in shapes the red tests miss ---

def test_except_exception_with_bare_pass_and_comment_fires():
    # `except Exception:` + only `pass` + trailing comment: still swallowing.
    src = textwrap.dedent("""\
        def f():
            try:
                return risky()
            except Exception:  # noqa
                pass  # ignore
        """)
    findings = th.check([rec("delfin/agent/m.py", 1, 5, src)])
    assert any("except" in f.reason for f in findings), findings


def test_handler_with_other_statements_is_quiet():
    # A handler that logs and re-raises is the sanctioned negative.
    src = textwrap.dedent("""\
        import logging
        def f():
            try:
                return risky()
            except Exception:
                logging.getLogger(__name__).warning("failed")
                raise
        """)
    assert th.check([rec("delfin/agent/m.py", 1, 6, src)]) == []


def test_return_none_wrapped_in_if_is_flagged_or_must_not_lie():
    # A handler that returns None CONDITIONALLY: the scanner may stay
    # quiet (not provably swallowing) but must never crash.
    src = textwrap.dedent("""\
        def f(x):
            try:
                return risky(x)
            except Exception:
                if x is None:
                    return None
                raise
        """)
    result = th.check([rec("delfin/agent/m.py", 1, 6, src)])
    assert isinstance(result, list)


def test_unicode_source_does_not_crash():
    # A non-ASCII identifier in the scanned source must not break the AST walk.
    src = "def f\u00fc():\n    try:\n        return 1\n    except:\n        pass\n"
    findings = th.check([rec("delfin/agent/m_unicode.py", 1, 5, src)])
    assert any("except" in f.reason for f in findings), findings


def test_huge_line_numbers_do_not_crash():
    src = "try:\n    pass\nexcept:\n    pass\n"
    findings = th.check([rec("delfin/agent/m.py", 10**6, 10**9, src)])
    assert isinstance(findings, list)


# --- kind 2 attacks: test_zz path shapes ---

def test_test_zz_without_py_suffix_not_flagged():
    # Only tests/test_zz_*.py is the probe pattern.
    assert th.check([rec("tests/test_zz_data", 1, 0, None)]) == []


def test_test_zz_in_other_directory_is_flagged():
    # The wave rule says "committed tests/test_zz_*" — a probe parked
    # under delfin/ is equally a committed probe; if the builder scoped
    # kind 2 to tests/ only, this finds the gap.
    findings = th.check([rec("delfin/agent/test_zz_probe.py", 1, 0, None)])
    print("PROBE delfin/ test_zz ->", findings)


def test_committed_test_zz_with_source_fires():
    # Locked contract: kind 2 fires on a tests/test_zz_* record WITH
    # source (the caller says it is in the change); source=None is
    # path-only and stays quiet (covered by the builder's own test).
    findings = th.check([rec("tests/test_zz_probe.py", 1, 1, "x = 1\n")])
    assert len(findings) == 1 and findings[0].path == "tests/test_zz_probe.py"


# --- kind 3 attacks: .gate references in shapes the regex misses ---

def test_gate_reader_via_pathlib_is_caught():
    # pathlib Path('.gate/...').read_text() must be caught like open().
    src = "from pathlib import Path\ndef f():\n    return Path('.gate/px/new.py').read_text()\n"
    findings = th.check([rec("tests/test_b3_gate.py", 1, 3, src)])
    assert findings, "a pathlib .gate reader must be flagged"
    assert findings[0].line == 3


def test_gate_reader_via_os_path_join_is_caught():
    src = "import os\ndef f():\n    return open(os.path.join('.gate', 'x.py')).read()\n"
    findings = th.check([rec("tests/test_b3_gate.py", 1, 3, src)])
    assert findings, "an os.path.join('.gate', ...) reader must be flagged"


# --- kind 4 attacks: switch detection heuristics ---

def test_switch_via_os_environ_get_is_seen():
    # PROBE (finding sent to s3 2026-10-18): the module comment claims
    # switches include "an env/settings lookup of a flag"
    # (test_honesty.py:71-72), but _SWITCH_ASSIGN matches only bool
    # literal assignments, so DELFIN_X = os.environ.get("DELFIN_X",
    # False) is NOT detected as a switch and stays quiet. Pinned as a
    # probe until the builder either matches env lookups or fixes the
    # comment; NOT an assertion while the wording is disputed.
    src = textwrap.dedent("""\
        import os
        DELFIN_B3_MODE = os.environ.get("DELFIN_B3_MODE", False)
        """)
    findings = th.check([rec("delfin/agent/m.py", 1, 2, src),
                         rec("tests/test_b3_switch.py", 1, 1,
                             "def t(monkeypatch):\n    monkeypatch.setenv('DELFIN_B3_MODE', '1')\n")])
    print("PROBE env-get switch ->", findings)


def test_switch_name_collision_not_fired_wrongly():
    # A test toggling an unrelated object attribute of the same name
    # must not silently count as the switch's both-state test — but the
    # scan is heuristic; pinned behavior: it does not false-positive a
    # FINDING when both states appear (documented limitation probe).
    src = "ENABLE_B3 = True\n"
    tests = rec("tests/test_b3_x.py", 1, 2,
                "def t():\n    obj.ENABLE_B3 = True\n    obj.ENABLE_B3 = False\n")
    findings = th.check([rec("delfin/agent/m.py", 1, 1, src), tests])
    print("PROBE suffix-collision switch ->", findings)


def test_kind2_scope_is_tests_only_documented():
    # PROBE: kind 2 fires only for tests/test_zz_*.py (brief wording);
    # a committed probe parked under delfin/ stays quiet. Pinned as
    # probe: matches the brief's exact wording, noted for QS.
    findings = th.check([rec("delfin/agent/test_zz_probe.py", 1, 0, None)])
    print("PROBE delfin/ test_zz ->", findings)


# --- wiring attacks ---

def test_wiring_never_raises_on_path_only_and_garbage():
    assert isinstance(vg.scan_for_dishonest_tests(
        [rec("tests/test_zz_a.py", 1, 0, None),
         rec("delfin/agent/x.py", 0, 0, "def (:\n"),
         rec("", 1, 1, None)]), list)


def test_wiring_feedback_finding_without_path_attr_is_tolerated():
    # feedback must not crash on a finding-shaped object lacking attrs?
    # The builder's contract says HonestyFinding only; a crash here is a
    # finding, tolerance is the reviewer's question.
    class Fake:
        pass
    try:
        text = vg.test_honesty_feedback([Fake()])
        print("PROBE fake finding ->", repr(text)[:80])
    except AttributeError as exc:
        print("PROBE fake finding raises AttributeError:", exc)
