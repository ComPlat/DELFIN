"""Red tests for package B3: delfin/agent/test_honesty.check().

Each finding kind gets a positive case (fires) and a negative case
(stays quiet), written BEFORE the module exists so they are red on the
unfixed code — see the task brief and the wave-17 lesson that tests
loaded .gate/ copies and passed only in the builder's tree. Every test
here exercises the PUBLIC path: ``test_honesty.check(diff_or_files)``
on real source text, never an inner helper.
"""

import textwrap

import delfin.agent.test_honesty as th


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Kind 1: bare `except:` / `except Exception: pass` (or `return None`)
#          around new code
# ---------------------------------------------------------------------------

NEW_CODE = (
    "delfin/agent/demo_module.py",
    40, 45,
    textwrap.dedent("""\
        def load(path):
            try:
                return open(path).read()
            except:
                pass
        """),
)

NARROW_EXCEPT = (
    "delfin/agent/demo_module.py",
    40, 48,
    textwrap.dedent("""\
        import logging
        _log = logging.getLogger(__name__)

        def load(path):
            try:
                return open(path).read()
            except OSError:
                _log.exception("cannot read %s", path)
                raise
        """),
)


def test_new_bare_except_is_flagged():
    findings = th.check([NEW_CODE])
    hits = [f for f in findings if f.path == NEW_CODE[0]]
    assert hits, "a new bare except: pass must be flagged"
    assert any("bare" in f.reason or "pass" in f.reason
               or "except" in f.reason for f in hits)
    assert all(f.line for f in hits), "every finding carries file:line"


def test_narrow_except_that_reraises_is_quiet():
    findings = th.check([NARROW_EXCEPT])
    assert [f for f in findings if f.path == NARROW_EXCEPT[0]] == []


def test_return_none_instead_of_result_is_flagged():
    case = (
        "delfin/agent/demo_module.py",
        50, 58,
        textwrap.dedent("""\
            def parse(text):
                try:
                    return json.loads(text)
                except Exception:
                    return None
            """),
    )
    findings = th.check([case])
    assert any(f.path == case[0] for f in findings)


def test_pre_existing_bare_except_outside_new_code_is_quiet():
    # The scanner judges only NEW code: a handler in the file's text
    # that lies OUTSIDE the record's hunk range must stay silent (the
    # file may be full of old bare excepts; we flag what THIS change
    # adds). Full-file text, hunk = lines 100-200: the handler at text
    # line 4 maps to file line 4, below the hunk start → quiet.
    full_file = (
        "delfin/agent/demo_module.py",
        100, 200,
        textwrap.dedent("""\
            def load(path):
                try:
                    return open(path).read()
                except:
                    pass

            x = 1
            """) + "y = 2\n" * 200,
    )
    findings = th.check([full_file])
    assert not [f for f in findings if f.path == full_file[0]], (
        "a handler outside the hunk range is pre-existing history")


# ---------------------------------------------------------------------------
# Kind 2: a committed tests/test_zz_* file
# ---------------------------------------------------------------------------

def test_committed_test_zz_file_is_flagged():
    case = ("tests/test_zz_probe_something.py", 1, 5,
            "def test_probe():\n    assert True\n")
    findings = th.check([case])
    assert any(f.path == case[0] for f in findings), (
        "a committed tests/test_zz_* file must be flagged")
    assert any("test_zz" in f.reason for f in findings)


def test_untracked_test_zz_file_is_quiet():
    # The locked API: tracked-ness is the CALLER's determination — the
    # caller hands over only committed files, so an untracked probe
    # never appears as a record with source. check() never touches the
    # disk or git to re-derive it (purity). A path-only mention
    # (source=None) is also never judged for kind 2.
    records = [
        ("tests/test_zz_b3_quiet.py", 1, 0, None),   # path-only mention
        ("tests/test_b3_present.py", 1, 2, "x = 1\n"),
    ]
    findings = th.check(records)
    assert not [f for f in findings if "test_zz" in f.path]


# ---------------------------------------------------------------------------
# Kind 3: a test that reads from .gate/
# ---------------------------------------------------------------------------

def test_test_reading_gate_dir_is_flagged():
    case = (
        "tests/test_b3_gate_reader.py",
        1, 10,
        textwrap.dedent("""\
            def test_something():
                data = open(".gate/px/new_module.py").read()
                assert data
            """),
    )
    findings = th.check([case])
    assert any(f.path == case[0] for f in findings), (
        "a test reading from .gate/ must be flagged")
    assert any(".gate" in f.reason for f in findings)


def test_test_reading_tmp_data_is_quiet():
    case = (
        "tests/test_b3_tmp_user.py",
        1, 10,
        textwrap.dedent("""\
            def test_something(tmp_path):
                data = open(tmp_path / "in.txt").read()
                assert data
            """),
    )
    findings = th.check([case])
    assert not [f for f in findings if f.path == case[0]]


# ---------------------------------------------------------------------------
# Kind 4: a settings switch introduced by the change with no test that
#          sets it True AND one that leaves it off
# ---------------------------------------------------------------------------

SWITCH_GOOD = [
    ("delfin/agent/feature.py", 1, 20,
     'FEATURE_SWITCH = True\n'),
    ("tests/test_b3_feature.py", 1, 30,
     "def test_on():\n"
     "    mod.FEATURE_SWITCH = True\n"
     "    assert run()\n\n"
     "def test_off():\n"
     "    mod.FEATURE_SWITCH = False\n"
     "    assert not run()\n"),
]

SWITCH_ON_ONLY = [
    ("delfin/agent/feature.py", 1, 20,
     'FEATURE_SWITCH = True\n'),
    ("tests/test_b3_feature.py", 1, 30,
     "def test_on():\n"
     "    mod.FEATURE_SWITCH = True\n"
     "    assert run()\n"),
]


def test_switch_with_both_states_tested_is_quiet():
    findings = th.check(SWITCH_GOOD)
    assert not [f for f in findings
                if "switch" in f.reason.lower()]


def test_switch_tested_only_on_is_flagged():
    findings = th.check(SWITCH_ON_ONLY)
    assert any("switch" in f.reason.lower() for f in findings), (
        "a settings switch tested only in one state must be flagged")


# ---------------------------------------------------------------------------
# Envelope: purity, envelope shape, and the real-diff path
# ---------------------------------------------------------------------------

def test_finding_envelope_carries_file_line_reason():
    findings = th.check([NEW_CODE])
    for f in findings:
        assert f.path
        assert f.line
        assert f.reason


def test_empty_input_returns_empty_list():
    assert th.check([]) == []


def test_real_change_record_with_source_scans_clean():
    """A record with REAL source text (like the wiring's caller hands
    over) that follows the rules stays quiet — replaces an earlier
    git-dependent version (reviewer finding on 20ebdab2: a clean tree
    made it fail, and source=None records measured nothing)."""
    clean_real = _record_full(
        "delfin/agent/demo_module.py",
        textwrap.dedent("""\
            import logging
            _log = logging.getLogger(__name__)


            def load(path):
                try:
                    return open(path).read()
                except OSError:
                    _log.exception("cannot read %s", path)
                    raise
            """),
        switch=("DEMO_ON", "tests/test_b3_demo.py",
                "mod.DEMO_ON = True\n\n\ndef test_off():\n"
                "    mod.DEMO_ON = False\n"
                "    assert not run()\n"),
    )
    assert th.check(clean_real) == []


def _record_full(path, text, switch=None):
    """A full-file record; optionally pairs a switch with a test file
    that exercises both states (kind-4 quiet shape)."""
    if switch:
        name, tpath, ttext = switch
        return [
            (path, 1, text.count("\n") + 1, text),
            (tpath, 1, ttext.count("\n") + 1, ttext
             + f"\ndef test_on():\n    mod.{name} = True\n    assert run()\n"),
        ]
    return [(path, 1, text.count("\n") + 1, text)]
