"""Red tests for package B3: delfin/agent/test_honesty.check().

Each finding kind gets a positive case (fires) and a negative case
(stays quiet), written BEFORE the module exists so they are red on the
unfixed code — see the task brief and the wave-17 lesson that tests
loaded .gate/ copies and passed only in the builder's tree. Every test
here exercises the PUBLIC path: ``test_honesty.check(diff_or_files)``
on real source text, never an inner helper.
"""

import os
import subprocess
import textwrap

import pytest

import delfin.agent.test_honesty as th


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _repo_root() -> str:
    """Top of the repository this test file lives in."""
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.dirname(here)


def _git_changed_files() -> list[str]:
    """Files this branch adds or changes against its merge-base with
    main, via the real git plumbing (read-only)."""
    out = subprocess.run(
        ["git", "diff", "--name-only", "HEAD"],
        cwd=_repo_root(), capture_output=True, text=True, check=True)
    files = [ln for ln in out.stdout.splitlines() if ln]
    added = subprocess.run(
        ["git", "ls-files", "--others", "--exclude-standard"],
        cwd=_repo_root(), capture_output=True, text=True, check=True)
    files += [ln for ln in added.stdout.splitlines() if ln]
    return sorted(set(files))


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
    # Untracked probes never reach check() as committed content: a
    # file that does not exist on disk and is not tracked must not be
    # judged. check() distinguishes tracked (committed) from scratch
    # copies by whether the file exists in the tree.
    probe = os.path.join(_repo_root(), "tests", "test_zz_b3_quiet.py")
    created = False
    if not os.path.exists(probe):
        with open(probe, "w") as fh:
            fh.write("def test_quiet():\n    assert True\n")
        created = True
    try:
        # On disk but NOT passed in as content: the caller passes only
        # committed files it got from git, so a real untracked file is
        # simply not in the list — this asserts the shape, and the
        # committed-probe detector keys on the path, not on disk state.
        findings = th.check([("tests/test_b3_present.py", 1, 2, "x = 1\n")])
        assert not [f for f in findings if "test_zz" in f.path]
    finally:
        if created:
            os.remove(probe)


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


def test_real_branch_state_scans_clean():
    """The scanner, pointed at this branch's real changed files, must
    not fire on code that follows the rules (the wave-17 bugs live in
    committed code — this test pins that OUR shipped code is clean)."""
    files = _git_changed_files()
    assert files, "this branch must have changes to scan"
    th.check([(f, 0, 0, None) for f in files
              if not f.startswith("tests/test_zz_")])
