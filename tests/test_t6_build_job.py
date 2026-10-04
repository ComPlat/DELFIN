"""T6 phase 2: ``slurm_tests.build_job`` renders a node-local test job.

The wave grows only when sessions' tests stop running on the shared login
node. ``build_job`` turns a git ref plus a set of test paths into the text of
a SLURM job script that re-runs exactly those tests on a compute node from a
node-local copy:

* the repo is copied at the given ref onto node-local disk (``$TMPDIR``) via
  ``git archive`` -- never run from HOME / the shared tree;
* the interpreter's venv is staged node-local too;
* only the named test paths run, not a whole suite;
* only the JSON summary and the job's own log are written back;
* ``--export=ALL`` exactly -- the incident rule that ``--export=ALL,VARS``
  leaks HOME state into compute jobs.

``build_job`` is pure: same arguments, identical text, no I/O. The tests
inspect the returned script text; nothing reads the network or the
scheduler.
"""

from __future__ import annotations

import re

import pytest

from delfin.agent.slurm_tests import build_job


def test_returns_a_bash_script_text():
    text = build_job("delfin", "main", ["tests/test_x.py"], "cpu", 10)
    assert text.lstrip().startswith("#!/bin/bash")
    assert "#SBATCH" in text


def test_export_is_always_exactly_all():
    """The incident rule: --export=ALL, never --export=ALL,VARS or a bare
    environment-variable ambiguity that pulls login-node state along."""
    text = build_job("delfin", "main", ["tests/test_x.py"], "cpu", 10)
    for m in re.finditer(r"--export=[^\s]+", text):
        assert m.group(0) == "--export=ALL", m.group(0)


def test_copies_the_ref_node_local_via_git_archive():
    text = build_job("delfin", "myref123", ["tests/test_x.py"], "cpu", 10)
    # The script must put a copy of the repo at the given ref onto node-local
    # disk (TMPDIR), not touch HOME or the shared checkout.
    assert "git archive" in text
    assert "myref123" in text
    assert "$TMPDIR" in text or "${TMPDIR" in text
    assert re.search(r"cd\s+(\$HOME\b|~/|\$HOME/)", text, re.MULTILINE) is None


def test_runs_only_the_named_test_paths():
    text = build_job("delfin", "main", ["tests/a.py", "tests/b.py::C::t"], "cpu", 10)
    assert "tests/a.py" in text
    assert "tests/b.py::C::t" in text


def test_runs_pytest_from_the_node_local_copy_not_home():
    """No pytest may execute against the shared checkout / HOME paths."""
    text = build_job("delfin", "main", ["tests/x.py"], "cpu", 10)
    # pytest is invoked on the staged copy; a pytest of a HOME path is absent.
    assert re.search(r"pytest(\s[^#]*)?\s/(root|home|users)?/?~?", text) is None
    assert "~/tests" not in text


def test_stages_the_venv_node_local():
    text = build_job("delfin", "main", ["tests/x.py"], "cpu", 10)
    # The interpreter is unpacked onto node-local disk, not used from the
    # shared file system path in its original location. The node-local root
    # derives from TMPDIR, and the staged interpreter is found under it.
    assert re.search(r'LOCAL="\$\{TMPDIR:-/tmp\}/delfin_tests_', text) is not None
    assert re.search(r'STAGE_PYTHON="\$LOCAL/.*/bin/python"', text) is not None
    assert 'STAGE_PYTHON="$LOCAL/$(basename "$VENV_DIR")/bin/python"' in text
    # The venv is unpacked into $LOCAL, and a working staged interpreter is
    # preferred over the shared file-system python.
    assert 'tar -x -C "$LOCAL"' in text


def test_writes_a_json_summary():
    text = build_job("delfin", "main", ["tests/x.py"], "cpu", 10)
    assert ".json" in text
    assert "summary" in text.lower() or "result" in text.lower()


def test_partition_and_minutes_appear_in_the_header():
    text = build_job("delfin", "main", ["tests/x.py"], "gpu_4", 25)
    assert "--partition=gpu_4" in text
    assert re.search(r"--time=\S*25", text) is not None


def test_deterministic_pure():
    a = build_job("delfin", "main", ["tests/x.py"], "cpu", 10)
    b = build_job("delfin", "main", ["tests/x.py"], "cpu", 10)
    assert a == b


def test_rejects_bad_test_paths():
    """Test paths that are not under tests/ are refused -- a malicious or
    mistaken path must not become a pytest target outside the tree."""
    with pytest.raises(ValueError):
        build_job("delfin", "main", ["/etc/passwd"], "cpu", 10)
    with pytest.raises(ValueError):
        build_job("delfin", "main", ["../outside.py"], "cpu", 10)
