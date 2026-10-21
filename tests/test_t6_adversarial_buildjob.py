"""T6 adversarial review of ``slurm_tests.build_job`` (phase 2, nacht-s24).

build_job renders a Bash script that runs on a SLURM compute node. Everything
the caller supplies crosses into shell text, so every interpolated value must
stay inside a shell context that cannot be re-parsed as a command. The tests
here push on the ALLOW surface the builder's own suite did not cover:

* command substitution smuggled through ``ref`` / ``repo`` (double-quoting of
  an already-shlex-quoted value turns a payload into an executed ``$(...)``);
* the ``--export=ALL`` incident rule, exact string and no attached ``,VAR``;
* the node-local invariant (pytest runs under the staged python in $TMPDIR);
* the test-path ALLOW-list (unicode look-alikes, backslash, bare ``tests/``,
  node ids and ``::`` selectors that MUST keep working).
"""
from __future__ import annotations

import json
import os
import re
import sys
import tempfile
from pathlib import Path

import pytest

from delfin.agent.slurm_tests import _SAFE_TEST_PATH, build_job


def test_ref_command_substitution_stays_out_of_double_quotes():
    """A ref carrying ``$(...)`` must never reach the shell: build_job refuses it.

    Old decision (F1): the payload was single-quote-contained in the rendered
    argv. F4 hardened the contract: command substitution / metacharacter
    payloads in ref are now refused outright (ValueError), which enforces the
    same "never executed" intent more strongly.
    """
    for payload in (
        "$(touch /tmp/x)",
        "$(id)",
        "`touch /tmp/x`",
        "main$(echo hi)",
    ):
        with pytest.raises(ValueError):
            build_job("delfin", payload, ["tests/a.py"], "cpu", 10)


def test_repo_metacharacters_are_single_quoted_for_shell_commands():
    """repo must never let a metachar payload reach the archival line: refused.

    Old decision (F1): repo="repo$(id)" was single-quoted as -C 'repo$(id)'.
    F4 hardened the contract: shell metacharacters in repo are now refused
    (ValueError) before any rendering, which is the stronger form of the same
    "no extra command on the archival line" intent.
    """
    with pytest.raises(ValueError):
        build_job("repo$(id)", "main", ["tests/a.py"], "cpu", 10)


def test_export_all_exactly_once_and_no_attached_vars():
    """Incident rule: the script must carry --export=ALL and never a
    --export=ALL,VAR variant or a repeat."""
    text = build_job("delfin", "main", ["tests/a.py"], "cpu", 10)
    occurrences = re.findall(r"--export=[^\s]+", text)
    assert occurrences == ["--export=ALL"], occurrences


def test_no_export_all_vars_variant_present():
    text = build_job("delfin", "main", ["tests/a.py"], "cpu", 10)
    # --export=ALL,VARS or --export=ALL <space> VAR (an env leak) is banned.
    assert not re.search(r"--export=ALL,[A-Za-z0-9_]+", text)
    assert re.search(r"--export=ALL\b", text) is not None


def test_pytest_runs_under_the_staged_node_local_python():
    """The pytest invocation must use the node-local staged interpreter
    ($STAGE_PYTHON), never a HOME / shared-file-system python path."""
    text = build_job("delfin", "main", ["tests/x.py"], "cpu", 10)
    # pytest is invoked via "$STAGE_PYTHON"; no absolute HOME/shared python.
    m = re.search(r'"(?:\$STAGE_PYTHON|PYTHON_PLACEHOLDER)"\s+-m pytest', text)
    assert m is not None, "pytest must run through the staged interpreter"
    assert re.search(r'["\']/(?:home|users|root)?[^\s"\']*bin/python"\s+-m pytest', text) is None


def test_tests_run_from_node_local_copy_not_home():
    """The whole pytest block must execute with cwd in the TMPDIR copy."""
    text = build_job("delfin", "main", ["tests/x.py"], "cpu", 10)
    body = text[text.index('cd "$LOCAL"'):]
    # after cd $LOCAL, the only cd is into $LOCAL; no HOME cd.
    assert re.search(r"cd\s+(\$HOME|\$HOME/|~)", body) is None
    # test paths are relative to the node-local root, not absolute/HOME.
    assert "TESTS=(tests/x.py)" in body


def test_path_allowlist_rejects_unicode_look_alikes():
    """A fullwidth slash or look-alike must not pass the allow-list: the
    rendering is consumed by a real shell later."""
    for hostile in ("tests/é.py", "tests/／.py", "tests/a\u00b7py", "tests/a b.py"):
        assert not _SAFE_TEST_PATH.match(hostile), hostile
        with pytest.raises(ValueError):
            build_job("delfin", "main", [hostile], "cpu", 10)


def test_path_allowlist_rejects_backslash_and_bare_tests_dir():
    """Backslash separators and a bare ``tests/`` (scope the whole suite)
    are refused; a single test file under tests/ passes."""
    for hostile in ("tests\\a.py", "tests/"):
        assert not _SAFE_TEST_PATH.match(hostile), hostile
        with pytest.raises(ValueError):
            build_job("delfin", "main", [hostile], "cpu", 10)
    assert _SAFE_TEST_PATH.match("tests/a.py")
    assert build_job("delfin", "main", ["tests/a.py"], "cpu", 10) is not None


def test_node_id_selectors_and_colons_still_accepted():
    """Legitimate pytest selectors (file::Class::method, [param]) must survive
    the allow-list — they are the documented <paths> input."""
    good = "tests/test_x.py::TestY::test_z"
    assert _SAFE_TEST_PATH.match(good)
    text = build_job("delfin", "main", [good], "cpu", 10)
    assert good in text


# --- node-side JSON summary accuracy (run the REAL embedded logic) ---

def _extract_summary_py(script: str) -> str:
    """The python source of the job's <<PY heredoc, verbatim."""
    start = script.index("<<'PY'") + len("<<'PY'\n")
    end = script.index("\nPY\n", start)
    return script[start:end]


def test_summary_counts_match_the_pytest_log():
    """The JSON the job emits must report the real pass/fail numbers from the
    pytest log, not a count of how often the words 'passed'/'failed' appear."""
    script = build_job("delfin", "main", ["tests/x.py"], "cpu", 10)
    inner = _extract_summary_py(script)
    # the SHELL invokes the summary writer as: python - <rc> <ref> <out>
    # (slurm_tests.py:151), so argv[1]=rc, argv[2]=ref, argv[3]=out. Reproduce
    # that exact order; the code must bind rc/out/ref accordingly and write a
    # valid summary to <out>.
    for log_text, exp_passed, exp_failed in [
        ("........................\n25 passed in 2.50s\n", 25, 0),
        ("...............F..........\n1 failed, 24 passed in 2.60s\n", 24, 1),
        ("F.F.F..F.F.F\n8 failed, 4 passed in 1.10s\n", 4, 8),
    ]:
        ns: dict = {}
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            logp = td / "pytest.log"
            logp.write_text(log_text)
            outp = td / "summary.json"
            sys.argv = ["py", "1", "main", str(outp)]
            os.environ["DELFIN_TJ_LOG"] = str(logp)
            try:
                exec(compile(inner, "<summar>", "exec"), ns)
            finally:
                os.environ.pop("DELFIN_TJ_LOG", None)
            d = json.loads((td / "summary.json").read_text())
        assert (d["passed_tests"], d["failed_tests"]) == (exp_passed, exp_failed), (
            f"log {log_text.strip()!r} -> summary {d['passed_tests']}/{d['failed_tests']}, "
            f"expected {exp_passed}/{exp_failed}"
        )
        assert d["ref"] == "main"

