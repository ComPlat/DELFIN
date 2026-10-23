"""run_tests with several files, with a timeout, and with a target it cannot find.

Driven 2026-09-16 in a dashboard session: run_tests(target="tests/a.py
tests/b.py ...") reported "failed" while the same files passed from bash --
the string reached pytest as ONE missing path (exit 4). A second call timed
out and the tool crashed with "sequence item 0: expected str instance, bytes
found", hiding the timeout.
"""
import shutil
import tempfile
import textwrap
from pathlib import Path

from delfin.agent import test_runner as TR


def _repo():
    d = Path(tempfile.mkdtemp(prefix="rt_"))
    (d / "test_a.py").write_text("def test_a():\n    assert True\n")
    (d / "test_b.py").write_text("def test_b():\n    assert True\n")
    return d


def test_several_space_separated_files_are_several_arguments():
    repo = _repo()
    try:
        r = TR.run_tests(repo, target="test_a.py test_b.py")
        assert r["status"] == "ok", r
        assert r["summary"]["passed"] == 2
    finally:
        shutil.rmtree(repo, ignore_errors=True)


def test_a_node_id_with_spaces_stays_one_argument(tmp_path):
    (tmp_path / "test_p.py").write_text("def test_x():\n    pass\n")
    assert TR.split_target("test_p.py::test_x[a b]", tmp_path) == ["test_p.py::test_x[a b]"]
    assert TR.split_target("test_p.py missing.py", tmp_path) == ["test_p.py missing.py"]


def test_a_target_pytest_cannot_find_is_an_error_not_a_failure():
    repo = _repo()
    try:
        r = TR.run_tests(repo, target="no_such_test.py")
        assert r["status"] == "error", r
    finally:
        shutil.rmtree(repo, ignore_errors=True)


def test_a_timeout_with_long_output_is_reported_not_crashed():
    repo = Path(tempfile.mkdtemp(prefix="rt_"))
    try:
        (repo / "test_slow.py").write_text(textwrap.dedent("""\
            import sys, time
            def test_slow():
                for i in range(80):
                    print("line", i, flush=True)
                    sys.stdout.flush()
                time.sleep(30)
            """))
        r = TR.run_tests(repo, target="test_slow.py", pytest_args=["-s"], timeout_s=5)
        assert r["status"] == "timeout", r
        assert isinstance(r["raw_stdout_tail"], str)
    finally:
        shutil.rmtree(repo, ignore_errors=True)


def test_tail_takes_bytes():
    assert TR._tail(b"a\nb\n", n=1).endswith("b")
