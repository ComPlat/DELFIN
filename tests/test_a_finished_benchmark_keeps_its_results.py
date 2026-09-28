"""A completed benchmark run does not lose its scores to the filing.

Input: the results of one run. Output: a JSONL file, and the path it was
actually written to.

Measured 2026-09-28, on a login-node trial reported by the session that
built the bench suite:

    OSError: [Errno 30] Read-only file system:
    '.../.delfin/benchmark_runs/....jsonl'

Every task had run and been scored; the run died at the line that files
the scores, and the measurement was gone. A benchmark exists to produce
that file, so a directory that refuses it is stepped over rather than
raised on -- and the step is announced, because a reader looking where
the files have always been has to be told once where this one went.

An explicitly named ``runs_dir`` is NOT stepped over: a caller that named
a directory wants that directory, and writing somewhere else quietly
would hide the fault from whoever collects the files.
"""

from __future__ import annotations

import os
import pathlib

import pytest

from delfin.agent import benchmark as B


def _results():
    return []


def _readonly(tmp_path, name="ro"):
    d = tmp_path / name
    d.mkdir()
    os.chmod(d, 0o500)
    return d


@pytest.fixture(autouse=True)
def _no_scratch(monkeypatch):
    monkeypatch.delenv("DELFIN_SCRATCH_STATE", raising=False)


def test_the_normal_case_is_unchanged(tmp_path, monkeypatch):
    monkeypatch.setattr(B, "_DEFAULT_RUNS_DIR", tmp_path / "runs")
    out = B.write_run(_results(), model="m")
    assert out.parent == tmp_path / "runs"
    assert out.exists()


def test_a_refused_default_steps_to_the_scratch_state(tmp_path, monkeypatch,
                                                      capsys):
    monkeypatch.setattr(B, "_DEFAULT_RUNS_DIR",
                        _readonly(tmp_path) / "benchmark_runs")
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("DELFIN_SCRATCH_STATE", str(scratch))
    out = B.write_run(_results(), model="m")
    assert out.exists()
    assert scratch in out.parents, out
    said = capsys.readouterr().err
    assert "results written to" in said and "could not use" in said, (
        "the step was taken silently; a reader would look in the old place")


def test_with_no_scratch_it_still_lands_somewhere(tmp_path, monkeypatch):
    monkeypatch.setattr(B, "_DEFAULT_RUNS_DIR",
                        _readonly(tmp_path) / "benchmark_runs")
    out = B.write_run(_results(), model="m")
    assert out.exists(), "a completed run still lost its results"


def test_an_explicit_directory_is_not_stepped_over(tmp_path, monkeypatch):
    """A caller that named a directory gets that directory or an error."""
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("DELFIN_SCRATCH_STATE", str(scratch))
    with pytest.raises(OSError):
        B.write_run(_results(), model="m",
                    runs_dir=_readonly(tmp_path) / "named")


def test_when_nothing_is_writable_the_reasons_are_named(tmp_path,
                                                        monkeypatch):
    monkeypatch.setattr(B, "_DEFAULT_RUNS_DIR",
                        _readonly(tmp_path, "a") / "runs")
    monkeypatch.setenv("DELFIN_SCRATCH_STATE", str(_readonly(tmp_path, "b")))
    monkeypatch.setattr(B.tempfile, "gettempdir",
                        lambda: str(_readonly(tmp_path, "c")))
    with pytest.raises(OSError) as caught:
        B.write_run(_results(), model="m")
    assert "tried" in str(caught.value), (
        "the failure does not say which places refused it")


def test_the_rows_are_written_once_not_per_attempt(tmp_path, monkeypatch):
    """The results are serialised before the first attempt: a generator
    consumed by a failed write would leave the retry with nothing."""
    import inspect

    src = inspect.getsource(B.write_run)
    at_rows = src.index("rows = [")
    at_loop = src.index("for d in candidates")
    assert at_rows < at_loop, (
        "results are serialised inside the retry loop; an Iterable would "
        "be empty by the second attempt")
