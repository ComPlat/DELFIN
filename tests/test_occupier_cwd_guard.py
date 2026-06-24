"""CWD-race guard for the parallel OCCUPIER execution.

Regression for the MarcDy3 bug. OCCUPIER FoB jobs run in a single-process
``ThreadPoolExecutor``, where ``os.chdir`` mutates the *process*-wide working
directory. Several independent CWD locks let one thread ``chdir`` while another
performed CWD-relative ``CONTROL.txt`` writes — leaking an OCCUPIER sub-folder's
CONTROL (``input_file=input.xyz`` + the auto-sequence block) into the main run
directory. A later recalc then read the leaked ``input.xyz`` reference, could not
find the (never-staged) file, and aborted at the electron-count step.

The fix routes every ``os.chdir`` through a single process-wide guard
(``delfin.common.paths.pushd`` / ``GLOBAL_CWD_LOCK``) and anchors
``run_OCCUPIER``'s CONTROL read/write to an explicit ``work_dir``.
"""

import os
import threading

from delfin.common.paths import GLOBAL_CWD_LOCK, pushd


def test_all_occupier_cwd_locks_are_unified():
    """Every os.chdir site must serialize on the *same* lock object; otherwise
    two locks provide no mutual exclusion and the CWD can still be raced."""
    import delfin.occupier_flat_extraction as fe
    import delfin.workflows.engine.occupier as eo

    assert fe._cwd_lock is GLOBAL_CWD_LOCK
    assert eo._cwd_lock is GLOBAL_CWD_LOCK


def test_run_occupier_accepts_work_dir():
    """run_OCCUPIER must expose an optional, backward-compatible work_dir."""
    import inspect

    from delfin.occupier import run_OCCUPIER

    params = inspect.signature(run_OCCUPIER).parameters
    assert "work_dir" in params
    assert params["work_dir"].default is None  # default == historical behaviour


def test_pushd_restores_cwd_and_is_reentrant(tmp_path):
    start = os.getcwd()
    sub = tmp_path / "a"
    sub.mkdir()
    try:
        with pushd(tmp_path):
            assert os.path.realpath(os.getcwd()) == os.path.realpath(str(tmp_path))
            with pushd(sub):  # re-entrant (RLock) — nested chdir must work
                assert os.path.realpath(os.getcwd()) == os.path.realpath(str(sub))
            assert os.path.realpath(os.getcwd()) == os.path.realpath(str(tmp_path))
        assert os.getcwd() == start
    finally:
        os.chdir(start)


def test_pushd_serializes_concurrent_chdir(tmp_path):
    """Two threads each ``pushd`` into their own directory and do a CWD-relative
    write — the essence of the leaked CONTROL write. Under the global guard the
    write+read always stays consistent with the thread's own directory, even
    though ``os.chdir`` is process-global. Without the guard this interleaves and
    the file content/location leaks across threads.
    """
    dir_a = tmp_path / "A"
    dir_a.mkdir()
    dir_b = tmp_path / "B"
    dir_b.mkdir()
    errors: list = []
    barrier = threading.Barrier(2)

    def worker(target):
        try:
            barrier.wait()
            for _ in range(200):
                with pushd(target):
                    with open("marker.txt", "w") as fh:
                        fh.write(target.name)
                    with open("marker.txt") as fh:
                        assert fh.read() == target.name
        except Exception as exc:  # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=worker, args=(d,)) for d in (dir_a, dir_b)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors, errors
    assert (dir_a / "marker.txt").read_text() == "A"
    assert (dir_b / "marker.txt").read_text() == "B"
