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


# ---------------------------------------------------------------------------
# Defense-in-depth: the run-root CONTROL.txt must never receive auto-sequence
# rewrites, even if an unanchored/raced write ever targets it again. The leak
# above is fixed at the cause (work_dir anchoring); this is the structural
# backstop so the master CONTROL can never be silently corrupted.
# ---------------------------------------------------------------------------

_ROOT_CONTROL = (
    "input_file=input.xyz\n"
    "charge=+2\n"
    "OCCUPIER_method=auto\n"
    "OCCUPIER_tree=own\n"
    "OCCUPIER_sequence_profiles:\n"
    "-3,-2,-1,0,+1,+2,+3=[\n"
    "even electron number:\n"
    "even_seq = [\n"
    '  {"index": 1, "m": 1, "BS": "", "from": 0}\n'
    "]\n"
    "]\n"
    "INFOS:\n"
    "some trailing template\n"
)

_BUNDLE = {"even_seq": [{"index": 1, "m": 5, "BS": "", "from": 0}]}


def test_append_sequence_overrides_refuses_run_root(tmp_path, monkeypatch):
    """A run root (a dir owning *_OCCUPIER sub-folders) must not get the auto block."""
    monkeypatch.delenv("DELFIN_OCC_ROOT", raising=False)
    from delfin.occupier_sequences import append_sequence_overrides

    (tmp_path / "initial_OCCUPIER").mkdir()  # makes tmp_path a run root
    control = tmp_path / "CONTROL.txt"
    control.write_text(_ROOT_CONTROL, encoding="utf-8")

    append_sequence_overrides(control, _BUNDLE)

    assert "# AUTO sequence overrides" not in control.read_text(encoding="utf-8")
    assert "OCCUPIER_sequence_profiles:" in control.read_text(encoding="utf-8")


def test_remove_existing_sequence_blocks_refuses_run_root(tmp_path, monkeypatch):
    """force=True strip must be refused on the run-root CONTROL (would discard the
    user's OCCUPIER_sequence_profiles block)."""
    monkeypatch.delenv("DELFIN_OCC_ROOT", raising=False)
    from delfin.occupier_sequences import remove_existing_sequence_blocks

    (tmp_path / "ox_step_1_OCCUPIER").mkdir()  # makes tmp_path a run root
    control = tmp_path / "CONTROL.txt"
    control.write_text(_ROOT_CONTROL, encoding="utf-8")

    result = remove_existing_sequence_blocks(control, force=True)

    assert control.read_text(encoding="utf-8") == _ROOT_CONTROL  # unchanged on disk
    assert "OCCUPIER_sequence_profiles:" in result


def test_stage_control_is_still_rewritten(tmp_path, monkeypatch):
    """Positive control: a genuine *_OCCUPIER stage CONTROL must still be rewritten."""
    monkeypatch.delenv("DELFIN_OCC_ROOT", raising=False)
    from delfin.occupier_sequences import (
        append_sequence_overrides,
        remove_existing_sequence_blocks,
    )

    stage = tmp_path / "initial_OCCUPIER"
    stage.mkdir()
    control = stage / "CONTROL.txt"
    control.write_text(_ROOT_CONTROL, encoding="utf-8")

    remove_existing_sequence_blocks(control, force=True)
    append_sequence_overrides(control, _BUNDLE)
    text = control.read_text(encoding="utf-8")

    assert "OCCUPIER_sequence_profiles:" not in text  # template stripped
    assert "# AUTO sequence overrides" in text        # auto block written


def test_process_jobs_path_anchors_work_dir():
    """build_occupier_process_jobs must no longer call a bare run_OCCUPIER() under a
    local lock; it must anchor via pushd + work_dir like the other chdir sites."""
    import inspect
    import delfin.workflows.engine.occupier as eo

    src = inspect.getsource(eo.build_occupier_process_jobs)
    assert "run_OCCUPIER(work_dir=" in src
    assert "pushd(stage_dir)" in src
    assert "getattr(build_occupier_process_jobs, '_cwd_lock'" not in src
