import importlib.util
import signal

import pytest
import sys
import types
from pathlib import Path


_ROOT = Path(__file__).resolve().parents[1] / "delfin"
_MODULE_PATH = _ROOT / "dashboard" / "local_runner.py"

#: The three ``main`` installs; see the fixture below.
_ENTRY_POINT_SIGNALS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)

if "delfin" not in sys.modules:
    _PKG = types.ModuleType("delfin")
    _PKG.__path__ = [str(_ROOT)]
    _PKG.__version__ = "test"
    sys.modules["delfin"] = _PKG

if "delfin.dashboard" not in sys.modules:
    _DASHBOARD_PKG = types.ModuleType("delfin.dashboard")
    _DASHBOARD_PKG.__path__ = [str(_ROOT / "dashboard")]
    sys.modules["delfin.dashboard"] = _DASHBOARD_PKG

_SPEC = importlib.util.spec_from_file_location("delfin.dashboard.local_runner", _MODULE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Could not load local_runner module from {_MODULE_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)


def test_run_command_starts_child_in_own_session(monkeypatch, tmp_path):
    calls = {}

    class FakeProc:
        pid = 4321

        def wait(self):
            return 7

    def fake_popen(cmd, cwd=None, start_new_session=None):
        calls["cmd"] = cmd
        calls["cwd"] = cwd
        calls["start_new_session"] = start_new_session
        return FakeProc()

    monkeypatch.setattr(_MODULE.subprocess, "Popen", fake_popen)

    code = _MODULE._run_command(["echo", "test"], cwd=tmp_path)

    assert code == 7
    assert calls["cwd"] == str(tmp_path)
    assert calls["start_new_session"] is True


def test_main_preserves_signal_exit_code(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(_MODULE, "_configure_environment", lambda: None)
    monkeypatch.setattr(_MODULE, "_print_job_banner", lambda *args, **kwargs: None)
    monkeypatch.setattr(_MODULE, "_run_mode", lambda mode: (_ for _ in ()).throw(SystemExit(124)))

    written = []
    monkeypatch.setattr(_MODULE, "_write_exit_code", lambda code: written.append(code))

    try:
        _MODULE.main([])
    except SystemExit as exc:
        assert exc.code == 124
    else:
        raise AssertionError("Expected SystemExit")

    assert written[-1] == 124


def test_handle_termination_forwards_signal_and_exits(monkeypatch):
    forwarded = []
    written = []
    monkeypatch.setattr(_MODULE, "_signal_active_child", lambda signum: forwarded.append(signum))
    monkeypatch.setattr(_MODULE, "_write_exit_code", lambda code: written.append(code))

    try:
        _MODULE._handle_termination(signal.SIGTERM, None)
    except SystemExit as exc:
        assert exc.code == 124
    else:
        raise AssertionError("Expected SystemExit")

    assert forwarded == [signal.SIGTERM]
    assert written == [124]


def test_run_mode_orca_triggers_nmr_postprocess(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    inp_path = tmp_path / "mol_NMR.inp"
    inp_path.write_text("! test\n%EPRNMR\nEND\n", encoding="utf-8")

    monkeypatch.setenv("DELFIN_INP_FILE", inp_path.name)
    monkeypatch.setattr(_MODULE, "_resolve_orca_bin", lambda: "/opt/orca/orca")
    run_calls = []
    monkeypatch.setattr(
        "delfin.orca.run_orca",
        lambda inp_file, out_file: run_calls.append((inp_file, out_file)) or True,
    )

    rc = _MODULE._run_mode("orca")

    assert rc == 0
    assert run_calls == [(inp_path.name, "mol_NMR.out")]


def test_run_mode_orca_skips_nmr_postprocess_for_non_nmr_input(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    inp_path = tmp_path / "plain.inp"
    inp_path.write_text("! test\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")

    monkeypatch.setenv("DELFIN_INP_FILE", inp_path.name)
    monkeypatch.setattr(_MODULE, "_resolve_orca_bin", lambda: "/opt/orca/orca")
    run_calls = []
    monkeypatch.setattr(
        "delfin.orca.run_orca",
        lambda inp_file, out_file: run_calls.append((inp_file, out_file)) or True,
    )

    rc = _MODULE._run_mode("orca")

    assert rc == 0
    assert run_calls == [(inp_path.name, "plain.out")]


def test_run_mode_censo_anmr_passes_resume_flag(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    xyz_path = tmp_path / "mol.xyz"
    xyz_path.write_text("2\nx\nH 0 0 0\nH 0 0 1\n", encoding="utf-8")

    monkeypatch.setenv("DELFIN_XYZ_FILE", str(xyz_path))
    monkeypatch.setenv("DELFIN_WORKFLOW_LABEL", "mol_CENSO_ANMR")
    monkeypatch.setenv("DELFIN_CENSO_NMR_SOLVENT", "chcl3")
    monkeypatch.setenv("DELFIN_CENSO_NMR_CHARGE", "0")
    monkeypatch.setenv("DELFIN_CENSO_NMR_MULTIPLICITY", "1")
    monkeypatch.setenv("DELFIN_CENSO_NMR_MHZ", "400")
    monkeypatch.setenv("DELFIN_CENSO_NMR_RESUME", "1")

    captured = {}

    def fake_run_and_tee(cmd, path):
        captured["cmd"] = cmd
        captured["path"] = path
        return 0

    monkeypatch.setattr(_MODULE, "_run_and_tee", fake_run_and_tee)

    rc = _MODULE._run_mode("censo_anmr")

    assert rc == 0
    assert "--resume" in captured["cmd"]


@pytest.fixture(autouse=True)
def _the_handlers_this_file_installs_do_not_outlive_it():
    """Put SIGINT, SIGTERM and SIGHUP back after every test here.

    Input: the three dispositions before the test. Output: the same three
    after it. Semantics: nothing this file calls changes what a signal
    does to the rest of the run.

    `main` is a process entry point and installs
    `_handle_termination` for all three, which is right for
    `python -m delfin.dashboard.local_runner` and wrong inside pytest:
    `monkeypatch` restores attributes, not `signal.signal`, so calling
    `main` here left the handler installed for every test that followed.
    Measured: after this file, `signal.getsignal` returned
    `local_runner._handle_termination` for SIGINT, SIGTERM and SIGHUP;
    with a control file that does not call `main`, SIG_DFL.

    What it cost: the handler writes `.exit_code_<job id>` into
    `Path.cwd()` and raises SystemExit(124). In a suite run the cwd is
    the checkout, so one signal 11000 tests later put `.exit_code_0`
    there and the run failed on the guard against writing into the
    checkout -- blamed on whichever test happened to be last. Ctrl-C
    after this file did the same instead of interrupting pytest.
    """
    before = {num: signal.getsignal(num) for num in _ENTRY_POINT_SIGNALS}
    try:
        yield
    finally:
        for num, handler in before.items():
            if handler is not None:
                signal.signal(num, handler)


def test_main_installs_the_handlers_and_this_file_gives_them_back():
    """The leak, asserted from both ends.

    `main` must install them -- a job runner that ignores SIGTERM cannot
    stop its child -- and this file must not pass them on. Written out
    because the symptom was a file in the checkout and a failure
    attributed to an unrelated test, which is not a thing anyone traces
    back to a signal disposition.
    """
    seen = [signal.getsignal(num) for num in _ENTRY_POINT_SIGNALS]
    assert all(h is not _MODULE._handle_termination for h in seen), (
        "the fixture did not put the handlers back before this test")

    # Nothing of the run itself: `main` resolves a mode and would start a
    # real job. The handlers are installed before any of that, which is
    # the whole point -- a runner must be stoppable from its first
    # instruction -- so stubbing the body is enough to reach them.
    import io
    import contextlib
    monkeypatch = pytest.MonkeyPatch()
    try:
        monkeypatch.setattr(_MODULE, "_run_mode", lambda mode: 0)
        monkeypatch.setattr(_MODULE, "_configure_environment", lambda: None)
        monkeypatch.setattr(_MODULE, "_write_exit_code", lambda code: None)
        with contextlib.redirect_stdout(io.StringIO()):
            assert _MODULE.main([]) == 0
    finally:
        monkeypatch.undo()
    for num in _ENTRY_POINT_SIGNALS:
        assert signal.getsignal(num) is _MODULE._handle_termination, (
            f"main no longer installs a handler for {num}; a job runner "
            "that cannot be told to stop leaves its child behind")
    # Left installed on purpose: the fixture takes them away again, and
    # test_the_handlers_do_not_reach_the_next_test is what says so.


def test_the_handlers_do_not_reach_the_next_test():
    """Runs after the test above and sees SIG_DFL, not the runner's handler.

    The order is the assertion: pytest runs the tests in a file in the
    order they are written, so this one is the next test, and it is the
    one that would have seen the leak.
    """
    for num in _ENTRY_POINT_SIGNALS:
        handler = signal.getsignal(num)
        assert handler is not _MODULE._handle_termination, (
            f"{num} still runs local_runner._handle_termination, which "
            "writes .exit_code_* into the working directory and exits 124")


def test_the_exit_marker_lands_where_the_reader_looks(tmp_path, monkeypatch):
    """`_write_exit_code` writes into the process's working directory.

    That is the job directory when the runner is launched the way
    `backend_local` launches it, which is the only way it is launched --
    and `backend_local` reads the marker as
    `Path(job_dir) / f".exit_code_{job_id}"`. The two agree only through
    the cwd, so the cwd is what this pins: called from somewhere else,
    it writes somewhere else, and that is how a marker reached the
    checkout.
    """
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("DELFIN_JOB_ID", "77")
    _MODULE._write_exit_code(3)
    assert (tmp_path / ".exit_code_77").read_text().strip() == "3"
    assert sorted(p.name for p in tmp_path.iterdir()) == [".exit_code_77"]
