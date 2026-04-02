import importlib.util
import signal
import sys
import types
from pathlib import Path


_ROOT = Path(__file__).resolve().parents[1] / "delfin"
_MODULE_PATH = _ROOT / "dashboard" / "local_runner.py"

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
