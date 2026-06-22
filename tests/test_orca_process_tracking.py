from delfin import orca


class _FakeManager:
    def __init__(self):
        self.registered = []
        self.unregistered = []
        self._shutdown_requested = type("_ShutdownFlag", (), {"is_set": lambda self: False})()

    def register_subprocess(self, process, *, label="", cwd=None):
        token = f"token-{len(self.registered) + 1}"
        self.registered.append((token, process.pid, label, cwd))
        return token

    def unregister_subprocess(self, token):
        self.unregistered.append(token)


class _FakeProcess:
    def __init__(self, pid=4242, returncode=0):
        self.pid = pid
        self.returncode = returncode

    def wait(self, timeout=None):
        return self.returncode


def test_run_orca_subprocess_unregisters_tracked_process_on_success(monkeypatch, tmp_path):
    output_path = tmp_path / "job.out"
    input_path = tmp_path / "job.inp"
    input_path.write_text("! SP\n", encoding="utf-8")

    manager = _FakeManager()

    monkeypatch.setattr(orca, "get_global_manager", lambda: manager)
    monkeypatch.setattr(orca, "_should_monitor_orca_progress", lambda: False)
    monkeypatch.setattr(orca, "_prepare_orca_environment", lambda scratch_subdir: {})
    monkeypatch.setattr(orca, "_check_orca_success", lambda output_file: True)
    monkeypatch.setattr(orca, "_copy_densitiesinfo", lambda *args, **kwargs: None)
    monkeypatch.setattr(orca, "_ensure_process_group_terminated", lambda process, grace_timeout=2.0: None)
    monkeypatch.setattr(orca.subprocess, "Popen", lambda *args, **kwargs: _FakeProcess())

    assert orca._run_orca_subprocess("orca", str(input_path), str(output_path)) is True
    assert [token for token, *_ in manager.registered] == ["token-1"]
    assert manager.unregistered == ["token-1"]


def test_run_orca_subprocess_unregisters_tracked_process_on_failure(monkeypatch, tmp_path):
    """The zombie bug: ORCA fails (returncode!=0) but unregister must still be called."""
    output_path = tmp_path / "job.out"
    input_path = tmp_path / "job.inp"
    input_path.write_text("! SP\n", encoding="utf-8")

    manager = _FakeManager()

    monkeypatch.setattr(orca, "get_global_manager", lambda: manager)
    monkeypatch.setattr(orca, "_should_monitor_orca_progress", lambda: False)
    monkeypatch.setattr(orca, "_prepare_orca_environment", lambda scratch_subdir: {})
    monkeypatch.setattr(orca, "_ensure_process_group_terminated", lambda process, grace_timeout=2.0: None)
    monkeypatch.setattr(orca.subprocess, "Popen", lambda *args, **kwargs: _FakeProcess(returncode=1))

    assert orca._run_orca_subprocess("orca", str(input_path), str(output_path)) is False
    assert manager.unregistered == ["token-1"], "unregister must be called even when ORCA fails"


def test_run_orca_subprocess_unregisters_on_exception(monkeypatch, tmp_path):
    """If Popen itself raises, unregister must still be called for any prior registration."""
    input_path = tmp_path / "job.inp"
    input_path.write_text("! SP\n", encoding="utf-8")
    output_path = tmp_path / "job.out"

    manager = _FakeManager()

    def _boom(*args, **kwargs):
        raise OSError("orca binary not found")

    monkeypatch.setattr(orca, "get_global_manager", lambda: manager)
    monkeypatch.setattr(orca, "_should_monitor_orca_progress", lambda: False)
    monkeypatch.setattr(orca, "_prepare_orca_environment", lambda scratch_subdir: {})
    monkeypatch.setattr(orca.subprocess, "Popen", _boom)

    # Popen raises before register_subprocess is ever called, so nothing to unregister.
    # The important thing is that _run_orca_subprocess doesn't hang or leak.
    assert orca._run_orca_subprocess("orca", str(input_path), str(output_path)) is False
    assert manager.registered == []
    assert manager.unregistered == []


def test_copy_densitiesinfo_swallows_internal_errors(monkeypatch, tmp_path):
    monkeypatch.setattr(orca, "get_orca_scratch_dir", lambda: (_ for _ in ()).throw(RuntimeError("boom")))

    orca._copy_densitiesinfo(str(tmp_path / "job.inp"), None, None)
