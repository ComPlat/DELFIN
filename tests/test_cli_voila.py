import sys
from pathlib import Path

from delfin import cli_voila


def test_stage_notebook_keeps_paths_inside_root(tmp_path):
    root_dir = tmp_path / "home"
    root_dir.mkdir()
    notebook = root_dir / "delfin_dashboard.ipynb"
    notebook.write_text("{}", encoding="utf-8")

    staged = cli_voila._stage_notebook_under_root(str(notebook), str(root_dir))

    assert staged == str(notebook.resolve())


def test_stage_notebook_copies_packaged_notebook_under_root(tmp_path):
    root_dir = tmp_path / "home"
    source_dir = tmp_path / "package"
    root_dir.mkdir()
    source_dir.mkdir()

    notebook = source_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    staged = Path(
        cli_voila._stage_notebook_under_root(str(notebook), str(root_dir))
    )

    assert staged == root_dir / "delfin_voila_runtime" / notebook.name
    assert staged.read_text(encoding="utf-8") == '{"cells":[]}'


def test_main_stages_out_of_root_notebook_before_launch(monkeypatch, tmp_path, capsys):
    root_dir = tmp_path / "home"
    package_dir = tmp_path / "package"
    root_dir.mkdir()
    package_dir.mkdir()

    notebook = package_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    captured = {}

    def fake_run(cmd, env, check):
        captured["cmd"] = cmd
        captured["env"] = env
        captured["check"] = check

        class Result:
            returncode = 0

        return Result()

    class FakeProc:
        def wait(self):
            return 0

    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_prepare_voila_env", lambda open_browser: {})
    monkeypatch.setattr(cli_voila, "_get_voila_static_root", lambda: "/tmp/voila-static")
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila, "_wait_for_port", lambda host, port, timeout=10.0: False)
    monkeypatch.setattr(cli_voila.subprocess, "run", fake_run)
    def fake_popen(cmd, env):
        captured["cmd"] = cmd
        captured["env"] = env
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main(["--no-browser", "--port", "9001"])
    except SystemExit as exc:
        assert exc.code == 0

    staged_path = root_dir / "delfin_voila_runtime" / notebook.name
    assert staged_path.exists()
    assert captured["cmd"][:4] == [
        sys.executable,
        "-m",
        "voila",
        str(staged_path),
    ]
    assert "--no-browser" in captured["cmd"]
    assert f"--port=9001" in captured["cmd"]
    assert f"--Voila.root_dir={root_dir.resolve()}" in captured["cmd"]
    assert (
        "--VoilaConfiguration.file_allowlist=.*\\.(png|jpg|gif|svg|js|css|html|ico)"
        in captured["cmd"]
    )
    assert "--ServerApp.disable_check_xsrf=True" in captured["cmd"]
    # Token is passed via env var (not CLI) for security — /proc/PID/cmdline
    # is world-readable, so the token must not appear in the command line.
    assert "--ServerApp.token=" not in captured["cmd"]
    assert "JUPYTER_TOKEN" in captured["env"]
    assert len(captured["env"]["JUPYTER_TOKEN"]) > 0
    assert "--ServerApp.websocket_ping_interval=30000" in captured["cmd"]
    assert "--ServerApp.websocket_ping_timeout=30000" in captured["cmd"]
    assert "--Voila.static_root=/tmp/voila-static" in captured["cmd"]
    assert captured["env"]["DELFIN_VOILA_ROOT_DIR"] == str(root_dir.resolve())
    assert "--VoilaConfiguration.preheat_kernel=False" in captured["cmd"]

    stdout = capsys.readouterr().out
    assert "Starting DELFIN Dashboard on http://localhost:9001" in stdout
    assert "?token=" in stdout  # URL includes the token for convenience


def test_main_reports_missing_voila_module(monkeypatch, capsys):
    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: False)

    try:
        cli_voila.main([])
    except SystemExit as exc:
        assert exc.code == 1

    stderr = capsys.readouterr().err
    assert "Error: voila is not installed." in stderr
    assert f"{sys.executable} -m pip install voila" in stderr


def test_main_defaults_to_no_browser(monkeypatch, tmp_path):
    root_dir = tmp_path / "home"
    package_dir = tmp_path / "package"
    root_dir.mkdir()
    package_dir.mkdir()

    notebook = package_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    captured = {}

    def fake_run(cmd, env, check):
        captured["cmd"] = cmd

        class Result:
            returncode = 0

        return Result()

    class FakeProc:
        def wait(self):
            return 0

    monkeypatch.delenv("TERM_PROGRAM", raising=False)
    monkeypatch.delenv("BROWSER", raising=False)
    monkeypatch.delenv("VSCODE_IPC_HOOK_CLI", raising=False)
    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_prepare_voila_env", lambda open_browser: {})
    monkeypatch.setattr(cli_voila, "_get_voila_static_root", lambda: "/tmp/voila-static")
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila, "_wait_for_port", lambda host, port, timeout=10.0: False)
    monkeypatch.setattr(cli_voila.subprocess, "run", fake_run)

    def fake_popen(cmd, env):
        captured["cmd"] = cmd
        captured["env"] = env
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main([])
    except SystemExit as exc:
        assert exc.code == 0

    assert "--no-browser" in captured["cmd"]
    assert "--Voila.ip=127.0.0.1" in captured["cmd"]
    assert "--ServerApp.disable_check_xsrf=True" in captured["cmd"]
    # Token via env var, not CLI
    assert "--ServerApp.token=" not in captured["cmd"]
    assert "JUPYTER_TOKEN" in captured["env"]
    assert "--ServerApp.websocket_ping_interval=30000" in captured["cmd"]
    assert "--ServerApp.websocket_ping_timeout=30000" in captured["cmd"]
    assert "--Voila.static_root=/tmp/voila-static" in captured["cmd"]
    assert "--Voila.open_browser=True" not in captured["cmd"]


def test_main_defaults_to_open_browser_in_vscode(monkeypatch, tmp_path, capsys):
    root_dir = tmp_path / "home"
    package_dir = tmp_path / "package"
    root_dir.mkdir()
    package_dir.mkdir()

    notebook = package_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    captured = {}
    observed = {}

    def fake_run(cmd, env, check):
        captured["cmd"] = cmd

        class Result:
            returncode = 0

        return Result()

    class FakeProc:
        def wait(self):
            return 0

    monkeypatch.setenv("TERM_PROGRAM", "vscode")
    monkeypatch.delenv("BROWSER", raising=False)
    monkeypatch.delenv("VSCODE_IPC_HOOK_CLI", raising=False)
    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_get_voila_static_root", lambda: "/tmp/voila-static")
    prepared_env = {"TERM_PROGRAM": "vscode"}

    def fake_prepare_voila_env(open_browser):
        observed["open_browser"] = open_browser
        return prepared_env

    browser_urls = []

    monkeypatch.setattr(cli_voila, "_prepare_voila_env", fake_prepare_voila_env)
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila, "_wait_for_port", lambda host, port, timeout=10.0: True)
    monkeypatch.setattr(cli_voila, "_open_browser_url", lambda url, env: browser_urls.append(url))
    monkeypatch.setattr(cli_voila.subprocess, "run", fake_run)

    def fake_popen(cmd, env):
        captured["cmd"] = cmd
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main([])
    except SystemExit as exc:
        assert exc.code == 0

    assert observed["open_browser"] is True
    assert browser_urls == ["http://localhost:8866"]
    assert "--Voila.open_browser=True" in captured["cmd"]
    assert "--no-browser" not in captured["cmd"]
    assert "--Voila.ip=127.0.0.1" in captured["cmd"]

    stdout = capsys.readouterr().out
    assert "Starting DELFIN Dashboard on http://localhost:8866" in stdout


def test_main_can_explicitly_open_browser(monkeypatch, tmp_path):
    root_dir = tmp_path / "home"
    package_dir = tmp_path / "package"
    root_dir.mkdir()
    package_dir.mkdir()

    notebook = package_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    captured = {}

    def fake_run(cmd, env, check):
        captured["cmd"] = cmd

        class Result:
            returncode = 0

        return Result()

    class FakeProc:
        def wait(self):
            return 0

    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_prepare_voila_env", lambda open_browser: {})
    monkeypatch.setattr(cli_voila, "_get_voila_static_root", lambda: "/tmp/voila-static")
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila.subprocess, "run", fake_run)
    def fake_popen(cmd, env):
        captured["cmd"] = cmd
        captured["env"] = env
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main(["--open-browser"])
    except SystemExit as exc:
        assert exc.code == 0

    assert "--Voila.open_browser=True" in captured["cmd"]
    assert "--ServerApp.disable_check_xsrf=True" in captured["cmd"]
    # Token via env var, not CLI
    assert "--ServerApp.token=" not in captured["cmd"]
    assert "JUPYTER_TOKEN" in captured["env"]
    assert "--ServerApp.websocket_ping_interval=30000" in captured["cmd"]
    assert "--ServerApp.websocket_ping_timeout=30000" in captured["cmd"]
    assert "--Voila.static_root=/tmp/voila-static" in captured["cmd"]
    assert "--no-browser" not in captured["cmd"]


def test_no_token_disables_auth(monkeypatch, tmp_path, capsys):
    """--no-token should pass --ServerApp.token= in CLI (empty = disabled)."""
    root_dir = tmp_path / "home"
    package_dir = tmp_path / "package"
    root_dir.mkdir()
    package_dir.mkdir()

    notebook = package_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    captured = {}

    class FakeProc:
        def wait(self):
            return 0

    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_prepare_voila_env", lambda open_browser: {})
    monkeypatch.setattr(cli_voila, "_get_voila_static_root", lambda: None)
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila, "_wait_for_port", lambda host, port, timeout=10.0: False)

    def fake_popen(cmd, env):
        captured["cmd"] = cmd
        captured["env"] = env
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main(["--no-token", "--no-browser"])
    except SystemExit as exc:
        assert exc.code == 0

    # With --no-token: empty token in CLI (auth disabled), NOT in env var
    assert "--ServerApp.token=" in captured["cmd"]
    assert "JUPYTER_TOKEN" not in captured["env"]

    stderr = capsys.readouterr().err
    assert "WARNING" in stderr
