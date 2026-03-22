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
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
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
    assert captured["env"]["DELFIN_VOILA_ROOT_DIR"] == str(root_dir.resolve())

    stdout = capsys.readouterr().out
    assert "Starting DELFIN Dashboard on http://0.0.0.0:9001" in stdout


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

    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_prepare_voila_env", lambda open_browser: {})
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
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

    assert "--no-browser" in captured["cmd"]
    assert "--ServerApp.disable_check_xsrf=True" in captured["cmd"]
    assert "--Voila.open_browser=True" not in captured["cmd"]


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
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila.subprocess, "run", fake_run)
    def fake_popen(cmd, env):
        captured["cmd"] = cmd
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main(["--open-browser"])
    except SystemExit as exc:
        assert exc.code == 0

    assert "--Voila.open_browser=True" in captured["cmd"]
    assert "--ServerApp.disable_check_xsrf=True" in captured["cmd"]
    assert "--no-browser" not in captured["cmd"]
