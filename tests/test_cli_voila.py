import os
import sys
from pathlib import Path

import pytest

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
    def fake_popen(cmd, env=None, **kwargs):
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
    # SECURITY: launch under jupyter-server + the voila extension so the token
    # is actually ENFORCED. Standalone `voila` ignores the token entirely
    # (dashboard, kernel API and kernel websocket are all reachable with no
    # auth). A refactor back to standalone `voila` must fail this test.
    assert captured["cmd"][:4] == [
        sys.executable,
        "-m",
        "jupyter",
        "server",
    ]
    # voila MUST be enabled under jupyter-server; the other auto-loading
    # extensions (jupyterlab/notebook/lsp/terminals) are explicitly DISABLED to
    # shrink the exposed HTTP surface. A refactor that drops voila or re-enables
    # the extras must fail here.
    _ext_arg = next(
        a for a in captured["cmd"] if a.startswith("--ServerApp.jpserver_extensions=")
    )
    assert "'voila': True" in _ext_arg
    for _disabled in ("jupyterlab", "notebook", "jupyter_lsp", "jupyter_server_terminals"):
        assert f"'{_disabled}': False" in _ext_arg
    staged_rel = staged_path.resolve().relative_to(root_dir.resolve()).as_posix()
    assert f"--ServerApp.default_url=/voila/render/{staged_rel}" in captured["cmd"]
    assert "--no-browser" in captured["cmd"]
    assert "--port=9001" in captured["cmd"]
    assert f"--ServerApp.root_dir={root_dir.resolve()}" in captured["cmd"]
    assert (
        "--VoilaConfiguration.file_allowlist=.*\\.(png|jpg|gif|svg|js|css|html|ico|pdf)"
        in captured["cmd"]
    )
    assert "--ServerApp.disable_check_xsrf=True" not in captured["cmd"]
    # Tracebacks OFF (info disclosure) and server terminals disabled (RCE surface).
    assert "--VoilaConfiguration.show_tracebacks=False" in captured["cmd"]
    assert "--show_tracebacks=True" not in captured["cmd"]
    assert "--ServerApp.terminals_enabled=False" in captured["cmd"]
    # Token is passed via env var (not CLI) for security — /proc/PID/cmdline
    # is world-readable, so the token must not appear in the command line.
    assert "--ServerApp.token=" not in captured["cmd"]
    assert "JUPYTER_TOKEN" in captured["env"]
    assert len(captured["env"]["JUPYTER_TOKEN"]) > 0
    assert "--ServerApp.websocket_ping_interval=30000" in captured["cmd"]
    assert "--ServerApp.websocket_ping_timeout=30000" in captured["cmd"]
    assert captured["env"]["DELFIN_VOILA_ROOT_DIR"] == str(root_dir.resolve())
    assert "--VoilaConfiguration.preheat_kernel=False" in captured["cmd"]

    stdout = capsys.readouterr().out
    assert "Starting DELFIN Dashboard on http://localhost:9001" in stdout
    # SECURITY: the raw token value MUST NOT appear in stdout — stdout may be
    # captured by journald/log files readable by other users on shared hosts.
    # In non-TTY mode, only the token-file path is printed; users read the
    # value with `cat $(file)`.
    _real_token = captured["env"]["JUPYTER_TOKEN"]
    assert _real_token not in stdout, "token value leaked to stdout"
    assert "voila-token-" in stdout, "token-file hint missing from stdout"


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

    def fake_popen(cmd, env=None, **kwargs):
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
    assert "--ServerApp.ip=127.0.0.1" in captured["cmd"]
    assert "--ServerApp.disable_check_xsrf=True" not in captured["cmd"]
    # Token via env var, not CLI
    assert "--ServerApp.token=" not in captured["cmd"]
    assert "JUPYTER_TOKEN" in captured["env"]
    assert "--ServerApp.websocket_ping_interval=30000" in captured["cmd"]
    assert "--ServerApp.websocket_ping_timeout=30000" in captured["cmd"]
    assert "--ServerApp.terminals_enabled=False" in captured["cmd"]
    assert "--ServerApp.open_browser=True" not in captured["cmd"]


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

    def fake_popen(cmd, env=None, **kwargs):
        captured["cmd"] = cmd
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main([])
    except SystemExit as exc:
        assert exc.code == 0

    assert observed["open_browser"] is True
    # Auto-open must land DIRECTLY in the dashboard: the render URL WITH the
    # token embedded, never the bare "/" (which would show the login page).
    assert len(browser_urls) == 1
    opened = browser_urls[0]
    assert opened.startswith("http://localhost:8866/voila/render/")
    assert "?token=" in opened
    # jupyter must NOT open its own (token-less) tab; the launcher controls it.
    assert "--no-browser" in captured["cmd"]
    assert "--ServerApp.open_browser=True" not in captured["cmd"]
    assert "--ServerApp.ip=127.0.0.1" in captured["cmd"]

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
    monkeypatch.setattr(cli_voila, "_wait_for_port", lambda host, port, timeout=10.0: True)
    browser_urls = []
    monkeypatch.setattr(cli_voila, "_open_browser_url", lambda url, env: browser_urls.append(url))
    monkeypatch.setattr(cli_voila.subprocess, "run", fake_run)
    def fake_popen(cmd, env=None, **kwargs):
        captured["cmd"] = cmd
        captured["env"] = env
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main(["--open-browser"])
    except SystemExit as exc:
        assert exc.code == 0

    # jupyter is told NOT to open (the launcher opens the token URL itself).
    assert "--no-browser" in captured["cmd"]
    assert "--ServerApp.open_browser=True" not in captured["cmd"]
    # The opened URL is the dashboard render URL WITH the token (direct login).
    assert browser_urls and "/voila/render/" in browser_urls[0]
    assert "?token=" in browser_urls[0]
    assert "--ServerApp.disable_check_xsrf=True" not in captured["cmd"]
    # Token via env var, not CLI
    assert "--ServerApp.token=" not in captured["cmd"]
    assert "JUPYTER_TOKEN" in captured["env"]
    assert "--ServerApp.websocket_ping_interval=30000" in captured["cmd"]
    assert "--ServerApp.websocket_ping_timeout=30000" in captured["cmd"]


def test_no_token_is_refused(monkeypatch, tmp_path, capsys):
    """--no-token is refused because Voilà exposes kernel/agent execution."""
    root_dir = tmp_path / "home"
    package_dir = tmp_path / "package"
    root_dir.mkdir()
    package_dir.mkdir()

    notebook = package_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main(["--no-token", "--no-browser"])
    except SystemExit as exc:
        assert exc.code == 2

    stderr = capsys.readouterr().err
    assert "--no-token is no longer allowed" in stderr


def test_short_custom_token_is_refused(monkeypatch, tmp_path, capsys):
    root_dir = tmp_path / "home"
    package_dir = tmp_path / "package"
    root_dir.mkdir()
    package_dir.mkdir()
    notebook = package_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main(["--token", "test", "--no-browser"])
    except SystemExit as exc:
        assert exc.code == 2

    stderr = capsys.readouterr().err
    assert "insecure --token" in stderr


def test_generated_token_is_never_printed(monkeypatch, tmp_path, capsys):
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
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    def fake_popen(cmd, env=None, **kwargs):
        captured["cmd"] = cmd
        captured["env"] = env
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)

    try:
        cli_voila.main(["--no-browser"])
    except SystemExit as exc:
        assert exc.code == 0

    token = captured["env"]["JUPYTER_TOKEN"]
    stdout = capsys.readouterr().out
    assert token not in stdout
    assert f"?token={token}" not in stdout
    assert "Append ?token=$(cat " in stdout


def test_remote_bind_is_refused_without_override(monkeypatch, tmp_path, capsys):
    root_dir = tmp_path / "home"
    package_dir = tmp_path / "package"
    root_dir.mkdir()
    package_dir.mkdir()
    notebook = package_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(notebook))
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    try:
        cli_voila.main(["--ip", "0.0.0.0", "--no-browser"])
    except SystemExit as exc:
        assert exc.code == 2

    stderr = capsys.readouterr().err
    assert "refusing to bind delfin-voila to a non-loopback address" in stderr


def test_remote_bind_requires_explicit_override(monkeypatch, tmp_path):
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
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))

    def fake_popen(cmd, env=None, **kwargs):
        captured["cmd"] = cmd
        captured["env"] = env
        return FakeProc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", fake_popen)

    try:
        cli_voila.main([
            "--ip", "0.0.0.0", "--allow-remote-bind", "--no-browser",
        ])
    except SystemExit as exc:
        assert exc.code == 0

    assert "--ServerApp.ip=0.0.0.0" in captured["cmd"]
    assert "JUPYTER_TOKEN" in captured["env"]


# ---------------------------------------------------------------------------
# Security hardening helpers
# ---------------------------------------------------------------------------

def test_dashboard_default_url_points_at_rendered_notebook(tmp_path):
    root = tmp_path / "root"
    (root / "delfin_voila_runtime").mkdir(parents=True)
    nb = root / "delfin_voila_runtime" / "delfin_dashboard.ipynb"
    nb.write_text("{}", encoding="utf-8")

    url = cli_voila._dashboard_default_url(str(nb), str(root))
    assert url == "/voila/render/delfin_voila_runtime/delfin_dashboard.ipynb"


def test_dashboard_default_url_encodes_spaces(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    nb = root / "my dashboard.ipynb"
    nb.write_text("{}", encoding="utf-8")

    url = cli_voila._dashboard_default_url(str(nb), str(root))
    assert url == "/voila/render/my%20dashboard.ipynb"


@pytest.mark.skipif(
    not hasattr(os, "geteuid"), reason="POSIX ownership check only"
)
def test_world_writable_notebook_source_is_rejected(tmp_path):
    """A world-writable notebook must be refused (search-path hijack guard)."""
    nb = tmp_path / "delfin_dashboard.ipynb"
    nb.write_text("{}", encoding="utf-8")

    os.chmod(nb, 0o644)
    assert cli_voila._is_safe_notebook_source(nb) is True

    os.chmod(nb, 0o666)  # world-writable
    assert cli_voila._is_safe_notebook_source(nb) is False


def test_find_notebook_ignores_unsafe_cwd_notebook(monkeypatch, tmp_path):
    """A world-writable cwd notebook must NOT be picked up; fall back to packaged."""
    if not hasattr(os, "geteuid"):
        pytest.skip("POSIX ownership check only")

    work = tmp_path / "work"
    work.mkdir()
    planted = work / "delfin_dashboard.ipynb"
    planted.write_text("{}", encoding="utf-8")
    os.chmod(planted, 0o666)

    monkeypatch.chdir(work)
    found = cli_voila._find_notebook()
    assert Path(found).resolve() != planted.resolve()
