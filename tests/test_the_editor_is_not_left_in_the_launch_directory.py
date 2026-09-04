"""Thirty megabytes of Ketcher were found sitting in somebody's archive folder.

The editor is a browser application, so it has to be fetched over HTTP, so it
has to sit somewhere the server hands out.  Until there was a route of its own
the only such place was Voila's root directory -- and that root is whatever
directory the dashboard was started in.  Start it in ``archive`` and the
editor is copied into ``archive/.delfin/ketcher``; start it somewhere else
tomorrow and there are two.  Nothing ever took one back.

``delfin-voila`` now publishes one directory instead, ``~/.delfin/served``,
and the editor is loaded from where it is kept.  Nothing is written into the
launch directory at all -- and a copy an older DELFIN left in one is taken
back, once the same editor is kept somewhere that does not move.

The published directory holds the editor and nothing else on purpose: the
route jupyter_server serves it on carries no token (it is the one the login
page loads its own assets from), and ``~/.delfin`` is where the credentials,
the settings and the agent's memory live.
"""

from __future__ import annotations

import inspect
import pathlib

import pytest

from delfin.dashboard import ketcher


@pytest.fixture
def home(tmp_path, monkeypatch):
    """A home of this test's own, so the machine's real editor is not found."""
    where = tmp_path / "home"
    where.mkdir()
    monkeypatch.setenv("HOME", str(where))
    monkeypatch.setattr(pathlib.Path, "home", lambda: where)
    monkeypatch.delenv(ketcher._URL_ENV, raising=False)
    monkeypatch.delenv("DELFIN_VOILA_ROOT_DIR", raising=False)
    return where


def _an_editor(folder: pathlib.Path, version: str = "3.17.0") -> pathlib.Path:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "index.html").write_text("<html>ketcher</html>", encoding="utf-8")
    (folder / ketcher._STAMP).write_text(version, encoding="utf-8")
    return folder


# ---------------------------------------------------------------------------
# where it is loaded from
# ---------------------------------------------------------------------------
def test_it_is_published_from_where_it_is_kept(home, tmp_path, monkeypatch):
    """With a route of its own there is nothing to copy and nowhere to copy
    it to: the directory the browser loads is the directory it lives in."""
    root = tmp_path / "archive"
    root.mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    monkeypatch.setenv(ketcher._URL_ENV, "/static/ketcher")
    _an_editor(ketcher.stored_directory())

    assert ketcher.app_directory() == ketcher.stored_directory()
    assert ketcher.installed_version() == "3.17.0"
    assert ketcher.app_url() == "/static/ketcher/index.html?v=3.17.0"
    assert not (root / ".delfin").exists(), (
        "the launch directory was written into anyway"
    )


def test_the_url_does_not_change_when_the_dashboard_is_started_elsewhere(
    home, tmp_path, monkeypatch
):
    """It used to carry the path it was copied to, which is the launch
    directory -- so the same install answered on two different URLs."""
    monkeypatch.setenv(ketcher._URL_ENV, "/static/ketcher")
    _an_editor(ketcher.stored_directory())

    urls = set()
    for name in ("first", "second"):
        root = tmp_path / name
        root.mkdir()
        monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
        urls.add(ketcher.app_url())
    assert urls == {"/static/ketcher/index.html?v=3.17.0"}


def test_without_a_route_it_still_goes_under_the_voila_root(
    home, tmp_path, monkeypatch
):
    """Only the launcher can add a route, so a dashboard started any other
    way -- a bare ``voila``, a notebook server somebody already had -- must
    still find its editor."""
    root = tmp_path / "wherever it was started"
    root.mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    _an_editor(ketcher.stored_directory())

    assert ketcher.app_directory() == root / ".delfin" / "ketcher"
    assert ketcher.installed_version() == "3.17.0"
    assert ketcher.app_url().startswith("/voila/files/.delfin/ketcher/index.html")


def test_a_route_is_taken_only_when_it_looks_like_one(home, monkeypatch):
    monkeypatch.setenv(ketcher._URL_ENV, "static/ketcher")
    assert ketcher.url_prefix() is None, "a URL path starts with a slash"
    monkeypatch.setenv(ketcher._URL_ENV, "   ")
    assert ketcher.url_prefix() is None
    monkeypatch.setenv(ketcher._URL_ENV, "/static/ketcher/")
    assert ketcher.url_prefix() == "/static/ketcher", "no doubled slash"


# ---------------------------------------------------------------------------
# what was left behind
# ---------------------------------------------------------------------------
def test_a_copy_left_in_a_launch_directory_is_taken_back(home, tmp_path, monkeypatch):
    """The one this whole change is about: an archive folder with thirty
    megabytes of JavaScript in it, beside the calculations."""
    root = tmp_path / "archive"
    left = _an_editor(root / ".delfin" / "ketcher")
    _an_editor(ketcher.stored_directory())
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    monkeypatch.setenv(ketcher._URL_ENV, "/static/ketcher")

    assert ketcher.installed_version() == "3.17.0"
    assert not left.exists()
    assert not (root / ".delfin").exists(), "the folder that held it stayed"
    assert root.exists(), "only what DELFIN put there"


def test_the_folder_that_held_it_stays_if_anything_else_is_in_it(
    home, tmp_path, monkeypatch
):
    root = tmp_path / "archive"
    _an_editor(root / ".delfin" / "ketcher")
    kept_by_somebody = root / ".delfin" / "settings.json"
    kept_by_somebody.write_text("{}", encoding="utf-8")
    _an_editor(ketcher.stored_directory())
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    monkeypatch.setenv(ketcher._URL_ENV, "/static/ketcher")

    ketcher.installed_version()

    assert not (root / ".delfin" / "ketcher").exists()
    assert kept_by_somebody.is_file()


def test_a_folder_that_is_not_ours_is_not_touched(home, tmp_path, monkeypatch):
    """No stamp, no page: not something this ever put there."""
    root = tmp_path / "archive"
    theirs = root / ".delfin" / "ketcher"
    theirs.mkdir(parents=True)
    (theirs / "notes.txt").write_text("mine", encoding="utf-8")
    _an_editor(ketcher.stored_directory())
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    monkeypatch.setenv(ketcher._URL_ENV, "/static/ketcher")

    ketcher.installed_version()

    assert (theirs / "notes.txt").is_file()


def test_the_only_copy_there_is_is_never_the_one_removed(home, tmp_path, monkeypatch):
    """Somebody who fetched it before the route existed must not be asked to
    fetch thirty megabytes again: it is taken in first, and taken back only
    once it is safely kept."""
    root = tmp_path / "archive"
    left = _an_editor(root / ".delfin" / "ketcher", "3.16.0")
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    monkeypatch.setenv(ketcher._URL_ENV, "/static/ketcher")

    assert ketcher.installed_version() == "3.16.0", "it asked to fetch again"
    assert (ketcher.stored_directory() / "index.html").is_file()
    assert not left.exists()


def test_a_store_from_before_the_route_is_moved_rather_than_fetched_again(
    home, monkeypatch
):
    """It was kept at ``~/.delfin/ketcher``; the published directory is a
    level down from there, because everything in that one is handed out."""
    monkeypatch.setenv(ketcher._URL_ENV, "/static/ketcher")
    was = _an_editor(home / ".delfin" / "ketcher", "3.15.0")

    assert ketcher.installed_version() == "3.15.0"
    assert (ketcher.stored_directory() / "index.html").is_file()
    assert not was.exists()


# ---------------------------------------------------------------------------
# what is in the published directory
# ---------------------------------------------------------------------------
def test_the_editor_is_the_only_thing_published(home):
    """The route carries no token -- jupyter_server serves the login page's
    own assets on it -- so the directory it points at may hold nothing that
    is not meant for anybody who can reach the port.  ``~/.delfin`` holds the
    credentials, the settings and the agent's memory."""
    published = ketcher.published_root()

    assert published != home / ".delfin", "that would publish the credentials"
    assert published.parent == home / ".delfin"
    assert ketcher.stored_directory().parent == published

    published.mkdir(parents=True)
    _an_editor(ketcher.stored_directory())
    assert [p.name for p in published.iterdir()] == ["ketcher"]


def test_a_download_is_never_unpacked_into_the_published_directory():
    """Half an unpacked zip under a route that hands out whatever is there."""
    body = inspect.getsource(ketcher.install)

    assert "staging = Path.home() / '.delfin'" in body
    assert "folder.parent / ('ketcher-download" not in body


def test_the_editor_is_not_copied_onto_itself():
    """Served where it is kept, the old last step -- remove the served copy,
    then copy the kept one into it -- would delete the editor and then copy
    it from nowhere."""
    body = inspect.getsource(ketcher.install)

    assert "if folder != kept:" in body
    assert body.index("shutil.move(str(unpacked), str(kept))") < body.index(
        "shutil.copytree(kept, folder)"
    )


# ---------------------------------------------------------------------------
# the launcher
# ---------------------------------------------------------------------------
def test_the_launcher_publishes_that_directory_and_says_so(home):
    from delfin import cli_voila

    assert cli_voila._prepare_published_dir() == str(ketcher.published_root())
    assert ketcher.published_root().is_dir()

    source = pathlib.Path(cli_voila.__file__).read_text(encoding="utf-8")
    launch = source.split("_published_dir = _prepare_published_dir()")[1][:600]
    assert "--ServerApp.extra_static_paths=" in launch
    assert '"/static/ketcher"' in launch
    assert ketcher._URL_ENV in launch


def test_the_launcher_really_puts_it_on_the_command_line(home, tmp_path, monkeypatch):
    """The whole thing turns on two strings reaching the server and the
    kernel: the route, and the URL the dashboard is to ask on."""
    from delfin import cli_voila

    root_dir = tmp_path / "started here"
    package_dir = tmp_path / "package"
    root_dir.mkdir()
    package_dir.mkdir()
    notebook = package_dir / "delfin_dashboard.ipynb"
    notebook.write_text('{"cells":[]}', encoding="utf-8")

    seen = {}

    class Proc:
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
    monkeypatch.setattr(
        cli_voila, "_wait_for_port", lambda host, port, timeout=10.0: False
    )
    monkeypatch.setattr(
        cli_voila.subprocess, "run", lambda cmd, env, check: type("R", (), {"returncode": 0})()
    )

    def popen(cmd, env=None, **kwargs):
        seen["cmd"], seen["env"] = cmd, env
        return Proc()

    monkeypatch.setattr(cli_voila.subprocess, "Popen", popen)
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root_dir))
    # main() remembers the launch directory in the environment; set it here so
    # it is put back with the rest of what this test changed.
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path))

    try:
        cli_voila.main([])
    except SystemExit as exc:
        assert exc.code == 0

    published = str(ketcher.published_root())
    assert f'--ServerApp.extra_static_paths=["{published}"]' in seen["cmd"]
    assert seen["env"][ketcher._URL_ENV] == "/static/ketcher"
    assert not (root_dir / ".delfin").exists()


def test_the_variable_is_set_only_when_the_route_was_actually_added(home):
    """Told there is a route when there is none, the editor would be asked
    for on a URL that answers 404 and no other route would be tried."""
    from delfin import cli_voila

    source = pathlib.Path(cli_voila.__file__).read_text(encoding="utf-8")
    launch = source.split("_published_dir = _prepare_published_dir()")[1][:600]
    assert launch.lstrip().startswith("if _published_dir:")
