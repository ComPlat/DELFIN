"""The working-directory box of a new session offers the folders a person
can pick, the way an editor's open dialog does.

Asked for on 2026-09-16: the box listed four fixed paths; everything else
had to be typed in full. It now lists the folders of the home directory,
and as a path is typed, the folders of that path or of its parent.
"""
import inspect
from types import SimpleNamespace

from delfin.dashboard import agent_sessions as AS


def _ctx(tmp_path):
    return SimpleNamespace(agent_dir=str(tmp_path / "agent"), calc_dir=str(tmp_path / "calc"),
                           repo_dir=str(tmp_path / "repo"))


def _home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    for name in ("DELFIN", "calc", "projects", ".hidden", "notes.txt"):
        (home / name).mkdir(parents=True) if name != "notes.txt" else None
    (home / "notes.txt").write_text("x")
    monkeypatch.setattr(AS.Path, "home", classmethod(lambda cls: home))
    monkeypatch.setattr(AS, "default_workspace", lambda ctx: str(home / "DELFIN"))
    return home


def test_nothing_typed_lists_the_folders_of_home(tmp_path, monkeypatch):
    home = _home(tmp_path, monkeypatch)
    got = AS.workspace_choices(_ctx(tmp_path), "")
    assert got[:4] == [str(home / "DELFIN"), str(tmp_path / "agent"), str(tmp_path / "calc"), str(tmp_path / "repo")]
    assert str(home / "projects") in got and str(home / "calc") in got
    assert str(home / ".hidden") not in got
    assert str(home / "notes.txt") not in got


def test_a_typed_prefix_completes_like_an_editor(tmp_path, monkeypatch):
    home = _home(tmp_path, monkeypatch)
    got = AS._folders_like(str(home / "pro"))
    assert got == [str(home / "projects")]
    (home / "projects" / "nmr").mkdir()
    assert AS._folders_like(str(home / "projects") + "/") == [str(home / "projects" / "nmr")]
    assert AS._folders_like(str(home / "projects")) == [str(home / "projects" / "nmr")]


def test_a_dot_shows_hidden_folders_and_nonsense_shows_nothing(tmp_path, monkeypatch):
    home = _home(tmp_path, monkeypatch)
    assert AS._folders_like(str(home) + "/.") == [str(home / ".hidden")]
    assert AS._folders_like("/no/such/dir/at/all") == []


def test_the_list_is_bounded(tmp_path, monkeypatch):
    home = _home(tmp_path, monkeypatch)
    for i in range(AS._PICKER_MAX + 50):
        (home / f"d{i:04d}").mkdir()
    assert len(AS.workspace_choices(_ctx(tmp_path), "")) <= AS._PICKER_MAX


def test_the_box_follows_what_is_typed():
    src = inspect.getsource(AS.create_tab)
    assert "workdir_box.observe(_refresh_folder_choices, names=\"value\")" in src
    assert "workspace_choices(ctx, typed)" in src


def test_the_session_bar_stays_on_screen_while_the_chat_scrolls():
    """Asked for on 2026-09-16: scrolling a long chat scrolled the session
    list away. The sidebar is sticky to the top of the page; the shell must
    not clip with `overflow: hidden`, which would make it the scroll box
    and end the stickiness."""
    css = AS._SIDEBAR_CSS
    i = css.index(".delfin-sessions {")
    block = css[i:css.index("}", i)]
    assert "position: sticky" in block and "top:" in block
    shell = css[css.index(".delfin-session-shell {"):]
    shell = shell[:shell.index("}")]
    assert "overflow-x: hidden" not in shell
