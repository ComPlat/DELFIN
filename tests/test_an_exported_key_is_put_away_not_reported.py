"""A key exported in the shell is put away, not reported every start.

DELFIN found the key, said so at every start, and left four manual steps:
store it, find the export, delete the line, open a new shell. On a
cluster the key sat in a group-readable ``~/.bashrc`` in the meantime
(2026-09-17). Finding a problem and handing it back is not the same as
fixing it.

  the key is taken into the store  0600, where DELFIN reads it anyway
  the export is commented out      with the reason, in the file itself
  a copy is kept                   and the copy is owner-only: it still
                                   holds the key
  a different stored value stops   which of the two is current is a
                                   question, not a chore
  nothing exported, nothing done   silent
"""

from __future__ import annotations

import os

import pytest

from delfin.agent import credentials as cred
from delfin.agent import process_guard


@pytest.fixture
def home(tmp_path):
    rc = tmp_path / ".bashrc"
    rc.write_text('# my shell\nexport KIT_TOOLBOX_API_KEY="sk-live-1234"\n'
                  'echo hello\n', encoding="utf-8")
    os.chmod(rc, 0o644)
    return tmp_path


def _secure(home, env=None, **kw):
    return cred.secure_exported_keys(
        ["KIT_TOOLBOX_API_KEY"],
        path=home / ".delfin" / "credentials.json",
        env=env if env is not None else {"KIT_TOOLBOX_API_KEY": "sk-live-1234"},
        home=home, **kw)


def test_the_key_lands_in_the_store(home):
    rows = _secure(home)
    assert rows[0]["action"] == "stored"
    store = home / ".delfin" / "credentials.json"
    assert store.exists()
    assert os.stat(store).st_mode & 0o777 == 0o600
    assert cred.load_credential("KIT_TOOLBOX_API_KEY", path=store) == "sk-live-1234"


def test_the_export_is_commented_out_with_its_reason(home):
    _secure(home)
    text = (home / ".bashrc").read_text(encoding="utf-8")
    assert '# export KIT_TOOLBOX_API_KEY="sk-live-1234"' in text
    assert "disabled by DELFIN" in text
    assert "echo hello" in text, "the rest of the file is untouched"


def test_the_copy_it_keeps_is_owner_only(home):
    """The copy still holds the key; inheriting a group-readable mode
    would move the exposure rather than end it."""
    _secure(home)
    backup = home / ".bashrc.delfin-backup"
    assert backup.exists()
    assert "sk-live-1234" in backup.read_text(encoding="utf-8")
    assert os.stat(backup).st_mode & 0o777 == 0o600


def test_a_stored_value_that_differs_is_left_alone(home):
    store = home / ".delfin" / "credentials.json"
    cred.set_credential("KIT_TOOLBOX_API_KEY", "sk-older-9999", path=store)
    rows = _secure(home)
    assert rows[0]["action"] == "differs"
    assert cred.load_credential("KIT_TOOLBOX_API_KEY", path=store) == "sk-older-9999"
    assert "export KIT_TOOLBOX_API_KEY" in (home / ".bashrc").read_text()
    assert rows[0]["cleaned"] == []


def test_nothing_exported_means_nothing_done(home):
    rows = _secure(home, env={})
    assert rows[0]["action"] == "absent"
    assert "export KIT_TOOLBOX_API_KEY" in (home / ".bashrc").read_text()


def test_a_commented_export_is_not_found_again(home):
    _secure(home)
    again = cred.exported_in_shell_files("KIT_TOOLBOX_API_KEY", home=home)
    assert again == [], "a line already commented out is not an export"


def test_the_line_it_reports_says_what_happened(home, monkeypatch):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-live-1234")
    monkeypatch.setattr(cred, "_DEFAULT_PATH", home / ".delfin" / "credentials.json")
    monkeypatch.setattr(cred.Path, "home", staticmethod(lambda: home))
    said = process_guard.put_exported_keys_away(["KIT_TOOLBOX_API_KEY"])
    assert "taken into" in said and "commented out" in said
    assert "Open a new shell" in said
    assert "sk-live-1234" not in said, "never the value itself"


def test_it_says_nothing_when_there_is_nothing_to_do(monkeypatch, tmp_path):
    monkeypatch.delenv("KIT_TOOLBOX_API_KEY", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    monkeypatch.setattr(cred.Path, "home", staticmethod(lambda: tmp_path))
    assert process_guard.put_exported_keys_away() == ""


# -- where a login session, not a shell file, hands it over -----------------

def test_the_session_files_are_searched_too(tmp_path):
    """No rc file named the key on the cluster, and every process still
    had it: a login session is given variables in places that are not
    shells."""
    (tmp_path / ".ssh").mkdir()
    (tmp_path / ".ssh" / "environment").write_text(
        "KIT_TOOLBOX_API_KEY=sk-from-ssh\n", encoding="utf-8")
    (tmp_path / ".config" / "environment.d").mkdir(parents=True)
    (tmp_path / ".config" / "environment.d" / "keys.conf").write_text(
        "KIT_TOOLBOX_API_KEY=sk-from-systemd\n", encoding="utf-8")

    found = {str(f) for f, _n, _l in
             cred.exported_in_shell_files("KIT_TOOLBOX_API_KEY", home=tmp_path)}
    assert str(tmp_path / ".ssh" / "environment") in found
    assert str(tmp_path / ".config" / "environment.d" / "keys.conf") in found


def test_a_session_file_is_commented_out_like_any_other(tmp_path):
    (tmp_path / ".ssh").mkdir()
    rc = tmp_path / ".ssh" / "environment"
    rc.write_text("KIT_TOOLBOX_API_KEY=sk-from-ssh\n", encoding="utf-8")
    _secure(tmp_path, env={"KIT_TOOLBOX_API_KEY": "sk-from-ssh"})
    text = rc.read_text(encoding="utf-8")
    assert text.startswith("# disabled by DELFIN")
    assert "# KIT_TOOLBOX_API_KEY=sk-from-ssh" in text


def test_an_export_from_nowhere_in_the_home_is_said_so(tmp_path, monkeypatch):
    """A site profile or a job script is outside reach: the key is put
    in the store and the remaining export is named, not invented away."""
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-live-1234")
    monkeypatch.setattr(cred, "_DEFAULT_PATH", tmp_path / ".delfin" / "c.json")
    monkeypatch.setattr(cred.Path, "home", staticmethod(lambda: tmp_path))
    monkeypatch.setattr(cred, "unset_in_systemd_user_environment",
                        lambda name: False)
    said = process_guard.put_exported_keys_away(["KIT_TOOLBOX_API_KEY"])
    assert "taken into" in said
    assert "still exported by something outside the home directory" in said
