"""An API key can be set from the Settings tab.

The agent does not start without one. Until now a key could be set in
two places, neither of which a person reaching their dashboard through
the Windows launcher has:

  * an export in the shell before the dashboard process started, and
  * ``delfin-agent credentials set NAME`` on a terminal -- the launcher's
    working shell is off by default for new connections.

The dashboard READ keys (``tab_agent._provider_key`` ->
``credentials.load_credential``) but could not write one:
``set_credential`` had exactly one caller in the whole tree, the CLI.
What the page offered instead was a sentence telling the person to type
``os.environ['KIT_TOOLBOX_API_KEY']='...'`` into a notebook cell, which
puts a credential in the notebook and in the kernel's history.

Measured 2026-11-19. The field is in the Settings tab because that is
where settings and repairs belong.

What these tests hold onto, beyond "it stores": the value never comes
back out. A status line, an error message or a log line that quotes a
key is the failure that matters here, and it is the kind that no one
notices until the key is somewhere else.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from delfin.dashboard import tab_settings

KEY = "KIT_TOOLBOX_API_KEY"
SECRET = "stand-in-value-never-a-real-key-7f3a9c1e"


@pytest.fixture
def own_store(monkeypatch, tmp_path):
    """A credential store of this test's own; never the person's."""
    home = tmp_path / "home"
    (home / ".delfin").mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    from delfin.agent import credentials as cred
    monkeypatch.setattr(cred, "_DEFAULT_PATH",
                        home / ".delfin" / "credentials.json")
    monkeypatch.setattr(cred, "configured_store", lambda: "file")
    for name, _label in tab_settings.API_KEY_NAMES:
        monkeypatch.delenv(name, raising=False)
    return home


def test_a_key_set_here_is_the_key_the_agent_finds(own_store):
    """The whole point: not "a write returned", but "the agent's own
    lookup now answers"."""
    stored, message = tab_settings.store_api_key(KEY, SECRET)

    assert stored, message
    from delfin.agent.credentials import load_credential
    assert load_credential(KEY) == SECRET


def test_the_key_works_in_this_dashboard_without_a_restart(own_store):
    """A person who just pasted their key expects the next New Cycle to
    work, not the next login."""
    tab_settings.store_api_key(KEY, SECRET)

    from delfin.dashboard.tab_agent import _provider_key
    assert _provider_key(KEY) == SECRET


def test_the_message_never_carries_the_key(own_store):
    stored, message = tab_settings.store_api_key(KEY, SECRET)

    assert stored
    assert SECRET not in message
    assert SECRET[:8] not in message


def test_a_refusal_never_carries_the_key_either(own_store):
    """The error path is where a value usually leaks: it is the one that
    gets written while thinking about something else."""
    _ok, message = tab_settings.store_api_key("NOT_A_DELFIN_KEY", SECRET)

    assert not _ok
    assert SECRET not in message


def test_a_store_that_raises_does_not_quote_what_it_was_given(own_store,
                                                              monkeypatch):
    from delfin.agent import credentials as cred

    def boom(name, value, **kw):
        raise RuntimeError(f"disk full while writing {value}")

    monkeypatch.setattr(cred, "set_credential", boom)

    stored, message = tab_settings.store_api_key(KEY, SECRET)

    assert not stored
    assert SECRET not in message, (
        "the store's own exception text quoted the key, and this message "
        "goes straight onto the page")


def test_nothing_about_the_key_reaches_the_log(own_store, caplog):
    with caplog.at_level(logging.DEBUG):
        tab_settings.store_api_key(KEY, SECRET)

    assert SECRET not in caplog.text


def test_an_empty_field_stores_nothing_and_says_so(own_store):
    stored, message = tab_settings.store_api_key(KEY, "   ")

    assert not stored and "nothing" in message.lower()
    from delfin.agent.credentials import load_credential
    assert load_credential(KEY) == ""


def test_the_status_says_which_keys_are_missing(own_store):
    """A person whose agent will not start has to see WHICH key it is."""
    rows = dict(tab_settings.api_key_status_rows())

    assert set(rows) == {n for n, _ in tab_settings.API_KEY_NAMES}
    assert rows[KEY] == "not set"

    tab_settings.store_api_key(KEY, SECRET)

    assert dict(tab_settings.api_key_status_rows())[KEY].startswith("set")


def test_the_status_never_carries_a_value(own_store):
    tab_settings.store_api_key(KEY, SECRET)

    for name, state in tab_settings.api_key_status_rows():
        assert SECRET not in state
        assert SECRET[:8] not in state


def _settings_widgets(tmp_path):
    """The real Settings tab, built. Not a source check: the field has to
    exist on the page, and the button has to be wired to it."""
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda _s: None
    tab, _ = tab_settings.create_tab(ctx)

    def walk(widget):
        yield widget
        for child in getattr(widget, "children", ()) or ():
            yield from walk(child)

    return list(walk(tab))


def test_the_field_and_its_button_are_on_the_page(own_store, tmp_path):
    parts = _settings_widgets(tmp_path)

    buttons = [w for w in parts if type(w).__name__ == "Button"
               and "store key" in str(getattr(w, "description", "")).lower()]
    fields = [w for w in parts if type(w).__name__ == "Password"
              and "key" in str(getattr(w, "placeholder", "")).lower()]
    choices = [w for w in parts if type(w).__name__ == "Dropdown"
               and "api key" in str(getattr(w, "description", "")).lower()]

    assert buttons and fields and choices, (
        "the Settings tab must offer somewhere to put a key; without it "
        "the only routes are a shell export and the CLI")
    assert KEY in [v for _l, v in choices[0].options]


def test_clicking_store_on_the_built_page_stores_the_key(own_store,
                                                         tmp_path):
    """The click itself, through the real widgets."""
    parts = _settings_widgets(tmp_path)
    field = [w for w in parts if type(w).__name__ == "Password"
             and "key" in str(getattr(w, "placeholder", "")).lower()][0]
    button = [w for w in parts if type(w).__name__ == "Button"
              and "store key" in str(getattr(w, "description", "")).lower()][0]

    field.value = SECRET
    for handler in button._click_handlers.callbacks:
        handler(button)

    from delfin.agent.credentials import load_credential
    assert load_credential(KEY) == SECRET
    assert field.value == "", (
        "the field must be cleared after storing; a key left in a widget "
        "is a key in the page's state for the rest of the session")


def test_only_keys_delfin_reads_are_offered(own_store):
    """An arbitrary name would be stored and never read, which looks from
    the page exactly like success."""
    names = {n for n, _ in tab_settings.API_KEY_NAMES}

    assert KEY in names
    stored, _msg = tab_settings.store_api_key("SOME_OTHER_TOKEN", SECRET)
    assert not stored

    from delfin.agent.credentials import load_credential
    assert load_credential("SOME_OTHER_TOKEN") == ""


def test_one_persons_key_does_not_replace_the_accounts(tmp_path,
                                                        monkeypatch):
    """On a shared account the person's store is a LINK to the account's.

    ``remote_bootstrap.sh`` gives each person a home under ``~/people``
    and links the account's files into it, ``.delfin/credentials.json``
    among them, so that someone who has set no key of their own still
    finds the account's. Measured on the lab Raspberry Pi on 2026-11-19:
    every person's credentials file is such a link, and the account's
    file held the one shared KIT key.

    A write that FOLLOWED that link would replace the key for everyone
    on the box -- the other people and the machine's own agents -- from
    a button in the Settings tab. It does not: the store writes a
    private temporary file and renames it over the path, which replaces
    the LINK rather than its target. The person ends up with a file of
    their own holding what the account's held plus their key, and the
    account's file is untouched.

    Pinned here because the whole arrangement is load-bearing for that
    machine and nothing else holds it.
    """
    from delfin.agent import credentials as cred

    account = tmp_path / "account" / ".delfin"
    account.mkdir(parents=True)
    account.chmod(0o700)
    shared = account / "credentials.json"
    shared.write_text('{"KIT_TOOLBOX_API_KEY": "the-accounts-key"}',
                      encoding="utf-8")
    shared.chmod(0o600)

    person = tmp_path / "account" / "people" / "simone" / ".delfin"
    person.mkdir(parents=True)
    person.chmod(0o700)
    link = person / "credentials.json"
    link.symlink_to(shared)

    monkeypatch.setenv("HOME", str(person.parent))
    monkeypatch.setattr(cred, "_DEFAULT_PATH", link)
    monkeypatch.setattr(cred, "configured_store", lambda: "file")

    stored, _message = tab_settings.store_api_key(KEY, SECRET)

    assert stored
    assert not link.is_symlink(), (
        "the person's store must become a file of their own; while it is "
        "a link, every write lands in the account's file")
    assert json.loads(shared.read_text())[KEY] == "the-accounts-key", (
        "one person setting a key in Settings replaced the key the whole "
        "account uses")
    assert json.loads(link.read_text())[KEY] == SECRET
