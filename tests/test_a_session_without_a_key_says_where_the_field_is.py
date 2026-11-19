"""A session with no key says so, and says where the field is.

The Settings tab has had somewhere to put an API key since 2026-11-17
(section "Diagnosis & API key"). What it did not have was anything
pointing at it from where the trouble actually shows up.

A session that could not start showed a symptom instead: a 401 from the
model endpoint, or the model list quietly falling back to its hard-coded
entries under "KIT_TOOLBOX_API_KEY is not set in the dashboard process".
Both are true, and neither says "you have not entered a key yet". The
remedy on offer was to type ``os.environ['...']='...'`` into a notebook
cell -- which puts a credential in the notebook and in the kernel's
history, and assumes the reader knows what a notebook cell is.

The second half of this file is about the store rather than the message,
and is here because the Raspberry Pi arrangement makes it load-bearing:
each person's credentials file is a symlink to the account's, so a write
that followed the link would replace the key for everyone on the box.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path

from delfin.dashboard import tab_agent, tab_settings

KEY = "KIT_TOOLBOX_API_KEY"
SECRET = "stand-in-value-never-a-real-key-7f3a9c1e"


# ------------------- what the person is told when there is no key yet

def test_a_session_without_a_key_says_so_and_says_where():
    said = tab_agent.missing_key_notice("kit", lambda _name: "")

    assert "No API key yet" in said
    assert KEY in said, "name the key, so the dropdown choice is obvious"
    assert "Settings tab" in said and "Diagnosis & API key" in said
    assert "Store key" in said, "name the button, not a vague 'save'"
    assert "New Cycle" in said, "say how to start once it is stored"
    assert "notebook" not in said and "os.environ" not in said, (
        "a credential typed into a notebook cell lands in the notebook "
        "and in the kernel history")


def test_nothing_is_said_once_a_key_is_there():
    assert tab_agent.missing_key_notice("kit", lambda _n: "a-key") == ""


def test_a_provider_that_needs_no_key_is_not_nagged():
    """The CLI backend signs in through its own OAuth."""
    assert tab_agent.missing_key_notice("claude", lambda _n: "") == ""
    assert tab_agent.missing_key_notice("", lambda _n: "") == ""


def test_a_lookup_that_fails_does_not_cry_wolf():
    """Telling someone their key is missing when the lookup merely broke
    sends them to re-enter a key that is already there."""
    def boom(_name):
        raise RuntimeError("store unreadable")

    assert tab_agent.missing_key_notice("kit", boom) == ""


def test_the_notice_sends_people_to_a_field_that_can_take_their_key():
    """The notice names a key and sends the person to set it. If the
    Settings dropdown had no such entry, that would be a dead end."""
    src = inspect.getsource(tab_settings.create_tab)
    i = src.index("key_name = widgets.Dropdown(")
    offered = src[i:i + 300]

    for name in tab_agent.PROVIDER_KEY_NAMES.values():
        assert f"'{name}'" in offered, (
            f"the agent asks for {name} and the Settings tab cannot "
            "store it")


def test_no_route_still_points_at_a_notebook_cell():
    """Both model-fallback hints used to."""
    src = inspect.getsource(tab_agent.create_tab)

    assert "os.environ['" not in src
    assert "notebook cell" not in src


def test_the_backend_hint_no_longer_claims_keys_are_never_stored():
    """It said so while the same tab stored them."""
    src = inspect.getsource(tab_settings.create_tab)

    assert "Keys are never stored on disk" not in src
    assert "Diagnosis &amp; API key" in src


# --------------------------- the store, on a shared account

def _linked_person_store(tmp_path, monkeypatch, shared_json: str):
    """The Pi's layout: the person's store is a link to the account's."""
    from delfin.agent import credentials as cred

    account = tmp_path / "account" / ".delfin"
    account.mkdir(parents=True)
    account.chmod(0o700)
    shared = account / "credentials.json"
    shared.write_text(shared_json, encoding="utf-8")
    shared.chmod(0o600)

    person = tmp_path / "account" / "people" / "simone" / ".delfin"
    person.mkdir(parents=True)
    person.chmod(0o700)
    link = person / "credentials.json"
    link.symlink_to(shared)

    monkeypatch.setenv("HOME", str(person.parent))
    monkeypatch.setattr(cred, "_DEFAULT_PATH", link)
    monkeypatch.setattr(cred, "configured_store", lambda: "file")
    return cred, shared, link


def test_one_persons_key_does_not_replace_the_accounts(tmp_path,
                                                        monkeypatch):
    """``remote_bootstrap.sh`` gives each person a home under ``~/people``
    and links the account's files into it, ``.delfin/credentials.json``
    among them, so someone who has set no key of their own still finds
    the account's. Measured on the lab Raspberry Pi on 2026-11-19: every
    person's credentials file is such a link, and the account's file held
    the one shared KIT key.

    A write that FOLLOWED that link would replace the key for everyone on
    the box -- the other people and the machine's own agents -- from a
    button in the Settings tab. It does not: the store writes a private
    temporary file and renames it over the path, which replaces the LINK
    rather than its target.

    Pinned here because the arrangement is load-bearing for that machine
    and nothing else holds it.
    """
    cred, shared, link = _linked_person_store(
        tmp_path, monkeypatch, '{"KIT_TOOLBOX_API_KEY": "the-accounts-key"}')

    cred.set_credential(KEY, SECRET)

    assert not link.is_symlink(), (
        "the person's store must become a file of their own; while it is "
        "a link, every write lands in the account's file")
    assert json.loads(shared.read_text())[KEY] == "the-accounts-key", (
        "one person setting a key replaced the key the whole account uses")
    assert json.loads(link.read_text())[KEY] == SECRET


def test_what_was_in_the_accounts_file_is_not_lost(tmp_path, monkeypatch):
    """Breaking the link must carry the account's other credentials over,
    or a person setting one key loses every key they had through it."""
    cred, _shared, link = _linked_person_store(
        tmp_path, monkeypatch,
        '{"KIT_TOOLBOX_API_KEY": "shared", "OPENAI_API_KEY": "also-shared"}')

    cred.set_credential(KEY, SECRET)

    mine = json.loads(link.read_text())
    assert mine["OPENAI_API_KEY"] == "also-shared", (
        "the account's other keys were reachable through the link and "
        "must survive it being replaced")
    assert mine[KEY] == SECRET
