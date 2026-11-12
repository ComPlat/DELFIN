"""Exported provider keys are moved and shell files edited only with consent.

The hunt (12.11., finding 17): any chat/run with an exported key moved it
into the store and commented out the export line in the shell file, without
asking. gh, aws and kubectl only advise. Now: asked on a terminal (default
yes), never done without one, and DELFIN_KEEP_EXPORTED_KEYS=1 never asks.
"""
from delfin.agent import process_guard as pg


def test_without_a_terminal_nothing_is_changed(monkeypatch):
    monkeypatch.delenv("DELFIN_KEEP_EXPORTED_KEYS", raising=False)
    assert pg._may_put_keys_away(["KIT_TOOLBOX_API_KEY"], isatty=False) is False


def test_on_a_terminal_enter_means_yes_and_n_means_no(monkeypatch):
    monkeypatch.delenv("DELFIN_KEEP_EXPORTED_KEYS", raising=False)
    assert pg._may_put_keys_away(["K"], isatty=True, ask=lambda q: "") is True
    assert pg._may_put_keys_away(["K"], isatty=True, ask=lambda q: "n") is False


def test_the_opt_out_never_asks(monkeypatch):
    monkeypatch.setenv("DELFIN_KEEP_EXPORTED_KEYS", "1")
    asked = []
    assert pg._may_put_keys_away(["K"], isatty=True,
                                 ask=lambda q: asked.append(q) or "") is False
    assert asked == []


def test_declined_gives_advice_and_changes_nothing(monkeypatch, tmp_path):
    monkeypatch.setattr(pg, "_may_put_keys_away", lambda names, **k: False)
    called = []
    from delfin.agent import credentials as cred
    monkeypatch.setattr(cred, "secure_exported_keys",
                        lambda *a, **k: called.append(a) or [])
    said = pg.put_exported_keys_away(["KIT_TOOLBOX_API_KEY"])
    assert called == [] and said.startswith("Note:")
