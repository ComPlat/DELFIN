"""A key that is present and a key that is accepted are two questions.

Input: the stored KIT-Toolbox credential. Output: two checks -- one that
says it is there, one that asks whether it still works.

Measured 2026-09-28: a benchmark trial was blocked for a whole session by
a stored key the server had expired. The doctor called it "configured",
because presence was all it looked at, and the run died much later at
engine start with a message about no key at all. The presence check now
says what it did not check, and the live one asks.

The live check is deliberately outside run_all: it reaches a network, so
it would make every doctor run wait on a remote host and could hang one
that only wanted to know what is installed.
"""

from __future__ import annotations

import urllib.error

import pytest

from delfin import doctor


def _no_stored_key(monkeypatch):
    monkeypatch.setattr("delfin.agent.credentials.load_credential",
                        lambda *a, **k: "")
    monkeypatch.delenv("KIT_TOOLBOX_API_KEY", raising=False)


def test_the_presence_check_says_what_it_did_not_check(monkeypatch):
    monkeypatch.setattr("delfin.agent.credentials.load_credential",
                        lambda *a, **k: "sk-whatever")
    out = doctor.check_kit_toolbox_key()
    assert out.status == doctor.OK
    assert "not checked against the server" in out.detail, (
        "a present key still reads as a working one")


def test_a_refused_key_is_reported_as_refused(monkeypatch):
    monkeypatch.setattr("delfin.agent.credentials.load_credential",
                        lambda *a, **k: "sk-expired")

    def _401(*a, **k):
        raise urllib.error.HTTPError("u", 401, "Unauthorized", {}, None)

    monkeypatch.setattr("urllib.request.urlopen", _401)
    out = doctor.check_kit_toolbox_key_live()
    assert out.status != doctor.OK
    assert "401" in out.detail and "refused" in out.detail
    assert "credentials set" in (out.fix_hint or ""), "no way out is named"


def test_unreachable_is_not_reported_as_refused(monkeypatch):
    """Replacing a key that was never the problem costs the user a trip to
    a web interface for nothing."""
    monkeypatch.setattr("delfin.agent.credentials.load_credential",
                        lambda *a, **k: "sk-fine")

    def _down(*a, **k):
        raise OSError("network unreachable")

    monkeypatch.setattr("urllib.request.urlopen", _down)
    out = doctor.check_kit_toolbox_key_live()
    assert "refused" not in out.detail
    assert "says nothing about the key" in (out.fix_hint or "")


def test_an_accepted_key_says_how_many_models(monkeypatch):
    monkeypatch.setattr("delfin.agent.credentials.load_credential",
                        lambda *a, **k: "sk-fine")

    class _Answer:
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def read(self):
            return b'{"data": [{"id": "a"}, {"id": "b"}]}'

    monkeypatch.setattr("urllib.request.urlopen", lambda *a, **k: _Answer())
    out = doctor.check_kit_toolbox_key_live()
    assert out.status == doctor.OK and "2 models" in out.detail


def test_no_key_is_not_a_server_question(monkeypatch):
    _no_stored_key(monkeypatch)
    out = doctor.check_kit_toolbox_key_live()
    assert "no KIT_TOOLBOX_API_KEY" in out.detail


def test_the_key_never_reaches_the_result(monkeypatch):
    """An expired key and a wrong key are the same answer here; printing
    either would put a secret in whatever reads this."""
    monkeypatch.setattr("delfin.agent.credentials.load_credential",
                        lambda *a, **k: "sk-SECRETVALUE-123")

    def _401(*a, **k):
        raise urllib.error.HTTPError("u", 401, "Unauthorized", {}, None)

    monkeypatch.setattr("urllib.request.urlopen", _401)
    out = doctor.check_kit_toolbox_key_live()
    assert "SECRETVALUE" not in (out.detail + (out.fix_hint or ""))


def test_the_live_check_stays_out_of_the_default_run():
    """run_all must not wait on a remote host."""
    import inspect

    src = inspect.getsource(doctor.run_all)
    assert "check_kit_toolbox_key_live" not in src, (
        "the network probe is in the default doctor run")
    assert "check_kit_toolbox_key" in src


def test_it_is_exported_so_something_can_call_it():
    assert "check_kit_toolbox_key_live" in doctor.__all__
