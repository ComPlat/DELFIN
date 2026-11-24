"""A key check would have passed on the evening the models did not answer.

MEASURED 2026-11-23. Twenty-two swarm agents on Blablador ended within
minutes on a reply no model had produced ("There's an issue with the
selected model ..."), with 0 tokens. An authenticated request at the
same time listed 76 models, ours among them: the KEY WAS ACCEPTED. So a
live check of the kind that existed -- list the models, report the key
accepted -- would have said all was well. The question that tells a
working provider from a broken one is whether the MODEL ANSWERS, read
from the usage it reports, since a model that answered has read the
prompt.

Three more gaps, in the same place:

  * ``check_kit_toolbox_key_live`` was named in ``__all__`` and called by
    nothing in the product: no command could run it;
  * the doctor's whole network path (``fast=False``) was reached by none
    of its six callers -- CLI, REPL, agent tab, settings tab,
    prerequisites, the agent's own tool -- so MCP reachability was never
    probed either;
  * the credentials row listed KIT, Anthropic and OpenAI, and not
    Blablador, so a stored Blablador key was invisible to it.

``delfin-agent doctor --live`` opens the network path, and with it the
key checks and the new one. No test here reaches the network: urlopen
and the credential store are replaced, as in the KIT key test.
"""
from __future__ import annotations

import io
import json
import urllib.error

import pytest

from delfin import doctor


class _Answer(io.BytesIO):
    status = 200

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _key(monkeypatch, value="stand-in-key"):
    monkeypatch.setattr("delfin.agent.credentials.load_credential",
                        lambda name: value)


def _serve(monkeypatch, body=None, *, code=200, exc=None, seen=None):
    def _urlopen(request, timeout=None):
        if seen is not None:
            seen.append(request)
        if exc is not None:
            raise exc
        if code != 200:
            raise urllib.error.HTTPError(
                request.full_url, code, "no", {}, io.BytesIO(b"{}"))
        return _Answer(json.dumps(body or {}).encode())
    monkeypatch.setattr("urllib.request.urlopen", _urlopen)


GATEWAY = ("There's an issue with the selected model (alias-glm-5.3-flash). "
           "It may not exist or you may not have access to it.")


# --- the Blablador key, as KIT's is checked -------------------------------

def test_a_blablador_key_the_server_accepts_is_ok(monkeypatch):
    _key(monkeypatch)
    _serve(monkeypatch, {"data": [{"id": "a"}, {"id": "b"}]})
    out = doctor.check_blablador_key_live()
    assert out.status == doctor.OK and "2 models" in out.detail


def test_a_refused_blablador_key_is_refused(monkeypatch):
    _key(monkeypatch)
    _serve(monkeypatch, code=401)
    out = doctor.check_blablador_key_live()
    assert out.status != doctor.OK and "refused" in out.detail


def test_the_key_never_appears_in_the_result(monkeypatch):
    _key(monkeypatch, "glpat-must-never-be-shown")
    _serve(monkeypatch, code=401)
    out = doctor.check_blablador_key_live()
    assert "must-never-be-shown" not in f"{out.detail} {out.fix_hint}"


# --- whether the model answers: the check that would have caught it -------

def test_a_model_that_answers_is_ok(monkeypatch):
    _key(monkeypatch)
    _serve(monkeypatch, {"choices": [{"message": {"content": "OK"}}],
                         "usage": {"prompt_tokens": 14,
                                   "completion_tokens": 1}})
    out = doctor.check_model_answers_live("blablador", "alias-glm-5.3-flash")
    assert out.status == doctor.OK, out.detail


def test_a_reply_no_model_produced_is_not_ok(monkeypatch):
    """The measured evening: accepted, listed, and not answering."""
    _key(monkeypatch)
    _serve(monkeypatch, {"choices": [{"message": {"content": GATEWAY}}],
                         "usage": {"prompt_tokens": 0,
                                   "completion_tokens": 0}})
    out = doctor.check_model_answers_live("blablador", "alias-glm-5.3-flash")
    assert out.status != doctor.OK
    assert "without running the model" in out.detail
    assert "issue with the selected model" in out.detail


def test_a_reply_with_no_usage_at_all_is_not_ok(monkeypatch):
    _key(monkeypatch)
    _serve(monkeypatch, {"choices": [{"message": {"content": GATEWAY}}]})
    out = doctor.check_model_answers_live("blablador", "alias-glm-5.3-flash")
    assert out.status != doctor.OK


def test_a_reasoning_model_with_no_visible_text_still_answered(monkeypatch):
    """A thinking model can spend a small budget before writing a word;
    the prompt tokens say it ran."""
    _key(monkeypatch)
    _serve(monkeypatch, {"choices": [{"message": {"content": None}}],
                         "usage": {"prompt_tokens": 14,
                                   "completion_tokens": 64}})
    out = doctor.check_model_answers_live("blablador", "alias-glm-5.3-flash")
    assert out.status == doctor.OK, out.detail


def test_unreachable_is_not_refused(monkeypatch):
    """Telling somebody to replace a key that was never the problem is
    the mistake the first diagnosis of that evening made."""
    _key(monkeypatch)
    _serve(monkeypatch, exc=urllib.error.URLError("timed out"))
    out = doctor.check_model_answers_live("kit", "kit.glm-5.3-flash")
    assert out.status != doctor.OK
    assert "refused" not in out.detail
    assert "could not be checked" in out.detail


def test_the_request_goes_to_the_provider_it_names(monkeypatch):
    _key(monkeypatch)
    seen = []
    _serve(monkeypatch, {"usage": {"prompt_tokens": 3}}, seen=seen)
    doctor.check_model_answers_live("blablador", "m")
    assert seen and "blablador" in seen[0].full_url
    assert seen[0].full_url.endswith("/chat/completions")


def test_an_unknown_provider_is_named_not_guessed():
    out = doctor.check_model_answers_live("nobody-knows", "m")
    assert out.status != doctor.OK and "nobody-knows" in out.detail


# --- and the doctor can now be asked ---------------------------------------

def test_doctor_live_runs_the_network_path(monkeypatch):
    from delfin.agent import cli as agent_cli
    from delfin.agent import doctor as agent_doctor

    seen = {}

    def _run(workspace=None, *, settings=None, fast=True):
        seen["fast"] = fast
        return []

    monkeypatch.setattr(agent_doctor, "run_doctor", _run)
    agent_cli.main(["doctor", "--live"])
    assert seen.get("fast") is False


def test_without_live_the_doctor_stays_off_the_network(monkeypatch):
    from delfin.agent import cli as agent_cli
    from delfin.agent import doctor as agent_doctor

    seen = {}

    def _run(workspace=None, *, settings=None, fast=True):
        seen["fast"] = fast
        return []

    monkeypatch.setattr(agent_doctor, "run_doctor", _run)
    agent_cli.main(["doctor"])
    assert seen.get("fast", True) is True


def test_the_fast_report_names_the_live_check():
    """A check nobody is told about is a check nobody runs."""
    from delfin.agent import doctor as agent_doctor

    rows = agent_doctor._check_provider_live({"fast": True})
    text = json.dumps(rows)
    assert "doctor --live" in text


def test_the_credentials_row_knows_blablador():
    from delfin.agent import doctor as agent_doctor

    assert any(name == "BLABLADOR_API_KEY"
               for _label, name in agent_doctor._PROVIDER_KEYS)


@pytest.mark.parametrize("fast", [True, False])
def test_the_live_group_is_in_the_report(fast, monkeypatch):
    from delfin.agent import doctor as agent_doctor

    assert ("provider live", "_check_provider_live") in \
        agent_doctor._CHECK_ATTRS
