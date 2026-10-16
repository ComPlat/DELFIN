"""``bug triage`` prints no credential: not from the settings, not from
an exception, not from a report's own text on the withheld line.

Input: settings that carry a fake API key in the places a key is kept or
mistyped (``agent.api_key``, a provider entry, ``bug_watcher.model`` and
``bug_watcher.provider`` set to the key), and an engine whose failure
message carries the key.

Output: :func:`bug_watcher.destination`, :func:`bug_watcher.send_line`,
the consent question and every line ``bug triage`` prints are free of the
key. A provider or model name is shown only when it is shaped like one
(``[A-Za-z0-9._:/@+-]{1,80}``, no 24-character alphanumeric run, nothing
the secret scrubber would change); an error message only when the
scrubber would leave it unchanged, otherwise its exception type alone.
The ``[WITHHELD]`` line prints the refusal reason, one of the fixed
sentences of :func:`bug_watcher.send_refusal`, never ``Triage.summary``
(which holds the model's analysis on the other paths; CodeQL alert 72,
py/clear-text-logging-sensitive-data, traced that flow to this line).
"""
from __future__ import annotations

import argparse
import json

import pytest

from delfin.agent import bug_watcher as bw

FAKE = "sk-" + "Zq7XvB2mN9pL4kR8tW1yC6dF3hJ5sA0e"   # shape of a provider key


def _settings(policy="ask"):
    return {
        "agent": {
            "api_key": FAKE,
            "provider": "kit",
            "model": "glm-4.6",
            "providers": {"kit": {"api_key": FAKE}},
            "bug_watcher": {"send_reports": policy, "model": FAKE,
                            "provider": FAKE},
        },
    }


def _report(archive, name="20260101-000000_u_sess_ab12"):
    d = archive / name
    d.mkdir(parents=True)
    (d / "report.json").write_text(json.dumps({
        "schema": "delfin-bug-report/1", "description": "it crashed",
        "model": "m", "mode": "solo"}), encoding="utf-8")
    (d / "report.md").write_text("# Bug report\n\ntranscript\n",
                                 encoding="utf-8")
    return d


def test_destination_never_names_a_key():
    s = _settings()
    assert FAKE not in bw.destination(s)
    for policy in bw.SEND_POLICIES:
        assert FAKE not in bw.send_line(_settings(policy))


def test_a_model_name_is_still_shown():
    s = {"agent": {"provider": "kit", "model": "openai/gpt-oss-120b",
                   "bug_watcher": {}}}
    assert bw.destination(s) == "kit / openai/gpt-oss-120b"


@pytest.mark.parametrize("value", [
    FAKE, "a" * 81, "model name with spaces", "m\nINJECTED",
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ012345",
])
def test_a_value_not_shaped_like_a_name_is_not_shown(value):
    s = {"agent": {"model": value, "bug_watcher": {}}}
    assert value not in bw.destination(s)


def test_the_consent_question_never_names_a_key():
    asked = []
    bw.ask_consent(1, _settings(), ask=lambda q: asked.append(q) or "n")
    assert asked and FAKE not in asked[0]


def test_an_error_carrying_a_key_is_printed_by_its_type_only():
    assert FAKE not in bw.printable_error(
        f"AuthenticationError: invalid key {FAKE}")
    assert bw.printable_error(
        f"AuthenticationError: invalid key {FAKE}").startswith(
            "AuthenticationError")
    plain = "TimeoutError: read timed out after 600 s"
    assert bw.printable_error(plain) == plain


def _cli(tmp_path, monkeypatch, capsys, *, policy, answer, engine):
    import delfin.user_settings as us
    from delfin.agent import cli
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: _settings(policy))
    monkeypatch.setattr(bw, "_default_engine_factory", engine)
    monkeypatch.setattr("sys.stdin.isatty", lambda: True, raising=False)
    monkeypatch.setattr("builtins.input", lambda _q="": answer)
    cli.cmd_bug(argparse.Namespace(bug_action="triage", archive=str(tmp_path),
                                   name=""))
    return capsys.readouterr()


def _engine(fail: bool):
    def factory(repo_root, settings=None):
        class _E:
            def stream_response(self, user_message="", on_token=None, **kw):
                if fail:
                    raise RuntimeError(f"401 Unauthorized for key {FAKE}")
                if on_token:
                    on_token("ROOT CAUSE: x")
        return _E()
    return factory


@pytest.mark.parametrize("policy,answer,fail", [
    ("ask", "n", False),       # withheld
    ("never", "", False),      # withheld
    ("ask", "y", False),       # analysed
    ("always", "", True),      # engine failure
])
def test_bug_triage_output_carries_no_key(tmp_path, monkeypatch, capsys,
                                          policy, answer, fail):
    _report(tmp_path)
    out = _cli(tmp_path, monkeypatch, capsys, policy=policy, answer=answer,
               engine=_engine(fail))
    assert FAKE not in out.out and FAKE not in out.err


def test_the_withheld_line_prints_the_reason_not_the_summary(
        tmp_path, monkeypatch, capsys):
    _report(tmp_path)
    reason = bw.send_refusal(consent=False, settings=_settings("never"))

    def fake_triage(d, consent=False, **kw):
        t = bw.Triage(report_name=d.name, report_path=str(d),
                      summary=f"model analysis quoting {FAKE}")
        t.withheld = True
        t.withheld_reason = reason
        return t
    monkeypatch.setattr(bw, "triage_report", fake_triage)
    out = _cli(tmp_path, monkeypatch, capsys, policy="never", answer="",
               engine=_engine(False))
    assert f"[WITHHELD] ({reason})" in out.out
    assert FAKE not in out.out


def test_a_withheld_triage_names_its_reason(tmp_path):
    d = _report(tmp_path)
    t = bw.triage_report(d, settings=_settings("never"))
    assert t.withheld
    assert t.withheld_reason == bw.send_refusal(
        consent=False, settings=_settings("never"))
