"""The step that sends a bug report out has one name.

Input: a written report directory, the ``agent.bug_upload`` policy and an
``ask``. Output: the line ``bug_filing.upload`` hands back, and exactly the
transfers the consent gate allowed -- ``push_report_to_remote`` is counted,
its rsync never runs (the stand-in in its place).

Semantics: ``upload`` is the only spelling of "copy it to the remote
archive"; the earlier ``push_if_configured`` alias is deleted. Rejected
alternative: keeping the alias so old callers keep working -- it was never in
the module's ``__all__``, so there was no published surface to keep working.

The invariant test fails at the baseline (the second spelling is still
defined); the two behaviour tests pass at both arms -- they are the
never-worse guard that the one remaining spelling still sends after consent
and sends nothing without it.
"""
from __future__ import annotations

from delfin.agent import bug_filing
from delfin.agent import bug_report as br

TRANSFER = {"host": "login.cluster.example", "user": "someone",
            "remote_path": "/archive", "port": 22}


def _settings(policy=None):
    s = {"transfer": dict(TRANSFER)}
    if policy is not None:
        s["agent"] = {"bug_upload": policy}
    return s


def test_a_remote_copy_has_one_spelling_so_an_audit_finds_every_door():
    assert not hasattr(bug_filing, "push_if_configured")
    assert callable(bug_filing.upload)


def test_the_one_spelling_sends_only_after_consent(monkeypatch, tmp_path):
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: _settings("ask"))
    sent = []

    def stand_in(report_dir, **kw):
        sent.append(kw["consent"])
        return True, (f"{TRANSFER['user']}@{TRANSFER['host']}:"
                      f"{TRANSFER['remote_path']}/AGENT_BUGS/"
                      f"{report_dir.name}")

    monkeypatch.setattr(br, "push_report_to_remote", stand_in)
    line = bug_filing.upload(tmp_path, ask=lambda _q: "y")
    assert sent == [True]
    assert "copied to the remote archive" in line


def test_without_consent_the_one_spelling_sends_nothing(monkeypatch, tmp_path):
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: _settings("ask"))

    def boom(*a, **kw):
        raise AssertionError("a report left without consent")

    monkeypatch.setattr(br, "push_report_to_remote", boom)
    line = bug_filing.upload(tmp_path, ask=None)
    assert "nobody was asked" in line and "always" in line
    line = bug_filing.upload(tmp_path, ask=lambda _q: "n")
    assert "stays local" in line