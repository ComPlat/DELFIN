"""A bug report reaches the shared archive only with consent.

The report is the session's transcript: the system prompt with the
person's memories, tool outputs, copies of touched files. Writing it
locally is the person's act (the button, /bug, ``delfin-agent bug
report``); copying it to ``<transfer.remote_path>/AGENT_BUGS``, where the
archive's maintainers read it, is a second disclosure. Before, it
happened whenever a transfer host was configured -- a setting made for
moving job files -- including under ``--yes``, which nobody reads.

``agent.bug_upload``: ``ask`` (default) copies only after a yes for this
report, ``never`` keeps every report local, ``always`` is the standing
opt-in for unattended runs. The gate sits in ``push_report_to_remote``,
the one function that runs rsync, so a caller that forgets to ask
uploads nothing.
"""
from __future__ import annotations

import pytest

from delfin.agent import bug_filing
from delfin.agent import bug_report as br

TRANSFER = {"host": "login.cluster.example", "user": "someone",
            "remote_path": "/archive", "port": 22}


@pytest.fixture
def report(tmp_path):
    return br.write_bug_report(chat_messages=[], settings={},
                               archive_dir=tmp_path / "archive",
                               session_id="sess1234")


@pytest.fixture
def runs(monkeypatch):
    calls = []

    class _OK:
        returncode = 0
        stderr = ""

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return _OK()

    monkeypatch.setattr(br.subprocess, "run", fake_run)
    return calls


def _settings(policy=None):
    out = {"transfer": dict(TRANSFER), "agent": {}}
    if policy is not None:
        out["agent"]["bug_upload"] = policy
    return out


def _push(report, **kw):
    return br.push_report_to_remote(
        report, host=TRANSFER["host"], user=TRANSFER["user"],
        remote_path=TRANSFER["remote_path"], port=22, **kw)


def test_a_configured_transfer_alone_uploads_nothing(
        report, runs, monkeypatch):
    """Called as every caller on main calls it: no consent argument, the
    settings read from the settings file."""
    _load(monkeypatch, _settings())
    ok, why = _push(report)
    assert runs == [] and ok is False
    assert "agent.bug_upload" in why and "local" in why


def test_never_keeps_the_report_local_even_with_a_yes(report, runs):
    ok, why = _push(report, settings=_settings("never"), consent=True)
    assert ok is False and runs == []
    assert "never" in why


@pytest.mark.parametrize("policy", [None, "ask", "ASK ", "sometimes"])
def test_ask_uploads_after_a_yes_for_this_report(report, runs, policy):
    ok, where = _push(report, settings=_settings(policy), consent=True)
    assert ok is True and len(runs) == 2
    assert where.endswith(report.name)


def test_always_is_the_standing_opt_in_for_unattended_runs(report, runs):
    ok, _ = _push(report, settings=_settings("always"))
    assert ok is True and len(runs) == 2


def test_an_unknown_value_is_ask_not_always(report, runs):
    assert br.upload_policy(_settings("yes please")) == "ask"
    assert br.upload_policy({}) == "ask"
    ok, _ = _push(report, settings=_settings("yes please"))
    assert ok is False and runs == []


def _load(monkeypatch, settings):
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: settings)


def test_the_filing_flow_asks_and_a_no_keeps_it_local(
        report, runs, monkeypatch):
    _load(monkeypatch, _settings())
    asked = []

    def ask(question):
        asked.append(question)
        return "n"

    line = bug_filing.upload(report, ask=ask)
    assert runs == [] and len(asked) == 1
    assert TRANSFER["host"] in asked[0] and "[y/N]" in asked[0]
    assert "local" in line


def test_the_filing_flow_uploads_after_a_yes(report, runs, monkeypatch):
    _load(monkeypatch, _settings())
    line = bug_filing.upload(report, ask=lambda _q: "y")
    assert len(runs) == 2 and "copied" in line


@pytest.mark.parametrize("err", [EOFError, KeyboardInterrupt])
def test_no_answer_is_a_no(report, runs, monkeypatch, err):
    _load(monkeypatch, _settings())

    def ask(_q):
        raise err

    bug_filing.upload(report, ask=ask)
    assert runs == []


def test_without_anyone_to_ask_it_stays_local_and_says_so(
        report, runs, monkeypatch):
    """--yes and every unattended caller: no question can be asked, so
    nothing leaves unless agent.bug_upload is 'always'."""
    _load(monkeypatch, _settings())
    line = bug_filing.upload(report, ask=None)
    assert runs == []
    assert "always" in line and str(report) in line
    _load(monkeypatch, _settings("always"))
    bug_filing.upload(report, ask=None)
    assert len(runs) == 2


def test_the_preview_says_where_the_report_would_go(monkeypatch, tmp_path):
    _load(monkeypatch, _settings())
    text = bug_filing.preview({}, tmp_path)
    assert TRANSFER["host"] in text and "ask" in text
    _load(monkeypatch, _settings("never"))
    assert "stays on this machine" in bug_filing.preview({}, tmp_path)


def _settings_tab(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    monkeypatch.setenv("HOME", str(tmp_path))
    from delfin.dashboard import tab_settings
    from delfin.dashboard.context import DashboardContext
    for n in ("calc", "archive", "office"):
        (tmp_path / n).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda _s: None
    tab, _ = tab_settings.create_tab(ctx)

    def walk(w):
        yield w
        for c in getattr(w, "children", ()) or ():
            yield from walk(c)
    return list(walk(tab))


def test_the_settings_tab_sets_bug_upload(tmp_path, monkeypatch):
    """The setting that decides whether a report leaves the machine is
    chosen where the archive and the transfer host are set, from the
    three values, and is saved into agent.bug_upload."""
    from delfin import user_settings as us
    data = us.load_settings()
    data.setdefault("agent", {})["bug_upload"] = "never"
    us.save_settings(data)
    parts = _settings_tab(tmp_path, monkeypatch)
    (dd,) = [w for w in parts if type(w).__name__ == "Dropdown"
             and "upload" in str(getattr(w, "description", "")).lower()]
    assert {v for _label, v in dd.options} == set(br.UPLOAD_POLICIES)
    assert dd.value == "never"
    dd.value = "always"
    (save,) = [w for w in parts if type(w).__name__ == "Button"
               and w.description == "Save Settings"]
    assert not save.disabled
    save.click()
    assert us.load_settings()["agent"]["bug_upload"] == "always"
