"""The header's Dashboards menu: every kept dashboard, Open and End, New.

Switching between kept dashboards is one click, and so is starting
another one. The menu is built from the kept-session records
(``session.live_records``, the registry a resume follows) and, for the
windows and last activity, from the server's own kernel list; there is
no second registry.

  a row per kept dashboard   name, open windows, last activity
  Open                       the resume address, relative to this page
                             (whatever host and port the browser used),
                             in a new tab, without the token
  End                        after one confirmation: the record is
                             dropped and the server shuts the kernel down
                             (its kept agent sessions go to the service)
  this window                marked, no Open, no End (Keep dashboard is
                             the switch for it)
  another login node         named, nothing to click here
  New dashboard              the dashboard's own address, in a new tab

The heartbeat of a kept dashboard refreshes its record and never
re-creates one that was dropped: otherwise an End from another window
could be undone by the next beat.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timedelta, timezone

import pytest

pytest.importorskip("ipywidgets")

from delfin.dashboard import dashboards_menu as M  # noqa: E402
from delfin.dashboard import server_record as SR  # noqa: E402
from delfin.dashboard import session as S  # noqa: E402

RESUME = "/voila/render/delfin_voila_runtime/delfin_resume.ipynb"
DASH = "/voila/render/delfin_voila_runtime/delfin_dashboard.ipynb"


@pytest.fixture
def kept(monkeypatch, tmp_path):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "kept"))
    monkeypatch.setattr(S, "_server_stdout", lambda: None)
    monkeypatch.setenv(S.RESUME_PATH_ENV, RESUME)
    monkeypatch.setenv(S.DASHBOARD_PATH_ENV, DASH)
    monkeypatch.setenv("JUPYTER_TOKEN", "k" * 43)
    monkeypatch.setenv("DELFIN_VOILA_PORT", "8899")
    S.write_record("node-aaaa", kid="kid-a")
    S.write_record("node-bbbb", kid="kid-b")
    path = os.path.join(S.RECORD_DIR, "node-cccc.json")
    with open(path, "w", encoding="utf-8") as fh:
        import time
        json.dump({"session_name": "node-cccc", "kernel_id": "kid-c",
                   "host": "elsewhere", "pid": 1, "started_at": time.time(),
                   "updated_at": time.time()}, fh)
    recent = (datetime.now(timezone.utc) - timedelta(minutes=4)).isoformat()
    monkeypatch.setattr(SR, "kernels", lambda record, timeout=5.0: [
        {"id": "kid-a", "connections": 0, "last_activity": recent},
        {"id": "kid-b", "connections": 2, "last_activity": recent},
        {"id": "kid-now", "connections": 1, "last_activity": recent}])
    monkeypatch.setattr(S, "kernel_id", lambda: "kid-b")


def test_one_row_per_kept_dashboard(kept):
    rows = {r["name"]: r for r in M.entries()}
    assert set(rows) == {"node-aaaa", "node-bbbb", "node-cccc"}
    assert rows["node-bbbb"]["current"] is True
    assert rows["node-aaaa"]["windows"] == 0
    assert rows["node-bbbb"]["windows"] == 2
    assert rows["node-aaaa"]["last_activity"] == pytest.approx(240, abs=30)
    assert rows["node-cccc"]["here"] is False


def test_open_is_relative_in_a_new_tab_and_carries_no_token(kept):
    menu = M.DashboardsMenu()
    menu.toggle()
    assert menu.pop.layout.display != "none"
    markup = " ".join(getattr(w, "value", "") for w in menu.all_widgets())
    assert f'href="{RESUME}?session=node-aaaa"' in markup
    assert 'target="_blank"' in markup
    assert "k" * 43 not in markup
    assert "token" not in markup.lower()
    assert f'href="{DASH}"' in markup and "New dashboard" in markup
    assert "runs on" in markup and "elsewhere" in markup
    assert "this window" in markup
    assert f'?session=node-bbbb"' not in markup, "no Open for this window"


def test_end_asks_once_then_ends_that_dashboard(kept, monkeypatch):
    ended = []
    monkeypatch.setattr(SR, "end_kernel",
                        lambda record, kid, timeout=30.0: ended.append(kid) or True)
    menu = M.DashboardsMenu()
    menu.toggle()
    ends = menu.end_buttons()
    assert sorted(ends) == ["node-aaaa"], "End only for another kept dashboard here"
    ends["node-aaaa"].click()
    assert ended == [], "the first click only asks"
    menu.confirm_button().click()
    assert ended == ["kid-a"]
    assert "node-aaaa" not in {r["session_name"] for r in S.list_records()}


def test_end_reports_a_server_that_cannot_be_asked(kept, monkeypatch):
    monkeypatch.delenv("JUPYTER_TOKEN")
    ok, text = M.end("node-aaaa")
    assert ok is False
    assert "delfin-voila stop" in text
    assert "node-aaaa" in {r["session_name"] for r in S.list_records()}, \
        "nothing is dropped when the kernel cannot be ended"


def test_the_heartbeat_never_recreates_a_dropped_record(kept):
    assert S.refresh_record("node-aaaa") is True
    S.drop_record("node-aaaa", reason="test")
    assert S.refresh_record("node-aaaa") is False
    assert not os.path.exists(S.record_path("node-aaaa"))


def test_the_dashboard_header_shows_the_menu():
    from pathlib import Path

    import delfin.dashboard as D

    source = Path(D.__file__).read_text(encoding="utf-8")
    assert "dashboards_menu.DashboardsMenu(" in source
