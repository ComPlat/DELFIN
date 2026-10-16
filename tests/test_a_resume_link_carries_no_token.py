"""A link back into a kept session does not carry the server token.

The browser that opened the dashboard holds the login cookie, so the link
works without the token there; any other browser lands on the login page,
which takes the token and then continues to the same link (``next=``).
A token inside a URL ends up in browser history, terminal scrollback that
is copied into tickets, and proxy logs; ``delfin-agent where --token``
prints it on its own when the login page asks for it.
"""

from __future__ import annotations

import json

from delfin.agent import where
from delfin.dashboard import session as S

_RESUME_PATH = "/voila/render/delfin_voila_runtime/delfin_resume.ipynb"


def test_the_session_resume_url_drops_the_token():
    url = S.resume_url(
        "login01-ab12",
        request_url=("http://login01:8866/voila/render/"
                     "delfin_voila_runtime/delfin_dashboard.ipynb?token=abc"),
        resume_path=_RESUME_PATH)
    assert url == "http://login01:8866" + _RESUME_PATH + "?session=login01-ab12"


def test_where_prints_resume_links_without_the_token(monkeypatch):
    monkeypatch.setattr(where, "kept_session_names", lambda: ["kept-1"])
    record = {"port": 8866, "token": "SECRET-TOKEN-7f3a",
              "resume_path": _RESUME_PATH}
    links = where.resume_links(record, [])
    assert links == [("kept-1",
                      "http://localhost:8866" + _RESUME_PATH + "?session=kept-1")]
    text = where.format_text(dict(record, host="h", running=True), [])
    assert "SECRET-TOKEN-7f3a" not in text
    assert "where --token" in text


def test_where_token_prints_only_the_token(tmp_path, monkeypatch, capsys):
    note = tmp_path / "dashboard_here.json"
    note.write_text(json.dumps({"port": 8866, "token": "SECRET-TOKEN-7f3a"}),
                    encoding="utf-8")
    monkeypatch.setattr(where, "record_path", lambda: note)
    assert where.main(["--token"]) == 0
    assert capsys.readouterr().out.strip() == "SECRET-TOKEN-7f3a"
