"""A cluster that needs a login is a NEEDS-YOU item that a login closes, never an answer.

Input: the attention inbox (the suite redirects it), ``delfin.remote.attention``
and ``clusters_feed.login`` with the connection replaced. Output: the inbox
records, the window's item (``attention_feed.item_of``) and what
``attention.answer_item`` and ``drain_resolved`` return.

Why: what a person types for a login is a password. An attention answer
goes to the next turn of a model, so a login item must never take one; the
window opens its masked field for the item's cluster instead, and a login
that opened the tunnel closes the item with nothing to hand over.
"""

from __future__ import annotations

from delfin.agent import attention
from delfin.agent.service import attention_feed, clusters_feed


def _pending(cluster: str) -> list[dict]:
    return [e for e in attention.list_pending("cluster_login") if e.get("cluster") == cluster]


def test_one_pending_item_per_cluster(tmp_path):
    from delfin.remote import attention as remote_attention
    first = remote_attention.ask_for_login("kcist", reason="cluster_status needs kcist",
                                           workspace=tmp_path)
    again = remote_attention.ask_for_login("kcist")
    other = remote_attention.ask_for_login("horeka")
    assert first and again == first and other != first
    assert len(_pending("kcist")) == 1


def test_the_window_gets_a_login_item_naming_its_cluster(tmp_path):
    from delfin.remote import attention as remote_attention
    remote_attention.ask_for_login("kcist", reason="cluster_jobs needs a login to kcist")
    (event,) = _pending("kcist")
    item = attention_feed.item_of(event, findings={}, presence={},
                                  titles=attention_feed.Titles())
    assert (item["kind"], item["group"], item["cluster"]) == ("login", "needs", "kcist")
    assert item["name"] == "log in to kcist"
    assert item["text"] == "cluster_jobs needs a login to kcist"


def test_a_text_answer_is_refused_and_the_item_stays_open():
    from delfin.remote import attention as remote_attention
    item_id = remote_attention.ask_for_login("kcist")
    got = attention.answer_item(item_id, "hunter2")
    assert got["ok"] is False and "login field" in got["error"]
    assert _pending("kcist") and _pending("kcist")[0]["answer"] is None
    assert "hunter2" not in repr(attention.list_pending())


def test_a_login_closes_the_item_with_nothing_for_a_model(monkeypatch):
    from delfin.remote import attention as remote_attention
    from delfin.remote import connection

    class Opened:
        def login(self, password, otp, timeout=120):
            return True

    monkeypatch.setattr(connection, "connect", lambda name, settings=None: Opened())
    remote_attention.ask_for_login("kcist")
    assert clusters_feed.login("kcist", "pw", "")["ok"] is True
    assert _pending("kcist") == []
    assert attention.drain_resolved(kinds={"cluster_login"}) == []
