"""The CLI with its network edge replaced by recorders.

Run as ``python -P child.py <cli arguments>`` with ``PYTHONPATH`` set to the
tree under test. Input (environment): ``INJECTION_SUITE_FIXTURE``, a JSON
file ``{"pages": {url: text}, "search": {query-word: [results]}}``, and
``INJECTION_SUITE_CANARY``, a JSONL file this process appends one record to
for every call that reached an edge. Output: whatever ``delfin.agent.cli``
prints; the exit code is the CLI's.

What is replaced, and only that: the function that downloads a page
(``web_tools.web_fetch``), the one that runs a web search
(``web_tools.web_search``), the desktop/phone notification
(``notify.send_notification``, ``messenger.activation.post_note``; the
harness's own alerts go through it too, recorded as ``desktop``), the
``push_notification`` and ``session_message`` executors (recorded under
their names, so only the tool's own calls count). Everything in front of them -- the
permission gate, the classifier, the fences, the dispatcher -- is the
tree's own code, so a record means the harness let the call through.

Rejected: patching the gate or the dispatcher (the thing measured); a real
network target (nothing may leave this machine during a measurement).
"""

from __future__ import annotations

import json
import os
import sys
import time
from urllib.parse import urlsplit

_FIXTURE = os.environ.get("INJECTION_SUITE_FIXTURE", "")
_CANARY = os.environ.get("INJECTION_SUITE_CANARY", "")


def _record(kind: str, **fields) -> None:
    if not _CANARY:
        return
    with open(_CANARY, "a", encoding="utf-8") as fh:
        fh.write(json.dumps({"edge": kind, "ts": time.time(), **fields}) + "\n")


def _fixture() -> dict:
    try:
        with open(_FIXTURE, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def install() -> None:
    fixture = _fixture()
    pages = fixture.get("pages") or {}
    search = fixture.get("search") or {}

    from delfin.agent import web_tools

    def web_fetch(url: str, *, timeout_s: int = 15) -> dict:
        _record("web_fetch", url=url, host=urlsplit(url).hostname or "")
        text = pages.get(url)
        if text is None:
            text = pages.get(url.rstrip("/"))
        if text is None:
            host = urlsplit(url).hostname or ""
            if "attacker" in host:
                return {"url": url, "status": "ok", "content_type":
                        "text/plain", "text": "ok", "truncated": False,
                        "bytes": 2}
            return {"error": "HTTP 404: Not Found", "url": url}
        return {"url": url, "status": "ok", "content_type": "text/html",
                "text": text, "truncated": False, "bytes": len(text)}

    def web_search(query: str, *, max_results: int = 8,
                   timeout_s: int = 15) -> dict:
        _record("web_search", query=query)
        hits = []
        for word, results in search.items():
            if word.lower() in str(query).lower():
                hits.extend(results)
        return {"query": query, "results": hits[:max_results]}

    web_tools.web_fetch = web_fetch
    web_tools.web_search = web_search

    from delfin.agent import notify

    def send_notification(title: str, body: str = "", *,
                          urgency: str = "normal") -> bool:
        _record("desktop", title=title, body=body)
        return True

    notify.send_notification = send_notification
    try:
        from delfin.messenger import activation

        def post_note(sid, title, body):
            _record("phone", title=title, body=body)
            return True

        activation.post_note = post_note
    except Exception:
        pass

    from delfin.agent import api_client

    def _session_message(self, arguments: dict, perms=None) -> str:
        to = str(arguments.get("to") or "").strip()
        if not to:
            return json.dumps({"sessions": [{"key": "ops-bot",
                                             "title": "ops"}]})
        _record("session_message", to=to,
                message=str(arguments.get("message") or "")[:500])
        return json.dumps({"status": "queued", "to": to,
                           "message_id": "m-1"})

    api_client._DocToolExecutor._execute_session_message = _session_message

    def _push_notification(self, arguments: dict, perms=None) -> str:
        _record("push_notification", title=str(arguments.get("title", "")),
                body=str(arguments.get("body", ""))[:500])
        return json.dumps({"status": "ok", "sent": True, "phone": True})

    api_client._DocToolExecutor._execute_push_notification = (
        _push_notification)
    import delfin
    _record("loaded", path=os.path.dirname(delfin.__file__))


def main() -> int:
    install()
    from delfin.agent.cli import main as cli_main
    return int(cli_main(sys.argv[1:]) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
