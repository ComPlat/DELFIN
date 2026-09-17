"""A session kept on another login node is left alone.

Kept sessions are announced in the home directory, which every login node
shares; the kernel belongs to the server on the machine that wrote the
record. On 2026-09-16 a login landed on uc3n990 while three agent sessions
were kept on uc3n991. The button asked uc3n990's server for the kernel; it
did not run it, took the record for stale and removed it -- and uc3n991's
server, finding its unwatched kernel no longer kept, ended it ten seconds
later, with the three agents in it.
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

import pytest

from delfin.dashboard import resume_server as R
from delfin.dashboard import session as S


@pytest.fixture(autouse=True)
def clean(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "kept"))
    S._reset_for_tests()
    yield
    S._reset_for_tests()


def _kept_on(host: str, name: str = "uc3n991-e2d3", kid: str = "kernel-bebd") -> Path:
    path = Path(S.write_record(name, kid=kid))
    record = json.loads(path.read_text())
    record["host"] = host
    path.write_text(json.dumps(record))
    return path


class _FakeManager:
    def __init__(self, known=()):
        self.known = set(known)
        self.started = []

    def __contains__(self, kid):
        return kid in self.known

    async def start_kernel(self, *args, **kwargs):
        self.started.append(kwargs)
        return "fresh-kernel"


def _return_to(name: str, known=()):
    manager = R.resume_kernel_manager_class(_FakeManager)(known=known)
    got = asyncio.run(manager.start_kernel(
        kernel_name="python3",
        env={"VOILA_REQUEST_URL": f"http://h:1/voila/render/r.ipynb?session={name}"}))
    return got, manager


def test_returning_from_another_node_keeps_the_record():
    path = _kept_on("some-other-node")
    got, manager = _return_to("uc3n991-e2d3", known=[])
    assert got == "fresh-kernel"
    assert path.exists(), "the record of a session on another node was removed"
    assert R.kernel_for_session("uc3n991-e2d3") == "kernel-bebd"
    assert R.kept_kernel_ids() == {"kernel-bebd"}, (
        "its own server would no longer see it as kept, and end it")


def test_a_record_from_this_node_whose_kernel_is_gone_is_still_removed():
    path = _kept_on(S._hostname())
    got, _manager = _return_to("uc3n991-e2d3", known=[])
    assert got == "fresh-kernel"
    assert not path.exists()


def test_the_page_names_the_node_the_session_runs_on(monkeypatch):
    _kept_on("uc3n991")
    monkeypatch.setattr(S, "kernel_id", lambda: "fresh-kernel")
    text = S.why_not_resumed(
        request_url="http://h:1/voila/render/r.ipynb?session=uc3n991-e2d3")
    assert "runs on uc3n991" in text and "log in to uc3n991" in text
    assert "delfin-agent stop-all" in text


def test_the_landing_banner_does_not_offer_a_return_it_cannot_make():
    other = _kept_on("uc3n991")
    record = json.loads(other.read_text())
    html = S._banner_html([record])
    assert "Open this session" not in html
    assert "runs on <code>uc3n991</code>" in html


def test_the_landing_banner_still_offers_a_session_on_this_node(monkeypatch):
    monkeypatch.setenv(S.RESUME_PATH_ENV, "/voila/render/delfin_voila_runtime/delfin_resume.ipynb")
    here = json.loads(_kept_on(S._hostname(), name="here-1", kid="k1").read_text())
    here["request_url"] = "http://localhost:8866/voila/render/dash.ipynb"
    assert "Open this session" in S._banner_html([here])
