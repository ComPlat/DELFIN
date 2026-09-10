"""The dashboard's own URL cannot bring you back, and that is the point.

It goes through Voila's renderer, which EXECUTES the notebook. Opening
it against a live kernel would run all nineteen tabs a second time in a
kernel that already has them — the opposite of coming back.

So a resume renders a different notebook, one cell long, into the kernel
the session is ALREADY running in. That cell re-displays the widget
objects that never left. Nothing is serialised and nothing is restored,
which is why this is four small pieces instead of a persistence layer.

The half that needs the most care is refusal: a record can outlive its
kernel — a crash, a machine restart — and handing Voila a dead id would
render the resume into nothing at all.
"""

from __future__ import annotations

import asyncio
import json
import tempfile
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


# ---------------------------------------------------------------------------
# The notebook a resume renders
# ---------------------------------------------------------------------------

def test_the_resume_notebook_has_exactly_one_cell():
    """It runs in somebody's live kernel. Anything more than the
    re-display would be executed against state that is already there."""
    nb = R.resume_notebook_source()
    assert len(nb["cells"]) == 1
    assert nb["cells"][0]["cell_type"] == "code"


def test_the_cell_only_re_displays():
    src = "".join(R.resume_notebook_source()["cells"][0]["source"])
    assert "session" in src and "resume()" in src
    assert "create_dashboard" not in src, (
        "building a second dashboard in a kernel that has one would "
        "double every tab")


def test_it_says_so_when_the_session_is_gone():
    src = "".join(R.resume_notebook_source()["cells"][0]["source"])
    assert "if not _s.resume():" in src


def test_staging_overwrites_a_stale_copy(tmp_path):
    """It is generated and one cell long; a copy from an older DELFIN
    would run against a kernel from a newer one."""
    first = Path(R.stage_resume_notebook(tmp_path))
    first.write_text('{"cells": [], "nbformat": 4, "nbformat_minor": 5}')
    second = Path(R.stage_resume_notebook(tmp_path))
    assert first == second
    assert len(json.loads(second.read_text())["cells"]) == 1


def test_the_staged_notebook_is_valid_json_nbformat(tmp_path):
    data = json.loads(Path(R.stage_resume_notebook(tmp_path)).read_text())
    assert data["nbformat"] == 4
    assert data["metadata"]["kernelspec"]["name"] == "python3"


# ---------------------------------------------------------------------------
# Recognising a resume
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("url,want", [
    ("http://h:1/voila/render/x.ipynb?session=abc&token=t", "abc"),
    ("http://h:1/voila/render/x.ipynb?token=t&session=abc", "abc"),
    ("http://h:1/voila/render/x.ipynb?token=t", ""),
    ("http://h:1/voila/render/x.ipynb", ""),
    ("", ""),
    ("not a url", ""),
    ("http://h:1/x?session=", ""),
    ("http://h:1/x?session=%20%20", ""),
])
def test_a_request_names_its_session_or_it_does_not(url, want):
    assert R.requested_session(url) == want


def test_the_render_url_carries_the_name(tmp_path):
    nb = R.stage_resume_notebook(tmp_path)
    url = R.resume_render_url(tmp_path, nb, "uc3n990-ab12")
    assert url.startswith("/voila/render/")
    assert f"{R.SESSION_QUERY_KEY}=uc3n990-ab12" in url


def test_a_name_that_needs_quoting_survives(tmp_path):
    nb = R.stage_resume_notebook(tmp_path)
    url = R.resume_render_url(tmp_path, nb, "host with space")
    assert R.requested_session("http://h:1" + url) == "host with space"


# ---------------------------------------------------------------------------
# Finding the kernel — and refusing when it is gone
# ---------------------------------------------------------------------------

def test_a_kept_session_is_found_by_name():
    S.write_record("probe-a", kid="kernel-aaa")
    assert R.kernel_for_session("probe-a") == "kernel-aaa"


def test_an_unknown_name_finds_nothing():
    S.write_record("probe-a", kid="kernel-aaa")
    assert R.kernel_for_session("probe-b") == ""
    assert R.kernel_for_session("") == ""


class _FakeManager:
    """Stands in for the configured kernel manager."""

    def __init__(self, known=()):
        self.known = set(known)
        self.started = []

    def __contains__(self, kid):
        return kid in self.known

    async def start_kernel(self, *args, **kwargs):
        self.started.append(kwargs)
        return "fresh-kernel"


def _managed(known=()):
    cls = R.resume_kernel_manager_class(_FakeManager)
    return cls(known=known)


def _start(manager, url):
    return asyncio.run(manager.start_kernel(
        kernel_name="python3", env={"VOILA_REQUEST_URL": url}))


def test_a_resume_reuses_the_kernel_the_session_runs_in():
    S.write_record("probe-a", kid="kernel-aaa")
    mgr = _managed(known=["kernel-aaa"])
    got = _start(mgr, "http://h:1/voila/render/r.ipynb?session=probe-a")
    assert got == "kernel-aaa"
    assert mgr.started == [], "a resume must not start a kernel"


def test_an_ordinary_request_still_starts_a_kernel():
    """The default path must be untouched — this class wraps whatever is
    configured and has to stay invisible to everything else."""
    mgr = _managed(known=["kernel-aaa"])
    got = _start(mgr, "http://h:1/voila/render/dash.ipynb?token=t")
    assert got == "fresh-kernel"
    assert len(mgr.started) == 1


def test_a_record_that_outlived_its_kernel_is_refused_and_removed():
    """A crash or a machine restart leaves the record behind. Handing
    Voila a dead id would render the resume into nothing at all, and
    leaving the record would mislead the next visitor too."""
    S.write_record("probe-a", kid="kernel-gone")
    mgr = _managed(known=[])                 # the kernel is not there
    got = _start(mgr, "http://h:1/voila/render/r.ipynb?session=probe-a")
    assert got == "fresh-kernel"
    assert R.kernel_for_session("probe-a") == "", (
        "the stale record was left to mislead the next visitor")


def test_a_resume_for_an_unknown_session_falls_through():
    mgr = _managed(known=["kernel-aaa"])
    got = _start(mgr, "http://h:1/voila/render/r.ipynb?session=never-existed")
    assert got == "fresh-kernel"


def test_the_wrapper_keeps_the_wrapped_class_recognisable():
    """Voila builds its own manager on top of the server's; a wrapper
    that hid which one it was would make a misconfiguration unreadable."""
    cls = R.resume_kernel_manager_class(_FakeManager)
    assert "_FakeManager" in cls.__name__
    assert issubclass(cls, _FakeManager)
