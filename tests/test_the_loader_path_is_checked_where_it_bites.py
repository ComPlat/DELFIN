"""A shared library taken from the wrong directory silently removed a capability.

The dynamic loader resolves a shared library once per process and the
first resolution wins for everything after it. On a host where a
directory on LD_LIBRARY_PATH ships an older libstdc++ than the
interpreter's own libraries need, importing pymupdf (the agent's PDF
reader) loads that older one, and every later import needing the newer
C++ ABI is fatal -- here sqlite3, through libicu.

Nothing reported it, and two layers of optional-import fallback turned it
into a capability that was absent rather than an error: stk reaches
sqlite3 through atomlite, and the module that imports stk swallows the
failure and records "not available". A build then runs without it and the
result depends on which import happened first.

It had already been measured once on this installation, in September, and
worked around by dropping LD_LIBRARY_PATH in the gate wrapper. The
workaround is correct for a gate that wants to measure a branch, and it
is also why nobody saw the condition again -- including the module that
was degrading under it. A workaround nothing reports is a fault that
stops being looked for.

So: the doctor asks, in a subprocess, with this process's environment.

Universal by construction. Nothing here names ORCA or any site: it
imports two modules DELFIN itself needs and reports what the loader did.
On a host whose library order is sound both imports succeed and the row
passes, which is what makes it safe to ship.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

from delfin.agent import doctor as D


# ---------------------------------------------------------------------------
# It is asked, and it is asked in the right place
# ---------------------------------------------------------------------------

def test_the_check_is_registered():
    assert any(attr == "_check_loader_path" for _, attr in D._CHECK_ATTRS)


def test_it_asks_a_subprocess_not_this_interpreter():
    """The loader caches LD_LIBRARY_PATH at process start, so an
    in-process answer is already fixed -- and importing the pair here
    would poison the interpreter doing the asking."""
    import inspect

    body = inspect.getsource(D._check_loader_path)
    assert "subprocess.run" in body
    assert "sys.executable" in body
    i_code = body.index("code = (")
    i_run = body.index("subprocess.run")
    assert i_code < i_run


def test_the_pair_is_imported_in_a_fixed_order():
    """The point is the ORDER: the first resolution wins. A check that
    imported them the other way round would pass on a broken host."""
    import inspect

    body = inspect.getsource(D._check_loader_path)
    assert "for first, second, what in _IMPORT_ORDER_PAIRS" in body
    first, second, _what = D._IMPORT_ORDER_PAIRS[0]
    assert (first, second) == ("pymupdf", "sqlite3")


def test_each_pair_names_the_subsystem_it_belongs_to():
    """A row saying "an import failed" sends the reader into Python, and
    the fault is three layers below it."""
    for first, second, what in D._IMPORT_ORDER_PAIRS:
        assert first and second
        assert len(what) > 10


# ---------------------------------------------------------------------------
# What it reports, driven both ways
# ---------------------------------------------------------------------------

def _rows(monkeypatch, *, returncode, stderr=""):
    def fake_run(argv, **kw):
        return subprocess.CompletedProcess(argv, returncode, "", stderr)

    monkeypatch.setattr(D.subprocess, "run", fake_run)
    return D._check_loader_path({})


def test_a_sound_host_passes(monkeypatch):
    rows = _rows(monkeypatch, returncode=0)
    assert rows
    assert all(r["status"] == D.PASS for r in rows)


def test_a_broken_order_warns_and_names_the_library(monkeypatch):
    stderr = (
        "Traceback (most recent call last):\n"
        "ImportError: /opt/example/lib/libstdc++.so.6: version "
        "`CXXABI_1.3.15' not found (required by /env/lib/libicui18n.so.78)\n")
    rows = _rows(monkeypatch, returncode=1, stderr=stderr)
    row = rows[0]
    assert row["status"] != D.PASS
    assert "CXXABI_1.3.15" in row["detail"]
    assert "/opt/example/lib/libstdc++.so.6" in row["fix"], (
        "the fix must name the file the loader took, not the exception")
    assert "ImportError" not in row["fix"], (
        "splitting the message on ':' took the exception name as the library")
    assert "LD_LIBRARY_PATH" in row["fix"]


def test_a_failure_with_no_library_named_still_reports(monkeypatch):
    rows = _rows(monkeypatch, returncode=1,
                 stderr="ModuleNotFoundError: No module named 'pymupdf'\n")
    row = rows[0]
    assert row["status"] != D.PASS
    assert row["detail"]
    assert row["fix"]


def test_an_empty_stderr_does_not_produce_an_empty_row(monkeypatch):
    rows = _rows(monkeypatch, returncode=1, stderr="")
    assert rows[0]["detail"].strip()


def test_an_interpreter_that_cannot_be_asked_is_reported_not_assumed(
        monkeypatch):
    def boom(*a, **k):
        raise OSError("no exec")

    monkeypatch.setattr(D.subprocess, "run", boom)
    rows = D._check_loader_path({})
    assert rows[0]["status"] != D.PASS
    assert "could not ask" in rows[0]["detail"]
    assert rows[0]["fix"]


def test_a_timeout_does_not_escape(monkeypatch):
    def slow(*a, **k):
        raise subprocess.TimeoutExpired(cmd="python", timeout=120)

    monkeypatch.setattr(D.subprocess, "run", slow)
    rows = D._check_loader_path({})
    assert rows[0]["status"] != D.PASS


def test_every_non_passing_row_carries_a_remedy(monkeypatch):
    """The module's convention, and the reason this check exists: the
    condition is three layers below Python and nobody finds it by
    reading a status."""
    for kwargs in ({"returncode": 1, "stderr": "ImportError: boom\n"},
                   {"returncode": 2, "stderr": ""}):
        for row in _rows(monkeypatch, **kwargs):
            if row["status"] != D.PASS:
                assert row["fix"], row


def test_the_doctor_never_raises_with_this_check(monkeypatch, tmp_path):
    """The module's stated contract."""
    monkeypatch.setattr(D.subprocess, "run", lambda *a, **k: (_ for _ in ())
                        .throw(OSError("nope")))
    assert D._check_loader_path({})
