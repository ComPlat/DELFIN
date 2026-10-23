"""A kept session tells the server's terminal where to come back, also
when the server is protected.

The kernel reached that terminal by opening its parent's stdout through
/proc. The dashboard server is protected (process_guard: non-dumpable, its
/proc belongs to root), so on Linux the open was refused and the
"Return:" line never appeared. The kernel's own stdout is the server's:
ipykernel keeps a copy of what it redirected, and a stream it did not
redirect left fd 1 as it was. Those are written to now.
"""

from __future__ import annotations

import os
import sys

from delfin.dashboard import session as S


class _Redirected:
    """ipykernel's OutStream, as far as this is concerned."""

    def __init__(self, fd: int):
        self._original_stdstream_copy = fd

    def write(self, text):              # what reaches the page
        raise AssertionError("the address must not go to the page")


def test_the_line_goes_to_the_stdout_the_kernel_was_started_with(tmp_path, monkeypatch):
    terminal = tmp_path / "server.log"
    fd = os.open(terminal, os.O_WRONLY | os.O_CREAT | os.O_APPEND)
    try:
        monkeypatch.setattr(sys, "stdout", _Redirected(fd))
        # The parent's /proc is refused, as for a protected server.
        real_open = open

        def no_proc(path, *a, **kw):
            if str(path).startswith("/proc/"):
                raise PermissionError(13, "Permission denied", path)
            return real_open(path, *a, **kw)

        monkeypatch.setattr("builtins.open", no_proc)
        monkeypatch.setattr(S, "kernel_id", lambda: "k" * 8)
        monkeypatch.setattr(S, "resume_url", lambda who: "http://host:1/delfin/resume?session=x")

        S.announce("x")
    finally:
        monkeypatch.undo()
    os.write(fd, b"")                   # ipykernel's copy is still open
    os.close(fd)

    told = terminal.read_text()
    assert "Return:  http://host:1/delfin/resume?session=x" in told


def test_without_a_kernel_stream_the_parent_route_is_still_taken(monkeypatch):
    monkeypatch.setattr(sys, "stdout", object())
    monkeypatch.setattr(sys, "stderr", object())
    opened: list = []

    def fake_open(path, *a, **kw):
        opened.append(str(path))
        raise PermissionError(13, "Permission denied", path)

    monkeypatch.setattr("builtins.open", fake_open)
    assert S._server_stdout() is None
    assert opened == [f"/proc/{os.getppid()}/fd/1", f"/proc/{os.getppid()}/fd/2"]


class _Unwatched:
    """An OutStream that left its descriptor alone (ipykernel under pytest)."""

    _should_watch = False


def test_a_kernel_stream_that_was_never_redirected_writes_to_its_fd(tmp_path, monkeypatch):
    terminal = tmp_path / "server.log"
    fd = os.open(terminal, os.O_WRONLY | os.O_CREAT | os.O_APPEND)
    saved = os.dup(1)
    try:
        os.dup2(fd, 1)
        monkeypatch.setattr(sys, "stdout", _Unwatched())
        out = S._server_stdout()
        with out:
            out.write("Return:  here\n")
        os.write(1, b"")                # fd 1 is still open
    finally:
        os.dup2(saved, 1)
        os.close(saved)
        os.close(fd)
    assert terminal.read_text() == "Return:  here\n"
