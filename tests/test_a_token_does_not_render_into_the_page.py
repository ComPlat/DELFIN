"""The address a kept session comes back to is spoken on the terminal.

It is written for the terminal that runs the server, which is where
somebody looks tomorrow after `tmux attach`. When that terminal cannot
be reached, the code fell back to ``print`` — and in a kernel a print is
not a terminal: ipykernel forwards it to the frontend, so the block
rendered IN THE DASHBOARD. Reported from a live session on a cluster:

    [delfin] Session "uc3n991-f3c3" is kept.
             Return:  http://localhost:8866/voila/render/...?token=ZK2d...
             End:     in the dashboard, or Ctrl+C here

That line carries the server token in a URL. The page is the one surface
that does not need it — the status strip already offers the same address
behind a button — and the one where it is worth least, because whoever
reads it is already authenticated. So a kernel that cannot reach the
server's terminal now says nothing at all.

Outside a kernel — a CLI, a test — ``print`` is still exactly right, and
that half must not be lost with the fix.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import session as S


@pytest.fixture()
def no_server_terminal(monkeypatch):
    """The condition that produced the report: a kernel whose parent's
    stdout cannot be opened."""
    monkeypatch.setattr(S, "_server_stdout", lambda: None)
    monkeypatch.setattr(S, "session_name", lambda: "uc3n991-f3c3")
    monkeypatch.setattr(
        S, "resume_url",
        lambda *a, **k: "http://localhost:8866/voila/render/x.ipynb"
                        "?token=SECRET-TOKEN-VALUE&session=uc3n991-f3c3")


def test_a_kernel_says_nothing_it_cannot_say_on_the_terminal(
        no_server_terminal, monkeypatch, capsys):
    monkeypatch.setattr(S, "kernel_id", lambda: "8e7532e0")

    said = S.announce()

    printed = capsys.readouterr()
    assert "SECRET-TOKEN-VALUE" not in printed.out, (
        "the token rendered into the page: a kernel's stdout is the "
        "frontend, not a terminal")
    assert "SECRET-TOKEN-VALUE" not in printed.err
    assert printed.out == "" and printed.err == ""
    # The caller still gets the line; what changed is where it goes.
    assert "is kept" in said and "SECRET-TOKEN-VALUE" in said


def test_outside_a_kernel_it_is_still_printed(no_server_terminal,
                                              monkeypatch, capsys):
    """A CLI or a test has no frontend to leak into, and a caller that
    sees nothing at all would be the other kind of wrong."""
    monkeypatch.setattr(S, "kernel_id", lambda: "")

    S.announce()

    printed = capsys.readouterr().out
    assert "is kept" in printed
    assert "SECRET-TOKEN-VALUE" in printed


def test_a_kernel_with_a_terminal_writes_there(monkeypatch):
    """The normal path is untouched: the line reaches the server's own
    stdout, which is the terminal inside tmux."""
    written: list[str] = []

    class _Sink:
        def write(self, text):
            written.append(text)

        def flush(self):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(S, "kernel_id", lambda: "8e7532e0")
    monkeypatch.setattr(S, "_server_stdout", lambda: _Sink())
    monkeypatch.setattr(S, "session_name", lambda: "uc3n991-f3c3")
    monkeypatch.setattr(S, "resume_url", lambda *a, **k: "http://x/?token=T")

    S.announce()

    assert any("is kept" in w for w in written)
