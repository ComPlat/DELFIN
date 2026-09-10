"""Answering the address a kept session left behind.

The dashboard's own URL cannot do this. It goes through Voila's
renderer, which EXECUTES the notebook — so opening it against a live
kernel would run all nineteen tabs a second time in a kernel that
already has them, which is the opposite of coming back.

So a resume takes a different route and a different notebook. The
notebook has one cell:

    from delfin.dashboard import session; session.resume()

Voila renders that the way it renders anything, and the one thing that
changes is WHICH kernel it renders into: a kernel manager that,
recognising a resume request, hands back the kernel the session is
already running in instead of starting a fresh one. The cell then
re-displays the widget objects that never left, and the page is the
dashboard again — every value, every callback, the agent's thread still
going.

Nothing is serialised and nothing is restored. That is why this is four
small pieces rather than a persistence layer.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Optional
from urllib.parse import parse_qs, quote, urlsplit

from delfin.dashboard import session as _session


#: The query key that marks a request as a resume, and names which
#: session to resume. It rides on the URL because that is the only thing
#: Voila hands the kernel manager -- the request URL reaches it through
#: the kernel environment it builds.
SESSION_QUERY_KEY = "session"

RESUME_NOTEBOOK_NAME = "delfin_resume.ipynb"


# ---------------------------------------------------------------------------
# The notebook a resume renders
# ---------------------------------------------------------------------------

def resume_notebook_source() -> dict:
    """One cell. Anything more would run in somebody's live kernel."""
    return {
        "cells": [{
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "# Re-display the dashboard this kernel already holds.\n",
                "# Nothing is rebuilt: the widget objects, their values and\n",
                "# their callbacks never went away.\n",
                "from delfin.dashboard import session as _s\n",
                "if not _s.resume():\n",
                "    print('Diese Sitzung ist nicht mehr da.')\n",
            ],
        }],
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3", "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def stage_resume_notebook(root_dir: str | os.PathLike) -> str:
    """Write the resume notebook beside the dashboard, and return it.

    Rewritten every start rather than only when missing: it is generated,
    it is one cell, and a stale copy from an older DELFIN would run
    against a kernel from a newer one.
    """
    directory = Path(root_dir).resolve() / "delfin_voila_runtime"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / RESUME_NOTEBOOK_NAME
    path.write_text(json.dumps(resume_notebook_source(), indent=1),
                    encoding="utf-8")
    return str(path)


def resume_render_url(root_dir: str | os.PathLike, notebook: str,
                      name: str) -> str:
    """The Voila render URL for the resume notebook, naming the session."""
    rel = Path(notebook).resolve().relative_to(
        Path(root_dir).resolve()).as_posix()
    return f"/voila/render/{quote(rel)}?{SESSION_QUERY_KEY}={quote(name)}"


# ---------------------------------------------------------------------------
# Recognising a resume, from the only thing the kernel manager is given
# ---------------------------------------------------------------------------

def requested_session(request_url: str) -> str:
    """The session a request asks to resume, or "" for an ordinary one."""
    if not request_url:
        return ""
    try:
        query = urlsplit(request_url).query
    except ValueError:
        return ""
    values = parse_qs(query).get(SESSION_QUERY_KEY) or []
    return (values[0] or "").strip() if values else ""


def kernel_for_session(name: str, *, root: str = "") -> str:
    """The kernel id a session announced, or "" if it announced none."""
    if not name:
        return ""
    for record in _session.list_records(root=root):
        if record.get("session_name") == name:
            return str(record.get("kernel_id") or "")
    return ""


def resume_kernel_manager_class(base: type) -> type:
    """A kernel manager that reuses a kept session's kernel.

    Wrapping whatever class is configured rather than a fixed one: Voila
    already builds its own manager on top of the server's, and replacing
    that would take the preheat and pooling behaviour with it.

    The reuse is refused unless the named kernel is still known to this
    manager. A record can outlive its kernel -- a crash, a machine
    restart -- and handing Voila a dead id would render the resume
    notebook into nothing at all.
    """

    class ResumeAwareKernelManager(base):                # type: ignore[misc]

        async def _delfin_existing(self, env: Optional[dict]) -> str:
            name = requested_session((env or {}).get("VOILA_REQUEST_URL", ""))
            if not name:
                return ""
            kid = kernel_for_session(name)
            if not kid:
                return ""
            try:
                known = kid in self
            except Exception:
                known = False
            if not known:
                # The session announced a kernel that is gone. Drop the
                # record here rather than leaving it to mislead the next
                # visitor.
                _session.drop_record(name)
                return ""
            return kid

        async def start_kernel(self, *args: Any, **kwargs: Any):
            existing = await self._delfin_existing(kwargs.get("env"))
            if existing:
                return existing
            return await super().start_kernel(*args, **kwargs)

    ResumeAwareKernelManager.__name__ = f"ResumeAware{base.__name__}"
    return ResumeAwareKernelManager
