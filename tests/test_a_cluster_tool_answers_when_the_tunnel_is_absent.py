"""A cluster tool on a machine without the tunnel: an answer, never a raise.

``cluster_status``, ``cluster_read``, ``cluster_fetch``, ``cluster_upload``,
``cluster_submit``, ``cluster_cancel``, ``cluster_run`` and ``cluster_propose``
are dispatched through ``api_client._DocToolExecutor.execute`` -- the path a
model's call takes. Each one, directly or through ``connection``, imports
``delfin.remote``, whose ``clusters`` module imports the cluster tunnel at
module level (``delfin/remote/clusters.py:47``). Where that package is not
installed the import raised ``ModuleNotFoundError`` out of ``execute``.

The executor's own contract says the opposite: ``_execute_pinned`` re-raises
a dispatch exception unchanged on purpose (``api_client.py:17165-17174``), so
the tool body is the thing that has to answer. The agent loop catches it once
per call and hands the model ``"tool 'cluster_status' failed: No module named
'cluster_tunnel'"`` -- a Python package name where the documented answer is a
refusal code. And the agent is told to fix the machine, not to tell the user
what to install.

``delfin/doctor.py:487`` already settled the same question for the check:
``_clusters_without_the_tunnel`` reports a missing package as missing, never
as "report this as a DELFIN bug". The tool path never got the same treatment.

The absence is forced here through the one seam each path imports behind, so
this measures the same thing on a machine that has the package and on one
that does not -- it must not depend on what this host happens to have
installed.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A

#: Every cluster tool, with the arguments that reach the dispatch of the tool
#: itself -- past the argument check, at the point the remote layer is needed.
_CALLS = [
    ("cluster_status", {}),
    ("cluster_status", {"cluster": "demo"}),
    ("cluster_read", {"cluster": "demo", "path": "run.out"}),
    ("cluster_fetch", {"cluster": "demo", "path": "run.out"}),
    ("cluster_upload", {"cluster": "demo", "src": "input.inp"}),
    ("cluster_submit", {"cluster": "demo", "jobs": [{"folder": "job1"}]}),
    ("cluster_cancel", {"cluster": "demo", "job_ids": ["4242"]}),
    ("cluster_run", {"cluster": "demo", "command": "ls"}),
    ("cluster_propose", {"name": "demo", "host": "demo.example"}),
]


@pytest.fixture(name="absent_tunnel")
def absent_tunnel_fixture(monkeypatch, tmp_path):
    """Make the remote layer see no cluster tunnel, on any machine.

    Forced through ``sys.modules`` before ``delfin.remote`` is imported, and
    re-forced afterwards because ``delfin.remote.clusters`` reads the two
    modules at its own import: an entry equal to None makes ``import
    cluster_tunnel`` raise ImportError even once the real package is in
    memory. The settings name no cluster, so nothing here reads a real one.
    """
    import sys

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    for mod in ("delfin.remote.clusters", "delfin.remote.connection",
                "delfin.remote.setup", "delfin.remote.jobs",
                "delfin.remote.install", "delfin.remote"):
        monkeypatch.delitem(sys.modules, mod, raising=False)
    monkeypatch.setitem(sys.modules, "cluster_tunnel", None)
    monkeypatch.setitem(sys.modules, "cluster_tunnel.api", None)
    monkeypatch.setitem(sys.modules, "cluster_tunnel.config", None)
    monkeypatch.setitem(sys.modules, "cluster_tunnel.budget", None)
    return tmp_path


def _call(tmp_path, tool: str, args: dict) -> dict:
    """One tool call through the executor's dispatch, as a model makes it."""
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    perms = A.KitToolPermissions(workspace=ws)
    perms.mode = "bypassPermissions"
    perms.confirm_callback = lambda *a, **k: True
    perms.task_session_id = "s1"
    executor = A._DocToolExecutor.__new__(A._DocToolExecutor)
    raw = executor.execute(tool, args, perms)
    out = json.loads(raw)
    assert isinstance(out, dict), out
    return out


@pytest.mark.parametrize("tool,args", _CALLS)
def test_the_tool_answers_instead_of_raising(absent_tunnel, tool, args):
    """The call returns a JSON object; it does not raise out of execute."""
    out = _call(absent_tunnel, tool, args)
    assert isinstance(out, dict), out


@pytest.mark.parametrize("tool,args", _CALLS)
def test_the_answer_names_the_missing_package(absent_tunnel, tool, args):
    """A refusal code, and what to install -- not a bare ImportError text."""
    out = _call(absent_tunnel, tool, args)
    assert "error" in out, out
    assert out.get("code") == "tunnel_not_installed", out
    # the name a person can install, not the module the traceback named
    assert "cluster-tunnel" in out["error"], out
    assert "pip install" in out["error"], out


def test_a_call_that_is_refused_on_its_arguments_still_says_so(absent_tunnel):
    """Control: the tunnel answer is not a blanket for every call.

    A tool that returned ``tunnel_not_installed`` for every call would pass
    the tests above while having stopped reading its arguments. ``cluster_run``
    refuses a command longer than 4000 characters in its own body, before it
    needs the remote layer (``cluster_tools.py:571``), so with the tunnel
    absent the answer must still be its ``bad_argument`` -- and it is on the
    baseline too, which is what makes this a control rather than a second
    copy of the fix. (A REQUIRED argument that is missing or blank never
    reaches the tool at all: ``execute`` refuses that before dispatch.)
    """
    too_long = _call(absent_tunnel, "cluster_run",
                     {"cluster": "demo", "command": "x" * 5000})
    assert too_long.get("code") == "bad_argument", too_long
