"""The Domain State a session's turn carries describes that session's workspace.

Input: a dashboard session (Code mode) whose workspace is a git worktree,
launched from a checkout that holds a task file the worktree does not
have. Output: the Domain State block and the session primer on the wire of
the first turn name the worktree's files and branch, never the launch
checkout's.

The defect this pins: the block listed ``ctx.agent_dir`` -- a directory
shared by every session of the dashboard -- under the label ``workspace``,
and the primer read branch and commits from ``ctx.repo_dir``. A session in
a worktree was told about a task file that did not exist where its tools
resolve paths, and looked for it. The workspace the tools resolve against
is the engine's permission workspace (``kit_permissions.workspace``), so
that is the one the block reads.

Driven through the real tab: create_tab, Send, and the messages the stub
transport received.
"""
from __future__ import annotations

import json
import subprocess
import time

import pytest

from tests.conftest import child_env

_SID = "d0d1d2d3d4d5d6d7d8d9dadbdcdddedf"


def _git(env, cwd, *args):
    subprocess.run(["git", *args], cwd=str(cwd), env=env, check=True,
                   capture_output=True, text=True, timeout=30)


def _widgets(tab):
    import ipywidgets as widgets
    found, stack = {}, [tab]
    while stack:
        node = stack.pop()
        if isinstance(node, (widgets.Dropdown, widgets.Button,
                             widgets.Textarea)):
            found.setdefault((type(node).__name__, node.description), node)
        stack.extend(getattr(node, "children", ()) or ())
    return found


def _wire_text(messages):
    out = []
    for m in messages or []:
        c = m.get("content")
        if isinstance(c, list):
            c = " ".join(str(p.get("text", "")) for p in c
                         if isinstance(p, dict))
        out.append(str(c or ""))
    return "\n".join(out)


def test_the_domain_state_names_the_worktree_not_the_launch_checkout(
        tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.agent import api_client as A
    from delfin.agent import scheduler as S
    from delfin import user_settings
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    env = child_env(tmp_path)
    launch = tmp_path / "launch"
    launch.mkdir()
    _git(env, launch, "init", "-q", "-b", "launch-branch")
    (launch / "TASK_FROM_LAUNCH.md").write_text("the launch checkout's task\n")
    (launch / "shared.txt").write_text("both\n")
    _git(env, launch, "add", "shared.txt")
    _git(env, launch, "-c", "user.name=t", "-c", "user.email=t@example.invalid",
         "commit", "-q", "-m", "launch commit")
    wt = tmp_path / "wt"
    _git(env, launch, "worktree", "add", "-q", "-b", "session-branch", str(wt))
    (wt / "ONLY_IN_WORKTREE.py").write_text("x = 1\n")

    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-test")
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "kit", "model": "kit.deepseek-v4-flash",
        "permission_profile": "all_free", "probe_endpoint": False}}))

    sent = []

    def _answer(self, system="", messages=None, **kw):
        sent.append(list(messages or []))
        return
        yield

    monkeypatch.setattr(A.OpenAIClient, "stream_message", _answer)

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office",
                           agent_dir=launch, repo_dir=launch)
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(wt)
    ctx.initial_session_id = _SID
    ctx.sessions_managed = True
    ctx.presence_key = "domain-state"
    tab, refs = tab_agent.create_tab(ctx)
    w = _widgets(tab)
    try:
        mode = next(d for (kind, _), d in w.items() if kind == "Dropdown"
                    and "solo" in [v for _, v in d.options])
        mode.value = "solo"
        w[("Textarea", "")].value = "Where am I?"
        w[("Button", "Send")].click()
        deadline = time.monotonic() + 30
        while not sent and time.monotonic() < deadline:
            time.sleep(0.05)
        assert sent, "the turn never reached the transport"
        engine = refs["state"]["engine"]
        assert str(engine.kit_permissions.workspace) == str(wt.resolve()) or \
            str(engine.kit_permissions.workspace) == str(wt)
        text = _wire_text(sent[0])
        assert "--- Domain State ---" in text
        assert "ONLY_IN_WORKTREE.py" in text, text[-2000:]
        assert "TASK_FROM_LAUNCH.md" not in text, text[-2000:]
        assert "launch-branch" not in text
        assert "session-branch" in text
    finally:
        refs["shutdown"](save=False)
