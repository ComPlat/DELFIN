"""End-to-end behavioural test of the KIT-Toolbox coding agent.

Simulates a typical workflow for any user who wants the agent to
build software in their own project: starting in DELFIN's repo,
granting access to a separate generic project via the
``project_dev`` bundle, reading + editing source there, creating a
venv, installing dependencies, running the script + pytest, observing
a bug, fixing the source, re-running pytest until green. Along the
way also editing a Jupyter analysis notebook, tracking progress via
the persistent task list, and verifying the audit log captured every
code-modifying action with the correct decision label.

The point is to confirm the agent BEHAVES correctly — not just that
each tool returns without crashing. Every step asserts something
about the expected outcome:

* writes that target Self-Modification-Guard files are blocked even
  with the workspace and a confirm-callback in place
* the auto-test hint after an edit names the *correct* test file
* background bash output is readable while still running and after
  exit, with smart head+tail truncation when output is long
* persistent state (tasks + bundle settings) survives a simulated
  "restart" (drop in-memory caches, re-instantiate)
* the audit log lists every recorded action with the right
  ``decision`` (``ok`` for successful work, ``denied`` for blocked
  attempts) and never lets a denied action slip through as ``ok``
* destructive bash patterns and secret reads stay blocked across
  all modes — including ``bypassPermissions``

If any assertion here fails, the agent has regressed in a way that
would visibly hurt the user's day-to-day project work.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from delfin.agent import audit_log as al
from delfin.agent import agent_tasks
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def workspace(tmp_path) -> Path:
    """Simulates the DELFIN repo where the dashboard runs."""
    ws = tmp_path / "delfin"
    ws.mkdir()
    return ws


@pytest.fixture
def user_project(tmp_path) -> Path:
    """Simulates a generic user project granted to the agent via
    extra_workspace_dir / remember_permission_bundle. The contents
    are intentionally a small Python package with one module + one
    test so the build/run/test loop has something real to act on."""
    p = tmp_path / "user_project"
    p.mkdir()
    (p / "mylib").mkdir()
    (p / "mylib" / "core.py").write_text(
        "def compute(n_trials=100):\n"
        "    return [i * 0.5 for i in range(n_trials)]\n"
    )
    (p / "tests").mkdir()
    (p / "tests" / "test_core.py").write_text(
        "from mylib.core import compute\n"
        "def test_compute():\n"
        "    assert len(compute(5)) == 5\n"
    )
    return p


@pytest.fixture
def audit_to_tmp(tmp_path, monkeypatch) -> Path:
    """Redirect ~/.delfin/audit.log to a tmp file for the test."""
    log = tmp_path / "audit.log"
    monkeypatch.setattr(al, "_default_log_path", lambda: log)
    return log


@pytest.fixture(autouse=True)
def reset_task_store_cache():
    """Drop in-memory store cache so each test starts fresh."""
    agent_tasks._STORES.clear()
    yield
    agent_tasks._STORES.clear()


# ---------------------------------------------------------------------------
# The actual workflow simulation
# ---------------------------------------------------------------------------


def test_full_agent_session(
    workspace, user_project, audit_to_tmp, tmp_path, monkeypatch,
):
    # ---- 0. Engine boots in DELFIN's repo. Initial perms have only
    #         DELFIN as workspace; user-project is unknown to the agent. ---
    seen_confirms: list[dict] = []

    def confirm_cb(name: str, args: dict, preview: str) -> bool:
        # User clicks "Erlauben" for any genuine bundle / extra-dir
        # confirm; this matches the dashboard's mode chip flow.
        seen_confirms.append({"name": name, "args": dict(args)})
        return True

    perms = KitToolPermissions(
        workspace=workspace,
        mode="default",
        confirm_callback=confirm_cb,
    )
    # Redirect KIT-settings persistence to tmp so the test never
    # touches the user's real ~/.delfin.
    from delfin.agent import kit_settings as ks
    monkeypatch.setattr(ks, "USER_SETTINGS_PATH", tmp_path / "user_kit.json")

    # ---- 1. the user says "arbeite in /home/user/user-project". The agent
    #         calls remember_permission_bundle to set up the project. -
    bundle_out = json.loads(_doc_executor.execute(
        "remember_permission_bundle",
        {"profile": "project_dev",
         "directory": str(user_project),
         "scope": "repo",
         "rationale": "Bayesian-opt project"},
        perms))
    assert bundle_out["status"] == "persisted", bundle_out
    assert seen_confirms[-1]["name"] == "remember_permission_bundle"
    # user-project is now a workspace root and venv patterns are live.
    roots = [Path(p) for p in perms.all_workspace_roots()]
    assert user_project.resolve() in roots
    assert any("pip" in p for p in perms.bash_auto_allow_patterns)

    # ---- 2. Read + edit the optimizer source. -------------------------
    src = user_project / "mylib" / "core.py"
    read_out = _doc_executor.execute(
        "read_file", {"path": str(src)}, perms)
    # read_file returns formatted text with line-number prefix, not JSON.
    assert "n_trials=100" in read_out

    edit_out = _doc_executor.execute("edit_file", {
        "path": str(src),
        "old_string": "n_trials=100",
        "new_string": "n_trials=200",
    }, perms)
    assert "Edited" in edit_out
    # The auto-test hint must point at the EXISTING test module.
    assert "tests/test_core.py" in edit_out, (
        "auto-test hint missing or wrong target:\n" + edit_out
    )
    # And the file actually changed.
    assert "n_trials=200" in src.read_text()

    # ---- 3. Self-Mod-Guard: editing api_client.py STILL refused even
    #         though we're in a fully-confirmed flow. -------------------
    delfin_self = workspace / "delfin" / "agent" / "api_client.py"
    delfin_self.parent.mkdir(parents=True)
    delfin_self.write_text("# fake stub\n")
    self_mod_perms = KitToolPermissions(
        workspace=workspace, mode="bypassPermissions",
        # Auto-DENY the Self-Mod-Guard prompt — that's what a
        # vigilant user would do for a surprise rewrite.
        confirm_callback=lambda *a: False,
    )
    _doc_executor.execute("read_file", {"path": str(delfin_self)},
                           self_mod_perms)
    smg_out = json.loads(_doc_executor.execute("edit_file", {
        "path": str(delfin_self),
        "old_string": "# fake stub",
        "new_string": "# evil rewrite",
    }, self_mod_perms))
    assert "error" in smg_out
    assert "denied" in smg_out["error"].lower() or \
           "self" in smg_out["error"].lower() or \
           "protected" in smg_out["error"].lower()
    # File must NOT have been modified.
    assert "evil rewrite" not in delfin_self.read_text()

    # ---- 4. Background optimization — typical long Bayesian-opt run. -
    perms.mode = "bypassPermissions"  # equiv. to user clicking the chip
    bg = json.loads(_doc_executor.execute("bash_background", {
        "command": (
            "for i in 1 2 3 4 5; do "
            "echo 'iter '$i' / 5'; "
            "sleep 0.15; "
            "done; echo 'OPTIMIZATION DONE'"
        ),
        "description": "AX trials",
        "cwd": str(user_project),
    }, perms))
    assert bg["status"] == "started"
    job_id = bg["job_id"]
    # Read partial output while still running.
    time.sleep(0.25)
    partial = json.loads(_doc_executor.execute("bash_output",
        {"job_id": job_id}, perms))
    assert partial["running"] is True
    assert "iter 1" in partial["stdout"]
    assert "OPTIMIZATION DONE" not in partial["stdout"]
    # Wait for completion + verify final output.
    for _ in range(50):
        s = json.loads(_doc_executor.execute("bash_status",
            {"job_id": job_id}, perms))
        if not s["running"]:
            break
        time.sleep(0.1)
    assert s["exit_code"] == 0, s
    final = json.loads(_doc_executor.execute("bash_output",
        {"job_id": job_id}, perms))
    assert "OPTIMIZATION DONE" in final["stdout"]
    assert final["stdout_total_lines"] >= 6

    # ---- 5. Notebook edit — analysis report. -------------------------
    nb_path = user_project / "compare.ipynb"
    nb_initial = {
        "cells": [
            {"cell_type": "markdown", "source": "# Optimizer Comparison\n",
             "metadata": {}},
            {"cell_type": "code", "source": "results = {}\n",
             "outputs": [], "execution_count": None, "metadata": {}},
        ],
        "metadata": {}, "nbformat": 4, "nbformat_minor": 5,
    }
    nb_path.write_text(json.dumps(nb_initial))

    nb_read = json.loads(_doc_executor.execute("notebook_read",
        {"path": str(nb_path)}, perms))
    assert nb_read["cell_count"] == 2

    nb_edit = json.loads(_doc_executor.execute("notebook_edit", {
        "path": str(nb_path), "cell_idx": 1, "mode": "insert_after",
        "source": "import matplotlib.pyplot as plt\n", "cell_type": "code",
    }, perms))
    assert nb_edit["status"] == "ok"
    assert nb_edit["cells_after"] == 3
    saved = json.loads(nb_path.read_text())
    assert len(saved["cells"]) == 3

    # ---- 6. Task tracking across a "restart". ------------------------
    t1 = json.loads(_doc_executor.execute("task_create", {
        "subject": "Wire AX optimizer",
        "description": "Add mylib/optimizers/wrapper.py",
        "active_form": "Wiring AX optimizer",
    }, perms))
    t2 = json.loads(_doc_executor.execute("task_create", {
        "subject": "Compare optimizers in notebook",
    }, perms))
    assert t1["task"]["id"] != t2["task"]["id"]

    _doc_executor.execute("task_update", {
        "task_id": t1["task"]["id"], "status": "in_progress",
    }, perms)
    _doc_executor.execute("task_update", {
        "task_id": t1["task"]["id"], "status": "completed",
    }, perms)

    # Simulate "restart": drop the in-memory cache, then re-list.
    agent_tasks._STORES.clear()
    listed = json.loads(_doc_executor.execute("task_list", {}, perms))
    assert listed["count"] == 2
    statuses = {t["id"]: t["status"] for t in listed["tasks"]}
    assert statuses[t1["task"]["id"]] == "completed"
    assert statuses[t2["task"]["id"]] == "pending"
    # And the file is in the right place (project-scoped).
    tasks_file = workspace / ".delfin" / "session_tasks.json"
    assert tasks_file.exists(), (
        f"task store should live under workspace, not user home; "
        f"contents of workspace: {list(workspace.rglob('*'))}"
    )

    # ---- 7. Audit-log inspection. ------------------------------------
    log_text = audit_to_tmp.read_text()
    log_lines = [json.loads(l) for l in log_text.splitlines() if l]
    by_tool: dict[str, list[dict]] = {}
    for r in log_lines:
        by_tool.setdefault(r["tool"], []).append(r)

    # The Self-Mod-Guard denial must be in the log as 'denied'.
    smg_records = [r for r in by_tool.get("edit_file", [])
                   if r["decision"] == "denied"]
    assert smg_records, (
        f"audit log missing the Self-Mod-Guard denial; "
        f"edit_file records were: {by_tool.get('edit_file', [])}"
    )
    # The successful user-project edit must be 'ok'.
    ok_records = [r for r in by_tool.get("edit_file", [])
                  if r["decision"] == "ok"
                  and r.get("path", "").endswith("core.py")]
    assert ok_records, "audit log missing the successful user-project edit"

    # bash_background and notebook_edit are also tracked.
    assert any(r["decision"] == "ok"
               for r in by_tool.get("bash_background", []))
    assert any(r["decision"] == "ok"
               for r in by_tool.get("notebook_edit", []))


def test_safety_invariants_after_full_setup(
    workspace, user_project, audit_to_tmp, tmp_path, monkeypatch,
):
    """Even after the project_dev bundle is in place, the destructive
    pattern deny-list, secret scanner, and sandbox boundary still hold.

    These are the non-negotiable safety guarantees the user reaffirmed:
    'sicherheit ist ganz wichtig'. If any of these regress, the agent
    can do real damage; this test is the trip-wire.
    """
    from delfin.agent import kit_settings as ks
    monkeypatch.setattr(ks, "USER_SETTINGS_PATH", tmp_path / "u.json")

    perms = KitToolPermissions(
        workspace=workspace, mode="bypassPermissions",
        confirm_callback=lambda *a: True,
    )
    # Apply the bundle (most permissive realistic state).
    _doc_executor.execute("remember_permission_bundle", {
        "profile": "project_dev",
        "directory": str(user_project),
        "scope": "repo",
        "rationale": "test",
    }, perms)

    # 1. rm -rf STILL blocked.
    out = json.loads(_doc_executor.execute("bash", {
        "command": "rm -rf /tmp/something",
        "description": "evil",
    }, perms))
    assert "error" in out
    assert "deny" in out["error"].lower()

    # 2. SSH key access STILL blocked.
    out = json.loads(_doc_executor.execute("bash", {
        "command": "cat /home/user/.ssh/id_rsa",
        "description": "leak",
    }, perms))
    assert "error" in out
    assert "secret" in out["error"].lower() or "deny" in out["error"].lower()

    # 3. Force-push STILL blocked.
    out = json.loads(_doc_executor.execute("bash", {
        "command": "git push --force origin main",
        "description": "rewrite history",
    }, perms))
    assert "error" in out
    assert "deny" in out["error"].lower()

    # 4. Branch deletion STILL blocked.
    out = json.loads(_doc_executor.execute("bash", {
        "command": "git branch -D main",
        "description": "drop branch",
    }, perms))
    assert "error" in out

    # 5. Background bash inherits ALL of these checks.
    out = json.loads(_doc_executor.execute("bash_background", {
        "command": "rm -rf $HOME/.config",
        "description": "evil bg",
    }, perms))
    assert "error" in out

    # 6. Sandbox: writing to a path outside both DELFIN and user-project fails.
    foreign = tmp_path / "foreign"
    foreign.mkdir()
    out = json.loads(_doc_executor.execute("write_file", {
        "path": str(foreign / "bad.py"),
        "content": "x = 1",
    }, perms))
    assert "error" in out

    # 7. Network deny-list still applies to web tools.
    from delfin.agent import web_tools as wt
    bad_urls = [
        "http://localhost:8080/x",
        "http://169.254.169.254/latest/",
        "http://192.168.1.1/",
        "http://kit.internal/secrets",
    ]
    for url in bad_urls:
        r = wt.web_fetch(url)
        assert "error" in r, f"{url} should be blocked"


def test_full_software_build_lifecycle(
    workspace, user_project, audit_to_tmp, tmp_path, monkeypatch,
):
    """End-to-end: agent builds the user's program from scratch.

    Replays the actual software-build steps a coding agent should
    handle without manual approval once the project_dev bundle is in
    place: create venv, install dependencies, run the script, run
    pytest, observe failure, fix the source, re-run pytest, see green.

    The whole cycle uses ONLY the auto-allow patterns persisted by
    the bundle — if any of those steps unexpectedly hits the gate,
    the user's flow breaks and this test fails with a useful error.
    """
    from delfin.agent import kit_settings as ks
    monkeypatch.setattr(ks, "USER_SETTINGS_PATH", tmp_path / "u.json")

    perms = KitToolPermissions(
        workspace=workspace, mode="default",
        confirm_callback=lambda *a: True,
    )
    # Step 0: agent registers the project. -----------------------------
    bundle = json.loads(_doc_executor.execute("remember_permission_bundle", {
        "profile": "project_dev",
        "directory": str(user_project),
        "scope": "repo",
        "rationale": "Bayesian-opt project setup",
    }, perms))
    assert bundle["status"] == "persisted"

    # Step 1: write a minimal pyproject so pip install -e . works. ----
    _doc_executor.execute("read_file", {
        "path": str(user_project / "mylib" / "core.py"),
    }, perms)
    pyproj = user_project / "pyproject.toml"
    pyproj.write_text(
        "[project]\n"
        'name = "testopt"\n'
        'version = "0.1.0"\n'
        'requires-python = ">=3.9"\n'
        "[tool.setuptools.packages.find]\n"
        'where = ["."]\n'
        'include = ["mylib*"]\n'
    )
    (user_project / "mylib" / "__init__.py").write_text("")

    # Step 2: create the venv. The bundle's first pattern matches
    # 'python -m venv <path>'. Use bash with cwd=project. -------------
    out = json.loads(_doc_executor.execute("bash", {
        "command": "python -m venv .venv-user",
        "description": "create project venv",
        "cwd": str(user_project),
        "timeout_s": 60,
    }, perms))
    assert out.get("exit_code") == 0, (
        f"venv creation failed; gate or runtime issue: {out}"
    )
    venv_python = user_project / ".venv-user" / "bin" / "python"
    assert venv_python.exists(), (
        f"venv python missing at {venv_python}; "
        f"venv contents: {list((user_project / '.venv-user').rglob('*'))[:20]}"
    )

    # Step 3: install pytest + the project itself (-e .) into the venv.
    # The bundle's pip pattern matches both .venv-*/bin/pip install
    # forms; the project install puts mylib on the venv's
    # sys.path so pytest can import it from the test module. -------
    out = json.loads(_doc_executor.execute("bash", {
        "command": (
            ".venv-user/bin/pip install --quiet --disable-pip-version-check "
            "-e . pytest"
        ),
        "description": "install project + pytest",
        "cwd": str(user_project),
        "timeout_s": 180,
    }, perms))
    assert out.get("exit_code") == 0, (
        f"pip install failed (gate/runtime): {out}"
    )

    # Step 4: run the existing script through the venv python. -------
    out = json.loads(_doc_executor.execute("bash", {
        "command": (
            ".venv-user/bin/python -c "
            "'from mylib.core import compute; "
            "print(len(compute(7)))'"
        ),
        "description": "smoke run",
        "cwd": str(user_project),
        "timeout_s": 30,
    }, perms))
    assert out.get("exit_code") == 0, f"script run failed: {out}"
    assert "7" in (out.get("stdout") or ""), (
        f"expected '7' in stdout; got: {out.get('stdout')}"
    )

    # Step 5: run pytest on the existing test, expect green. ---------
    out = json.loads(_doc_executor.execute("bash", {
        "command": ".venv-user/bin/pytest -q",
        "description": "run tests",
        "cwd": str(user_project),
        "timeout_s": 60,
    }, perms))
    assert out.get("exit_code") == 0, (
        f"initial pytest failed: stdout={out.get('stdout', '')[:300]} "
        f"stderr={out.get('stderr', '')[:300]}"
    )

    # Step 6: agent INTRODUCES a bug, sees pytest fail. --------------
    src = user_project / "mylib" / "core.py"
    _doc_executor.execute("read_file", {"path": str(src)}, perms)
    edit = _doc_executor.execute("edit_file", {
        "path": str(src),
        "old_string": "for i in range(n_trials)",
        "new_string": "for i in range(n_trials - 1)",  # off-by-one bug
    }, perms)
    assert "Edited" in edit
    # The auto-test hint should still point at the right test.
    assert "tests/test_core.py" in edit

    out = json.loads(_doc_executor.execute("bash", {
        "command": ".venv-user/bin/pytest -q",
        "description": "verify bug",
        "cwd": str(user_project),
        "timeout_s": 60,
    }, perms))
    assert out.get("exit_code") != 0, (
        "test should have failed after the off-by-one bug, but it "
        f"passed: {out}"
    )
    # The traceback contains 'AssertionError' — agent can see it.
    combined = (out.get("stdout") or "") + (out.get("stderr") or "")
    assert "AssertionError" in combined or "assert" in combined.lower(), (
        f"expected an assertion error in pytest output; got: {combined[:300]}"
    )

    # Step 7: agent fixes the bug, re-runs, sees green again. --------
    _doc_executor.execute("read_file", {"path": str(src)}, perms)
    fix = _doc_executor.execute("edit_file", {
        "path": str(src),
        "old_string": "for i in range(n_trials - 1)",
        "new_string": "for i in range(n_trials)",
    }, perms)
    assert "Edited" in fix

    out = json.loads(_doc_executor.execute("bash", {
        "command": ".venv-user/bin/pytest -q",
        "description": "verify fix",
        "cwd": str(user_project),
        "timeout_s": 60,
    }, perms))
    assert out.get("exit_code") == 0, (
        f"tests should be green after the fix; got: {out}"
    )

    # Step 8: long-running 'optimization' kicked off in the background
    # while the agent moves on. Verifies the bash_background tool is
    # also gated by the project_dev bundle, not just sync bash. ------
    bg = json.loads(_doc_executor.execute("bash_background", {
        "command": (
            ".venv-user/bin/python -c "
            "'import time; "
            "[print(f\"trial {i}\", flush=True) or time.sleep(0.05) "
            "for i in range(5)]'"
        ),
        "description": "background opt",
        "cwd": str(user_project),
    }, perms))
    assert bg["status"] == "started", bg
    job_id = bg["job_id"]
    # Wait for it.
    for _ in range(50):
        s = json.loads(_doc_executor.execute("bash_status",
            {"job_id": job_id}, perms))
        if not s["running"]:
            break
        time.sleep(0.1)
    out = json.loads(_doc_executor.execute("bash_output",
        {"job_id": job_id}, perms))
    assert out["exit_code"] == 0
    assert "trial 0" in out["stdout"]
    assert "trial 4" in out["stdout"]


def test_audit_log_does_not_break_tool_path(
    workspace, monkeypatch, tmp_path,
):
    """If the audit log can't be written (disk full, permission denied),
    the tool call must still succeed — audit is observability, not
    a gate. Simulate by pointing the log at an unwritable path."""
    # Read-only directory for the log target.
    bad_dir = tmp_path / "ro"
    bad_dir.mkdir()
    bad_dir.chmod(0o500)
    bad_log = bad_dir / "subdir" / "audit.log"
    monkeypatch.setattr(al, "_default_log_path", lambda: bad_log)

    try:
        perms = KitToolPermissions(workspace=workspace, mode="default")
        f = workspace / "x.py"
        f.write_text("a = 1\n")
        _doc_executor.execute("read_file", {"path": str(f)}, perms)
        out = _doc_executor.execute("edit_file", {
            "path": str(f), "old_string": "a = 1", "new_string": "a = 2",
        }, perms)
        # Edit must have succeeded despite the audit failure.
        assert "Edited" in out
        assert "a = 2" in f.read_text()
    finally:
        bad_dir.chmod(0o700)
