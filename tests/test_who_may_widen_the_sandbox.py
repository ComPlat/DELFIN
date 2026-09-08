"""Every place that makes a new directory writable, enumerated.

The workspace boundary is the agent's main containment. Widening it is a
decision with consequences that outlast the turn, so the ways to do it
are few and each is either the user's own act or gated on their explicit
confirmation:

  * ``remember_permission(kind='extra_dir')`` and
    ``remember_permission_bundle`` — the sanctioned grants, both refused
    unless ``confirm_callback`` returns true. Without a callback there is
    nobody to ask and they refuse.
  * ``enter_worktree`` — derives a writable root from a repository the
    session already works in, and refuses any other (see
    test_a_worktree_is_not_a_way_out_of_the_workspace.py).
  * ``AgentEngine.grant_directory`` — the dashboard's "Erlaubte
    Verzeichnisse" panel. Not a tool: no model can call it.

A fourth appeared once and nobody noticed: enter_worktree took a
repo_dir the model named, made a worktree of it and registered that as
writable, turning a refused path into a writable copy in one call. This
file is the counting rule that would have caught it — the same shape as
test_a_checked_out_repository_cannot_run_commands.py's enumeration of who
may grant trust.
"""

from __future__ import annotations

import pathlib
import re

_ROOT = pathlib.Path(__file__).resolve().parents[1]

# Where the call is allowed to appear, and what makes each one safe.
_ALLOWED = {
    # The definition itself.
    "delfin/agent/api_client.py",
    # The engine method the dashboard panel calls. Not a tool.
    "delfin/agent/engine.py",
    # The panel.
    "delfin/dashboard/tab_agent.py",
}


def _call_sites(name: str) -> set[str]:
    hits = set()
    for path in sorted(_ROOT.joinpath("delfin").rglob("*.py")):
        for line in path.read_text(encoding="utf-8",
                                   errors="replace").splitlines():
            stripped = line.lstrip()
            if stripped.startswith("#") or f"{name}(" not in line:
                continue
            if re.search(rf"def\s+{name}\s*\(", line):
                continue
            hits.add(str(path.relative_to(_ROOT)))
    return hits


def test_nothing_new_widens_the_sandbox():
    assert _call_sites("add_extra_dir") <= _ALLOWED, (
        f"a new place makes a directory writable: "
        f"{sorted(_call_sites('add_extra_dir') - _ALLOWED)}")


def test_the_two_sanctioned_tools_ask_first():
    """Both grants sit behind confirm_callback in the same function, so
    a refusal or a missing callback stops them before anything is
    written."""
    src = (_ROOT / "delfin" / "agent" / "api_client.py").read_text(
        encoding="utf-8")
    for anchor in ("_execute_remember_permission",
                   "_execute_remember_permission_bundle"):
        i = src.find(f"def {anchor}")
        assert i > 0, anchor
        body = src[i:i + 6000]
        add = body.find("add_extra_dir(")
        confirm = body.find("confirm_callback(")
        assert confirm > 0, f"{anchor} does not ask"
        assert add < 0 or confirm < add, (
            f"{anchor} widens the sandbox before asking")


def test_a_headless_session_cannot_be_talked_into_a_grant():
    """No callback means nobody to ask, and the answer is no rather than
    yes — which is what a benchmark run and the scheduler are."""
    import json
    import tempfile

    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    ws = pathlib.Path(tempfile.mkdtemp())
    target = ws.parent / "somewhere-else"
    target.mkdir(exist_ok=True)
    perms = KitToolPermissions(workspace=ws)
    assert perms.confirm_callback is None
    out = json.loads(_DocToolExecutor()._execute_remember_permission(
        {"kind": "extra_dir", "value": str(target), "scope": "user"}, perms))
    assert out.get("status") != "ok", out
    assert perms.extra_workspace_dirs == ()
