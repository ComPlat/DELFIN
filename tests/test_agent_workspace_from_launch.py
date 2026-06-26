"""The agent must build where the user LAUNCHED delfin-voila, not where the
notebook lives. Voila pins the kernel cwd to the notebook dir (inside the
delfin checkout), so the agent reads the real launch dir from DELFIN_LAUNCH_CWD.
Bug 2026-06-26: launched from ~/agent_workspace, the agent still built in
software/delfin. ONLY the agent's workspace changes — ctx.repo_dir is untouched."""
from pathlib import Path
from delfin.dashboard.tab_agent import _agent_workspace_from_launch as R


def _mk_delfin_tree(root: Path) -> Path:
    (root / "delfin").mkdir(parents=True)
    (root / "delfin" / "__init__.py").write_text("")
    return root


def test_launch_inside_delfin_tree_works_on_delfin(tmp_path):
    repo = _mk_delfin_tree(tmp_path / "software" / "delfin")
    assert R(str(repo), tmp_path) == repo
    # a SUBDIR of the checkout still resolves up to the checkout root
    sub = repo / "delfin" / "dashboard"
    sub.mkdir(parents=True, exist_ok=True)
    assert R(str(sub), tmp_path) == repo


def test_launch_in_plain_project_builds_there(tmp_path):
    ws = tmp_path / "agent_workspace"
    ws.mkdir()
    # NOT a delfin tree → the agent builds right where it was launched
    assert R(str(ws), tmp_path / "fallback") == ws


def test_empty_launch_cwd_uses_fallback(tmp_path):
    fb = tmp_path / "fb"
    assert R("", fb) == fb
    assert R("   ", fb) == fb
