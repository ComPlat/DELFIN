"""pre_pr_check routes every git call through delfin.agent.safe_git.

RED CONTROL, written before the fix: on the pre-fix module the three
call sites (_git, _rev_exists, run's merge-base probe) call
subprocess.run(["git", ...]) directly — exactly what the operator's
invariant test tests/test_delfins_own_git_reads_go_through_safe_git.py
names. This test pins the DELFIN-owned route, so it must pass only
after the routing fix.
"""

import ast
from pathlib import Path


def _module_source() -> str:
    path = Path(__file__).resolve().parent.parent / (
        "delfin/agent/pre_pr_check.py")
    return path.read_text(encoding="utf-8")


def test_no_subprocess_run_git_site_remains_in_pre_pr_check():
    tree = ast.parse(_module_source())
    direct = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "run"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "subprocess"):
            for arg in node.args:
                if (isinstance(arg, (ast.List, ast.Tuple))
                        and arg.elts
                        and isinstance(arg.elts[0], ast.Constant)
                        and arg.elts[0].value == "git"):
                    direct.append(node.lineno)
    assert not direct, (
        "pre_pr_check runs git directly (route through "
        f"delfin.agent.safe_git): lines {direct}")


def test_pre_pr_check_imports_safe_git():
    source = _module_source()
    assert "from delfin.agent import safe_git" in source, (
        "pre_pr_check must use the same safe-git route as the rest of "
        "DELFIN")
