"""U8 (operator build): what a git/patch command writes, case by case."""

import pytest

from delfin.agent import patch_targets as pt

_DIFF = ("diff --git a/{n} b/{n}\n--- a/{n}\n+++ b/{n}\n"
         "@@ -1 +1 @@\n-x\n+y\n")


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "a.diff").write_text(_DIFF.format(n="a.py"))
    (tmp_path / "b.diff").write_text(_DIFF.format(n="b.py"))
    (tmp_path / "afile").write_text("x")
    return tmp_path


def test_every_patch_file_named_is_read(ws):
    assert pt.targets("git apply a.diff b.diff", ws) == [
        str(ws / "a.py"), str(ws / "b.py")]


def test_a_patch_file_is_read_where_the_command_runs_not_where_we_stand(
        ws, tmp_path_factory, monkeypatch):
    monkeypatch.chdir(tmp_path_factory.mktemp("elsewhere"))
    assert pt.targets("git apply a.diff", ws) == [str(ws / "a.py")]
    assert pt.targets("patch -p1 -i a.diff", ws) == [str(ws / "a.py")]


def test_a_patch_that_cannot_be_read_is_unknowable(ws):
    assert pt.targets("git apply missing.diff", ws) is pt.UNKNOWABLE
    assert pt.targets("git apply a.diff missing.diff", ws) is pt.UNKNOWABLE


@pytest.mark.parametrize("cmd", ["git checkout main", "git checkout -b new",
                                 "git switch main"])
def test_a_branch_switch_writes_no_named_file(ws, cmd):
    assert pt.targets(cmd, ws) == []


def test_after_the_double_dash_every_operand_is_a_path(ws):
    assert pt.targets("git checkout -- main", ws) == [str(ws / "main")]
    assert pt.targets("git checkout HEAD~1 -- a b", ws) == [
        str(ws / "a"), str(ws / "b")]


def test_without_the_dash_an_existing_or_path_shaped_operand_is_a_path(ws):
    assert pt.targets("git checkout afile", ws) == [str(ws / "afile")]
    assert pt.targets("git checkout src/x.py", ws) == [str(ws / "src/x.py")]


@pytest.mark.parametrize("cmd", ["git --work-tree=/x checkout -- f",
                                 "git -C /x checkout -- f",
                                 "git -C /x restore f"])
def test_a_named_tree_is_where_the_paths_resolve(ws, cmd):
    assert pt.targets(cmd, ws) == ["/x/f"]


def test_index_only_forms_write_no_work_tree_file(ws):
    assert pt.targets("git apply --cached a.diff", ws) == []
    assert pt.targets("git apply --check a.diff", ws) == []
    assert pt.targets("git restore --staged afile", ws) == []


@pytest.mark.parametrize("flag", ["--3way", "--index", "-3"])
def test_three_way_and_index_still_write_the_work_tree(ws, flag):
    assert pt.targets(f"git apply {flag} a.diff", ws) == [str(ws / "a.py")]


def test_the_gate_refuses_a_restore_outside_the_workspace(tmp_path):
    from delfin.agent.api_client import KitToolPermissions, _doc_executor
    work = tmp_path / "ws"
    work.mkdir()
    perms = KitToolPermissions(workspace=work, mode="bypassPermissions",
                               confirm_callback=None)

    def gate(cmd):
        return _doc_executor._gate_bash_write_targets(
            cmd, {"command": cmd}, perms)

    outside = tmp_path / "outside.txt"
    assert gate(f"git checkout -- {outside}") is not None
    assert gate(f"git -C {tmp_path} restore outside.txt") is not None
    assert gate("git checkout -- inside.txt") is None
    assert gate("git checkout main") is None
