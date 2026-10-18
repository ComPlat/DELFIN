"""Phase-2 red control for delfin.agent.refusal_memory.same_change.

same_change(refusal, tool, args, workspace, write_targets=None) is a new
pure module function: "does this tool call CHANGE (write) the file the
person refused?". Built on _targets_of (file tools) and the caller-injected
write-target tokenizer (bash -- the gate's own _bash_write_targets, so no
second parser). It never allows; outputs are None ('no') or the target it
would change.

Control: this file collects as ERROR on the previous commit, because
same_change does not exist yet (ImportError at the from-import). Every case
below is RED with the module unchanged.
"""
from pathlib import Path

from delfin.agent.api_client import _bash_write_targets
from delfin.agent.refusal_memory import Refusal, same_change

WS = "/home/u/workspace"


def _refused_write(path, tool="bash", reason="the user said no to this"):
    return Refusal(
        tool=tool,
        target=str(Path(path).resolve()),
        reason=reason,
        time="2026-10-16T00:00:00+00:00",
    )


# ---- bash write to a refused file, another spelling -------------------


class TestBashWriteSameFile:
    def test_redirect_changes_refused_file(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "bash",
            {"command": f"echo x > {WS}/f.txt"}, WS,
            write_targets=_bash_write_targets,
        ) == str(Path(WS + "/f.txt").resolve())

    def test_sed_in_place_changes_refused_file(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "bash",
            {"command": f"sed -i s/a/b/ {WS}/f.txt"}, WS,
            write_targets=_bash_write_targets,
        ) == str(Path(WS + "/f.txt").resolve())

    def test_tee_changes_refused_file(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "bash",
            {"command": f"echo x | tee {WS}/f.txt"}, WS,
            write_targets=_bash_write_targets,
        ) == str(Path(WS + "/f.txt").resolve())

    def test_python_open_w_changes_refused_file(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "bash",
            {"command": f"python3 -c 'open(\"{WS}/f.txt\",\"w\").write(\"x\")'"},
            WS, write_targets=_bash_write_targets,
        ) == str(Path(WS + "/f.txt").resolve())

    def test_heredoc_changes_refused_file(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "bash",
            {"command": f"cat > {WS}/f.txt <<EOF\ncontent\nEOF"}, WS,
            write_targets=_bash_write_targets,
        ) == str(Path(WS + "/f.txt").resolve())

    def test_append_redirect_changes_refused_file(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "bash",
            {"command": f"echo x >> {WS}/f.txt"}, WS,
            write_targets=_bash_write_targets,
        ) == str(Path(WS + "/f.txt").resolve())


# ---- a DIFFERENT file is not a change of the refused one -------------


class TestDifferentFile:
    def test_write_to_a_different_file_is_no_change(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "bash",
            {"command": f"echo x > {WS}/other.txt"}, WS,
            write_targets=_bash_write_targets,
        ) is None

    def test_write_file_to_a_different_file_is_no_change(self):
        refused = _refused_write(WS + "/f.txt", tool="write_file")
        assert (
            same_change(refused, "write_file",
                        {"path": WS + "/other.txt"}, WS)
            is None
        )


# ---- reading is not changing ------------------------------------------


class TestReadingIsNotChanging:
    def test_read_bash_is_no_change(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "bash",
            {"command": f"cat {WS}/f.txt"}, WS,
            write_targets=_bash_write_targets,
        ) is None

    def test_read_file_is_no_change(self):
        refused = _refused_write(WS + "/f.txt", tool="write_file")
        assert (
            same_change(refused, "read_file", {"path": WS + "/f.txt"}, WS)
            is None
        )


# ---- write_file / edit_file to the refused file -----------------------


class TestDocToolWrite:
    def test_write_file_to_refused_file_is_a_change(self):
        refused = _refused_write(WS + "/f.txt", tool="write_file")
        assert (
            same_change(refused, "write_file",
                        {"path": WS + "/f.txt", "content": "new"}, WS)
            == str(Path(WS + "/f.txt").resolve())
        )

    def test_edit_file_to_refused_file_is_a_change(self):
        refused = _refused_write(WS + "/f.txt", tool="edit_file")
        assert (
            same_change(refused, "edit_file",
                        {"path": WS + "/f.txt", "old_string": "a"}, WS)
            == str(Path(WS + "/f.txt").resolve())
        )


# ---- directory refusal containment ------------------------------------


class TestDirRefusal:
    def test_write_inside_refused_dir_is_a_change(self):
        refused = _refused_write(WS + "/secrets", tool="bash")
        refused.is_dir = True
        assert same_change(
            refused, "bash",
            {"command": f"echo x > {WS}/secrets/k.txt"}, WS,
            write_targets=_bash_write_targets,
        ) == str(Path(WS + "/secrets/k.txt").resolve())

    def test_write_outside_refused_dir_is_no_change(self):
        refused = _refused_write(WS + "/secrets", tool="bash")
        refused.is_dir = True
        assert same_change(
            refused, "bash",
            {"command": f"echo x > {WS}/public/k.txt"}, WS,
            write_targets=_bash_write_targets,
        ) is None


# ---- path normalisation: relative vs absolute vs dot-dot -------------


class TestNormalisation:
    def test_relative_path_resolves_to_refused_absolute(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "write_file", {"path": "f.txt"}, WS,
        ) == str(Path(WS + "/f.txt").resolve())

    def test_dotdot_normalises_to_refused_absolute(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "write_file",
            {"path": WS + "/sub/../f.txt"}, WS,
        ) == str(Path(WS + "/f.txt").resolve())


# ---- write_targets is None: a bash attempt cannot be judged -----------


class TestNoExtractor:
    def test_bash_without_extractor_returns_none(self):
        # Caller always passes the real tokenizer; without it we cannot know
        # the write targets, so we return None rather than guess.
        refused = _refused_write(WS + "/f.txt")
        assert same_change(
            refused, "bash", {"command": f"echo x > {WS}/f.txt"}, WS,
        ) is None

    def test_file_tool_needs_no_extractor(self):
        refused = _refused_write(WS + "/f.txt", tool="write_file")
        assert (
            same_change(refused, "write_file", {"path": WS + "/f.txt"}, WS)
            == str(Path(WS + "/f.txt").resolve())
        )


# ---- odd / hostile args never raise ------------------------------------


class TestNeverRaises:
    def test_empty_args(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(refused, "bash", {}, WS) is None
        assert same_change(refused, "write_file", {}, WS) is None

    def test_non_dict_args(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(refused, "bash", None, WS) is None
        assert same_change(refused, "write_file", "not-a-dict", WS) is None

    def test_no_such_tool(self):
        refused = _refused_write(WS + "/f.txt")
        assert same_change(refused, "janitor_service", {}, WS) is None
