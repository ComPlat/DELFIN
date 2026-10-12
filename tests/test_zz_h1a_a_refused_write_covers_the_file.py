"""H1A (operator build): a refused write covers the file, not the spelling.

Through the real gate, with a dialog that records what it was asked.
"""

import pytest

from delfin.agent import api_client as A


class Dialog:
    last_timed_out = False

    def __init__(self, answer):
        self.answer = answer
        self.asked = []

    def ask(self, name, args, preview):
        self.asked.append((name, preview))
        return self.answer


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "f.txt").write_text("old\n")
    (tmp_path / "a").write_text("src\n")
    (tmp_path / "sub").mkdir()
    return tmp_path


def _session(ws, mode="default", answer=False):
    dialog = Dialog(answer)
    perms = A.KitToolPermissions(
        workspace=ws, mode=mode, confirm_callback=dialog.ask)
    return A._DocToolExecutor(), perms, dialog


def _bash(ex, perms, cmd, **args):
    return ex._run_permission_gate("bash", dict(command=cmd, **args), perms)


RESPELLED = [
    "tee f.txt", "cat a | tee -a f.txt", "cp a f.txt", "mv a f.txt",
    "install -m 644 a f.txt", "truncate -s 0 f.txt",
    "curl -o f.txt http://example.org/x", "ln -sf a f.txt",
    "./tool --run > f.txt", "./tool --run >> ./f.txt",
    "./tool --run > sub/../f.txt", "rm f.txt",
]


@pytest.mark.parametrize("second", RESPELLED)
def test_the_same_file_in_another_spelling_is_refused_without_a_dialog(ws, second):
    ex, perms, dialog = _session(ws)
    assert "user denied" in _bash(ex, perms, "./tool --run > f.txt.tmp && mv f.txt.tmp f.txt")
    asked = len(dialog.asked)
    out = _bash(ex, perms, second)
    assert out and "already refused" in out, second
    assert "f.txt" in out
    assert len(dialog.asked) == asked, second


def test_it_holds_in_the_modes_that_ask_nothing(ws):
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool > f.txt") is not None
    for mode in ("acceptEdits", "bypassPermissions"):
        perms.mode = mode
        out = _bash(ex, perms, "cp a f.txt")
        assert out and "already refused" in out, mode
    assert len(dialog.asked) == 1


def test_reading_the_file_is_not_what_was_refused(ws):
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool > f.txt") is not None
    assert _bash(ex, perms, "cat f.txt") is None
    assert _bash(ex, perms, "grep old f.txt") is None
    assert len(dialog.asked) == 1


def test_another_file_is_asked_about_as_before(ws):
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool > f.txt") is not None
    dialog.answer = True
    assert _bash(ex, perms, "./tool > g.txt") is None
    assert len(dialog.asked) == 2
    assert _bash(ex, perms, "./tool > sub/f.txt") is None
    assert len(dialog.asked) == 3


def test_a_refused_command_that_writes_nothing_closes_no_file(ws):
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool --run f.txt") is not None
    dialog.answer = True
    # An arrow inside a commit message is text, not a write.
    assert _bash(ex, perms, "git commit -m 'a > f.txt'") is None
    assert _bash(ex, perms, "./tool --other > f.txt") is None
    assert len(dialog.asked) == 2


def test_the_directory_the_command_ran_in_decides_which_file(ws):
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool > f.txt", cwd="sub") is not None
    dialog.answer = True
    # sub/f.txt was refused, the workspace's own f.txt was not.
    assert _bash(ex, perms, "./tool2 > f.txt") is None
    out = _bash(ex, perms, "./tool2 > sub/f.txt")
    assert out and "already refused" in out


def test_a_write_tool_after_a_refused_shell_write_is_asked_never_silent(ws):
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool > f.txt") is not None
    perms.mode = "acceptEdits"
    args = {"path": "f.txt", "content": "new\n"}
    out = ex._run_permission_gate("write_file", args, perms)
    assert out is not None and "user denied" in out
    assert len(dialog.asked) == 2
    assert "REFUSED EARLIER" in dialog.asked[1][1]
    # ... and with nobody to ask it is refused, not written.
    perms.confirm_callback = None
    perms.denied_actions.pop(next(
        k for k in perms.denied_actions if k.startswith("write_change:")))
    assert ex._run_permission_gate("write_file", args, perms) is not None


def test_a_different_change_is_shown_and_can_be_approved(ws):
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool > f.txt") is not None
    perms.mode = "acceptEdits"
    dialog.answer = True
    assert ex._run_permission_gate(
        "write_file", {"path": "f.txt", "content": "agreed\n"}, perms) is None
    assert len(dialog.asked) == 2
    # An untouched file needs no dialog in acceptEdits, as before.
    assert ex._run_permission_gate(
        "write_file", {"path": "g.txt", "content": "x"}, perms) is None
    assert len(dialog.asked) == 2


def test_the_exact_command_stays_refused_as_before(ws):
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool > f.txt") is not None
    out = _bash(ex, perms, "./tool > f.txt")
    assert out and "this exact command" in out
    assert len(dialog.asked) == 1


def test_a_timeout_is_absence_and_closes_nothing(ws):
    ex, perms, dialog = _session(ws)
    dialog.last_timed_out = True
    perms.confirm_callback = dialog.ask
    assert "TIMED OUT" in _bash(ex, perms, "./tool > f.txt")
    dialog.last_timed_out = False
    dialog.answer = True
    assert _bash(ex, perms, "cp a f.txt") is None


def test_a_cd_in_the_command_decides_which_file(ws):
    """`cd sub && ... > f` writes sub/f: that is the file refused, not f."""
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "cd sub && ./tool > f.txt") is not None
    dialog.answer = True
    assert _bash(ex, perms, "./tool2 > f.txt") is None
    out = _bash(ex, perms, "./tool2 > sub/f.txt")
    assert out and "already refused" in out
    out = _bash(ex, perms, "cd sub && ./tool2 > ../sub/f.txt")
    assert out and "already refused" in out


def test_a_git_restore_of_the_file_is_the_same_write(ws):
    """patch_targets names the paths of git checkout/restore; the write
    gate gates them, so the refusal reads them too -- both ways."""
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool > f.txt") is not None
    out = _bash(ex, perms, "git checkout -- f.txt")
    assert out and "already refused" in out
    assert len(dialog.asked) == 1
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "git restore f.txt") is not None
    out = _bash(ex, perms, "./tool > f.txt")
    assert out and "already refused" in out
    assert len(dialog.asked) == 1


def test_a_refusal_in_the_classifiers_dialog_closes_the_file_too(
        ws, monkeypatch):
    monkeypatch.setattr(
        A, "_classifier_settings",
        lambda: {"model": "m", "rules": [], "broad": True})
    ex, perms, dialog = _session(ws, mode="bypassPermissions")
    perms.classifier_runner = (
        lambda prompt, model:
        "DECISION: ask\nREASON: the request does not mention this.")
    perms.classifier_messages = [{"role": "user", "content": "look around"}]
    assert _bash(ex, perms, "./tool > f.txt") is not None
    assert len(dialog.asked) == 1
    out = _bash(ex, perms, "cp a f.txt")
    assert out and "already refused" in out
    assert len(dialog.asked) == 1


def test_a_refusal_never_closes_the_workspace_itself(ws):
    """`cp a .` names the workspace as its destination; refusing it must
    not refuse every later copy into the workspace without a dialog."""
    ex, perms, dialog = _session(ws)
    assert _bash(ex, perms, "./tool && cp ../x .") is not None
    dialog.answer = True
    assert _bash(ex, perms, "./tool && cp a .") is None
    assert len(dialog.asked) == 2
