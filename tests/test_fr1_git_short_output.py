"""Control for the FR1 parser baseline (reviewer nacht-s17's finding,
verified by s10): git show/log/diff accept the short form ``-o FILE``
of ``--output=FILE``. _GIT_WRITING_OPTION_PREFIXES only knows the long
form, so a file-writing git command passes as a plain read."""

from delfin.agent.action_classifier import is_plain_read


def test_git_show_short_output_writes_a_file():
    action = {"tool": "bash",
              "command": "git show -o /tmp/evil HEAD --stat"}
    assert not is_plain_read(action)


def test_git_log_short_output_writes_a_file():
    action = {"tool": "bash", "command": "git log -o out.patch -1"}
    assert not is_plain_read(action)


def test_git_diff_attached_short_output_writes_a_file():
    action = {"tool": "bash", "command": "git diff -oout.patch main"}
    assert not is_plain_read(action)


def test_long_output_form_still_not_a_plain_read():
    # The long form is already covered; pinned so the fix keeps it.
    action = {"tool": "bash",
              "command": "git show --output=out.patch HEAD --stat"}
    assert not is_plain_read(action)


def test_git_show_without_output_still_a_plain_read():
    action = {"tool": "bash", "command": "git show HEAD --stat"}
    assert is_plain_read(action)
