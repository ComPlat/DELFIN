"""Heredoc tags with '.' or '-' break the bash read gate.

The read scan blanks heredoc bodies before scanning (`_strip_heredoc_bodies`,
api_client.py). Its tag regex only accepted ``\\w+`` tags, so real bash tags
like ``<<'PY.SCRIPT'`` or ``<<EOF-X`` were not recognised, the body was not
blanked, and an absolute path inside the body fired the read gate -- a false
block of a command the shell would treat as data-only.

Tested through the public path: the read-scan collector ``_bash_outside_reads``
(what ``_gate_bash_read_paths`` consumes). A command outside the heredoc that
reads an absolute path must still be gated.

Licence: GPL-3.0-or-later
"""
import pytest

from delfin.agent import api_client as ac


def _reads(cmd: str) -> list[str]:
    return ac._bash_outside_reads(cmd, cwd=None)


# --- the defect: dotted / dashed tags, quoted and unquoted, <<- with tabs ---

@pytest.mark.parametrize("cmd", [
    # quoted dotted tag, body holds an absolute path as DATA
    "cat <<'PY.SCRIPT'\nprint('x')\nls /etc/passwd\nPY.SCRIPT\n",
    # quoted dashed tag
    "cat <<'EOF-FILE'\ncat /etc/passwd\nEOF-FILE\n",
    # <<- dashes tag, body indented with TABS, closing line tab-prefixed
    "cat <<-'CONF.BASE'\n\tread /etc/passwd\n\tCONF.BASE\n",
    # unquoted dashed tag
    "cat <<EOF-X\nls /etc/passwd\nEOF-X\n",
    # dotted unquoted tag
    "cat <<CONF.B\ncat /etc/passwd\nCONF.B\n",
])
def test_body_with_absolute_path_is_data_not_a_read(cmd):
    assert _reads(cmd) == [], (
        f"heredoc body leaked into the read scan for: {cmd!r}")


@pytest.mark.parametrize("cmd", [
    # unquoted dashed tag whose body holds a line equal to the WRONG short
    # tag the old regex would match: blanking must stop at the true tag
    # EOF-X, and the body after `EOF` must still be blanked (old code
    # matched tag `EOF`, stopped there, and leaked the rest)
    "cat <<EOF-X\nEOF\nls /etc/passwd\nEOF-X\n",
    "cat <<CONF.B\nCONF\nls /etc/passwd\nCONF.B\n",
])
def test_unquoted_dashed_tag_closes_at_the_true_tag(cmd):
    assert _reads(cmd) == [], (
        f"body after the true tag's delimiter leaked: {cmd!r}")


# --- negative control: outside the heredoc the path must still be gated ---

def test_absolute_path_outside_heredoc_is_still_a_read():
    assert _reads("cat /etc/passwd") == ["/etc/passwd"]


def test_path_after_redirect_line_is_still_a_read():
    # the rest of the redirect's own line is KEPT by design
    # (api_client.py:9219-9226): a second command after the heredoc must
    # still be scanned -- `python - <<'PY' ... PY; cat /outside/secret`
    cmd = "python - <<'PY'\nprint('x')\nPY\ncat /outside/secret\n"
    assert _reads(cmd) == ["/outside/secret"]


def test_redirect_target_after_tag_not_a_read_target_but_write_gated():
    # `cat <<EOF > /outside/x` -- the read scan skips redirect targets by
    # design (_segment_reads, api_client.py:2044-2049); the WRITE gate owns
    # that path. The body below the redirect line must stay blanked so the
    # body's words are never scanned as reads.
    cmd = "cat <<EOF > /outside/x\n/etc/passwd\nEOF\n"
    assert _reads(cmd) == []


# --- closing-line rule: the body ends only at a line equal to the tag ---

def test_unrelated_word_does_not_close_the_body():
    # 'PY.SCRIPT' appears as a substring word but is not the whole line:
    # the body continues, so its path stays blanked.
    cmd = ("cat <<'PY.SCRIPT'\n/etc/passwd\nPY.SCRIPTS\nPY.SCRIPT\n")
    assert _reads(cmd) == []


def test_quoted_closing_line_closes_the_body():
    # bash accepts a closing line that merely matches the tag's word
    cmd = "cat <<'EOF-FILE'\n/etc/passwd\n'EOF-FILE'\nEOF-FILE\n"
    assert _reads(cmd) == []
