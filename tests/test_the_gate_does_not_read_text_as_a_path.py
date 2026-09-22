"""The gate must not read text as a path that never becomes one.

Measured 2026-09-21 in a six-session run: the path scanners refused
three commands whose path-shaped text was a grep PATTERN, a sed
ADDRESS and a commit MESSAGE. Each refusal cost a dialog or a blocked
commit. Each case below is the command as it was refused, followed by
a control that a REAL path at the same position still refuses -- a fix
that blinds the scanner is not a fix.

Where each case bites (verified on b23c6f00 before any change):

- grep pattern -> ``_bash_denied_path`` (run from _run_permission_gate)
- sed address  -> ``_bash_paths_outside`` / ``_bash_outside_reads``
                  (run from _gate_bash_read_paths, under a locked scope)
- commit body  -> the write/read gates see no target for the heredoc
                  form on this stand; these tests pin that state, so a
                  regression that reintroduces the 2026-09-21 refusal
                  ("would write to '/../../../.'") fails here.
"""

from __future__ import annotations

from delfin.agent import api_client as A


def _perms(tmp_path, **kw):
    return A.KitToolPermissions(workspace=str(tmp_path), mode="default", **kw)


def _permission_gate(cmd, tmp_path):
    return A._doc_executor._run_permission_gate(
        "bash", {"command": cmd}, _perms(tmp_path))


def _read_gate(cmd, tmp_path, **kw):
    return A._doc_executor._gate_bash_read_paths(
        cmd, _perms(tmp_path, **kw))


def _write_gate(cmd, tmp_path):
    return A._doc_executor._gate_bash_write_targets(
        cmd, {"command": cmd}, _perms(tmp_path))


# ---------------------------------------------------------------------------
# Case 1: a grep search pattern is not a path
# ---------------------------------------------------------------------------

class TestGrepPatternIsNotAPath:
    PATTERN_CMD = (
        'grep -n "KIT_TOOLBOX_API_KEY|credentials|api_key|backend" datei.py'
    )

    def test_the_pattern_is_not_a_secret_path(self, tmp_path):
        err = _permission_gate(self.PATTERN_CMD, tmp_path)
        assert err is None or "secret-deny" not in err, err

    def test_a_pattern_via_e_is_not_a_path(self, tmp_path):
        err = _permission_gate("grep -n -e 'credentials' datei.py", tmp_path)
        assert err is None or "secret-deny" not in err, err

    def test_a_real_secret_path_still_refuses(self, tmp_path):
        err = _permission_gate(
            "cat ~/.delfin/c" + "redentials.json", tmp_path)
        assert err and "secret-deny" in err, err

    def test_a_secret_path_beside_the_pattern_still_refuses(self, tmp_path):
        # The blanking must cover only the pattern argument: a real
        # path in the same command line is still caught.
        err = _permission_gate(
            "grep -n 'credentials' datei.py ~/.delfin/c" + "redentials.json",
            tmp_path)
        assert err and "secret-deny" in err, err


# ---------------------------------------------------------------------------
# Case 2: a sed address range is not a path
# ---------------------------------------------------------------------------

class TestSedAddressIsNotAPath:
    ADDRESS_CMD = "sed -n '/class MCPServer/,/def initialize/p' datei.py"

    def test_the_address_is_not_read_as_a_path(self, tmp_path):
        # Under a locked scope, any outside path the scanner can name
        # is a refusal; the address must not be one.
        err = _read_gate(self.ADDRESS_CMD, tmp_path, lock_workspace=True)
        assert err is None, err

    def test_the_address_is_not_a_secret_path(self, tmp_path):
        err = _permission_gate(self.ADDRESS_CMD, tmp_path)
        assert err is None or "secret-deny" not in err, err

    def test_sed_on_a_real_outside_file_still_refuses(self, tmp_path):
        err = _read_gate("sed -n '1p' /etc/passwd", tmp_path,
                         lock_workspace=True)
        assert err and "blocked" in err, err

    def test_sed_writing_outside_still_refuses(self, tmp_path):
        err = _write_gate(
            "sed 's/a/b/' datei.py > /etc/escape.txt", tmp_path)
        assert err and "would write to" in err, err


# ---------------------------------------------------------------------------
# Case 3: a commit message body is prose, not a write target
# ---------------------------------------------------------------------------

HEREDOC_COMMIT = (
    'git commit -m "$(cat <<\'EOF\'\n'
    'sys.executable arrived as ".../worktrees/name/../../../.venv/bin/python"\n'
    'EOF\n'
    ')"'
)


class TestCommitMessageIsProse:
    def test_a_heredoc_message_body_is_not_a_write_target(self, tmp_path):
        err = _write_gate(HEREDOC_COMMIT, tmp_path)
        assert err is None, err

    def test_a_heredoc_message_body_is_not_an_outside_read(self, tmp_path):
        err = _read_gate(HEREDOC_COMMIT, tmp_path, lock_workspace=True)
        assert err is None, err

    def test_commit_dash_F_names_no_outside_target(self, tmp_path):
        (tmp_path / ".gate").mkdir(exist_ok=True)
        (tmp_path / ".gate" / "commitmsg.txt").write_text("x\n")
        err = _write_gate("git commit -F .gate/commitmsg.txt", tmp_path)
        assert err is None, err

    def test_a_heredoc_that_really_writes_outside_still_refuses(
            self, tmp_path):
        # The commit message body is prose, but a heredoc that really
        # writes outside is not.
        err = _write_gate(
            "cat > /etc/escape.txt <<'EOF'\n"
            "../../outside\n"
            "EOF", tmp_path)
        assert err and "would write to" in err, err

    def test_a_secret_read_inside_a_substitution_still_refuses(
            self, tmp_path):
        # The prose rule must not blanket a substitution that names a
        # REAL secret path: _prose_blanked keeps $()-substitutions for
        # exactly this.
        err = _permission_gate(
            'git commit -m "$(cat ~/.delfin/c' + 'redentials.json)"',
            tmp_path)
        assert err and "secret-deny" in err, err


# ---------------------------------------------------------------------------
# The blanking never hides a path: every shape it does not fully know
# stays visible. Each command below was a hole in an earlier draft.
# ---------------------------------------------------------------------------

SECRET = "~/.delfin/c" + "redentials.json"


class TestTheBlankingNeverHidesAPath:
    def _refused(self, cmd, tmp_path):
        err = _permission_gate(cmd, tmp_path)
        assert err and "secret-deny" in err, (cmd, err)

    def test_after_e_every_free_word_is_a_file(self, tmp_path):
        self._refused(f"grep -e x {SECRET}", tmp_path)

    def test_a_pattern_file_is_a_path(self, tmp_path):
        self._refused(f"grep -f {SECRET} datei.py", tmp_path)
        self._refused(f"grep --file={SECRET} datei.py", tmp_path)
        self._refused(f"grep -rf {SECRET} .", tmp_path)

    def test_bundles_and_markers_it_cannot_read_leave_the_line_alone(
            self, tmp_path):
        self._refused(f"grep -ie x {SECRET}", tmp_path)
        self._refused(f"grep -ex {SECRET}", tmp_path)
        self._refused(f"grep -- -x {SECRET}", tmp_path)
        self._refused(f"grep --include=*.py x {SECRET}", tmp_path)

    def test_a_substitution_in_the_pattern_still_runs(self, tmp_path):
        self._refused(f'grep "$(cat {SECRET})" datei.py', tmp_path)
        self._refused(f'grep -e "$(cat {SECRET})" datei.py', tmp_path)
        self._refused(f"grep `cat {SECRET}` datei.py", tmp_path)

    def test_a_sed_script_that_can_write_or_run_stays_visible(
            self, tmp_path):
        self._refused(f"sed 's/x/y/w {SECRET}' datei.py", tmp_path)
        self._refused(f"sed '/x/w {SECRET}' datei.py", tmp_path)
        self._refused(f"sed '1e cat {SECRET}' datei.py", tmp_path)
        self._refused(f"sed -e 's/x/y/' {SECRET}", tmp_path)

    def test_a_file_argument_of_sed_is_a_path(self, tmp_path):
        self._refused(f"sed -n '/x/p' {SECRET}", tmp_path)

    def test_outside_files_still_refuse_under_a_locked_scope(self, tmp_path):
        for cmd in ("sed -n '/a/,/b/p' /etc/passwd",
                    "grep -n x /etc/passwd",
                    "grep -e /etc/hosts /etc/passwd"):
            err = _read_gate(cmd, tmp_path, lock_workspace=True)
            assert err and "blocked" in err, (cmd, err)
