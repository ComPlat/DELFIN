# Package H1F — adversarial tests (nacht-s27, reviewer).
#
# These target `delfin.agent.read_only_shapes.matches(command, workspace) ->
# shape name or ""`. The contract (package task): a PARSED command (shlex, one
# simple command, no operator/substitution/newline/redirect), a closed list of
# shapes, every path argument inside the workspace, none of the writing options.
# Anything not a clean read-only shape must return "" (so the normal gate,
# deny-list, secret scan and read-gate keep their say).
#
# Test discipline (ALLOW-lists, not deny-lists): for every shape, each escape
# that could turn the read into a write, an execution or an out-of-workspace
# read MUST return "". The safe canonical form MUST return a non-empty name.
#
# Written against the contract before the builder's module merged, so this file
# is RED (uses read_only_shapes which does not exist yet) until phase 2 lands.

import pytest


def _matches(command, workspace):
    from delfin.agent import read_only_shapes
    return read_only_shapes.matches(command, workspace)


WS = "/srv/delfin/ws"


class TestGitShow:
    """git show <rev>:<path> — safe; path escapes / --output / -o / exec veto."""


    def test_path_outside_workspace(self):
        # absolute path outside the workspace -> must not be a shape
        assert _matches("git show HEAD:config.py /etc/passwd", WS) == ""

    def test_relative_slash_escape(self):
        # ../../ escapes upward even when a leading token is inside
        assert _matches("git show HEAD:../../etc/passwd %s/readme.md" % WS, WS) == ""

    def test_output_option_writes(self):
        assert _matches("git show HEAD:readme.md --output=/tmp/x %s" % WS, WS) == ""

    def test_short_o_writes(self):
        assert _matches("git show HEAD:readme.md -o /tmp/x %s" % WS, WS) == ""

    def test_ext_diff_runs(self):
        assert _matches("git show --ext-diff HEAD:readme.md %s" % WS, WS) == ""

    def test_open_files_in_pager_runs(self):
        assert _matches("git show --open-files-in-pager=x HEAD %s" % WS, WS) == ""

    def test_command_operator_rejected(self):
        # a `;` makes it two commands — must not be a single shape
        assert _matches("git show HEAD:readme.md; rm -rf / %s" % WS, WS) == ""

    def test_substitution_rejected(self):
        assert _matches("git show HEAD:$(id) %s" % WS, WS) == ""


class TestGitGrep:
    """git grep — safe; -O / --open-files-in-pager runs a pager/program."""


    def test_open_files_in_pager_runs(self):
        assert _matches("git grep -O sh TODO %s" % WS, WS) == ""

    def test_ext_diff_runs(self):
        assert _matches("git grep --ext-diff TODO %s" % WS, WS) == ""



class TestGitApplyCheck:
    """git apply --check — safe dry-run; --index / --cached mutate the index."""

    def test_safe(self):
        assert _matches("git apply --check %s/patch.patch" % WS, WS)

    def test_index_mutates(self):
        assert _matches("git apply --check --index %s/patch.patch" % WS, WS) == ""

    def test_cached_mutates(self):
        assert _matches("git apply --check --cached %s/patch.patch" % WS, WS) == ""

    def test_apply_without_check_writes(self):
        assert _matches("git apply %s/patch.patch" % WS, WS) == ""


class TestGitHashObject:
    """git hash-object <rel> — safe unless -w writes the object DB."""

    def test_safe(self):
        assert _matches("git hash-object %s/readme.md" % WS, WS)


    def test_write_flag_rejected(self):
        assert _matches("git hash-object -w %s/readme.md" % WS, WS) == ""


class TestGhVersion:
    """gh --version — pure; any other gh subcommand is a write/action."""

    def test_safe(self):
        # gh --version takes no path argument — the pure two-token command
        assert _matches("gh --version", WS)

    def test_gh_not_just_version_rejected(self):
        assert _matches("gh pr create --base main", WS) == ""

    def test_gh_version_extra_arg_rejected(self):
        # a trailing junk arg is not the shape — must stay ""
        assert _matches("gh --version %s" % WS, WS) == ""


class TestGitSwitch:
    """git switch -c pr/<name> origin/main — the approved creation shape."""

    def test_safe_branch_create(self):
        # git switch -c pr/<name> origin/main — takes no workspace path arg
        assert _matches("git switch -c pr/fix origin/main", WS)

    def test_trailing_arg_rejected(self):
        # a junk trailing argument is not the exact shape -> must stay ""
        assert _matches("git switch -c pr/fix origin/main %s" % WS, WS) == ""

    def test_core_pager_global_runs(self):
        assert _matches("git -c core.pager=x log %s" % WS, WS) == ""


class TestCrossCutting:
    """Cases that must be "" for EVERY shape: operators, redirects, env-pre."""

    @pytest.mark.parametrize("cmd", [
        "git show HEAD:readme.md %s | sh" % WS,            # pipe
        "git grep TODO %s > /tmp/o" % WS,                 # redirect
        "git show HEAD:readme.md %s &" % WS,              # background
        "X=y git show HEAD:readme.md %s" % WS,            # env prefix
        "git show HEAD:readme.md %s && true" % WS,        # && operator
        "git show HEAD:readme.md %s || true" % WS,        # || operator
    ])
    def test_operator_redirect_env_all_rejected(self, cmd):
        assert _matches(cmd, WS) == ""
