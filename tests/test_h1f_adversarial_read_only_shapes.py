"""Adversarial tests for H1F read_only_shapes — s28 reviewer (wave 15).

These go beyond tests/test_h1f_shapes.py (the builder's control, which covers
the *named* task gaps). Here I attack the ALLOW-list surface a shape matcher
is most likely to widen: attached/combined short values, `-f FILE` reading
its argument, git-apply index writes, `switch` extra flags, rev:path escapes,
and option-look-alike punctuation. Every case asserts the shape does NOT fire
(green == no widening). A red here is a finding for nacht-s24.
"""

import pytest

from delfin.agent.read_only_shapes import SHAPE_NAMES, matches, matches_with

WS = "/pfs/data6/home/ka/example/delfin"


# --- git grep: reading a pattern from an arbitrary FILE ---------------------

def test_grep_dash_f_absolute_file_does_not_fire():
    # git grep -f FILE reads patterns from FILE -> an arbitrary FILE read that
    # is NOT inside the workspace. Must not be approved.
    assert matches("git grep -f /etc/shadow TODO", WS) == ""


def test_grep_dash_f_writes_none_but_reads_outside_does_not_fire():
    # attached value form -f/relative
    assert matches("git grep TODO -f ../../../etc/passwd", WS) == ""


def test_grep_dash_f_attached_does_not_fire():
    assert matches("git grep -f/etc/passwd TODO", WS) == ""


# --- git grep: combined / attached short options ----------------------------

def test_grep_combined_short_o_does_not_fire():
    # -o combined with another short value flag must not pass as pure grep
    assert matches("git grep -e x -o", WS) == ""


def test_grep_count_c_does_not_fire():
    # `git grep -c` is read-only for grep BUT -c is a write/exec global for
    # git; over-strict is safe, so assert "".
    assert matches("git grep -c TODO", WS) == ""


def test_grep_open_files_attached_value_does_not_fire():
    assert matches("git grep -Ovim TODO", WS) == ""


# --- git apply: index / cached writes ---------------------------------------

def test_apply_check_no_fire_with_worktree_variant():
    assert matches("git apply --check -3 x.patch", WS) == ""


def test_apply_check_tool_does_not_fire():
    assert matches("git apply --check --tool=x x.patch", WS) == ""


# --- git switch: extra flags after the branch --------------------------------

def test_switch_extra_flag_does_not_fire():
    assert matches("git switch -c pr/x origin/main --orphan", WS) == ""


def test_switch_double_dash_does_not_fire():
    assert matches("git switch -c pr/x -- origin/main", WS) == ""


# --- git show: rev:path escapes ---------------------------------------------

def test_show_absolute_path_inside_colon_does_not_fire():
    # path must be repo-relative; an absolute-looking path must not clear
    assert matches("git show HEAD:/etc/passwd", WS) == ""


def test_show_empty_revkolon_reads_index():
    # `git show :path` reads a path from the index -- read-only and inside
    # the workspace. NOT a finding: withdrawn by s28 (my earlier assumption
    # that it must be rejected was wrong). Assert it fires (is a valid shape).
    assert matches("git show :delfin/api.py", WS) == "git_show_path"


def test_show_colon_path_dotdot_does_not_fire():
    assert matches("git show HEAD:../delfin/api.py", WS) == ""


# --- generic: option-look-alike punctuation ----------------------------------

def test_semicolon_in_quoted_pattern_still_does_not_fire():
    # the matcher is conservative: even a quoted ';' is not a shape
    assert matches("git grep 'foo;bar'", WS) == ""


def test_error_redirect_look_alike_does_not_fire():
    assert matches("git show 'HEAD:x>&2'", WS) == ""


def test_compound_on_other_shape_does_not_fire():
    assert matches("gh --version; id", WS) == ""
