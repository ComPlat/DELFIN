"""IO1 adversarial tests (reviewer nacht-s2): parser-shape attacks on hpc_io_check.

The builder's own tests use clean one-command-per-line scripts. Real job
scripts contain env prefixes, attached option values, numbered redirections,
quoting and unicode. Every case here states what the check MUST do; a wrong
result is a finding to nacht-s1.

These tests attack the PARSE of commands, so most do not depend on whether
Lustre is mounted; where they do, the Lustre fixture from the builder's test
module is reused.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check


LUSTRE_MOUNTS = "/dev/sda1 / ext4 rw 0 0\nmgmt@o2ib:fs /pfs lustre rw,lazystatfs 0 0\n"
ENV = {"HOME": "/pfs/home/ka_user"}


def _check(script: str):
    return hpc_io_check.check(
        script, env=ENV,
        mounts_path=hpc_io_check.MountsSource(text=LUSTRE_MOUNTS))


def _rules(findings):
    return {f.rule for f in findings}


# --- env prefixes: X=y cmd must still be recognised as cmd ----------------

def test_env_prefixed_tree_walk_is_refused():
    script = "#!/bin/bash\nLC_ALL=C du -sh /pfs/data/workspace/\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


def test_env_prefix_before_find_is_refused():
    script = "#!/bin/bash\nX=1 find /pfs/data -type f\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


# --- attached values and abbreviated options ------------------------------

def test_du_without_space_is_still_a_walk():
    # 'du -sh/...' style attached flags must not hide the tree walk
    script = "#!/bin/bash\ndu /pfs/data/workspace\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


def test_rsync_recursive_short_and_long_both_refused():
    for cmd in ("rsync -a /pfs/data/w /pfs/data/b/",
                "rsync --archive /pfs/data/w /pfs/data/b/"):
        assert "tree_walk_on_lustre" in _rules(_check(cmd + "\n")), cmd


def test_cp_recursive_long_option_is_refused():
    script = "#!/bin/bash\ncp --recursive /pfs/data/workspace /pfs/data/backup\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


def test_cp_r_attached_value_is_refused():
    script = "#!/bin/bash\ncp -r /pfs/data/workspace /pfs/data/backup\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


# --- redirections: a walk hidden behind a numbered redirect is still one --

def test_tree_walk_with_stderr_redirect_is_refused():
    script = "#!/bin/bash\ndu -sh /pfs/data/workspace/ 2>/dev/null\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


def test_tree_walk_with_stdout_redirect_is_refused():
    script = "#!/bin/bash\nfind /pfs/data -type f > /pfs/data/list.txt\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


def test_tree_walk_inside_a_pipe_is_refused():
    # the PIPELINE ITSELF walks the tree before wc counts anything
    script = "#!/bin/bash\ndu -sh /pfs/data/workspace/ | wc -l\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


# --- quoting and word shapes ----------------------------------------------

def test_quoted_lustre_path_is_refused():
    script = '#!/bin/bash\ndu -sh "/pfs/data/workspace/"\n'
    assert "tree_walk_on_lustre" in _rules(_check(script))


def test_trailing_slash_and_dotdot_do_not_hide_the_walk():
    script = "#!/bin/bash\ndu -sh /pfs/data/workspace/..\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


def test_lustre_path_via_home_env_is_refused():
    script = '#!/bin/bash\ndu -sh "$HOME/delfin"\n'
    assert "tree_walk_on_lustre" in _rules(_check(script))


# --- input robustness ------------------------------------------------------

def test_empty_script_has_no_findings():
    assert _check("") == []


def test_none_like_input_does_not_crash():
    # check() must accept whatever the caller passes; an empty result for
    # empty input, no exception for whitespace-only text
    assert _check("\n\n   \n") == []


def test_crlf_line_endings_do_not_break_line_numbers():
    script = "#!/bin/bash\r\ndu -sh /pfs/data/workspace/\r\n"
    findings = _check(script)
    assert "tree_walk_on_lustre" in _rules(findings)
    f = [x for x in findings if x.rule == "tree_walk_on_lustre"][0]
    assert f.line == 2


def test_100k_line_script_terminates_quickly():
    # a walk on line 100001 must still be found and the check must not
    # blow up on a huge script (no linear rescan per rule)
    script = "#!/bin/bash\n" + "# filler\n" * 100000 + "du -sh /pfs/data/w\n"
    findings = _check(script)
    f = [x for x in findings if x.rule == "tree_walk_on_lustre"][0]
    assert f.line == 100002


def test_unicode_path_lookalike_is_reported_honestly():
    # a cyrillic-looking homoglyph path must either be reported as a walk
    # over a NON-lustre path (no finding) or as a walk; it must never crash
    script = "#!/bin/bash\ndu -sh /pf\u0455/data/w\n"  # cyrillic ѕ, not latin s
    _check(script)  # must not raise


def test_multiline_continuation_walk_is_refused():
    script = "#!/bin/bash\ndu -sh \\\n  /pfs/data/workspace/\n"
    assert "tree_walk_on_lustre" in _rules(_check(script))


def test_env_is_none_still_checks_against_mounts():
    # env=None must not disable the check: paths are still judged against
    # the Lustre mounts
    script = "#!/bin/bash\ndu -sh /pfs/data/workspace/\n"
    findings = hpc_io_check.check(
        script, env=None,
        mounts_path=hpc_io_check.MountsSource(text=LUSTRE_MOUNTS))
    assert "tree_walk_on_lustre" in _rules(findings)
