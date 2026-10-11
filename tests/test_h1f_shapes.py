"""RED control for package H1F: read-only command shapes approved by rule.

delfin.agent.read_only_shapes is a NEW module (does not exist on main). These
tests exercise its public API. Against the unchanged tree the import fails and
every case is RED — that is the control. After the module lands, every case
must be GREEN and every shape is tested in BOTH directions (fires when its
safe read-only form matches, does not fire when a write/execute option, an
out-of-workspace path, a compound ';', a substitution, a redirect or a
disabled switch is present).
"""

from delfin.agent.read_only_shapes import matches, matches_with, SHAPE_NAMES

WORKSPACE = "/pfs/data6/home/ka/example/delfin"


def test_import_lands():
    assert callable(matches)


# --- shape: git_show_path --------------------------------------------------



def test_git_show_path_outside_workspace_does_not_fire():
    assert matches("git show HEAD:../../etc/passwd", WORKSPACE) == ""


def test_git_show_path_output_option_does_not_fire():
    assert matches("git show HEAD:delfin/api.py --output=/tmp/x", WORKSPACE) == ""


def test_git_show_path_config_evil_pager_does_not_fire():
    assert (
        matches("git -c core.pager=x log", WORKSPACE) == ""
    )


# --- git_show bare-arg containment (reviewer-s27) ---------------------------

def test_git_show_bare_absolute_path_does_not_fire():
    # a bare absolute path arg (no colon) must not ride along unvalidated
    assert matches("git show HEAD:config.py /etc/passwd", WORKSPACE) == ""




def test_git_show_bare_parent_escape_does_not_fire():
    assert matches("git show HEAD ../config.py", WORKSPACE) == ""


# --- shape: git_grep -------------------------------------------------------



def test_git_grep_open_in_pager_does_not_fire():
    assert (
        matches("git grep -O sh TODO delfin/", WORKSPACE) == ""
    )


def test_git_grep_textconv_does_not_fire():
    assert (
        matches("git grep --textconv TODO delfin/", WORKSPACE) == ""
    )


# --- git_grep pattern-file reads (reviewer-s28 finding) ---------------------

def test_git_grep_dash_f_absolute_file_does_not_fire():
    # -f FILE reads the pattern list from FILE, which may be outside the
    # workspace; the shape must not approve the read.
    assert matches("git grep -f /etc/shadow TODO", WORKSPACE) == ""


def test_git_grep_dash_f_attached_does_not_fire():
    assert matches("git grep -f/etc/passwd TODO", WORKSPACE) == ""


def test_git_grep_long_file_option_does_not_fire():
    assert matches("git grep --file /etc/shadow TODO", WORKSPACE) == ""


def test_git_grep_long_file_equals_does_not_fire():
    assert matches("git grep --file=/etc/passwd TODO", WORKSPACE) == ""




# --- shape: git_apply_check ------------------------------------------------

def test_git_apply_check_fires():
    assert (
        matches("git apply --check .gate/h1f.patch", WORKSPACE) == "git_apply_check"
    )


def test_git_apply_without_check_does_not_fire():
    assert matches("git apply .gate/h1f.patch", WORKSPACE) == ""


def test_git_apply_check_index_does_not_fire():
    assert matches("git apply --check --index .gate/h1f.patch", WORKSPACE) == ""


def test_git_apply_check_write_option_does_not_fire():
    assert (
        matches("git apply --check .gate/h1f.patch --output=/etc/out", WORKSPACE)
        == ""
    )


def test_git_apply_check_outside_workspace_does_not_fire():
    assert (
        matches("git apply --check /etc/hostname", WORKSPACE) == ""
    )


def test_git_apply_check_unknown_option_does_not_fire():
    # --tool is not a git apply option (that is difftool/mergetool); an
    # unknown option means we are not looking at a known read-only form.
    assert (
        matches("git apply --check --tool=x .gate/h1f.patch", WORKSPACE) == ""
    )


# --- shape: git_hash_object ------------------------------------------------

def test_git_hash_object_fires():
    assert (
        matches("git hash-object delfin/api.py", WORKSPACE) == "git_hash_object"
    )


def test_git_hash_object_write_does_not_fire():
    assert (
        matches("git hash-object -w delfin/api.py", WORKSPACE) == ""
    )


def test_git_hash_object_outside_workspace_does_not_fire():
    assert matches("git hash-object /etc/bashrc", WORKSPACE) == ""


# --- shape: git_switch_pr --------------------------------------------------

def test_git_switch_pr_fires():
    assert (
        matches("git switch -c pr/h1f origin/main", WORKSPACE)
        == "git_switch_pr"
    )


def test_git_switch_existing_does_not_fire():
    assert matches("git switch main", WORKSPACE) == ""


def test_git_switch_pr_name_with_dotdot_does_not_fire():
    # the branch NAME after pr/ is bounded: `..` must be rejected
    assert matches("git switch -c pr/../evil origin/main", WORKSPACE) == ""


def test_git_switch_pr_bare_dotdot_does_not_fire():
    # reviewer-s28 probe spelling: a lone ".." name is a path special
    assert matches("git switch -c pr/.. origin/main", WORKSPACE) == ""


def test_git_switch_pr_bare_dot_does_not_fire():
    # reviewer-s28 probe spelling: a lone "." name is a path special
    assert matches("git switch -c pr/. origin/main", WORKSPACE) == ""


def test_git_switch_pr_name_with_inner_dotdot_does_not_fire():
    assert matches("git switch -c pr/a..b origin/main", WORKSPACE) == ""


def test_git_switch_pr_name_leading_dash_does_not_fire():
    assert matches("git switch -c pr/-dash origin/main", WORKSPACE) == ""


def test_git_switch_pr_extra_token_does_not_fire():
    assert matches("git switch -c pr/x origin/main extra", WORKSPACE) == ""


def test_git_switch_pr_valid_named_branch_fires():
    assert (
        matches("git switch -c pr/h1f-name.2 origin/main", WORKSPACE)
        == "git_switch_pr"
    )


# --- shape: gh_version -----------------------------------------------------

def test_gh_version_fires():
    assert matches("gh --version", WORKSPACE) == "gh_version"


def test_gh_action_does_not_fire():
    assert matches("gh pr create", WORKSPACE) == ""


# --- shape: gate_probe was removed (operator finding 1) ---------------------
# A program run is never a read-only shape: a command whose program is named
# `gate` (or a workspace file named `gate-tools/gate`) must still be asked
# about, not silently approved by rule.

def test_gate_named_program_is_not_a_shape():
    assert matches("gate tests/test_zz_h1f_probe.py", WORKSPACE) == ""


def test_workspace_gate_tools_gate_file_is_not_a_shape():
    # a session could write ./tools/gate-tools/gate itself and then have that
    # program run -- that must still be asked about, never auto-approved
    assert matches("tools/gate-tools/gate tests/test_zz_h1f_probe.py", WORKSPACE) == ""


# --- compounds / structure are never a shape -------------------------------

def test_semicolon_compound_does_not_fire():
    assert matches("git show HEAD:delfin/api.py; echo hi", WORKSPACE) == ""


def test_command_substitution_does_not_fire():
    assert matches("git show HEAD:$(ls)", WORKSPACE) == ""


def test_redirect_does_not_fire():
    assert matches("git show HEAD:delfin/api.py > /tmp/x", WORKSPACE) == ""


def test_source_does_not_fire():
    assert matches("source .envrc", WORKSPACE) == ""


# --- per-shape enablement --------------------------------------------------

def test_matches_with_disabled_shape_does_not_fire():
    enabled = {s for s in SHAPE_NAMES if s != "git_show_path"}
    assert (
        matches_with("git show HEAD:delfin/api.py", WORKSPACE, enabled) == ""
    )






