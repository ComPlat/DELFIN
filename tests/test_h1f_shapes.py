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

def test_git_show_path_fires():
    assert matches("git show HEAD:delfin/api.py", WORKSPACE) == "git_show_path"


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


def test_git_show_bare_relative_rev_inside_fires():
    # a bare relative rev/ref token resolves inside the workspace and is fine
    assert matches("git show v1.2.3 delfin/api.py", WORKSPACE) == "git_show_path"


def test_git_show_bare_parent_escape_does_not_fire():
    assert matches("git show HEAD ../config.py", WORKSPACE) == ""


# --- shape: git_grep -------------------------------------------------------

def test_git_grep_fires():
    assert matches("git grep -n TODO delfin/", WORKSPACE) == "git_grep"


def test_git_grep_open_in_pager_does_not_fire():
    assert (
        matches("git grep -O sh TODO delfin/", WORKSPACE) == ""
    )


def test_git_grep_textconv_does_not_fire():
    assert (
        matches("git grep --textconv TODO delfin/", WORKSPACE) == ""
    )


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


# --- shape: gh_version -----------------------------------------------------

def test_gh_version_fires():
    assert matches("gh --version", WORKSPACE) == "gh_version"


def test_gh_action_does_not_fire():
    assert matches("gh pr create", WORKSPACE) == ""


# --- shape: gate_probe -----------------------------------------------------

def test_gate_probe_fires():
    assert (
        matches("gate tests/test_zz_h1f_probe.py", WORKSPACE) == "gate_probe"
    )


def test_gate_probe_outside_workspace_does_not_fire():
    assert matches("gate /tmp/x.py", WORKSPACE) == ""


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


def test_matches_with_enabled_shape_fires():
    enabled = set(SHAPE_NAMES)
    assert (
        matches_with("git show HEAD:delfin/api.py", WORKSPACE, enabled)
        == "git_show_path"
    )


def test_shape_names_are_known():
    assert "git_show_path" in SHAPE_NAMES
    assert "git_grep" in SHAPE_NAMES
    assert "git_apply_check" in SHAPE_NAMES
    assert "git_hash_object" in SHAPE_NAMES
    assert "git_switch_pr" in SHAPE_NAMES
    assert "gh_version" in SHAPE_NAMES
    assert "gate_probe" in SHAPE_NAMES


def test_no_write_option_in_head_matches():
    # every shape's safety depends on the write/exec-option veto; the set
    # must include the task's named options
    from delfin.agent.read_only_shapes import _WRITE_OR_EXEC_OPTS
    for o in ("--output", "-o", "--ext-diff", "--textconv", "-c",
              "--exec-path", "--upload-pack", "--open-files-in-pager"):
        assert o in _WRITE_OR_EXEC_OPTS
