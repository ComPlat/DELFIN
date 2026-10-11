"""U8: patch_targets.targets() names the files a git/patch command writes.

Field gap (wave 14): ``git apply``, ``git am``, ``patch`` and
``git checkout -- paths`` / ``git restore`` write files the write gate never
saw, so write_allow_globs and the workspace boundary did not hold there --
only a wave-level deny pattern stopped them. ``_bash_write_targets`` returned
[] for every spelling. ``patch_targets.targets(command, cwd)`` turns each
spelling into the absolute paths it would write, or the UNKNOWABLE sentinel
for a patch whose targets cannot be read from the command text (stdin, pipe,
mailbox), which the caller must refuse.

These tests drive the public path: import ``targets`` and ``UNKNOWABLE`` from
delfin.agent.patch_targets.
"""

from __future__ import annotations

from delfin.agent.patch_targets import UNKNOWABLE, targets


# --- assisted fixtures ----------------------------------------------------

def _patch(tmp_path, text: str):
    """Write *text* into a patch file under *tmp_path* and return its path."""
    p = tmp_path / "changes.diff"
    p.write_text(text, encoding="utf-8")
    return str(p)


GIT_DIFF = "\n".join([
    "diff --git a/delfin/agent/a.py b/delfin/agent/a.py",
    "index 1111111..2222222 100644",
    "--- a/delfin/agent/a.py",
    "+++ b/delfin/agent/a.py",
    "@@ -1 +1 @@",
    "-old",
    "+new",
    "diff --git a/delfin/agent/b.py b/delfin/agent/b.py",
    "index 3333333..4444444 100644",
    "--- a/delfin/agent/b.py",
    "+++ b/delfin/agent/b.py",
    "@@ -1 +1 @@",
    "-x",
    "+y",
    "",
])


# --- git apply: named patch file ------------------------------------------

def test_git_apply_named_patch_lists_b_side(tmp_path):
    patch = _patch(tmp_path, GIT_DIFF)
    assert targets(f"git apply {patch}", cwd=str(tmp_path)) == [
        f"{tmp_path}/delfin/agent/a.py",
        f"{tmp_path}/delfin/agent/b.py",
    ]


def test_git_apply_p1_applies_strip(tmp_path):
    # Default strip=1 drops the leading a//b/ prefix; the paths are unchanged
    # because the a/ b/ forms are prefixes, not directory components.
    patch = _patch(tmp_path, GIT_DIFF)
    got = targets(f"git apply -p1 {patch}", cwd=str(tmp_path))
    assert f"{tmp_path}/delfin/agent/a.py" in got


def test_git_apply_p0_keeps_a_prefix(tmp_path):
    # -p0 keeps the literal a/ b/ prefix components.
    patch = _patch(tmp_path, GIT_DIFF)
    got = targets(f"git apply -p0 {patch}", cwd=str(tmp_path))
    assert f"{tmp_path}/a/delfin/agent/a.py" in got


def test_git_apply_directory_prefixes_files_not_phantom_dir(tmp_path):
    # The --directory VALUE must prefix the patched files, not be reported as
    # a standalone write target (reviewer note, nacht-s14).
    patch = _patch(tmp_path, GIT_DIFF)
    got = targets(
        f"git apply --directory=sub/dir {patch}", cwd=str(tmp_path))
    assert f"{tmp_path}/sub/dir/delfin/agent/a.py" in got
    assert f"{tmp_path}/sub/dir" not in got


def test_git_apply_reverse_is_same_target_set(tmp_path):
    patch = _patch(tmp_path, GIT_DIFF)
    fwd = targets(f"git apply {patch}", cwd=str(tmp_path))
    rev = targets(f"git apply -R {patch}", cwd=str(tmp_path))
    assert sorted(rev) == sorted(fwd)


def test_git_apply_check_is_not_a_write(tmp_path):
    # --check only validates; it must not be treated as writing.
    patch = _patch(tmp_path, GIT_DIFF)
    assert targets(f"git apply --check {patch}", cwd=str(tmp_path)) == []


# --- git apply: unknowable patch sources -----------------------------------

def test_git_apply_stdin_is_unknowable(tmp_path):
    assert targets("git apply", cwd=str(tmp_path)) is UNKNOWABLE


def test_git_apply_dash_is_unknowable(tmp_path):
    assert targets("git apply -", cwd=str(tmp_path)) is UNKNOWABLE


def test_git_apply_piped_is_unknowable(tmp_path):
    assert targets("git show abc | git apply", cwd=str(tmp_path)) \
        is UNKNOWABLE


def test_git_apply_read_redirect_is_unknowable(tmp_path):
    assert targets("git apply < changes.diff", cwd=str(tmp_path)) \
        is UNKNOWABLE


def test_git_apply_missing_patch_file_is_unknowable(tmp_path):
    # The patch file is named but absent: we cannot read what it would write.
    assert targets("git apply nope.diff", cwd=str(tmp_path)) is UNKNOWABLE


def test_git_am_is_unknowable(tmp_path):
    assert targets("git am 0001-fix.patch", cwd=str(tmp_path)) is UNKNOWABLE


# --- patch ----------------------------------------------------------------

def test_gnupatch_stdin_is_unknowable(tmp_path):
    assert targets("patch -p1 < fix.diff", cwd=str(tmp_path)) is UNKNOWABLE


def test_gnupatch_i_reads_header(tmp_path):
    patch = _patch(tmp_path, "+++ b/var/www/f.py\n@@ -1 +1 @@\n-x\n+y\n")
    got = targets(f"patch -i {patch}", cwd=str(tmp_path))
    assert f"{tmp_path}/var/www/f.py" in got


def test_gnupatch_positional_origplus_patchfile(tmp_path):
    # patch orig patchfile: the SECOND positional is the patch reading source,
    # the first the original file.
    patch = _patch(tmp_path, "+++ b/var/www/f.py\n@@ -1 +1 @@\n-x\n+y\n")
    got = targets(f"patch fix.txt {patch}", cwd=str(tmp_path))
    assert f"{tmp_path}/var/www/f.py" in got


def test_gnupatch_dry_run_is_not_a_write(tmp_path):
    assert targets("patch --dry-run -i x.diff", cwd=str(tmp_path)) == []


# --- checkout / restore pathspecs ------------------------------------------

def test_git_checkout_double_dash_paths(tmp_path):
    got = targets("git checkout other -- delfin/x.py tests/y.py",
                  cwd=str(tmp_path))
    assert got == [f"{tmp_path}/delfin/x.py", f"{tmp_path}/tests/y.py"]


def test_git_checkout_branchless_path(tmp_path):
    # git checkout <path> (path-mode checkout) writes the path.
    got = targets("git checkout delfin/x.py", cwd=str(tmp_path))
    assert f"{tmp_path}/delfin/x.py" in got


def test_git_restore_source(tmp_path):
    got = targets("git restore --source HEAD delfin/x.py",
                  cwd=str(tmp_path))
    assert got == [f"{tmp_path}/delfin/x.py"]


def test_git_restore_staged(tmp_path):
    got = targets("git restore --staged -- delfin/x.py",
                  cwd=str(tmp_path))
    assert got == [f"{tmp_path}/delfin/x.py"]


# --- non-write / read-only git must yield nothing --------------------------

def test_git_read_only_commands_yield_nothing(tmp_path):
    assert targets("git log --oneline", cwd=str(tmp_path)) == []
    assert targets("git status", cwd=str(tmp_path)) == []
    assert targets("git diff --stat", cwd=str(tmp_path)) == []


def test_git_commit_is_not_a_worktree_write(tmp_path):
    # A commit writes git objects, not working-tree files; and a `>` inside a
    # quoted -m message must never surface here (the phase-3 fix handles the
    # redirection; this module simply does not see git commit as a file write).
    assert targets('git commit -m "use -> apply_proposal now"',
                   cwd=str(tmp_path)) == []


def test_non_patch_command_yields_nothing(tmp_path):
    # Left for the existing _bash_write_targets; this module is scoped to
    # git/patch file-writes.
    assert targets("echo hi > out.txt", cwd=str(tmp_path)) == []


def test_env_var_prefix_is_skipped(tmp_path):
    patch = _patch(tmp_path, GIT_DIFF)
    got = targets(f"GIT_PAGER=cat git apply {patch}", cwd=str(tmp_path))
    assert f"{tmp_path}/delfin/agent/a.py" in got
