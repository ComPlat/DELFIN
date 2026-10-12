"""H1B reviewer (nacht-s21) GREEN adversarial tests for self_written.verdict.

These lock in behaviours that were NOT in the builder's own test_h1b_written.py
-- the confirming direction (fires when it should, does not fire when it
should not) for spellings the module is expected to get right:

- shell-quoted script path (`python3 'x.py'`) still resolves to the written
  file and fires;
- interpreter flags before the path (`-O`, `-u`, `-B`) do not swallow the
  path (the gate's interpreter branch skips `-*` tokens);
- multiple spaces between interpreter and path still fire (shlex collapses);
- `$(cat …)` / `$(python3 …)` command substitution fires in any position
  (head and mid-command), matching the builder's read-then-run extension;
- a unicode look-alike slash (U+2215) is NOT a real path separator, so it
  must NOT accidentally resolve to the written file (differs by name, and
  the danger is a false positive in the wrong direction, not a miss);
- data-position paths (grep -f x.py, cat x.py, git diff x.py) do NOT fire:
  reading the file is not executing it.

These pass on the current phase-2 module (28001521). The RED adversarial
findings live in .gate/h1b_red_s21.py and were reported to the builder.
"""

import delfin.agent.self_written as sw


def _ws(tmp_path) -> str:
    return str(tmp_path)


def _written(tmp_path, *names) -> list:
    paths = []
    for n in names:
        p = tmp_path / n
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("echo payload\n")
        paths.append(str(p))
    return paths


def test_shell_quoted_script_path_fires(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    v = sw.verdict("python3 'x.py'", written, workspace)
    assert v.matched
    assert v.key == written[0]


def test_interp_flag_before_path_does_not_swallow(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    for cmd in ("python3 -O x.py", "python3 -u x.py", "python3 -B x.py"):
        v = sw.verdict(cmd, written, workspace)
        assert v.matched, cmd
        assert v.key == written[0], cmd


def test_multi_space_between_interp_and_path(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    v = sw.verdict("python3   x.py", written, workspace)
    assert v.matched
    assert v.key == written[0]


def test_cmdsub_cat_at_segment_head_fires(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.sh")
    v = sw.verdict("$(cat x.sh)", written, workspace)
    assert v.matched


def test_cmdsub_python_mid_command_fires(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    v = sw.verdict("echo $(python3 x.py)", written, workspace)
    assert v.matched
    assert v.key == written[0]


def test_cmdsub_cat_absolute_path_fires(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.sh")
    v = sw.verdict('sh -c "$(cat %s)"' % written[0], written, workspace)
    assert v.matched
    assert v.key == written[0]


def test_unicode_slash_lookalike_does_not_fire(tmp_path):
    """U+2215 (⁄) is not the ASCII '/', so x⁄py.py is a DIFFERENT file name.
    Whatever verdict does with it, it must NOT resolve to the written x/py.py
    (a false positive). We only pin the "does not fire" end here."""
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    v = sw.verdict("python3 x\u2215.py", written, workspace)
    assert not v.matched


def test_data_position_paths_do_not_fire(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    for cmd in ("grep -f x.py", "cat x.py", "git diff x.py", "sed -n 1,5p x.py"):
        v = sw.verdict(cmd, written, workspace)
        assert not v.matched, cmd
        assert v.hint == "", cmd


def test_script_not_in_written_paths_does_not_fire(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "mine.py")
    other = tmp_path / "other.py"
    other.write_text("print('x')\n")
    v = sw.verdict("python3 other.py", written, workspace)
    assert not v.matched
    assert v.hint == ""


def test_empty_command_or_empty_written_returns_not_matched(tmp_path):
    workspace = _ws(tmp_path)
    assert not sw.verdict("", [], workspace).matched
    assert not sw.verdict("python3 x.py", [], workspace).matched
    assert not sw.verdict("", [str(tmp_path / "x.py")], workspace).matched


# --- folded-in findings, GREEN after builder amend fe8cfa76 ---

def test_nested_sh_c_python3_fires(tmp_path):
    """sh -c 'python3 x.py' (the -c string descends with the same tokenizer)."""
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    v = sw.verdict('sh -c "python3 x.py"', written, workspace)
    assert v.matched


def test_nested_sh_c_bash_fires(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.sh")
    v = sw.verdict('sh -c "bash x.sh"', written, workspace)
    assert v.matched


def test_stdin_redir_bash_fires(tmp_path):
    """bash < x.sh feeds the file as the program -- execution."""
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.sh")
    v = sw.verdict("bash < x.sh", written, workspace)
    assert v.matched


def test_cmdsub_bash_fires(tmp_path):
    """$(bash x.sh) reads AND runs x.sh in a subshell."""
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.sh")
    v = sw.verdict('sh -c "$(bash x.sh)"', written, workspace)
    assert v.matched


def test_cmdsub_redir_read_fires(tmp_path):
    """$(< x.sh) is bash's builtin read-without-cat."""
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.sh")
    v = sw.verdict('sh -c "$(< x.sh)"', written, workspace)
    assert v.matched


def test_python_c_runpy_fires(tmp_path):
    """python3 -c 'import runpy; runpy.run_path(x)' executes the file inside -c."""
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    v = sw.verdict("python3 -c \"import runpy; runpy.run_path('x.py')\"",
                   written, workspace)
    assert v.matched


def test_cd_relative_then_python3_fires(tmp_path):
    """cd d && python3 x.py resolves x.py against d (F1 cd-then-run, relative)."""
    workspace = _ws(tmp_path)
    (tmp_path / "d").mkdir(parents=True, exist_ok=True)
    (tmp_path / "d" / "x.py").write_text("import sys\n")
    written = [str(tmp_path / "d" / "x.py")]
    v = sw.verdict("cd d && python3 x.py", written, workspace)
    assert v.matched


def test_heredoc_inline_program_not_fp(tmp_path):
    """A heredoc feeding an inline program (not a self-written file) does NOT fire."""
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    v = sw.verdict("python3 - <<'EOF'\nprint('hi')\nEOF", written, workspace)
    assert not v.matched
