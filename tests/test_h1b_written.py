"""H1B: a command that EXECUTES a file the session itself wrote.

self_written.verdict(command, written_paths, workspace) is a pure
detect-and-hint: it answers whether `command` runs one of `written_paths`
and, if so, returns a hint naming the sanctioned route plus a *key* (the
normalised absolute path of the executed self-written file) so the gate can
record the repeat-stop on the SAME store it already uses (denied_actions,
key "self_written:<path>"). It never allows anything: matched=False, or
matched=True with key + hint.

Controls:
- a written script executed via interpreter (`python3 scripts/x.py`) -> matched
- via `bash f.sh`, `./gen.py`, `source setup.sh`, `. setup.sh` -> matched
- inline execution `sh -c "$(cat payload.py)"` -> matched
- a script the session did NOT write -> not matched (judged as before)
- READING the file (`cat`, `git diff`, `sed -n`) -> not execution -> not matched
- the SAME file under a second spelling -> the SAME returning key (so the
  gate's existing already-refused tier can stop the turn and ask the user).
"""

import delfin.agent.self_written as sw
import os

from pathlib import Path


def _ws(tmp_path) -> str:
    return str(tmp_path)


def _written(tmp_path, *names) -> list:
    paths = []
    for n in names:
        p = tmp_path / n
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("print('payload')\n")
        paths.append(str(p))
    return paths


def test_interpreter_exec_written_script_matches(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "scripts/x.py")
    v = sw.verdict("python3 scripts/x.py", written, workspace)
    assert v.matched
    assert v.key == written[0]
    assert v.hint


def test_bash_script_matches(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "run.sh")
    v = sw.verdict("bash run.sh", written, workspace)
    assert v.matched
    assert v.key == written[0]


def test_dot_slash_matches(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "gen.py")
    v = sw.verdict("./gen.py", written, workspace)
    assert v.matched
    assert v.key == written[0]


def test_source_matches(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "env.sh")
    assert sw.verdict("source env.sh", written, workspace).matched
    assert sw.verdict(". env.sh", written, workspace).matched


def test_inline_cat_execution_matches(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "payload.py")
    v = sw.verdict('sh -c "$(cat payload.py)"', written, workspace)
    assert v.matched


def test_script_not_written_does_not_match(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "mine.py")
    # a different file the session never wrote
    other = tmp_path / "other.py"
    other.write_text("print('x')\n")
    v = sw.verdict("python3 other.py", written, workspace)
    assert not v.matched
    assert v.hint == ""


def test_reading_is_not_execution(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "x.py")
    for cmd in ("cat x.py", "git diff x.py", "sed -n 1,5p x.py"):
        v = sw.verdict(cmd, written, workspace)
        assert not v.matched, cmd
        assert v.hint == "", cmd


def test_second_spelling_same_key(tmp_path):
    workspace = _ws(tmp_path)
    written = _written(tmp_path, "go.sh")
    first = sw.verdict("bash go.sh", written, workspace)
    second = sw.verdict("sh go.sh", written, workspace)
    assert first.matched and second.matched
    assert first.key == second.key == written[0]


def test_relative_written_path_matches_abs_run(tmp_path):
    """The change journal stores each written path in the CALLER'S spelling
    (change_journal.py:398, Path(str(path))) — a relative write stays
    relative in list_changes. verdict() must absolutize the WRITTEN side
    against the workspace too, or a script written "x.py" and run
    "./x.py" is missed — the exact case the four waves reported.
    """
    workspace = _ws(tmp_path)
    # written_paths as the journal echoes them: RELATIVE spellings "x.py"
    # and "scripts/setup.sh" (no leading workspace, no "./").
    written = ["x.py", "scripts/setup.sh"]
    assert sw.verdict("python3 ./x.py -q", written, workspace).matched
    assert sw.verdict("bash scripts/setup.sh", written, workspace).matched
    # and the repeat-stop key is the abs path, the same as an abs write:
    key_rel = sw.verdict("python3 x.py", written, workspace).key
    key_abs = sw.verdict("bash scripts/setup.sh", written, workspace).key
    assert key_rel == _ws(tmp_path) + "/x.py"
    assert key_abs == _ws(tmp_path) + "/scripts/setup.sh"

def test_real_journal_record_shape_is_accepted(tmp_path, monkeypatch):
    """The change journal stores records as {path, ts, created} (list_changes
    echoes them) and _journal_changes(sid) passes those to the gate, which
    hands verdict() the raw records. A record must be matched by its "path"
    key, not crash. We drive the REAL writer (change_journal.record_change is
    what the executor's write path calls) into a redirected home so no real
    user state is touched (state_paths rule: never touch the real home)."""
    import delfin.agent.change_journal as cj

    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    workspace = _ws(tmp_path)
    import os
    x = os.path.join(workspace, "x.py")   # absolute, as record_change stores it
    (tmp_path / "x.py").write_text("print(1)")
    cj.record_change("sid-h1b", tool="write_file", path=x,
                     old_text=None, new_text="print(1)")
    # the shape _journal_changes feeds the gate:
    records = [{"path": r["path"], "ts": r["ts"], "created": True}
               for r in (cj.list_changes("sid-h1b") or []) if r.get("path")]
    assert any(r["path"] == x for r in records)          # real record present
    v = sw.verdict("python3 x.py", records, workspace)   # real records, no crash
    assert v.matched
    assert v.key == x


def test_f1_cd_then_run_resolves_against_cd(tmp_path):
    """cd d && python3 x.py — x.py resolves under d/, not the workspace root.
    (reviewer s19 finding 1 / operator probe)."""
    workspace = _ws(tmp_path)
    (tmp_path / "d").mkdir()
    written = [os.path.join(workspace, "d", "x.py")]
    assert sw.verdict("cd d && python3 x.py", written, workspace).matched
    assert sw.verdict("cd d; python3 x.py", written, workspace).matched


def test_f2_right_of_pipe_and_is_not_dropped(tmp_path):
    """echo go | ./gen.py / cd d && python3 x.py — a leading-space segment
    after |/&&/; is kept (reviewer s19 finding 2)."""
    workspace = _ws(tmp_path)
    written = [os.path.join(workspace, "gen.py")]
    assert sw.verdict("echo go | ./gen.py", written, workspace).matched


def test_f4_minus_c_descent(tmp_path):
    """sh -c "python3 x.py" descends into the -c string (s21 finding 1)."""
    workspace = _ws(tmp_path)
    written = [os.path.join(workspace, "x.py")]
    assert sw.verdict('sh -c "python3 x.py"', written, workspace).matched


def test_f4b_stdin_redirect(tmp_path):
    """python3 < x.sh feeds the file to the interpreter (s21 finding 2)."""
    workspace = _ws(tmp_path)
    written = [os.path.join(workspace, "x.sh")]
    assert sw.verdict("bash < x.sh", written, workspace).matched
    assert sw.verdict("python3 < x.sh", written, workspace).matched


def test_f5_cmdsub_runs_reader(tmp_path):
    """$(bash x.sh) and $(< x.sh) read-and-run (s21 findings 3-4)."""
    workspace = _ws(tmp_path)
    written = [os.path.join(workspace, "x.sh")]
    assert sw.verdict('sh -c "$(bash x.sh)"', written, workspace).matched
    assert sw.verdict('sh -c "$(< x.sh)"', written, workspace).matched
