"""H1B reviewer adversarial tests for self_written.verdict/executed_script_paths.

The builder's tests (test_h1b_written.py) cover the happy paths. These probe
the boundary of what the tokenizer treats as EXECUTION of a file the session
wrote.

QS operator rule (this run): probes use the REAL shape, not hand-made dicts.
``written_paths`` is therefore taken from real change_journal records: each
case writes the file, records it through ``change_journal.record_change``
(the real writer the write executor calls), reads the records back verbatim
with ``change_journal.list_changes``, and passes those real ``path`` values
to ``verdict``. A unique session id keeps test journals out of real ones.

Two cases here are RED probes (they assert the correctly-desired behaviour
and FAIL against the current implementation): ``test_pipe_segment_exec`` and
``test_bare_name_via_env``. They are the findings the reviewer reports to the
builder, not fixes in this file.
"""

import uuid

import delfin.agent.self_written as sw
from delfin.agent import change_journal


def _ws(tmp_path):
    return str(tmp_path)


def _record_written(tmp_path, name, text="print('x')\n"):
    """Real-shape written_paths: write file, journal it via the real writer,
    read the real records back, return the real list of written paths."""
    sid = "h1b-adv-" + uuid.uuid4().hex[:12]
    p = tmp_path / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    change_journal.record_change(
        sid,
        tool="write_file",
        path=str(p),
        old_text=None,  # created
        new_text=text,
    )
    records = change_journal.list_changes(sid)
    assert records, f"no real record for {name!r} under sid {sid}"
    # Return the FULL record dicts exactly as the P3 gate's _journal_changes
    # wraps them ({path, ts, created}), so verdict() is probed with the true
    # gate shape (builder's record-shape amendment accepts these).
    return records


# --------------------------------------------------------------------------
# PER-MATCH: interpreters / wrappers / path-position run a self-written file.
# --------------------------------------------------------------------------

def test_interpreter_abs_exec(tmp_path):
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "x.py")
    # use the REAL absolute path the journal recorded, run under python3
    real_abs = str(tmp_path / "x.py")
    v = sw.verdict("python3 %s" % real_abs, written, ws)
    assert v.matched, "python3 <abs self-written path> must be caught"


def test_interpreter_relative_exec_matches_real_record(tmp_path):
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "scripts/x.py")
    v = sw.verdict("python3 scripts/x.py", written, ws)
    assert v.matched and v.recognised.endswith("scripts/x.py")


def test_timeout_env_nohup_wrappers(tmp_path):
    """('timeout','env','nohup') run the following command; a ./path target
    in execution position is caught."""
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "x.py")
    for cmd in (
        "timeout 30 ./x.py",
        "env VAR=1 ./x.py",
        "nohup ./x.py > /dev/null 2>&1",
    ):
        v = sw.verdict(cmd, written, ws)
        assert v.matched, f"{cmd!r} (self-written ./x.py) must be caught"


def test_dot_slash_exec(tmp_path):
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "gen.py")
    v = sw.verdict("./gen.py", written, ws)
    assert v.matched and v.key and v.hint, "./ self-written script is execution"


def test_sh_c_dollar_cat_exec(tmp_path):
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "payload.py")
    v = sw.verdict('sh -c "$(cat payload.py)"', written, ws)
    assert v.matched, "sh -c \"$(cat self-written)\" reads THEN runs; it is execution"


def test_source_and_dot_exec(tmp_path):
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "setup.sh")
    for cmd in ("source setup.sh", ". setup.sh"):
        v = sw.verdict(cmd, written, ws)
        assert v.matched, f"{cmd!r} runs the self-written file in-shell"


def test_cd_dir_then_python(tmp_path):
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "d/x.py")
    v = sw.verdict("cd d && python3 x.py", written, ws)
    assert v.matched, "cd to the dir then python3 the self-written file"


def test_heredoc_self_written_not_exec(tmp_path):
    """A heredoc body is not a path; verdict must not crash and not match."""
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "x.py")
    v = sw.verdict("python3 - <<EOF\nx.py\nEOF", written, ws)
    assert not v.matched, "heredoc body is data, not an executed script path"


def test_python_dash_m_not_a_script(tmp_path):
    """python3 -m <module> is module execution, not an executed script path."""
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "x.py")
    v = sw.verdict("python3 -m pytest x.py", written, ws)
    assert not v.matched, "-m runs a module; x.py is an argument, not executed"


def test_not_self_written_not_matched(tmp_path):
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "x.py")
    v = sw.verdict("python3 /other/y.py", written, ws)
    assert not v.matched, "a script the session did NOT write is judged as before"


def test_reading_not_execution(tmp_path):
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "x.py")
    for cmd in ("cat x.py", "git diff x.py", "sed -n 1,5p x.py", "ls x.py"):
        v = sw.verdict(cmd, written, ws)
        assert not v.matched, f"{cmd!r} READS, not executes"


def test_second_spelling_same_key(tmp_path):
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "go.sh")
    a = sw.verdict("bash go.sh", written, ws)
    b = sw.verdict("sh go.sh", written, ws)
    assert a.matched and b.matched and a.key == b.key, (
        "same file under a second spelling must yield the SAME key "
        "(so the gate already-refused tier can repeat-stop)"
    )


def test_relative_write_abs_run_same_key(tmp_path):
    """Written 'scripts/x.py' (real record is absolute under tmp_path), run
    './scripts/x.py' — both absolutize to the same key."""
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "scripts/x.py")
    v = sw.verdict("./scripts/x.py", written, ws)
    assert v.matched and v.key.endswith("scripts/x.py"), (
        "a ./ run of a written file must match via absolutisation"
    )


# --------------------------------------------------------------------------
# RED PROBES — correctly-desired behaviour that FAILS today (the findings).
# These are reported to the builder; they are not fixed in this file.
# --------------------------------------------------------------------------

def test_pipe_segment_exec(tmp_path):
    """RED: a ./script on the right of a pipe is execution and must be caught."""
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "gen.py")
    v = sw.verdict("echo go | ./gen.py", written, ws)
    assert v.matched, "a ./script on the right of a pipe is execution"


def test_bare_name_via_env(tmp_path):
    """RED: a bare self-written executable run via ``env`` (no ./ or extension)
    is execution of a written file and must be caught."""
    ws = _ws(tmp_path)
    written = _record_written(tmp_path, "runme", "#! /bin/sh\nexit 0\n")
    v = sw.verdict("env VAR=1 runme", written, ws)
    assert v.matched, "env-prefixed bare self-written exec must be caught"
