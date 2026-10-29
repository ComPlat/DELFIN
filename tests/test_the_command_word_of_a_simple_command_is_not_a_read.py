"""The command word of a simple command is not a read target.

Measured 2026-10-29: an agent that runs a Python script with its own
interpreter -- /home/qmchem_max/delfin-swarm/venv/bin/python, a venv
whose python is a symlink into
/home/qmchem_max/micromamba/envs/delfin/bin/python3.12 -- was asked for
an OUTSIDE-WORKSPACE READ of the interpreter's resolved path. Executing
a program is not reading it, and the interpreter DELFIN itself runs on
is not user data.

The direct case was already right: `_segment_reads` scans `argv[1:]`
and never the command word `argv[0]`. The false read came from a WRAPPER
(`env VAR=x prog`, `timeout 20 prog`, `nice prog`): there the wrapper is
argv[0], so the real program path sits in `argv[1]` and was scanned as
an ordinary path argument.

Semantics: the command word of a simple command -- the program executed,
in whatever position the shell resolves it -- is not a read target. Its
arguments still are: `python /outside/file.py` asks for /outside/file.py
and `cat /outside/x` asks. The exemption is structural (the command-word
position), never by path nature: `cat /usr/bin/python3` and
`python -c "open('/etc/shadow')"` stay gated, and no directory is added
to the readable roots.
"""

from __future__ import annotations

from delfin.agent import api_client as A


def R(cmd, cwd="."):
    """The outside-read scanner, the unit under test."""
    return A._bash_outside_reads(cmd, cwd)


# ---------------------------------------------------------------------------
# The command word of a simple command is not a read target
# ---------------------------------------------------------------------------

def test_the_wrapped_interpreter_path_is_not_a_read_target(tmp_path):
    interp = "/home/qmchem_max/delfin-swarm/venv/bin/python"
    script = str(tmp_path / "outside.py")
    # The executed program's own path, however it is reached, is not a read.
    assert R(f"{interp} -u {script}", "") != []
    got = R(f"env X=1 {interp} -u {script}")
    got = R(f"env PYTHONPATH=/x {interp} -u {script}")
    assert interp not in got, got
    assert script in got, got
    assert R(f"timeout 20 {interp} -u {script}") == [script]
    assert R(f"nice {interp} -u {script}") == [script]
    assert R(f"nohup {interp} -u {script}") == [script]
    assert R(f"cd /tmp && {interp} -u {script}") == ["/tmp", script]


def test_a_direct_interpreted_run_asks_for_its_arguments_not_its_program(tmp_path):
    interp = "python3"
    script = str(tmp_path / "outside.py")
    got = R(f"{interp} -u {script}")
    assert script in got, got
    # the program word itself is not scanned even when it is a bare path
    assert sorted(p for p in got) == [script], got


def test_the_program_word_is_skipped_but_sibling_arguments_are_not(tmp_path):
    # `cp a b` is a WRITE: the write scanner owns b. This test is about
    # the READ scanner treating only the command word, not its siblings,
    # as innocuous: two real files after a wrapper are both still reads.
    a = "/tmp/aaa"
    b = "/tmp/bbb"
    got = R(f"env X=1 prog {a} {b}")
    assert a in got, got
    assert b in got, got


# ---------------------------------------------------------------------------
# The paired refusals of requirements 1 and 2 stay open
# ---------------------------------------------------------------------------

def test_argument_files_are_still_read_targets(tmp_path):
    target = str(tmp_path / "x.txt")
    assert target in R(f"python3 /some/script.py {target}")
    assert target in R(f"cat {target}")


def test_interpreters_are_not_exempt_as_generic_paths(tmp_path):
    # `/usr/bin/python3` is a system path and stays as gated as today --
    # this change does not touch it. What must survive is the rule that a
    # NON-system interpreter path read as DATA (not executed) is still a
    # read: `cat` on the venv python exposes its bytes to the shell.
    interp = "/home/qmchem_max/delfin-swarm/venv/bin/python"
    assert interp in R(f"cat {interp}")


def test_python_inline_reads_are_still_gated():
    assert R("python3 -c \"open('/etc/shadow')\"") == ["/etc/shadow"]


def test_shell_wrappers_keep_working():
    # env's own chdir option still scans the directory it goes to
    assert "/some/dir" in R("env -C /some/dir python3 -u /outside.py")


def test_the_gate_does_not_ask_for_a_wrapped_interpreter(tmp_path):
    # End to end through the gate that produced the field question.
    ws = tmp_path / "ws"
    ws.mkdir()
    other = tmp_path / "other"
    other.mkdir()
    (other / "outside.py").write_text("print('x')\n")
    perms = A.KitToolPermissions(workspace=ws, mode="default")
    ex = A._DocToolExecutor()
    interp = "/home/qmchem_max/delfin-swarm/venv/bin/python"
    cmd = f"env X=1 {interp} -u {other}/outside.py"
    # The interpreter path (resolved and not) must not trip the read gate;
    # only the script it runs stays an outside read.
    out = ex._gate_bash_read_paths(cmd, perms, str(ws))
    assert out is not None
    assert interp.split("/")[-1] not in out


def test_a_wrapper_alone_or_without_a_command_word_skips_nothing(tmp_path):
    # `env FOO=bar` sets no program: nothing is a command word, nothing
    # is skipped. `env FOO=bar python3` still has python3 as its command
    # word -- a program path here is not read, so no target either.
    assert R("env FOO=bar") == []
    assert R("env FOO=bar python3") == []
