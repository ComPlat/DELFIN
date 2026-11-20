"""CM1: the documented API route is reachable from the REAL entry point.

The Code Mode skill (delfin/agent/pack/skills/code-mode.md) tells the model
two commands:

    python -m delfin.agent.code_mode stub > tools_api.py
    python -m delfin.agent.code_mode run my_program.py

The first was missing from ``main()``'s dispatch: ``main(["stub"])`` printed
the ``run`` usage and exited 2, so ``generate_stub()`` — the typed API the
skill documents — could not be obtained at all, and the generated stub's
``make_client`` import never happened. These tests drive ``main()`` (the
process the bash tool actually starts) and a real ``RpcServer``, never
``generate_stub`` directly, so they measure the route and not the function.

Socket paths use a short mkdtemp base: pytest's tmp_path can exceed the
108-byte AF_UNIX sockaddr limit and then bind() fails for a reason that has
nothing to do with the route.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading

from delfin.agent import code_mode


class _RecordingExecutor:
    """Stands in for KitToolExecutor at the public entry the RPC calls."""

    def __init__(self):
        self.calls = []

    def execute(self, name, arguments, permissions=None):
        self.calls.append((name, arguments))
        return json.dumps({"tool": name, "arguments": arguments})


def _short_base():
    """A directory short enough for an AF_UNIX socket path.

    The RPC socket is ``<base>/.code_mode_rpc/rpc-<32-hex>.sock``, and an
    AF_UNIX sockaddr tops out at 108 bytes. ``tempfile.gettempdir()`` honours
    TMPDIR, which a sandboxed or CI host can make long enough to overflow
    that on its own; bind() would then fail for a reason unrelated to the
    route under test. Take the shortest candidate that is a writable
    directory, so this holds wherever the suite runs instead of assuming
    /tmp exists.
    """
    budget = 108 - len(os.path.join("", ".code_mode_rpc", "rpc-" + "0" * 32 + ".sock"))
    candidates = [tempfile.gettempdir(), "/tmp"]
    for base in sorted({c for c in candidates if c and os.path.isdir(c)},
                       key=len):
        if len(base) + 1 + len("cm1x_") + 8 <= budget:
            return tempfile.mkdtemp(prefix="cm1x_", dir=base)
    raise RuntimeError("no temp directory short enough for an AF_UNIX socket")


def test_the_documented_stub_command_prints_the_api(capsys):
    """``main(["stub"])`` writes the typed API to stdout and succeeds.

    Red before this route existed: the same call printed the ``run`` usage
    and returned 2, so the documented command produced no API module.
    """
    rc = code_mode.main(["stub"])
    out = capsys.readouterr().out
    assert rc == 0
    for name in code_mode.READ_ONLY_TOOLS:
        assert "def %s(" % name in out
    assert "write_file" not in out


def test_unknown_and_missing_subcommands_still_refuse(capsys):
    """Adding the route must not turn a bad invocation into a success."""
    assert code_mode.main([]) == 2
    assert code_mode.main(["nonsense"]) == 2
    assert code_mode.main(["run"]) == 2
    err = capsys.readouterr().err
    assert "stub" in err, "the usage line must name the route it now has"


def test_a_generated_stub_reaches_a_tool_through_the_real_route():
    """End to end: the documented stub, imported by a program the real
    runner starts, reaches an executor call over a real RpcServer.

    Neither ``generate_stub`` nor ``make_client`` is called by this test:
    the stub it exercises came from ``main(["stub"])`` and the client is the
    one the stub's own text imports.
    """
    base = _short_base()
    try:
        sock = code_mode.rpc_socket_path(base, "e2e")
        executor = _RecordingExecutor()
        server = code_mode.RpcServer(executor, permissions=None)
        server.bind(sock)
        threading.Thread(
            target=lambda: server.serve_one(idle_timeout_s=30),
            daemon=True).start()

        stub_path = os.path.join(base, "tools_api.py")
        written = subprocess.run(
            [sys.executable, "-m", "delfin.agent.code_mode", "stub"],
            capture_output=True, text=True, timeout=120,
            cwd=base,
            env=dict(os.environ,
                     PYTHONPATH=os.path.dirname(
                         os.path.dirname(os.path.dirname(
                             os.path.abspath(code_mode.__file__))), )))
        assert written.returncode == 0, written.stderr
        with open(stub_path, "w") as fh:
            fh.write(written.stdout)

        program = os.path.join(base, "census.py")
        with open(program, "w") as fh:
            fh.write(
                "import json, tools_api as api\n"
                "print(json.dumps(api.get_calc_info(calc_id='c1')))\n")

        os.environ["DELFIN_CODE_MODE_RPC"] = sock
        try:
            rc = code_mode.main(["run", program])
        finally:
            os.environ.pop("DELFIN_CODE_MODE_RPC", None)
        assert rc == 0
        assert executor.calls == [("get_calc_info", {"calc_id": "c1"})]
    finally:
        try:
            server.close()
        except Exception:
            pass
        shutil.rmtree(base, ignore_errors=True)


def test_a_stub_without_the_env_var_refuses_instead_of_pretending():
    """Without the socket path the stub raises and names the runner.

    This is the property ``make_client`` exists for: no connection means no
    call, never a silent success.
    """
    base = _short_base()
    try:
        stub_path = os.path.join(base, "tools_api.py")
        written = subprocess.run(
            [sys.executable, "-m", "delfin.agent.code_mode", "stub"],
            capture_output=True, text=True, timeout=120, cwd=base,
            env=dict(os.environ,
                     PYTHONPATH=os.path.dirname(
                         os.path.dirname(os.path.dirname(
                             os.path.abspath(code_mode.__file__))), )))
        assert written.returncode == 0, written.stderr
        with open(stub_path, "w") as fh:
            fh.write(written.stdout)

        probe = os.path.join(base, "probe.py")
        with open(probe, "w") as fh:
            fh.write(
                "import tools_api as api\n"
                "try:\n"
                "    api.get_calc_info(calc_id='c1')\n"
                "except RuntimeError as exc:\n"
                "    print('REFUSED:', exc)\n")
        env = dict(os.environ)
        env.pop("DELFIN_CODE_MODE_RPC", None)
        env["PYTHONPATH"] = os.path.dirname(
            os.path.dirname(os.path.dirname(os.path.abspath(
                code_mode.__file__))))
        run = subprocess.run([sys.executable, probe], capture_output=True,
                             text=True, timeout=120, cwd=base, env=env)
        assert run.returncode == 0, run.stderr
        assert "REFUSED:" in run.stdout
        assert "DELFIN_CODE_MODE_RPC" in run.stdout
    finally:
        shutil.rmtree(base, ignore_errors=True)
