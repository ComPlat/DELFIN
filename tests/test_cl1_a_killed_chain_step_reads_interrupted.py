"""A killed chain step reads "interrupted", not "running".

Input: a calculation folder that is a DELFIN chain step, with the files a
killed ORCA optimisation leaves (a fresh output without a termination
line, temporary files), and a chain run ended by SIGTERM.

Before: calc_status judged by the files and the clock alone. In JC3's live
run 6 (2026-11-15) it called the killed step "running" -- the files had
been written minutes before -- and a run ended by ``scancel`` read "still
running" in the chain until its heartbeat was two minutes old, so the
agent's resume in that window was refused.
Now calc_status and extract_energy_table ask the chain's runner record
when the folder is a step's (or an earlier attempt's): "interrupted" when
the run verifiably ended, "running" when it verifiably runs; an exit code
or a termination line still has the last word. A chain run ended by
SIGTERM/SIGHUP records that at once and then ends exactly as before; a
signal some handler catches and survives records nothing.
"""
from __future__ import annotations

import json
import os
import signal
import socket
import subprocess
import sys
import textwrap
import threading
import time
from pathlib import Path

import pytest

from delfin import api
from delfin.agent import proc_identity
from delfin.tools import chain_liveness as cl
from delfin.tools import chain_state as cs

REPO = Path(__file__).resolve().parents[1]


def _ended_runner() -> dict:
    p = subprocess.Popen([sys.executable, "-c", "pass"])
    start = proc_identity.process_start(p.pid)
    p.wait()
    return {"pid": p.pid, "start": start, "host": socket.gethostname(),
            "pid_ns": cl._pid_namespace()}


def _killed_orca_dir(d: Path) -> Path:
    """What a killed optimisation leaves: no exit code, no termination."""
    d.mkdir(parents=True, exist_ok=True)
    (d / "calc.inp").write_text("! PBE0 def2-SVP Opt\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n")
    (d / "calc.out").write_text("synthetic stand-in output\nGEOMETRY OPTIMIZATION CYCLE 2\n")
    (d / "calc.densities.tmp").write_text("x")
    return d


def _chain(root: Path, runner) -> tuple[cs.ChainState, Path]:
    chain = cs.ChainState.create(name="demo", work_dir=root)
    chain.add_step("xtb_opt", "xtb_opt", params={"charge": 0})
    chain.add_step("orca_opt", "orca_opt", params={"charge": 0})
    step_dir = _killed_orca_dir(root / "01_orca_opt")
    chain.steps[1].mark_submitted("local", runner=runner)
    chain.steps[1].work_dir = str(step_dir)
    chain.save()
    return chain, step_dir


def test_a_killed_step_is_interrupted_with_the_chain_as_evidence(tmp_path):
    _, step_dir = _chain(tmp_path / "c", _ended_runner())
    st = api.calculation_status(str(step_dir))
    assert st.state == "interrupted", st
    says = [e["says"] for e in st.evidence if e["source"] == "chain_state.json"]
    assert says and "orca_opt" in says[0] and "not running" in says[0]


def test_a_step_whose_runner_lives_is_running(tmp_path):
    _, step_dir = _chain(tmp_path / "c", cl.this_runner())
    assert api.calculation_status(str(step_dir)).state == "running"


def test_an_earlier_attempt_is_interrupted(tmp_path):
    chain, step_dir = _chain(tmp_path / "c", _ended_runner())
    step = chain.steps[1]
    step.mark_failed("killed")
    step.reset_for_retry(attempt_dir=str(tmp_path / "c" / "orca_opt" / "attempt_2"),
                         end="killed/job gone")
    chain.save()
    st = api.calculation_status(str(step_dir))
    assert st.state == "interrupted"
    assert any("attempt 1 ended" in e["says"] for e in st.evidence)


def test_the_files_keep_the_last_word(tmp_path):
    _, step_dir = _chain(tmp_path / "c", _ended_runner())
    with (step_dir / "calc.out").open("a") as fh:
        fh.write("****ORCA TERMINATED NORMALLY****\n")
    assert api.calculation_status(str(step_dir)).state not in ("interrupted", "running")


def test_a_folder_outside_any_chain_is_unchanged(tmp_path):
    d = _killed_orca_dir(tmp_path / "plain")
    st = api.calculation_status(str(d))
    assert st.state != "interrupted"
    assert not [e for e in st.evidence if e["source"] == "chain_state.json"]


def test_the_energy_table_says_interrupted_too(tmp_path):
    _, step_dir = _chain(tmp_path / "c", _ended_runner())
    rows = api.extract_energy_table([str(step_dir)], ["single_point"])
    assert rows[0]["state"] == "interrupted"


_CHILD = textwrap.dedent('''
    import sys, time
    from pathlib import Path
    import delfin.tools._runner as runner
    wd, flag = Path(sys.argv[1]), Path(sys.argv[2])
    def slow(step_name, **kw):
        flag.write_text("started")
        time.sleep(120)
    runner.run_step = slow
    {prelude}
    from delfin.tools.pipeline import Pipeline
    geom = wd.parent / "start.xyz"
    geom.parent.mkdir(parents=True, exist_ok=True)
    geom.write_text("2\\n\\nH 0 0 0\\nH 0 0 0.74\\n")
    p = Pipeline(name="demo")
    p.add("xtb_opt", charge=0)
    p.run(work_dir=wd, geometry=geom, cores=1)
''')


def _run_child(tmp_path, sig, prelude=""):
    wd, flag = tmp_path / "runs" / "a", tmp_path / "started"
    env = dict(os.environ, HOME=str(tmp_path / "home"), PYTHONPATH=str(REPO))
    env.pop("DELFIN_STATE", None)
    child = subprocess.Popen(
        [sys.executable, "-c", _CHILD.replace("{prelude}", prelude), str(wd), str(flag)],
        env=env)
    deadline = time.time() + 60
    while not flag.exists() and time.time() < deadline and child.poll() is None:
        time.sleep(0.1)
    assert flag.exists(), "the child never reached the step"
    child.send_signal(sig)
    return child, wd


@pytest.mark.parametrize("sig", [signal.SIGTERM, signal.SIGHUP])
def test_a_signal_records_the_end_and_still_ends_the_run(tmp_path, sig):
    child, wd = _run_child(tmp_path, sig)
    assert child.wait(timeout=30) == -sig          # ended exactly as before
    step = cs.ChainState.load(wd / "chain_state.json").steps[0]
    assert step.state == "submitted"
    assert step.runner["ended"] == sig.name
    alive, why = cl.runner_liveness(step)
    assert alive == cl.GONE and sig.name in why


def test_a_signal_a_handler_survives_records_nothing(tmp_path):
    prelude = "import signal; signal.signal(signal.SIGHUP, lambda *a: None)"
    child, wd = _run_child(tmp_path, signal.SIGHUP, prelude)
    time.sleep(1.0)
    try:
        assert child.poll() is None                # still running
        data = json.loads((wd / "chain_state.json").read_text())
        assert "ended" not in (data["steps"][0].get("runner") or {})
    finally:
        child.kill()
        child.wait(timeout=30)


def test_outside_the_main_thread_nothing_is_installed():
    seen = []

    def body():
        before = signal.getsignal(signal.SIGTERM)
        with cl.ended_on_signal(seen.append):
            seen.append(signal.getsignal(signal.SIGTERM) is before)

    t = threading.Thread(target=body)
    t.start()
    t.join()
    assert seen == [True]


def test_the_handlers_are_restored_after_the_step():
    before = {s: signal.getsignal(s) for s in (signal.SIGTERM, signal.SIGHUP)}
    with cl.ended_on_signal(lambda name: None):
        assert signal.getsignal(signal.SIGTERM) is not before[signal.SIGTERM]
    assert {s: signal.getsignal(s) for s in before} == before
