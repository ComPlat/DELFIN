"""Adversarial reviewer tests for package F, part 2: the USR1 flush path.

Reviewer s4, wave 11.  Phase 3 of the task: the USR1 trap must still leave
finished systems resumable.  These tests execute the RENDERED trap line under
the sbatch's own `set -euo pipefail` and drive the runner's flush handshake
with bounded, event-ordered waits (release-file handshakes, never machine
speed: a fast and a slow machine see the same events in the same order).
"""
from __future__ import annotations

import json
import subprocess
import sys
import threading
import time
from pathlib import Path

from delfin.cluster_bench import prepare as cp
from delfin.cluster_bench import provenance as cprov
from delfin.cluster_bench import runner as crun
from delfin.cluster_bench.slurm_script import cbatch_render_sbatch

CISPLATIN = "[Cl][Pt-2]([Cl])([NH3+])[NH3+]"


def _man(tool_python=None) -> dict:
    return {"tool": "architector", "label": "adv",
            "settings": {"tool_python": tool_python, "timeout_base_s": 21600,
                         "speed_factor": "1.0", "workers": 48, "threads": 1},
            "sets": {"main": {"n_shards": 2}}}


def _trap_line(text: str) -> str:
    lines = text.splitlines()
    return next(ln for ln in lines if ln.startswith("trap "))


def _run_trap(tmp_path: Path, local_root: str) -> subprocess.CompletedProcess:
    """Run the rendered trap line under `set -euo pipefail`, fire USR1 at the script
    itself, and report what happened.  `local_root` plays DELFIN_CLUSTER_LOCAL_ROOT."""
    text = cbatch_render_sbatch(tmp_path / "run", _man(tool_python="/ws/venv/bin/python"))
    script = tmp_path / "trap.sh"
    script.write_text("#!/bin/bash\nset -euo pipefail\n"
                      f"export DELFIN_CLUSTER_LOCAL_ROOT={local_root!r}\n"
                      + _trap_line(text) + "\n"
                      "sleep 0.3 &\nPID=$!\nkill -USR1 $$\n"
                      "wait $PID || true\necho SURVIVED\n")
    return subprocess.run(["bash", str(script)], capture_output=True, text=True, timeout=60)


def test_the_usr1_trap_touches_the_flush_flag_and_survives(tmp_path):
    """With a local root the trap must create the flag file the runner polls for, and the
    script must keep running afterwards (the runner wait loop follows the trap)."""
    local = tmp_path / "local"
    local.mkdir()
    done = _run_trap(tmp_path, str(local))
    assert done.returncode == 0, done.stderr
    assert "SURVIVED" in done.stdout
    assert (local / ".usr1").exists(), "the trap did not create the flush flag"


def test_the_usr1_trap_is_harmless_in_fallback_mode(tmp_path):
    """Fallback (no local root): the trap must not kill the script and must not write a
    flag into the current directory -- with an empty root `touch /.usr1` would fail or
    scatter a file; the guard must hold under `set -euo pipefail`."""
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    done = _run_trap(tmp_path, "")
    assert done.returncode == 0, done.stderr
    assert "SURVIVED" in done.stdout
    assert not (tmp_path / ".usr1").exists() and not (run_dir / ".usr1").exists(), \
        "the trap wrote a flag although no local root was staged"


def _wait_until(what, budget_s: float) -> bool:
    """Bounded, event-ordered wait: poll a cheap condition, never a machine-speed race."""
    end = time.monotonic() + budget_s
    while time.monotonic() < end:
        if what():
            return True
        time.sleep(0.02)
    return False


def _fake_workers_dir(tmp_path: Path) -> Path:
    """A stdlib worker that writes its result, then waits on a release file the test
    controls before exiting: the runner only sees the completion AFTER the test has set
    up the situation under test (e.g. the USR1 flag).  CB_HANG_IDS chooses the waiters."""
    wd = tmp_path / "workers"
    wd.mkdir(exist_ok=True)
    (wd / "env_probe.py").write_text("import json\nprint(json.dumps({'packages': {}}))\n")
    (wd / "bench_worker.py").write_text(
        "import json, os, sys, time\n"
        "tool, ref, specs, archive, work = sys.argv[1:6]\n"
        "os.makedirs(os.path.join(archive, '_meta'), exist_ok=True)\n"
        "open(os.path.join(archive, ref + '.xyz'), 'w').write('1\\n%s frame0\\nFe 0 0 0\\n' % ref)\n"
        "json.dump({'refcode': ref, 'tool': tool, 'status': 'ok', 'n_frames': 1},\n"
        "          open(os.path.join(archive, '_meta', ref + '.json'), 'w'))\n"
        "ran = os.environ.get('CB_RAN')\n"
        "if ran:\n"
        "    open(ran, 'a').write(ref + '\\n')\n"
        "if ref in os.environ.get('CB_HANG_IDS', '').split(','):\n"
        "    while not os.path.exists(os.environ['CB_RELEASE']):\n"
        "        time.sleep(0.02)\n")
    return wd


def _prepared_run(tmp_path, monkeypatch, ids):
    wd = _fake_workers_dir(tmp_path)
    monkeypatch.setattr(cprov, "WORKERS_DIR", wd)
    monkeypatch.setattr(crun, "WORKERS_DIR", wd)
    tool_python = tmp_path / "manifest_venv" / "bin" / "python"
    tool_python.parent.mkdir(parents=True, exist_ok=True)
    tool_python.write_text("#!/bin/sh\nexec \"$CB_REALPY\" \"$@\"\n")
    tool_python.chmod(0o755)
    inp = tmp_path / "in.txt"
    inp.write_text("".join(f"{rid};{CISPLATIN}\n" for rid in ids))
    specs = tmp_path / "specs.jsonl"
    specs.write_text("".join(json.dumps({"refcode": rid, "status": "ok", "cn": 4}) + "\n"
                             for rid in ids))
    cp.cbatch_prepare(tool="molsimplify", input_list=inp, run_dir=tmp_path / "run",
                      specs_file=specs, tool_python=str(tool_python), timeout_base=60,
                      size=50, check_tool=False)
    return tmp_path / "run"


def test_the_usr1_flag_flushes_finished_systems_without_a_further_completion(
        tmp_path, monkeypatch):
    """The trap's promise is 'flushing finished systems'.  The runner consults the flag
    only when another system COMPLETES (runner.py:338, single call site) -- on a real
    node builds take hours, so systems finished BEFORE USR1 would never reach the
    workspace before the node is reaped at wall time.  The flag must be honoured by
    itself: finished-but-pending systems must land on the workspace archive even when
    no further completion follows."""
    monkeypatch.setenv("CB_INTERP_MARKER", str(tmp_path / "interp"))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_PATHS", str(tmp_path / "paths"))
    monkeypatch.setenv("CB_RELEASE", str(tmp_path / "release"))
    run_dir = _prepared_run(tmp_path, monkeypatch, ["U1", "U2"])
    monkeypatch.setenv("CB_HANG_IDS", "U1,U2")   # both workers hang after writing results
    local = tmp_path / "local"
    monkeypatch.setenv("DELFIN_CLUSTER_LOCAL_ROOT", str(local))
    monkeypatch.setenv("DELFIN_CLUSTER_SYNC_SIZE", "32")

    # run the shard in a thread; U1+U2 write results and then hang on the release file
    result = {}

    def _run():
        result["rc"] = crun.cbatch_run_shard(run_dir, 0, log=lambda m: None)

    th = threading.Thread(target=_run)
    th.start()
    try:
        assert _wait_until(lambda: (Path(run_dir) / "out" / "main_main" / "chunk_0000"
                                    / "run_summary.jsonl").exists(), 10.0), \
            "the shard never started"
        local.mkdir(parents=True, exist_ok=True)
        (local / ".usr1").write_text("")           # the trap fired BEFORE any exit
        # no child completes after the flag: still, the finished systems must be flushed
        assert _wait_until(lambda: ((Path(run_dir) / "out" / "main_main" / "chunk_0000"
                                     / "archive" / "_meta" / "U1.json").exists()
                                    and (Path(run_dir) / "out" / "main_main" / "chunk_0000"
                                         / "archive" / "_meta" / "U2.json").exists()), 10.0), \
            "USR1 flag set, no further completion, yet finished systems were not flushed"
    finally:
        (tmp_path / "release").write_text("")      # let the workers exit
        th.join(timeout=60)
    assert result.get("rc") == 0


def test_a_system_with_a_final_record_is_skipped_so_resubmitting_the_shard_resumes(
        tmp_path, monkeypatch):
    """Phase-3 contract: the USR1 kill must never cost finished work.  A system whose
    workspace _meta/<rid>.json carries final:true is skipped on a re-run -- after the
    trap's flush and re-reap, resubmitting the array index rebuilds nothing finished."""
    monkeypatch.setenv("CB_INTERP_MARKER", str(tmp_path / "interp"))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_RAN", str(tmp_path / "ran.txt"))
    run_dir = _prepared_run(tmp_path, monkeypatch, ["R1", "R2"])
    archive_ws = (Path(run_dir) / "out" / "main_main" / "chunk_0000" / "archive")
    (archive_ws / "_meta").mkdir(parents=True, exist_ok=True)
    (archive_ws / "_meta" / "R1.json").write_text(
        json.dumps({"refcode": "R1", "tool": "molsimplify", "status": "ok",
                    "n_frames": 1, "final": True}))
    rc = crun.cbatch_run_shard(run_dir, 0, log=lambda m: None)
    assert rc == 0, "pre-finished R1 must complete the shard"
    ran = (tmp_path / "ran.txt").read_text().split()
    assert ran == ["R2"], f"worker ran {ran}, final-recorded R1 was rebuilt"
    assert (archive_ws / "_meta" / "R2.json").exists()


def test_a_crash_residue_meta_without_final_does_NOT_skip_the_system(tmp_path, monkeypatch):
    """The negative case: a meta lacking final:true (a partial write, a crash between
    xyz and meta, or last run's non-final residue) must NOT mark the system finished --
    otherwise a half-baked result survives a resubmit and masquerades as a success."""
    monkeypatch.setenv("CB_INTERP_MARKER", str(tmp_path / "interp"))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_RAN", str(tmp_path / "ran.txt"))
    run_dir = _prepared_run(tmp_path, monkeypatch, ["N1", "N2"])
    archive_ws = (Path(run_dir) / "out" / "main_main" / "chunk_0000" / "archive")
    (archive_ws / "_meta").mkdir(parents=True, exist_ok=True)
    (archive_ws / "_meta" / "N1.json").write_text(
        json.dumps({"refcode": "N1", "status": "ok", "n_frames": 1}))   # no final key
    rc = crun.cbatch_run_shard(run_dir, 0, log=lambda m: None)
    assert rc == 0
    ran = (tmp_path / "ran.txt").read_text().split()
    assert "N1" in ran, "system without a final record was skipped: stale residue read as done"
    assert (archive_ws / "_meta" / "N1.json").exists()
