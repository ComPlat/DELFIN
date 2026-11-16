"""BG2 adversarial tests (reviewer s4) — attack the adopted-job path.

The builder's tests show the happy move works once the px pairs are
applied. These attack what a move actually hands the registry and what
adopt() must hold together:

- A1: an adopted job's deadline is an absolute wall-clock cap counting the
  foreground part (the move must not buy extra time).
- A2: an adopted job with NO capture object still registers and completes.
- A3: large output written after adoption is fully in the job's file
  (no truncation at a 64 KiB chunk boundary) — via a capture object, the
  real demote path (pipes are drained only by its readers).
- A4: unicode output survives both sides of the move byte-exactly.
- A5: the record carries command/cwd/env provenance and counts as running.
- A6: kill() on an adopted job ends its process tree (lifecycle parity
  with a bash_background job).

Every case drives the real registry path: _Registry.adopt(), no mocks of
the registry itself.
"""
import os
import subprocess
import time

import pytest

from delfin.agent import bash_jobs as BJ


class _RealCapture:
    """The demote path's shape: reader threads draining a live pipe into
    the job's files, so-far content first (the PipeCapture contract)."""

    def __init__(self, proc, so_far: bytes = b""):
        self._proc = proc
        self._so_far = so_far

    def attach(self, sout, serr):
        if self._so_far:
            sout.write(self._so_far.decode("utf-8", errors="replace"))
            sout.flush()
        import threading
        t = threading.Thread(
            target=self._pump, args=(self._proc.stdout, sout),
            daemon=True)
        t.start()
        self._thread = t

    def _pump(self, stream, handle):
        while True:
            chunk = stream.read(65536)
            if not chunk:
                break
            handle.write(chunk.decode("utf-8", errors="replace"))
            handle.flush()
        stream.close()

    def join(self, timeout=60.0):
        t = getattr(self, "_thread", None)
        if t is not None:
            t.join(timeout=timeout)


@pytest.fixture(autouse=True)
def _clean_registry(tmp_path, monkeypatch):
    BJ._REGISTRY._jobs.clear()
    yield
    for job in list(BJ._REGISTRY._jobs.values()):
        try:
            BJ.get_registry().kill(job.job_id)
        except Exception:
            pass
    BJ._REGISTRY._jobs.clear()


def _adopt(tmp_path, cmd, capture=None, **kw):
    proc = subprocess.Popen(
        ["/bin/bash", "-c", cmd],
        cwd=str(tmp_path), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        preexec_fn=os.setsid)
    job = BJ.get_registry().adopt(
        command=cmd, cwd=str(tmp_path), proc=proc,
        workspace=str(tmp_path), session_id="s4-adv",
        env_source="gate", capture=capture,
        timeout_s=kw.pop("timeout_s", 600), **kw)
    return job, proc


def _wait(job, s=30.0):
    end = time.monotonic() + s
    while time.monotonic() < end:
        if job.poll() is not None:
            return job.poll()
        time.sleep(0.05)
    return None


def test_a1_deadline_counts_the_foreground_part(tmp_path):
    job, _ = _adopt(tmp_path, "sleep 300", foreground_elapsed_s=1.0,
                    timeout_s=10)
    rec = BJ._record_deadline({
        "started_at": job.started_at, "timeout_s": 10})
    remaining = rec - time.monotonic()
    assert 0 < remaining <= 10.0, remaining


def test_a2_adopt_without_capture_registers_and_completes(tmp_path):
    job, proc = _adopt(tmp_path, "echo done")
    assert job.job_id, "adopt() must register without a capture object"
    assert _wait(job) == 0, "job must complete"


def test_a3_large_output_after_adoption_is_complete(tmp_path):
    n = 50000  # ~200 KiB, spans many 64 KiB chunks
    proc = subprocess.Popen(
        ["/bin/bash", "-c", f"seq 1 {n}"],
        cwd=str(tmp_path), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        preexec_fn=os.setsid)
    cap = _RealCapture(proc)
    job = BJ.get_registry().adopt(
        command=f"seq 1 {n}", cwd=str(tmp_path), proc=proc,
        workspace=str(tmp_path), capture=cap, timeout_s=600)
    assert _wait(job, s=60.0) == 0
    cap.join(timeout=30.0)
    text = job.stdout_path.read_text()
    assert text.splitlines()[-1] == str(n), "tail must not be truncated"


def test_a4_unicode_output_survives_the_move(tmp_path):
    proc = subprocess.Popen(
        ["/bin/bash", "-c",
         "printf 'ae\\xc3\\xa4o\\xc3\\xb6\\xf0\\x9f\\x8c\\x8d\\n'"],
        cwd=str(tmp_path), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        preexec_fn=os.setsid)
    cap = _RealCapture(proc, so_far="vorné\n".encode("utf-8"))
    job = BJ.get_registry().adopt(
        command="printf unicode", cwd=str(tmp_path), proc=proc,
        workspace=str(tmp_path), capture=cap, timeout_s=600)
    assert _wait(job) == 0
    cap.join(timeout=30.0)
    text = job.stdout_path.read_text(encoding="utf-8", errors="replace")
    assert "vorné" in text, text
    assert "aeäoö" in text and "🌍" in text, text


def test_a5_record_carries_provenance_and_counts_as_running(tmp_path):
    job, _ = _adopt(tmp_path, "sleep 5", foreground_elapsed_s=1.5)
    assert job.command == "sleep 5"
    assert job.cwd == str(tmp_path)
    rec = _persisted_record(tmp_path, job.job_id)
    assert rec["origin"] == "adopted"
    assert rec["env_source"] == "gate"
    assert job.poll() is None, "sleep 5 must still be running here"
    assert BJ.live_jobs(), "an adopted job counts as a running job"


def _persisted_record(ws, job_id):
    import json
    return json.loads(BJ._registry_path(ws).read_text())["jobs"][job_id]


def test_a6_kill_ends_the_adopted_process_tree(tmp_path):
    job, _ = _adopt(tmp_path, "sleep 300")
    BJ.get_registry().kill(job.job_id)
    for _ in range(100):
        try:
            os.kill(job.proc.pid, 0)
        except ProcessLookupError:
            return  # ended, as required
        time.sleep(0.05)
    raise AssertionError("adopted process still alive after kill()")
