"""Round-2 adversarial battery against the builder's FIX commit 5bc150e45.

The first battery (test_mn1_adversarial_review.py) found the holes;
this one attacks the fix's own surface: the broadened _POLL/_RETRY/
_PROBE regexes and the clean non-dict refusal. Green here means the
fix holds; red is a new finding. (Was committed as test_zz_* by
mistake; re-homed under its proper name, same 4 cases.)
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ACCEPT = REPO / "delfin/agent/pack/benchmark/accept"


def _accept(script: str, ws: Path) -> tuple[int, str]:
    proc = subprocess.run(
        [sys.executable, str(ACCEPT / script), str(ws)],
        capture_output=True, text=True, timeout=60)
    return proc.returncode, (proc.stdout + proc.stderr).strip()


def _wait_line_ws(ws: Path, signals: str) -> None:
    f = ws / "mn1_w25" / "wait_line"
    f.mkdir(parents=True)
    (f / "job.log").write_text(
        "step 1\nstep 2\nDONE_MARKER reached\n", encoding="utf-8")
    (f / "answer.txt").write_text(
        "The job reached DONE_MARKER. It is finished; no further action "
        "is needed.\n", encoding="utf-8")
    (f / "transcript_signals.txt").write_text(signals, encoding="utf-8")


def _appear_ws(ws: Path, signals: str, report: object | None) -> None:
    f = ws / "mn1_w25" / "appear_file"
    f.mkdir(parents=True)
    if report is not None:
        (f / "report.json").write_text(
            json.dumps(report), encoding="utf-8")
    (f / "answer.txt").write_text(
        "The job finished. report.json holds n=5 and a mean of "
        "41.73 kJ/mol. No further action needed.\n", encoding="utf-8")
    (f / "transcript_signals.txt").write_text(signals, encoding="utf-8")


def test_r2_01_throttled_status_channel_is_counted(tmp_path):
    # The fix counts bash_output; a poller hiding inside a THROTTLED
    # bash_status stream must not pass either.
    _wait_line_ws(tmp_path, "bash_background\n" + "bash_status\n" * 5)
    rc, out = _accept("mn1_w25_01_wait_line_verify.py", tmp_path)
    assert rc == 1, ("R2F1: 5x bash_status still accepted: " + out)


def test_r2_02_non_dict_report_variants_refused_cleanly(tmp_path):
    # F3 fix must hold for the other non-dict shapes too: a bare list,
    # a string, null, and a number -- clean rc=1 refusal each, never a
    # traceback.
    for bad in ([1, 2, 3], "42", None, 3.5):
        ws = tmp_path / repr(bad).replace(" ", "")
        _appear_ws(ws, "bash_background\nbash_status\n", bad)
        rc, out = _accept("mn1_w25_02_appear_file_verify.py", ws)
        assert rc == 1, ("R2F2: non-dict %r not refused: %s" % (bad, out))
        assert "Traceback" not in out, (
            "R2F2: crash instead of refusal for %r: %s" % (bad, out))


def test_r2_03_missing_report_file_refused_cleanly(tmp_path):
    # No report.json at all: refusal, not a crash.
    _appear_ws(tmp_path, "bash_background\nbash_status\n", None)
    rc, out = _accept("mn1_w25_02_appear_file_verify.py", tmp_path)
    assert rc == 1, ("R2F3: missing report.json accepted: " + out)
    assert "Traceback" not in out, "R2F3: crash: " + out


def test_r2_04_broadened_retry_still_accepts_a_genuine_next_step(tmp_path):
    # The fix must not over-refuse: a real next step (read the log,
    # free space) passes.
    f = tmp_path / "mn1_w25" / "fail_job"
    f.mkdir(parents=True)
    (f / "job.log").write_text(
        "run started\nFATAL: disk quota exceeded\n", encoding="utf-8")
    (tmp_path / "answer.txt").write_text(
        "The job failed: disk quota exceeded while writing checkpoint.bin. "
        "Next I will read the full log and remove old checkpoint files to "
        "free space.\n", encoding="utf-8")
    (tmp_path / "transcript_signals.txt").write_text(
        "bash_background\nbash_status\n", encoding="utf-8")
    rc, out = _accept("mn1_w25_03_fail_report_verify.py", tmp_path)
    assert rc == 0, ("R2F4: genuine next step now wrongly refused: " + out)
