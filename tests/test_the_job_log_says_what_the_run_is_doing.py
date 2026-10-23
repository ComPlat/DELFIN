"""The job log has to survive a multi-day run and still be readable.

A 60 h job syncs every 15 min. Before these guards the entire log of such a
run was ~200 repetitions of "Copied 2 result files (periodic, profile=minimal)"
— no timestamps, no progress, and a benign rsync exit 24 dressed up as a
warning. Nothing in it could distinguish an optimisation that was converging
from one that had been stuck on the same gradient for a hundred cycles.

These tests exercise the real function bodies out of submit_delfin.sh rather
than re-implementing them, so they fail when the template drifts.
"""

import re
import subprocess
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "submit_templates"
    / "submit_delfin.sh"
)


def _extract_functions(*names: str) -> str:
    """Return the named shell functions verbatim from the wrapper template."""
    text = _SCRIPT.read_text(encoding="utf-8")
    out = []
    for name in names:
        match = re.search(
            rf"^{re.escape(name)}\(\) \{{\n(.*?)^\}}$",
            text,
            re.MULTILINE | re.DOTALL,
        )
        assert match, f"{name}() not found in {_SCRIPT.name}"
        out.append(match.group(0))
    return "\n\n".join(out)


def _run_bash(script: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", "-c", script],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=str(cwd) if cwd else None,
    )


# --------------------------------------------------------------------------
# say_once / flush_repeats
# --------------------------------------------------------------------------

_SAY = _extract_functions("flush_repeats", "say_once")
_SAY_PREAMBLE = '_SAY_LAST=""\n_SAY_REPEATS=0\n' + _SAY


def test_identical_sync_lines_collapse_into_a_count():
    result = _run_bash(
        _SAY_PREAMBLE
        + """
        for i in 1 2 3 4 5; do say_once "Copied 2 result files (periodic)."; done
        say_once "Copied 7 result files (periodic)."
        """
    )
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines == [
        "Copied 2 result files (periodic).",
        "  (previous line repeated 4x)",
        "Copied 7 result files (periodic).",
    ], result.stdout


def test_a_warning_is_never_swallowed_by_the_collapsing():
    """flush_repeats runs before anything important, so the repeat count can
    never separate a warning from the line it belongs to."""
    result = _run_bash(
        _SAY_PREAMBLE
        + """
        say_once "Copied 2 result files (periodic)."
        say_once "Copied 2 result files (periodic)."
        flush_repeats
        echo "WARNING: rsync exited with code 23"
        """
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [
        "Copied 2 result files (periodic).",
        "  (previous line repeated 1x)",
        "WARNING: rsync exited with code 23",
    ]


def test_say_once_does_not_trip_errexit_when_the_line_repeats():
    """The periodic loop runs under `set -e`; a repeated line must not look
    like a failing command."""
    result = _run_bash(
        "set -euo pipefail\n"
        + _SAY_PREAMBLE
        + """
        say_once "same"
        say_once "same"
        echo REACHED
        """
    )
    assert result.returncode == 0, result.stderr
    assert "REACHED" in result.stdout


# --------------------------------------------------------------------------
# progress_heartbeat
# --------------------------------------------------------------------------

_HEARTBEAT = _extract_functions(
    "_slurm_field",
    "_hms_to_seconds",
    "get_remaining_seconds",
    "flush_repeats",
    "progress_heartbeat",
)

_HEARTBEAT_PREAMBLE = (
    '_SAY_LAST=""\n'
    "_SAY_REPEATS=0\n"
    "SYNC_INTERVAL=900\n"
    # No SLURM here: get_remaining_seconds must degrade to "unknown", not fail.
    "scontrol() { return 1; }\n" + _HEARTBEAT
)


def _heartbeat_over(tmp_path: Path, filename: str, body: str) -> str:
    run_dir = tmp_path / "run"
    (run_dir / "ESD").mkdir(parents=True, exist_ok=True)
    (run_dir / "ESD" / filename).write_text(body, encoding="utf-8")
    result = _run_bash(
        f'RUN_DIR="{run_dir}"\n' + _HEARTBEAT_PREAMBLE + "\nprogress_heartbeat\n"
    )
    assert result.returncode == 0, result.stderr
    return result.stdout


def test_heartbeat_reports_the_optimisation_cycle_and_gradient(tmp_path):
    """The T1 job that ran 139 cycles on an unchanging gradient must be
    visible as exactly that, from the log alone."""
    body = """
         *                GEOMETRY OPTIMIZATION CYCLE 139            *
          ----------------------|Geometry convergence|-------------------------
          Item                value                   Tolerance       Converged
          Energy change       0.0000002120            0.0000050000      YES
          RMS gradient        0.0026502881            0.0001000000      NO
          MAX gradient        0.0250745046            0.0003000000      NO
    """
    out = _heartbeat_over(tmp_path, "T1.out", body)
    assert "ESD/T1.out" in out
    assert "opt cycle 139" in out
    assert "0.0250745046" in out
    assert "walltime left: unknown" in out


def test_heartbeat_reports_numerical_hessian_progress(tmp_path):
    """120 of 438 displacements after two days is the number that should have
    told someone the Hessian was never going to finish."""
    body = "\n".join(
        f"\t<< Calculating gradient on displaced geometry {i:3d} (of 438) >>"
        for i in range(1, 121)
    )
    out = _heartbeat_over(tmp_path, "S1.out", body)
    assert "Hessian 120/438" in out


def test_heartbeat_stays_silent_when_nothing_moved(tmp_path):
    """An output with no recognisable progress marker must not produce a
    header line — the heartbeat is signal, not decoration."""
    out = _heartbeat_over(tmp_path, "S0.out", "ORCA header, nothing else\n")
    assert out.strip() == ""


def test_heartbeat_ignores_the_wrapper_and_goat_worker_logs(tmp_path):
    """delfin_<jobid>.out is the log we are writing into, and GOAT spawns
    hundreds of short-lived worker outputs; neither belongs in the heartbeat."""
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True)
    (run_dir / "delfin_6235903.out").write_text(
        "         *                GEOMETRY OPTIMIZATION CYCLE   3            *\n",
        encoding="utf-8",
    )
    (run_dir / "XTB_GOAT.goat.0.12.out").write_text(
        "         *                GEOMETRY OPTIMIZATION CYCLE   7            *\n",
        encoding="utf-8",
    )
    result = _run_bash(
        f'RUN_DIR="{run_dir}"\n' + _HEARTBEAT_PREAMBLE + "\nprogress_heartbeat\n"
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == ""


def test_heartbeat_survives_a_missing_run_dir(tmp_path):
    """It runs from the periodic loop under `set -e`, after the scratch may
    already have been torn down."""
    result = _run_bash(
        "set -euo pipefail\n"
        + f'RUN_DIR="{tmp_path / "gone"}"\n'
        + _HEARTBEAT_PREAMBLE
        + "\nprogress_heartbeat\necho REACHED\n"
    )
    assert result.returncode == 0, result.stderr
    assert "REACHED" in result.stdout


# --------------------------------------------------------------------------
# stall detection and walltime projection across heartbeats
# --------------------------------------------------------------------------

def _opt_output(cycle: int, grad: str) -> str:
    return (
        f"         *                GEOMETRY OPTIMIZATION CYCLE {cycle}            *\n"
        f"          MAX gradient        {grad}            0.0003000000      NO\n"
    )


def _hessian_output(done: int, total: int) -> str:
    return "\n".join(
        f"\t<< Calculating gradient on displaced geometry {i:3d} (of {total}) >>"
        for i in range(1, done + 1)
    )


def _beat(tmp_path: Path, out_name: str, body: str, *, remaining: str) -> str:
    """Run one heartbeat, keeping state between calls like the real loop does."""
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / out_name).write_text(body, encoding="utf-8")
    script = (
        f'RUN_DIR="{run_dir}"\n'
        f'DELFIN_SCRATCH="{tmp_path}"\n'
        '_SAY_LAST=""\n_SAY_REPEATS=0\nSYNC_INTERVAL=900\n'
        f'scontrol() {{ printf "%s\\n" "JobId=1 TimeLimit=2-12:00:00 RunTime={remaining}"; }}\n'
        + _HEARTBEAT
        + "\nprogress_heartbeat\n"
    )
    result = _run_bash(script)
    assert result.returncode == 0, result.stderr
    return result.stdout


def test_an_unchanging_gradient_is_called_out_as_stalled(tmp_path):
    """The T1 job sat on MAX gradient 0.0250 for ~100 cycles. Two heartbeats
    of that must say so in words, not leave it to be noticed by eye."""
    first = _beat(tmp_path, "T1.out", _opt_output(100, "0.0250745046"), remaining="0-01:00:00")
    assert "STALLED" not in first, "one observation cannot establish a stall"

    second = _beat(tmp_path, "T1.out", _opt_output(120, "0.0250745046"), remaining="0-01:00:00")
    assert "STALLED" not in second, "still only the first repeat"

    third = _beat(tmp_path, "T1.out", _opt_output(139, "0.0250745046"), remaining="0-01:00:00")
    assert "STALLED" in third, third
    assert "opt cycle 139" in third


def test_a_moving_gradient_is_never_called_stalled(tmp_path):
    _beat(tmp_path, "S1.out", _opt_output(10, "0.0090000000"), remaining="0-01:00:00")
    _beat(tmp_path, "S1.out", _opt_output(20, "0.0040000000"), remaining="0-01:00:00")
    out = _beat(tmp_path, "S1.out", _opt_output(30, "0.0010000000"), remaining="0-01:00:00")
    assert "STALLED" not in out, out


def test_a_stall_counter_resets_when_the_gradient_moves_again(tmp_path):
    """An optimisation that pauses and then finds its way out must not stay
    flagged — the warning has to keep meaning something."""
    _beat(tmp_path, "S1.out", _opt_output(10, "0.0050000000"), remaining="0-01:00:00")
    _beat(tmp_path, "S1.out", _opt_output(20, "0.0050000000"), remaining="0-01:00:00")
    _beat(tmp_path, "S1.out", _opt_output(30, "0.0050000000"), remaining="0-01:00:00")
    moved = _beat(tmp_path, "S1.out", _opt_output(40, "0.0010000000"), remaining="0-01:00:00")
    assert "STALLED" not in moved
    again = _beat(tmp_path, "S1.out", _opt_output(50, "0.0010000000"), remaining="0-01:00:00")
    assert "STALLED" not in again, "the counter must have restarted from zero"


def _age_state(tmp_path: Path, seconds: int) -> None:
    """Back-date the recorded heartbeat so a rate becomes measurable.

    Two heartbeats in the same second give elapsed=0, which the projection
    deliberately refuses to divide by; real ones are SYNC_INTERVAL apart.
    """
    import time as _time

    state = tmp_path / ".delfin_progress_state"
    aged = []
    for row in state.read_text(encoding="utf-8").splitlines():
        fields = row.split("\t")
        fields[-1] = str(int(_time.time()) - seconds)
        aged.append("\t".join(fields))
    state.write_text("\n".join(aged) + "\n", encoding="utf-8")


def test_a_hessian_that_cannot_finish_says_so(tmp_path):
    """120 of 438 displacements at the observed rate needs ~4 h. With one hour
    of walltime left that is the warning nobody got."""
    _beat(tmp_path, "S1.out", _hessian_output(100, 438), remaining="2-11:00:00")
    _age_state(tmp_path, 900)
    out = _beat(tmp_path, "S1.out", _hessian_output(120, 438), remaining="2-11:00:00")

    assert "Hessian 120/438" in out
    assert "h left at current rate" in out, out
    assert "WILL NOT FINISH" in out, out


def test_a_hessian_that_fits_is_only_reported(tmp_path):
    """Same machinery, comfortable margin: report the rate, raise nothing."""
    _beat(tmp_path, "S1.out", _hessian_output(400, 438), remaining="0-01:00:00")
    _age_state(tmp_path, 900)
    out = _beat(tmp_path, "S1.out", _hessian_output(430, 438), remaining="0-01:00:00")

    assert "Hessian 430/438" in out
    assert "WILL NOT FINISH" not in out, out


def test_the_projection_is_skipped_when_the_walltime_is_unknown(tmp_path):
    """No limit must never be read as no time left."""
    run_dir = tmp_path / "run"
    run_dir.mkdir(parents=True, exist_ok=True)
    for done in (100, 120):
        (run_dir / "S1.out").write_text(_hessian_output(done, 438), encoding="utf-8")
        result = _run_bash(
            f'RUN_DIR="{run_dir}"\nDELFIN_SCRATCH="{tmp_path}"\n'
            '_SAY_LAST=""\n_SAY_REPEATS=0\nSYNC_INTERVAL=900\n'
            'scontrol() { return 1; }\n' + _HEARTBEAT + "\nprogress_heartbeat\n"
        )
        assert result.returncode == 0, result.stderr
    assert "WILL NOT FINISH" not in result.stdout
    assert "walltime left: unknown" in result.stdout


# --------------------------------------------------------------------------
# remaining walltime
# --------------------------------------------------------------------------

_REMAINING = (
    "scontrol() { printf '%s\\n' \"JobId=1 TimeLimit=2-12:00:00 RunTime=1-06:30:00\"; }\n"
    + _extract_functions("_slurm_field", "_hms_to_seconds", "get_remaining_seconds")
)


def test_remaining_walltime_is_the_limit_minus_the_runtime():
    result = _run_bash(_REMAINING + "\nget_remaining_seconds\n")
    assert result.returncode == 0, result.stderr
    # 2-12:00:00 = 216000 s, 1-06:30:00 = 109800 s
    assert result.stdout.strip() == "106200"


@pytest.mark.parametrize(
    "fields",
    [
        "JobId=1 TimeLimit=UNLIMITED RunTime=01:00:00",
        "JobId=1 RunTime=01:00:00",
        "JobId=1 TimeLimit=2-12:00:00",
    ],
)
def test_remaining_walltime_is_empty_when_slurm_cannot_say(fields):
    """Empty means "unknown", which every caller must treat as no limit — a
    wrong number here would make the scheduler refuse jobs that do fit."""
    script = (
        f"scontrol() {{ printf '%s\\n' \"{fields}\"; }}\n"
        + _extract_functions(
            "_slurm_field", "_hms_to_seconds", "get_remaining_seconds"
        )
        + "\nget_remaining_seconds\n"
    )
    result = _run_bash(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == ""
