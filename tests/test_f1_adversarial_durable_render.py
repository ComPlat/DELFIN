"""Adversarial battery against delfin.agent.service.durable_host.render_durable_job.

Package F1 (Targets 2-4). The renderer builds a SLURM script that starts the
agent supervisor on a compute node. A defect here is a job that silently fails
to checkpoint before the walltime end, or a job script whose header a hostile
input can extend (injection into the SLURM options).

Rule under test (Task target 2): "a walltime with a --signal before the end so
the supervisor checkpoints, and a requeue/resubmit of itself". SLURM accepts
``--signal=<sig>[@<time>]`` where ``time`` (per the SLURM ``sbatch`` manual)
**must have an integer value** — i.e. seconds. Rendered 5:00 (MM:SS) is neither
an integer number of seconds nor a well-formed HH:MM:SS timestamp, so SLURM
cannot checkpoint the supervisor "x minutes before the end".

Structure of the battery:
* test_signal_* — the checkpoint-window contract; expected RED until the
  builder renders an integer-seconds --signal time.
* test_*_refused — hostile inputs (newline, unicode look-alike, leading dash)
  that must never extend a rendered SBATCH option line; the renderer refuses
  them with ValueError. Expected GREEN (the allow-lists already hold).
* test_export_all_exact / test_supervisor_command_is_immutable — Target 4
  invariants (same launcher, same settings, no wider permission).
"""

from __future__ import annotations

import pytest

from delfin.agent.service.durable_host import render_durable_job

_SIGNAL_LINE = "#SBATCH --signal="
_MMSS_SIGNAL = f"{_SIGNAL_LINE}INT@5:00"
_SECONDS_SIGNAL = f"{_SIGNAL_LINE}INT@300"


def _has(line: str, needle: str) -> bool:
    return line == needle or needle in line


def test_signal_time_is_integer_seconds_default() -> None:
    """SLURM sig_time must be integer seconds; 5 min => INT@300, not INT@5:00.

    Expected RED: the builder renders ``#SBATCH --signal=INT@5:00``.
    """
    script = render_durable_job()
    assert _has(script, _SECONDS_SIGNAL), (
        "--signal time must be integer seconds (SLURM sbatch: 'sig_time must "
        "have an integer value'); got the MM:SS form"
    )
    assert not _has(script, _MMSS_SIGNAL), (
        "MM:SS --signal time ('INT@5:00') is not interpreted by SLURM as "
        "'5 minutes before end'; the supervisor would not checkpoint"
    )


def test_signal_zero_has_a_nonzero_checkpoint_window() -> None:
    """signal_minutes_before=0 must still produce a sane --signal time.

    Expected RED: the current code renders ``INT@0`` (signal at the walltime
    end, no checkpoint window) instead of refusing or computing a real window.
    """
    script = render_durable_job(signal_minutes_before=0)
    assert _SIGNAL_LINE in script
    # @0 fires at the deadline — contradicts "signal before the end"
    assert f"{_SIGNAL_LINE}INT@0" not in script


def test_export_all_exact() -> None:
    """--export=ALL exactly, never a named variable list (Target 4)."""
    script = render_durable_job()
    assert "#SBATCH --export=ALL\n" in script
    # no explicit variable export, e.g. --export=ALL,PATH=... or ~= separator
    assert "\n#SBATCH --export=ALL," not in script


def test_supervisor_command_is_immutable() -> None:
    """The exec line starts the same supervisor; no caller-visible command."""
    script = render_durable_job()
    assert "exec " in script
    assert " -m delfin.agent.service.supervisor" in script
    # no shell metacharacter injection surface in the exec line
    for bad in (";", "$(", "`", "|", "&&"):
        line = next((l for l in script.splitlines() if "exec " in l), "")
        assert bad not in line, f"supervisor exec line must not allow {bad!r}"


def test_newline_in_partition_refused() -> None:
    with pytest.raises(ValueError):
        render_durable_job(partition="cpu\n#SBATCH --time=00:00:01")


def test_newline_in_qos_refused() -> None:
    with pytest.raises(ValueError):
        render_durable_job(qos="normal\n#SBATCH --time=00:00:01")


def test_newline_in_walltime_refused() -> None:
    with pytest.raises(ValueError):
        render_durable_job(walltime="12:00:00\n#SBATCH --partition=gpu")


def test_unicode_lookalike_in_partition_refused() -> None:
    # full-width dash (U+FF0D) is not an ASCII name char
    with pytest.raises(ValueError):
        render_durable_job(partition="cpu\uFF0Dx")


def test_leading_dash_partition_refused() -> None:
    with pytest.raises(ValueError):
        render_durable_job(partition="-o/tmp/pwn")


def test_short_walltime_for_signal_is_rejected() -> None:
    """A walltime that cannot hold a signal window must raise, not render."""
    with pytest.raises(ValueError):
        render_durable_job(walltime="00:00:05", signal_minutes_before=5)


def test_default_walltime_valid() -> None:
    script = render_durable_job()
    assert "#SBATCH --time=" in script
