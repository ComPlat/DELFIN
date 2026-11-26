"""`delfin stop --cleanup` must keep what a resumed stage reads, and remove only scratch.

The ``*_OCCUPIER`` stage directories are resume state, not scratch. A resumed
OCCUPIER stage sets itself up from them:

* ``read_occ`` (delfin/workflows/engine/occupier.py:646) raises ``RuntimeError``
  when ``<stage>/OCCUPIER.txt`` is missing, and its caller at :738 then falls
  back to ``read_occ_from_control`` -- a different electron configuration than
  the one the stage chose.
* ``extract_preferred_spin`` (delfin/copy_helpers.py:153) reads
  ``<stage>/OCCUPIER.txt`` for the preferred multiplicity; with the folder gone
  it returns ``(None, None)`` and ``cli._sync_occupier_spin_metadata``
  (delfin/cli.py:1470) silently keeps CONTROL's multiplicity instead.
* the per-stage ``CONTROL.txt`` (delfin/occupier_sequences.py:305) and the
  ``.gbw``/``.xyz`` a stage restarts from live there too.

``delfin stop --cleanup`` and ``delfin cleanup --orca`` both reach
``cleanup_orca`` -> ``_remove_occuper_directories``, which used to delete those
folders outright. These tests drive the real CLI entry point ``cli.main``.

What must still go is true scratch: the ``DEFAULT_PATTERNS`` temporary files
(``*.cpcm``, ``*.tmp*`` ...) and ``.orca_scratch`` -- including the copies that
sit *inside* a kept stage directory, because ``cleanup_all`` recurses.
"""
from pathlib import Path

import pytest

# A real OCCUPIER.txt: the winning FSPE entry is tagged so that
# extract_preferred_spin and read_occupier_file find the preferred spin.
OCCUPIER_TXT = """\
Method: PBE0 def2-SVP
Charge: 0
-------------
FINAL SINGLE POINT ENERGY (1)       = -2609.340 (H)
multiplicity 1
----------------------------------------------------------------
FINAL SINGLE POINT ENERGY (2)       = -2609.352 (H) <-- PREFERRED VALUE
multiplicity 3, BrokenSym 1,1
----------------------------------------------------------------
Preferred Index: 2
"""


def _fake_proc(workspace: Path):
    """A recorded-but-long-gone DELFIN entry, so `stop` has something to signal."""
    return [{
        "pid": 4_999_999,
        "pgid": 4_999_999,
        "cmd": "delfin pipeline",
        "cwd": str(workspace),
    }]


@pytest.fixture
def run_cli(monkeypatch, tmp_path):
    """Run delfin.cli.main() against tmp, with the process scans faked."""
    def _run(argv):
        monkeypatch.chdir(tmp_path)
        monkeypatch.setattr(
            "delfin.cli._iter_delfin_procs", lambda ws: _fake_proc(Path(ws)))
        # No `ps` dependency and no real signal: the ORCA process scan finds
        # nothing, and os.kill on the fake pid raises ProcessLookupError, which
        # `stop` already handles as "already gone".
        import subprocess

        monkeypatch.setattr(
            "delfin.cleanup.subprocess.run",
            lambda *a, **k: subprocess.CompletedProcess(
                args=[], returncode=1, stdout="", stderr=""))
        from delfin import cli
        return cli.main(argv)
    return _run


def make_workspace(root: Path) -> Path:
    """The workspace a stopped OCCUPIER run leaves: resume state + true scratch."""
    ws = root / "ws"
    ws.mkdir()
    (ws / "CONTROL.txt").write_text("STAGE=0\n", encoding="utf-8")

    # Resume state: three stage folders, each with what the stage reads.
    initial = ws / "initial_OCCUPIER"
    initial.mkdir()
    (initial / "CONTROL.txt").write_text("STAGE=initial\n", encoding="utf-8")
    (initial / "OCCUPIER.txt").write_text(OCCUPIER_TXT, encoding="utf-8")
    (initial / "initial.xyz").write_text("2\nfrom opt\nC 0 0 0\nH 0 0 1\n",
                                        encoding="utf-8")

    ox1 = ws / "ox_step_1_OCCUPIER"
    ox1.mkdir()
    (ox1 / "OCCUPIER.txt").write_text(OCCUPIER_TXT, encoding="utf-8")
    (ox1 / "input_ox_step_1_OCCUPIER.gbw").write_bytes(b"fake gbw")
    (ox1 / "input_ox_step_1_OCCUPIER.xyz").write_text(
        "2\nfrom opt\nC 0 0 0\nH 0 0 1\n", encoding="utf-8")

    red2 = ws / "red_step_2_OCCUPIER"
    red2.mkdir()
    (red2 / "OCCUPIER.txt").write_text(OCCUPIER_TXT, encoding="utf-8")
    (red2 / "input_red_step_2_OCCUPIER.gbw").write_bytes(b"fake gbw")

    # True scratch at the workspace root.
    (ws / "leftover.cpcm").write_text("scratch cpcm\n", encoding="utf-8")
    (ws / ".orca_scratch").mkdir()
    (ws / ".orca_scratch" / "job.gbw").write_bytes(b"scratch")

    # True scratch *inside* a kept stage: cleanup_all recurses, so a fix that
    # only spares the directory but leaves the recursion unguarded deletes these.
    (ox1 / "output.gbw.tmp").write_text("tmp\n", encoding="utf-8")
    (initial / "run.cpcm").write_text("cpcm\n", encoding="utf-8")

    return ws


def test_stop_cleanup_keeps_what_a_resume_reads(run_cli, tmp_path):
    """Red on main: stop --cleanup deletes the stage folders a resume reads."""
    ws = make_workspace(tmp_path)
    run_cli(["stop", "--workspace", str(ws), "--cleanup"])

    # The stage folder and the OCCUPIER.txt a resumed stage sets itself up from.
    for stage, needed in (
        ("initial_OCCUPIER", ["OCCUPIER.txt", "CONTROL.txt", "initial.xyz"]),
        ("ox_step_1_OCCUPIER",
         ["OCCUPIER.txt", "input_ox_step_1_OCCUPIER.gbw",
          "input_ox_step_1_OCCUPIER.xyz"]),
        ("red_step_2_OCCUPIER", ["OCCUPIER.txt", "input_red_step_2_OCCUPIER.gbw"]),
    ):
        assert (ws / stage).is_dir(), f"{stage} is resume state, not scratch"
        for name in needed:
            assert (ws / stage / name).is_file(), (
                f"{stage}/{name} is what a resumed stage reads")


def test_the_preferred_spin_is_readable_after_stop_cleanup(run_cli, tmp_path):
    """Red on main: with the folder gone the preferred spin cannot be read.

    This is the consequence the finding names: ``extract_preferred_spin``
    answers ``(None, None)`` and ``_sync_occupier_spin_metadata`` keeps
    CONTROL's multiplicity instead of the one OCCUPIER chose.
    """
    from delfin.copy_helpers import extract_preferred_spin

    ws = make_workspace(tmp_path)
    mult_before, bs_before = extract_preferred_spin(ws / "initial_OCCUPIER")
    assert (mult_before, bs_before) == (3, "1,1"), "fixture must be readable"

    run_cli(["stop", "--workspace", str(ws), "--cleanup"])

    mult_after, bs_after = extract_preferred_spin(ws / "initial_OCCUPIER")
    assert (mult_after, bs_after) == (3, "1,1"), (
        "the preferred spin OCCUPIER chose must survive stop --cleanup")


def test_cleanup_orca_dry_run_reports_no_occupier_removals(run_cli, tmp_path, capsys):
    """Red on main: the dry run lists the stage folders as removed."""
    ws = make_workspace(tmp_path)
    run_cli(["cleanup", "--workspace", str(ws), "--orca", "--dry-run"])
    out = capsys.readouterr().out
    assert "Removed OCCUPIER folders" not in out, (
        "a dry run must not promise to delete resume state:\n" + out)
    # The stage folders may be NAMED (as kept) but never as removed.
    assert "Kept OCCUPIER folders" in out, (
        "the report should say what it spared, so the operator knows why:\n" + out)


def test_a_stage_survives_cleanup_as_resumable_state(run_cli, tmp_path):
    """Red on main: after stop --cleanup the pipeline cannot reuse the stage.

    ``_occupier_results_fresh`` (delfin/workflows/pipeline.py:139) is the gate
    that decides whether a recalc continues from the existing OCCUPIER results
    or starts the sequence over. It fails closed to a rerun -- a new set of
    FSPE calculations -- as soon as ``OCCUPIER.txt`` or the control
    fingerprint is missing, which is exactly what the cleanup left behind.
    """
    from delfin.workflows.pipeline import _stamp_occupier_control, _occupier_results_fresh

    ws = make_workspace(tmp_path)
    # The state a finished stage leaves: results stamped with the run CONTROL.
    _stamp_occupier_control(ws / "CONTROL.txt", ws / "initial_OCCUPIER")
    assert _occupier_results_fresh(ws / "CONTROL.txt", ws / "initial_OCCUPIER"), (
        "fixture must be reusable before the cleanup")

    run_cli(["stop", "--workspace", str(ws), "--cleanup"])

    assert _occupier_results_fresh(ws / "CONTROL.txt", ws / "initial_OCCUPIER"), (
        "stop --cleanup must not force the OCCUPIER sequence to start over")


def test_stop_cleanup_still_removes_true_scratch(run_cli, tmp_path):
    """Control, green on both: genuine scratch is still removed."""
    ws = make_workspace(tmp_path)
    run_cli(["stop", "--workspace", str(ws), "--cleanup"])

    assert not (ws / ".orca_scratch").exists(), ".orca_scratch is true scratch"
    assert not (ws / "leftover.cpcm").exists(), "*.cpcm is true scratch"
    # the copies inside a kept stage folder
    assert not (ws / "ox_step_1_OCCUPIER" / "output.gbw.tmp").exists(), (
        "*.tmp inside a kept stage is still true scratch")
    assert not (ws / "initial_OCCUPIER" / "run.cpcm").exists(), (
        "*.cpcm inside a kept stage is still true scratch")
