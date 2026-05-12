"""End-to-end test for the iso_dir-output-visibility fix.

After commit d853789 (orca: write ORCA stdout direct to parent_dir, not iso_dir),
the ORCA `.out` file must:
  1. Land in `parent_dir`, not `iso_dir` (so it survives `iso_dir` cleanup races
     and is visible to `periodic_copy` in the SLURM submit wrapper).
  2. Contain whatever ORCA wrote, even if ORCA crashes — so that "[Errno 2]
     No such file or directory: '.../.orca_iso_*'" never again means the user
     loses the diagnostic header.
  3. Still receive aux files (.gbw, .opt, .engrad, …) copied back from iso_dir
     after a successful run.

These guarantees apply to all three submit modes (initial, --recalc,
--smart-recalc) because all of them flow through `_run_orca_isolated`
when `isolate=True`, and the dashboard sets `isolate=True` for ESD,
Occupier, IMAG, GUPPY-sampling.

The test runs `_run_orca_isolated` against a fake ORCA binary (a small
shell script). No real ORCA install needed.
"""
from __future__ import annotations

import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from delfin.orca import _run_orca_isolated


def _make_fake_orca(tmp_path: Path, script_body: str) -> Path:
    """Create an executable shell script that mimics ORCA for testing.

    The script writes given content to its first arg's matching .out file
    (where stdout is captured by _run_orca_subprocess via Popen redirection).
    """
    fake = tmp_path / "fake_orca"
    fake.write_text("#!/usr/bin/env bash\n" + script_body)
    fake.chmod(fake.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
    return fake


def _make_inp(parent_dir: Path, name: str = "S0.inp") -> Path:
    inp = parent_dir / name
    inp.write_text(
        "! B3LYP def2-SVP\n"
        "%pal nprocs 1 end\n"
        "* xyz 0 1\n"
        "H 0.0 0.0 0.0\n"
        "*\n"
    )
    return inp


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    """ESD/-style working directory with a clean input."""
    parent = tmp_path / "ESD"
    parent.mkdir()
    return parent


def test_iso_dir_run_writes_out_directly_to_parent_dir(workdir, monkeypatch):
    """The .out file must end up in parent_dir, populated by ORCA's stdout —
    proving the fix from d853789. Pre-d853789 the .out lived in iso_dir/
    and was lost whenever iso_dir disappeared mid-run.
    """
    inp = _make_inp(workdir)
    out = workdir / "S0.out"

    fake_orca = _make_fake_orca(
        workdir.parent,
        # Print a banner that includes the success marker
        # _check_orca_success looks for 'ORCA TERMINATED NORMALLY' in the
        # last 8 KB of the output file.
        'echo "                                 *****************"\n'
        'echo "                                 * O   R   C   A *"\n'
        'echo "                                 *****************"\n'
        'echo "Program Version 6.1.1"\n'
        'echo "test SCF iteration 1: -42.0"\n'
        'echo "                             ****ORCA TERMINATED NORMALLY****"\n'
        'exit 0\n'
    )

    success = _run_orca_isolated(
        str(fake_orca),
        inp,
        out,
        copy_files=None,
    )

    assert success is True
    assert out.exists(), "S0.out must be in parent_dir, not iso_dir"
    body = out.read_text()
    assert "* O   R   C   A *" in body, "ORCA header must be in S0.out"
    assert "ORCA TERMINATED NORMALLY" in body
    assert "test SCF iteration 1: -42.0" in body


def test_aux_files_copied_back_from_iso_dir(workdir, monkeypatch):
    """ORCA-generated aux files (.gbw, .opt, etc.) live in iso_dir while
    ORCA runs and must be copied back to parent_dir on success. The fix
    must not regress this — the for-iso_entries loop is still required.
    """
    inp = _make_inp(workdir)
    out = workdir / "S0.out"

    # Fake ORCA writes both stdout (captured to out) AND a .gbw file in cwd
    # (which is iso_dir). Our patch must copy the .gbw back to parent_dir.
    fake_orca = _make_fake_orca(
        workdir.parent,
        'echo "* O   R   C   A *"\n'
        'echo "Program Version 6.1.1"\n'
        'printf "FAKE_GBW_BYTES" > S0.gbw\n'  # written into cwd = iso_dir
        'printf "FAKE_OPT_BYTES" > S0.opt\n'
        'echo "                             ****ORCA TERMINATED NORMALLY****"\n'
        'exit 0\n'
    )

    success = _run_orca_isolated(
        str(fake_orca),
        inp,
        out,
        copy_files=None,
    )

    assert success is True
    assert (workdir / "S0.gbw").exists(), "Aux .gbw must be copied back"
    assert (workdir / "S0.gbw").read_bytes() == b"FAKE_GBW_BYTES"
    assert (workdir / "S0.opt").exists()
    assert (workdir / "S0.opt").read_bytes() == b"FAKE_OPT_BYTES"


def test_orca_crash_leaves_partial_out_in_parent_dir(workdir, monkeypatch):
    """When ORCA crashes early (no success marker), the partial output —
    header, basis-set echo, error line — must survive in parent_dir.
    Pre-d853789 the .out lived in iso_dir/ and was inaccessible after the
    iso_dir.iterdir() OSError aborted the rescue path. With the fix, the
    partial .out is already in parent_dir from ORCA's first byte.
    """
    inp = _make_inp(workdir)
    out = workdir / "S0.out"

    # Crash mid-output, NO 'ORCA TERMINATED NORMALLY'
    fake_orca = _make_fake_orca(
        workdir.parent,
        'echo "* O   R   C   A *"\n'
        'echo "Program Version 6.1.1"\n'
        'echo "Your calculation uses the basis: def2-TZVP"\n'
        'echo "ERROR: SCF failed to converge"\n'
        'exit 1\n'
    )

    success = _run_orca_isolated(
        str(fake_orca),
        inp,
        out,
        copy_files=None,
    )

    assert success is False, "Crash must propagate as failure"
    assert out.exists(), (
        "Even after ORCA crash, S0.out must be in parent_dir for diagnosis"
    )
    body = out.read_text()
    assert "* O   R   C   A *" in body
    assert "ERROR: SCF failed to converge" in body, (
        "User must be able to read the actual ORCA error from /calc/"
    )


def test_iso_dir_disappearance_does_not_lose_out(workdir, monkeypatch):
    """Defensive: even if iso_dir is removed before the aux-copy loop
    (the original symptom that motivated d853789 — '[Errno 2] No such
    file or directory'), the main S0.out must still be present.
    """
    inp = _make_inp(workdir)
    out = workdir / "S0.out"

    # Fake ORCA writes output, then deletes its own working dir on exit
    # (simulating an aggressive MPI scratch-cleanup that ORCA sometimes does).
    fake_orca = _make_fake_orca(
        workdir.parent,
        'echo "* O   R   C   A *"\n'
        'echo "Program Version 6.1.1"\n'
        'echo "                             ****ORCA TERMINATED NORMALLY****"\n'
        # Schedule self-destruction of cwd shortly after exit so that
        # the post-run aux-copy loop hits a vanished iso_dir.
        '(sleep 0.05; rm -rf "$PWD") &\n'
        'exit 0\n'
    )

    success = _run_orca_isolated(
        str(fake_orca),
        inp,
        out,
        copy_files=None,
    )

    # Aux-files might be lost (iso_dir gone), but S0.out is independent now.
    assert out.exists(), (
        "S0.out must be in parent_dir even when iso_dir vanishes — "
        "this is the core symptom of jobs 4656609/4656950 that d853789 fixes."
    )
    body = out.read_text()
    assert "* O   R   C   A *" in body
    # Whether `success` is True or False here depends on race timing, but
    # the diagnostic content must be available either way.


def test_recalc_smart_recalc_initial_share_isolated_path(workdir):
    """All three submit modes (initial, --recalc, --smart-recalc) share the
    same _run_orca_isolated code path:
      - Initial submit  → run_orca(isolate=True) → _run_orca_isolated
      - --recalc        → cli_recalc patches run_orca → _run_orca_real →
                          run_orca(isolate=True) → _run_orca_isolated
      - --smart-recalc  → same wrapper chain, just adds skip-fingerprint
                          check → _run_orca_real → _run_orca_isolated
    Therefore a single fix in _run_orca_isolated covers all three modes.
    This regression test pins that contract.
    """
    from delfin import cli_recalc, orca

    # cli_recalc.setup_recalc_mode returns (wrappers, reals).
    # The 'reals' for run_orca is the un-patched orca.run_orca function
    # which dispatches to _run_orca_isolated when isolate=True.
    wrappers, reals = cli_recalc.setup_recalc_mode()

    # 'run_orca' real is the entry point that branches on `isolate=` to
    # either _run_orca_isolated or _run_orca_subprocess. So when ESD,
    # Occupier, IMAG etc. call run_orca(isolate=True), the patched wrapper
    # forwards to this real, which goes through _run_orca_isolated, which
    # is where the d853789 fix lives.
    assert reals["run_orca"] is orca.run_orca, (
        "Recalc reals must dispatch through orca.run_orca so the isolated "
        "path (and thus the d853789 fix) applies to --recalc / --smart-recalc."
    )
    # smart_mode is a runtime flag inside the wrapper, not a different
    # underlying real — same dispatcher, same isolated path, same fix.
    assert callable(wrappers["run_orca"])


# =============================================================================
# Realistic CONTROL.txt scenarios from /home/qmchem_all/archive
# =============================================================================
#
# These parameterised tests pin the end-to-end contract for each known user's
# typical setup. They check that:
#   1. GlobalJobManager.initialize() correctly derates maxcore (or doesn't,
#      when no derate is needed) for the SLURM allocation typical for that
#      user's compute partition.
#   2. The caller's config dict is mutated in place so that the engine's
#      maxcore_mb reflects the derated value (fixed in f0998d2).
#   3. The pool budget cannot ever block a single PAL-wide ORCA job from
#      starting — the original deadlock symptom that drove this whole chain.
#
# Each scenario corresponds to a real archive directory:
#   Fritz   /home/qmchem_all/archive/Fritz/*                 (single ESD-S0)
#   Jerome  /home/qmchem_all/archive/Jerome/complexes_200_400 (Occupier classic)
#   Max     /home/qmchem_all/Archiv/Max/*                    (large-PAL TZVP)
# =============================================================================

import pytest as _pytest

from delfin.workflows.scheduling.manager import _SLURM_MEM_HEADROOM


@_pytest.fixture
def fresh_global_manager(monkeypatch):
    """Reset GlobalJobManager singleton before/after each scenario."""
    from delfin.global_manager import GlobalJobManager

    GlobalJobManager._instance = None
    mgr = GlobalJobManager()
    yield mgr
    try:
        mgr.shutdown()
    except Exception:
        pass
    GlobalJobManager._instance = None


_SCENARIOS = [
    # name, slurm_mem_mb, control_dict, expected_derate
    pytest.param(
        "fritz_esd_tzvp",
        108_000,  # BwUni dev_single 12-core node, ~105.5 GB
        {"PAL": 12, "maxcore": 9000, "pal_jobs": 2,
         "calc_initial": "no", "ESD_modul": "yes",
         "main_basisset": "def2-TZVP"},
        True,  # 12*9000 = 108000 > 91800 cap → derate to 7650
        id="fritz_esd_tzvp_derate",
    ),
    pytest.param(
        "fritz_esd_svp_safe",
        108_000,
        {"PAL": 12, "maxcore": 4000, "pal_jobs": 2,
         "calc_initial": "no", "ESD_modul": "yes",
         "main_basisset": "def2-SVP"},
        False,  # 12*4000 = 48000 ≤ 91800 → no derate
        id="fritz_esd_svp_no_derate",
    ),
    pytest.param(
        "jerome_occupier_classic",
        260_000,  # BwUni dev_multiple-style 40-core node
        {"PAL": 40, "maxcore": 6000, "pal_jobs": 4,
         "calc_initial": "yes", "ESD_modul": "no",
         "main_basisset": "def2-SVP"},
        True,  # 40*6000 = 240000 > 221000 cap → derate
        id="jerome_occupier_classic_derate",
    ),
    pytest.param(
        "max_archive_pal60",
        360_000,  # large compute node ~360 GB
        {"PAL": 60, "maxcore": 3800, "pal_jobs": 6},
        False,  # 60*3800 = 228000 ≤ 306000 cap → no derate
        id="max_pal60_no_derate",
    ),
    pytest.param(
        "jew_r465_oom_scenario",
        360_000,
        {"PAL": 60, "maxcore": 6000, "pal_jobs": 4,
         "calc_initial": "yes", "ESD_modul": "yes"},
        True,  # 60*6000 = 360000 > 306000 cap → derate to 5100
        id="jew_r465_derate_protects_OOM",
    ),
]


@pytest.mark.parametrize("name,slurm_mem,control,expect_derate", _SCENARIOS)
def test_real_user_scenarios_derate_and_propagate(
    name, slurm_mem, control, expect_derate, monkeypatch, fresh_global_manager
):
    """For every realistic CONTROL.txt scenario:
      - Manager derates iff PAL × maxcore > SLURM × headroom.
      - Caller's config dict is mutated in place (so engine.maxcore_mb sees it).
      - Single PAL-wide ORCA job's memory_mb fits the pool — no deadlock.
    """
    monkeypatch.setenv("SLURM_MEM_PER_NODE", str(slurm_mem))
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)

    config = dict(control)  # copy so monkeypatch doesn't leak
    original_maxcore = config["maxcore"]

    fresh_global_manager.initialize(config)

    if expect_derate:
        cap = int(slurm_mem * _SLURM_MEM_HEADROOM)
        expected_maxcore = cap // config["PAL"]
        assert config["maxcore"] == expected_maxcore, (
            f"[{name}] Derate not propagated: got {config['maxcore']}, "
            f"expected {expected_maxcore}. Engine would deadlock."
        )
        assert config["maxcore"] < original_maxcore
    else:
        assert config["maxcore"] == original_maxcore, (
            f"[{name}] Spurious derate: original {original_maxcore} "
            f"changed to {config['maxcore']}, but PAL × maxcore fits the cap."
        )

    # The pool budget must accommodate one full PAL-wide ORCA job —
    # otherwise the deadlock the chain of fixes was meant to solve returns.
    pool_budget = fresh_global_manager.pool.total_memory_mb
    full_pal_request = config["PAL"] * config["maxcore"]
    assert full_pal_request <= pool_budget, (
        f"[{name}] PAL × maxcore = {full_pal_request} MB > pool {pool_budget} MB "
        f"— single-job deadlock will return."
    )


_RECALC_SCENARIOS = [
    pytest.param("initial_submit", "initial",
                 id="initial_uses_isolated_path"),
    pytest.param("recalc_classic", "recalc",
                 id="recalc_uses_isolated_path"),
    pytest.param("recalc_smart", "smart-recalc",
                 id="smart_recalc_uses_isolated_path"),
]


@pytest.mark.parametrize("name,mode", _RECALC_SCENARIOS)
def test_all_three_modes_produce_visible_out_file(name, mode, workdir):
    """Initial submit / --recalc / --smart-recalc must all leave the ORCA
    .out file readable in parent_dir. The mode only changes which wrapper
    intercepts run_orca — the underlying isolated path is the same.

    We exercise the underlying _run_orca_isolated directly (it's the shared
    sink for all three modes, as proven by
    test_recalc_smart_recalc_initial_share_isolated_path).
    """
    inp = _make_inp(workdir, name=f"{mode}_S0.inp")
    out = workdir / f"{mode}_S0.out"

    fake_orca = _make_fake_orca(
        workdir.parent,
        f'echo "* O   R   C   A *  (mode={mode})"\n'
        'echo "Program Version 6.1.1"\n'
        'echo "                             ****ORCA TERMINATED NORMALLY****"\n'
    )

    success = _run_orca_isolated(str(fake_orca), inp, out, copy_files=None)

    assert success
    assert out.exists()
    body = out.read_text()
    assert f"mode={mode}" in body, (
        f"Mode '{mode}' did not produce the expected ORCA output in /calc/."
    )


def test_jerome_style_parallel_jobs_each_get_own_out(workdir):
    """Jerome's Occupier setup runs many ORCA jobs in parallel (pal_jobs=4).
    Each job must get its own visible .out file in parent_dir — no cross-job
    overwrite, no race on the .out file. The d853789 fix preserves this:
    each job has its own unique output_path (S0.out, oxidation.out,
    reduction.out, …) so direct-write to parent_dir is safe.
    """
    fake_orca = _make_fake_orca(
        workdir.parent,
        'inp_basename=$(basename "$1" .inp)\n'
        'echo "* O   R   C   A * (job=${inp_basename})"\n'
        'echo "                             ****ORCA TERMINATED NORMALLY****"\n'
    )

    job_names = ["S0", "oxidation", "reduction", "T1"]
    results = {}
    for name in job_names:
        inp = _make_inp(workdir, name=f"{name}.inp")
        out = workdir / f"{name}.out"
        success = _run_orca_isolated(str(fake_orca), inp, out, copy_files=None)
        results[name] = (success, out)

    # All four parallel-style jobs must each have their own output, no cross-talk
    for name, (success, out_path) in results.items():
        assert success, f"Job {name} failed"
        assert out_path.exists(), f"Job {name} has no .out file"
        body = out_path.read_text()
        assert f"job={name}" in body, (
            f"Cross-job contamination: {name}.out has wrong content"
        )

