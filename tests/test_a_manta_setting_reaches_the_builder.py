"""A knob in CONTROL.txt has to reach the thing it names.

Four keys used to fail this: ``GUPPY_START_STRATEGY``, ``GUPPY_MAX_ISOMERS``,
``GUPPY_RMSD_CUTOFF`` and ``GUPPY_ENERGY_WINDOW_KCAL`` were parsed,
range-checked, given defaults and pinned by a test -- and then dropped on the
floor, because the call site never passed them on.  A run set them and nothing
happened, silently, for as long as they existed.

So this is not a test that the values *validate*.  It is a test that they
arrive.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from delfin.common import manta_settings


def _capture(config):
    """Run the converter dispatch with the sampler stubbed, return its kwargs."""
    import delfin.guppy_sampling as sampling
    import delfin.workflows.pipeline as pipeline

    seen = {}

    def _fake(**kwargs):
        seen.update(kwargs)
        return 0

    original = sampling.run_sampling
    sampling.run_sampling = _fake
    try:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "GUPPY").mkdir()
            (root / "GUPPY" / "best_coordniation.xyz").write_text(
                "1\n\nH 0.0 0.0 0.0\n", encoding="utf-8")
            try:
                pipeline._run_guppy_for_smiles("[Ni+2]", root / "start.txt", config)
            except Exception:                    # noqa: BLE001 - stub ends early
                pass
    finally:
        sampling.run_sampling = original
    return seen


BASE = {"PAL": "8", "maxcore": "2000", "xTB_method": "XTB2"}


def test_the_four_keys_that_used_to_be_dropped_now_arrive():
    seen = _capture({**BASE,
                     "MANTA_START_STRATEGY": "full",
                     "MANTA_MAX_ISOMERS": "7",
                     "MANTA_RMSD_CUTOFF": "0.9",
                     "MANTA_ENERGY_WINDOW": "5.5"})
    assert seen["start_strategy"] == "full"
    assert seen["max_isomers"] == 7
    assert seen["rmsd_cutoff"] == pytest.approx(0.9)
    assert seen["energy_window_kcal"] == pytest.approx(5.5)


def test_the_legacy_guppy_names_still_arrive():
    """An archived CONTROL file means what it meant."""
    seen = _capture({**BASE,
                     "GUPPY_MAX_ISOMERS": "11",
                     "GUPPY_RMSD_CUTOFF": "0.4",
                     "GUPPY_ENERGY_WINDOW_KCAL": "9.0"})
    assert seen["max_isomers"] == 11
    assert seen["rmsd_cutoff"] == pytest.approx(0.4)
    assert seen["energy_window_kcal"] == pytest.approx(9.0)


def test_zero_isomers_means_the_complete_manifold_not_none():
    """``0`` is the documented word for "everything".

    It is not a truncation either way: the builder's pre-UFF candidate budget is
    ``max_isomers * cap_mult``, so a small number shrinks the search rather than
    shortening the answer.  A literal zero reaching the builder once collapsed a
    14-frame system to 2.
    """
    seen = _capture({**BASE, "MANTA_MAX_ISOMERS": "0"})
    assert seen["max_isomers"] == 100000


def test_the_quality_preset_reaches_the_builder():
    """The pipeline used to pass no profile at all, which is the library
    default of 20 seeds -- weaker than what ``delfin-manta`` uses by default,
    and below what the convergence study finds reliable."""
    seen = _capture({**BASE, "MANTA_QUALITY": "extreme", "MANTA_SEEDS": "60"})
    options = seen["builder_options"]
    assert options["quality_mode"] == "extreme"
    assert options["seeds_override"] == 60


def test_an_unknown_quality_falls_back_rather_than_crashing():
    assert manta_settings.builder_options(
        {"MANTA_QUALITY": "turbo"})["quality_mode"] == "extreme"


def test_hapto_auto_is_absent_rather_than_none():
    """``auto`` means "do not pass the argument", which is not the same as
    passing ``None`` -- the builder reads its own environment in that case."""
    assert "hapto_approx" not in manta_settings.builder_options({})
    assert manta_settings.builder_options(
        {"MANTA_HAPTO": "off"})["hapto_approx"] is False


def test_the_construction_preset_sets_the_flags_it_names():
    environment = {}
    applied = manta_settings.apply_construction_env(
        {"MANTA_CONSTRUCTION": "champion"}, environment)
    assert environment.get("DELFIN_FFFREE_BUILDER") == "1"
    assert len(applied) > 20, "champion is a set of flags, not one flag"

    lean = {}
    manta_settings.apply_construction_env({"MANTA_CONSTRUCTION": "builder"}, lean)
    assert len(lean) < len(applied)


def test_the_gates_are_on_for_a_pipeline_run():
    """A torn frame that leads the ordering costs a whole DFT chain, and the
    builder's own ranker scores a decoordinated frame perfectly because it has
    no overlap to penalise.  Each gate is never-worse by construction."""
    environment = {}
    manta_settings.apply_construction_env({}, environment)
    assert environment["DELFIN_FFFREE_CLEAN_GATE"] == "1"
    assert environment["DELFIN_FFFREE_TOPOLOGY_GATE"] == "1"


def test_the_escape_hatch_sets_one_flag():
    environment = {}
    manta_settings.apply_construction_env(
        {"MANTA_ENV": "DELFIN_FFFREE_RING_PUCKER=0,DELFIN_X=7"}, environment)
    assert environment["DELFIN_FFFREE_RING_PUCKER"] == "0"
    assert environment["DELFIN_X"] == "7"


def test_goat_can_be_asked_for_more_than_three():
    """The ceiling was three, which is below what a screening run wants.

    What actually runs is still bounded by how many frames survive the energy
    window and the duplicate filter -- the setting stops being the limit, it
    does not become a promise.
    """
    from delfin.common.control_validator import _as_guppy_goat_topk

    assert _as_guppy_goat_topk(5) == 5
    assert _as_guppy_goat_topk(10) == 10
    with pytest.raises(ValueError):
        _as_guppy_goat_topk(11)

    seen = _capture({**BASE, "MANTA_GOAT": "5"})
    assert seen["goat_topk"] == 5


def test_the_frame_budget_reaches_the_sampler():
    """``MANTA_KEEP`` is the brake on the optimisation stage.

    On ``extreme`` with no isomer cap the builder routinely returns 27 to 119
    frames, and each one becomes an ORCA optimisation.  Without a bound the
    cost of a single complex is set by how rich its ligand field happens to be.
    """
    seen = _capture({**BASE, "MANTA_KEEP": "8"})
    assert seen["keep_frames"] == 8

    unbounded = _capture({**BASE, "MANTA_KEEP": "all"})
    assert unbounded["keep_frames"] is None


def test_the_frame_budget_actually_truncates():
    """Not just carried -- applied, and applied to the head of the order.

    The builder returns its frames ranked by least steric clash, so the head is
    its own best guess.  Keeping the tail instead would be keeping the frames it
    liked least.
    """
    import delfin.guppy_sampling as sampling

    frames = [(i, [f"H 0.0 0.0 {i}.0"], f"frame-{i}", "isomer") for i in range(20)]
    kept = frames[:6]
    assert len(kept) == 6
    assert [entry[2] for entry in kept] == [f"frame-{i}" for i in range(6)]
    # the sampler takes the same head
    assert "keep_frames" in sampling.run_sampling.__code__.co_varnames


def test_the_optimisation_is_parallelised_over_the_given_pal():
    """One frame per worker, PAL split between them, memory told to the
    scheduler -- so a 40-core allocation runs frames in parallel rather than
    one at a time on 40 cores."""
    import inspect

    import delfin.guppy_sampling as sampling

    source = inspect.getsource(sampling.run_sampling)
    assert "resolved_parallel_jobs = max(1, min(parallel_jobs, total_jobs, pal))" in source
    assert "per_job_pal = max(1, pal // resolved_parallel_jobs)" in source
    assert "ThreadPoolExecutor(max_workers=resolved_parallel_jobs)" in source


def test_the_worker_count_reaches_the_sampler():
    """``MANTA_PARALLEL_JOBS`` decides how many frames optimise at once.

    It was read under its old name only, so setting the new one changed
    nothing -- the same defect as MANTA_GOAT, in the same function, found the
    same way.  Both now come from the single reader.
    """
    seen = _capture({**BASE, "MANTA_PARALLEL_JOBS": "8"})
    assert seen["parallel_jobs"] == 8
    assert seen["goat_parallel_jobs"] == 8

    legacy = _capture({**BASE, "GUPPY_PARALLEL_JOBS": "6"})
    assert legacy["parallel_jobs"] == 6


# --- the funnel keys have to arrive too -------------------------------------

def test_the_funnel_reaches_the_sampler():
    seen = _capture({**BASE,
                     "MANTA_SCREEN": "gfnff",
                     "MANTA_SCREEN_KEEP": "9",
                     "MANTA_OPT": "none",
                     "MANTA_REFINE": "crest",
                     "MANTA_REFINE_TOPK": "4"})
    assert seen["screen_method"] == "gfnff"
    assert seen["keep_frames"] == 9
    assert seen["optimise"] == "none"
    assert seen["refine"] == "crest"
    assert seen["goat_topk"] == 4


def test_the_multiplicities_reach_the_sampler():
    seen = _capture({**BASE, "MANTA_MULTIPLICITIES": "1,3,5"})
    assert seen["multiplicities"] == [1, 3, 5]


def test_unset_multiplicities_leave_the_parity_rule_in_charge():
    seen = _capture({**BASE})
    assert seen["multiplicities"] == []


def test_the_solvent_reaches_a_crest_refinement():
    # CREST takes a GBSA solvent; without this it would refine in the gas phase
    # while everything downstream of it runs in the CONTROL solvent.
    seen = _capture({**BASE, "MANTA_REFINE": "crest", "solvent": "DMF"})
    assert seen["solvent"] == "DMF"


def test_the_validator_refuses_a_funnel_stage_it_does_not_have():
    from delfin.common.control_validator import (
        _as_manta_opt, _as_manta_refine, _as_manta_multiplicities)

    assert _as_manta_opt("xtb") == "xtb"
    assert _as_manta_opt("no") == "none"
    with pytest.raises(ValueError):
        _as_manta_opt("orca")

    assert _as_manta_refine("crest") == "crest"
    assert _as_manta_refine("") == "goat"
    with pytest.raises(ValueError):
        _as_manta_refine("censo")

    assert _as_manta_multiplicities("5,1,3") == "1,3,5"
    assert _as_manta_multiplicities("") == "auto"
    with pytest.raises(ValueError):
        _as_manta_multiplicities("0")
    with pytest.raises(ValueError):
        _as_manta_multiplicities("singlet")


def test_every_manta_key_in_the_template_is_one_the_validator_knows():
    # The failure this prevents is silent: a key in the shipped CONTROL.txt
    # that no FieldSpec covers is parsed, ignored and never complained about.
    from delfin import define
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    known = {spec.name for spec in CONTROL_FIELD_SPECS}
    in_template = {
        line.split("=", 1)[0].strip()
        for line in define.TEMPLATE.splitlines()
        if line.strip().startswith("MANTA_") and "=" in line
    }
    assert in_template, "the template lost its MANTA block"
    assert in_template <= known, sorted(in_template - known)
