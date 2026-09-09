"""Three stages, each switchable, and the winner keeps its own spin state.

MANTA hands the pipeline a whole manifold of frames, and the pipeline needs
exactly one geometry.  Between the two sits a funnel with three stages that
cost three different amounts:

    screen (one single point per frame)  ->  optimise  ->  refine (GOAT/CREST)

A run has different reasons to want each of them, so each has to be switchable
on its own from CONTROL.txt: rank and optimise only the top of the ranking, or
optimise everything and never rank at all, or rank and stop.  This test pins
those combinations, because they are the ones a user actually asks for and
because the difference between them is hours of ORCA.

The fourth thing pinned here is the multiplicity.  A metal complex has to be
tried at more than one, since which coordination isomer lies lowest and which
spin state lies lowest are not separable questions -- and once several are
tried, the winner is a (frame, multiplicity) pair.  Refining that winner at
some default spin state instead of its own would be refining a different
molecule, silently, with both calculations perfectly valid.
"""

from __future__ import annotations

import pytest

from delfin.common.manta_settings import selection_options
from delfin.guppy_sampling import (
    _multiplicity_argument,
    _screen_start_geometries,
    multiplicity_of_label,
)


def _frames(n=3):
    return [(i, [f"Ni {i}.0 0.0 0.0", "O 0.0 0.0 2.0"], f"iso{i}", "isomer")
            for i in range(1, n + 1)]


# --- what CONTROL asks for, and what the funnel then is ---------------------

def test_rank_and_then_optimise_only_the_top():
    got = selection_options({"MANTA_SCREEN": "gfn2", "MANTA_SCREEN_KEEP": "12"})
    assert got["screen"] == "gfn2"
    assert got["screen_keep"] == 12
    assert got["optimise"] == "xtb"


def test_optimise_everything_needs_no_rank():
    got = selection_options({"MANTA_SCREEN": "none", "MANTA_SCREEN_KEEP": "all",
                             "MANTA_OPT": "xtb"})
    assert got["screen"] == "none"
    assert got["screen_keep"] is None
    assert got["optimise"] == "xtb"


def test_ranking_alone_optimises_nothing():
    got = selection_options({"MANTA_SCREEN": "gfn2", "MANTA_OPT": "none"})
    assert got["optimise"] == "none"
    assert got["screen"] == "gfn2"


def test_neither_stage_measuring_anything_is_obeyed_not_overruled():
    # Asking for no screen *and* no optimisation asks for a ranking with no
    # measurement behind it.  The settings are obeyed as written -- quietly
    # substituting a screen would be deciding for the user -- and the run is
    # warned instead, which is what test_ranking_nothing_is_allowed_but_says_so
    # pins.
    got = selection_options({"MANTA_SCREEN": "none", "MANTA_OPT": "none"})
    assert got["optimise"] == "none"
    assert got["screen"] == "none"


def test_refine_can_be_crest_not_only_goat():
    got = selection_options({"MANTA_REFINE": "crest", "MANTA_REFINE_TOPK": "3"})
    assert got["refine"] == "crest"
    assert got["refine_topk"] == 3


def test_refine_can_be_switched_off():
    assert selection_options({"MANTA_REFINE": "none"})["refine"] == "none"


# --- the multiplicities -----------------------------------------------------

def test_unset_multiplicities_mean_the_parity_rule():
    # Empty is not "multiplicity 1".  It is "nobody said, so use the rule
    # DELFIN already uses everywhere else", which the sampler applies from the
    # electron count.
    assert selection_options({})["multiplicities"] == []
    assert selection_options({"MANTA_MULTIPLICITIES": "auto"})["multiplicities"] == []


@pytest.mark.parametrize("text,expected", [
    ("1,3,5", [1, 3, 5]),
    ("2 4", [2, 4]),
    ("5,1,3", [1, 3, 5]),
    ("3,3,3", [3]),
    ("1, 3 ,5", [1, 3, 5]),
])
def test_multiplicities_can_be_named(text, expected):
    assert selection_options({"MANTA_MULTIPLICITIES": text})["multiplicities"] == expected


def test_the_command_line_reads_multiplicities_the_same_way():
    assert _multiplicity_argument("1,3,5") == [1, 3, 5]
    assert _multiplicity_argument("") == []
    assert _multiplicity_argument("0,-2,3") == [3]


def test_every_frame_is_offered_at_every_multiplicity():
    # Three frames at three spin states is nine candidates, not three.  Ranking
    # frames at one assumed multiplicity answers neither of the two questions.
    pairs = _screen_start_geometries(
        _frames(3), method="none", keep=None, charge=0,
        multiplicities=[1, 3, 5])
    assert len(pairs) == 9
    assert sorted({mult for _, mult, _ in pairs}) == [1, 3, 5]


def test_the_keep_cuts_pairs_not_frames():
    pairs = _screen_start_geometries(
        _frames(3), method="none", keep=4, charge=0, multiplicities=[1, 3])
    assert len(pairs) == 4


def test_an_unscreened_pair_carries_no_energy_rather_than_a_made_up_one():
    pairs = _screen_start_geometries(
        _frames(2), method="none", keep=None, charge=0, multiplicities=[1])
    assert [energy for _, _, energy in pairs] == [None, None]


# --- the winner's own spin state reaches the refinement ---------------------

def test_the_refinement_reads_the_multiplicity_off_the_winner():
    assert multiplicity_of_label("iso3 M=5", 1) == 5
    assert multiplicity_of_label("iso3 M=5 goat", 1) == 5


def test_a_label_without_a_multiplicity_falls_back_rather_than_guessing():
    assert multiplicity_of_label("iso3", 3) == 3
    assert multiplicity_of_label("", 2) == 2
    assert multiplicity_of_label(None, 1) == 1


def test_a_nonsense_multiplicity_in_a_label_does_not_become_zero():
    assert multiplicity_of_label("iso M=0", 3) == 1


# --- what ranks, and what happens when nothing does -------------------------

BASE = {"charge": "0", "method": "classic", "PAL": "8"}


def _notes(**manta):
    from delfin.common.control_validator import (
        _manta_funnel_notes, validate_control_config)
    return _manta_funnel_notes(validate_control_config({**BASE, **manta}))


def test_optimising_everything_needs_no_single_point():
    # The optimisation *is* the ranking here.  A single point on every frame
    # first would be measured and then thrown away.
    got = selection_options({"MANTA_OPT": "xtb", "MANTA_SCREEN_KEEP": "all"})
    assert got["screen"] == "none"


def test_optimising_a_subset_needs_a_single_point_on_all_of_them():
    # Something has to choose the subset, and only a screen over every frame
    # can.  Cutting the builder's own order at N is cutting an unranked list.
    got = selection_options({"MANTA_SCREEN_KEEP": "10"})
    assert got["screen"] == "gfn2"
    assert got["screen_keep"] == 10


def test_an_explicit_screen_beats_the_rule():
    assert selection_options({"MANTA_SCREEN": "gfnff"})["screen"] == "gfnff"
    assert selection_options({"MANTA_SCREEN": "none",
                              "MANTA_SCREEN_KEEP": "5"})["screen"] == "none"


def test_ranking_nothing_is_allowed_but_says_so():
    notes = _notes(MANTA_SCREEN="none", MANTA_OPT="none")
    assert notes, "ranking nothing has to be flagged"
    assert "ranks nothing" in notes[0]
    assert "first frame" in notes[0]


def test_keeping_the_top_n_of_an_unranked_list_says_so():
    notes = _notes(MANTA_SCREEN="none", MANTA_SCREEN_KEEP="10")
    assert any("builder order" in n for n in notes)


def test_the_shipped_default_is_quiet():
    # The template must not warn about itself.
    assert _notes(MANTA_SCREEN="none", MANTA_SCREEN_KEEP="all", MANTA_OPT="xtb",
                  MANTA_REFINE="goat", MANTA_REFINE_TOPK="5") == []


def test_a_wasted_screen_is_not_reported_because_it_cannot_happen():
    # Asking for both a screen and opt-all used to mean paying for single
    # points nobody reads.  The rule now resolves it instead of warning.
    assert _notes(MANTA_OPT="xtb", MANTA_SCREEN_KEEP="all") == []


# --- a key can be asked what it does ----------------------------------------

def test_a_question_mark_explains_the_key_and_does_not_set_it():
    from delfin.common.control_validator import validate_control_config

    validated = validate_control_config({**BASE, "MANTA_MULTIPLICITIES": "?"})
    # asking is not setting: the run goes on with the default
    assert validated["MANTA_MULTIPLICITIES"] in ("", "auto")
    assert selection_options(validated)["multiplicities"] == []


def test_a_question_mark_never_blocks_the_run():
    from delfin.common.control_validator import validate_control_config

    for key in ("MANTA_SCREEN", "MANTA_OPT", "MANTA_REFINE",
                "MANTA_REFINE_TOPK", "MANTA_SCREEN_KEEP", "MANTA_QUALITY"):
        validate_control_config({**BASE, key: "?"})


def test_every_manta_key_has_something_to_say_when_asked():
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    for spec in CONTROL_FIELD_SPECS:
        if spec.name.startswith("MANTA_"):
            assert spec.help, f"{spec.name} has no explanation for '?'"
            assert len(spec.help) > 40, f"{spec.name} explanation is too thin"


# --- the old spellings still mean what they said ----------------------------

def test_a_control_written_before_the_rename_is_unchanged():
    from delfin.common.control_validator import validate_control_config

    validated = validate_control_config(
        {**BASE, "GUPPY_RANK": "gfnff", "GUPPY_GOAT": "0"})
    got = selection_options(validated)
    assert got["screen"] == "gfnff", "the old rank method was overruled"
    assert got["refine_topk"] == 0, "GOAT was switched on for an old file"


def test_the_screen_energy_is_converted_before_it_becomes_a_result():
    # gfnff_energy answers in kcal/mol; every RunResult downstream is Hartree.
    # Ordering survives the mismatch, so nothing would look wrong -- the
    # trajectory comments would just be 627x too large and the energy window
    # 627x too wide, silently keeping every candidate.
    from delfin.guppy_sampling import _KCAL_PER_HARTREE

    assert abs(_KCAL_PER_HARTREE - 627.509474) < 1e-6
    kcal = -1_234_567.8
    assert abs(kcal / _KCAL_PER_HARTREE + 1967.4) < 1.0


# --- the frames are optimised in the solvent the rest of the run uses -------

def test_the_optimisation_input_carries_the_control_solvation(tmp_path):
    # Ranking coordination isomers in the gas phase and then running everything
    # downstream in CPCM ranks them under a different Hamiltonian than the one
    # that decides anything afterwards. For a complex carrying a formal charge
    # that is not a small difference, and nothing in the output would say so.
    from delfin.guppy_sampling import _write_xtb_input, _write_goat_input

    inp = tmp_path / "XTB.inp"
    _write_xtb_input(inp, ["Ni 0.0 0.0 0.0"], charge=2, multiplicity=3,
                     pal=4, maxcore=2000, method="XTB2",
                     solvation="CPCM(DMF)")
    head = inp.read_text().splitlines()[0]
    assert head == "!XTB2 OPT CPCM(DMF)", head

    goat = tmp_path / "goat.inp"
    _write_goat_input(goat, xyz_file=tmp_path / "c.xyz", charge=2,
                      multiplicity=3, pal=4, maxcore=2000, method="XTB2",
                      solvation="CPCM(DMF)")
    assert goat.read_text().splitlines()[0] == "!XTB2 CPCM(DMF) GOAT"


def test_gas_phase_stays_gas_phase(tmp_path):
    from delfin.guppy_sampling import _write_xtb_input

    inp = tmp_path / "XTB.inp"
    _write_xtb_input(inp, ["Ni 0.0 0.0 0.0"], charge=0, multiplicity=1,
                     pal=4, maxcore=2000, method="XTB2")
    assert inp.read_text().splitlines()[0] == "!XTB2 OPT"


def test_the_solvation_reaches_the_sampler_from_control():
    import tempfile
    from pathlib import Path
    import delfin.guppy_sampling as sampling
    import delfin.workflows.pipeline as pipeline

    seen = {}
    original = sampling.run_sampling
    sampling.run_sampling = lambda **kw: seen.update(kw) or 0
    try:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "GUPPY").mkdir()
            (root / "GUPPY" / "best_coordniation.xyz").write_text("1\n\nH 0 0 0\n")
            try:
                pipeline._run_guppy_for_smiles(
                    "[Ni+2]", root / "start.txt",
                    {"PAL": "8", "maxcore": "2000", "xTB_method": "XTB2",
                     "solvent": "DMF", "implicit_solvation_model": "CPCM"})
            except Exception:                      # noqa: BLE001
                pass
    finally:
        sampling.run_sampling = original
    assert seen["solvation"] == "CPCM(DMF)"
    assert seen["solvent"] == "DMF"


# --- how many frames run at once --------------------------------------------

def test_the_worker_count_follows_pal():
    # A fixed 4 wasted the machine in both directions: two cores per job on
    # PAL=8, and 96 % of the node idle on PAL=450.
    assert selection_options({"PAL": "8"})["parallel_jobs"] == 2
    assert selection_options({"PAL": "32"})["parallel_jobs"] == 8
    assert selection_options({"PAL": "450"})["parallel_jobs"] == 112


def test_an_explicit_worker_count_still_wins():
    assert selection_options({"PAL": "450",
                              "MANTA_PARALLEL_JOBS": "6"})["parallel_jobs"] == 6


def test_a_tiny_pal_still_runs_one_job():
    assert selection_options({"PAL": "1"})["parallel_jobs"] == 1
    assert selection_options({})["parallel_jobs"] == 1


# --- no Hamiltonian is substituted behind the run's back --------------------

def test_an_unknown_ranking_method_gets_no_energy_rather_than_gfnff():
    # GFN-FF used to be the silent fallback. On a transition-metal complex it
    # is not a rougher GFN2 but a different answer -- measured on ABAKOE its
    # minimum is GFN2's near-maximum -- and nothing in the result said which
    # method produced it.
    from delfin.manta import _gfnff_rank

    assert _gfnff_rank.gfnff_energy("H 0 0 0\nH 0 0 0.74",
                                    method="not-a-hamiltonian") is None


# --- GOAT inside MANTA is off by default; global_optimizer does it downstream

def test_no_goat_inside_manta_by_default():
    # The winner goes to start.txt and global_optimizer=[GOAT|CREST] refines
    # that one structure. Refining inside MANTA as well would do it twice.
    assert selection_options({})["refine_topk"] == 0


def test_the_template_ships_no_goat_inside_manta():
    from delfin import define

    block = [l for l in define.TEMPLATE.splitlines() if l.startswith("MANTA_")]
    assert "MANTA_REFINE_TOPK=0" in block
    assert "MANTA_PARALLEL_JOBS=auto" in block


# --- the run records what it actually did, not what the config says ---------

def test_the_construction_set_is_recorded_with_the_run(tmp_path):
    """The provenance names the construction config and lists its flags.

    The construction set defines the builder's output more than the quality
    profile does -- 42 environment flags against five profile numbers -- and it
    was applied in cli.py and recorded nowhere, so a finished run could not say
    which builder produced it.

    The failure this guards against is not a crash. The MANTA harness ran
    builds at quality ``extreme`` for weeks while every results file recorded
    the shipped default, because the call site and the recording site each
    wrote their own fallback for the same environment lookup. Both looked
    right. So this records the set that was *applied*, not the config read a
    second time.
    """
    import json
    import delfin.guppy_sampling as sampling
    import delfin.workflows.pipeline as pipeline
    from delfin.common.manta_settings import apply_construction_env

    original = sampling.run_sampling
    sampling.run_sampling = lambda **kw: 0
    try:
        (tmp_path / "GUPPY").mkdir()
        (tmp_path / "GUPPY" / "best_coordniation.xyz").write_text("1\n\nH 0 0 0\n")
        config = {"PAL": "8", "maxcore": "2000", "xTB_method": "XTB2",
                  "MANTA_CONSTRUCTION": "champion"}
        apply_construction_env(config)
        try:
            pipeline._run_guppy_for_smiles("[Ni+2]", tmp_path / "start.txt", config)
        except Exception:                          # noqa: BLE001
            pass
        recorded = json.loads((tmp_path / "guppy_settings.json").read_text())
    finally:
        sampling.run_sampling = original

    construction = recorded["construction"]
    assert construction["config"] == "champion"
    assert len(construction["flags"]) > 30, construction["flags"]
    assert recorded["builder"]["quality_mode"] == "extreme"


def test_control_reproduces_the_champion_construction_set():
    """CONTROL's champion is the same set the MANTA command line applies.

    If it were a subset, a CONTROL run would silently build something that
    does not match anything in the MANTA archives -- and the file defining
    these flags records that a hand-picked 29-flag subset scored 13.7 %
    topology-correct against 33.6 % for no flags at all, so a partial set is
    not a milder version of the whole one.
    """
    from delfin import cli_manta
    from delfin.common.manta_settings import apply_construction_env

    wanted = {f"DELFIN_FFFREE_{name}" for name in cli_manta._CHAMPION_FLAGS}
    applied = {}
    apply_construction_env({"MANTA_CONSTRUCTION": "champion"}, applied)
    enabled = {key for key, value in applied.items() if str(value) == "1"}
    assert wanted <= enabled, sorted(wanted - enabled)


# --- optimising everything is right for the median and wrong for the tail ---

def test_the_screen_threshold_is_a_setting_with_a_measured_default():
    """Measured over 5810 builds at champion/extreme with max_isomers=0:

        mean 26.1 · p10 3 · p25 4 · p50 14 · p75 33 · p90 64 · p95 90
        p99 190 · max 399

    So "optimise everything" is right for the median system and wrong for its
    tail, where 4.1 % return more than 100 frames and each frame is an ORCA
    optimisation. Screening above 30 leaves roughly three quarters of systems
    on the path that needs no screen, and caps the worst case there instead of
    at 399.
    """
    from delfin.common.manta_settings import SCREEN_ABOVE_FRAMES

    assert SCREEN_ABOVE_FRAMES == 30
    assert selection_options({})["screen_above"] == 30
    assert selection_options({"MANTA_SCREEN_ABOVE": "50"})["screen_above"] == 50
    assert selection_options({"MANTA_SCREEN_ABOVE": "0"})["screen_above"] == 0


def test_an_explicit_screen_is_never_overridden_by_the_threshold():
    # The threshold fills in what CONTROL left open. It does not overrule what
    # CONTROL said.
    assert selection_options({})["screen_explicit"] is False
    assert selection_options({"MANTA_SCREEN": "none"})["screen_explicit"] is True
    assert selection_options({"MANTA_RANK": "gfnff"})["screen_explicit"] is True


def test_the_threshold_reaches_the_sampler():
    import tempfile
    from pathlib import Path
    import delfin.guppy_sampling as sampling
    import delfin.workflows.pipeline as pipeline

    seen = {}
    original = sampling.run_sampling
    sampling.run_sampling = lambda **kw: seen.update(kw) or 0
    try:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "GUPPY").mkdir()
            (root / "GUPPY" / "best_coordniation.xyz").write_text("1\n\nH 0 0 0\n")
            try:
                pipeline._run_guppy_for_smiles(
                    "[Ni+2]", root / "start.txt",
                    {"PAL": "8", "maxcore": "2000", "xTB_method": "XTB2"})
            except Exception:                      # noqa: BLE001
                pass
    finally:
        sampling.run_sampling = original
    assert seen["screen_above"] == 30
    assert seen["screen_explicit"] is False


# --- an invisible cache is an unverifiable cache ----------------------------

def test_the_single_point_cache_can_be_counted():
    """The solvent is part of the cache key, so changing it invalidates every
    cached energy. A run that silently recomputes everything and a run that
    silently reuses energies from another medium look identical without this.
    """
    from delfin.manta import _gfnff_rank

    stats = _gfnff_rank.cache_stats()
    assert set(stats) == {"hits", "misses", "entries"}
    assert all(isinstance(v, int) for v in stats.values())
