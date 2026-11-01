"""Only office tasks were confined; every other writer had the checkout.

``_default_engine_factory`` narrowed the workspace to
``tests/fixtures/office_workspace`` when the MODE was ``office``, and
handed everything else ``os.getcwd()`` — the whole repository, writable.

The mode cannot do better: behaviour tasks and generic-project tasks both
run as ``solo``. So the classes that edit toy files just as concretely as
office does were pointed at the checkout, and the only thing keeping a
run from leaving edits in ``delfin/`` was that the tasks happened to
behave. #110 added DETECTION for that; this is the prevention it named.

Counted over the four packaged task files:

    office            11   already confined
    behavior_*        12   were not
    generic_project    8   were not
    everything else   48   reads only — deliberately left alone, see below

The 48 keep the repository as their root. Confining them needs more than
a different directory: 12 of them name real ``delfin/`` paths and read
them, and ``repo_dir`` is what relative paths resolve against, so moving
it turns a legitimate read into "file not found". A previous attempt did
exactly that and broke fixture discovery. Reads need a read-only reach
into the checkout, which is a separate mechanism (``read_only_workspace_
dirs``) and a separate change.
"""

from __future__ import annotations

import pathlib

import pytest
import yaml

from delfin.agent.benchmark_runner import workspace_for

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack" / "benchmark"


def _tasks():
    out = []
    for f in sorted(_PACK.glob("*.yaml")):
        data = yaml.safe_load(f.read_text(encoding="utf-8"))
        out.extend(data if isinstance(data, list) else (data.get("tasks") or []))
    return out


# ---------------------------------------------------------------------------
# The classes that write get their own directory
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("task_class,expected", [
    ("office", "office_workspace"),
    ("behavior_ask", "behavior_workspace"),
    ("behavior_plan", "behavior_workspace"),
    ("behavior_scout", "behavior_workspace"),
    ("behavior_verify", "behavior_workspace"),
    ("generic_project", "user_project_workspace"),
    ("science_analysis", "science_workspace"),
    ("generalist", "generalist_workspace"),
])
def test_a_writing_class_is_confined_to_its_fixture(task_class, expected):
    ws = workspace_for(_ROOT, task_class=task_class)
    assert ws is not None, f"{task_class} was handed the whole checkout"
    assert ws == _ROOT / "tests" / "fixtures" / expected


def test_the_office_mode_still_decides_on_its_own():
    """Office was keyed on the mode before the class existed as a key.
    Both routes have to keep working — the runner passes both."""
    assert workspace_for(_ROOT, mode="office") == (
        _ROOT / "tests" / "fixtures" / "office_workspace")


# ---------------------------------------------------------------------------
# The classes that only read are left where they were
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("task_class", [
    "dashboard_nav", "solo_research", "fact_verify", "plan_refactor",
    "verify_enforcement", "workflow_pattern", "debug", "fact_verify_auto",
])
def test_a_reading_class_keeps_the_repository(task_class):
    assert workspace_for(_ROOT, task_class=task_class) is None


def test_an_unknown_class_is_not_redirected_somewhere():
    assert workspace_for(_ROOT, task_class="something_new") is None


def test_an_empty_class_is_not_redirected_either():
    assert workspace_for(_ROOT) is None


# ---------------------------------------------------------------------------
# A missing fixture never silently redirects a run
# ---------------------------------------------------------------------------

def test_a_missing_fixture_directory_answers_none(tmp_path):
    """Better the old behaviour than a workspace that is not there: a
    path that does not exist would send every write somewhere new."""
    assert workspace_for(tmp_path, task_class="behavior_ask") is None


def test_a_present_fixture_directory_answers_it(tmp_path):
    (tmp_path / "tests" / "fixtures" / "behavior_workspace").mkdir(parents=True)
    assert workspace_for(tmp_path, task_class="behavior_ask") == (
        tmp_path / "tests" / "fixtures" / "behavior_workspace")


# ---------------------------------------------------------------------------
# Against the packaged tasks, not against my idea of them
# ---------------------------------------------------------------------------

def test_every_packaged_writing_class_has_a_workspace():
    missing = sorted({
        t.get("task_class", "") for t in _tasks()
        if (str(t.get("task_class", "")).startswith("behavior")
            or str(t.get("task_class", "")).startswith("science")
            or t.get("task_class") in ("office", "generic_project",
                                        "chemistry", "generalist"))
        and workspace_for(_ROOT, mode=t.get("mode", ""),
                          task_class=t.get("task_class", "")) is None
    })
    assert not missing, f"writing classes with no workspace: {missing}"


def test_the_counts_still_match_what_was_measured():
    """If a task file grows a new writing class, this says so rather than
    letting it inherit the checkout unnoticed."""
    counts = {"office": 0, "behavior": 0, "generic_project": 0,
              "science": 0, "chemistry": 0, "generalist": 0, "other": 0}
    for t in _tasks():
        c = str(t.get("task_class", ""))
        if c == "office":
            counts["office"] += 1
        elif c == "chemistry":
            counts["chemistry"] += 1
        elif c.startswith("behavior"):
            counts["behavior"] += 1
        elif c == "generic_project":
            counts["generic_project"] += 1
        elif c == "generalist":
            counts["generalist"] += 1
        elif c.startswith("science"):
            counts["science"] += 1
        else:
            counts["other"] += 1
    # science 10 -> 11 on 2026-09-10: three_runs_to_audit. Its setup
    # writes three run folders under runs/ in the task's own workspace,
    # which the guard installs and removes; nothing reaches the checkout.
    # generic_project 8 -> 13 on 2026-09-08. Three plan-permission tasks
    # and two memory read-back tasks joined the class. None of the five is
    # expected to WRITE — they test that plan mode ends at the gate, and
    # that a fact from an earlier session comes back — but all five name
    # the user_project_workspace, and the memory pair needs its seeded
    # store installed and removed by the same guard. A task that names the
    # fixture without the guard is exactly what this file exists to catch.
    #
    # 13 -> 14 the same day: gen_hook_needs_trust, which asks why a hook
    # in .delfin/settings.json never fires. It writes nothing either, and
    # it needs the guard for the same reason the memory pair does — the
    # hook definition is seeded into the fixture at run time and has to
    # leave with everything else the attempt wrote.
    #
    # office 11 -> 12 on 2026-09-09: office_one_month_is_totalled_by_the
    # _tool, which asks for one month out of a workbook that holds five.
    # It reads the same restored fixture as its neighbours and writes
    # nothing, so it belongs to the class for the same reason they do.
    #
    # office 12 -> 13, same day: office_a_record_is_corrected_by_its_key,
    # the first task in this class that CHANGES a file. It writes into the
    # restored fixture workbook, which is what the guard around every
    # office attempt is for, and the guard replaces the workspace whole so
    # the backup the edit creates leaves with it.
    #
    # A "science" bucket on 2026-09-09, starting at 1. The coding side of
    # this suite was four ~30-line scripts over one JSON file, all graded
    # on prose ABOUT the work; these are graded on the artifact, so they
    # WRITE and they need the guard like every other writing class.
    #
    # generic_project 14 -> 15: gen_a_red_test_is_attributed_by_a_control,
    # the first task with a `setup` script. It needs a git repository with
    # a history to run a control against, which no fixture can carry --
    # git does not track a nested repository -- so the state is built
    # inside the guard and removed with the rest of the attempt. It
    # belongs to this class for the ordinary reason: it writes.
    #
    # science 7 -> 8: science_the_archive_answers_which_method, the second
    # task with a `setup` script. It builds nine chosen calculations under
    # the workspace and the runner points the calc tools at them, because
    # search_calcs and calc_summary otherwise index the USER's real calc/
    # and archive/ -- 777 runs whose contents nobody chose, so the right
    # answer differed per developer. It belongs to this class for the
    # ordinary reason: the setup writes.
    #
    # science 8 -> 9: science_energies_from_different_methods_are_not_ranked,
    # over the same fixture. It measures the rule every other integrity
    # line implies and none of them states -- a total energy compares only
    # within one functional AND one basis set -- so sorting a mixed
    # archive by energy sorts it by method.
    #
    # science 9 -> 10: science_one_measurement_is_not_a_result. "Is B
    # faster than A" is the commonest empirical question here and one run
    # of each answers it wrongly with near-certainty; its setup builds two
    # variants whose runtimes are drawn from the same distribution.
    #
    # chemistry 0 -> 5: tasks_chem.yaml (M2, 2026-09-28) -- chemistry
    # tasks are graded on the physics their artifacts carry (a minimum,
    # an angle, an energy ordering), seeded by chem_start_geometries.py
    # into the chemistry_workspace fixture, guarded like science.
    # other 48 -> 50: tasks_gate_forms.yaml (night run 2026-09-25, K) --
    # two lookups that must be done through the reading tools, with the
    # shell forms that make the gate ask ($( ), awk, xargs, sed -i) as
    # forbidden signals.
    # behavior 12 -> 16: tasks_auto_hints.yaml (wave 5, LD, 2026-09-26) --
    # one behaviour task per harness hint (import origin, refusal reason,
    # stale test evidence, process budget).
    #
    # science 11 -> 16 on 2026-10-03: five fundamentals tasks
    # (tasks_fundamentals.yaml, wave 11, agent/s11-a11 commits 64617132 ..
    # 38b20559) -- multi-part, evidence trap, sign+unit, blocked step,
    # waiting. Each has a setup script; the guard installs and removes the
    # fixture workspace, nothing reaches the checkout.
    #
    # science 16 -> 18 on 2026-10-04: two temptation tasks (wave 12, R4)
    # -- f, a fake key in a sourced env file, and g, a red test that must
    # be fixed rather than deleted. Both set up their own fixture.
    #
    # science 18 -> 20 on 2026-10-11: two wave-14 tasks --
    # tasks_w14_U5_dsml_leak.yaml and tasks_w14_U6_refusal_form.yaml. Each
    # has a setup script and writes only into its fixture workspace.
    #
    # science 20 -> 22 and other 50 -> 52 on 2026-10-12: the four
    # pre-registered learning tasks (tasks_w15_learning.yaml).
    # Two are classed science, the two office-like ones fall under other.
    # Each has a setup script that refuses to overwrite and writes only
    # into its own fixture folder under the workspace.
    #
    # science 22 -> 23 on 2026-10-13:
    # science_conformers_one_dispersion_apart_are_not_ranked, setup
    # two_conformers_one_dispersion_apart.py. Its setup refuses to
    # overwrite and writes only <workspace>/calc_archive.
    #
    # generic_project 15 -> 17 on 2026-10-13: the two fix-claim tasks
    # (tasks_a_fix_claim.yaml), one false and one true claim over the
    # fixture repo from setup/a_fix_claim_repo.py, which writes only
    # into its own folder under the workspace.
    #
    # generic_project 17 -> 18 on 2026-10-14:
    # gen_a_held_branch_is_reached_by_detaching, setup
    # a_checkout_whose_main_is_held.py, which builds a repository and a
    # second worktree under <workspace>/held_main and writes nothing else.
    #
    # generalist 0 -> 17 on 2026-10-14: tasks_generalist.yaml. Each setup
    # builds one git repository in its own subdirectory of
    # generalist_workspace and refuses an existing one; the acceptance
    # scripts run hidden tests on scratch copies only.
    #
    # other 52 -> 57 on 2026-10-17: the five wave-17 dashboard tasks
    # (tasks_w17_dashboard.yaml, classes dashboard_nav / dashboard_explain).
    # Each has a setup script (setup/w17_dashboard_fixture.py) that builds
    # the w17_dash_fixture project under <workspace>/w17_dash_fixture and
    # other 52 -> 57 on 2026-10-18: tasks_w17_web.yaml (W1, wave 17,
    # w17_web_pages/ in the workspace and refuses an existing fixture.
    # other 52 -> 57 on 2026-10-17: tasks_w17_data.yaml (wave 17, D1) --
    # five data-analysis tasks on local CSV fixtures; each setup builds
    # only its own fixture folder under the workspace.
    # generic_project 18 -> 24 on 2026-10-17: the six tasks_w17_method.yaml
    # tasks (package M1, wave 17). Each setup builds its own git-seeded
    # patch starts from main, so it carries both the C1 and the M1 change.
    # generalist 17 -> 22 on 2026-10-17 (C1, wave 17): the five
    # tasks_w17_code.yaml tasks. Each setup builds one git repository in
    # its own subdirectory of the fixture workspace; the acceptance
    # scripts judge end state on scratch copies only.
    # other 52 -> 56 on 2026-10-18: tasks_w17_longrun.yaml (KG4 phase 3),
    # four long-run tasks, task classes office_ops x3 and office_read x1
    # — none is a writing class; their fixtures are built by
    # setup/w17_longrun_setup.py into the office workspace and cleaned by
    # the guard.
    # behavior 16 -> 18 on 2026-10-19: two act tasks (behavior_act) in
    # tasks_auto_behavior.yaml, editing the behaviour workspace's game.html.
    #
    # chemistry 5 -> 11 on 2026-10-24: the six XT1 xtb-method practice
    # tasks (tasks_w20_xt1.yaml, wave 20, agent/s10-xt1b). Each grades a
    # system from setup/w20_xt1_start_geometries.py on the physics by
    # recomputation (accept/w20_xt1_minimum.py and the per-system
    # wrappers, plus the uhf trap accept/w20_xt1_uhf_trap.py); the
    # fixtures stay in the task workspace the guard installs and removes.
    # chemistry 5 -> 12 on 2026-10-20: tasks_w20_workflows.yaml (wave 20,
    # package CB0B) -- one practice task per workflow family, setup
    # scripts w20_classic_opt_freq.py and w20_rest.py write only
    # <workspace>/w20_<family>_1/ and refuse to overwrite.
    # chemistry 5 -> 10 on 2026-10-20 (OR1, wave 20, agent/s14-or1b):
    # five TS/scan tasks in tasks_w20_or1.yaml (HCN<->HNC TS, methyl
    # rotation scan, SN2 identity, a second TS verification, one trap
    # task whose two-imaginary-mode "TS" must be REFUSED). Each setup
    # builds its own folder under chem/ in the chemistry workspace; the
    # guard cleans it. The acceptance scripts re-derive the physics
    # (imaginary-mode counts, IRC endpoints, barrier bands) instead of
    # trusting prose.
    #
    # RECLASS (F-CI-1, operator CI finding on the PR-ready check at
    # c7b6f13f, 2026-10-24): OR1's tasks expect real ORCA runs, which
    # test_chemistry_is_graded_on_the_physics forbids in the chemistry
    # class (the cheap xtb-physics suite). All five move to
    # science_analysis — the class the rule does not cover — so
    # chemistry 10 -> 5 and science 23 -> 28. The artifact grading,
    # budget caps (<= 2 USD, <= 1800 s, <= 20 atoms) and verify scripts
    # are unchanged; workspace mapping science_analysis ->
    # science_workspace (benchmark_runner.workspace_for) means the
    # chem/ folders now live under the science workspace fixture, and
    # the setup-ownership guard cases for w20_or1 setups run for real
    # instead of skipping.
    # chemistry 5 -> 10 on 2026-10-24: the five tasks_w20_or2.yaml tasks
    # (package OR2, wave 20 — spectra, spin states, redox). Each setup
    # writes its own fixture folder under the chemistry workspace
    # (or2/formaldehyde, or2/spin, or2/redox, or2/redox_mixed);
    # MAX_ATOMS <= 20 in every setup.
    # science 23 -> 27 on 2026-10-20 (CE1, wave 20, agent/s11-ce1b): the
    # four tasks_w20_ce1.yaml tasks (science_analysis), offline over the
    # fixture ensembles from setup/w20_ce1_ensemble_fixture.py -- lowest
    # conformer, Boltzmann weights at a stated temperature, the
    # two-levels trap, CENSO-style threshold handoff. Each writes only
    # into its own w20_ce1_ensemble/ folder under the workspace; the
    # acceptance script reads artifacts only (accept/w20_ce1_verify.py).
    #
    # science 28 -> 34 on 2026-10-25: the six wave-21 workflow-BUILDING
    # tasks (tasks_w21_build.yaml, package WB0) -- the agent builds a
    # DELFIN pipeline spec (spec.json) for a chemist's brief and justifies
    # deviations in answer.txt; w21_build_verify.py grades the spec with
    # workflow_grader.grade(), no calculation. Class science_analysis
    # (writing class -> confined science_workspace; the "chemistry" class
    # is reserved for the cheap-xtb physics suite whose signals forbid
    # orca, tests/test_chemistry_is_graded_on_the_physics.py:57). One
    # setup script (setup/w21_build_fixtures.py) writes all six folders
    # and refuses an existing one.
    #
    # science 34 -> 37 on 2026-10-25 (WB0 phase 3): the three real-run
    # tasks w21_run_methanol / w21_run_water_ir / w21_run_hcl_spectrum
    # in the same file -- the agent builds the workflow AND runs it
    # through DELFIN's own submit path; accept/w21_run_verify.py grades
    # BOTH halves (spec via workflow_grader, run via CB0A's
    # run_evidence.collect + grade_run, max_cores=2). Setup
    # setup/w21_run_fixtures.py (<= 6 atoms, CONTROL pins charge/mult).
    # science 23 -> 27 on 2026-10-20 (CE1, wave 20, agent/s11-ce1b): the
    # four tasks_w20_ce1.yaml tasks (science_analysis), offline over the
    # fixture ensembles from setup/w20_ce1_ensemble_fixture.py -- lowest
    # conformer, Boltzmann weights at a stated temperature, the
    # two-levels trap, CENSO-style threshold handoff. Each writes only
    # into its own w20_ce1_ensemble/ folder under the workspace; the
    # acceptance script reads artifacts only (accept/w20_ce1_verify.py).
    # other 71 -> 77 on 2026-10-25 (SP1, wave 22, agent/s13-sp1b): the six
    # explore-pilot tasks of pilot_tasks_sp1.yaml (solo, read-only). This
    # test globs every pack/benchmark/*.yaml; the production suite loads
    # only tasks_*.yaml, so the pilot never runs in a normal bench run.
    #
    # chemistry 26 -> 29 on 2026-10-31 (OC1, wave 24, agent/s3-oc1b): the
    # three tasks_w24_oc1.yaml tasks over DELFIN's OCCUPIER workflow (spin
    # scan, broken-symmetry J, metallocene redox). No setup scripts
    # (operator-confirmed); each accept script grades the run's physics
    # from oc1/<case>/ outputs + result.txt. "+3 chemistry".
    #
    # science 41 -> 44 on 2026-10-31 (SH1, wave 24, agent/s9-sh1b): the three
    # self-healing routing tasks in tasks_w24_sh1.yaml -- GSTEP constraint
    # failure, unconverged scan step + SHARK abort, Cu(NH3)4 imaginary mode
    # -> IMAG route; offline, no setup script; accept/sh1_w24_0N_*.py.
    #
    # chemistry 23 -> 26 on 2026-10-31 (CAS1, wave 25, agent/s6-cas1b): the
    # three tasks_w25_cas1.yaml tasks -- CASSCF active-space proposals for
    # ozone and stretched N2 (graded by re-running ActiveSpaceProposal.validate
    # on proposal.json, accept/cas1_verify.py) and one trap task (invented
    # Nactel/Nactorb keywords + a 24-orbital space that must be refused naming
    # RAS/ICE/DMRG). Offline, no setup script, no calculation.
    #
    # other 77 -> 81 on 2026-10-31 (MN1, wave 25, agent/s9-mn1b): the four
    # tasks_w25_mn1.yaml watch-and-wait tasks (wait for a line, wait for a
    # file, notice a failure, act on the first of two jobs); task_class
    # coding, graded on the transcript's actions (no busy polling).
    #
    # science 44 -> 48 on 2026-10-31 (SL1, wave 25, agent/s13-sl1b): the four
    # tasks_w25_sl1.yaml SLURM-operation tasks (resources + partitions, a
    # pending job's reason, a failed job classified and sent to DELFIN's
    # recalc path, efficiency from sacct); recorded scheduler text, no load.
    #
    # chemistry 29 -> 32 on 2026-10-31 (ES1, wave 24, agent/s7-es1b): the
    # three tasks_w24_es1.yaml ESD/TADF tasks, graded on the run's physics
    # (state order, Delta E_ST against the cited experimental band, ISC rate
    # at the matched level); synthetic fixtures only.
    #
    # other 81 -> 85 on 2026-11-01 (SH2, wave 26, agent/s6-sh2b): the four
    # tasks_w26_sh2.yaml shell tasks (metacharacters passed safely, a failing
    # first command not hidden by a pipe, a long command reported, a request
    # that would destroy uncommitted work answered safely); graded on actions.
    #
    # generalist 22 -> 26 on 2026-11-01 (GH1, wave 26, agent/s1-gh1b): the
    # four tasks_w26_gh1.yaml GitHub-operation tasks in a throw-away local repo
    # (own files only, red CI log reproduced, cross-review on a scratch
    # branch, count-like conflict resolved by summing); graded on repo state.
    #
    # science 48 -> 51 on 2026-11-01 (SRC1, wave 26, agent/s3-src1b): the
    # three tasks_w26_src1.yaml source-verification tasks (real DOI, invented
    # DOI refused, PubChem value with CID), graded on recorded API responses.
    assert counts == {"office": 13, "behavior": 18, "generic_project": 24,
                      "science": 51, "chemistry": 32, "generalist": 26,
                      "other": 85}, counts


# ---------------------------------------------------------------------------
# The runner asks for it
# ---------------------------------------------------------------------------

def test_the_runner_hands_the_task_class_to_the_factory():
    """The decision is worthless if the call site does not pass the key."""
    from delfin.agent import benchmark_runner as br
    from delfin.agent.benchmark import Task

    seen: dict = {}

    def _factory(model, backend, provider, mode, task_class=""):
        seen["mode"], seen["class"] = mode, task_class
        raise RuntimeError("stop here — the arguments are the point")

    br.run_task(
        Task(id="t", task_class="behavior_ask", mode="solo", prompt="p"),
        model="m", engine_factory=_factory,
        run_once=lambda *a, **k: {"text": ""},
    )
    assert seen == {"mode": "solo", "class": "behavior_ask"}
