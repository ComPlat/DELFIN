"""Self-healing ORCA OPI block: classify, escalate, never change chemistry.

The block (``delfin.tools.adapters.orca_selfheal_opi``) classifies an ORCA
failure from the ``.out`` text (SCF not converged, geometry not converged,
frequency failure, disk/memory) and proposes an escalating retry plan.  Every
plan step is a *catalog-known* keyword / ``%scf`` / ``%geom`` block option (the
manual sections are cited in the module docstring) and — by construction and by
an explicit guard — never alters method, basis, charge, multiplicity or
electronic state.  Self-healing is opt-in on each OPI adapter
(``self_heal``, off by default).

These are unit tests: the classifiers run on hand-written synthetic ORCA-format
text, never real ORCA output (fixtures are lines of the shape ORCA prints).
The probe test hands OPI-built inputs to ORCA's own parser when the binary is
present and skips otherwise.  Runs under OPI (py3.12 job); skipped on the
no-OPI interpreter.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.tools.adapters import orca_selfheal_opi as heal  # noqa: E402
from delfin.orca_knowledge import known_block_option, known_keyword  # noqa: E402

# ---------------------------------------------------------------------------
# Fixtures: synthetic ORCA-format failure tails (hand-written, not ORCA output)
# ---------------------------------------------------------------------------

_SCF_NOT_CONVERGED = """\
SCF ITERATIONS
------------------
ORCA SCF NOT CONVERGED AFTER 128 ITERATIONS

****ORCA TERMINATED NORMALLY****
"""

_GEOM_NOT_CONVERGED = """\
GEOMETRY OPTIMIZATION CYCLE 150
        **CONVERGENCE FAILURE**
        MAX Gradient:   0.00004200
        RMS Gradient:   0.00002100
   The optimization has not yet converged
   *** MAXIMUM NUMBER OF OPTIMIZATION CYCLES REACHED ***
"""

# The exact sentence ORCA 6.1.1 prints when a distorted geometry exceeds the
# optimization budget yet ORCA still exits normally (see REAL-RUN FAIL #597).
# Isolated on purpose: it must be classified by the gave-up sentence alone, not
# by a co-present "CONVERGENCE FAILURE" line.
_GEOM_GAVE_UP_REAL_611 = """\
GEOMETRY OPTIMIZATION CYCLE 2
The optimization did not converge but reached the maximum number of optimization cycles.
****ORCA TERMINATED NORMALLY****
"""

_FREQ_FAILURE = """\
FREQUENCY CALCULATION
  ... ORCA frequency calculation aborted in module ORCA_FREQ
  ****ORCA TERMINATED WITH ERROR****
CARTESIAN FORCE CONSTANT MATRIX EVALUATION...
"""

_MEMORY_ERROR = """\
Error  (ORCA_SCF): Not enough memory available!
      ====>        Please increase MaxCore to more than:  800.0 MB
CARTESIAN FORCE CONSTANT MATRIX EVALUATION...
"""

_CLEAN = """\
SCF CONVERGED AFTER 11 ITERATIONS
****ORCA TERMINATED NORMALLY****
"""


@pytest.fixture
def elements() -> list[str]:
    return ["O", "H", "H"]


@pytest.fixture
def coords() -> list[list[float]]:
    return [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.95],
        [0.9, 0.0, -0.3],
    ]


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


class TestClassifyFailure:
    def test_scf_not_converged(self):
        assert heal.classify_failure(_SCF_NOT_CONVERGED) is heal.SelfHealFailure.SCF_NOT_CONVERGED

    def test_geometry_not_converged(self):
        assert heal.classify_failure(_GEOM_NOT_CONVERGED) is heal.SelfHealFailure.GEOMETRY_NOT_CONVERGED

    def test_geometry_gave_up_real_611_sentence(self):
        # ORCA 6.1.1 prints "The optimization did not converge but reached the
        # maximum number of optimization cycles." and still terminates normally;
        # the recovery module detects this via _OPT_GAVE_UP and so must we.
        assert heal.classify_failure(_GEOM_GAVE_UP_REAL_611) \
            is heal.SelfHealFailure.GEOMETRY_NOT_CONVERGED

    def test_memory_error_takes_precedence_over_scf(self):
        # SCF text present, but the aborting cause is memory: Memory > SCF.
        assert heal.classify_failure(_MEMORY_ERROR + "\n" + _SCF_NOT_CONVERGED) \
            is heal.SelfHealFailure.MEMORY_ERROR

    def test_frequency_failure(self):
        assert heal.classify_failure(_FREQ_FAILURE) is heal.SelfHealFailure.FREQUENCY_FAILURE

    def test_clean_output_is_none(self):
        assert heal.classify_failure(_CLEAN) is heal.SelfHealFailure.NONE

    def test_empty_output_is_none(self):
        assert heal.classify_failure("") is heal.SelfHealFailure.NONE


# ---------------------------------------------------------------------------
# Escalation ladder: every emitted keyword / block option is catalog-known and
# the chemistry surface is never touched.
# ---------------------------------------------------------------------------

_CHEM_KEYS = {"functional", "basis", "ri", "aux_basis", "dispersion",
              "charge", "mult", "method_keywords"}


class TestPlanEscalation:
    def test_scf_ladder_produces_distinct_steps(self, elements, coords):
        plan1 = heal.plan_escalation(heal.SelfHealFailure.SCF_NOT_CONVERGED, 1, job_kind="SP")
        plan2 = heal.plan_escalation(heal.SelfHealFailure.SCF_NOT_CONVERGED, 2, job_kind="SP")
        plan3 = heal.plan_escalation(heal.SelfHealFailure.SCF_NOT_CONVERGED, 3, job_kind="SP")
        assert plan1 != plan2 and plan2 != plan3  # each rung is distinct
        # Distinct progress: strictly weaker at step 1 than at step 3.
        assert plan1.step == 1 and plan3.step == 3
        assert plan1 != plan3
        assert plan3.scf_options.get("MaxIter", plan3.scf_options.get("maxiter", 0)) >= \
            plan1.scf_options.get("MaxIter", plan1.scf_options.get("maxiter", 0))

    def test_scf_ladder_keywords_are_catalog_known(self):
        for step in (1, 2, 3):
            plan = heal.plan_escalation(heal.SelfHealFailure.SCF_NOT_CONVERGED, step, job_kind="SP")
            for kw in plan.keywords:
                assert known_keyword(kw), f"{kw!r} (step {step}) not in catalog"
            for key, _val in plan.scf_options.items():
                assert known_block_option("scf", key), \
                    f"%scf {key!r} (step {step}) not in catalog"

    def test_geometry_ladder_uses_known_geom_options(self):
        for step in (1, 2, 3):
            plan = heal.plan_escalation(
                heal.SelfHealFailure.GEOMETRY_NOT_CONVERGED, step, job_kind="OPT")
            for key, _val in plan.geom_options.items():
                assert known_block_option("geom", key), \
                    f"%geom {key!r} (step {step}) not in catalog"

    def test_scf_frequency_rungs_pin_identical_thresh(self):
        # A healed run must reproduce the base energy. Different convergence
        # levels set different %scf Thresh (manual table 2.9), so every ladder
        # rung pins the SAME tight DIIS threshold to converge to one energy.
        seen = set()
        for failure in (heal.SelfHealFailure.SCF_NOT_CONVERGED,
                        heal.SelfHealFailure.FREQUENCY_FAILURE):
            for step in (1, 2, 3):
                plan = heal.plan_escalation(failure, step, job_kind="SP")
                assert known_block_option("scf", "thresh")
                assert "thresh" in plan.scf_options, \
                    f"rung {failure.value}/{step} missing identical thresh"
                assert isinstance(plan.scf_options["thresh"], float)
                seen.add(plan.scf_options["thresh"])
        assert len(seen) == 1, f"rungs must all pin the same Thresh, got {seen}"

    def test_sp_geometry_escalation_is_refused(self):
        # A single-point has no geometry optimization: escalation is refused.
        with pytest.raises(ValueError, match="geometry|single.point|SP"):
            heal.plan_escalation(heal.SelfHealFailure.GEOMETRY_NOT_CONVERGED, 1, job_kind="SP")

    def test_plan_never_touches_chemistry(self):
        plan = heal.plan_escalation(heal.SelfHealFailure.SCF_NOT_CONVERGED, 2, job_kind="OPT_FREQ")
        for surface in (_CHEM_KEYS,):
            assert not (surface & set(map(str, plan.scf_options)))
            assert not (surface & set(map(str, plan.geom_options)))
            assert not (surface & set(plan.keywords))
        # Never rewrites functional/basis/charge/mult: no such key appears.
        planned_keys = {str(k).lower() for k in plan.keywords}
        planned_keys |= {str(k).lower() for k in plan.scf_options}
        planned_keys |= {str(k).lower() for k in plan.geom_options}
        assert not planned_keys & {"functional", "basis", "ri", "aux_basis",
                                   "dispersion", "charge", "mult", "method_keywords"}

    def test_out_of_range_step_raises(self):
        with pytest.raises(ValueError, match="step"):
            heal.plan_escalation(heal.SelfHealFailure.SCF_NOT_CONVERGED, 0)

    def test_beyond_max_steps_raises(self):
        with pytest.raises(ValueError, match="step"):
            heal.plan_escalation(heal.SelfHealFailure.SCF_NOT_CONVERGED, 99)


# ---------------------------------------------------------------------------
# compose_overrides — merges onto user options without mutating them
# ---------------------------------------------------------------------------


class TestComposeOverrides:
    def test_returns_new_dict_and_keeps_user_values(self):
        user = {"MaxIter": 200}
        added = {"MaxIter": 500, "DampFac": 0.9}
        merged = heal.compose_overrides(user, added)
        assert merged["MaxIter"] == 500            # escalation wins
        assert merged["DampFac"] == 0.9
        assert user == {"MaxIter": 200}            # unchanged
        assert merged is not user

    def test_none_user_options(self):
        assert heal.compose_overrides(None, {"MaxIter": 500}) == {"MaxIter": 500}
        assert heal.compose_overrides({"MaxIter": 500}, None) == {"MaxIter": 500}

    def test_does_not_merge_chemistry_keys(self):
        added = {"Functional": "B3LYP"}  # a rogue key that belongs in the chemistry surface
        merged = heal.compose_overrides({"MaxIter": 200}, added)
        assert "Functional" not in merged
        assert merged == {"MaxIter": 200}


# ---------------------------------------------------------------------------
# guard_unchanged_chemistry — refuse phantom chemistry changes
# ---------------------------------------------------------------------------


class TestGuardUnchangedChemistry:
    def test_identical_surfaces_pass(self):
        base = {"functional": "B3LYP", "basis": "def2-SVP", "charge": 0, "mult": 1}
        esc = {"scf_options": {"MaxIter": 500}, "keywords": ["SlowConv"]}
        heal.guard_unchanged_chemistry(base=base, **esc)  # no raise

    def test_changed_method_raises(self):
        base = {"functional": "B3LYP", "basis": "def2-SVP", "charge": 0, "mult": 1}
        esc = {"functional": "PBE0", "scf_options": {"MaxIter": 500}}
        with pytest.raises(ValueError, match="functional"):
            heal.guard_unchanged_chemistry(base=base, **esc)

    def test_changed_charge_raises(self):
        base = {"functional": "B3LYP", "basis": "def2-SVP", "charge": 0, "mult": 1}
        esc = {"charge": -1}
        with pytest.raises(ValueError, match="charge"):
            heal.guard_unchanged_chemistry(base=base, **esc)

    def test_changed_mult_raises(self):
        base = {"functional": "B3LYP", "basis": "def2-SVP", "charge": 0, "mult": 1}
        esc = {"mult": 3}
        with pytest.raises(ValueError, match="mult"):
            heal.guard_unchanged_chemistry(base=base, **esc)

    def test_defaults_match_when_absent(self):
        # Absent keys use the documented defaults, so a "different" key that is
        # actually the default is fine.
        base = {"functional": "B3LYP", "basis": "def2-SVP", "charge": 0, "mult": 1}
        heal.guard_unchanged_chemistry(base=base, scf_options={"MaxIter": 500})


# ---------------------------------------------------------------------------
# Registration — the OPI adapters expose an opt-in, off-by-default knob
# ---------------------------------------------------------------------------


class TestRegistration:
    def test_sp_adapter_exposes_self_heal_off(self):
        from delfin.tools.adapters.orca_sp_opi import OrcaSpOpiAdapter

        names = {p.name for p in OrcaSpOpiAdapter.params}
        assert "self_heal" in names, "orca_sp_opi must expose self_heal"
        spec = next(p for p in OrcaSpOpiAdapter.params if p.name == "self_heal")
        assert spec.default == 0, "self_heal must be off by default"

    def test_opt_adapter_exposes_self_heal_off(self):
        from delfin.tools.adapters.orca_opt_freq_opi import OrcaOptOpiAdapter

        spec = next(p for p in OrcaOptOpiAdapter.params if p.name == "self_heal")
        assert spec.default == 0

    def test_opt_freq_adapter_exposes_self_heal_off(self):
        from delfin.tools.adapters.orca_opt_freq_opi import OrcaOptFreqOpiAdapter

        spec = next(p for p in OrcaOptFreqOpiAdapter.params if p.name == "self_heal")
        assert spec.default == 0


# ---------------------------------------------------------------------------
# Result reading — the retry log surfaces typed per-attempt records.
# ---------------------------------------------------------------------------


class TestResultReading:
    def test_retry_log_records_attempts(self):
        attempts = heal.build_attempt_log([
            heal.plan_escalation(heal.SelfHealFailure.SCF_NOT_CONVERGED, 1, job_kind="SP"),
            heal.plan_escalation(heal.SelfHealFailure.SCF_NOT_CONVERGED, 2, job_kind="SP"),
        ])
        assert len(attempts) == 2
        assert attempts[0]["failure"] == "scf_not_converged"
        assert attempts[0]["step"] == 1
        assert attempts[0]["keywords"]  # non-empty escalation
        assert attempts[1]["step"] == 2
        # Second is a strictly later rung than the first by its MaxIter.
        assert (attempts[1]["scf_options"].get("MaxIter")
                or attempts[1]["scf_options"].get("maxiter")) > \
            (attempts[0]["scf_options"].get("MaxIter")
             or attempts[0]["scf_options"].get("maxiter"))


# ---------------------------------------------------------------------------
# Probe — ORCA's own parser accepts the escalated input (skips when absent)
# ---------------------------------------------------------------------------


class TestProbe:
    def test_probe_accepts_escalated_opt_input(self, tmp_path, elements, coords):
        from delfin.orca_knowledge.probe import probe_input
        from delfin.tools.opi_chain import build_orca_input

        plan = heal.plan_escalation(
            heal.SelfHealFailure.GEOMETRY_NOT_CONVERGED, 2, job_kind="OPT")
        geom = heal.compose_overrides(None, plan.geom_options)
        inp = build_orca_input(
            tmp_path,
            job_kind="OPT",
            elements=elements,
            coordinates=coords,
            charge=0,
            mult=1,
            functional="B3LYP",
            basis="def2-SVP",
            scf_options=heal.compose_overrides(None, plan.scf_options),
            geom_options=geom,
            moread=None,
        )
        result = probe_input(inp.read_text())
        if result["verdict"] == "unavailable":
            pytest.skip("ORCA binary not available for probe")
        assert result["verdict"] in ("parsed", "finished"), result


# ---------------------------------------------------------------------------
# Input build — typical + edge
# ---------------------------------------------------------------------------


class TestInputBuild:
    def test_mathced_opt_geometry_block_emitted(self, tmp_path, elements, coords):
        from delfin.tools.opi_chain import build_orca_input

        plan = heal.plan_escalation(
            heal.SelfHealFailure.GEOMETRY_NOT_CONVERGED, 1, job_kind="OPT")
        inp = build_orca_input(
            tmp_path,
            job_kind="OPT",
            elements=elements,
            coordinates=coords,
            charge=0,
            mult=1,
            functional="B3LYP",
            basis="def2-SVP",
            scf_options=heal.compose_overrides(None, plan.scf_options),
            geom_options=heal.compose_overrides(None, plan.geom_options),
        )
        text = inp.read_text()
        # %geom block with the trust radius and the escalated MaxIter.
        assert "%geom" in text.lower()
        assert "maxiter" in text.lower()

    def test_sp_escalated_scf_block_emitted(self, tmp_path, elements, coords):
        from delfin.tools.opi_chain import build_orca_input

        plan = heal.plan_escalation(
            heal.SelfHealFailure.SCF_NOT_CONVERGED, 1, job_kind="SP")
        inp = build_orca_input(
            tmp_path,
            job_kind="OPT",  # build-only; the ladder's SCF surface is what matters
            elements=elements,
            coordinates=coords,
            charge=0,
            mult=1,
            functional="B3LYP",
            basis="def2-SVP",
            scf_options=plan.scf_options,
        )
        text = inp.read_text()
        assert "%scf" in text.lower()
        assert "maxiter" in text.lower()
