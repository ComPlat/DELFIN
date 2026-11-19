"""Tests for the relativistic / spin-property ORCA building block.

The building block (``delfin.tools.adapters.orca_relspin_opi``) builds three
single-job ORCA inputs entirely with OPI's typed objects:

* ``rel_sp``    — a scalar-relativistic single point (ZORA / DKH1 / DKH2)
  with the matching recontracted basis,
* ``zfs``       — the zero-field splitting via the ``%eprnmr`` ``DTensor``
  block (axial ``D``, rhombic ``E``, and their ratio),
* ``soc_tddft`` — SOC-corrected excited states via ``TDDFT`` with the
  ``%tddft`` ``soc`` flag.

Every refusal is asserted with its reason.  Result reading goes through
``delfin.orca_read`` (extended with ``read_zfs_parameters`` and
``read_soc_states``) and is tested against synthetic ORCA-format text written
here in the test -- never real ORCA output.
"""
import json

import pytest

# The main py3.11 CI jobs run ALL tests without OPI installed; skip there.
# The OPI job (DELFIN_OPI_REQUIRED=1, tests/conftest.py) turns that skip
# into a failure, so the OPI-required job still gates on these tests.
pytest.importorskip("opi")

from delfin.tools.adapters.orca_relspin_opi import (  # noqa: E402
    read_rel_sp_result,
    read_zfs_result,
    build_rel_sp_input,
    build_soc_tddft_input,
    build_zfs_input,
    OrcaRelSpinOpiAdapter,
    validate_relspin,
)

# A trivial geometry: one He atom at the origin.  Nothing CSD-derived, no
# crystal data -- just a placeholder point for an input-format test.
ELS = ["He"]
XYZ = [[0.0, 0.0, 0.0]]
# An open-shell placeholder for the ZFS tests: one H atom carries 1 electron,
# so the doublet (mult 2) is the physically allowed spin state.
ZFS_ELS = ["H"]


def _read_inp(path):
    return path.read_text()


# ---------------------------------------------------------------------------
# validate_relspin: keyword / block-option validation against the knowledge
# store, plus the physical checks (multiplicity, electron count, basis match).
# ---------------------------------------------------------------------------

class TestValidateRelspin:
    def test_rel_sp_zora_with_recontracted_basis_is_accepted(self):
        errors, warnings = validate_relspin(
            functional="B3LYP", job="rel_sp", hamiltonian="ZORA", basis="ZORA-def2-TZVP",
        )
        assert errors == []

    def test_rel_sp_dkh2_with_recontracted_basis_is_accepted(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="rel_sp", hamiltonian="DKH2", basis="DKH-def2-TZVP",
        )
        assert errors == []

    def test_unknown_hamiltonian_is_an_error(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="rel_sp", hamiltonian="ZORAX", basis="def2-TZVP",
        )
        assert any("ZORAX" in e for e in errors)

    def test_unknown_basis_is_an_error(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="rel_sp", hamiltonian="ZORA", basis="not-a-basis",
        )
        assert any("not-a-basis" in e for e in errors)

    def test_nonrel_basis_with_relativistic_hamiltonian_is_refused(self):
        """Manual 6.1 sec 2.12.6: ZORA/DKH need recontracted (ZORA-/DKH-)
        basis sets; a plain def2 basis silently changes the Hamiltonian."""
        errors, _ = validate_relspin(
            functional="B3LYP", job="rel_sp", hamiltonian="ZORA", basis="def2-TZVP",
        )
        assert errors, "non-relativistic basis with ZORA must be refused"
        assert any("ZORA" in e for e in errors)

    def test_nonrel_basis_allowed_when_explicitly_requested(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="rel_sp", hamiltonian="ZORA", basis="def2-TZVP",
            allow_nonrel_basis=True,
        )
        assert errors == []

    def test_unknown_job_is_an_error(self):
        errors, _ = validate_relspin(job="not-a-job", functional="B3LYP", basis="def2-TZVP")
        assert any("not-a-job" in e for e in errors)

    def test_zfs_with_singlet_is_refused(self):
        """A closed shell has no zero-field splitting (manual epr.txt: the
        ZFS tensor needs at least one unpaired electron)."""
        errors, _ = validate_relspin(
            functional="B3LYP", job="zfs", basis="def2-TZVP", mult=1,
        )
        assert any("multiplicity" in e.lower() for e in errors)

    def test_zfs_with_triplet_is_accepted(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="zfs", basis="def2-TZVP", mult=2, elements=ZFS_ELS,
        )
        assert errors == []

    def test_unknown_epr_option_is_an_error(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="zfs", basis="def2-TZVP", mult=2,
            epr_options={"dtensor": "nonsense"},
        )
        assert any("dtensor" in e.lower() for e in errors)

    def test_unknown_rel_option_is_an_error(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="zfs", basis="def2-TZVP", mult=2, elements=ZFS_ELS,
            rel_options={"nonsenserel": 1},
        )
        assert any("nonsenserel" in e.lower() for e in errors)

    def test_soc_tddft_with_zero_nroots_is_refused(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="soc_tddft", basis="def2-TZVP", mult=1, nroots=0,
        )
        assert any("nroots" in e.lower() for e in errors)

    def test_soc_tddft_typical_is_accepted(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="soc_tddft", basis="def2-TZVP", mult=1, nroots=8,
        )
        assert errors == []

    def test_even_electrons_with_even_multiplicity_is_refused(self):
        """He has 2 electrons; multiplicity 2 would be an even-electron
        doublet, which is impossible (parity of N-1 must match 2S)."""
        errors, _ = validate_relspin(
            functional="B3LYP", job="rel_sp", hamiltonian="ZORA", basis="ZORA-def2-TZVP",
            charge=0, mult=2, elements=ELS,
        )
        assert any("multiplicity" in e.lower() for e in errors)

    def test_charge_scales_the_electron_count(self):
        """He+ has 1 electron; multiplicity 2 is the physically possible
        doublet and must be accepted."""
        errors, _ = validate_relspin(
            functional="B3LYP", job="rel_sp", hamiltonian="ZORA", basis="ZORA-def2-TZVP",
            charge=1, mult=2,
        )
        assert errors == []

    def test_unknown_dsoc_value_is_refused(self):
        """OPI's BlockEprnmr model allows only qro|pk|cp|pcp|cvw."""
        errors, _ = validate_relspin(
            functional="B3LYP", job="zfs", basis="def2-TZVP", mult=2,
            epr_options={"dsoc": "qrx"},
        )
        assert any("dsoc" in e.lower() for e in errors)

    def test_dsoc_cpcm_with_triplet_is_accepted(self):
        errors, _ = validate_relspin(
            functional="B3LYP", job="zfs", basis="def2-TZVP", mult=2,
            epr_options={"dsoc": "cp"},
        )
        assert errors == []


# ---------------------------------------------------------------------------
# Input build: the written .inp must contain the keywords the manual names.
# ---------------------------------------------------------------------------

class TestBuildRelSpInput:
    def test_zora_builds_recontracted_basis_and_rel_block(self, tmp_path):
        path = build_rel_sp_input(
            tmp_path, functional="B3LYP", hamiltonian="ZORA", basis="ZORA-def2-TZVP",
            elements=ELS, coordinates=XYZ, charge=0, mult=1,
            maxcore=1000, ncores=2, base_name="zora_he",
        )
        text = _read_inp(path)
        assert "!zora-def2-tzvp" in text
        assert "!zora" in text
        assert "%rel" in text and "method zora" in text
        assert "nprocs 2" in text

    def test_dkh2_builds_dkh_basis_and_rel_block(self, tmp_path):
        path = build_rel_sp_input(
            tmp_path, functional="B3LYP", hamiltonian="DKH2", basis="DKH-def2-TZVP",
            elements=ELS, coordinates=XYZ, charge=0, mult=1,
            maxcore=1000, ncores=1, base_name="dkh_he",
        )
        text = _read_inp(path)
        assert "!dkh-def2-tzvp" in text
        assert "!dkh2" in text
        assert "method dkh" in text

    def test_picture_change_and_order_reach_the_rel_block(self, tmp_path):
        path = build_rel_sp_input(
            tmp_path, functional="B3LYP", hamiltonian="DKH2", basis="DKH-def2-TZVP",
            elements=ELS, coordinates=XYZ, charge=0, mult=1,
            rel_options={"picturechange": 2, "order": 4},
            maxcore=1000, ncores=1, base_name="dkh_pc",
        )
        text = _read_inp(path)
        assert "picturechange 2" in text
        assert "order 4" in text

    def test_nonrel_basis_raises(self, tmp_path):
        with pytest.raises(ValueError, match="def2-TZVP"):
            build_rel_sp_input(
                tmp_path, functional="B3LYP", hamiltonian="ZORA", basis="def2-TZVP",
                elements=ELS, coordinates=XYZ, charge=0, mult=1,
                maxcore=1000, ncores=1, base_name="bad",
            )

    def test_impossible_multiplicity_raises(self, tmp_path):
        with pytest.raises(ValueError, match="multiplicity"):
            build_rel_sp_input(
                tmp_path, functional="B3LYP", hamiltonian="ZORA", basis="ZORA-def2-TZVP",
                elements=ELS, coordinates=XYZ, charge=0, mult=2,
                maxcore=1000, ncores=1, base_name="bad",
            )

    def test_extra_keywords_are_validated_too(self, tmp_path):
        with pytest.raises(ValueError, match="NOTREAL"):
            build_rel_sp_input(
                tmp_path, functional="B3LYP", hamiltonian="ZORA", basis="ZORA-def2-TZVP",
                elements=ELS, coordinates=XYZ, charge=0, mult=1,
                extra_keywords=("NOTREAL",),
                maxcore=1000, ncores=1, base_name="bad",
            )


class TestBuildZfsInput:
    def test_triplet_builds_eprnmr_block(self, tmp_path):
        path = build_zfs_input(
            tmp_path, functional="PBE", basis="def2-TZVP", elements=ZFS_ELS, coordinates=XYZ,
            charge=0, mult=2, maxcore=1000, ncores=1, base_name="zfs_h",
        )
        text = _read_inp(path)
        assert "%eprnmr" in text
        assert "dtensor ssandso" in text
        assert "nprocs 1" in text

    def test_epr_options_reach_the_block(self, tmp_path):
        path = build_zfs_input(
            tmp_path, functional="PBE", basis="def2-TZVP", elements=ZFS_ELS, coordinates=XYZ,
            charge=0, mult=2, epr_options={"dsoc": "cp"},
            maxcore=1000, ncores=1, base_name="zfs_h",
        )
        text = _read_inp(path)
        assert "dsoc cp" in text

    def test_zora_hamiltonian_adds_rel_block(self, tmp_path):
        path = build_zfs_input(
            tmp_path, functional="PBE", basis="ZORA-def2-TZVP", hamiltonian="ZORA",
            elements=ZFS_ELS, coordinates=XYZ, charge=0, mult=2,
            maxcore=1000, ncores=1, base_name="zfs_zora",
        )
        text = _read_inp(path)
        assert "!zora-def2-tzvp" in text
        assert "!zora" in text
        assert "%rel" in text and "method zora" in text

    def test_singlet_raises(self, tmp_path):
        with pytest.raises(ValueError, match="multiplicity"):
            build_zfs_input(
                tmp_path, functional="PBE", basis="def2-TZVP", elements=ZFS_ELS,
                coordinates=XYZ,
                charge=0, mult=1, maxcore=1000, ncores=1, base_name="bad",
            )

    def test_unknown_epr_option_raises(self, tmp_path):
        with pytest.raises(ValueError, match="dtensor"):
            build_zfs_input(
                tmp_path, functional="PBE", basis="def2-TZVP", elements=ZFS_ELS,
                coordinates=XYZ,
                charge=0, mult=2, epr_options={"dtensor": "nonsense"},
                maxcore=1000, ncores=1, base_name="bad",
            )


class TestBuildSocTddftInput:
    def test_typical_builds_tddft_block_with_soc(self, tmp_path):
        path = build_soc_tddft_input(
            tmp_path, functional="PBE", basis="def2-TZVP", nroots=8, elements=ELS,
            coordinates=XYZ, charge=0, mult=1,
            maxcore=1000, ncores=1, base_name="soc_he",
        )
        text = _read_inp(path)
        # TDDFT is selected by the %tddft block only: a "!tddft" simple
        # keyword is rejected by ORCA 6.1.1 ("UNRECOGNIZED OR DUPLICATED
        # KEYWORD(S) IN SIMPLE INPUT LINE TDDFT").
        assert "%tddft" in text and "nroots 8" in text
        assert not any(
            ln.startswith("!") and "tddft" in ln.lower()
            for ln in text.splitlines()
        )
        assert "%tddft" in text
        assert "nroots 8" in text
        assert "soc True" in text

    def test_nroots_zero_raises(self, tmp_path):
        with pytest.raises(ValueError, match="nroots"):
            build_soc_tddft_input(
                tmp_path, functional="PBE", basis="def2-TZVP", nroots=0, elements=ELS,
                coordinates=XYZ, charge=0, mult=1,
                maxcore=1000, ncores=1, base_name="bad",
            )

    def test_dkh2_adds_rel_block(self, tmp_path):
        path = build_soc_tddft_input(
            tmp_path, functional="PBE", basis="DKH-def2-TZVP", hamiltonian="DKH2", nroots=6,
            elements=ELS, coordinates=XYZ, charge=0, mult=1,
            maxcore=1000, ncores=1, base_name="soc_dkh",
        )
        text = _read_inp(path)
        assert "!dkh-def2-tzvp" in text
        assert "!dkh2" in text
        assert "%rel" in text and "method dkh" in text


# ---------------------------------------------------------------------------
# Result reading: synthetic ORCA-format text, never real output.
# ---------------------------------------------------------------------------

def _write_zfs_property(tmp_path, d, e):
    """A synthetic ORCA property-JSON: the D_Tensor record the ORCA manual
    documents in property_file.txt (D axial, E rhombic, cm^-1)."""
    prop = {"D_Tensor": {"D": d, "E": e, "d_raw": [d, 0.0, 0.0]}}
    (tmp_path / "zfs_he.property.json").write_text(json.dumps(prop))
    return tmp_path


class TestReadZfsParameters:
    def test_d_e_and_ratio_are_read(self, tmp_path):
        from delfin.orca_read import read_zfs_parameters

        _write_zfs_property(tmp_path, 1.23, 0.456)
        zfs = read_zfs_parameters(tmp_path, "zfs_he")
        assert zfs is not None
        assert zfs.d_cm1 == pytest.approx(1.23)
        assert zfs.e_cm1 == pytest.approx(0.456)
        assert zfs.e_over_d == pytest.approx(0.456 / 1.23)

    def test_missing_property_file_gives_none(self, tmp_path):
        from delfin.orca_read import read_zfs_parameters

        assert read_zfs_parameters(tmp_path, "nothere") is None


# Layout of a real ORCA 6.1.1 ``%tddft soc true`` run; numbers invented.
SOC_OUT = """# SYNTHETIC test fixture written for DELFIN -- not ORCA output

      SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
      Transition         Energy     Energy  Wavelength fosc(D2)      D2       |DX|      |DY|      |DZ|
                          (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)
  0-1.0A  ->  1-3.0A    1.530000   12340.3   810.4   0.000000010   0.00000   0.00000   0.00010   0.00000
  0-1.0A  ->  2-1.0A    2.909000   23462.6   426.2   0.020000000   0.10000   0.00000   0.31623   0.00000

                             ****ORCA TERMINATED NORMALLY****
"""


class TestReadSocStates:
    def test_states_are_read_from_the_spectrum(self, tmp_path):
        from delfin.orca_read import read_soc_states

        (tmp_path / "soc_he.out").write_text(SOC_OUT)
        states = read_soc_states(tmp_path, "soc_he")
        assert len(states) == 2
        assert states[0].e_cm1 == pytest.approx(12340.3)
        assert states[0].e_eV == pytest.approx(1.530)
        assert states[1].e_cm1 == pytest.approx(23462.6)
        assert states[1].fosc == pytest.approx(0.02)

    def test_missing_out_gives_empty_list(self, tmp_path):
        from delfin.orca_read import read_soc_states

        assert read_soc_states(tmp_path, "nothere") == []


# ---------------------------------------------------------------------------
# The adapter: typed result, registration in the tool registry.
# ---------------------------------------------------------------------------

class TestAdapterRegistration:
    def test_adapter_is_registered(self):
        from delfin.tools._registry import get

        assert get("orca_relspin_opi") is not None

    def test_adapter_builds_result_dataclass(self, tmp_path):
        from delfin.tools.adapters.orca_relspin_opi import OrcaRelSpinResult

        adapter = OrcaRelSpinOpiAdapter()
        assert adapter.name == "orca_relspin_opi"
        # The adapter's typed result carries the requested quantities.
        r = OrcaRelSpinResult(job="zfs", work_dir=tmp_path, ok=False)
        assert r.zfs_D_cm1 is None and r.zfs_E_cm1 is None


# ---------------------------------------------------------------------------
# Probe: ORCA's own parser must accept the input.  Parser check only, <= 2
# cores, never a full run; skipped with a reason where ORCA is absent
# (same shape as tests/test_opi_esd_probe.py).
# ---------------------------------------------------------------------------

def test_probe_orca_parser_accepts_the_relspin_input(tmp_path):
    from delfin.orca_knowledge.probe import probe_input

    path = build_rel_sp_input(
        tmp_path, functional="PBE", hamiltonian="ZORA", basis="ZORA-def2-TZVP",
        elements=ELS, coordinates=XYZ, charge=0, mult=1,
        maxcore=400, ncores=1, base_name="probe_rel",
    )
    verdict = probe_input(path.read_text(), timeout_s=300)
    if verdict["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert verdict["verdict"] in ("parsed", "finished"), (
        f"ORCA rejected the OPI-built relspin input; messages: "
        f"{verdict['messages']}"
    )


def test_probe_orca_parser_accepts_the_zfs_input(tmp_path):
    from delfin.orca_knowledge.probe import probe_input

    path = build_zfs_input(
        tmp_path, functional="PBE", basis="def2-TZVP", elements=ZFS_ELS, coordinates=XYZ,
        charge=0, mult=2, maxcore=400, ncores=1, base_name="probe_zfs",
    )
    verdict = probe_input(path.read_text(), timeout_s=300)
    if verdict["verdict"] == "unavailable":
        pytest.skip("ORCA binary not available for probe")
    assert verdict["verdict"] in ("parsed", "finished"), (
        f"ORCA rejected the OPI-built ZFS input; messages: "
        f"{verdict['messages']}"
    )


# ---------------------------------------------------------------------------
# Real-run findings (operator report): method required, no !tddft simple
# keyword, launch ignores the NMR-report False, OPI Output signature,
# probe inputs must not land in the checkout's scratch/.
# ---------------------------------------------------------------------------

class TestRealRunFindings:
    def test_rel_sp_requires_a_functional(self, tmp_path):
        with pytest.raises((TypeError, ValueError), match="functional"):
            build_rel_sp_input(
                tmp_path, hamiltonian="ZORA", basis="ZORA-def2-TZVP",
                elements=ELS, coordinates=XYZ, charge=0, mult=1,
                maxcore=1000, ncores=1, base_name="x",
            )

    def test_rel_sp_writes_the_functional(self, tmp_path):
        path = build_rel_sp_input(
            tmp_path, functional="B3LYP", hamiltonian="ZORA", basis="ZORA-def2-TZVP",
            elements=ELS, coordinates=XYZ, charge=0,
            mult=1, maxcore=1000, ncores=1, base_name="fn",
        )
        assert "!b3lyp" in _read_inp(path).lower()

    def test_zfs_writes_the_functional(self, tmp_path):
        path = build_zfs_input(
            tmp_path, basis="def2-TZVP", functional="PBE",
            elements=ZFS_ELS, coordinates=XYZ, charge=0, mult=2,
            maxcore=1000, ncores=1, base_name="zf",
        )
        assert "!pbe" in _read_inp(path).lower()

    def test_soc_tddft_writes_the_functional(self, tmp_path):
        path = build_soc_tddft_input(
            tmp_path, basis="def2-TZVP", functional="PBE", nroots=5,
            elements=ZFS_ELS, coordinates=XYZ, charge=0, mult=2,
            maxcore=1000, ncores=1, base_name="sc",
        )
        assert "!pbe" in _read_inp(path).lower()

    def test_soc_tddft_input_has_no_tddft_simple_keyword(self, tmp_path):
        # ORCA: "UNRECOGNIZED OR DUPLICATED KEYWORD(S) IN SIMPLE INPUT LINE
        # TDDFT" -- TDDFT is selected by the %tddft block only.
        path = build_soc_tddft_input(
            tmp_path, basis="def2-TZVP", functional="PBE", nroots=5,
            elements=ZFS_ELS, coordinates=XYZ, charge=0, mult=2,
            maxcore=1000, ncores=1, base_name="sc",
        )
        lines = [
            ln for ln in _read_inp(path).splitlines()
            if ln.startswith("!") and "tddft" in ln.lower()
        ]
        assert lines == []

    def test_launch_tolerates_a_false_return_from_the_runner(self, tmp_path, monkeypatch):
        # A %eprnmr input goes through run_orca's NMR-report path, which
        # returns False when no H shieldings exist -- the run itself finished
        # normally.  _launch must error only on a raised exception; the
        # readers gate on terminated_normally().
        import delfin.tools.adapters.orca_relspin_opi as mod

        def fake_orca(inp, out, working_dir=None, **_kw):
            return False
        monkeypatch.setattr("delfin.orca.run_orca", fake_orca)
        assert mod._launch(tmp_path, "j") is None

    def test_launch_reports_a_raised_exception(self, tmp_path, monkeypatch):
        import delfin.tools.adapters.orca_relspin_opi as mod

        def boom(inp, out, working_dir=None, **_kw):
            raise RuntimeError("orca gone")
        monkeypatch.setattr("delfin.orca.run_orca", boom)
        err = mod._launch(tmp_path, "j")
        assert err is not None and "orca gone" in err

    def test_read_rel_sp_result_never_imports_a_missing_reader(self, tmp_path):
        # The reader must not depend on a non-existent orca_read symbol;
        # a bogus .out yields ok=False, not an ImportError.
        (tmp_path / "j.out").write_text("# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n")
        res = read_rel_sp_result(tmp_path, "j")
        assert res.ok is False

    def test_read_zfs_result_never_imports_a_missing_reader(self, tmp_path):
        (tmp_path / "z.out").write_text(
            "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n")
        res = read_zfs_result(tmp_path, "z")
        assert res.ok is False


# ---------------------------------------------------------------------------
# Live ORCA 6.1.1 findings (real runs through the adapter's execute path).
# ---------------------------------------------------------------------------

_JOB_KWARGS = {
    "rel_sp": dict(hamiltonian="ZORA", basis="ZORA-def2-TZVP",
                   elements=ELS, coordinates=XYZ, charge=0, mult=1),
    "zfs": dict(basis="def2-TZVP", elements=ZFS_ELS, coordinates=XYZ,
                charge=0, mult=2),
    "soc_tddft": dict(basis="def2-TZVP", nroots=5, elements=ELS,
                      coordinates=XYZ, charge=0, mult=1),
}


class TestLiveRunFindings:
    @pytest.mark.parametrize("job", sorted(_JOB_KWARGS))
    def test_execute_passes_the_functional_to_the_input(
            self, job, tmp_path, monkeypatch):
        # The builders require a functional, so the step path must accept
        # one and hand it on; before, run_* had no such parameter and every
        # execute() call failed with an unexpected-keyword TypeError.
        import delfin.tools.adapters.orca_relspin_opi as mod

        monkeypatch.setattr(mod, "_launch", lambda wdir, base: "stop: no ORCA")
        res = OrcaRelSpinOpiAdapter().execute(
            tmp_path, cores=1, job=job, functional="B3LYP",
            **_JOB_KWARGS[job])
        assert "stop: no ORCA" in (res.error or ""), res.error
        assert "!b3lyp" in (tmp_path / f"{job}.inp").read_text().lower()

    @staticmethod
    def _opi_rejected_jsons(tmp_path, base):
        # A property JSON OPI 2.0.0's typed model rejects, as it rejects the
        # real ORCA 6.1.1 SOC-TDDFT one (SOC states carry multiplicity 0),
        # plus a gbw JSON so OPI never shells out to orca_2json.
        prop = {"Geometries": [{"Energy": [
            {"Method": "TDA/CIS", "Mult": [1, 3, 0]}]}]}
        (tmp_path / f"{base}.property.json").write_text(json.dumps(prop))
        (tmp_path / f"{base}.json").write_text("{}")

    def test_soc_tddft_result_does_not_depend_on_opis_property_model(
            self, tmp_path):
        from delfin.tools.adapters.orca_relspin_opi import read_soc_tddft_result

        (tmp_path / "soc.out").write_text(SOC_OUT)
        self._opi_rejected_jsons(tmp_path, "soc")
        res = read_soc_tddft_result(tmp_path, "soc")
        assert res.ok and res.terminated_normally and res.meaningful, res
        assert [s.states for s in res.soc_states] == [1, 2]

    def test_zfs_result_does_not_depend_on_opis_property_model(self, tmp_path):
        (tmp_path / "z.out").write_text(
            "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
            "                             ****ORCA TERMINATED NORMALLY****\n")
        self._opi_rejected_jsons(tmp_path, "z")
        prop = json.loads((tmp_path / "z.property.json").read_text())
        prop["Geometries"][0]["D_Tensor"] = [{"D": 2.0, "E": 0.2, "Mult": 3}]
        (tmp_path / "z.property.json").write_text(json.dumps(prop))
        res = read_zfs_result(tmp_path, "z")
        assert res.ok and res.meaningful, res
        assert res.zfs_D_cm1 == pytest.approx(2.0)
        assert res.zfs_E_over_D == pytest.approx(0.1)

    def test_a_failed_rel_sp_read_says_why(self, tmp_path):
        (tmp_path / "j.out").write_text(
            "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n")
        self._opi_rejected_jsons(tmp_path, "j")
        res = read_rel_sp_result(tmp_path, "j")
        assert res.ok is False
        assert res.error, "a failed read must carry its reason"

    def test_execute_reports_fosc_and_spin_label_of_each_soc_state(
            self, tmp_path, monkeypatch):
        import delfin.tools.adapters.orca_relspin_opi as mod

        def fake_launch(wdir, base):
            (wdir / f"{base}.out").write_text(SOC_OUT)
            TestLiveRunFindings._opi_rejected_jsons(wdir, base)
            return None
        monkeypatch.setattr(mod, "_launch", fake_launch)
        res = OrcaRelSpinOpiAdapter().execute(
            tmp_path, cores=1, job="soc_tddft", functional="B3LYP",
            **_JOB_KWARGS["soc_tddft"])
        assert res.ok, res.error
        first = res.data["soc_states"][1]
        assert first["fosc"] == pytest.approx(0.02)
        assert first["mult"] == pytest.approx(1.0)

    def test_soc_tddft_without_the_soc_table_is_a_failed_step(
            self, tmp_path, monkeypatch):
        # Fail closed: ORCA ended normally but no SOC-corrected states were
        # printed -> FAILED with the reason, never SUCCESS with no data.
        import delfin.tools.adapters.orca_relspin_opi as mod

        def fake_launch(wdir, base):
            (wdir / f"{base}.out").write_text(
                "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
                "      SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION "
                "ELECTRIC DIPOLE MOMENTS\n\n"
                "                             ****ORCA TERMINATED NORMALLY****\n")
            TestLiveRunFindings._opi_rejected_jsons(wdir, base)
            return None
        monkeypatch.setattr(mod, "_launch", fake_launch)
        res = OrcaRelSpinOpiAdapter().execute(
            tmp_path, cores=1, job="soc_tddft", functional="B3LYP",
            **_JOB_KWARGS["soc_tddft"])
        assert not res.ok
        assert res.status.value == "failed"
        assert "SOC-corrected states" in (res.error or ""), res.error

    def test_rel_sp_with_an_unconverged_scf_is_refused(self, tmp_path,
                                                       monkeypatch):
        import delfin.tools.adapters.orca_relspin_opi as mod

        class _Out:
            def terminated_normally(self):
                return True

            def get_final_energy(self):
                return -1.0

        monkeypatch.setattr(mod, "_read_output_with_reason",
                            lambda wdir, base: (_Out(), None))
        (tmp_path / "u.out").write_text(
            "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
            "                   SCF NOT CONVERGED AFTER 125 CYCLES\n"
            "                             ****ORCA TERMINATED NORMALLY****\n")
        res = read_rel_sp_result(tmp_path, "u")
        assert res.ok is False and res.energy_Eh is None
        assert "SCF did not converge" in (res.error or ""), res.error

    def test_rel_sp_with_a_converged_scf_reports_the_energy(self, tmp_path,
                                                            monkeypatch):
        import delfin.tools.adapters.orca_relspin_opi as mod

        class _Out:
            def terminated_normally(self):
                return True

            def get_final_energy(self):
                return -1.0

        monkeypatch.setattr(mod, "_read_output_with_reason",
                            lambda wdir, base: (_Out(), None))
        (tmp_path / "c.out").write_text(
            "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
            "               *           SCF CONVERGED AFTER  10 CYCLES          *\n"
            "                             ****ORCA TERMINATED NORMALLY****\n")
        res = read_rel_sp_result(tmp_path, "c")
        assert res.ok and res.energy_Eh == pytest.approx(-1.0)
