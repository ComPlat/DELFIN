"""OPI building block: geometry optimisation with an MLIP via ORCA's ExtOpt.

The block (``delfin.tools.adapters.orca_mlip_opt_opi``) drives an ORCA
``! ExtOpt`` optimisation whose energy / gradients come from an MLIP through
the orca-mlips command, reusing the shared input builder in
``delfin.mlp_tools.orca_mlips.build_orca_mlips_input`` (the same route the GOAT
block uses) and reading results back through ``delfin.orca_read``.

Manual ground truth (ORCA 6.1): 4.1.10 "Optimizing with External Methods" --
``!ExtOpt`` "must be combined with other keywords like !Opt or !GOAT"; the
``%method`` block carries ``ProgExt`` and ``Ext_Params``.  The ExtOpt energy is
the MLIP's "current total energy in Eh" so the FINAL SINGLE POINT ENERGY line
in the ``.out`` is reported in Hartree (``Eh``).

Build-only tests import OPI lazily so they run in every py3.11 CI job without
OPI; the read-back and registration tests run wherever OPI is installed.
"""

from __future__ import annotations

import pytest

WATER_XYZ = [
    ("O", 0.0, 0.0, 0.0),
    ("H", 0.757, 0.586, 0.0),
    ("H", -0.757, 0.586, 0.0),
]
WATER_ELEMENTS = [a[0] for a in WATER_XYZ]
WATER_COORDS = [a[1:] for a in WATER_XYZ]


# ---------------------------------------------------------------------------
# Build + refusal tests: never launch ORCA.  OPI imported lazily so they pass
# in the py3.11 no-OPI CI jobs too.
# ---------------------------------------------------------------------------


class TestBuildValidation:
    """The block refuses bad inputs with a named reason; never fixes chemistry."""

    def test_builds_plain_opt_extopt_input(self, tmp_path):
        from delfin.tools.adapters.orca_mlip_opt_opi import build_mlip_opt_input

        path = build_mlip_opt_input(
            tmp_path,
            run_type="opt",
            elements=WATER_ELEMENTS,
            coordinates=WATER_COORDS,
        )
        text = path.read_text().lower()
        # the shared builder writes simple keywords on separate lines when OPI
        # is available ("!opt"/"!extopt") and on one line in the plain-text
        # fallback writer ("! opt extopt"); tokenise on '!' and whitespace so
        # the "opt" and "extopt" run keywords hold for both forms
        tokens = set(text.replace("!", " ").split())
        assert "opt" in tokens
        assert "extopt" in tokens
        assert 'progext "orca-mlips-uma"' in text
        assert "--backend uma --device cpu" in text
        assert "*xyzfile 0 1 extopt.xyz" in text

    def test_build_sets_charge_mult(self, tmp_path):
        from delfin.tools.adapters.orca_mlip_opt_opi import build_mlip_opt_input

        path = build_mlip_opt_input(
            tmp_path,
            run_type="opt",
            elements=["O", "H", "H"],
            coordinates=WATER_COORDS,
            charge=-1,
            mult=2,
        )
        assert "*xyzfile -1 2 extopt.xyz" in path.read_text()

    def test_model_passes_through(self, tmp_path):
        from delfin.tools.adapters.orca_mlip_opt_opi import build_mlip_opt_input

        path = build_mlip_opt_input(
            tmp_path,
            run_type="opt",
            elements=WATER_ELEMENTS,
            coordinates=WATER_COORDS,
            model="a_model",
        )
        assert "--model a_model" in path.read_text()

    def test_refuses_unknown_run_type(self, tmp_path):
        from delfin.tools.adapters.orca_mlip_opt_opi import (
            build_mlip_opt_input as b,
        )

        with pytest.raises(ValueError, match="run_type"):
            b(tmp_path, run_type="md", elements=WATER_ELEMENTS,
              coordinates=WATER_COORDS)

    def test_refuses_unknown_backend(self, tmp_path):
        from delfin.tools.adapters.orca_mlip_opt_opi import (
            build_mlip_opt_input as b,
        )

        with pytest.raises(ValueError, match="not an orca-mlips plugin"):
            b(tmp_path, run_type="opt", elements=WATER_ELEMENTS,
              coordinates=WATER_COORDS, backend="nonesuch")

    def test_refuses_mlip_not_installed_with_hint(
        self, monkeypatch
    ):
        import delfin.mlp_tools.orca_mlips as om
        from delfin.tools.adapters.orca_mlip_opt_opi import validate_mlip_opt

        monkeypatch.setattr(
            om, "detect_orca_mlips",
            lambda: {
                "available": False,
                "progext_names": [],
                "version": "",
                "install_hint": ("pip install orca-mlips[uma] "
                                 "(or an email install hint)"),
                "xtb_required_for_solvent": True,
            },
        )
        # The not-installed refusal belongs at run time (check_installed=True),
        # never at input-build time: building the .inp needs no MLIP present.
        errors = validate_mlip_opt(
            run_type="opt", backend="uma", device="cpu",
            elements=WATER_ELEMENTS,
        )
        assert any("orca-mlips" in e for e in errors), errors
        assert any("pip install orca-mlips[uma]" in e for e in errors), errors
        assert validate_mlip_opt(
            run_type="opt", backend="uma", device="cpu",
            elements=WATER_ELEMENTS, check_installed=False,
        ) == []

    def test_refuses_element_not_covered(self, tmp_path):
        # aimnet2 is the one orca-mlips backend whose registry element list is
        # enumerated (the others are model-dependent / universal): Fe is not in
        # the AIMNet2 element set.
        from delfin.tools.adapters.orca_mlip_opt_opi import (
            build_mlip_opt_input as b,
        )

        with pytest.raises(ValueError, match="Fe"):
            b(
                tmp_path,
                run_type="opt",
                elements=["Fe", "O"],
                coordinates=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                backend="aimnet2",
            )

    def test_refuses_gpu_local(self, tmp_path):
        from delfin.tools.adapters.orca_mlip_opt_opi import (
            build_mlip_opt_input as b,
        )

        with pytest.raises(ValueError, match="GPU|cluster"):
            b(
                tmp_path,
                run_type="opt",
                elements=WATER_ELEMENTS,
                coordinates=WATER_COORDS,
                device="cuda",
            )

    def test_refuses_impossible_mult(self, tmp_path):
        from delfin.tools.adapters.orca_mlip_opt_opi import (
            build_mlip_opt_input as b,
        )

        with pytest.raises(ValueError, match="multiplicity"):
            b(
                tmp_path,
                run_type="opt",
                elements=WATER_ELEMENTS,
                coordinates=WATER_COORDS,
                mult=0,
            )


# ---------------------------------------------------------------------------
# Result reading: synthetic ORCA-format text written in the test (never real
# ORCA output).  No OPI needed -- reads through delfin.orca_read.
# ---------------------------------------------------------------------------


class TestResultReading:
    SYNTH_HEADER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output"
    INP = (
        "! Opt ExtOpt\n"
        "%method\n"
        '  ProgExt "orca-mlips-uma"\n'
        '  Ext_Params "--backend uma --device cpu"\n'
        "end\n"
        "* xyzfile 0 1 extopt.xyz\n"
    )
    CONVERGED_OUT = (
        "FINAL SINGLE POINT ENERGY       -76.026632737078\n"
        "\n"
        "        **********************\n"
        "        ** OPTIMIZATION RUN **\n"
        "        **********************\n"
        "\n"
        "        THE OPTIMIZATION HAS CONVERGED\n"
        "\n"
        "****ORCA TERMINATED NORMALLY****\n"
    )

    def test_reads_final_energy_geometry_convergence(self, tmp_path):
        from delfin.tools.adapters.orca_mlip_opt_opi import read_mlip_opt

        (tmp_path / "extopt.out").write_text(self.SYNTH_HEADER + "\n" + self.CONVERGED_OUT)
        (tmp_path / "extopt.inp").write_text(self.SYNTH_HEADER + "\n" + self.INP)
        (tmp_path / "extopt.xyz").write_text(
            "3\nwater\nO 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n"
        )
        res = read_mlip_opt(tmp_path)
        assert res["energy_Eh"] == pytest.approx(-76.026632737078)
        assert res["energy_unit"] == "Eh"
        assert res["geometry"] == tmp_path / "extopt.xyz"
        assert res["geometry_converged"] is True
        assert res["terminated_normally"] is True

    def test_reads_unconverged_as_not_converged(self, tmp_path):
        from delfin.tools.adapters.orca_mlip_opt_opi import read_mlip_opt

        (tmp_path / "extopt.inp").write_text(self.SYNTH_HEADER + "\n" + self.INP)
        (tmp_path / "extopt.out").write_text(
            self.SYNTH_HEADER + "\n"
            "FINAL SINGLE POINT ENERGY       -76.0\n"
            "        THE OPTIMIZATION HAS CONVERGED\n"
        )
        (tmp_path / "extopt.xyz").write_text("1\nx\nO 0.0 0.0 0.0\n")
        res = read_mlip_opt(tmp_path)
        assert res["terminated_normally"] is False
        assert res["error"] is not None


# ---------------------------------------------------------------------------
# ORCA probe.  The shared probe helper guards against external-program
# constructs (delfin/orca_knowledge/probe.py:184-187 forbids extopt, progsext,
# ext_params) because a probe must never let ORCA invoke an external program.
# So an ExtOpt input always yields verdict "refused" *before* ORCA launches;
# the parser-acceptance check cannot be exercised through the probe by design.
# This is the same limitation test_opi_orca_mlips_goat.py documents for the
# GOAT ExtOpt input (measured once manually against ORCA 6.1.1).  The test
# below exercises the probe plumbing and asserts the guard behaviour; it is
# skipped with a reason where ORCA is absent.
# ---------------------------------------------------------------------------


class TestOrcaProbe:
    def test_probe_guards_extopt_with_reason(self):
        import tempfile

        pytest.importorskip("opi")
        from delfin.orca_knowledge.probe import probe_input
        from delfin.tools.adapters.orca_mlip_opt_opi import build_mlip_opt_input

        with tempfile.TemporaryDirectory() as wd:
            inp = build_mlip_opt_input(
                wd,
                run_type="opt",
                elements=WATER_ELEMENTS,
                coordinates=WATER_COORDS,
            )
            text = inp.read_text()
        result = probe_input(text, timeout_s=300)
        if result.get("verdict") == "unavailable":
            pytest.skip("ORCA binary not available for probe")
        # ExtOpt is a guarded external-program construct: the probe must refuse
        # it before ORCA runs, never invoke the external MLIP.
        assert result.get("verdict") in ("refused",)


# ---------------------------------------------------------------------------
# Registration in the step-adapter registry.
# ---------------------------------------------------------------------------


class TestRegistration:
    def test_block_is_registered(self):
        from delfin.tools._registry import get as get_adapter

        adapter = get_adapter("orca_mlip_opt_opi")
        assert adapter is not None
        assert adapter.name == "orca_mlip_opt_opi"
        assert "opi" in adapter.requires_python


def test_a_refused_request_is_a_failed_step_not_a_crash(tmp_path):
    """Real-run finding: the refusal path used ErrorKind.INVALID_INPUT, which
    does not exist, so every refusal raised AttributeError."""
    from delfin.tools.adapters.orca_mlip_opt_opi import OrcaMlipOptOpiAdapter

    res = OrcaMlipOptOpiAdapter().execute(
        str(tmp_path), elements=["O", "H", "H"],
        coordinates=[[0, 0, 0.117], [0, 0.757, -0.469], [0, -0.757, -0.469]],
        backend="xyzzy",
    )
    assert res.status.name == "FAILED"
    assert "not an orca-mlips plugin" in res.error
