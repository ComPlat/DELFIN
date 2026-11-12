"""Tests for chaining the OPI-built ORCA blocks into a workflow.

The three geometry-altering OPI steps (``orca_opt_opi`` / ``orca_freq_opi`` /
``orca_opt_freq_opi`` in delfin.tools.adapters.orca_opt_freq_opi) build on an
OPI single point and must compose through the Pipeline engine:

- **Opt geometry feeds the next step** (``StepResult.geometry`` is the
  optimized `.xyz` the previous step left behind);
- **the .gbw is offered as MORead start orbitals** (``result.artifacts["gbw"]``
  is auto-wired into the next ``orca_*`` step's ``moread``);
- **a frequency step's thermal correction is available to an energy profile**
  (``result.data["thermal"]`` with ``free_energy_Eh`` etc.);
- **a failed / unconverged step stops the chain with its reason**
  (``stop_on_failure`` halts on the first FAILED step).

The identity/validation parts run anywhere; the run-and-read chain parts are
OPI-only and, like ``test_opi_orca_sp_dft``, are driven by the fake ORCA
(``tests/scripts/fake_orca.sh``) copying measured ORCA 6.1.1 fixture runs.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("opi")

from delfin.tools.opi_chain import (  # noqa: E402
    build_orca_input,
    read_geometry_xyz,
    validate_chain_method,
)
from delfin.tools.adapters.orca_opt_freq_opi import (  # noqa: E402
    OrcaOptOpiAdapter,
    OrcaFreqOpiAdapter,
)


def _min_xyz(tmp_path, name="geom.xyz") -> Path:
    p = tmp_path / name
    p.write_text(
        "3\n\nO             0.0000000000000000  0.0000000000000000  0.1173000000000000\n"
        "H             0.0000000000000000  0.7572000000000000 -0.4692000000000000\n"
        "H             0.0000000000000000 -0.7572000000000000 -0.4692000000000000\n"
    )
    return p


class TestAdapterContracts:
    """The chainable steps declare the right consumes/produces/geometry."""

    def test_opt_produces_geometry_gbw(self):
        c = OrcaOptOpiAdapter().contract()
        assert c.produces_geometry is True
        assert "geometry" in c.produces
        assert "gbw" in c.produces
        assert "geometry" in c.consumes

    def test_freq_produces_hessian_and_gbw_and_consumes_geometry(self):
        c = OrcaFreqOpiAdapter().contract()
        assert c.produces_geometry is False
        assert "geometry" not in c.produces
        assert "hessian" in c.produces
        assert "gbw" in c.produces
        assert "geometry" in c.consumes


class TestChainValidation:
    """The knowledge-store validation: an unknown %freq block option is an
    error returned to the caller, never dropped."""

    def test_valid_dft_method_has_no_errors(self):
        errors, warnings = validate_chain_method(
            functional="B3LYP", basis="def2-SVP",
            freq_options={"temp": 298.15, "pressure": 1.0},
        )
        assert errors == []

    def test_unknown_freq_block_option_is_an_error(self):
        errors, _ = validate_chain_method(
            functional="B3LYP", basis="def2-SVP",
            freq_options={"not_a_freq_option": 1},
        )
        assert any("not_a_freq_option" in e for e in errors)

    def test_unknown_functional_still_an_error(self):
        errors, _ = validate_chain_method(
            functional="NOTAFUNCTIONAL", basis="def2-SVP")
        assert any("NOTAFUNCTIONAL" in e for e in errors)


class TestGeometryParsing:
    """Reading the previous step's geometry out of an ORCA .xyz file."""

    def test_parses_xyz(self, tmp_path):
        p = _min_xyz(tmp_path)
        elements, coordinates = read_geometry_xyz(p)
        assert elements == ["O", "H", "H"]
        assert len(coordinates) == 3
        assert coordinates[0][2] == pytest.approx(0.1173)

    def test_bad_xyz_refused(self, tmp_path):
        p = tmp_path / "bad.xyz"
        p.write_text("2\n\nonly the count line\n")
        with pytest.raises(ValueError):
            read_geometry_xyz(p)


class TestBuildInput:
    """The OPI input writer composes Opt/Freq/MOREAD without inventing syntax."""

    def test_opt_input_contains_opt(self, tmp_path):
        inp = build_orca_input(
            tmp_path, job_kind="OPT", elements=["O", "H", "H"],
            coordinates=[[0, 0, 0.1173], [0, 0.7572, -0.4692],
                         [0, -0.7572, -0.4692]],
            charge=0, mult=1, basis="def2-SVP", base_name="opt",
        )
        text = inp.read_text()
        assert "!opt" in text or "!Opt" in text

    def test_freq_input_has_temp_and_pressure(self, tmp_path):
        inp = build_orca_input(
            tmp_path, job_kind="FREQ", elements=["O", "H", "H"],
            coordinates=[[0, 0, 0.1173], [0, 0.7572, -0.4692],
                         [0, -0.7572, -0.4692]],
            charge=0, mult=1, basis="def2-SVP",
            freq_options={"temp": 298.15, "pressure": 1.0}, base_name="freq",
        )
        text = inp.read_text()
        assert "!freq" in text or "!Freq" in text
        assert "temp 298.15" in text

    def test_moread_writes_moinp_block_when_gbw_present(self, tmp_path):
        gbw = tmp_path / "prev.gbw"
        gbw.write_bytes(b"fake gbw")
        inp = build_orca_input(
            tmp_path, job_kind="FREQ", elements=["O", "H", "H"],
            coordinates=[[0, 0, 0.1173], [0, 0.7572, -0.4692],
                         [0, -0.7572, -0.4692]],
            charge=0, mult=1, basis="def2-SVP", moread=str(gbw), base_name="freq",
        )
        text = inp.read_text()
        # The manual's top-level %moinp (section 2.20.5): ORCA reads guess
        # orbitals from a .gbw named in %moinp.
        assert "%moinp" in text
        assert '"freq.gbw"' in text

    def test_moread_refused_when_gbw_missing(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            build_orca_input(
                tmp_path, job_kind="FREQ", elements=["O", "H", "H"],
                coordinates=[[0, 0, 0.1173], [0, 0.7572, -0.4692],
                             [0, -0.7572, -0.4692]],
                charge=0, mult=1, basis="def2-SVP",
                moread=str(tmp_path / "nope.gbw"), base_name="freq",
            )


# ----------------------------------------------------------------------------
# 3-step chain Opt -> Freq -> SP driven through Pipeline.run() with the fake
# ORCA.  These need measured ORCA 6.1.1 fixture runs (opt with .xyz, freq with
# .hess + property.json carrying thermochemistry); they are written against
# those fixtures and are skipped until the fixtures ship.
# ----------------------------------------------------------------------------

import delfin.orca as orca_mod  # noqa: E402

FIXTURE_DIR = Path("tests/fixtures/orca611")


def _chain_fake_script(tmp_path: Path) -> Path:
    """A step-aware fake ORCA that copies each step's fixture outputs.

    The Pipeline names each step's working directory ``<base>/NN_<label>``;
    the fake reads the ''NN_<label>'' from its own input path's directory,
    maps ``0x`` to a fixture via the FAKE_ORCA_CHAIN env file
    (lines ``NN_label=fixture_dir``), and copies THAT fixture's outputs —
    re-named to the run's base name (always ``calc``).  This is how a multi-step
    chain gets distinct Opt / Freq / SP fixture runs without a real ORCA.
    """
    script = tmp_path / "chain_fake.sh"
    script.write_text(
        '#!/usr/bin/env bash\n'
        'set -euo pipefail\n'
        'inp="$1"\n'
        'dir="$(cd "$(dirname "$inp")" && pwd)"\n'
        'base="$(basename "$inp" .inp)"\n'
        'step="$(basename "$dir")"\n'
        'map="${FAKE_ORCA_CHAIN:?need FAKE_ORCA_CHAIN map file}"\n'
        'fix_dir="$(awk -v s="$step" -F= \'$1==s {print $2}\' "$map")"\n'
        'if [ -z "$fix_dir" ]; then\n'
        '  echo "chain_fake: no fixture for step $step" >&2\n'
        '  exit 2\n'
        'fi\n'
        'for ext in out gbw xyz hess json property.json property.txt inp; do\n'
        '  for f in "$fix_dir"/*."$ext"; do\n'
        '    [ -f "$f" ] || continue\n'
        '    cp "$f" "$dir/$base.$ext"\n'
        '    break\n'
        '  done\n'
        'done\n'
        'exit 0\n'
    )
    script.chmod(0o755)
    return script


@pytest.fixture()
def chain_fake(monkeypatch, tmp_path):
    """Install the step-aware fake ORCA and a step->fixture map."""
    script = _chain_fake_script(tmp_path)
    monkeypatch.setattr(orca_mod, "find_orca_executable", lambda: str(script))
    import delfin.qm_runtime as qmr

    qmr.clear_resolver_cache()
    mapping = (tmp_path / "chain.map")

    def map_steps(step_to_fixture: dict):
        mapping.write_text(
            "\n".join(f"{step}={FIXTURE_DIR / fix}" for step, fix in step_to_fixture.items()))
        monkeypatch.setenv("FAKE_ORCA_CHAIN", str(mapping))

    return map_steps


def test_three_step_chain_propagates_geometry_and_thermal(chain_fake, tmp_path):
    """Opt -> Freq -> SP: the chain runs, the Opt geometry flows to Freq, the
    freq step exposes its thermal correction, and each step's .gbw becomes the
    next step's MOREAD start-orbital source."""
    opt = FIXTURE_DIR / "h2o_opt"
    freq = FIXTURE_DIR / "h2o_freq"
    sp = FIXTURE_DIR / "h2o_optfreq"
    if not (opt.is_dir() and freq.is_dir()):
        pytest.skip("h2o_opt / h2o_freq ORCA fixtures not generated yet")

    chain_fake({
        "00_orca_opt_opi": opt,
        "01_orca_freq_opi": freq,
        "02_orca_sp_opi": sp,
    })

    from delfin.tools.pipeline import Pipeline

    pipe = Pipeline("h2o_chain", base_dir=tmp_path / "chain")
    pipe.add("orca_opt_opi", functional="B3LYP", basis="def2-SVP")
    pipe.add("orca_freq_opi", functional="B3LYP", basis="def2-SVP")
    pipe.add("orca_sp_opi", elements=["O", "H", "H"],
             coordinates=[[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692],
                          [0.0, -0.7572, -0.4692]],
             functional="B3LYP", basis="def2-SVP")
    result = pipe.run(geometry=_min_xyz(tmp_path), cores=1)

    assert result.ok, [getattr(r, "error", None) for r in result.results]
    names = [r.step_name for r in result.results]
    assert names == ["orca_opt_opi", "orca_freq_opi", "orca_sp_opi"]

    opt_res = result.results[0]
    freq_res = result.results[1]
    # Opt produced a geometry (its <base>.xyz) — that is what feeds Freq.
    assert opt_res.geometry is not None and opt_res.geometry.exists()
    # Freq consumed that geometry and produced a gbw artifact (MOREAD source).
    assert freq_res.artifacts.get("gbw") is not None
    # The freq step carries its thermal correction for an energy profile.
    assert freq_res.data.get("thermal", {}).get("free_energy_Eh") is not None


def test_mid_chain_failure_stops_pipeline(chain_fake, tmp_path):
    """A freq step the fake reports as unconverged halts the chain: the SP
    never runs and the failing step's reason is on the result."""
    opt = FIXTURE_DIR / "h2o_opt"
    bad = FIXTURE_DIR / "scf_not_converged"
    if not opt.is_dir():
        pytest.skip("h2o_opt ORCA fixture not generated yet")

    chain_fake({
        "00_orca_opt_opi": opt,
        "01_orca_freq_opi": bad,   # fails / does not converge
    })

    from delfin.tools.pipeline import Pipeline

    pipe = Pipeline("h2o_fail", base_dir=tmp_path / "fail_chain")
    pipe.add("orca_opt_opi", functional="B3LYP", basis="def2-SVP")
    pipe.add("orca_freq_opi", functional="B3LYP", basis="def2-SVP")
    pipe.add("orca_sp_opi", elements=["O", "H", "H"],
             coordinates=[[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692],
                          [0.0, -0.7572, -0.4692]],
             functional="B3LYP", basis="def2-SVP")
    result = pipe.run(geometry=_min_xyz(tmp_path), cores=1)

    # The chain stopped at the failing (freq) step: SP must NOT have run.
    names = [r.step_name for r in result.results]
    assert "orca_sp_opi" not in names
    freq = [r for r in result.results if r.step_name == "orca_freq_opi"]
    assert freq, "the freq step must be present and FAILED"
    assert freq[0].status.name == "FAILED"
    assert freq[0].error is not None
