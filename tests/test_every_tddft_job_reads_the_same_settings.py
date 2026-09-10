"""Every %tddft block DELFIN writes reads CONTROL through one place.

Before delfin.common.tddft_settings, eleven writers read the TD-DFT keys
five different ways.  A user who set TDDFT_nroots=25 and TDDFT_TDA=FALSE got
them in the TDDFT-mode state optimisations and nowhere else: the deltaSCF
check jobs, the standalone S0 check and all ISC/IC/FLUOR/PHOSP rate jobs kept
15 roots and TDA, and TDDFT_TDDFT_maxiter -- the template's own spelling --
reached no job at all.  These tests generate every ESD input for every mode
and read the blocks back.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin.common import tddft_settings
from delfin.config import (
    get_esd_hints,
    parse_control_text,
    read_control_file,
    set_control_value,
    validate_control_text,
)
from delfin.define import TEMPLATE

_XYZ = "3\nwater\nO 0.000 0.000 0.117\nH 0.000 0.757 -0.467\nH 0.000 -0.757 -0.467\n"

_BASE = [
    ("charge", "0"), ("solvent", "water"), ("method", "classic"),
    ("ESD_modul", "yes"), ("ESD_modus", "TDDFT"), ("ESD_T1_opt", "uks"),
    ("states", "[S1,T1,S2,T2]"), ("ISCs", "[S1>T1,T1>S1]"), ("ICs", "[S2>S1,T2>T1]"),
    ("emission_rates", "[f,p]"),
]


def _control(**keys) -> str:
    text = TEMPLATE
    for key, value in _BASE + list(keys.items()):
        text = set_control_value(text, key, value)
    return text


def _generate_all(tmp_path: Path, monkeypatch, mode: str, **keys) -> dict:
    """Every ESD input one run can write, keyed by file name."""
    import delfin.esd_input_generator as gen
    import delfin.esd_module as esd_module

    work = tmp_path / mode
    esd = work / "ESD"
    esd.mkdir(parents=True)
    monkeypatch.chdir(work)
    (work / "initial.xyz").write_text(_XYZ)
    for stem in ("S0", "S1", "S2", "T1", "T2", "S1_second_deltaSCF", "S2_second_deltaSCF",
                 "T2_second_deltaSCF", "S1_second", "S2_second", "T2_second"):
        (esd / f"{stem}.xyz").write_text(_XYZ)
        (esd / f"{stem}.out").write_text("FINAL SINGLE POINT ENERGY      -76.0\n")
    control = work / "CONTROL.txt"
    control.write_text(_control(ESD_modus=mode, **keys))
    config = read_control_file(str(control))
    common = dict(esd_dir=esd, charge=0, solvent="water", metals=[],
                  main_basisset="def2-SVP", metal_basisset="def2-TZVP", config=config)
    for state in ("S0", "S1", "S2", "T1", "T2"):
        gen.create_state_input(state=state, **common)
    for pair in ("S1>T1", "T1>S1"):
        gen.create_isc_input(isc_pair=pair, trootssl=0, **common)
    for pair in ("S2>S1", "T2>T1"):
        gen.create_ic_input(ic_pair=pair, **common)
    gen.create_fluor_input(**common)
    gen.create_phosp_input(**common)
    esd_module._create_s0_tddft_check_input(output_path=esd / "S0_TDDFT_standalone.inp", **common)
    return {f.name: f.read_text() for f in esd.glob("*.inp")}


def _blocks(text: str) -> list:
    """Each %tddft block as {keyword: value}, lower-cased."""
    out = []
    for body in re.findall(r"%tddft\n(.*?)\nend\b", text, flags=re.I | re.S):
        out.append({line.split()[0].lower(): line.split()[1].lower()
                    for line in body.splitlines() if len(line.split()) >= 2})
    return out


@pytest.mark.parametrize("mode", ["TDDFT", "deltaSCF", "hybrid1"])
def test_every_writer_honours_the_tddft_keys(tmp_path, monkeypatch, mode):
    inputs = _generate_all(
        tmp_path, monkeypatch, mode,
        TDDFT_nroots="25", TDDFT_maxdim="8", TDDFT_maxiter="321",
        TDDFT_TDA="FALSE", TDDFT_SOC="true",
    )
    seen = 0
    for name, text in inputs.items():
        for block in _blocks(text):
            seen += 1
            assert block.get("nroots") == "25", (name, block)
            assert block.get("maxdim") == "8", (name, block)
            assert block.get("maxiter") == "321", (name, block)
            assert block.get("tda") == "false", (name, block)
            if "ISC" in name or "PHOSP" in name:
                assert block.get("dosoc") == "true", (name, block)  # the rate needs SOC
            elif "FLUOR" in name or "_IC" in name:
                assert "dosoc" not in block, (name, block)
            else:
                assert block.get("dosoc") == "true", (name, block)  # TDDFT_SOC
    # every kind of job wrote at least one block
    assert seen >= 12
    assert any("ISC" in n for n in inputs) and any("PHOSP" in n for n in inputs)


def test_maxdim_auto_is_orcas_own_multiplier(tmp_path, monkeypatch):
    inputs = _generate_all(tmp_path, monkeypatch, "TDDFT")
    maxdims = {block["maxdim"] for text in inputs.values() for block in _blocks(text)}
    # ORCA 6.0.1 and 6.1.1 use 10 x NRoots when the input names no MaxDim
    assert maxdims == {"10"}
    assert tddft_settings.ORCA_DEFAULT_MAXDIM == 10


def test_template_maxiter_reaches_the_jobs(tmp_path, monkeypatch):
    inputs = _generate_all(tmp_path, monkeypatch, "deltaSCF")
    values = {block.get("maxiter") for text in inputs.values() for block in _blocks(text)}
    assert values == {"500"}


def test_followiroot_only_where_a_root_is_optimised(tmp_path, monkeypatch):
    inputs = _generate_all(tmp_path, monkeypatch, "TDDFT")
    assert _blocks(inputs["S1.inp"])[0].get("followiroot") == "true"
    for name in ("S1_T1_ISC_ms0.inp", "S1_S0_FLUOR.inp", "S0_TDDFT_standalone.inp"):
        assert all("followiroot" not in b for b in _blocks(inputs[name])), name


def test_additions_reach_every_block(tmp_path, monkeypatch):
    inputs = _generate_all(
        tmp_path, monkeypatch, "hybrid1",
        TDDFT_additions="DoNTO true; OrbWin[0] 2,-1,-1,14; ETol 1e-7",
    )
    for name, text in inputs.items():
        for body in re.findall(r"%tddft\n(.*?)\nend\b", text, flags=re.I | re.S):
            assert "  DoNTO true" in body, name
            # the comma the CONTROL reader split at is back
            assert "  OrbWin[0] 2,-1,-1,14" in body, name
            assert "  ETol 1e-7" in body, name


@pytest.mark.parametrize("line, fragment", [
    ("nroots 30", "use TDDFT_nroots"),
    ("MaxDim 5", "use TDDFT_maxdim"),
    ("DoSOC true", "use TDDFT_SOC"),
    ("iroot 2", "per job"),
    ("triplets true", "per job"),
    ("%scf maxiter 300 end", "not blocks"),
])
def test_additions_refuse_what_has_another_spelling(line, fragment):
    errors = validate_control_text(_control(TDDFT_additions=line))
    assert any(fragment in err for err in errors), errors


def test_validator_takes_auto_and_refuses_typos():
    assert validate_control_text(_control(TDDFT_maxdim="auto", TDDFT_maxiter="auto")) == []
    errors = validate_control_text(_control(TDDFT_maxdim="lots"))
    assert any("TDDFT_maxdim must be auto or a positive whole number" in e for e in errors)
    # a yes/no reading would make this "false" and switch on full TD-DFT silently
    errors = validate_control_text(_control(TDDFT_TDA="TURE"))
    assert any("TDDFT_TDA must be true or false" in e for e in errors)


@pytest.mark.parametrize("spelling, value, canonical, expected", [
    ("TDDFT_TDDFT_maxiter", "400", "TDDFT_maxiter", 400),
    ("ESD_TDDFT_maxiter", "400", "TDDFT_maxiter", 400),
    ("ESD_nroots", "20", "TDDFT_nroots", 20),
    ("ESD_maxdim", "6", "TDDFT_maxdim", 6),
    ("ESD_SOC", "yes", "TDDFT_SOC", "true"),
    ("tddft_nroots", "20", "TDDFT_nroots", 20),
])
def test_old_spellings_land_on_the_one_key(tmp_path, spelling, value, canonical, expected):
    text = _control(ESD_modus="TDDFT") + f"\n{spelling}={value}\n"
    # the template line for the canonical key would win the tie; drop it
    text = re.sub(rf"(?m)^{canonical}=.*\n", "", text)
    control = tmp_path / "CONTROL.txt"
    control.write_text(text)
    config = read_control_file(str(control))
    assert config[canonical] == expected
    assert spelling not in config or spelling == canonical


def test_an_explicit_maxdim_outside_orcas_range_is_explained():
    hints = get_esd_hints(_control(TDDFT_maxdim="30"))
    assert any("450 vectors" in h and "auto" in h for h in hints), hints
    assert not any("TDDFT_maxdim" in h for h in get_esd_hints(_control(TDDFT_maxdim="auto")))
    assert not any("TDDFT_maxdim" in h for h in get_esd_hints(_control(TDDFT_maxdim="7")))


def test_recovery_asks_for_a_multiplier_not_a_vector_count():
    from delfin.orca_recovery import OrcaErrorType, RecoveryStrategy

    strategy = RecoveryStrategy(OrcaErrorType.TDDFT_ROOT_COLLAPSE, attempt=3, config={})
    fixes = strategy._tddft_instability_fixes()
    assert fixes["tddft_block"]["maxdim"] == 2 * tddft_settings.ORCA_DEFAULT_MAXDIM


def test_hand_built_configs_fall_back_to_the_old_names():
    config = {"ESD_nroots": 12, "ESD_maxdim": 6, "TDA": "FALSE", "ESD_SOC": "yes",
              "ESD_TDDFT_maxiter": 250}
    settings = tddft_settings.read_settings(config)
    assert (settings.nroots, settings.maxdim, settings.tda, settings.soc, settings.maxiter) == \
        (12, 6, False, True, 250)
    parsed = parse_control_text("TDDFT_nroots=\nESD_nroots=9\n")
    # a blank canonical key must not hide the value under the old name
    assert tddft_settings.read_settings(parsed).nroots == 9
