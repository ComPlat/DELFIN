from delfin.analysis_tools import collect_analysis_summary
from delfin.ensemble_nmr import (
    build_anmrrc_text,
    build_censo_anmr_rc,
    build_orca_reference_input,
    normalize_anmr_solvent_name,
    normalize_cpcm_solvent_name,
    xyz_body_to_coord_text,
)
from delfin.qm_runtime import settings_selectable_tools


def test_xyz_body_to_coord_text_converts_to_turbomole_coord():
    coord_text = xyz_body_to_coord_text(
        [
            "H 0.000000 0.000000 0.000000",
            "C 1.000000 0.000000 0.000000",
        ]
    )
    assert coord_text.startswith("$coord\n")
    assert "  h" in coord_text
    assert "1.88972599" in coord_text
    assert coord_text.rstrip().endswith("$end")


def test_build_censo_anmr_rc_contains_orca_nmr_setup():
    text = build_censo_anmr_rc(
        solvent="chcl3",
        resonance_frequency=500.0,
        active_nuclei=("h",),
        orca_path="/opt/orca/orca",
        xtb_path="/opt/xtb/xtb",
    )
    assert "[nmr]" in text
    assert "prog = orca" in text
    assert "run = True" in text
    assert "sm_s = cpcm" in text
    assert "sm_j = cpcm" in text
    assert "resonance_frequency = 500.0" in text
    assert "h_active = True" in text
    assert "c_active = False" in text
    assert "orcapath = /opt/orca/orca" in text
    assert "xtbpath = /opt/xtb/xtb" in text
    assert "solvent = chloroform" in text
    assert "[refinement]" in text
    assert "threshold = 0.0" in text


def test_build_censo_anmr_rc_emits_orca_version_when_parseable():
    text = build_censo_anmr_rc(
        solvent="chcl3",
        resonance_frequency=400.0,
        active_nuclei=("h",),
        orca_path="/opt/orca_6_1_1/orca",
        xtb_path="/opt/xtb/xtb",
    )
    assert "orcaversion = 6.1.1" in text


def test_build_anmrrc_text_uses_computed_reference_values():
    text = build_anmrrc_text(
        solvent="chcl3",
        resonance_frequency=400.0,
        shielding_ref_h=31.234567,
        shielding_ref_c=184.876543,
    )
    assert "ENSO qm= ORCA mf= 400.0" in text
    assert "TMS[chcl3]" in text
    assert "1  31.234567    0.0    1" in text
    assert "6  184.876543    0.0    0" in text


def test_build_orca_reference_input_uses_cpcm_and_nmr():
    text = build_orca_reference_input(
        "Si 0.0 0.0 0.0\nC 1.0 0.0 0.0\nH 1.0 1.0 0.0",
        solvent="chcl3",
        pal=8,
        maxcore=2000,
    )
    assert "! PBE0 D4 DEF2-TZVP NMR CPCM(chloroform)" in text
    assert "nprocs 8" in text
    assert "%maxcore 2000" in text
    assert "* xyz 0 1" in text


def test_solvent_normalization_uses_full_names_for_cpcm_and_short_names_for_anmr():
    assert normalize_cpcm_solvent_name("chcl3") == "chloroform"
    assert normalize_cpcm_solvent_name("h2o") == "water"
    assert normalize_cpcm_solvent_name("ch2cl2") == "dichloromethane"
    assert normalize_anmr_solvent_name("chloroform") == "chcl3"


def test_settings_selectable_tools_include_censo_anmr_helpers():
    selectable = settings_selectable_tools()
    assert "censo" in selectable
    assert "anmr" in selectable
    assert "c2anmr" in selectable
    assert "nmrplot" in selectable


def test_analysis_summary_lists_censo_anmr_helpers():
    summary = collect_analysis_summary()
    tool_names = {item["name"] for item in summary["tools"]}
    assert "CENSO" in tool_names
    assert "ANMR" in tool_names
    assert "c2anmr" in tool_names
    assert "nmrplot" in tool_names
