"""Unit tests for delfin.cli_fukui — argparse + helpers + dispatch logic."""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

from delfin import cli_fukui


# ---------------------------------------------------------------------------
# input-type detection
# ---------------------------------------------------------------------------

def test_looks_like_xyz_path_true(tmp_path):
    p = tmp_path / "mol.xyz"
    p.write_text("1\ncomment\nC 0 0 0\n")
    assert cli_fukui._looks_like_xyz_path(str(p))


def test_looks_like_xyz_path_false_for_missing_file(tmp_path):
    assert not cli_fukui._looks_like_xyz_path(str(tmp_path / "nope.xyz"))


def test_looks_like_smiles_basic():
    assert cli_fukui._looks_like_smiles("C=O")
    assert cli_fukui._looks_like_smiles("CCO")
    assert cli_fukui._looks_like_smiles("C1=CC=CC=C1")


def test_looks_like_smiles_rejects_path_like(tmp_path):
    existing = tmp_path / "real.xyz"
    existing.write_text("1\n\nC 0 0 0\n")
    assert not cli_fukui._looks_like_smiles(str(existing))
    assert not cli_fukui._looks_like_smiles("path/to/file.xyz")


def test_looks_like_smiles_rejects_empty():
    assert not cli_fukui._looks_like_smiles("")


# ---------------------------------------------------------------------------
# XYZ header helpers
# ---------------------------------------------------------------------------

def test_ensure_xyz_header_adds_when_missing():
    headerless = "C 0.0 0.0 0.0\nO 1.2 0.0 0.0\n"
    out = cli_fukui._ensure_xyz_header(headerless)
    lines = out.splitlines()
    assert lines[0] == "2"
    assert "fukui" in lines[1].lower()
    assert "C" in lines[2]
    assert "O" in lines[3]


def test_ensure_xyz_header_keeps_existing():
    with_header = "2\ntitle\nC 0 0 0\nO 1.2 0 0\n"
    assert cli_fukui._ensure_xyz_header(with_header).startswith("2\n")


def test_strip_xyz_header_removes_header():
    with_header = "2\ntitle\nC 0 0 0\nO 1.2 0 0\n"
    stripped = cli_fukui._strip_xyz_header(with_header)
    assert stripped == "C 0 0 0\nO 1.2 0 0"


def test_strip_xyz_header_handles_headerless():
    headerless = "C 0 0 0\nO 1.2 0 0"
    assert cli_fukui._strip_xyz_header(headerless) == "C 0 0 0\nO 1.2 0 0"


# ---------------------------------------------------------------------------
# Workdir + pre-opt resolution
# ---------------------------------------------------------------------------

def test_resolve_workdir_explicit(tmp_path):
    target = tmp_path / "explicit"
    result = cli_fukui._resolve_workdir("C=O", str(target))
    assert result == target.resolve()


def test_resolve_workdir_from_xyz(tmp_path, monkeypatch):
    p = tmp_path / "mymol.xyz"
    p.write_text("1\n\nC 0 0 0\n")
    monkeypatch.chdir(tmp_path)
    result = cli_fukui._resolve_workdir(str(p), None)
    assert result.name == "fukui_mymol"


def test_resolve_workdir_from_smiles_sanitizes(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = cli_fukui._resolve_workdir("C=O", None)
    assert result.name == "fukui_C_O"


def test_resolve_pre_opt_forced_for_smiles():
    assert cli_fukui._resolve_pre_opt("C=O", False) is True
    assert cli_fukui._resolve_pre_opt("C=O", None) is True
    assert cli_fukui._resolve_pre_opt("C=O", True) is True


def test_resolve_pre_opt_optional_for_xyz(tmp_path):
    p = tmp_path / "mol.xyz"
    p.write_text("1\n\nC 0 0 0\n")
    assert cli_fukui._resolve_pre_opt(str(p), False) is False
    assert cli_fukui._resolve_pre_opt(str(p), None) is False
    assert cli_fukui._resolve_pre_opt(str(p), True) is True


# ---------------------------------------------------------------------------
# Settings resolution
# ---------------------------------------------------------------------------

def test_resolve_setting_cli_wins():
    out = cli_fukui._resolve_setting("CLI", {"main_basisset": "control"}, "main_basisset", "default")
    assert out == "CLI"


def test_resolve_setting_control_when_no_cli():
    out = cli_fukui._resolve_setting(None, {"main_basisset": "control"}, "main_basisset", "default")
    assert out == "control"


def test_resolve_setting_default_when_neither():
    out = cli_fukui._resolve_setting(None, {}, "main_basisset", "default")
    assert out == "default"


def test_resolve_setting_empty_string_falls_back():
    out = cli_fukui._resolve_setting("", {"main_basisset": "control"}, "main_basisset", "default")
    assert out == "control"


# ---------------------------------------------------------------------------
# ORCA input assembly (uses DELFIN's canonical bang-line builder)
# ---------------------------------------------------------------------------

_BASE_SETTINGS = {
    "functional": "B3LYP",
    "basis": "def2-SVP",
    "dispersion": "D3BJ",
    "solvation": "",
    "solvent": "",
    "pal": "4",
    "maxcore": "3000",
}


def test_build_orca_input_minimal_sp(tmp_path):
    config = cli_fukui._resolve_fukui_config(tmp_path, dict(_BASE_SETTINGS))
    inp = cli_fukui._build_orca_input(
        job_type="SP", charge=0, mult=1,
        coords="C 0 0 0\nO 1.2 0 0",
        config=config,
        found_metals=[],
    )
    bang = inp.splitlines()[0]
    assert bang.startswith("!")
    assert "B3LYP" in bang
    assert "def2-SVP" in bang
    assert "UKS" not in bang
    # SP has neither FREQ nor geom_opt token
    assert "FREQ" not in bang.upper()
    assert "OPT" not in bang.upper()
    assert "%maxcore 3000" in inp
    assert "%pal nprocs 4 end" in inp
    assert "* xyz 0 1" in inp
    assert "C 0 0 0" in inp
    assert inp.rstrip().endswith("*")


def test_build_orca_input_anion_marks_uks(tmp_path):
    config = cli_fukui._resolve_fukui_config(tmp_path, dict(_BASE_SETTINGS))
    inp = cli_fukui._build_orca_input(
        job_type="SP", charge=-1, mult=2,
        coords="C 0 0 0",
        config=config,
        found_metals=[],
    )
    assert "UKS" in inp.splitlines()[0]
    assert "* xyz -1 2" in inp


def test_build_orca_input_solvation_passes_through(tmp_path):
    settings = dict(_BASE_SETTINGS, solvation="CPCM", solvent="water")
    config = cli_fukui._resolve_fukui_config(tmp_path, settings)
    inp = cli_fukui._build_orca_input(
        job_type="SP", charge=0, mult=1,
        coords="C 0 0 0",
        config=config,
        found_metals=[],
    )
    assert "CPCM(water)" in inp.splitlines()[0]


def test_build_orca_input_opt_includes_freq_and_geom_opt(tmp_path):
    settings = dict(_BASE_SETTINGS)
    config = cli_fukui._resolve_fukui_config(tmp_path, settings)
    config["geom_opt"] = "TIGHTOPT"
    inp = cli_fukui._build_orca_input(
        job_type="OPT", charge=0, mult=1,
        coords="C 0 0 0",
        config=config,
        found_metals=[],
    )
    bang = inp.splitlines()[0].upper()
    assert "FREQ" in bang
    assert "TIGHTOPT" in bang


def test_resolve_fukui_config_reads_control_txt(tmp_path):
    (tmp_path / "CONTROL.txt").write_text(
        "functional=PBE0\nmain_basisset=def2-TZVP\ndisp_corr=D4\n"
        "implicit_solvation_model=CPCM\nsolvent=acetonitrile\n"
        "PAL=8\nmaxcore=2500\n",
        encoding="utf-8",
    )
    # CLI provides no overrides → CONTROL.txt wins
    cli_settings = {
        "functional": "", "basis": "", "dispersion": "",
        "solvation": "", "solvent": "", "pal": "", "maxcore": "",
    }
    # Manually mimic the CLI's resolve step
    control = cli_fukui._maybe_load_control(tmp_path)
    merged = {
        "functional": cli_fukui._resolve_setting("", control, "functional", "B3LYP"),
        "basis":      cli_fukui._resolve_setting("", control, "main_basisset", "def2-SVP"),
        "dispersion": cli_fukui._resolve_setting("", control, "disp_corr", "D3BJ"),
        "solvation":  cli_fukui._resolve_setting("", control, "implicit_solvation_model", ""),
        "solvent":    cli_fukui._resolve_setting("", control, "solvent", ""),
        "pal":        cli_fukui._resolve_setting("", control, "PAL", "4"),
        "maxcore":    cli_fukui._resolve_setting("", control, "maxcore", "3000"),
    }
    assert merged["functional"] == "PBE0"
    assert merged["basis"] == "def2-TZVP"
    assert merged["dispersion"] == "D4"
    assert merged["solvation"] == "CPCM"
    assert merged["solvent"] == "acetonitrile"
    assert merged["pal"] == "8"
    assert merged["maxcore"] == "2500"


# ---------------------------------------------------------------------------
# argparse smoke
# ---------------------------------------------------------------------------

def test_parser_input_optional_now():
    """``--input`` is optional: dashboard pre-seeds input.txt in the workdir."""
    parser = cli_fukui._build_parser()
    args = parser.parse_args([])
    assert args.input is None


def test_parser_accepts_smiles_input():
    parser = cli_fukui._build_parser()
    args = parser.parse_args(["--input", "C=O"])
    assert args.input == "C=O"
    assert args.scheme == "mulliken"
    assert args.pre_opt is None  # tri-state default


def test_parser_pre_opt_toggle():
    parser = cli_fukui._build_parser()
    args = parser.parse_args(["--input", "C=O", "--pre-opt"])
    assert args.pre_opt is True
    args = parser.parse_args(["--input", "C=O", "--no-pre-opt"])
    assert args.pre_opt is False


def test_parser_scheme_validation():
    parser = cli_fukui._build_parser()
    args = parser.parse_args(["--input", "C=O", "--scheme", "loewdin"])
    assert args.scheme == "loewdin"
    with pytest.raises(SystemExit):
        parser.parse_args(["--input", "C=O", "--scheme", "mayer"])


# ---------------------------------------------------------------------------
# End-to-end with mocked ORCA — verifies dispatch + result writing
# ---------------------------------------------------------------------------

_FAKE_OUT = """\
-----------------------
MULLIKEN ATOMIC CHARGES
-----------------------
   0 C    :    {c:.6f}
   1 O    :    {o:.6f}
   2 H    :    {h1:.6f}
   3 H    :    {h2:.6f}
Sum of atomic charges:    0.0000000
"""


def _make_fake_run_orca(out_text_per_state: dict[str, str]):
    """Return a stub for delfin.orca.run_orca that writes a fake .out per state."""
    def stub(inp_path, out_path, *args, **kwargs):
        name = Path(out_path).stem
        Path(out_path).write_text(out_text_per_state.get(name, "no charges block\n"))
        return True
    return stub


def test_main_e2e_with_mocked_orca(tmp_path, monkeypatch):
    """Run main() against a pre-optimized .xyz, with ORCA mocked.

    Verifies the full pipeline assembles 3 jobs, parses their charges,
    and writes fukui_result.json + .fukui_job marker.
    """
    xyz = tmp_path / "mol.xyz"
    xyz.write_text("4\nformaldehyde\nC 0 0 0\nO 1.2 0 0\nH -0.5 0.9 0\nH -0.5 -0.9 0\n")
    workdir = tmp_path / "fukui_out"

    fakes = {
        "neutral": _FAKE_OUT.format(c=0.10, o=-0.30, h1=0.10, h2=0.10),
        "anion":   _FAKE_OUT.format(c=-0.20, o=-0.40, h1=0.05, h2=0.05),
        "cation":  _FAKE_OUT.format(c=0.20, o=0.50, h1=0.15, h2=0.15),
    }

    with mock.patch("delfin.orca.run_orca", side_effect=_make_fake_run_orca(fakes)):
        rc = cli_fukui.main([
            "--input", str(xyz),
            "--workdir", str(workdir),
            "--no-pre-opt",
            "--skip-cubes",
        ])

    assert rc == 0
    marker = workdir / ".fukui_job"
    assert marker.exists()
    result = workdir / "fukui_result.json"
    assert result.exists()

    import json
    payload = json.loads(result.read_text())
    assert payload["scheme"] == "mulliken"
    assert payload["atoms"] == ["C", "O", "H", "H"]
    # f_plus on C should be largest (constructed test data)
    assert payload["f_plus"][0] == max(payload["f_plus"])
    # f_minus on O should be largest
    assert payload["f_minus"][1] == max(payload["f_minus"])


def test_main_e2e_orca_failure_exits_nonzero(tmp_path):
    xyz = tmp_path / "mol.xyz"
    xyz.write_text("1\n\nC 0 0 0\n")

    with mock.patch("delfin.orca.run_orca", return_value=False):
        with pytest.raises(SystemExit) as exc:
            cli_fukui.main([
                "--input", str(xyz),
                "--workdir", str(tmp_path / "out"),
                "--no-pre-opt",
                "--skip-cubes",
            ])
        assert "failed" in str(exc.value).lower()
