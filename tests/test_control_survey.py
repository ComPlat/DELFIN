"""Tests for the CONTROL.txt survey (synthetic fixtures only).

No archive file is ever copied into the repository: every fixture here is
written by the tests themselves.
"""

import json
import os

import pytest

from delfin.orca_knowledge import control_survey as cs


@pytest.fixture()
def run_dir(tmp_path):
    """A synthetic DELFIN run: CONTROL.txt, xyz, ORCA input + outputs."""
    d = tmp_path / "run1"
    d.mkdir()
    control = d / "CONTROL.txt"
    control.write_text(
        "# synthetic fixture\n"
        "NAME=Test System\n"
        "charge=1\n"
        "multiplicity_global_opt=2\n"
        "method=classic\n"
        "functional=B3LYP\n"
        "main_basisset=def2-TZVP\n"
        "implicit_solvation_model=CPCM\n"
        "solvent=ACN\n"
        "PAL=4\n"
        "GUPPY=yes\n"
        "totally_unknown_key=somevalue\n"
    )
    # xyz with a metal (Ni), NO real coordinates, dummy values only.
    xyz = d / "input.xyz"
    xyz.write_text("5\nsynthetic fixture - not a real structure\n"
                   "Ni 0.0 0.0 0.0\nC 0.1 0.0 0.0\nN 0.0 0.1 0.0\n"
                   "O 0.0 0.0 0.1\nH 0.1 0.1 0.0\n")
    inp = d / "input.inp"
    inp.write_text("! B3LYP OPT FREQ\n* xyz 1 2\nNi 0 0 0\nC 0 0 1\nN 0 1 0\n"
                   "O 1 0 0\nH 1 1 0\n*\n")
    (d / "output.out").write_text("some orca output\nORCA TERMINATED NORMALLY\n")
    return d


def test_parse_returns_every_known_key_and_unknown_separately(run_dir):
    parsed = cs.parse(str(run_dir / "CONTROL.txt"))
    config = parsed["config"]
    unknown = parsed["unknown"]
    assert config["method"] == "classic"
    assert config["functional"] == "B3LYP"
    assert config["main_basisset"] == "def2-TZVP"
    assert config["implicit_solvation_model"] == "CPCM"
    assert config["solvent"] == "ACN"
    assert config["GUPPY"] == "yes"
    assert unknown == {"totally_unknown_key": "somevalue"}
    # Every key DELFIN knows must not leak into "unknown".
    assert "charge" in config and "PAL" in config
    # Aliases are canonicalised: NAME is in CONTROL_FIELD_SPECS as NAME.
    assert config["NAME"] == "Test System"


def test_parse_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        cs.parse(str(tmp_path / "nope.txt"))


def test_parse_refuses_oversized_control(tmp_path):
    big = tmp_path / "CONTROL.txt"
    big.write_bytes(b"x" * (cs.MAX_FILE_SIZE + 1))
    with pytest.raises(ValueError):
        cs.parse(str(big))


def test_system_counts_metals_and_uses_delfin_charge(run_dir):
    sysd = cs.system(str(run_dir))
    assert sysd["n_atoms"] == 5
    assert sysd["metals"] == {"Ni": 1}
    assert sysd["elements"] == ["C", "H", "N", "Ni", "O"]
    assert sysd["charge"] == 1
    assert sysd["charge_source"] == "CONTROL.txt"
    assert sysd["multiplicity"] == 2
    assert sysd["multiplicity_source"] == "CONTROL.txt:multiplicity_global_opt"


def test_system_falls_back_to_orca_recorded_charge(tmp_path):
    d = tmp_path / "run2"
    d.mkdir()
    (d / "CONTROL.txt").write_text("charge=0\nmethod=classic\n")
    # CONTROL's charge wins; multiplicity only from ORCA input.
    (d / "input.xyz").write_text("3\ncomment\nFe 0 0 0\nCl 1 0 0\nCl 0 1 0\n")
    (d / "input.inp").write_text("! PBE\n* xyz 0 5\nFe 0 0 0\nCl 1 0 0\nCl 0 1 0\n*\n")
    sysd = cs.system(str(d))
    assert sysd["charge"] == 0 and sysd["charge_source"] == "CONTROL.txt"
    assert sysd["multiplicity"] == 5
    assert sysd["multiplicity_source"] == "orca_input"
    assert sysd["metals"] == {"Fe": 1}


def test_system_without_structure_leaves_n_atoms_none(tmp_path):
    d = tmp_path / "run3"
    d.mkdir()
    (d / "CONTROL.txt").write_text("charge=-2\n")
    sysd = cs.system(str(d))
    assert sysd["n_atoms"] is None
    assert sysd["elements"] == []
    assert sysd["charge"] == -2


def test_outcome_finished_when_summary_and_all_normal(run_dir):
    (run_dir / "DELFIN.txt").write_text("final summary\n")
    oc = cs.outcome(str(run_dir))
    assert oc["finished"] is True
    assert oc["orca_steps_total"] == 1
    assert oc["failed_steps"] == []


def test_outcome_failed_step_listed(run_dir):
    (run_dir / "DELFIN.txt").write_text("summary\n")
    (run_dir / "output2.out").write_text("SCF not converged\n****Maximum Opt iterations reached****\n")
    oc = cs.outcome(str(run_dir))
    assert oc["finished"] is False
    assert "output2" in oc["failed_steps"]


def test_outcome_aborted_step_distinguished(run_dir):
    (run_dir / "DELFIN.txt").write_text("summary\n")
    (run_dir / "output3.out").write_text("something\nABORTING THE RUN\n")
    oc = cs.outcome(str(run_dir))
    assert "output3" in oc["aborted_steps"]
    assert "output3" not in oc["failed_steps"]


def test_outcome_no_summary_and_no_outputs_is_unverifiable(tmp_path):
    d = tmp_path / "run4"
    d.mkdir()
    (d / "CONTROL.txt").write_text("charge=0\n")
    oc = cs.outcome(str(d))
    assert oc["finished"] is False
    assert oc["unverifiable"] is True


def test_outcome_reads_subdirectory_outputs(run_dir):
    sub = run_dir / "IMAG"
    sub.mkdir()
    (run_dir / "DELFIN.txt").write_text("s\n")
    (sub / "output.out").write_text("x\nORCA TERMINATED NORMALLY\n")
    (sub / "output2.out").write_text("x\nbroken\n")
    oc = cs.outcome(str(run_dir))
    assert oc["orca_steps_total"] == 3
    assert oc["failed_steps"] == ["IMAG/output2"]


def test_survey_writes_jsonl_with_hash_only(tmp_path):
    root = tmp_path / "archive"
    d = root / "someone" / "runX"
    d.mkdir(parents=True)
    (d / "CONTROL.txt").write_text("charge=0\nmethod=classic\n")
    (d / "input.xyz").write_text("2\nc\nFe 0 0 0\nCl 1 0 0\n")
    (d / "output.out").write_text("o\nORCA TERMINATED NORMALLY\n")
    (d / "DELFIN.txt").write_text("summary\n")
    out = tmp_path / "survey.jsonl"
    n = cs.survey(str(root), str(out))
    assert n == 1
    rec = json.loads(out.read_text().splitlines()[0])
    assert set(rec) == {"hash", "system", "control", "unknown", "outcome"}
    # The path never appears; the hash is a sha256 hex digest.
    assert rec["hash"] and len(rec["hash"]) == 64
    assert all(c in "0123456789abcdef" for c in rec["hash"])
    blob = out.read_text()
    assert "runX" not in blob and "someone" not in blob
    # No coordinate-like number: xyz coordinates were 0/1 only, but the record
    # must not contain them at all.
    assert "0.0" not in blob


def test_survey_skips_symlink_escape_and_oversize(tmp_path):
    root = tmp_path / "archive"
    d1 = root / "ok"
    d1.mkdir(parents=True)
    (d1 / "CONTROL.txt").write_text("charge=0\n")
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "CONTROL.txt").write_text("charge=9\n")
    link = root / "link_to_outside"
    os.symlink(outside, link)
    big = root / "big" / "CONTROL.txt"
    big.parent.mkdir()
    big.write_bytes(b"x" * (cs.MAX_FILE_SIZE + 1))
    out = tmp_path / "survey.jsonl"
    n = cs.survey(str(root), str(out))
    assert n == 1  # only the ok dir; symlink and oversize skipped


def test_cli_main_writes_file(tmp_path, capsys):
    root = tmp_path / "archive"
    d = root / "r"
    d.mkdir(parents=True)
    (d / "CONTROL.txt").write_text("charge=0\n")
    out = tmp_path / "s.jsonl"
    rc = cs.main([str(root), "--out", str(out)])
    assert rc == 0
    assert out.is_file()
    assert "1" in capsys.readouterr().out


def test_outcome_ignores_slurm_wrapper_and_old_outputs(run_dir):
    (run_dir / "DELFIN.txt").write_text("summary\n")
    (run_dir / "delfin_12345.out").write_text("slurm wrapper stdout, no ORCA line\n")
    (run_dir / "output.old1.out").write_text("earlier attempt, no termination\n")
    (run_dir / "output2.out").write_text("x\nORCA TERMINATED NORMALLY\n")
    oc = cs.outcome(str(run_dir))
    assert oc["orca_steps_total"] == 2
    assert oc["orca_steps_normal"] == 2
    assert oc["failed_steps"] == []
    assert oc["finished"] is True
