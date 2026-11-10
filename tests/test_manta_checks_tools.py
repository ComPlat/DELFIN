"""Tests for the agent-side MANTA checks (delfin/tools/manta_checks.py).

Fakes only: every "build" is a fake CLI -- a tiny standalone python file that
writes synthetic manifest.json + frame .xyz files (first line of every
fixture: "# SYNTHETIC test fixture written for DELFIN -- not ORCA output").
No real MANTA build runs in this suite and no ORCA binary is needed.
"""

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from delfin.tools import manta_checks as mc

ROOT = Path(__file__).resolve().parent.parent

FAKE_CLI = '''\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
# Fake MANTA CLI: writes two deterministic synthetic isomer frames.
import json, sys
from pathlib import Path
smi, out = sys.argv[1], Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
lines = smi.splitlines()
frames = []
for i, label in enumerate(("A", "B")):
    name = "%03d__%s.xyz" % (i, label)
    body = "2\\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\\n" \\
           "C 0.0 0.0 0.0\\nO 1.3 0.0 0.0\\n"
    (out / name).write_text(body)
    frames.append({"file": name, "label": label})
(out / "manifest.json").write_text(json.dumps({"isomers": frames}))
'''

FAKE_SMILES = "C(=O)O"


def _fake_cmd(tmp: Path) -> list:
    script = tmp / "fake_cli.py"
    script.write_text(FAKE_CLI)
    return [sys.executable, str(script)]


# ---------------------------------------------------------------- compare
def _ref_archive(tmp: Path, frames):
    arch = tmp / "arch"
    arch.mkdir()
    body = str(len(frames)) + "\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
    body += "".join("C %d.0 0.0 0.0\nO 1.3 %d.0 0.0\n" % (i, i)
                    for i in range(len(frames)))
    (arch / "sys1.xyz").write_text(body)
    return arch


def _build(tmp: Path, rid="sys1"):
    out = tmp / "out" / rid
    script = tmp / "fake_cli.py"
    script.write_text(FAKE_CLI)
    subprocess.run(
        [sys.executable, str(script), FAKE_SMILES, "-o", str(out)],
        check=True, capture_output=True)
    return out


def test_blocks_of_parses_multiframe():
    text = ("2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
            "C 0.0  0.0 0.0\nO 1.3 0.0 0.0\n"
            "2\nsecond\nN 0.0 0.0 0.0\nH 1.0 0.0 0.0\n")
    b = mc.blocks_of(text)
    assert len(b) == 4
    assert b[0] == "C 0.0 0.0 0.0"          # whitespace joined
    assert b[-1] == "H 1.0 0.0 0.0"


def test_compare_identical_other_order_different_missing(tmp_path):
    out = _build(tmp_path)
    ref = tmp_path / "arch_ident"
    ref.mkdir()
    (ref / "sys1.xyz").write_text(
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.3 0.0 0.0\n"
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.3 0.0 0.0\n")
    row = mc.compare_build_to_archive(out, ref / "sys1.xyz")
    assert row["verdict"] == "identical"

    # Same frame set, other frame order: two single-atom frames reordered.
    out2 = tmp_path / "out2" / "s2"
    out2.mkdir(parents=True)
    (out2 / "000__a.xyz").write_text(
        "1\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\nC 0.0 0.0 0.0\n")
    (out2 / "001__b.xyz").write_text(
        "1\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\nO 1.3 0.0 0.0\n")
    (out2 / "manifest.json").write_text(json.dumps(
        {"isomers": [{"file": "000__a.xyz", "label": "a"},
                     {"file": "001__b.xyz", "label": "b"}]}))
    ref2 = tmp_path / "arch2"
    ref2.mkdir()
    (ref2 / "s2.xyz").write_text(
        "1\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\nO 1.3 0.0 0.0\n"
        "1\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\nC 0.0 0.0 0.0\n")
    row2 = mc.compare_build_to_archive(out2, ref2 / "s2.xyz")
    assert row2["verdict"] == "same_set_other_order"

    # Different: extra coordinate line in the reference.
    ref3 = tmp_path / "arch3"
    ref3.mkdir()
    (ref3 / "sys1.xyz").write_text(
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.4 0.0 0.0\n")
    row3 = mc.compare_build_to_archive(out, ref3 / "sys1.xyz")
    assert row3["verdict"] == "different"

    # Missing: no build output at all.
    row4 = mc.compare_build_to_archive(tmp_path / "nope" / "sys1",
                                       ref3 / "sys1.xyz")
    assert row4["verdict"] == "missing"


# ---------------------------------------------------------------- parity_env
def test_parity_env_runs_fake_cli_and_classifies(tmp_path):
    ids = tmp_path / "ids.txt"
    ids.write_text("sys1|" + FAKE_SMILES + "\n")
    arch = tmp_path / "arch"
    arch.mkdir()
    (arch / "sys1.xyz").write_text(
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.3 0.0 0.0\n"
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.3 0.0 0.0\n")
    work = tmp_path / "work"
    fake = _fake_cmd(tmp_path)
    res = mc.parity_env(ids, arch, work,
                        allow_builds=True, command=fake)
    assert res["summary"] == {"identical": 1, "same_set_other_order": 0,
                              "different": 0, "missing": 0}


def test_parity_env_dry_run_classifies_existing_only(tmp_path):
    ids = tmp_path / "ids.txt"
    ids.write_text("sys1|" + FAKE_SMILES + "\n")
    arch = tmp_path / "arch"
    arch.mkdir()
    (arch / "sys1.xyz").write_text(
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.3 0.0 0.0\n")
    res = mc.parity_env(ids, arch, tmp_path / "work", allow_builds=False)
    assert res["mode"] == "dry_run"
    assert res["summary"]["missing"] == 1


# ---------------------------------------------------------------- parity_gate
def test_parity_gate_dry_run_by_default(tmp_path, monkeypatch):
    ids = tmp_path / "ids.txt"
    ids.write_text("sys1|" + FAKE_SMILES + "\n")
    arch = tmp_path / "arch"
    arch.mkdir()
    (arch / "sys1.xyz").write_text(
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.3 0.0 0.0\n")
    ran = []
    monkeypatch.setattr(mc.subprocess, "run",
                        lambda *a, **k: ran.append(a) or (_ for _ in ()).throw(
                            AssertionError("no subprocess in dry run")))
    res = mc.parity_gate("ref0", ids, arch, tmp_path / "work")
    assert res["mode"] == "dry_run"
    assert not ran
    assert res["summary"]["missing"] == 1
    assert not (tmp_path / "work" / "clone").exists()


def test_parity_gate_executed_uses_fake_install_and_cli(tmp_path, monkeypatch):
    """Full gate path with a fake installer and a fake build CLI."""
    ids = tmp_path / "ids.txt"
    ids.write_text("sys1|" + FAKE_SMILES + "\n")
    arch = tmp_path / "arch"
    arch.mkdir()
    (arch / "sys1.xyz").write_text(
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.3 0.0 0.0\n"
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.3 0.0 0.0\n")
    work = tmp_path / "work"
    fake = _fake_cmd(tmp_path)

    # fake venv: skip real venv creation; pretend the venv python exists
    venv = work / "venv"
    (venv / "bin").mkdir(parents=True)
    (venv / "bin" / "python").write_text("# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n")
    calls = []
    real_run = subprocess.run

    def fake_git(args, cwd, **kw):
        # git goes through safe_git: no real clone, an empty repo dir instead
        s = ["git", *[str(a) for a in args]]
        calls.append(s)
        if args[0] == "clone":
            Path(s[-1]).mkdir(parents=True, exist_ok=True)
        return subprocess.CompletedProcess(s, 0, "", "")

    def fake_run(cmd, **kw):
        calls.append([str(c) for c in cmd])
        return real_run(cmd, **kw)   # the build itself really runs the fake CLI

    from delfin.agent import safe_git
    monkeypatch.setattr(safe_git, "run", fake_git)
    monkeypatch.setattr(mc.subprocess, "run", fake_run)
    monkeypatch.setattr(mc, "_ensure_venv",
                        lambda *a, **k: venv / "bin" / "python")
    res = mc.parity_gate("ref0", ids, arch, work, allow_install=True,
                         command=fake, install_cmd=["true", "install"])
    assert res["mode"] == "executed"
    assert res["summary"]["identical"] == 1
    assert any("clone" in " ".join(c) for c in calls)        # git clone into work/clone
    assert any("checkout" in " ".join(c) for c in calls)     # git checkout of ref
    assert any("fake_cli.py" in " ".join(c) for c in calls)  # builds via injected command


# ---------------------------------------------------------------- load_repro
def test_load_repro_counts_empty_and_distinct(tmp_path):
    fake = _fake_cmd(tmp_path)
    res = mc.load_repro(FAKE_SMILES, tmp_path / "lr", n=3,
                        command=fake)
    assert res["n"] == 3
    assert res["n_empty"] == 0
    assert res["n_distinct"] == 1          # deterministic fake -> one result
    assert res["n_failed"] == 0


def test_load_repro_flags_empty_frames(tmp_path):
    script = tmp_path / "fake_empty.py"
    script.write_text(textwrap.dedent(FAKE_CLI).replace(
        "frames.append(", "frames.append(  # unchanged") )
    # an empty-producing fake: writes no frames
    script = tmp_path / "fake_empty.py"
    script.write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "import json, sys\nfrom pathlib import Path\n"
        "out = Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)\n"
        "(out / 'manifest.json').write_text(json.dumps({'isomers': []}))\n")
    res = mc.load_repro("X", tmp_path / "lr2", n=2,
                        command=[sys.executable, str(script)])
    assert res["n_empty"] == 2


def test_load_repro_rejects_n_over_16(tmp_path):
    with pytest.raises(ValueError):
        mc.load_repro("X", tmp_path / "lr3", n=17)


# ---------------------------------------------------------------- geom_audit
def test_geom_audit_raw_measures_no_table(tmp_path):
    xyz = tmp_path / "f1.xyz"
    xyz.write_text(
        "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "3\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "Pt 0.0 0.0 0.0\nN 2.0 0.0 0.0\nH 0.5 0.5 0.0\n")
    syms, coords = mc.parse_xyz(xyz.read_text())
    rows = mc.geom_audit_frames([("f1", syms, coords)], use_rdkit=False)
    qs = {r["quantity"] for r in rows}
    assert {"bond_len", "angle"} <= qs
    # Pt-N at 2.0 A: normalized by rcov sum (~2.07) -> value < 1 (close contact)
    bn = [r["value"] for r in rows if r["quantity"] == "bond_len"][0]
    assert 0.9 < bn < 1.05
    # angle at N between Pt and H
    ang = [r for r in rows if r["quantity"] == "angle"]
    assert ang and all(0 <= r["value"] <= 180 for r in ang)


def test_geom_audit_runs_without_any_threshold_table(tmp_path):
    xyz = tmp_path / "f2.xyz"
    xyz.write_text(
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.4 0.0 0.0\n")
    syms, coords = mc.parse_xyz(xyz.read_text())
    rows = mc.geom_audit_frames([("f2", syms, coords)], use_rdkit=False)
    assert rows, "raw audit produced no units"
    # and no table file was ever needed:
    assert not list(tmp_path.glob("*table*"))


def test_geom_compare_table_optional_cod_table(tmp_path):
    xyz = tmp_path / "f3.xyz"
    xyz.write_text(
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.70 0.0 0.0\n")  # bonded but stretched: 1.197
    syms, coords = mc.parse_xyz(xyz.read_text())
    rows = mc.geom_audit_frames([("f3", syms, coords)], use_rdkit=False)
    tbl = tmp_path / "cod_table.json"
    tbl.write_text(json.dumps({"bond_len": {"bond_len": [1.0] * 400}}))
    res = mc.geom_compare_table(rows, tbl)
    assert res["n_frames_flagged"] == 1     # 1.197 > p999 = 1.0
    assert res["frames"]["f3"]["bond_len"] == "sev"


# ---------------------------------------------------------------- CLI
def test_cli_geom_audit_json(tmp_path, capsys):
    xyz = tmp_path / "f4.xyz"
    xyz.write_text(
        "2\nSYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
        "C 0.0 0.0 0.0\nO 1.4 0.0 0.0\n")
    rc = mc.main(["geom_audit", str(xyz), "--no-rdkit"])
    out = json.loads(capsys.readouterr().out)
    assert rc == 0
    assert out["n_units"] > 0


def test_cli_load_repro_smoke(tmp_path, capsys, monkeypatch):
    monkeypatch.setattr(mc, "load_repro",
                        lambda *a, **k: {"n": 2, "n_empty": 0,
                                         "n_distinct": 1})
    rc = mc.main(["load_repro", "--smiles", FAKE_SMILES,
                  "--work", str(tmp_path / "lrcli"),
                  "--n", "2"])
    capsys.readouterr()
    assert rc == 0


def test_no_delfin_env_set_in_process(tmp_path):
    """No DELFIN_* construction variable is ever set globally in-process."""
    import os
    before = {k for k in os.environ if k.startswith("DELFIN_")}
    fake = _fake_cmd(tmp_path)
    mc.load_repro(FAKE_SMILES, tmp_path / "lr4", n=1, command=fake)
    assert {k for k in os.environ if k.startswith("DELFIN_")} == before
