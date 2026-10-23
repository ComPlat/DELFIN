"""epic-MACE builds a SMILES wherever MANTA, Architector and molSimplify do.

epic-MACE (Chernyshov & Pidko, JCTC 2024; GPL-3.0) is an external program in
an environment of its own -- Python 3.7, RDKit 2020.09 -- that DELFIN starts
through ``delfin.common.external_builders`` and never imports.  The same build
(``build_frames('mace', smiles)`` with ``MACE_DEFAULTS``) is behind the MACE
button of the structure editor (Submit Job and ORCA Builder) and behind
``smiles_converter=MACE`` in CONTROL; ``python -m delfin.installer --install
epic-mace`` builds the environment, ``DELFIN_MACE_PYTHON`` points elsewhere.

The real builds run only where epic-MACE is: the interpreter
``DELFIN_MACE_PYTHON`` names, or the environment the installer built.
"""

from __future__ import annotations

import collections
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import child_env

pytest.importorskip("rdkit")

from delfin.common import external_builders as eb  # noqa: E402

REPO = Path(__file__).resolve().parents[1]

CISPLATIN = "[Pt](Cl)(Cl)([NH3])[NH3]"
PT_EN = "Cl[Pt]1(Cl)[NH2]CC[NH2]1"
CO_TRIAMMINE = "Cl[Co](Cl)(Cl)(N)(N)N"
# Zeise's anion: an eta2 ethylene, written with one bond from each carbon.
ZEISE = "[CH2]1=[CH2][Pt-]1(Cl)(Cl)Cl"
# Ferrocene: two eta5 cyclopentadienides, every ring carbon bound to iron.
FERROCENE = ("[Fe+2]12345%10%11%12%13%14.[CH-]1%20[CH]2=[CH]3[CH]4=[CH]5%20."
             "[CH-]%10%21[CH]%11=[CH]%12[CH]%13=[CH]%14%21")
# (eta6-benzene)Cr(CO)3, carbonyls dative.
ARENE_CR = ("[Cr]123456(<-[C-]#[O+])(<-[C-]#[O+])<-[C-]#[O+]."
            "[CH]1%20=[CH]2[CH]3=[CH]4[CH]5=[CH]6%20")


def _mace_formula_of_smiles(smiles):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    count = collections.Counter()
    for atom in mol.GetAtoms():
        count[atom.GetSymbol()] += 1
        count["H"] += atom.GetTotalNumHs()
    return count


def _mace_formula_of_xyz(text):
    count = collections.Counter()
    for row in text.splitlines():
        parts = row.split()
        if len(parts) == 4:
            try:
                [float(v) for v in parts[1:]]
            except ValueError:
                continue
            count[parts[0]] += 1
    return count


def _mace_python():
    """An interpreter with epic-MACE in it, or skip."""
    python = eb.tool_python("mace")
    probe = ("import sys; import mace; "
             "sys.exit(0 if hasattr(mace, 'ComplexFromLigands') else 1)")
    try:
        ok = subprocess.run([python, "-c", probe], timeout=120,
                            env=dict(os.environ, PYTHONNOUSERSITE="1"),
                            capture_output=True).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        ok = False
    if not ok:
        pytest.skip(f"epic-MACE is not installed in {python} (set DELFIN_MACE_PYTHON or "
                    f"run python -m delfin.installer --install epic-mace)")
    return python


# -- the input epic-MACE gets, which needs only RDKit -------------------------------

def test_monodentate_donors_are_mapped_and_the_metal_carries_its_oxidation_state():
    job = eb.mace_job(eb.split_complex_smiles(CISPLATIN), "extended")
    assert job["CA"] == "[Pt+2]"
    assert sorted(job["ligands"]) == ["[Cl-:1]", "[Cl-:1]", "[NH3:1]", "[NH3:1]"]
    assert job["geoms"] == ["SP", "TET"]
    assert eb.mace_job(eb.split_complex_smiles(CISPLATIN), "paper")["geoms"] == ["SP"]


def test_a_chelate_keeps_both_donors_in_one_ligand():
    job = eb.mace_job(eb.split_complex_smiles(PT_EN))
    en = [s for s in job["ligands"] if "Cl" not in s]
    assert len(en) == 1 and en[0].count(":1]") == 2


def test_a_hapto_ring_is_one_centroid_bound_to_one_ring_carbon():
    job = eb.mace_job(eb.split_complex_smiles(FERROCENE))
    assert job["CA"] == "[Fe+2]"
    assert job["geoms"] == ["SAN"]
    assert job["info"]["hapto"] == [(5, "anchor"), (5, "anchor")]
    for smi in job["ligands"]:
        assert smi.count("[*:1]") == 1

    job = eb.mace_job(eb.split_complex_smiles(ARENE_CR))
    assert job["CA"] == "[Cr]"
    assert (6, "anchor") in job["info"]["hapto"]
    assert job["info"]["n_sites"] == 4 and job["geoms"] == ["SP", "TET"]


def test_a_hapto_alkene_is_a_centroid_bound_to_both_carbons():
    from rdkit import Chem

    job = eb.mace_job(eb.split_complex_smiles(ZEISE))
    assert job["info"]["hapto"] == [(2, "star")]
    alkene = [s for s in job["ligands"] if "*" in s][0]
    mol = Chem.MolFromSmiles(alkene)
    star = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0][0]
    assert sorted(n.GetSymbol() for n in star.GetNeighbors()) == ["C", "C"]


def test_a_site_count_without_a_geometry_is_said():
    five = "Cl[Fe](Cl)(Cl)(Cl)Cl"
    with pytest.raises(eb.BuildError, match="no geometry for 5 donor sites"):
        eb.mace_job(eb.split_complex_smiles(five), "paper")
    assert eb.mace_job(eb.split_complex_smiles(five), "extended")["geoms"] == ["SPY", "TBP"]


# -- where epic-MACE is looked for -------------------------------------------------

def test_the_interpreter_is_the_variable_then_the_managed_environment_then_none(tmp_path):
    root = tmp_path / "ai_tools"
    env = {"DELFIN_AI_TOOLS_ROOT": str(root)}
    assert eb.tool_python("mace", env) == sys.executable

    python = root / ".mamba_env" / "epic_mace" / "bin" / "python"
    python.parent.mkdir(parents=True)
    python.write_text("#!/bin/sh\n")
    python.chmod(0o755)
    assert eb.managed_python("mace", env) == str(python)
    assert eb.tool_python("mace", env) == str(python)

    env["DELFIN_MACE_PYTHON"] = "/somewhere/else/python"
    assert eb.tool_python("mace", env) == "/somewhere/else/python"
    # Architector and molSimplify have no environment of their own.
    assert eb.managed_python("architector", env) is None


def test_a_missing_mace_names_the_installer(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(tmp_path))
    monkeypatch.delenv("DELFIN_MACE_PYTHON", raising=False)
    frames, error = eb.build_frames("mace", CISPLATIN, python=sys.executable)
    assert frames == []
    assert "python -m delfin.installer --install epic-mace" in error
    assert "DELFIN_MACE_PYTHON" in error


def test_the_mace_potential_is_not_taken_for_epic_mace(tmp_path):
    """``import mace`` is also mace-torch; the worker must refuse it."""
    fake = tmp_path / "site" / "mace"
    fake.mkdir(parents=True)
    (fake / "__init__.py").write_text("# the machine-learning potential\n")
    request = tmp_path / "request.json"
    result = tmp_path / "result.json"
    request.write_text(json.dumps({"smiles": CISPLATIN, "work": str(tmp_path), "options": {}}))
    env = {**child_env(tmp_path), "PYTHONPATH": os.pathsep.join([str(REPO), str(tmp_path / "site")])}
    subprocess.run([sys.executable, "-m", "delfin.common.external_builders", "mace",
                    str(request), str(result)], env=env, check=True, timeout=120)
    answer = json.loads(result.read_text())
    assert answer["frames"] == [] and answer.get("not_installed")


# -- one build, the same settings, from the dashboard and from CONTROL --------------

def _mace_recording_stub(calls):
    def stub(tool, smiles, **kwargs):
        calls.append((tool, smiles, kwargs))
        return [("Pt 0 0 0\nCl 2.3 0 0", 2, "epic-MACE SP-iso0-conf0")], None
    return stub


def test_dashboard_and_control_call_the_same_build_with_the_same_settings(monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin import smiles_converter
    from delfin.dashboard import structure_editor as se

    calls = []
    monkeypatch.setattr(eb, "build_frames", _mace_recording_stub(calls))
    se._run_external_build(CISPLATIN, "mace")
    xyz, error = smiles_converter.smiles_to_xyz_mace(CISPLATIN)
    assert error is None and xyz.startswith("Pt")
    assert [c[:2] for c in calls] == [("mace", CISPLATIN), ("mace", CISPLATIN)]
    # Neither side passes options of its own: both get MACE_DEFAULTS.
    assert calls[0][2] == calls[1][2] == {}
    assert eb.MACE_DEFAULTS == {"geometries": "extended", "num_confs": 10, "max_attempts": 10}


def test_control_accepts_mace_and_the_pipeline_builds_with_it(tmp_path, monkeypatch):
    from delfin.common.control_validator import _as_smiles_converter
    from delfin.define import TEMPLATE
    from delfin.workflows import pipeline

    assert _as_smiles_converter("mace") == "MACE"
    assert "smiles_converter=[QUICK|NORMAL|MANTA|ARCHITECTOR|MOLSIMPLIFY|MACE]" in TEMPLATE
    assert pipeline._resolve_smiles_converter({"smiles_converter": "MACE"}) == "MACE"

    calls = []
    monkeypatch.setattr(eb, "build_frames", _mace_recording_stub(calls))
    (tmp_path / "input.txt").write_text(CISPLATIN + "\n")
    control = tmp_path / "CONTROL.txt"
    control.write_text("smiles_converter=MACE\n")
    pipeline.normalize_input_file({"smiles_converter": "MACE", "input_file": "input.txt"}, control)
    assert [c[:2] for c in calls] == [("mace", CISPLATIN)]
    assert (tmp_path / "start.txt").read_text().split()[0] == "Pt"


# -- the button, in both tabs ---------------------------------------------------------

def _mace_widgets_under(widget):
    yield widget
    for child in getattr(widget, "children", ()) or ():
        yield from _mace_widgets_under(child)


@pytest.mark.parametrize("tab", ["submit", "orca_builder"])
def test_the_mace_button_is_shown_beside_the_other_constructors(tab, tmp_path):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_orca_builder, tab_submit
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc", archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda _script: None
    module = tab_submit if tab == "submit" else tab_orca_builder
    widget, refs = module.create_tab(ctx)
    shown = list(_mace_widgets_under(widget))
    for name in ("manta_button", "architector_button", "molsimplify_button", "mace_button"):
        assert any(w is refs[name] for w in shown), f"{name} is not in the {tab} tab"
    assert refs["mace_button"].description == "MACE"
    assert "epic-MACE" in refs["mace_button"].tooltip


# -- the real builds, where epic-MACE is ------------------------------------------------

@pytest.mark.parametrize("smiles, first", [
    (CISPLATIN, "SP"), (PT_EN, "SP"), (CO_TRIAMMINE, "OH"),
    (ZEISE, "SP"), (FERROCENE, "SAN"), (ARENE_CR, "SP"),
])
def test_every_frame_has_the_formula_of_the_smiles(smiles, first, monkeypatch):
    monkeypatch.setenv("DELFIN_MACE_PYTHON", _mace_python())
    frames, error = eb.build_frames("mace", smiles)
    assert error is None, error
    assert frames
    want = _mace_formula_of_smiles(smiles)
    for body, n_atoms, label in frames:
        # The hapto centroids (X) are not atoms and are not handed on.
        assert _mace_formula_of_xyz(body) == want, (label, _mace_formula_of_xyz(body), want)
        assert n_atoms == sum(want.values())
        assert label.startswith("epic-MACE ")
    assert frames[0][2].split()[1].startswith(first + "-")


def test_cis_and_trans_platin_are_both_built(monkeypatch):
    monkeypatch.setenv("DELFIN_MACE_PYTHON", _mace_python())
    frames, error = eb.build_frames("mace", CISPLATIN)
    assert error is None, error
    square = {label.split()[1].rsplit("-conf", 1)[0] for _b, _n, label in frames
              if label.split()[1].startswith("SP-")}
    assert square == {"SP-iso0", "SP-iso1"}


def test_the_control_converter_takes_the_first_frame(monkeypatch):
    monkeypatch.setenv("DELFIN_MACE_PYTHON", _mace_python())
    from delfin import smiles_converter

    xyz, error = smiles_converter.smiles_to_xyz_mace(PT_EN)
    assert error is None, error
    assert _mace_formula_of_xyz(xyz) == _mace_formula_of_smiles(PT_EN)


# -- installing it the usual DELFIN way ------------------------------------------------

def test_the_installer_knows_epic_mace_apart_from_the_mace_potential():
    from delfin import installer

    tool = installer.find("epic-mace")
    assert tool is installer.find("epic_mace") and tool.group == "ai"
    assert tool.switch == "INSTALL_EPIC_MACE" and tool.own_env == "mace"
    assert installer.find("mace").group == "mlp"
    assert installer.packages().get("mace", {}).get("tool", "mace") == "mace"
    assert "epic_mace" in installer.profile("all")


def test_present_looks_into_the_environment_of_its_own(tmp_path, monkeypatch):
    from delfin import installer

    monkeypatch.delenv("DELFIN_MACE_PYTHON", raising=False)
    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(tmp_path))
    tool = installer.find("epic-mace")
    assert not installer.present(tool)
    prefix = tmp_path / ".mamba_env" / "epic_mace"
    (prefix / "bin").mkdir(parents=True)
    (prefix / "bin" / "python").write_text("#!/bin/sh\n")
    (prefix / "bin" / "python").chmod(0o755)
    assert not installer.present(tool), "an environment without the package is not epic-MACE"
    site = prefix / "lib" / "python3.7" / "site-packages" / "mace"
    site.mkdir(parents=True)
    (site / "__init__.py").write_text("")
    assert installer.present(tool)


def test_the_ai_installer_builds_a_python37_environment_and_installs_the_commit(tmp_path):
    """With a fake micromamba: the environment goes where DELFIN looks for it."""
    fake = tmp_path / "fakebin"
    fake.mkdir()
    (fake / "micromamba").write_text(
        '#!/bin/sh\n'
        'echo "micromamba $*" >> "$LOGDIR/calls.log"\n'
        'prefix=""; while [ $# -gt 0 ]; do [ "$1" = "-p" ] && prefix="$2"; shift; done\n'
        'mkdir -p "$prefix/bin"\n'
        'cat > "$prefix/bin/python" <<EOS\n'
        '#!/bin/sh\n'
        'echo "python \\$*" >> "$LOGDIR/calls.log"\n'
        'case "\\$*" in *"pip install"*) touch "$prefix/installed" ;; esac\n'
        'case "\\$*" in *"import mace"*) [ -e "$prefix/installed" ] || exit 1 ;; esac\n'
        'exit 0\n'
        'EOS\n'
        'chmod 755 "$prefix/bin/python"\n')
    (fake / "micromamba").chmod(0o755)
    root = tmp_path / "ai_tools"
    root.mkdir()
    script = REPO / "delfin" / "ai_tools" / "install_ai_tools.sh"

    done = subprocess.run(
        ["bash", str(script)], capture_output=True, text=True, timeout=120,
        env={"PATH": os.pathsep.join([str(fake), "/usr/bin", "/bin"]), "HOME": str(tmp_path),
             "LOGDIR": str(tmp_path), "DELFIN_PYTHON": sys.executable,
             "DELFIN_AI_TOOLS_ROOT": str(root), "INSTALL_EPIC_MACE": "1"})

    said = done.stdout + done.stderr
    calls = (tmp_path / "calls.log").read_text()
    env_dir = root / ".mamba_env" / "epic_mace"
    assert (f"create -y -p {env_dir} -c conda-forge --override-channels "
            "python=3.7 rdkit=2020.09.5") in calls
    assert ("pip install --no-deps --no-cache-dir --force-reinstall https://github.com/EPiCs-group/"
            "epic-mace/archive/efb5778e715ea461f80cf3bbc752929101bd0bb3.tar.gz") in calls
    assert f"epic-MACE installed: {env_dir}/bin/python" in said, said[-2000:]
    assert eb.managed_python("mace", {"DELFIN_AI_TOOLS_ROOT": str(root)}) == str(env_dir / "bin" / "python")


def test_the_settings_install_button_runs_delfins_installer(tmp_path, monkeypatch):
    """epic-MACE cannot be pip-installed beside DELFIN; its row asks the installer."""
    pytest.importorskip("ipywidgets")
    import threading
    import time

    from delfin import ai_tools, installer
    from delfin.dashboard import tab_settings
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setenv("HOME", str(tmp_path))
    real_run = subprocess.run

    def no_pip(command, *args, **kwargs):          # no network for "pip list --outdated"
        if list(command[1:4]) == ["-m", "pip", "list"]:
            return subprocess.CompletedProcess(command, 1, "", "")
        return real_run(command, *args, **kwargs)

    monkeypatch.setattr(subprocess, "run", no_pip)
    monkeypatch.setattr(ai_tools, "collect_ai_summary", lambda: {"tools": [{
        "name": "epic-MACE", "category": "Metal Complex ML", "installed": False,
        "version": "", "description": "", "install_hint":
        "python -m delfin.installer --install epic-mace"}]})
    asked = []
    done = threading.Event()

    def fake_install(names, **_kw):
        asked.append(list(names))
        done.set()
        return {"ok": True, "results": []}

    monkeypatch.setattr(installer, "install", fake_install)
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc", archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda _script: None
    widget, _refs = tab_settings.create_tab(ctx)
    refresh = [w for w in _mace_widgets_under(widget)
               if getattr(w, "description", "") == "Refresh AI status"]
    assert refresh, "no Refresh AI status button"
    refresh[0].click()
    install = [w for w in _mace_widgets_under(widget)
               if getattr(w, "description", "") == "Install"
               and getattr(w, "tooltip", "") == "python -m delfin.installer --install epic-mace"]
    assert install, "the epic-MACE row has no installer button"
    install[0].click()
    assert done.wait(30)
    time.sleep(0.2)
    assert asked == [["epic-mace"]]
