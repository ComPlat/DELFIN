"""rxembed, Molassembler and MetalloGen build a SMILES wherever MANTA and the other constructors do.

Each runs in an environment of its own that DELFIN starts through
``delfin.common.external_builders`` and never imports (rxembed needs Python >= 3.12 and
RDKit >= 2026.9.1; Molassembler's environment has no RDKit; MetalloGen is run with its ORCA
calculator only, never Gaussian).  Their builds are defined once in
``delfin.cluster_bench.tool_workers.smiles_tools``; the RXEMBED / MOLASSEMBLER / METALLOGEN
buttons of the structure editor (Submit Job and ORCA Builder), ``smiles_converter=RXEMBED|
MOLASSEMBLER|METALLOGEN`` in CONTROL and ``delfin cluster --tool rxembed|molassembler|metallogen``
all run it.  ``python -m delfin.installer --install rxembed`` (``molassembler``, ``metallogen``)
builds the environment; ``DELFIN_<TOOL>_PYTHON`` points elsewhere.

The real builds run only where the tool is: the interpreter that variable names, or the
environment the installer built (``DELFIN_AI_TOOLS_ROOT``); MetalloGen also needs ORCA.
"""

from __future__ import annotations

import collections
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from delfin.common import external_builders as eb  # noqa: E402
from delfin.cluster_bench.tool_workers import smiles_tools as st  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
NEW = ("rxembed", "molassembler", "metallogen")

CISPLATIN = "[Pt](Cl)(Cl)([NH3])[NH3]"
CISPLATIN_DATIVE = "N->[Pt](<-N)(Cl)Cl"
CISPLATIN_CHARGED = "[NH3+][Pt-2]([NH3+])(Cl)Cl"
ZNCL4 = "Cl[Zn-2](Cl)(Cl)Cl"
PT_EN = "Cl[Pt]1(Cl)[NH2]CC[NH2]1"

_S = 1 / math.sqrt(3)
SQUARE = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]]
TETRAHEDRON = [[_S, _S, _S], [_S, -_S, -_S], [-_S, _S, -_S], [-_S, -_S, _S]]
OCTAHEDRON = [[0, 0, 1], [0, 0, -1], [0, 1, 0], [0, -1, 0], [1, 0, 0], [-1, 0, 0]]
VECTORS = {"4_square_planar": SQUARE, "4_tetrahedral": TETRAHEDRON, "6_octahedral": OCTAHEDRON}
KNOWN = {4: ["4_tetrahedral", "4_square_planar"], 6: ["6_octahedral"]}


def _formula_of_smiles(smiles):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    count = collections.Counter()
    for atom in mol.GetAtoms():
        count[atom.GetSymbol()] += 1
        count["H"] += atom.GetTotalNumHs()
    return count


def _formula_of_xyz(text):
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


def _python_with(tool):
    """An interpreter in which *tool* imports, or skip."""
    python = eb.tool_python(tool)
    probe = f"import {eb.TOOLS[tool]['module']}"
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env["PYTHONNOUSERSITE"] = "1"
    try:
        ok = subprocess.run([python, "-c", probe], timeout=300, env=env,
                            capture_output=True).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        ok = False
    if not ok:
        pytest.skip(f"{eb.TOOLS[tool]['display']} is not installed in {python} (set "
                    f"{eb.TOOLS[tool]['python_env']} or run {eb.TOOLS[tool]['install']})")
    return python


def _needs_orca():
    from delfin.qm_runtime import find_tool_executable

    if not find_tool_executable("orca"):
        pytest.skip("MetalloGen runs with its ORCA calculator, and no ORCA was found")


# -- the registry: optional, own environment, install hint -----------------------------

@pytest.mark.parametrize("tool", NEW)
def test_each_tool_has_an_environment_of_its_own_and_an_installer(tool, tmp_path):
    info = eb.TOOLS[tool]
    assert info["install"] == f"python -m delfin.installer --install {tool}"
    root = tmp_path / "ai_tools"
    env = {"DELFIN_AI_TOOLS_ROOT": str(root)}
    assert eb.tool_python(tool, env) == sys.executable
    python = root / ".mamba_env" / tool / "bin" / "python"
    python.parent.mkdir(parents=True)
    python.write_text("#!/bin/sh\n")
    python.chmod(0o755)
    assert eb.tool_python(tool, env) == str(python)
    env[info["python_env"]] = "/somewhere/else/python"
    assert eb.tool_python(tool, env) == "/somewhere/else/python"


@pytest.mark.parametrize("tool", NEW)
def test_a_missing_tool_names_the_installer_and_builds_nothing_else(tool, tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(tmp_path))
    monkeypatch.delenv(eb.TOOLS[tool]["python_env"], raising=False)
    frames, error = eb.build_frames(tool, CISPLATIN, python=sys.executable)
    assert frames == []
    assert eb.TOOLS[tool]["install"] in error and eb.TOOLS[tool]["python_env"] in error


def test_metallogen_without_orca_is_an_error_not_gaussian(tmp_path, monkeypatch):
    import delfin.qm_runtime as qr

    python = tmp_path / "python"
    python.write_text("#!/bin/sh\nexit 1\n")
    python.chmod(0o755)
    monkeypatch.setattr(qr, "find_tool_executable", lambda _name: None)
    frames, error = eb.build_frames("metallogen", CISPLATIN, python=str(python))
    assert frames == [] and "ORCA" in error and "Gaussian" in error


def test_importing_delfin_needs_none_of_them():
    code = ("import sys, delfin.common.external_builders, delfin.ai_tools, "
            "delfin.cluster_bench.tool_workers.smiles_tools; "
            "bad = [m for m in ('rxembed', 'scine_molassembler', 'MetalloGen') if m in sys.modules]; "
            "sys.exit(1 if bad else 0)")
    env = dict(os.environ, PYTHONPATH=str(REPO))
    assert subprocess.run([sys.executable, "-c", code], env=env, timeout=300).returncode == 0


def test_the_ai_registry_lists_them_with_their_install_hint(monkeypatch, tmp_path):
    from delfin import ai_tools

    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(tmp_path))
    for var in ("DELFIN_RXEMBED_PYTHON", "DELFIN_MOLASSEMBLER_PYTHON", "DELFIN_METALLOGEN_PYTHON"):
        monkeypatch.delenv(var, raising=False)
    rows = {r["name"]: r for r in ai_tools.collect_ai_summary()["tools"]}
    for name, tool in (("rxembed", "rxembed"), ("Molassembler", "molassembler"),
                       ("MetalloGen", "metallogen")):
        assert rows[name]["install_hint"] == f"python -m delfin.installer --install {tool}"
        assert rows[name]["installed"] is False and rows[name]["version"] == ""


def test_the_installer_knows_them_and_looks_into_their_environment(tmp_path, monkeypatch):
    from delfin import installer

    monkeypatch.setenv("DELFIN_AI_TOOLS_ROOT", str(tmp_path))
    for tool, module, switch in (("rxembed", "rxembed", "INSTALL_RXEMBED"),
                                 ("molassembler", "scine_molassembler", "INSTALL_MOLASSEMBLER"),
                                 ("metallogen", "MetalloGen", "INSTALL_METALLOGEN")):
        monkeypatch.delenv(eb.TOOLS[tool]["python_env"], raising=False)
        found = installer.find(tool)
        assert found.group == "ai" and found.own_env == tool and found.switch == switch
        assert tool in installer.profile("all") and module not in installer.packages()
        assert not installer.present(found)
        prefix = tmp_path / ".mamba_env" / tool
        (prefix / "bin").mkdir(parents=True)
        (prefix / "bin" / "python").write_text("#!/bin/sh\n")
        (prefix / "bin" / "python").chmod(0o755)
        assert not installer.present(found)
        site = prefix / "lib" / "python3.12" / "site-packages" / module
        site.mkdir(parents=True)
        (site / "__init__.py").write_text("")
        assert installer.present(found)
    assert installer.find("scine-molassembler") is installer.find("molassembler")


def test_the_ai_installer_builds_each_environment_with_its_pin(tmp_path):
    """With a fake micromamba: each environment goes where DELFIN looks for it."""
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
        'case "\\$*" in *"import "*) [ -e "$prefix/installed" ] || exit 1 ;; esac\n'
        'exit 0\n'
        'EOS\n'
        'chmod 755 "$prefix/bin/python"\n')
    (fake / "micromamba").chmod(0o755)
    # A stand-in download: the archive's bytes, recorded URL; the pinned digest is replaced
    # by the digest of these bytes for the archive that is to pass, and kept for the one
    # that is to be refused.
    (fake / "curl").write_text(
        '#!/bin/sh\n'
        'echo "curl $*" >> "$LOGDIR/calls.log"\n'
        'while [ $# -gt 0 ]; do [ "$1" = "-o" ] && dest="$2"; shift; done\n'
        'printf "fake-archive" > "$dest"\n')
    (fake / "curl").chmod(0o755)
    import hashlib

    good = hashlib.sha256(b"fake-archive").hexdigest()
    root = tmp_path / "ai_tools"
    root.mkdir()
    done = subprocess.run(
        ["bash", str(REPO / "delfin" / "ai_tools" / "install_ai_tools.sh")],
        capture_output=True, text=True, timeout=120,
        env={"PATH": os.pathsep.join([str(fake), "/usr/bin", "/bin"]), "HOME": str(tmp_path),
             "LOGDIR": str(tmp_path), "DELFIN_PYTHON": sys.executable,
             "DELFIN_AI_TOOLS_ROOT": str(root), "INSTALL_RXEMBED": "1",
             "INSTALL_MOLASSEMBLER": "1", "INSTALL_METALLOGEN": "1",
             "RXEMBED_SHA256": good, "MOLASSEMBLER_SHA256": good})
    said = done.stdout + done.stderr
    calls = (tmp_path / "calls.log").read_text()
    assert "curl -fsSL -o" in calls
    assert "aligfellow/rxembed/archive/9c0ba5a6048bdbcba94450df5d6a68067bd250c8.tar.gz" in calls
    assert ("pip install --no-cache-dir numpy==2.2.2 pillow==11.1.0 rdkit==2026.9.1 "
            "networkx==3.7 xyzgraph==1.7.0") in calls
    assert f"pip install --no-cache-dir --no-deps {root}/.mamba_env/rxembed/rxembed-" in calls
    assert ("pip install --no-cache-dir scine-utilities==10.1.0 numpy==2.4.6 scipy==1.17.1"
            in calls)
    # Molassembler: the pinned CPython 3.11 wheel, fetched and checked, installed under its own file name
    # (pip reads the Python and platform tags from it).
    assert "/scine_molassembler-3.0.1-cp311-cp311-manylinux_2_28_x86_64.whl" in calls
    assert (f"pip install --no-cache-dir --no-deps {root}/.mamba_env/molassembler/"
            "scine_molassembler-3.0.1-cp311-cp311-manylinux_2_28_x86_64.whl") in calls
    assert (f"create -y -p {root}/.mamba_env/rxembed -c conda-forge --override-channels "
            "python=3.12.16 pip") in calls
    assert ("python=3.11.18 rdkit=2026.09.2 numpy=2.4.6 scipy=1.17.1 cclib=1.8.1 pulp=3.3.2"
            in calls)
    for tool, label in (("rxembed", "rxembed"), ("molassembler", "Molassembler")):
        python = root / ".mamba_env" / tool / "bin" / "python"
        assert f"{label} installed: {python}" in said, said[-2000:]
        assert eb.managed_python(tool, {"DELFIN_AI_TOOLS_ROOT": str(root)}) == str(python)
    # The stand-in is not the pinned MetalloGen archive: refused, never installed.
    assert "is not the pinned 1a8a4f6a0223b7aa138b2eb0f57af4be6a84e056ed5c972bcbac2d6fdf3ab250" in said
    assert "MetalloGen installed" not in said
    assert "--no-deps https://github.com/kyunghoonlee777" not in calls


# -- what each tool can be given, without the tool ---------------------------------------

def test_molassembler_cannot_read_a_dative_bond_and_says_so():
    with pytest.raises(st.NotExpressible, match="dative"):
        st.run_molassembler(CISPLATIN_DATIVE, "/nonexistent")


def test_the_rotations_of_each_polyhedron():
    assert len(st.site_rotations(OCTAHEDRON)) == 24
    assert len(st.site_rotations(TETRAHEDRON)) == 12
    assert len(st.site_rotations(SQUARE)) == 8        # C4 about the normal and four C2 in plane


def _arrangements(smiles, geom):
    arrs, _sites, _ligs, truncated = st.metallogen_arrangements(
        eb.split_complex_smiles(smiles), VECTORS[geom])
    assert not truncated
    return len(arrs)


def test_trien_dichloro_has_five_arrangements():
    """[Co(trien)Cl2]+: cis-alpha, cis-beta (each Delta and Lambda) and trans.  The two
    terminal and the two inner N of trien are exchanged by the ligand's own symmetry; naming
    donors atom by atom counted 4."""
    assert _arrangements("Cl[Co+]123(Cl)[NH2]CC[NH]1CC[NH]2CC[NH2]3", "6_octahedral") == 5


def test_a_cap_on_the_arrangements_is_reported():
    arrs, _s, _l, truncated = st.metallogen_arrangements(
        eb.split_complex_smiles("Cl[Co+]123(Cl)[NH2]CC[NH]1CC[NH]2CC[NH2]3"), OCTAHEDRON, cap=2)
    assert len(arrs) == 2 and truncated
    report = {}
    st.metallogen_jobs(eb.split_complex_smiles(CISPLATIN), {"max_arrangements": 1}, KNOWN,
                       VECTORS, report)
    assert report == {"truncated": True, "truncated_geometries": ["4_square_planar"]}


@pytest.mark.parametrize("bad", ["/opt/my orca/orca", "/opt/orca;rm/orca", "/opt/$HOME/orca"])
def test_metallogen_refuses_a_path_the_shell_would_split(bad, tmp_path):
    with pytest.raises(st.NotExpressible, match="shell"):
        st._mg_check_path("ORCA path", bad)
    st._mg_check_path("working directory", str(tmp_path))      # an ordinary path passes


def test_frames_without_an_energy_go_last(monkeypatch, tmp_path):
    """The first frame is the lowest-energy one; a frame without an energy never comes first."""
    import shutil

    real = tmp_path / "result_source.json"
    real.write_text(json.dumps({"frames": [
        {"label": "a", "symbols": ["Zn"], "coords": [[0, 0, 0]], "energy": None},
        {"label": "b", "symbols": ["Zn"], "coords": [[0, 0, 1]], "energy": -2.0},
        {"label": "c", "symbols": ["Zn"], "coords": [[0, 0, 2]], "energy": -3.0}],
        "error": None}))
    fake = tmp_path / "python"
    fake.write_text(f'#!/bin/sh\ncp {real} "$5"\n')
    fake.chmod(0o755)
    frames, error = eb.build_frames("metallogen", ZNCL4, python=str(fake),
                                    options={"orca": shutil.which("sh")})
    assert error is None
    assert [label.split()[1] for _b, _n, label in frames] == ["c", "b", "a"]


def test_each_tool_has_its_own_time_limit_and_metallogen_keeps_what_it_finished(
        tmp_path, monkeypatch):
    import shutil

    monkeypatch.delenv("DELFIN_EXTERNAL_BUILDER_TIMEOUT", raising=False)
    assert eb.TOOLS["metallogen"]["timeout"] == 7200
    monkeypatch.setenv("DELFIN_EXTERNAL_BUILDER_KEEP", "0")
    fake = tmp_path / "python"
    # Writes one finished arrangement to the partial xyz it is told about, then hangs.
    # ($4 is the request; the partial xyz the request names sits beside it.)
    fake.write_text('#!/bin/sh\n'
                    'grep -q partial_xyz "$4" || exit 3\n'
                    'p="$(dirname "$4")/partial.xyz"\n'
                    'printf "1\\nx\\nZn 0 0 0\\n" > "$p"\nsleep 30\n')
    fake.chmod(0o755)
    frames, error = eb.build_frames("metallogen", ZNCL4, python=str(fake), timeout=2,
                                    workdir=str(tmp_path), options={"orca": shutil.which("sh")})
    assert frames == [] and "did not finish within 2 s" in error
    kept = error.split("kept in ", 1)[1].strip()
    assert Path(kept).read_text().startswith("1\nx\nZn")


def test_an_isolated_tool_gets_no_user_site_on_the_cluster():
    from delfin.cluster_bench.provenance import cbatch_tool_child_env, cbatch_tool_isolated

    for tool in NEW + ("mace",):
        assert cbatch_tool_isolated(tool)
    assert not cbatch_tool_isolated("architector") and not cbatch_tool_isolated("manta")
    assert cbatch_tool_child_env({}, nousersite=True)["PYTHONNOUSERSITE"] == "1"


def test_metallogen_seeds_each_m_smiles():
    assert st.METALLOGEN_DEFAULTS["seed"] == 42
    assert st.MOLASSEMBLER_MAX_LEAVES == 50000


@pytest.mark.parametrize("smiles, geom, n", [
    (CISPLATIN, "4_square_planar", 2),                       # cis, trans
    (CISPLATIN, "4_tetrahedral", 1),
    (PT_EN, "4_square_planar", 1),                           # en cannot span trans
    (ZNCL4, "4_tetrahedral", 1),
    ("[NH3][Co+2]([NH3])([NH3])([NH3])([NH3])Cl", "6_octahedral", 1),
    ("[NH3][Co+](Cl)([NH3])([NH3])([NH3])Cl", "6_octahedral", 2),    # cis, trans
    ("[NH3][Co](Cl)(Cl)([NH3])([NH3])Cl", "6_octahedral", 2),        # fac, mer
    ("[Co+3]123([NH2]CC[NH2]1)([NH2]CC[NH2]2)[NH2]CC[NH2]3", "6_octahedral", 2),  # Delta, Lambda
])
def test_metallogen_gets_every_distinct_arrangement_once(smiles, geom, n):
    assert _arrangements(smiles, geom) == n


def test_bis_en_dichloro_is_cis_delta_cis_lambda_and_trans():
    smiles = "Cl[Co+]1(Cl)([NH2]CC[NH2]1)1[NH2]CC[NH2]1"
    assert _arrangements(smiles, "6_octahedral") == 3


def test_the_m_smiles_carries_metal_charge_sites_and_geometry():
    jobs = st.metallogen_jobs(eb.split_complex_smiles(CISPLATIN), None, KNOWN, VECTORS)
    assert [(g, k) for g, k, _m in jobs] == [("4_square_planar", 0), ("4_square_planar", 1),
                                              ("4_tetrahedral", 0)]
    for _g, _k, msmiles in jobs:
        parts = msmiles.split("|")
        assert parts[0] == "[Pt+2]" and parts[-1] in VECTORS
        sites = sorted(int(x.split(":")[1].rstrip("]")) for x in parts[1:-1])
        assert sites == [1, 2, 3, 4]
        assert sorted(p.split(":")[0] for p in parts[1:-1]) == ["[Cl-", "[Cl-", "[NH3", "[NH3"]
    with pytest.raises(st.NotExpressible, match="no geometry for 5"):
        st.metallogen_jobs(eb.split_complex_smiles("Cl[Fe](Cl)(Cl)(Cl)Cl"), None, KNOWN, VECTORS)


# -- one build, the same settings, from the dashboard and from CONTROL --------------------

def _recording_stub(calls, tool_label):
    def stub(tool, smiles, **kwargs):
        calls.append((tool, smiles, kwargs))
        return [("Pt 0 0 0\nCl 2.3 0 0", 2, f"{tool_label} x")], None
    return stub


@pytest.mark.parametrize("tool", NEW)
def test_dashboard_and_control_call_the_same_build(tool, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin import smiles_converter
    from delfin.dashboard import structure_editor as se

    calls = []
    monkeypatch.setattr(eb, "build_frames", _recording_stub(calls, tool))
    se._run_external_build(CISPLATIN, tool)
    xyz, error = getattr(smiles_converter, f"smiles_to_xyz_{tool}")(CISPLATIN)
    assert error is None and xyz.startswith("Pt")
    assert [c[:2] for c in calls] == [(tool, CISPLATIN), (tool, CISPLATIN)]
    assert calls[0][2] == calls[1][2] == {}


@pytest.mark.parametrize("tool", NEW)
def test_control_accepts_the_converter_and_the_pipeline_builds_with_it(tool, tmp_path, monkeypatch):
    from delfin.common.control_validator import _as_smiles_converter
    from delfin.define import TEMPLATE
    from delfin.workflows import pipeline

    mode = tool.upper()
    assert _as_smiles_converter(tool) == mode
    assert mode in TEMPLATE.split("smiles_converter=[", 1)[1].split("]", 1)[0].split("|")
    assert pipeline._resolve_smiles_converter({"smiles_converter": mode}) == mode
    calls = []
    monkeypatch.setattr(eb, "build_frames", _recording_stub(calls, tool))
    (tmp_path / "input.txt").write_text(CISPLATIN + "\n")
    control = tmp_path / "CONTROL.txt"
    control.write_text(f"smiles_converter={mode}\n")
    pipeline.normalize_input_file({"smiles_converter": mode, "input_file": "input.txt"}, control)
    assert [c[:2] for c in calls] == [(tool, CISPLATIN)]
    assert (tmp_path / "start.txt").read_text().split()[0] == "Pt"


def test_the_old_template_placeholder_still_reads_as_normal():
    from delfin.common.control_validator import _as_smiles_converter

    assert _as_smiles_converter("[QUICK|NORMAL|MANTA|ARCHITECTOR|MOLSIMPLIFY|MACE]") == "NORMAL"
    assert _as_smiles_converter(
        "[QUICK|NORMAL|MANTA|ARCHITECTOR|MOLSIMPLIFY|MACE|RXEMBED|MOLASSEMBLER|METALLOGEN]") == "NORMAL"


# -- the buttons, in both tabs ----------------------------------------------------------

def _widgets_under(widget):
    yield widget
    for child in getattr(widget, "children", ()) or ():
        yield from _widgets_under(child)


@pytest.mark.parametrize("tab", ["submit", "orca_builder"])
def test_the_buttons_stand_beside_the_other_constructors(tab, tmp_path):
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
    shown = list(_widgets_under(widget))
    for tool in NEW:
        button = refs[f"{tool}_button"]
        assert any(w is button for w in shown), f"{tool}_button is not in the {tab} tab"
        assert button.description == tool.upper()
        assert eb.TOOLS[tool]["install"] in button.tooltip


def test_the_construction_batch_offers_them():
    pytest.importorskip("ipywidgets")
    from delfin.cluster_bench.provenance import MODES, TOOLS
    from delfin.dashboard.construction_batch import CB_MODE_HELP, CB_TOOL_CHOICES

    offered = dict(CB_TOOL_CHOICES)
    for tool in NEW:
        assert offered[tool.upper()] == tool
        assert TOOLS[tool]["worker"] == "smiles_worker.py" and TOOLS[tool]["needs_specs"]
        for mode in MODES[tool]:
            assert (tool, mode) in CB_MODE_HELP


# -- the real builds, where the tool is ----------------------------------------------------

@pytest.mark.parametrize("smiles", [CISPLATIN_DATIVE, ZNCL4])
def test_rxembed_builds_every_frame_with_the_formula_of_the_smiles(smiles, monkeypatch):
    monkeypatch.setenv("DELFIN_RXEMBED_PYTHON", _python_with("rxembed"))
    frames, error = eb.build_frames("rxembed", smiles)
    assert error is None, error
    want = _formula_of_smiles(smiles)
    for body, n_atoms, label in frames:
        assert _formula_of_xyz(body) == want and n_atoms == sum(want.values())
        assert label.startswith("rxembed ")
    if smiles == CISPLATIN_DATIVE:
        species = {label.split()[1].rsplit("-conf", 1)[0] for _b, _n, label in frames}
        assert any("cis" in s for s in species) and any("trans" in s for s in species)


@pytest.mark.parametrize("smiles", [CISPLATIN_CHARGED, ZNCL4])
def test_molassembler_builds_every_frame_with_the_formula_of_the_smiles(smiles, monkeypatch):
    monkeypatch.setenv("DELFIN_MOLASSEMBLER_PYTHON", _python_with("molassembler"))
    frames, error = eb.build_frames("molassembler", smiles)
    assert error is None, error
    assert len(frames) == 10           # one stereoisomer, ten conformers
    want = _formula_of_smiles(smiles)
    for body, n_atoms, label in frames:
        assert _formula_of_xyz(body) == want and n_atoms == sum(want.values())
        assert label.startswith("Molassembler iso0-conf")
    again, _ = eb.build_frames("molassembler", smiles)
    assert again == frames, "Molassembler is seeded: the same build twice is the same"


def test_molassembler_counts_a_meso_form_once(monkeypatch):
    """FC(Cl)C(F)Cl: (R,R), (S,S) and meso -- three, as the benchmark runner counts."""
    monkeypatch.setenv("DELFIN_MOLASSEMBLER_PYTHON", _python_with("molassembler"))
    frames, error = eb.build_frames("molassembler", "FC(Cl)C(F)Cl")
    assert error is None, error
    assert {label.split()[1].split("-")[0] for _b, _n, label in frames} == {"iso0", "iso1", "iso2"}


def test_molassembler_through_its_environment_refuses_a_dative_smiles(monkeypatch):
    monkeypatch.setenv("DELFIN_MOLASSEMBLER_PYTHON", _python_with("molassembler"))
    frames, error = eb.build_frames("molassembler", CISPLATIN_DATIVE)
    assert frames == [] and "dative" in error


def test_metallogen_builds_cis_and_trans_platin_with_orca(monkeypatch):
    monkeypatch.setenv("DELFIN_METALLOGEN_PYTHON", _python_with("metallogen"))
    _needs_orca()
    frames, error = eb.build_frames("metallogen", CISPLATIN)
    assert error is None, error
    want = _formula_of_smiles(CISPLATIN)
    built = set()
    for body, n_atoms, label in frames:
        assert _formula_of_xyz(body) == want and n_atoms == sum(want.values())
        assert "GFN2-xTB via ORCA" in label
        built.add(label.split()[1].rsplit("-conf", 1)[0])
    assert {"4_square_planar-arr0", "4_square_planar-arr1"} <= built


@pytest.mark.parametrize("tool, smiles", [("rxembed", ZNCL4), ("molassembler", ZNCL4)])
def test_the_control_converter_takes_the_first_frame(tool, smiles, monkeypatch):
    monkeypatch.setenv(eb.TOOLS[tool]["python_env"], _python_with(tool))
    from delfin import smiles_converter

    xyz, error = getattr(smiles_converter, f"smiles_to_xyz_{tool}")(smiles)
    assert error is None, error
    assert _formula_of_xyz(xyz) == _formula_of_smiles(smiles)


@pytest.mark.parametrize("tool, smiles", [("rxembed", ZNCL4), ("molassembler", ZNCL4)])
def test_the_cluster_worker_writes_the_archive(tool, smiles, tmp_path):
    python = _python_with(tool)
    from delfin.cluster_bench.prepare import cbatch_spec_of

    specs = tmp_path / "specs.jsonl"
    specs.write_text(json.dumps(cbatch_spec_of(("ZNCL4", smiles))) + "\n")
    archive = tmp_path / "archive"
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env["PYTHONNOUSERSITE"] = "1"
    subprocess.run([python, str(REPO / "delfin" / "cluster_bench" / "tool_workers" / "smiles_worker.py"),
                    tool, "ZNCL4", str(specs), str(archive), str(tmp_path / "work")],
                   env=env, check=True, timeout=600)
    meta = json.loads((archive / "_meta" / "ZNCL4.json").read_text())
    assert meta["status"] == "ok" and meta["n_frames"] > 0
    text = (archive / "ZNCL4.xyz").read_text()
    assert text.splitlines()[1].startswith("ZNCL4 frame0 ") and f"tool={tool}" in text
