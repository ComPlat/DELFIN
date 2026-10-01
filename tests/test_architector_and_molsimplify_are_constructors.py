"""Architector and molSimplify build a SMILES wherever MANTA does.

One definition per tool (``delfin.common.external_builders.build_frames``),
reached from the ARCHITECTOR / MOLSIMPLIFY buttons of the structure editor --
which the Submit Job and ORCA Builder tabs both embed -- and from
``smiles_converter=ARCHITECTOR|MOLSIMPLIFY`` in CONTROL.

The builds themselves run only where the tool is installed: in this
interpreter, or in the one ``DELFIN_ARCHITECTOR_PYTHON`` /
``DELFIN_MOLSIMPLIFY_PYTHON`` names.  Every frame must carry exactly the atoms
the SMILES says -- the old Architector path filled a four-coordinate platinum
up with two waters and failed on every cobalt ammine, because it never told
Architector the coordination number or the oxidation state.
"""

from __future__ import annotations

import collections
import subprocess
import time

import pytest

pytest.importorskip("rdkit")

from delfin.common import external_builders as eb  # noqa: E402

CISPLATIN = "[Pt](Cl)(Cl)([NH3])[NH3]"
COBALT_AMMINE = "[NH3][Co+2]([NH3])([NH3])([NH3])([NH3])Cl"
PT_EN = "Cl[Pt]1(Cl)[NH2]CC[NH2]1"
THREE = (CISPLATIN, COBALT_AMMINE, PT_EN)


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
    """Element count of the atom rows (symbol + three numbers) of *text*."""
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


def _installed(tool):
    python = eb.tool_python(tool)
    module = eb.TOOLS[tool]["module"]
    probe = ("import importlib.util, sys; "
             f"sys.exit(0 if importlib.util.find_spec({module!r}) else 1)")
    try:
        return subprocess.run([python, "-c", probe], timeout=120).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def _needs(tool):
    if not _installed(tool):
        pytest.skip(f"{tool} is not installed in {eb.tool_python(tool)} "
                    f"(set {eb.TOOLS[tool]['python_env']})")


# -- the split, which needs no tool ------------------------------------------

def test_the_split_reads_cn_and_oxidation_state_from_the_smiles():
    spec = eb.split_complex_smiles(CISPLATIN)
    assert (spec["metal"], spec["cn"], spec["metal_ox"]) == ("Pt", 4, 2)
    # An unbracketed Cl bound to the metal is chloride, not HCl.
    assert sorted(l["smiles"] for l in spec["ligands"]) == ["N", "N", "[Cl-]", "[Cl-]"]

    spec = eb.split_complex_smiles(COBALT_AMMINE)
    assert (spec["metal"], spec["cn"], spec["metal_ox"]) == ("Co", 6, 3)

    spec = eb.split_complex_smiles(PT_EN)
    assert (spec["cn"], spec["metal_ox"]) == (4, 2)
    en = [l for l in spec["ligands"] if l["denticity"] == 2]
    assert len(en) == 1 and en[0]["smiles"] == "NCCN" and en[0]["coordList"] == [0, 3]


@pytest.mark.parametrize("smiles, why", [
    ("CCO", "no metal"),
    ("Cl[Cu]Cl.Cl[Cu]Cl", "metal atoms"),
    ("[NH3][Pt+2]([NH3])([NH3])[NH3].[Cl-].[Cl-]", "not bound to the metal"),
])
def test_what_cannot_be_built_is_said(smiles, why):
    with pytest.raises(eb.BuildError, match=why):
        eb.split_complex_smiles(smiles)


def test_a_missing_tool_is_an_error_not_another_builder():
    import importlib.util
    import sys

    checked = 0
    for tool in eb.TOOLS:
        if importlib.util.find_spec(eb.TOOLS[tool]["module"]) is not None:
            continue
        frames, error = eb.build_frames(tool, CISPLATIN, python=sys.executable)
        assert frames == []
        assert "not installed" in error and eb.TOOLS[tool]["python_env"] in error
        if eb.TOOLS[tool].get("install"):
            assert eb.TOOLS[tool]["install"] in error
        checked += 1
    if not checked:
        pytest.skip("both tools are installed in this interpreter")


# -- the builds, where the tool is -----------------------------------------------

# epic-MACE has its own set (tests/test_mace_is_a_constructor.py): it does not
# embed [Co(NH3)5Cl]2+ within its ten attempts, and it builds hapto ligands.
@pytest.mark.parametrize("tool", sorted(set(eb.TOOLS) - {"mace"}))
@pytest.mark.parametrize("smiles", THREE)
def test_every_frame_has_the_formula_of_the_smiles(tool, smiles):
    """Through the dashboard's own build function, headless."""
    _needs(tool)
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import structure_editor as se

    result = se._run_external_build(smiles, tool)
    assert result["error"] is None, result["error"]
    frames = result["isomers"]
    assert frames, "no frame came back"
    want = _formula_of_smiles(smiles)
    for body, n_atoms, label in frames:
        assert _formula_of_xyz(body) == want, (label, _formula_of_xyz(body), want)
        assert n_atoms == sum(want.values())
        assert label.startswith(eb.TOOLS[tool]["display"])


@pytest.mark.parametrize("tool, convert", [
    ("architector", "smiles_to_xyz_architector"),
    ("molsimplify", "smiles_to_xyz_molsimplify"),
    ("mace", "smiles_to_xyz_mace"),
])
def test_the_control_converter_uses_the_same_build(tool, convert):
    _needs(tool)
    from delfin import smiles_converter

    xyz, error = getattr(smiles_converter, convert)(CISPLATIN)
    assert error is None, error
    assert _formula_of_xyz(xyz) == _formula_of_smiles(CISPLATIN)


# -- the buttons, in both tabs ---------------------------------------------------

def _a_dashboard_context(tmp_path):
    from delfin.dashboard.context import DashboardContext

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc", archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda _script: None
    return ctx


def _wait_until_the_build_is_over(state, budget=1800):
    began = time.time()
    time.sleep(0.2)
    while state.get("smiles_busy") and time.time() - began < budget:
        time.sleep(0.1)
    time.sleep(0.3)


def _press_in_tab(tab, tool, tmp_path):
    """Type the SMILES where the tab takes one and press the tool's button."""
    from delfin.dashboard import tab_orca_builder, tab_submit

    ctx = _a_dashboard_context(tmp_path)
    if tab == "submit":
        _widget, refs = tab_submit.create_tab(ctx)
        refs["coords_widget"].value = CISPLATIN
    else:
        _widget, refs = tab_orca_builder.create_tab(ctx)
        refs["orca_coords"].value = CISPLATIN
    button = refs.get(f"{tool}_button")
    assert button is not None, f"the {tab} tab hands out no {tool}_button"
    button.click()
    state = refs["editor_state"]
    _wait_until_the_build_is_over(state)
    return state


def _what_landed(tab, state):
    """The structures the tab now keeps: the isomer set, or the named blocks."""
    if tab == "submit":
        return [xyz for xyz, _n, _label in state.get("isomers") or []]
    return [xyz for _name, xyz in state.get("xyz_blocks") or []]


_TWO_FRAMES = [
    ("Pt 0 0 0\nCl 2.3 0 0\nCl 0 2.3 0\nN -2.0 0 0\nN 0 -2.0 0", 5, "Fake a"),
    ("Pt 0 0 0\nCl 2.3 0 0\nCl -2.3 0 0\nN 0 2.0 0\nN 0 -2.0 0", 5, "Fake b"),
]


@pytest.mark.parametrize("tab", ["submit", "orca_builder"])
@pytest.mark.parametrize("tool", sorted(eb.TOOLS))
def test_the_button_puts_every_frame_where_the_tab_keeps_structures(
        tab, tool, tmp_path, monkeypatch):
    """Wiring, with the tool stubbed: both tabs, both buttons, all frames."""
    pytest.importorskip("ipywidgets")
    asked = []

    def stub(tool_name, smiles, **_kw):
        asked.append((tool_name, smiles))
        return list(_TWO_FRAMES), None

    monkeypatch.setattr(eb, "build_frames", stub)
    state = _press_in_tab(tab, tool, tmp_path)
    assert asked == [(tool, CISPLATIN)]
    assert len(_what_landed(tab, state)) == len(_TWO_FRAMES)


@pytest.mark.parametrize("tab", ["submit", "orca_builder"])
@pytest.mark.parametrize("tool", sorted(eb.TOOLS))
def test_a_failed_build_lands_nothing(tab, tool, tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    monkeypatch.setattr(eb, "build_frames",
                        lambda *_a, **_k: ([], "the tool is not installed"))
    state = _press_in_tab(tab, tool, tmp_path)
    assert not state.get("isomers")
    assert not state.get("smiles_busy")


@pytest.mark.parametrize("tab", ["submit", "orca_builder"])
@pytest.mark.parametrize("tool", sorted(eb.TOOLS))
def test_the_real_button_builds_in_both_tabs(tab, tool, tmp_path):
    _needs(tool)
    pytest.importorskip("ipywidgets")
    state = _press_in_tab(tab, tool, tmp_path)
    landed = _what_landed(tab, state)
    assert landed, "no frame landed"
    want = _formula_of_smiles(CISPLATIN)
    for xyz in landed:
        assert _formula_of_xyz(xyz) == want
