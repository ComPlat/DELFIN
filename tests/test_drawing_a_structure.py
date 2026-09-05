"""Draw the structure instead of typing it.

The Submit tab takes a SMILES or a block of coordinates.  Typing a SMILES is
fine for ethanol and hopeless for anything a chemist would actually draw, so
Ketcher is offered beside the box: draw it, press TO SMILES, and it lands in
the field the rest of the tab already reads.
"""

from __future__ import annotations

import os
import pathlib

import pytest

from delfin.dashboard import ketcher
from delfin.dashboard.context import DashboardContext
from editor_source import SUBMIT_SOURCE


def molblock(smiles):
    rdkit = pytest.importorskip("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    return Chem.MolToMolBlock(mol)


# ---------------------------------------------------------------------------
# what was drawn, read the way the rest of DELFIN reads structures
# ---------------------------------------------------------------------------
def test_a_drawing_comes_back_as_a_smiles_rdkit_wrote():
    """The editor can write a SMILES itself, and it would be a different
    dialect.  Everything downstream here reads structures with RDKit, so the
    SMILES is made by RDKit from the drawing: one RDKit wrote is one RDKit will
    certainly read back."""
    for smiles in ("CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"):
        outcome = ketcher.smiles_from_molfile(molblock(smiles))
        assert outcome["ok"] is True, outcome["status"]
        assert outcome["smiles"] == smiles


def test_a_molfile_keeps_its_three_header_lines():
    """The first of them is usually blank.  Stripping the leading newline
    shifts every line up by one, the counts line is read as a title, and a
    perfectly good drawing comes back unreadable -- which is what happened."""
    body = molblock("CCO")
    assert body.startswith("\n"), "this molfile has the blank first line"

    assert ketcher.smiles_from_molfile(body)["ok"] is True

    source = open(ketcher.__file__, encoding="utf-8").read()
    reader = source.split("def smiles_from_molfile")[1]
    assert "text = str(molfile or '')" in reader
    assert "text.strip()" in reader, "emptiness is still checked"


def test_something_a_chemist_draws_that_no_valence_table_likes():
    """Drawn structures are often chemically loose -- an open valence, a metal
    with a bond count no table agrees with.  Refusing them would refuse the
    reason for having a drawing editor at all."""
    outcome = ketcher.smiles_from_molfile(molblock("C[Pt](Cl)(Cl)N"))

    assert outcome["ok"] is True, outcome["status"]
    assert "Pt" in outcome["smiles"]

    lone = ketcher.smiles_from_molfile(molblock("[Fe]"))
    assert lone["ok"] is True and "Fe" in lone["smiles"]


def test_nothing_drawn_is_said_plainly_not_raised():
    assert ketcher.smiles_from_molfile("")["ok"] is False
    assert "Nothing was drawn" in ketcher.smiles_from_molfile("")["status"]
    assert ketcher.smiles_from_molfile("not a molfile")["ok"] is False


# ---------------------------------------------------------------------------
# where it lives, and how the browser reaches it
# ---------------------------------------------------------------------------
def test_it_is_put_where_the_browser_can_reach_it(tmp_path, monkeypatch):
    """An absolute path on the machine the kernel runs on means nothing to a
    browser that may be on the other side of an SSH tunnel.  Voila serves
    everything below its root at /voila/files/, which is the same route the
    literature tab already uses for PDFs."""
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    # This is the route for a dashboard nobody launched with delfin-voila:
    # with one of its own the editor is loaded from where it is kept and
    # nothing is put under the root at all.
    monkeypatch.delenv(ketcher._URL_ENV, raising=False)
    # A home of its own: the editor is kept there now, and a copy on the
    # machine running the tests would otherwise be served into this root.
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(pathlib.Path, "home", lambda: home)

    folder = ketcher.app_directory()
    assert folder == tmp_path / ".delfin" / "ketcher"
    assert ketcher.installed_version() is None
    assert ketcher.app_url() is None, "no URL for something that is not there"

    folder.mkdir(parents=True)
    (folder / "index.html").write_text("<html></html>")
    (folder / ".delfin-ketcher-version").write_text("3.17.0")

    assert ketcher.installed_version() == "3.17.0"
    url = ketcher.app_url()
    assert url.startswith("/voila/files/.delfin/ketcher/index.html")
    assert "v=3.17.0" in url, (
        "an update that kept the same URL would be shown from the cache"
    )


def test_without_a_served_directory_it_says_so(monkeypatch):
    monkeypatch.delenv("DELFIN_VOILA_ROOT_DIR", raising=False)
    monkeypatch.delenv(ketcher._URL_ENV, raising=False)

    assert ketcher.app_directory() is None
    assert ketcher.app_url() is None
    outcome = ketcher.install()
    assert outcome["ok"] is False
    assert "not serving a directory" in outcome["status"]


def test_only_what_the_main_page_needs_is_kept():
    """Four bundles ship, one per demo page, each 29.5 MB of the same chemistry
    engine compiled in.  Keeping the one the main page loads is the difference
    between 115 MB on disk and 32."""
    page = ('<script defer="defer" src="./static/js/main.abc.js"></script>'
            '<link href="./static/css/main.def.css" rel="stylesheet">')
    names = [
        "standalone/",
        "standalone/index.html",
        "standalone/popup.html",
        "standalone/favicon.ico",
        "standalone/static/js/main.abc.js",
        "standalone/static/js/popup.xyz.js",
        "standalone/static/js/157.chunk.js",
        "standalone/static/css/main.def.css",
        "__MACOSX/standalone/._index.html",
    ]

    keep = ketcher._wanted(names, page)

    assert "standalone/index.html" in keep
    assert "standalone/static/js/main.abc.js" in keep
    assert "standalone/static/js/157.chunk.js" in keep
    assert "standalone/static/css/main.def.css" in keep
    assert "standalone/static/js/popup.xyz.js" not in keep, "29.5 MB unused"
    assert "standalone/popup.html" not in keep, "its bundle is not kept"
    assert "standalone/" not in keep, "a folder entry opened as a file is EISDIR"
    assert not any("__MACOSX" in name for name in keep)


def test_the_newest_build_is_asked_for_rather_than_pinned():
    """Keeping up to date is a button, not an edit to this file."""
    source = open(ketcher.__file__, encoding="utf-8").read()
    assert "releases/latest" in ketcher.RELEASES
    assert ketcher._ASSET.match("ketcher-standalone-3.17.0.zip")
    assert not ketcher._ASSET.match("ketcher-standalone.zip"), (
        "two files of the same content, and only one says which build it is"
    )
    assert "3.17" not in source.split('"""')[2], "no version is pinned in code"


# ---------------------------------------------------------------------------
# in the tab
# ---------------------------------------------------------------------------
@pytest.fixture
def editor(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = lambda _script: None
    _widget, refs = tab_submit.create_tab(ctx)
    return refs


def test_the_drawing_lands_in_the_box_the_rest_of_the_tab_reads(editor, tmp_path):
    """Convert, Build and everything downstream see it exactly as a typed
    SMILES, because it is in the same field."""
    folder = tmp_path / ".delfin" / "ketcher"
    folder.mkdir(parents=True)
    (folder / "index.html").write_text("<html></html>")
    (folder / ".delfin-ketcher-version").write_text("3.17.0")

    editor["submit_draw_open_btn"].value = True
    assert editor["submit_draw_frame"].layout.display == ""
    assert editor["submit_draw_get_btn"].layout.display == ""
    assert "/voila/files/.delfin/ketcher/index.html" in editor["submit_draw_frame"].value

    editor["submit_draw_sync"].value = "1\n" + molblock("CCO")
    assert editor["coords_widget"].value == "CCO"

    drawn_frame = editor["submit_draw_frame"].value
    editor["submit_draw_open_btn"].value = False
    assert editor["submit_draw_frame"].layout.display == "none"
    assert editor["submit_draw_frame"].value == drawn_frame, (
        "folding it away must not end the editor: what was drawn has to still "
        "be there when it is folded open again"
    )

    editor["submit_draw_open_btn"].value = True
    assert editor["submit_draw_frame"].value == drawn_frame, (
        "and opening it again must not reload it, which would empty it"
    )


def test_the_same_drawing_twice_reads_as_two_answers(editor, tmp_path):
    """A widget only reports a value that changed, so the same structure drawn
    again would look like an answer that never came."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    handler = source.split("def on_submit_draw_get")[1].split("\n    def ")[0]
    assert "Date.now()" in handler
    reader = source.split("def on_submit_draw_sync")[1].split("\n    def ")[0]
    assert "raw.split('\\n', 1)[1]" in reader


def test_the_editor_is_asked_for_a_molfile_not_for_its_own_smiles(editor):
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    handler = source.split("def on_submit_draw_get")[1].split("\n    def ")[0]
    assert "getMolfile()" in handler
    assert "getSmiles" not in handler
    assert "contentWindow.ketcher" in handler, "same origin, so it can be asked"


def test_an_editor_that_is_not_open_says_so_rather_than_nothing(editor):
    editor["submit_draw_sync"].value = "1\n!no-editor"
    assert "not open yet" in editor["mol_status"].value


def test_it_is_offered_rather_than_fetched(editor):
    """Thirty-odd megabytes over the network, and on a machine without one it
    is a wait that ends in nothing."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    opening = source.split("def on_submit_draw_open")[1].split("\n    def ")[0]
    assert "_ketcher.install(" not in opening, "opening it must not fetch it"
    assert "latest_release()" in opening
    assert "Press Fetch it" in opening

    fetching = source.split("def on_submit_draw_update")[1].split("\n    def ")[0]
    assert "_ketcher.install(" in fetching


def test_the_editor_stands_between_drawing_and_converting(tmp_path, monkeypatch):
    """The order on screen is the order of the work: draw it, hand it over as a
    SMILES, then turn that into coordinates.

    With the editor open, the buttons that act on what it produced are below
    it, where the eye arrives after the drawing rather than before it.
    """
    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets
    from delfin.dashboard import tab_submit

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(tmp_path))
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = lambda _script: None
    tab, refs = tab_submit.create_tab(ctx)

    def column(node):
        kids = getattr(node, "children", ())
        if any(kid is refs["coords_widget"] for kid in kids):
            return node
        for kid in kids:
            found = column(kid)
            if found is not None:
                return found
        return None

    order = []
    for child in column(tab).children:
        if child is refs["coords_widget"]:
            order.append("input")
        elif child is refs["submit_draw_frame"]:
            order.append("editor")
        elif isinstance(child, widgets.HBox):
            labels = [getattr(k, "description", "") for k in child.children]
            if "DRAW" in labels:
                order.append("draw")
            elif "CONVERT SMILES" in labels:
                order.append("convert")

    assert order.index("input") < order.index("draw") < order.index("editor")
    assert order.index("editor") < order.index("convert"), (
        "the convert buttons belong below the editor, not above it"
    )


# ---------------------------------------------------------------------------
# a complex has to come out in the form the rest of DELFIN reads
# ---------------------------------------------------------------------------
def _molfile(atoms, bonds):
    """A molfile written by hand, the way a drawing arrives.

    The hydrogen field counts from one: 1 means no hydrogen, 4 means three.
    """
    head = f'{len(atoms):3d}{len(bonds):3d}  0  0  0  0            999 V2000'
    rows = [f'    0.0000    0.0000    0.0000 {sym:<3s} 0  0  0  0  0  {h}'
            for sym, h in atoms]
    links = [f'{a:3d}{b:3d}{kind:3d}  0' for a, b, kind in bonds]
    return '\n'.join(['', '  Ketcher', '', head] + rows + links + ['M  END', ''])


def test_a_coordination_bond_comes_out_charge_separated():
    """RDKit writes a coordination bond as an arrow, [n]->[Cd].  The structures
    DELFIN is built around -- MANTA's inputs, the batch files, every internal
    conversion -- write it as a plain bond with the charge on both ends: the
    donor [N+] and the metal carrying one minus per bond it accepts.

    Both describe the same molecule and only the second is understood
    downstream, so a drawing is brought into it.
    """
    ammine = _molfile([("N", 4), ("Pt", 1), ("Cl", 1), ("Cl", 1)],
                      [(1, 2, 1), (2, 3, 1), (2, 4, 1)])

    outcome = ketcher.smiles_from_molfile(ammine)

    assert outcome["ok"] is True, outcome["status"]
    assert outcome["smiles"] == "[NH3+][Pt-]([Cl])[Cl]"
    assert "->" not in outcome["smiles"], "the arrow form is not what is read"
    assert outcome["dative"] == 1
    assert "charge on both ends" in outcome["status"]


def test_the_hydrogens_of_a_donor_survive_the_conversion():
    """[NH3]->[Pt] came back as [N+][Pt-]: the count was frozen before it had
    been written down, so the ammonia lost its hydrogens."""
    drawn = _molfile([("N", 0), ("Pt", 1)], [(1, 2, 9)])

    assert ketcher.smiles_from_molfile(drawn)["smiles"] == "[NH3+][Pt-]"


def test_what_was_drawn_is_what_comes_out():
    """An amine drawn with two hydrogens is an amide, not an ammine, and it
    keeps its plain bond: the nitrogen has room for it.  Guessing otherwise
    would be correcting the chemist rather than reading the drawing."""
    amide = _molfile([("N", 3), ("Pt", 1), ("Cl", 1), ("Cl", 1)],
                     [(1, 2, 1), (2, 3, 1), (2, 4, 1)])

    outcome = ketcher.smiles_from_molfile(amide)

    assert outcome["smiles"] == "[NH2][Pt]([Cl])[Cl]"
    assert outcome["dative"] == 0
    # and the halides stay plain covalent either way, as they are in the
    # batch files: [Cl][Cd-3] -- the metal's charge counts the dative bonds
    assert "[Cl-]" not in outcome["smiles"]


def test_a_complex_drawn_naively_matches_the_batch_entry_for_it():
    """The whole question in one test: draw a complex the natural way, without
    charges, and get the molecule DELFIN's own data says it is.

    ACUDOT in the batch file is
    CC1=[N+]2NC(C3=CC=CC=N3)=[O+][Cd-3]2([Br])([Br])[N+]2=CC=CC=C12 -- and
    drawing it with plain bonds gives the same molecule, 33 atoms either way.
    """
    from rdkit import Chem

    naive = "CC1=N2NC(C3=CC=CC=N3)=O[Cd]2([Br])([Br])N2=CC=CC=C12"
    entry = ("CC1=[N+]2NC(C3=CC=CC=N3)=[O+][Cd-3]2([Br])([Br])"
             "[N+]2=CC=CC=C12")

    drawn = ketcher.smiles_from_molfile(molblock(naive))

    assert drawn["ok"] is True, drawn["status"]
    assert drawn["dative"] == 3, "three donors, three minus on the metal"
    assert drawn["smiles"] == Chem.MolToSmiles(Chem.MolFromSmiles(entry)), (
        "a drawn complex must be the molecule the batch file says it is"
    )


def test_delfin_builds_coordinates_from_what_was_drawn():
    """The point of the whole path: what comes out of the editor goes into the
    converters the rest of DELFIN uses, and they make a structure of it."""
    pytest.importorskip("rdkit")
    from delfin.dashboard.input_processing import smiles_to_xyz_quick

    naive = "CC1=N2NC(C3=CC=CC=N3)=O[Cd]2([Br])([Br])N2=CC=CC=C12"
    drawn = ketcher.smiles_from_molfile(molblock(naive))

    _xyz, atoms, _method, error = smiles_to_xyz_quick(drawn["smiles"])

    assert error is None, error
    assert atoms == 33, f"the batch entry is 33 atoms, this gave {atoms}"


def test_a_molecule_without_a_metal_is_left_exactly_as_drawn():
    outcome = ketcher.smiles_from_molfile(molblock("CC(=O)Oc1ccccc1C(=O)O"))

    assert outcome["dative"] == 0
    assert outcome["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"


def test_the_viewers_own_draw_still_offers_an_element(editor):
    """Two functions of the same name in one scope: the later replaces the
    earlier.  The drawing panel's refresh was called _refresh_draw_controls,
    which is what the viewer's own Draw uses to show the element dropdown --
    so switching Draw on stopped offering an element to draw with.
    """
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    assert source.count("    def _refresh_draw_controls():") == 1, (
        "one name, one function"
    )
    assert "    def _refresh_ketcher_controls():" in source

    editor["coords_widget"].value = "3\nw\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n"
    assert editor["submit_element_dd"].layout.display == "none"

    editor["submit_draw_btn"].value = True
    assert editor["submit_element_dd"].layout.display == "", (
        "Draw has to offer an element to draw with"
    )

    editor["submit_draw_btn"].value = False
    assert editor["submit_element_dd"].layout.display == "none"


def test_a_converted_structure_can_be_worked_on_at_once(editor):
    """Showing a structure is what enables the editing toolbar -- there is no
    second step between converting and being able to take hold of it."""
    from delfin.dashboard import tab_submit

    source = SUBMIT_SOURCE
    shower = source.split("def _replace_mol_output_view")[1].split("\n    def ")[0]
    assert "_set_manip_toolbar_enabled(True)" in shower

    editor["coords_widget"].value = "3\nw\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n"
    assert editor["submit_ff_dd"].disabled is False
    assert editor["submit_draw_btn"].disabled is False


def test_the_editor_is_kept_somewhere_that_does_not_move(tmp_path, monkeypatch):
    """It was installed beside the launch directory, and so it was not
    installed at all.

    Voila serves whatever the dashboard was started in, and the editor has to
    live under that for a browser to load it -- but it was also *kept* there.
    Start the dashboard from somewhere else and a Ketcher that had been
    fetched looked absent, so the answer was to fetch thirty megabytes again,
    into a place that would lose it just as easily. One of those places is the
    cache directory under /tmp, which the system sweeps.
    """
    from delfin.dashboard import ketcher

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(pathlib.Path, "home", lambda: home)
    monkeypatch.delenv(ketcher._URL_ENV, raising=False)

    kept = ketcher.stored_directory()
    kept.mkdir(parents=True)
    (kept / "index.html").write_text("<html>ketcher</html>")
    (kept / ketcher._STAMP).write_text("2.31.0")

    for name in ("first", "second"):
        root = tmp_path / name
        root.mkdir()
        monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
        assert ketcher.installed_version() == "2.31.0", name
        assert (root / ".delfin" / "ketcher" / "index.html").is_file()


def test_an_editor_that_is_already_there_is_taken_in(tmp_path, monkeypatch):
    """Somebody who fetched it before there was a place to keep it must not
    have to fetch it again."""
    from delfin.dashboard import ketcher

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(pathlib.Path, "home", lambda: home)

    monkeypatch.delenv(ketcher._URL_ENV, raising=False)

    old_root = tmp_path / "where it was"
    served = old_root / ".delfin" / "ketcher"
    served.mkdir(parents=True)
    (served / "index.html").write_text("<html>ketcher</html>")
    (served / ketcher._STAMP).write_text("2.30.0")

    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(old_root))
    assert ketcher.installed_version() == "2.30.0"
    assert (ketcher.stored_directory() / ketcher._STAMP).is_file(), (
        "it was left where the next launch directory would lose it"
    )

    new_root = tmp_path / "where it is started next time"
    new_root.mkdir()
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(new_root))
    assert ketcher.installed_version() == "2.30.0", "it asked to fetch again"


def test_the_download_keeps_a_copy_that_outlives_the_launch_directory():
    import inspect

    from delfin.dashboard import ketcher

    body = inspect.getsource(ketcher.install)
    assert "stored_directory()" in body
    # Kept first, served from a copy of it -- not the other way round.
    assert body.index("shutil.move(str(unpacked), str(kept))") < body.index(
        "shutil.copytree(kept, folder)")
