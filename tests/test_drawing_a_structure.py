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

    editor["submit_draw_open_btn"].value = False
    assert editor["submit_draw_frame"].value == "", (
        "a closed editor must not go on running in a hidden frame"
    )


def test_the_same_drawing_twice_reads_as_two_answers(editor, tmp_path):
    """A widget only reports a value that changed, so the same structure drawn
    again would look like an answer that never came."""
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    handler = source.split("def on_submit_draw_get")[1].split("\n    def ")[0]
    assert "Date.now()" in handler
    reader = source.split("def on_submit_draw_sync")[1].split("\n    def ")[0]
    assert "raw.split('\\n', 1)[1]" in reader


def test_the_editor_is_asked_for_a_molfile_not_for_its_own_smiles(editor):
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
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

    source = open(tab_submit.__file__, encoding="utf-8").read()
    opening = source.split("def on_submit_draw_open")[1].split("\n    def ")[0]
    assert "_ketcher.install(" not in opening, "opening it must not fetch it"
    assert "latest_release()" in opening
    assert "Press Fetch it" in opening

    fetching = source.split("def on_submit_draw_update")[1].split("\n    def ")[0]
    assert "_ketcher.install(" in fetching
