"""Backspace goes up a folder; Delete deletes the selection, never the folder you are in.

Found in the bug hunt of 11.11. (job 8543260): in the Calculations browser
(and Archive, Home, Office, which are the same builder) Backspace and Delete
were both bound to "delete". With nothing selected -- after Esc, which the
list binds to "deselect" -- Backspace inside ``gamma/`` asked "Delete current
folder? gamma", and Yes ran ``shutil.rmtree`` on it.

Convention (Windows Explorer, Finder, VS Code explorer, JupyterLab):

* Backspace and Alt+Up go one folder up;
* Delete (Cmd+Backspace on a Mac) acts on the *selected* items, after a
  confirmation that names them, and on nothing else;
* no key press deletes the folder that is open.  Deleting it stays an
  explicit click on the Delete button.

Every test works in a ``tmp_path`` folder; nothing outside it is touched.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard import tab_remote_archive as remote
from delfin.dashboard.context import DashboardContext

REPO = Path(__file__).resolve().parents[1]
_NODE = shutil.which("node")


# -- the key map, as the browser runs it -------------------------------------------

def _key_map_js(source: str) -> str:
    """The lines of the list's keydown handler that pick the action."""
    match = re.search(
        r"var action = '';\n(.*?)\n\s*if \(!?action\)", source, re.S)
    assert match, "keydown action map not found"
    return match.group(1)


def _actions(key_map: str) -> dict:
    cases = {
        "Backspace": "{key:'Backspace'}",
        "Delete": "{key:'Delete'}",
        "Alt+ArrowUp": "{key:'ArrowUp', altKey:true}",
        "Meta+Backspace": "{key:'Backspace', metaKey:true}",
        "Ctrl+Backspace": "{key:'Backspace', ctrlKey:true}",
        "Escape": "{key:'Escape'}",
        "Enter": "{key:'Enter'}",
    }
    script = "var out = {};\n"
    for name, event in cases.items():
        script += (
            "(function(){ var e = Object.assign({ctrlKey:false, metaKey:false,"
            f" altKey:false, shiftKey:false}}, {event}); var action = '';\n"
            f"{key_map}\n out[{name!r}] = action; }})();\n")
    script += "console.log(JSON.stringify(out));\n"
    proc = subprocess.run([_NODE, "-e", script], capture_output=True,
                          text=True, timeout=30)
    assert proc.returncode == 0, proc.stderr
    import json
    return json.loads(proc.stdout)


@pytest.mark.parametrize("path", [
    REPO / "delfin" / "dashboard" / "tab_calculations_browser.py",
    REPO / "delfin" / "dashboard" / "tab_remote_archive.py",
])
def test_backspace_is_up_and_only_delete_deletes(path):
    key_map = _key_map_js(path.read_text(encoding="utf-8"))
    # Static floor, also where node is missing: Backspace alone is never delete.
    assert "e.key === 'Backspace') action = 'delete'" not in key_map
    if _NODE is None:
        pytest.skip("node is not installed")
    got = _actions(key_map)
    assert got["Backspace"] == "up"
    assert got["Alt+ArrowUp"] == "up"
    assert got["Delete"] == "delete"
    assert got["Meta+Backspace"] == "delete"
    assert got["Ctrl+Backspace"] == "delete"
    assert got["Escape"] == "deselect"
    assert got["Enter"] == "open"


# -- the Calculations browser (also Archive, Home, Office) -----------------------

@pytest.fixture
def calc(tmp_path):
    root = tmp_path / "calc"
    for name in ("alpha", "beta", "gamma"):
        (root / name).mkdir(parents=True)
        (root / name / "input.txt").write_text(f"input {name}\n", encoding="utf-8")
    (root / "gamma" / "orca.out").write_text("ORCA fake output\n", encoding="utf-8")
    (tmp_path / "archive").mkdir()
    (tmp_path / "office").mkdir()
    ctx = DashboardContext(calc_dir=root, archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    refs["calc_list_directory"]()

    class Calc:
        pass

    c = Calc()
    c.root, c.refs, c.state = root, refs, refs["xyz_batch_state"]
    c.confirm = refs["calc_delete_confirm"]
    c.yes = c.confirm.children[1]

    def key(action):
        refs["calc_keyboard_action_input"].value = action

    def enter(rel):
        refs["calc_path_input"].value = "/" + rel
        assert c.state["current_path"] == rel

    def select(*names):
        file_list = refs["calc_file_list"]
        values = []
        for name in names:
            hits = [o for o in file_list.options if str(o if not isinstance(o, tuple) else o[0]).endswith(name)]
            assert hits, f"{name} not in {file_list.options}"
            values.append(hits[0][1] if isinstance(hits[0], tuple) else hits[0])
        file_list.value = tuple(values)

    c.key, c.enter, c.select = key, enter, select
    c.confirm_open = lambda: c.confirm.layout.display != "none"
    return c


def test_backspace_with_nothing_selected_goes_up_and_deletes_nothing(calc):
    calc.enter("gamma")
    calc.key("deselect")
    calc.key("up")
    assert calc.state["current_path"] == ""
    assert not calc.confirm_open()
    assert (calc.root / "gamma" / "orca.out").exists()


def test_delete_key_with_nothing_selected_never_offers_the_open_folder(calc):
    calc.enter("gamma")
    calc.key("deselect")
    calc.key("delete")
    assert not calc.confirm_open()
    assert "current folder" not in calc.refs["calc_delete_label"].value
    assert "Select the items" in calc.refs["calc_delete_status"].value
    calc.yes.click()  # even a stray Yes has nothing to act on
    assert (calc.root / "gamma").is_dir()
    assert (calc.root / "gamma" / "orca.out").exists()


def test_delete_key_deletes_exactly_the_items_it_named(calc):
    calc.enter("gamma")
    calc.select("orca.out")
    calc.key("delete")
    assert calc.confirm_open()
    assert "orca.out" in calc.refs["calc_delete_label"].value
    # The selection moves on before Yes; Yes still means what it asked.
    calc.select("input.txt")
    calc.yes.click()
    assert not (calc.root / "gamma" / "orca.out").exists()
    assert (calc.root / "gamma" / "input.txt").exists()
    assert (calc.root / "gamma").is_dir()


def test_escape_answers_the_delete_question_with_no(calc):
    calc.enter("gamma")
    calc.select("orca.out")
    calc.key("delete")
    assert calc.confirm_open()
    calc.key("deselect")
    assert not calc.confirm_open()
    calc.yes.click()
    assert (calc.root / "gamma" / "orca.out").exists()


def test_a_folder_question_from_the_button_does_not_follow_the_reader_up(calc):
    """The Delete button keeps its "delete the open folder" (kept function),
    but a Backspace after it cannot turn Yes against the parent folder."""
    calc.enter("alpha")
    calc.key("deselect")
    calc.refs["calc_delete_btn"].click()
    assert calc.confirm_open()
    assert "current folder" in calc.refs["calc_delete_label"].value
    calc.key("up")
    assert calc.state["current_path"] == ""
    calc.yes.click()
    assert (calc.root / "alpha").is_dir()
    assert (calc.root / "beta").is_dir()
    assert (calc.root / "gamma").is_dir()
    assert calc.root.is_dir()


def test_the_delete_button_still_deletes_the_open_folder_on_yes(calc):
    calc.enter("beta")
    calc.key("deselect")
    calc.refs["calc_delete_btn"].click()
    calc.yes.click()
    assert not (calc.root / "beta").exists()
    assert calc.state["current_path"] == ""
    assert (calc.root / "alpha").is_dir() and (calc.root / "gamma").is_dir()


# -- the Remote Archive ---------------------------------------------------------

@pytest.fixture
def archive(monkeypatch):
    listing = {
        "": [{"name": "proj", "relative_path": "proj", "is_dir": True, "size": 0}],
        "proj": [
            {"name": n, "relative_path": f"proj/{n}", "is_dir": True, "size": 0}
            for n in ("a", "b", "c")
        ],
    }
    deleted = []
    config = {"host": "archive.example.org", "user": "someone",
              "remote_path": "/r", "port": 22}

    def fake_list(host, user, root, port, rel, sort_mode=None):
        rel = rel or ""
        return {"relative_path": rel, "entries": list(listing.get(rel, []))}

    def fake_delete(host, user, root, port, rel_paths):
        deleted.append(list(rel_paths))
        return {"deleted": list(rel_paths), "errors": []}

    monkeypatch.setattr(remote, "load_transfer_settings", lambda: config)
    monkeypatch.setattr(remote, "list_remote_entries", fake_list)
    monkeypatch.setattr(remote, "remote_delete", fake_delete)
    ctx = DashboardContext()
    ctx.run_js = lambda script: None
    _widget, refs = remote.create_tab(ctx)
    refs["path_input"].value = "/r/proj"
    assert {o[1] if isinstance(o, tuple) else o for o in refs["file_list"].options} >= {
        "proj/a", "proj/b", "proj/c"}
    refs["deleted"] = deleted
    return refs


def _yes(refs):
    panel = refs["confirm_panel"]
    assert panel.layout.display != "none"
    yes = [b for b in panel.children if getattr(b, "description", "") == "Yes, delete"]
    assert yes
    yes[0].click()


def test_archive_backspace_goes_up_and_delete_key_needs_a_selection(archive):
    archive["file_list"].value = ()
    archive["keyboard_action_input"].value = "delete"
    assert archive["confirm_panel"].layout.display == "none"
    archive["keyboard_action_input"].value = "up"
    options = {o[1] if isinstance(o, tuple) else o for o in archive["file_list"].options}
    assert "proj" in options
    assert archive["deleted"] == []


def test_archive_delete_names_the_items_and_deletes_only_those(archive):
    archive["file_list"].value = ("proj/a", "proj/b")
    archive["keyboard_action_input"].value = "delete"
    prompt = archive["confirm_panel"].children[0].value
    assert "a</code>" in prompt and "b</code>" in prompt
    archive["file_list"].value = ("proj/b", "proj/c")
    _yes(archive)
    assert archive["deleted"] == [["proj/a", "proj/b"]]


def test_archive_question_does_not_survive_going_up(archive):
    archive["file_list"].value = ("proj/a", "proj/b")
    archive["keyboard_action_input"].value = "delete"
    archive["keyboard_action_input"].value = "up"
    assert archive["confirm_panel"].layout.display == "none"
    assert archive["deleted"] == []
