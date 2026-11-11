"""A file list never sends the browser a selection that is not in its options.

The Home tab (the Calculations browser rooted at ~), the Remote Archive and
the Literature tab show their files in an ipywidgets ``SelectMultiple``.
Setting its ``options`` from Python sends two updates: the new labels, then
``index = ()``. In between the browser holds the old ``index`` against the
new labels; when the list shrank under a selected file the view's
``updateSelection`` touches an option that is not there and the page logs
``Error setting state: Cannot set properties of undefined (setting
'selected')``. ``select_state.replace_options`` sends labels and selection
as one update.

Input: a real ``SelectMultiple`` whose comm records the messages, and a
model of the browser view that applies each message as ipywidgets'
``SelectMultipleView.update`` does. Properties asserted:

* a plain ``options =`` under a selection sends an invalid state (the hazard
  is real in the installed ipywidgets);
* ``replace_options`` never does, whatever the old selection, and the kept
  selection arrives with its labels;
* Python observers see the same sequence of values either way, so no
  handler of the file lists runs differently;
* no file list in the three tabs sets ``.options`` directly.
"""

from __future__ import annotations

import ast
import uuid
from pathlib import Path

import pytest

widgets = pytest.importorskip("ipywidgets")

from delfin.dashboard.select_state import replace_options  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
DASHBOARD = REPO / "delfin" / "dashboard"
FILE_LIST_TABS = ("tab_calculations_browser.py", "tab_remote_archive.py",
                  "tab_literature.py")


class _RecordingComm:
    kernel = "recording"

    def __init__(self):
        self.comm_id = uuid.uuid4().hex
        self.messages: list[dict] = []

    def send(self, data=None, buffers=None, **_kw):
        self.messages.append(data)

    def on_msg(self, *_a, **_kw):
        pass

    def close(self, *_a, **_kw):
        pass


class _BrowserSelectMultiple:
    """What ``SelectMultipleView`` does with each ``update`` message.

    The browser's model applies the whole message, then the view updates
    once: rebuild the ``<option>`` list if the labels changed, then mark
    ``listbox.options[i].selected`` for every ``i`` in ``index``. An ``i``
    past the end is the "Cannot set properties of undefined" error.
    """

    def __init__(self, labels, index):
        self.labels = list(labels)
        self.index = list(index)
        self.errors: list[str] = []

    def apply(self, msg):
        if not msg or msg.get("method") != "update":
            return
        state = msg.get("state", {})
        if "_options_labels" in state:
            self.labels = list(state["_options_labels"])
        if "index" in state:
            self.index = list(state["index"] or [])
        if "_options_labels" in state or "index" in state:
            for i in self.index:
                if not 0 <= i < len(self.labels):
                    self.errors.append(
                        f"index {i} with {len(self.labels)} options "
                        f"(message keys {sorted(state)})")


def _file_list(options, value):
    w = widgets.SelectMultiple(options=options, value=value)
    comm = _RecordingComm()
    w.comm = comm
    browser = _BrowserSelectMultiple(w._options_labels, w.index)
    return w, comm, browser


def _replay(comm, browser):
    for msg in comm.messages:
        browser.apply(msg)
    return browser


def test_a_plain_options_assignment_under_a_selection_sends_an_invalid_state():
    w, comm, browser = _file_list(list("abcdef"), ("f",))
    w.options = ["x"]
    assert _replay(comm, browser).errors, (
        "the installed ipywidgets now sends options and index together; "
        "select_state.replace_options may no longer be needed")


@pytest.mark.parametrize("old_value", [(), ("a",), ("f",), ("b", "e", "f")])
@pytest.mark.parametrize("new_options", [[], ["x"], ["a", "b"], list("abcdefgh")])
def test_replace_options_never_sends_a_selection_outside_the_list(old_value, new_options):
    w, comm, browser = _file_list(list("abcdef"), old_value)
    replace_options(w, new_options)
    _replay(comm, browser)
    assert browser.errors == []
    assert browser.labels == list(new_options)
    assert browser.index == []
    assert w.value == ()


def test_a_kept_selection_arrives_with_its_labels():
    w, comm, browser = _file_list(
        [("one", "1"), ("two", "2"), ("three", "3"), ("four", "4")], ("4",))
    replace_options(w, [("four", "4"), ("five", "5")], ("4",))
    _replay(comm, browser)
    assert browser.errors == []
    assert browser.labels == ["four", "five"]
    assert browser.index == [0]
    assert w.value == ("4",)
    # one message: labels and selection are never apart
    assert len([m for m in comm.messages if m.get("method") == "update"]) == 1


def test_python_handlers_see_the_same_values_as_before():
    def observed(change_options):
        w = widgets.SelectMultiple(options=list("abcdef"), value=("f",))
        w.comm = _RecordingComm()
        seen = []
        w.observe(lambda ch: seen.append(ch["new"]), names="value")
        change_options(w)
        return seen

    def before(w):
        w.options = ["f", "g"]
        w.value = ("f",)

    def after(w):
        replace_options(w, ["f", "g"], ("f",))

    assert observed(before) == observed(after) == [(), ("f",)]


def test_a_forced_rerender_in_two_steps_stays_valid():
    """The XYZ batch list clears, then refills (to repaint changed marks)."""
    w, comm, browser = _file_list([("□ a", "a"), ("□ b", "b")], ("b",))
    replace_options(w, ())
    replace_options(w, [("✔ b", "b")], ("b",))
    _replay(comm, browser)
    assert browser.errors == []
    assert browser.labels == ["✔ b"] and browser.index == [0]


def _select_multiple_names(tree: ast.AST) -> set[str]:
    names = set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and isinstance(node.value, ast.Call)
                and isinstance(node.value.func, ast.Attribute)
                and node.value.func.attr == "SelectMultiple"):
            names.update(t.id for t in node.targets if isinstance(t, ast.Name))
    return names


@pytest.mark.parametrize("filename", FILE_LIST_TABS)
def test_no_file_list_sets_its_options_directly(filename):
    tree = ast.parse((DASHBOARD / filename).read_text(encoding="utf-8"))
    lists = _select_multiple_names(tree)
    assert lists, f"{filename}: no SelectMultiple found (renamed?)"
    direct = [
        f"{filename}:{node.lineno} {target.value.id}.options = ..."
        for node in ast.walk(tree) if isinstance(node, ast.Assign)
        for target in node.targets
        if isinstance(target, ast.Attribute) and target.attr == "options"
        and isinstance(target.value, ast.Name) and target.value.id in lists
    ]
    assert direct == [], (
        "set a file list's options with select_state.replace_options, so "
        "labels and selection reach the browser in one update: " + "; ".join(direct))
