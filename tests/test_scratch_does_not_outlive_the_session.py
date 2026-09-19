"""A scratch directory belongs to the process that made it.

Measured on the login node on 2026-09-19: ``/tmp`` held **1071**
directories named ``delfin-gfnff-topo-*``. They are the GFN-FF bonding
perception the structure editor keeps for one molecule — made with
``tempfile.mkdtemp``, removed by ``_drop_gfn_topology`` when the molecule
changes, and removed by nothing at all when the session ends. A browser
tab closed, a kernel killed, a machine rebooted: each leaves its folder
behind for good, and the editor makes a fresh one per molecule, per
charge, per session.

``atexit`` is half an answer. The way these sessions actually end is a
kernel that is killed, and a killed process runs no handler. So the
folder carries a stamp saying which process owns it, and a later run
sweeps the ones whose owner is gone.

"Whose owner is gone" is the question a pid alone cannot answer — the
number is handed to somebody else afterwards. It is asked here the same
way ``where.py`` asks it about a dashboard, with the process start time
beside the pid, and out of ONE implementation: two answers to one
question drift, and this codebase has paid for that twice.

What must never happen is the sweep removing a directory it does not own:
an unstamped one, one whose owner is alive, one stamped on another host.
"""

from __future__ import annotations

import json
import os
import socket

import pytest

from delfin.agent import scratch


PREFIX = "delfin-test-scratch-"


@pytest.fixture()
def base(tmp_path):
    """A base of our own — the real /tmp is not a test fixture."""
    room = tmp_path / "base"
    room.mkdir()
    return room


def _plant(base, *, stamp=None, name="planted"):
    folder = base / name
    folder.mkdir()
    if stamp is not None:
        (folder / scratch.OWNER_NAME).write_text(
            json.dumps(stamp), encoding="utf-8")
    return folder


# -- the stamp --------------------------------------------------------------

def test_a_folder_knows_who_made_it(base):
    folder = scratch.owned_dir(PREFIX, base=base)
    record = json.loads(
        (folder / scratch.OWNER_NAME).read_text(encoding="utf-8"))
    assert record["pid"] == os.getpid()
    assert record["host"] == socket.gethostname()
    assert record["proc_start"], (
        "without the start time a recycled pid keeps the folder alive "
        "forever")


def test_it_is_made_where_it_was_asked_for(base):
    assert scratch.owned_dir(PREFIX, base=base).parent == base


# -- the sweep --------------------------------------------------------------

def test_a_folder_whose_owner_is_gone_is_swept(base):
    gone = _plant(base, name=PREFIX + "gone",
                  stamp={"pid": 2 ** 22 - 1, "proc_start": "1",
                         "host": socket.gethostname()})
    scratch.sweep(PREFIX, base=base)
    assert not gone.exists(), "the leak this whole module exists to stop"


def test_a_folder_whose_owner_is_alive_survives(base):
    mine = scratch.owned_dir(PREFIX, base=base)
    scratch.sweep(PREFIX, base=base)
    assert mine.exists(), (
        "a sweep that takes a running session's scratch is worse than "
        "the leak")


def test_a_recycled_pid_does_not_pass_for_the_owner(base):
    """The number is in use — by somebody else. That is the case the
    start time is carried for."""
    stale = _plant(base, name=PREFIX + "recycled",
                   stamp={"pid": os.getpid(), "proc_start": "not-ours",
                          "host": socket.gethostname()})
    scratch.sweep(PREFIX, base=base)
    assert not stale.exists()


def test_an_unstamped_directory_is_left_alone(base):
    """We did not make it, so we do not know that it is finished."""
    theirs = _plant(base, name=PREFIX + "unstamped")
    scratch.sweep(PREFIX, base=base)
    assert theirs.exists()


def test_a_folder_from_another_host_is_left_alone(base):
    """A shared scratch filesystem carries other machines' work, and a
    pid there means nothing here."""
    elsewhere = _plant(base, name=PREFIX + "elsewhere",
                       stamp={"pid": 1, "proc_start": "1",
                              "host": socket.gethostname() + "-other"})
    scratch.sweep(PREFIX, base=base)
    assert elsewhere.exists()


def test_a_different_prefix_is_not_touched(base):
    other = _plant(base, name="someone-elses-",
                   stamp={"pid": 2 ** 22 - 1, "proc_start": "1",
                          "host": socket.gethostname()})
    scratch.sweep(PREFIX, base=base)
    assert other.exists()


def test_a_sweep_of_nothing_is_quiet(base):
    assert scratch.sweep(PREFIX, base=base) == []


def test_a_broken_stamp_is_left_alone(base):
    bad = base / (PREFIX + "broken")
    bad.mkdir()
    (bad / scratch.OWNER_NAME).write_text("{not json", encoding="utf-8")
    scratch.sweep(PREFIX, base=base)
    assert bad.exists(), "unreadable is not the same as finished"


def test_the_sweep_never_raises(base):
    assert scratch.sweep(PREFIX, base=base / "does-not-exist") == []


# -- released on the way out, when there is a way out -----------------------

def test_releasing_takes_it_now(base):
    folder = scratch.owned_dir(PREFIX, base=base)
    scratch.release(folder)
    assert not folder.exists()


def test_the_exit_handler_takes_what_this_process_made(base):
    scratch.owned_dir(PREFIX, base=base)
    scratch.owned_dir(PREFIX, base=base)
    scratch._release_mine()
    assert list(base.iterdir()) == [], (
        "a clean shutdown should not need the sweep to catch up later")


def test_releasing_something_twice_is_not_an_error(base):
    folder = scratch.owned_dir(PREFIX, base=base)
    scratch.release(folder)
    scratch.release(folder)


# -- and it is wired into the thing that leaked -----------------------------

def test_the_editors_topology_folder_is_stamped(tmp_path, monkeypatch):
    """A helper that is never called protects nothing: the folder the
    structure editor actually makes must carry the stamp."""
    pytest.importorskip("ipywidgets")
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    water = ("3\nwater\n"
             "O   0.000000   0.000000   0.000000\n"
             "H   0.758602   0.000000   0.504284\n"
             "H  -0.758602   0.000000   0.504284\n")

    room = tmp_path / "room"
    for name in ("calc", "archive", "office"):
        (room / name).mkdir(parents=True)
    ctx = DashboardContext(calc_dir=room / "calc", archive_dir=room / "archive",
                           office_dir=room / "office")
    ctx.run_js = lambda _script: None
    state: dict = {}
    part = structure_editor.build(
        ctx, state=state, coords_widget=widgets.Textarea(value=water),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    part.submit_ff_dd.value = "gfnff"

    # No seed: the atom count of the source has to match for the one-cycle
    # perception to run, and a real xtb call is not what this test is about.
    state["gfn_topology_source"] = "1\nnot the same molecule\nH 0 0 0\n"

    folder = part._gfn_topology_dir(water)
    assert (folder / scratch.OWNER_NAME).is_file(), (
        "the editor's own folder is the one that left 1071 copies in /tmp")
