"""``delfin <dir> --define`` leaves a workspace that can be run.

CONTROL.txt was written into the directory that was named, and input.txt into
whichever directory the command happened to be typed in -- the CONTROL path
was resolved against the workspace, the input path against the process's own
working directory. Nothing failed: the user was handed a workspace whose
CONTROL file names a geometry that is not there, and found out at the next
run. The program's own --help promised the opposite.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from delfin.define import create_control_file

GEOMETRY = """3
water
O 0.000000 0.000000 0.000000
H 0.758602 0.000000 0.504284
H -0.758602 0.000000 0.504284
"""


@pytest.fixture
def elsewhere(tmp_path, monkeypatch):
    """A working directory that is not the workspace."""
    caller = tmp_path / "wherever_the_user_stands"
    caller.mkdir()
    monkeypatch.chdir(caller)
    return caller


def test_both_files_land_in_the_workspace(tmp_path, elsewhere):
    workspace = tmp_path / "project"
    workspace.mkdir()

    create_control_file(filename=str(workspace / "CONTROL.txt"))

    assert (workspace / "CONTROL.txt").is_file()
    assert (workspace / "input.txt").is_file()          # this is the one that went missing
    assert not (elsewhere / "input.txt").exists()       # and this is where it used to go


def test_a_geometry_is_converted_into_the_workspace(tmp_path, elsewhere):
    workspace = tmp_path / "project"
    workspace.mkdir()
    xyz = elsewhere / "water.xyz"
    xyz.write_text(GEOMETRY, encoding="utf-8")

    create_control_file(filename=str(workspace / "CONTROL.txt"), input_file=str(xyz))

    written = (workspace / "input.txt").read_text(encoding="utf-8")
    assert written.splitlines()[0].startswith("O ")     # the two header lines are dropped
    assert len(written.splitlines()) == 3
    assert not (elsewhere / "input.txt").exists()


def test_the_workspace_is_created_when_it_is_not_there(tmp_path, elsewhere):
    workspace = tmp_path / "not_yet" / "project"

    create_control_file(filename=str(workspace / "CONTROL.txt"))

    assert (workspace / "input.txt").is_file()


def test_running_in_place_is_unchanged(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    create_control_file(filename="CONTROL.txt")

    assert (tmp_path / "CONTROL.txt").is_file()
    assert (tmp_path / "input.txt").is_file()


def test_an_existing_input_file_is_left_alone(tmp_path, elsewhere):
    workspace = tmp_path / "project"
    workspace.mkdir()
    (workspace / "input.txt").write_text("C 0.0 0.0 0.0\n", encoding="utf-8")

    create_control_file(filename=str(workspace / "CONTROL.txt"))

    assert (workspace / "input.txt").read_text(encoding="utf-8") == "C 0.0 0.0 0.0\n"


def test_the_help_and_the_behaviour_agree(tmp_path, elsewhere):
    """--help promises both files in the named directory; it now holds."""
    from delfin.cli_helpers import _build_parser

    help_text = _build_parser().format_help()
    assert "/path/to/project/CONTROL.txt" in help_text and "/path/to/project/input.txt" in help_text

    workspace = tmp_path / "path_to_project"
    workspace.mkdir()
    create_control_file(filename=str(workspace / "CONTROL.txt"))
    assert {"CONTROL.txt", "input.txt"} <= {p.name for p in workspace.iterdir()}
