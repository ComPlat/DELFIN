"""A pre-existing complete CONTROL must keep reading exactly as on origin/main.

Item (3) of the package: a new option defaults to the old behaviour.  The
hand-over resolves a CO2 field either from ``co2_<key>`` or from a bare
``<key>`` line, and that plain-key fallback is old behaviour every existing
CONTROL relies on -- so it must not become the route by which a NEW option
takes effect on a file nobody edited.

All three options here (the orientation-scan switch and the two ORCA keyword
overrides) have a name that used to be dropped by the hand-over.  Reaching
them through the bare name would therefore change what an existing file does
without its owner touching the CO2 CONTROL.  They are reachable only as
``co2_<key>``.

The CONTROL templates are synthetic: the keys and block structure
``delfin/define.py`` ships, with invented values.  No real workspace content.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import child_env

#: A complete CONTROL in the shape ``delfin/define.py`` ships.
COMPLETE_CONTROL = """\
input_file=input.txt
NAME=synthetic_complex
charge=+2
------------------------------------------------------
Solvation:
----------
implicit_solvation_model=CPCM
solvent=DMF
------------------------------------------------------
Redox steps:
------------------
method=OCCUPIER
oxidation_steps=1
reduction_steps=1
------------------------------------------------------
OCCUPIER-Settings:
----------
OCCUPIER_method=auto
occupier_selection=tolerance
------------------------------------------------------
CO2 Coordination:
----------
co2_coordination=on
co2_species_delta=0
{extra}\
"""

SPECIES_XYZ = """\
5

Ni   0.000000   0.000000   0.000000
N    2.000000   0.000000   0.000000
N   -2.000000   0.000000   0.000000
N    0.000000   2.000000   0.000000
N    0.000000  -2.000000   0.000000
"""

STATE = {"0": {"even": {"m": 1, "BS": "", "index": 1}}}


def _write_workspace(root: Path, extra: str) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "CONTROL.txt").write_text(
        COMPLETE_CONTROL.format(extra=extra), encoding="utf-8")
    (root / "initial.xyz").write_text(SPECIES_XYZ, encoding="utf-8")
    (root / ".delfin_occ_auto_state.json").write_text(
        json.dumps(STATE), encoding="utf-8")
    return root


def _generated_control(root: Path) -> str:
    """Run the hand-over through its real entry point; return the generated file."""
    repo_root = Path(__file__).resolve().parents[1]
    env = child_env(root)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(repo_root)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    proc = subprocess.run(
        [sys.executable, "-m", "delfin.co2.chain_setup", "0"],
        cwd=root, env=env, capture_output=True, text=True)
    out = root / "CO2_coordination" / "CONTROL.txt"
    assert out.is_file(), (
        f"the hand-over wrote no CO2 CONTROL; rc={proc.returncode}\n"
        f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")
    return out.read_text()


#: The three options that were read by the coordinator but not transferable,
#: each with a value that would visibly change the coordinator's behaviour.
NEW_OPTIONS = [
    ("perform_orientation_scan", "false"),
    ("orca_keywords", "GFN2-XTB OPT ALPB(DMF)"),
    ("rot_orca_keywords", "GFN2-XTB SP ALPB(DMF)"),
]


@pytest.mark.parametrize("key,value", NEW_OPTIONS)
def test_a_bare_line_for_a_new_option_does_not_change_an_existing_control(
        tmp_path, key, value):
    """The old, bare reading must still be the old reading.

    A CONTROL that carries a bare ``<key>`` line is what the hand-over used to
    DROP.  If adding the option to the transfer set made the bare line start
    reaching the coordinator, that file's meaning changed without its owner
    editing anything -- the exact behaviour item (3) forbids.

    This is the control case for the prefixed test below: same option, same
    value, only the prefix differs, and the expected generated CONTROL differs
    with it.
    """
    bare = _generated_control(_write_workspace(tmp_path / "bare", f"{key}={value}\n"))
    absent = _generated_control(_write_workspace(tmp_path / "absent", ""))

    assert f"{key}={value}" not in bare, (
        f"a bare {key}= line reached the generated CO2 CONTROL; it used to be "
        "dropped, so an existing CONTROL would change meaning without being "
        "edited. This option must be spelled co2_<key>.")
    assert bare == absent, (
        f"a bare {key}= line changed the generated CO2 CONTROL. The bare name "
        "is old behaviour and must keep its old (dropped) meaning.")


@pytest.mark.parametrize("key,value", NEW_OPTIONS)
def test_the_prefixed_form_sets_each_new_option(tmp_path, key, value):
    """The documented spelling reaches the coordinator with the user's value."""
    text = _generated_control(
        _write_workspace(tmp_path / "prefixed", f"co2_{key}={value}\n"))

    assert f"{key}={value}" in text, (
        f"co2_{key}={value} did not reach the generated CO2 CONTROL; the "
        "option the user set has no effect")


def test_the_scan_switch_keeps_the_coordinators_absent_meaning(tmp_path):
    """Unset stays absent, so the coordinator's own default still applies.

    ``CO2_Coordinator6.py:1427-1428`` reads an ABSENT ``perform_orientation_scan``
    as "run the scan".  A CONTROL that never mentions it must therefore keep
    producing a file without the line -- not one with a value of our choosing.
    """
    text = _generated_control(_write_workspace(tmp_path, ""))

    assert "perform_orientation_scan" not in text, (
        "the hand-over invented an orientation-scan switch the user never set")
    assert "orca_keywords" not in text, (
        "the hand-over invented an ORCA keyword override the user never set")
