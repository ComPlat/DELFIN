"""The OCCUPIER → CO2 coordinator hand-over, driven from the workflow's own
entry point and measured at each joint.

The group's central workflow is: one CONTROL.txt drives OCCUPIER (spin-state
screening) and then the CO2 coordinator (orientation single points at a fixed
distance, then a relaxed distance scan).  Between the two halves there is a
hand-over — ``python -m delfin.co2.chain_setup <delta>``, which
``delfin/dashboard/local_runner.py:797`` runs after DELFIN and before
``delfin co2`` — and that hand-over is where a CONTROL option can silently
stop meaning anything: it reads the main CONTROL.txt **raw** (no alias
normalisation, unlike ``delfin.config.OCCUPIER_parser``) and writes out only
the keys it knows about.

Every CONTROL template below is SYNTHETIC: invented values, but the real
*structure* — the key names, the block order, the OCCUPIER block and the
"CO2 Coordination" block as ``delfin/define.py`` ships them.  No content of
any real user's workspace is reproduced here.

What is asserted is the documented outcome for each case, never whether this
machine happens to have ORCA or OPI installed: the coordinator's ORCA launches
are replaced with recorders, so the tests run anywhere.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Synthetic CONTROL.txt — structure only
# ---------------------------------------------------------------------------

#: A complete CONTROL in the shape and key set ``delfin/define.py`` ships,
#: with invented chemistry values.  ``{extra}`` is where a test appends the
#: option under examination, so every other line stays byte-identical and a
#: difference can only come from that option.
COMPLETE_CONTROL = """\
input_file=input.txt
NAME=synthetic_complex
SMILES=[Ni]1(N)(N)N1
charge=+2
------------------------------------------------------
Solvation:
----------
implicit_solvation_model=CPCM
solvent=DMF
XTB_SOLVATOR=no
number_explicit_solv_molecules=2
------------------------------------------------------
SMILES conversion:
----------
smiles_converter=GUPPY
------------------------------------------------------
Global geometry optimisation:
----------
xTB_method=XTB2
XTB_preOPT=no
global_optimizer=GOAT
multiplicity_global_opt=
------------------------------------------------------
Imaginary mode elimination:
----------
IMAG=yes
IMAG_scope=initial
IMAG_option=2
allow_imaginary_freq=-50
------------------------------------------------------
Redox steps:
------------------
method=OCCUPIER
calc_initial=yes
oxidation_steps=1
reduction_steps=1
------------------------------------------------------
OCCUPIER-Settings:
----------
OCCUPIER_method=auto
OCCUPIER_compare=FSPE
occupier_selection=tolerance
maxiter_occupier=125
geom_opt_OCCUPIER=OPT
------------------------------------------------------
CO2 Coordination:
----------
co2_coordination=on
co2_species_delta=0
{extra}"""

#: The species geometry the hand-over copies as ``input.xyz``.
SPECIES_XYZ = """\
5

Ni   0.000000   0.000000   0.000000
N    2.000000   0.000000   0.000000
N   -2.000000   0.000000   0.000000
N    0.000000   2.000000   0.000000
N    0.000000  -2.000000   0.000000
"""


def _write_workspace(root: Path, *, extra: str = "", multiplicities=None) -> Path:
    """A finished-DELFIN workspace as the chain sees it after OCCUPIER.

    ``multiplicities`` is what the OCCUPIER state file records per redox
    delta; the default records the even-electron ladder's first rung.
    """
    root.mkdir(parents=True, exist_ok=True)
    (root / "CONTROL.txt").write_text(
        COMPLETE_CONTROL.format(extra=extra), encoding="utf-8")
    state = multiplicities if multiplicities is not None else {
        "0": {"even": {"m": 1, "BS": "", "index": 1}},
        "-1": {"odd": {"m": 2, "BS": "", "index": 1}},
        "1": {"odd": {"m": 2, "BS": "2,1", "index": 2}},
    }
    # One geometry per species, under the name the chain looks for
    # (``species_delta_to_name``: 0 -> initial, +1 -> ox_step_1, -1 ->
    # red_step_1). OCCUPIER leaves exactly these next to CONTROL.txt.
    from delfin.co2.chain_setup import species_delta_to_name

    for delta in state:
        (root / f"{species_delta_to_name(int(delta))}.xyz").write_text(
            SPECIES_XYZ, encoding="utf-8")
    (root / ".delfin_occ_auto_state.json").write_text(
        json.dumps(state), encoding="utf-8")
    return root


def _co2_control(root: Path) -> dict[str, str]:
    """The generated CO2 CONTROL.txt as the coordinator reads it.

    Parsed by the coordinator's OWN reader, so the values asserted below are
    the values the CO2 run would act on — not a re-parse that could differ.
    """
    from delfin.co2.CO2_Coordinator6 import _minimal_read_control_file

    co2_dir = root / "CO2_coordination"
    assert (co2_dir / "CONTROL.txt").is_file(), (
        "the hand-over wrote no CO2_coordination/CONTROL.txt")
    import os

    cwd = os.getcwd()
    os.chdir(co2_dir)
    try:
        return _minimal_read_control_file("CONTROL.txt")
    finally:
        os.chdir(cwd)


# ---------------------------------------------------------------------------
# The chain, through its real entry point
# ---------------------------------------------------------------------------


def _run_chain(root: Path, delta: int = 0):
    """Run ``python -m delfin.co2.chain_setup <delta>`` in *root*.

    This is verbatim the command ``delfin/dashboard/local_runner.py:797`` runs
    between DELFIN and ``delfin co2``, in the workspace it runs it in.  A
    subprocess, not an in-process call, because the module's own ``__main__``
    block — the one that parses ``sys.argv`` and chooses the job directory —
    is part of the path a user takes, and a test that skipped it would not
    notice that door breaking.

    ``PYTHONPATH`` names this checkout's parent, so the chain always runs the
    code under test rather than whatever ``delfin`` this machine has installed.
    """
    import os
    import subprocess
    import sys

    repo_root = Path(__file__).resolve().parents[1]
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(repo_root)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    proc = subprocess.run(
        [sys.executable, "-m", "delfin.co2.chain_setup", str(delta)],
        cwd=root, env=env, capture_output=True, text=True,
    )
    assert proc.returncode == 0, (
        f"the chain exited {proc.returncode}\n--- stdout ---\n{proc.stdout}\n"
        f"--- stderr ---\n{proc.stderr}")
    return root / "CO2_coordination"


def test_the_chain_writes_the_coordinator_workspace(tmp_path):
    """Control case: passes on the base commit too.

    The hand-over produces exactly the three files ``delfin co2`` then reads
    in its own directory, and the geometry it copies is the species' one.
    """
    root = _write_workspace(tmp_path)
    co2 = _run_chain(root, 0)

    assert (co2 / "CONTROL.txt").is_file()
    assert (co2 / "input.xyz").is_file()
    assert (co2 / "co2.xyz").is_file()
    assert (co2 / "input.xyz").read_bytes() == (root / "initial.xyz").read_bytes()
    # The CO2 file the coordinator places is O, C, O — its carbon is found by
    # symbol, which is why the adapter refuses to assume the first index.
    lines = (co2 / "co2.xyz").read_text().splitlines()[2:]
    assert [ln.split()[0] for ln in lines] == ["O", "C", "O"]


def test_the_chosen_spin_state_reaches_the_coordinator(tmp_path):
    """The heart of the hand-over: OCCUPIER's winner becomes the CO2 run's spin.

    For each redox species the multiplicity the state file records for that
    delta must appear as ``* xyz <charge> <multiplicity>`` in the CO2 CONTROL,
    with the charge shifted by the delta and BrokenSym carried when the state
    file has one.
    """
    root = _write_workspace(tmp_path)
    _run_chain(root, 0)
    ctrl = _co2_control(root)

    assert ctrl["multiplicity"] == 1, (
        f"the initial species' multiplicity is {ctrl['multiplicity']!r}; the "
        "OCCUPIER state file recorded m=1 for delta 0")
    assert ctrl["charge"] == 2, (
        "charge 2 + delta 0 must stay 2, got "
        f"{ctrl['charge']!r}")
    assert not ctrl["broken_sym"], (
        "a spin state without a BS label must not gain one; got "
        f"{ctrl['broken_sym']!r}. (The coordinator's reader maps the emitted "
        "empty value to None, so both '' and None are correct here.)")

    # The oxidised species: delta +1, m=2 with a BS label.
    root = _write_workspace(tmp_path / "ox")
    _run_chain(root, 1)
    ctrl = _co2_control(root)
    assert ctrl["multiplicity"] == 2
    assert ctrl["charge"] == 3, "charge 2 + delta 1 must be 3"
    assert re.search(r"BrokenSym\s+2,1", ctrl["broken_sym"]), (
        f"the BS label 2,1 from the state file did not reach the CO2 CONTROL: "
        f"{ctrl['broken_sym']!r}")


def test_an_absent_co2_option_is_absent_from_the_generated_control(
    tmp_path, monkeypatch
):
    """A new option the user did not write must not appear.

    This is the invariant that keeps every existing complete CONTROL reading
    exactly as it did on origin/main: the hand-over only emits a line for a
    key it can resolve, so an option nobody set is not invented.
    """
    root = _write_workspace(tmp_path)
    _run_chain(root, 0)
    text = (root / "CO2_coordination" / "CONTROL.txt").read_text()

    # Not written, so the coordinator's own default applies — unchanged.
    assert "perform_orientation_scan" not in text, (
        "an option the user did not write appeared in the generated CONTROL")
    assert re.search(r"^scan_steps=25$", text, re.MULTILINE), (
        "the shipped default scan_steps=25 must still be emitted")


# ---------------------------------------------------------------------------
# The options the coordinator reads but the hand-over could not pass on
# ---------------------------------------------------------------------------


def test_the_orientation_scan_switch_can_be_turned_off_from_the_main_control(
    tmp_path, monkeypatch
):
    """``co2_perform_orientation_scan=false`` must actually skip the scan.

    The coordinator reads ``perform_orientation_scan`` and treats an ABSENT
    key as "run the scan" (``CO2_Coordinator6.py:1427-1428``).  So an option
    the hand-over cannot pass on does not fall back to a safe default — it
    silently runs the expensive half of the workflow the user just turned off.

    Red on the base commit: the key is not in the transfer set, so no line is
    emitted, the coordinator reads ``None`` and scans anyway.
    """
    root = _write_workspace(tmp_path, extra="co2_perform_orientation_scan=false\n")
    _run_chain(root, 0)
    ctrl = _co2_control(root)

    assert "perform_orientation_scan" in ctrl, (
        "the switch never reached the CO2 CONTROL, so the coordinator takes "
        "the ABSENT branch and runs the orientation scan the user turned off")
    assert ctrl["perform_orientation_scan"] is False, (
        f"the switch reached the coordinator as "
        f"{ctrl['perform_orientation_scan']!r}, not False")


def test_the_two_keyword_overrides_reach_the_coordinator(tmp_path):
    """A custom ORCA keyword line must reach the job it names.

    ``rot_orca_keywords`` (orientation) and ``orca_keywords`` (distance scan)
    are the documented escape hatch from the per-keyword builder
    (``CO2_Coordinator6.py:1367-1377``), and both are read from the CO2
    CONTROL only.  A user setting them on the main CONTROL must get them —
    otherwise the only way to reach them is to hand-edit a file the hand-over
    regenerates.
    """
    root = _write_workspace(
        tmp_path,
        extra="co2_rot_orca_keywords=GFN2-XTB SP ALPB(DMF)\n"
              "co2_orca_keywords=GFN2-XTB OPT ALPB(DMF)\n",
    )
    _run_chain(root, 0)
    ctrl = _co2_control(root)

    assert ctrl.get("rot_orca_keywords") == "GFN2-XTB SP ALPB(DMF)", (
        f"orientation keyword override lost: {ctrl.get('rot_orca_keywords')!r}")
    assert ctrl.get("orca_keywords") == "GFN2-XTB OPT ALPB(DMF)", (
        f"scan keyword override lost: {ctrl.get('orca_keywords')!r}")
