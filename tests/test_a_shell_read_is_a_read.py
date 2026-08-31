"""Evidence gathered through the shell has to count as evidence.

Field report 20260831-083210_ka_ew7404: read the singlet-triplet gap out
of nine ORCA runs under /pfs/.../calc/TADFs-Jeneesh_*. The agent did it
correctly and the grounding guards fought it for two extra turns, because
the data lives OUTSIDE the workspace and is therefore reached with `grep`
through the shell rather than with `grep_file`.

Nothing recorded that. `_OBSERVATION_TOOLS` listed the typed readers and
not bash, so after 28 tool calls the observed-files ledger was empty while
the guards still ran against it. Three consequences, each reproduced
before it was fixed:

  * the answer cited "S1.out and T1.out in the respective ESD/ folders";
    the citation check resolved those basenames against the repo root,
    found nothing, and classified them "nonexistent" -- the hardest flag,
    which forces a correction turn. The absolute directory that clears
    both names stood in the same sentence.

  * that correction turn cost 91.8 s and six tool calls and produced a
    byte-identical table, because nothing had been wrong with the first.

  * the retry used grep too, so `new_files` stayed empty and the engine
    concluded it had "read nothing" -- and re-stated the ORIGINAL claims
    underneath the CORRECTED answer. The caveat ended up telling the user
    to distrust 0.858 eV and 2.390 eV, neither of which appears in the
    table it sits under.

A second, independent gap surfaced while measuring the first: the number
pool DID see the shell output, so the gaps in Hartree were grounded as
differences of observed values -- and the same gaps in eV were flagged,
because a unit conversion was not a derivation the guard knew. An answer
that reports 0.031522 Hartree and 0.858 eV side by side had one of them
called unsourced.

Both fixes are narrow on purpose. Recording every bash call as evidence
would re-open #124 (one call disarms the scanner for the whole turn);
accepting any scale factor would ground any number at all. So: only
readers, only when they exited 0 and printed something, only paths that
exist -- and only a closed table of physical conversion factors, at the
precision the answer itself printed.
"""

from __future__ import annotations

import json
import random
import tempfile
from pathlib import Path

import pytest

from delfin.agent import verify_guard as vg
from delfin.agent.api_client import _bash_read_targets, _observe_read_files


@pytest.fixture
def calc():
    """A calculation folder OUTSIDE the workspace — the reported shape."""
    with tempfile.TemporaryDirectory(prefix="ws-") as ws_dir, \
            tempfile.TemporaryDirectory(prefix="calc-") as calc_dir:
        ws = Path(ws_dir)
        esd = Path(calc_dir) / "TADFs-Jeneesh_3" / "ESD"
        esd.mkdir(parents=True)
        (esd / "S1.out").write_text(
            "FINAL SINGLE POINT ENERGY   -1705.197734757559\n")
        (esd / "T1.out").write_text(
            "FINAL SINGLE POINT ENERGY   -1705.219605478025\n")
        yield ws, esd


def _ran(cmd: str, *, exit_code: int = 0, stdout: str = "some output\n"):
    return json.dumps({"exit_code": exit_code, "stdout": stdout,
                       "stderr": "", "command": cmd})


# ---------------------------------------------------------------------------
# What a command demonstrably read
# ---------------------------------------------------------------------------

def test_the_files_a_grep_opened_are_recorded(calc):
    _ws, esd = calc
    got = _bash_read_targets(
        f'grep "FINAL SINGLE POINT ENERGY" {esd}/S1.out {esd}/T1.out')
    assert sorted(Path(p).name for p in got) == ["S1.out", "T1.out"]


def test_the_search_pattern_is_not_recorded_as_a_file(calc):
    """Existence does the disambiguation, so no per-tool flag grammar."""
    _ws, esd = calc
    got = _bash_read_targets(f'grep "S1.out" {esd}/T1.out')
    assert [Path(p).name for p in got] == ["T1.out"]


def test_a_relative_path_needs_the_cwd(calc):
    _ws, esd = calc
    assert _bash_read_targets("grep ENERGY S1.out") == []
    got = _bash_read_targets("grep ENERGY S1.out", esd)
    assert [Path(p).name for p in got] == ["S1.out"]


def test_a_glob_is_skipped_rather_than_expanded(calc):
    """Expanding it would record files the command may never have opened."""
    _ws, esd = calc
    assert _bash_read_targets(f'grep ENERGY {esd}/*.out') == []


def test_a_command_that_is_not_a_reader_grounds_nothing(calc):
    """Evidence must not be minted by a command whose effect nobody checked."""
    _ws, esd = calc
    for cmd in (f"cp {esd}/S1.out /tmp/x.out",
                f"rm {esd}/S1.out",
                f"python3 {esd}/S1.out",
                f"tar -cf /tmp/a.tar {esd}/S1.out"):
        assert _bash_read_targets(cmd) == [], cmd


def test_one_unknown_segment_does_not_poison_the_others(calc):
    """A compound command is only as trustworthy as its least understood
    part — which costs a missing entry here, never a wrong one."""
    _ws, esd = calc
    got = _bash_read_targets(f"grep ENERGY {esd}/S1.out && ./deploy.sh")
    assert [Path(p).name for p in got] == ["S1.out"]


# ---------------------------------------------------------------------------
# When a run counts as having shown something
# ---------------------------------------------------------------------------

def test_a_successful_grep_reaches_the_ledger(calc):
    ws, esd = calc
    cmd = f"grep ENERGY {esd}/S1.out"
    obs: set = set()
    _observe_read_files(obs, "bash", {"command": cmd}, _ran(cmd), workspace=ws)
    assert len(obs) == 1


def test_the_namespaced_mcp_form_reaches_it_too(calc):
    """The reported session ran every command as mcp__kit-coding__bash."""
    ws, esd = calc
    cmd = f"grep ENERGY {esd}/S1.out"
    obs: set = set()
    _observe_read_files(obs, "mcp__kit-coding__bash", {"command": cmd},
                        _ran(cmd), workspace=ws)
    assert len(obs) == 1


@pytest.mark.parametrize("kwargs", [
    {"exit_code": 1, "stdout": ""},          # no match
    {"exit_code": 0, "stdout": "   \n"},     # matched nothing visible
    {"exit_code": 2, "stdout": "x\n"},       # failed
])
def test_a_run_that_showed_nothing_grounds_nothing(calc, kwargs):
    """Same rule the typed readers answer to: one deliberately
    non-matching grep must not disarm the guard for a whole file."""
    ws, esd = calc
    cmd = f"grep NOPE {esd}/S1.out"
    obs: set = set()
    _observe_read_files(obs, "bash", {"command": cmd}, _ran(cmd, **kwargs),
                        workspace=ws)
    assert obs == set()


def test_a_background_job_has_not_read_anything_yet(calc):
    """bash_background returns a job id; the file is still unopened."""
    ws, esd = calc
    cmd = f"grep ENERGY {esd}/S1.out"
    obs: set = set()
    _observe_read_files(obs, "bash_background", {"command": cmd},
                        json.dumps({"job_id": "j1"}), workspace=ws)
    assert obs == set()


# ---------------------------------------------------------------------------
# The consequence the user actually saw
# ---------------------------------------------------------------------------

ANSWER = ("Alle Werte stammen aus den `FINAL SINGLE POINT ENERGY`-Zeilen "
          "der Dateien `S1.out` und `T1.out` in den jeweiligen `ESD/`-Ordnern.")


def test_a_cited_basename_is_not_called_invented(calc):
    ws, esd = calc
    cmd = f"grep ENERGY {esd}/S1.out {esd}/T1.out"
    obs: set = set()
    _observe_read_files(obs, "mcp__kit-coding__bash", {"command": cmd},
                        _ran(cmd), workspace=ws)
    flags = vg.scan_for_ungrounded_code_claims(
        ANSWER, repo_root=str(ws), observed_files=obs)
    assert [f.path for f in flags if f.kind == "nonexistent"] == []


def test_without_the_ledger_the_same_answer_is_still_judged(calc):
    """The guard must not have been switched off — only informed."""
    ws, _esd = calc
    flags = vg.scan_for_ungrounded_code_claims(
        ANSWER, repo_root=str(ws), observed_files=set())
    assert sorted(f.path for f in flags if f.kind == "nonexistent") == [
        "S1.out", "T1.out"]


# ---------------------------------------------------------------------------
# The same quantity in another unit is the same quantity
# ---------------------------------------------------------------------------

# The FSPE values the greps returned for two of the nine systems.
POOL = [-1705.188083533148, -1705.219605478025,
        -1694.015332057528, -1694.103175837060]


def test_a_gap_converted_to_ev_is_as_grounded_as_the_gap(calc):
    text = ("Kleinste ΔEST: 0.031522 Hartree = 0.858 eV. "
            "Größte ΔEST: 0.087844 Hartree = 2.390 eV.")
    assert vg.scan_for_unsourced_quantities(text, numbers=POOL) == []


def test_the_other_units_ride_the_same_table():
    gap = abs(POOL[0] - POOL[1])
    for factor, unit in ((627.5094740631, "kcal/mol"),
                         (2625.4996394799, "kJ/mol"),
                         (219474.6313632, "cm-1")):
        text = f"Die Lücke beträgt {gap * factor:.3f} {unit}."
        assert vg.scan_for_unsourced_quantities(text, numbers=POOL) == [], unit


def test_an_invented_number_is_still_flagged():
    flags = vg.scan_for_unsourced_quantities(
        "Die Lücke betrug 4.271 eV.", numbers=POOL)
    assert [f.quantity for f in flags] == ["4.271 eV"]


def test_a_wrong_conversion_of_a_real_gap_is_still_flagged():
    """The point is not to wave energies through — 0.031522 Ha is
    6918.4 cm-1, and an answer that prints 6917 has not converted it."""
    flags = vg.scan_for_unsourced_quantities(
        "Die Lücke beträgt 6917 cm-1.", numbers=POOL)
    assert [f.quantity for f in flags] == ["6917 cm-1"]


def test_the_conversion_costs_less_than_one_percent():
    """The COST side, pinned.

    Tying the tolerance to the claim's printed precision instead of the
    scanner's absolute floor is the whole reason this is affordable: with
    the floor, dividing a claim by 96.485 opened a ±0.48 eV window and
    29% of random eV-scale values came out "grounded". Measured here so a
    later widening of _UNIT_FACTORS or of the tolerance cannot pass
    unnoticed.
    """
    rng = random.Random(20260831)
    grounded = sum(
        bool(vg.scan_for_unsourced_quantities(
            f"Der Wert ist {rng.uniform(0.01, 10.0):.3f} eV.", numbers=POOL)
            == [])
        for _ in range(2000)
    )
    assert grounded / 2000 < 0.01, f"{grounded}/2000 wrongly grounded"
