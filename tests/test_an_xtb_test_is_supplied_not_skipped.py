"""xtb is supplied to a test, not a reason to skip it.

The rule is in test_a_test_is_bound_to_delfin_not_to_a_machine.py:
supply the binary, or record its output; skip last, with a named reason.
These tests were skipping.

Measured on 2026-10-05 across the 37 test files that gate on xtb:

    PATH                        passed  skipped  failed
    no xtb                        782     122       0
    real xtb                      899       5       0
    a stub that computes nothing  892       7       5

So 110 of the 122 skipped tests never consumed xtb output at all — they
gate on its presence and then exercise a command line, a parse, a guard
or an error path, which is DELFIN's own logic. And over a whole run with
real xtb the suite invokes it **15 times, in 10 distinct shapes**: the
gate was far more conservative than the need.

Five tests in four files do consume output and keep their skips. A
recorded energy replayed back to a test about that energy is not a
check.
"""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess

import pytest


_STUB = pathlib.Path(__file__).resolve().parent / "stubs" / "xtb"


def test_the_stand_in_ships_with_the_suite_and_is_executable():
    """"Universal" must not mean "install xtb first"."""
    assert _STUB.is_file(), _STUB
    assert os.access(_STUB, os.X_OK), "the stand-in is not executable"


def test_it_answers_the_two_questions_presence_gated_code_asks():
    out = subprocess.run([str(_STUB), "--version"], capture_output=True,
                         text=True, timeout=30)
    assert out.returncode == 0
    assert "version" in out.stdout.lower()


def test_it_refuses_loudly_instead_of_inventing_a_result():
    """A silent zero exit with empty output would make a parse return
    nothing and land the failure three layers from the cause."""
    out = subprocess.run([str(_STUB), "in.xyz", "--gfn", "2", "--opt"],
                         capture_output=True, text=True, timeout=30)
    assert out.returncode != 0
    assert not out.stdout.strip(), "a refusal must not look like output"
    assert "no output is recorded" in out.stderr
    assert "in.xyz --gfn 2 --opt" in out.stderr, (
        "the refusal must name the invocation, or the next person cannot "
        "tell which recording is missing")


def test_the_fixture_puts_it_on_path_when_the_machine_has_none(
        monkeypatch, request):
    monkeypatch.setenv("PATH", "/nonexistent-for-this-test")
    got = request.getfixturevalue("xtb_on_path")
    assert got is not None, "no xtb anywhere, and the fixture supplied none"
    assert shutil.which("xtb"), "the stand-in did not reach PATH"
    assert pathlib.Path(shutil.which("xtb")).resolve() == _STUB.resolve()


def test_the_fixture_leaves_a_real_xtb_alone(monkeypatch, request, tmp_path):
    """A developer's run must still exercise the real program; the
    stand-in is what CI and a bare machine get."""
    fake_real = tmp_path / "xtb"
    fake_real.write_text("#!/bin/sh\nexit 0\n")
    fake_real.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    got = request.getfixturevalue("xtb_on_path")
    assert got is None, "the fixture replaced an xtb that was already there"
    assert pathlib.Path(shutil.which("xtb")).resolve() == fake_real.resolve()


def test_the_converted_files_no_longer_skip_on_the_host():
    """The four files converted in this change must not reintroduce a
    host-dependent skip for the whole file. One NAMED skip on a single
    test is allowed and is what test_gfn_methods_in_the_viewer keeps."""
    tests_dir = pathlib.Path(__file__).resolve().parent
    for name in ("test_the_budget_prices_a_relaxed_path.py",
                 "test_gfn_methods_in_the_viewer.py",
                 "test_what_the_answer_already_computed.py",
                 "test_asking_what_a_structure_is.py"):
        text = (tests_dir / name).read_text(encoding="utf-8")
        assert '_needs_xtb = pytest.mark.usefixtures("xtb_on_path")' in text, (
            f"{name} gates on the host again")


@pytest.mark.parametrize("tool", ["xtb"])
def test_the_stub_directory_holds_only_stand_ins(tool):
    """A directory that lands on PATH must not grow anything unexpected:
    everything in it shadows a real program for every test that uses the
    fixture."""
    found = sorted(p.name for p in _STUB.parent.iterdir() if p.is_file())
    assert found == [tool], found


def test_a_supplied_test_does_not_assert_how_fast_the_machine_is():
    """The second kind of machine dependence, found by supplying xtb.

    `test_a_run_shorter_than_the_reading_interval_still_hands_its_path_over`
    asserted `spent < 0.2` — a wall-clock duration. It had never run in
    CI before this change, and the first time it did, the runner took
    0.323 s and the test reported the machine as a defect while the
    property it exists to protect held perfectly.

    A duration describes the host. The behaviour is asserted
    unconditionally now, and the timing only decides whether the
    STRONGER statement (one hand-over, because the loop never read the
    log) can be made at all.
    """
    text = (pathlib.Path(__file__).resolve().parent
            / "test_gfn_methods_in_the_viewer.py").read_text(encoding="utf-8")
    assert "assert spent < 0.2" not in text, (
        "a wall-clock assertion is back; it fails on a slow runner and "
        "says nothing about the code")
    assert "_READ_INTERVAL_S" in text, (
        "the interval must be named, not written as a bare 0.2 in an "
        "assertion")
