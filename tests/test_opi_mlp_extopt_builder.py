"""OPI-based tests for build_extopt_input (delfin/mlp_tools/orca_extopt.py).

Runs in the py3.12 OPI venv (CI sets DELFIN_OPI_REQUIRED=1); the input is
built through OPI (opi.core.Calculator, BlockMethod, Scan) and compared
against the ORCA 6.1.1 manual's external-tools input, chapter 4.1.10
(structurereactivity/optimizations.txt lines 1537-1701).  No ORCA run.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from delfin.mlp_tools.orca_extopt import (
    EXTOPT_RUN_TYPES,
    build_extopt_input,
)

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])


@pytest.fixture()
def work_dir(tmp_path):
    return tmp_path


def _read(work_dir, base_name="extopt"):
    return (Path(work_dir) / f"{base_name}.inp").read_text()


def test_run_types_complete():
    assert {"opt", "optts", "neb-ts", "neb", "irc", "goat", "scan"} <= set(EXTOPT_RUN_TYPES)


def test_opt_input_contains_extopt_and_progext(work_dir):
    elements, coords = WATER
    path = build_extopt_input(
        "opt", elements=elements, coordinates=coords, work_dir=work_dir
    )
    text = _read(work_dir)
    assert path.name == "extopt.inp"
    assert "extopt" in text.lower()
    assert "opt" in text.lower()
    # ProgExt in the %method block.  OPI writes the path unquoted via
    # InputFilePath; the manual's example quotes it, but quotes are only
    # required for paths containing whitespace.
    assert "%method" in text.lower()
    assert "progext delfin-mlp-extopt" in text.lower()
    assert "ext_params" in text.lower()


def test_progext_receives_absolute_wrapper_path(work_dir):
    elements, coords = WATER
    build_extopt_input(
        "opt", elements=elements, coordinates=coords,
        wrapper_path="/abs/wrapper", work_dir=work_dir,
    )
    text = _read(work_dir)
    assert "/abs/wrapper" in text


def test_ext_params_forwards_backend(work_dir):
    elements, coords = WATER
    build_extopt_input(
        "opt", elements=elements, coordinates=coords,
        backend="mace_off", work_dir=work_dir,
    )
    text = _read(work_dir)
    assert "--backend mace_off" in text


def test_neb_ts_irc_goat_keywords(work_dir):
    elements, coords = WATER
    for run_type, expected in (
        ("neb-ts", ["neb-ts", "extopt"]),
        ("irc", ["irc", "neb-ts", "extopt"]),
        ("goat", ["goat", "extopt"]),
    ):
        work = Path(work_dir) / run_type
        work.mkdir()
        build_extopt_input(
            run_type, elements=elements, coordinates=coords, work_dir=work
        )
        text = _read(work)
        for kw in expected:
            assert kw in text.lower(), f"{run_type}: {kw} missing"


def test_scan_block(work_dir):
    elements, coords = WATER
    build_extopt_input(
        "scan", elements=elements, coordinates=coords,
        scan={"mode": "B", "atom1": 1, "atom2": 2,
              "startvalue": 0.9, "endvalue": 1.4, "steps": 10},
        work_dir=work_dir,
    )
    text = _read(work_dir)
    assert "%geom" in text.lower()
    assert "scan" in text.lower()
    assert "0.9" in text and "1.4" in text and "10" in text


def test_scan_requires_scan_dict(work_dir):
    elements, coords = WATER
    with pytest.raises(ValueError, match="scan"):
        build_extopt_input("scan", elements=elements, coordinates=coords,
                           work_dir=work_dir)


def test_unknown_run_type_rejected(work_dir):
    elements, coords = WATER
    with pytest.raises(ValueError, match="run_type"):
        build_extopt_input("freq", elements=elements, coordinates=coords,
                           work_dir=work_dir)


def test_geometry_block(work_dir):
    elements, coords = WATER
    build_extopt_input(
        "opt", elements=elements, coordinates=coords,
        charge=0, mult=1, work_dir=work_dir,
    )
    text = _read(work_dir)
    assert "*xyzfile" in text.lower()
    assert "extopt.xyz" in text
