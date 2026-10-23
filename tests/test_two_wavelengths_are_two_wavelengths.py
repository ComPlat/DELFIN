"""``hyperpol_xTB_wavelengths=1064,532`` means two wavelengths, not a crash.

A CONTROL value holding a comma is a list by the time the workflow reads it --
``_parse_control_file`` splits it. The hyperpolarizability step put that list
through ``str()`` and split the result on commas again, so the manual's own
two-wavelength example produced ``float("['1064'")`` and ended the run with a
ValueError. One wavelength worked, which is why it survived.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from delfin.config import _parse_control_file


def _wavelengths_from(value) -> list[float]:
    """The reader as the workflow runs it (delfin/workflows/pipeline.py)."""
    parts = ([str(w) for w in value] if isinstance(value, (list, tuple))
             else str(value).split(','))
    parts = [w.strip() for w in parts if str(w).strip()]
    if parts and ' '.join(parts).lower() not in ('none', 'static'):
        return [float(w) for w in parts]
    return [math.inf]


def _control(tmp_path: Path, line: str):
    path = tmp_path / "CONTROL.txt"
    path.write_text(line + "\n", encoding="utf-8")
    return _parse_control_file(str(path), keep_steps_literal=False)


def test_a_comma_separated_value_arrives_as_a_list(tmp_path):
    # the premise of the bug: this is what the parser hands the workflow
    config = _control(tmp_path, "hyperpol_xTB_wavelengths=1064,532")
    assert config["hyperpol_xTB_wavelengths"] == ["1064", "532"]


def test_two_wavelengths_are_read_as_two_numbers(tmp_path):
    config = _control(tmp_path, "hyperpol_xTB_wavelengths=1064,532")
    assert _wavelengths_from(config["hyperpol_xTB_wavelengths"]) == [1064.0, 532.0]


def test_one_wavelength_still_works(tmp_path):
    config = _control(tmp_path, "hyperpol_xTB_wavelengths=1064")
    assert _wavelengths_from(config["hyperpol_xTB_wavelengths"]) == [1064.0]


@pytest.mark.parametrize("line", ["hyperpol_xTB_wavelengths=",
                                  "hyperpol_xTB_wavelengths=none",
                                  "hyperpol_xTB_wavelengths=static"])
def test_no_wavelength_means_the_static_calculation(tmp_path, line):
    config = _control(tmp_path, line)
    assert _wavelengths_from(config.get("hyperpol_xTB_wavelengths", "")) == [math.inf]


def test_the_workflow_reads_it_the_same_way():
    """The helper above must stay the code under test, not a copy of it."""
    source = Path(__file__).resolve().parents[1] / "delfin" / "workflows" / "pipeline.py"
    text = source.read_text(encoding="utf-8")
    assert "isinstance(raw_wl, (list, tuple))" in text
    assert "float(w) for w in parts" in text
