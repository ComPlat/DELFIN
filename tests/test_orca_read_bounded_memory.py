"""Every quantity read through :func:`delfin.orca_read.read` stays bounded.

``read()`` used to compute *everything* for a single requested quantity:
the energy finders, ``extract_thermochem`` and ``parse_orca_output`` each
call ``f.read()`` on the whole output, so one 12 MB ``.out`` peaked at
~200 MB and a realistic 100+ MB run loads the file many times over.

These tests pin that each quantity reads only what it needs from a bounded
region of the output (a tail, the input file, or its own section), measured
with tracemalloc on a generated multi-MB output. Same 8 MB budget as
``test_recovery_answers_the_error_orca_actually_reported.py``'s
``test_the_output_is_not_read_into_memory_whole``.
"""
from pathlib import Path

import pytest

from delfin.orca_read import Quantity, read

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir

FIXTURES = lc1_fixture_dir("tests/fixtures/orca611")


def _fixture_dir(name: str):
    """The specific orca611 subfolder: legacy path while its files are
    committed, else the private store (skip without the env var)."""
    return lc1_fixture_dir(f"tests/fixtures/orca611/{name}")

#: Same budget the recovery analyser's tail scan is held to.
PEAK_BUDGET = 8 * 1024 * 1024

#: Tail region of the output that carries termination, the final single-point
#: energy, the thermochemistry block and the geometry/outcome markers.
_TAIL = 2 * 1024 * 1024


def _generated_multimb_out(tmp_path: Path, source_name: str) -> Path:
    """A ~12 MB ``.out``: the fixture's own text with irrelevant dense padding
    inserted mid-file, keeping the head and the true tail (with the markers
    each quantity needs) intact."""
    src = _fixture_dir(source_name) / f"{source_name}.out"
    src_inp = _fixture_dir(source_name) / f"{source_name}.inp"
    out = tmp_path / f"{source_name}.out"
    inp = tmp_path / f"{source_name}.inp"
    text = src.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines(keepends=True)
    filler = ("x" * 4096 + "\n") * (12 * 1024 * 1024 // 4097)  # ~12 MB
    cut = len(lines) // 2
    out.write_text("".join(lines[:cut]) + filler + "".join(lines[cut:]),
                   encoding="utf-8")
    inp.write_text(src_inp.read_text(encoding="utf-8", errors="replace"),
                   encoding="utf-8")
    assert out.stat().st_size > 8 * 1024 * 1024
    return out


def _measure_peak(work_dir: Path, name: str, quantity: Quantity) -> float:
    import tracemalloc

    tracemalloc.start()
    read(work_dir, name, quantity)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak


@pytest.mark.parametrize("quantity", [
    Quantity.TERMINATED_NORMALLY,
    Quantity.FINAL_SINGLE_POINT_ENERGY,
])
def test_termination_and_final_energy_read_only_a_bounded_tail(tmp_path, quantity):
    """The two expensive tail-lived quantities must not load the whole file."""
    _generated_multimb_out(tmp_path, "h2o_optfreq")
    peak = _measure_peak(tmp_path, "h2o_optfreq", quantity)
    assert peak < PEAK_BUDGET, (
        f"read {peak / 1024 / 1024:.1f} MB for {quantity.name} on a "
        f"multi-MB output -- the whole file is being loaded"
    )


@pytest.mark.parametrize("quantity", [
    Quantity.GIBBS,
    Quantity.ZPE,
    Quantity.ENTHALPY,
    Quantity.INNER_ENERGY,
    Quantity.ENTROPY,
])
def test_thermochemistry_quantities_read_only_their_section(tmp_path, quantity):
    """G/ZPE/H/U/S all sit in the thermochemistry block near the end."""
    _generated_multimb_out(tmp_path, "h2o_optfreq")
    peak = _measure_peak(tmp_path, "h2o_optfreq", quantity)
    assert peak < PEAK_BUDGET, (
        f"read {peak / 1024 / 1024:.1f} MB for {quantity.name} -- "
        "the whole output is being loaded instead of the thermochemistry section"
    )


@pytest.mark.parametrize("quantity", [
    Quantity.CHARGE,
    Quantity.MULTIPLICITY,
    Quantity.NUMBER_OF_ELECTRONS,
])
def test_input_side_quantities_never_read_the_output(tmp_path, quantity):
    """Charge / multiplicity / electron count come from the .inp, never the .out."""
    _generated_multimb_out(tmp_path, "h2o_optfreq")
    peak = _measure_peak(tmp_path, "h2o_optfreq", quantity)
    assert peak < PEAK_BUDGET, (
        f"read {peak / 1024 / 1024:.1f} MB for {quantity.name} -- "
        "an input-side quantity must not touch the output"
    )


@pytest.mark.parametrize("quantity", [
    Quantity.GEOMETRY_CONVERGED,
    Quantity.FAILURE_REASON,
    Quantity.S2_DEVIATION,
])
def test_geometry_outcome_and_spin_read_only_a_bounded_tail(tmp_path, quantity):
    """Geometry convergence, the failure outcome and the UHF block all live
    near the end; a full-file parse is not required for one of them."""
    _generated_multimb_out(tmp_path, "h2o_optfreq")
    peak = _measure_peak(tmp_path, "h2o_optfreq", quantity)
    assert peak < PEAK_BUDGET, (
        f"read {peak / 1024 / 1024:.1f} MB for {quantity.name} -- "
        "only a bounded tail is needed"
    )
