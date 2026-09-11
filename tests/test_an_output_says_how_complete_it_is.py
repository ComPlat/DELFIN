"""Asked whether nine energies were reliable, an operator opened every
output by hand to learn they were eight lines each, and found
extract_optimization_trajectory reporting one optimization cycle for a
single point (2026-09-11). The parse result now says how many lines the
output has and which named blocks it holds; the trajectory tool calls a
single point what it is; and the archive fixture's CONTROL names the
molecule its geometry actually is.
"""

from __future__ import annotations

from delfin import api

_STUB = """
                  * O   R   C   A *

FINAL SINGLE POINT ENERGY       -113.302030000000

****ORCA TERMINATED NORMALLY****
"""

_OPT = """
*    Geometry Optimization Run   *
GEOMETRY OPTIMIZATION CYCLE   1
FINAL SINGLE POINT ENERGY       -113.300000000000
GEOMETRY OPTIMIZATION CYCLE   2
FINAL SINGLE POINT ENERGY       -113.302000000000
THE OPTIMIZATION HAS CONVERGED
VIBRATIONAL FREQUENCIES
   0:      0.00 cm**-1
THERMOCHEMISTRY AT 298.15K
Final Gibbs free energy         -113.290000000000 Eh
****ORCA TERMINATED NORMALLY****
"""


def _folder(tmp_path, name, text):
    d = tmp_path / name
    d.mkdir()
    (d / "run.out").write_text(text)
    return d


def test_the_parse_result_counts_lines_and_names_the_blocks(tmp_path):
    stub = api.parse_orca_output(str(_folder(tmp_path, "stub", _STUB)))
    assert stub.output_lines == 6
    assert stub.blocks_present == ["final_energy", "termination"]
    full = api.parse_orca_output(str(_folder(tmp_path, "full", _OPT)))
    assert "geometry_optimization" in full.blocks_present
    assert "frequencies" in full.blocks_present and "thermochemistry" in full.blocks_present
    assert full.blocks_present.index("final_energy") < full.blocks_present.index("termination")


def test_a_single_point_is_not_a_one_cycle_optimization(tmp_path):
    res = api.extract_optimization_trajectory(str(_folder(tmp_path, "sp", _STUB)))
    assert res.cycles == [] and res.n_cycles == 0
    assert res.error and "no geometry optimization" in res.error
    assert "parse_orca_output" in res.error, "the reader is told where the energy is"


def test_a_real_optimization_still_gives_its_cycles(tmp_path):
    res = api.extract_optimization_trajectory(str(_folder(tmp_path, "opt", _OPT)))
    assert res.error is None and res.n_cycles == 2 and res.converged is True
    assert res.final_energy_eh == -113.302


def test_the_tool_descriptions_say_so():
    from delfin.ops_server import server as ops
    assert "blocks_present" in (ops.tool_parse_orca_output.__doc__ or "")
    assert "single point" in (ops.tool_extract_optimization_trajectory.__doc__ or "").lower()


def test_the_archive_fixture_names_the_molecule_its_geometry_is(tmp_path):
    """CONTROL said benzene while every .inp held C=O; an operator read
    the mismatch as a reason to distrust the energies."""
    import subprocess, sys
    from pathlib import Path
    setup = Path(api.__file__).resolve().parent / "agent" / "pack" / "benchmark" / "setup" / "a_small_calc_archive.py"
    subprocess.run([sys.executable, str(setup), str(tmp_path)], check=True,
                   capture_output=True, timeout=120)
    controls = list(tmp_path.rglob("CONTROL.txt"))
    assert controls, "the fixture writes CONTROL files"
    for c in controls:
        assert "c1ccccc1" not in c.read_text()
        assert "[C-]#[O+]" in c.read_text()
