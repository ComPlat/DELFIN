"""The rotation-movie glue: ORCA .out -> multi-frame xyz.

A dihedral scan, an optimisation or an IRC each print one
``CARTESIAN COORDINATES (ANGSTROEM)`` block per geometry. The glue
collects every block into one multi-frame xyz so the agent's native
``show_molecule`` renders the rotation movie. No new renderer.

No OPI dependency; runs under the py3.11 CI interpreter.
"""
import pytest

from delfin.tools.geometry_edit import (
    GeometryError,
    parse_xyz,
    trajectory_xyz,
)

SYN = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"

# two-frame fake output: anti (phi=180) then gauche (phi=60). Hand-built
# arithmetic, no CSD data.
_AT = ("C   0.000000000    0.000000000    0.000000000\n"
       "C   0.000000000    0.000000000    1.526000000\n"
       "C   1.415000000    0.000000000    2.034000000\n"
       "C   2.830000000    0.000000000    1.012000000\n")
_GA = ("C   0.000000000    0.000000000    0.000000000\n"
       "C   0.000000000    0.000000000    1.526000000\n"
       "C   1.415000000    0.700000000    2.034000000\n"
       "C   2.830000000    1.400000000    1.012000000\n")


def _out(n_blocks: int = 2) -> str:
    parts = [SYN, "! r2SCAN-3c TightOpt TightSCF\n"]
    for i, body in enumerate([_AT, _GA][:n_blocks]):
        parts.append(f"CARTESIAN COORDINATES (ANGSTROEM)\n"
                     f"  0  1\n  0  1\n{body}")
    return "".join(parts)


def test_collects_every_cartesian_block_as_frames(tmp_path):
    out = tmp_path / "scan.out"
    out.write_text(_out())
    xyz = trajectory_xyz(str(out))
    # multi-frame: the file holds n_blocks count+comment headers
    assert xyz.count("frame") == 2
    g = parse_xyz("\n".join(xyz.splitlines()[:6]))
    assert g.symbols == ["C", "C", "C", "C"]


def test_returns_xyz_string_that_show_molecule_can_open(tmp_path):
    out = tmp_path / "scan.out"
    out.write_text(_out())
    xyz = trajectory_xyz(str(out))
    # multi-frame xyz: alternating count line and comment lines
    lines = xyz.splitlines()
    assert lines[0] == "4"
    assert lines[1].startswith("frame")
    assert lines[6] == "4"
    # second frame holds the gauche coords
    assert lines[10].split()[2].startswith("0.700")


def test_refuses_output_without_geometry_blocks(tmp_path):
    out = tmp_path / "sp.out"
    out.write_text(SYN + "FINAL SINGLE POINT ENERGY  -157.41\n")
    with pytest.raises(GeometryError):
        trajectory_xyz(str(out))


# --- agent-facing glue (api + ops tool + catalog) ---------------------------


def test_api_trajectory_xyz_returns_dict(tmp_path):
    from delfin import api as delfin_api
    d = tmp_path / "scan.out"
    d.write_text(_out())
    res = delfin_api.trajectory_xyz(str(d))
    assert res["n_frames"] == 2
    assert res["xyz"].count("frame") == 2


def test_ops_wrapper_and_catalog_list_trajectory_xyz(tmp_path):
    from delfin.ops_server import server as srv
    import json
    d = tmp_path / "scan.out"
    d.write_text(_out())
    out = json.loads(srv.tool_trajectory_xyz(str(d)))
    assert out["n_frames"] == 2
    catalog = delfin_catalog_names()
    assert "trajectory_xyz" in catalog


def delfin_catalog_names():
    from delfin import api as delfin_api
    return {e["name"] for e in delfin_api._TOOL_CATALOG}
