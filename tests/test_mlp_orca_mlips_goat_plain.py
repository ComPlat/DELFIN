"""orca-mlips builder without OPI: GOAT / NEB / IRC inputs from the fallback.

The plain-text fallback (py3.11 CI jobs, no ``opi``) must write the same
keywords and blocks as the OPI path (tests/test_opi_orca_mlips_goat.py).
The OPI path is forced off so the test means the same on every host.
"""

from __future__ import annotations

import pytest

from delfin.mlp_tools import orca_mlips

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])


@pytest.fixture(autouse=True)
def no_opi(monkeypatch):
    def _raise(*_a, **_k):
        raise ImportError("OPI forced off")

    monkeypatch.setattr(orca_mlips, "_build_with_opi", _raise)


def test_goat_fallback_text(tmp_path):
    inp = orca_mlips.build_orca_mlips_input(
        "goat", elements=WATER[0], coordinates=WATER[1], work_dir=tmp_path,
        ncores=8, goat_options={"workers": 8, "energy_window": 6.0},
        base_name="XTB_GOAT",
    )
    assert inp.read_text() == (
        "! GOAT ExtOpt\n"
        "\n"
        "%pal\n"
        "    nprocs 8\n"
        "end\n"
        "\n"
        "%method\n"
        '  ProgExt "orca-mlips-uma"\n'
        '  Ext_Params "--device cpu"\n'
        "end\n"
        "\n"
        "%goat\n"
        "    nworkers 8\n"
        "    maxen 6.0\n"
        "end\n"
        "\n"
        "*xyzfile 0 1 XTB_GOAT.xyz\n"
    )
    assert (tmp_path / "XTB_GOAT.xyz").read_text().splitlines()[0] == "3"


def test_neb_fallback_names_the_end_file(tmp_path):
    inp = orca_mlips.build_orca_mlips_input(
        "neb-ts", elements=WATER[0], coordinates=WATER[1], work_dir=tmp_path,
        neb_end_coordinates=WATER[1],
    )
    text = inp.read_text()
    assert text.startswith("! NEB-TS ExtOpt\n")
    assert '    neb_end_xyzfile "extopt_end.xyz"\n' in text
    assert (tmp_path / "extopt_end.xyz").is_file()


def test_irc_fallback(tmp_path):
    text = orca_mlips.build_orca_mlips_input(
        "irc", elements=WATER[0], coordinates=WATER[1], work_dir=tmp_path,
    ).read_text()
    assert "%irc\n    inithess calc_numfreq\nend\n" in text


def test_scan_fallback_with_pal(tmp_path):
    text = orca_mlips.build_orca_mlips_input(
        "scan", elements=WATER[0], coordinates=WATER[1], work_dir=tmp_path,
        ncores=2, scan={"mode": "B", "atom1": 0, "atom2": 1, "startvalue": 0.9,
                        "endvalue": 1.3, "steps": 5},
    ).read_text()
    assert text.startswith("! Opt ExtOpt\n\n%pal\n    nprocs 2\nend\n\n%method\n")
    assert "  Scan B 0 1 0.9 1.3 5\n" in text
