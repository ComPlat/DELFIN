"""orca-mlips builder: GOAT, NEB, NEB-TS, IRC and scan inputs through OPI.

Manual ground truth (ORCA 6.1): 4.1.10 "Optimizing with External Methods"
-- ``!ExtOpt`` "must be combined with other keywords like !Opt or !GOAT";
4.10 GOAT -- "The command to call the global optimizer is simply !GOAT
... its options can be given under the %GOAT block"; Table 4.8 (GOAT
simple keywords), Table 4.9 (%goat options; NWORKERS "any number multiple
of 4"); neb.txt (NEB_End_XYZFile); irc.txt (InitHess read |
calc_numfreq, Hess_Filename).

Measured once with ORCA 6.1.1 (parser check, not part of this suite): a
GOAT input written by this builder, with ProgExt pointed at a stub that
exits at the first call, was read without a refusal, entered GOAT and
called the stub as ``<stub> <base>_EXT.extinp.tmp --backend uma --device
cpu`` with NCores from %pal.

py3.12 OPI venv only (skipped without ``opi``).
"""

from __future__ import annotations

import re

import pytest

pytest.importorskip("opi")

from delfin.mlp_tools.orca_mlips import (  # noqa: E402
    GOAT_VARIANTS,
    build_orca_mlips_input,
    goat_output_files,
    normalize_goat_options,
)

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])
WATER_END = [[0.0, 0.0, 0.1], [0.757, 0.586, 0.1], [-0.757, 0.586, 0.1]]


def _goat(tmp_path, **kw):
    return build_orca_mlips_input(
        "goat", elements=WATER[0], coordinates=WATER[1], work_dir=tmp_path, **kw
    ).read_text()


def _block(text: str, name: str) -> dict[str, str]:
    m = re.search(rf"(?ims)^%{name}\s*\n(.*?)^end", text)
    assert m, f"no %{name} block in:\n{text}"
    return dict(line.split(None, 1) for line in m.group(1).strip().splitlines())


def test_goat_input_has_goat_extopt_and_progext(tmp_path):
    text = _goat(tmp_path)
    low = text.lower()
    assert "!goat" in low and "!extopt" in low
    assert 'progext "orca-mlips-uma"' in low
    assert 'Ext_Params "--backend uma --device cpu"' in text
    assert "*xyzfile 0 1 extopt.xyz" in text
    assert "%goat" not in low  # nothing given -> every option at ORCA's default


def test_goat_options_are_written_with_aliases(tmp_path):
    text = _goat(tmp_path, goat_options={
        "workers": 8, "energy_window": 6.0, "temperature": 298.15,
        "KEEPWORKERDATA": True, "maxoptiter": 64,
    })
    block = _block(text, "goat")
    assert block == {"nworkers": "8", "maxen": "6.0", "conftemp": "298.15",
                     "keepdata": "True", "maxoptiter": "64"}


def test_ncores_becomes_pal_nprocs(tmp_path):
    text = _goat(tmp_path, ncores=8)
    assert re.search(r"(?m)^%pal\s*\n\s*nprocs 8\s*\nend", text)


@pytest.mark.parametrize("variant", GOAT_VARIANTS)
def test_goat_variants(variant, tmp_path):
    assert f"!{variant.lower()}" in _goat(tmp_path, goat_variant=variant).lower()


def test_unknown_goat_option_refused_before_any_file(tmp_path):
    with pytest.raises(ValueError, match="not one ORCA knows"):
        _goat(tmp_path, goat_options={"definitelynotanoption": 1})
    assert not (tmp_path / "extopt.inp").exists()


@pytest.mark.parametrize("bad", [6, 0, -4, "8", True])
def test_nworkers_must_be_a_multiple_of_four(bad):
    with pytest.raises(ValueError, match="multiple of 4"):
        normalize_goat_options({"nworkers": bad})


def test_value_injection_refused():
    with pytest.raises(ValueError, match="value"):
        normalize_goat_options({"maxen": "6\nend\n%pal nprocs 64"})


def test_goat_variant_coarse_refused(tmp_path):
    with pytest.raises(ValueError, match="goat_variant"):
        _goat(tmp_path, goat_variant="GOAT-COARSE")


def test_goat_options_only_for_goat(tmp_path):
    with pytest.raises(ValueError, match="goat"):
        build_orca_mlips_input("opt", elements=WATER[0], coordinates=WATER[1],
                               work_dir=tmp_path, goat_options={"maxen": 6.0})


def test_two_step_hessian_refused_for_goat(tmp_path):
    with pytest.raises(ValueError, match="two_step"):
        _goat(tmp_path, hessian="two_step")


def test_goat_output_file_names():
    files = goat_output_files("XTB_GOAT")
    assert files["finalensemble"] == "XTB_GOAT.finalensemble.xyz"
    assert files["global_minimum"] == "XTB_GOAT.globalminimum.xyz"
    assert files["output"] == "XTB_GOAT.out"


# ------------------------------------------------------------- NEB / IRC / scan


@pytest.mark.parametrize("run_type,kw", [("neb", "!neb"), ("neb-ts", "!neb-ts")])
def test_neb_writes_end_geometry(run_type, kw, tmp_path):
    text = build_orca_mlips_input(
        run_type, elements=WATER[0], coordinates=WATER[1], work_dir=tmp_path,
        neb_end_coordinates=WATER_END,
    ).read_text()
    assert kw in text.lower()
    assert _block(text, "neb") == {"neb_end_xyzfile": '"extopt_end.xyz"'}
    end = (tmp_path / "extopt_end.xyz").read_text().splitlines()
    assert end[0] == "3" and end[2].split()[3] == "0.1000000000"


def test_neb_without_end_geometry_refused(tmp_path):
    with pytest.raises(ValueError, match="neb_end_coordinates"):
        build_orca_mlips_input("neb", elements=WATER[0], coordinates=WATER[1],
                               work_dir=tmp_path)


def test_irc_numerical_initial_hessian(tmp_path):
    text = build_orca_mlips_input("irc", elements=WATER[0], coordinates=WATER[1],
                                  work_dir=tmp_path).read_text()
    assert "!irc" in text.lower()
    assert _block(text, "irc") == {"inithess": "calc_numfreq"}


def test_irc_reads_the_dumped_hessian(tmp_path):
    text = build_orca_mlips_input("irc", elements=WATER[0], coordinates=WATER[1],
                                  work_dir=tmp_path, hessian="two_step").read_text()
    assert _block(text, "irc") == {"inithess": "read",
                                   "hess_filename": '"extopt.hess"'}
    assert "--dump-hessian extopt.hess" in text


def test_scan_writes_geom_scan(tmp_path):
    text = build_orca_mlips_input(
        "scan", elements=WATER[0], coordinates=WATER[1], work_dir=tmp_path,
        scan={"mode": "B", "atom1": 0, "atom2": 1, "startvalue": 0.9,
              "endvalue": 1.3, "steps": 5},
    ).read_text()
    assert "!opt" in text.lower()
    assert "scan" in text.lower()


def test_scan_without_dict_refused(tmp_path):
    with pytest.raises(ValueError, match="scan"):
        build_orca_mlips_input("scan", elements=WATER[0], coordinates=WATER[1],
                               work_dir=tmp_path)
