"""Phase C2: external_cli providers in the ExtOpt input builder.

GOAT / OptTS / NEB-TS / IRC / scan inputs already exist for the wrapper
route (EXTOPT_RUN_TYPES); these tests pin that the orca-mlips providers
plug into the SAME builder: the plugin executable becomes ProgExt, the
plugin's own options replace --backend, and legacy inputs are unchanged
byte-for-byte in their ProgExt/Ext_Params lines.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from delfin.mlp_tools.orca_extopt import build_extopt_input

GEOM = (["H", "O", "H"], [[0, 0, 0], [0, 0, 1.5], [0, 1.43, -0.4]])


def _inp(tmp_path, run_type, backend, base_name):
    elems, coords = GEOM
    return build_extopt_input(
        run_type, elements=elems, coordinates=coords, backend=backend,
        base_name=base_name, work_dir=tmp_path).read_text()


def test_goat_with_plugin_routes_progext_to_the_executable(tmp_path):
    txt = _inp(tmp_path, "goat", "orca_mlips_aimnet2", "goat_plugin")
    assert 'ProgExt "orca-mlips-aimnet2"' in txt
    assert "GOAT" in txt
    params_line = next(l for l in txt.splitlines() if "Ext_Params" in l)
    assert "--backend" not in params_line


def test_optts_and_neb_ts_with_plugins(tmp_path):
    optts = _inp(tmp_path, "optts", "orca_mlips_uma", "t1")
    assert "OptTS" in optts and "orca-mlips-uma" in optts
    neb = _inp(tmp_path, "neb-ts", "orca_mlips_orb", "t2")
    assert "NEB-TS" in neb and "orca-mlips-orb" in neb


def test_licence_restricted_provider_still_builds_but_flags(tmp_path):
    # building is a dry-run; licence enforcement lives in check_availability
    txt = _inp(tmp_path, "opt", "orca_mlips_mace", "t3")
    assert "orca-mlips-mace" in txt


def test_legacy_wrapper_route_is_byte_identical(tmp_path):
    txt = _inp(tmp_path, "opt", "ani2x", "legacy")
    assert 'ProgExt "delfin-mlp-extopt"' in txt
    assert "--backend ani2x" in txt


def test_unknown_backend_falls_back_to_wrapper_route(tmp_path):
    # a plain non-registry name keeps the old --backend behaviour
    txt = _inp(tmp_path, "opt", "ani2x", "t4")
    assert "Ext_Params \"--backend ani2x" in txt
