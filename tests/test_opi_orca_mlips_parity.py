"""Byte-level parity of the orca-mlips builder's opt/optts/freq inputs (OPI path).

Input: water, every opt/optts/freq variant the builder had before the
GOAT/NEB/IRC/scan run types were added.  Output: the exact ``.inp`` text
OPI writes.  The expected texts were produced by the builder at
origin/main 3315360cf, with one later change: OPI writes Ext_Params
quoted (tests/test_opi_ext_params_is_quoted.py).  A change to them is a
change of the ORCA input every orca-mlips job receives and must be made
here on purpose.

Runs in the py3.12 OPI venv only (skipped without ``opi``).
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.mlp_tools.orca_mlips import build_orca_mlips_input  # noqa: E402

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])

CASES = {
    "opt": ("opt", {}),
    "optts": ("optts", {}),
    "freq": ("freq", {}),
    "optts_two_step": ("optts", {"hessian": "two_step"}),
    "freq_two_step": ("freq", {"hessian": "two_step"}),
    "opt_mace_solvent": ("opt", {"solvent": "water", "model": "m1", "device": "cuda",
                                 "backend": "mace", "charge": 1, "mult": 2}),
}

EXPECTED_OPI = {
    'opt': (
        '!opt\n'
        '!extopt\n'
        '%maxcore 1000\n'
        '%pal\n'
        '    nprocs 1\n'
        'end\n'
        '\n'
        '%Method\n'
        '    ProgExt "orca-mlips-uma"\n'
        '    Ext_Params "--backend uma --device cpu"\n'
        'end\n'
        '\n'
        '%output\n'
        '    jsonpropfile True\n'
        '    jsongbwfile True\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
    ),
    'optts': (
        '!optts\n'
        '!extopt\n'
        '%maxcore 1000\n'
        '%pal\n'
        '    nprocs 1\n'
        'end\n'
        '\n'
        '%Method\n'
        '    ProgExt "orca-mlips-uma"\n'
        '    Ext_Params "--backend uma --device cpu"\n'
        'end\n'
        '\n'
        '%output\n'
        '    jsonpropfile True\n'
        '    jsongbwfile True\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
    ),
    'freq': (
        '!freq\n'
        '!extopt\n'
        '%maxcore 1000\n'
        '%pal\n'
        '    nprocs 1\n'
        'end\n'
        '\n'
        '%Method\n'
        '    ProgExt "orca-mlips-uma"\n'
        '    Ext_Params "--backend uma --device cpu --dump-hessian extopt.hess"\n'
        'end\n'
        '\n'
        '%output\n'
        '    jsonpropfile True\n'
        '    jsongbwfile True\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
    ),
    'optts_two_step': (
        '!optts\n'
        '!freq\n'
        '!extopt\n'
        '%maxcore 1000\n'
        '%pal\n'
        '    nprocs 1\n'
        'end\n'
        '\n'
        '%Method\n'
        '    ProgExt "orca-mlips-uma"\n'
        '    Ext_Params "--backend uma --device cpu --dump-hessian extopt.hess"\n'
        'end\n'
        '\n'
        '%output\n'
        '    jsonpropfile True\n'
        '    jsongbwfile True\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
        '\n'
        '%geom\n'
        '  inhess Read\n'
        '  InHessName "extopt.hess"\n'
        'end\n'
    ),
    'freq_two_step': (
        '!freq\n'
        '!extopt\n'
        '%maxcore 1000\n'
        '%pal\n'
        '    nprocs 1\n'
        'end\n'
        '\n'
        '%Method\n'
        '    ProgExt "orca-mlips-uma"\n'
        '    Ext_Params "--backend uma --device cpu --dump-hessian extopt.hess"\n'
        'end\n'
        '\n'
        '%output\n'
        '    jsonpropfile True\n'
        '    jsongbwfile True\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
        '\n'
        '%geom\n'
        '  inhess Read\n'
        '  InHessName "extopt.hess"\n'
        'end\n'
    ),
    'opt_mace_solvent': (
        '!opt\n'
        '!extopt\n'
        '%maxcore 1000\n'
        '%pal\n'
        '    nprocs 1\n'
        'end\n'
        '\n'
        '%Method\n'
        '    ProgExt "orca-mlips-mace"\n'
        '    Ext_Params "--backend mace --model m1 --device cuda --solvent water --solvent-model alpb"\n'
        'end\n'
        '\n'
        '%output\n'
        '    jsonpropfile True\n'
        '    jsongbwfile True\n'
        'end\n'
        '\n'
        '*xyzfile 1 2 extopt.xyz\n'
    ),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_opt_optts_freq_inputs_are_byte_identical(name, tmp_path):
    run_type, kwargs = CASES[name]
    inp = build_orca_mlips_input(
        run_type, elements=WATER[0], coordinates=WATER[1], work_dir=tmp_path, **kwargs
    )
    assert inp.read_text() == EXPECTED_OPI[name]
