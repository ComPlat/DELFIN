"""Byte-level parity of the orca-mlips builder's plain-text fallback.

The fallback writes the input when OPI is not importable (the py3.11 CI
jobs).  Input: water, every opt/optts/freq variant.  Output: the exact
``.inp`` and ``.xyz`` text, produced by the builder at origin/main
3315360cf.  The OPI path is forced off so the test means the same on a
host with and without OPI.
"""

from __future__ import annotations

import pytest

from delfin.mlp_tools import orca_mlips

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

EXPECTED_PLAIN = {
    'opt': (
        '! Opt ExtOpt\n'
        '\n'
        '%method\n'
        '  ProgExt "orca-mlips-uma"\n'
        '  Ext_Params "--device cpu"\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
    ),
    'optts': (
        '! OptTS ExtOpt\n'
        '\n'
        '%method\n'
        '  ProgExt "orca-mlips-uma"\n'
        '  Ext_Params "--device cpu"\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
    ),
    'freq': (
        '! Freq ExtOpt\n'
        '\n'
        '%method\n'
        '  ProgExt "orca-mlips-uma"\n'
        '  Ext_Params "--device cpu --dump-hessian extopt.hess"\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
    ),
    'optts_two_step': (
        '! OptTS Freq ExtOpt\n'
        '\n'
        '%method\n'
        '  ProgExt "orca-mlips-uma"\n'
        '  Ext_Params "--device cpu --dump-hessian extopt.hess"\n'
        'end\n'
        '\n'
        '%geom\n'
        '  inhess Read\n'
        '  InHessName "extopt.hess"\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
    ),
    'freq_two_step': (
        '! Freq ExtOpt\n'
        '\n'
        '%method\n'
        '  ProgExt "orca-mlips-uma"\n'
        '  Ext_Params "--device cpu --dump-hessian extopt.hess"\n'
        'end\n'
        '\n'
        '%geom\n'
        '  inhess Read\n'
        '  InHessName "extopt.hess"\n'
        'end\n'
        '\n'
        '*xyzfile 0 1 extopt.xyz\n'
    ),
    'opt_mace_solvent': (
        '! Opt ExtOpt\n'
        '\n'
        '%method\n'
        '  ProgExt "orca-mlips-mace"\n'
        '  Ext_Params "--model m1 --device cuda --solvent water --solvent-model alpb"\n'
        'end\n'
        '\n'
        '*xyzfile 1 2 extopt.xyz\n'
    ),
}
EXPECTED_XYZ = (
    '3\n'
    '\n'
    'O 0.0000000000 0.0000000000 0.0000000000\n'
    'H 0.7570000000 0.5860000000 0.0000000000\n'
    'H -0.7570000000 0.5860000000 0.0000000000\n'
)


@pytest.fixture()
def no_opi(monkeypatch):
    def _raise(*_a, **_k):
        raise ImportError("OPI forced off for the fallback parity test")

    monkeypatch.setattr(orca_mlips, "_build_with_opi", _raise)


@pytest.mark.parametrize("name", sorted(CASES))
def test_fallback_inputs_are_byte_identical(name, tmp_path, no_opi):
    run_type, kwargs = CASES[name]
    inp = orca_mlips.build_orca_mlips_input(
        run_type, elements=WATER[0], coordinates=WATER[1], work_dir=tmp_path, **kwargs
    )
    assert inp.read_text() == EXPECTED_PLAIN[name]
    assert (tmp_path / "extopt.xyz").read_text() == "".join(EXPECTED_XYZ)
