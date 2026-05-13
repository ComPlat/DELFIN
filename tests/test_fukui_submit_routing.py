"""Test that DELFIN_MODE=fukui in local_runner dispatches to cli_fukui.main()."""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest


def _set_env(monkeypatch, **kwargs):
    """Set env vars common to all fukui-mode dispatch tests."""
    monkeypatch.setenv('DELFIN_PAL', '4')
    monkeypatch.setenv('DELFIN_MAXCORE', '3000')
    for key in (
        'DELFIN_FUKUI_INPUT', 'DELFIN_FUKUI_SCHEME', 'DELFIN_FUKUI_PRE_OPT',
        'DELFIN_FUKUI_SKIP_CUBES', 'DELFIN_FUKUI_FUNCTIONAL',
        'DELFIN_FUKUI_BASIS', 'DELFIN_FUKUI_DISPERSION',
        'DELFIN_FUKUI_SOLVATION', 'DELFIN_FUKUI_SOLVENT',
    ):
        monkeypatch.delenv(key, raising=False)
    for key, value in kwargs.items():
        monkeypatch.setenv(key, value)


def test_fukui_mode_missing_input_returns_error(monkeypatch, tmp_path):
    from delfin.dashboard import local_runner

    _set_env(monkeypatch)
    monkeypatch.chdir(tmp_path)
    rc = local_runner._run_mode('fukui')
    assert rc == 1


def test_fukui_mode_dispatches_to_cli_with_smiles(monkeypatch, tmp_path):
    from delfin.dashboard import local_runner

    _set_env(
        monkeypatch,
        DELFIN_FUKUI_INPUT='C=O',
        DELFIN_FUKUI_SCHEME='loewdin',
        DELFIN_FUKUI_SKIP_CUBES='1',
    )
    monkeypatch.chdir(tmp_path)

    with mock.patch('delfin.cli_fukui.main', return_value=0) as fukui_main:
        rc = local_runner._run_mode('fukui')
        assert rc == 0
        args = fukui_main.call_args[0][0]
        assert '--input' in args and 'C=O' in args
        assert '--scheme' in args and 'loewdin' in args
        assert '--workdir' in args
        assert '--skip-cubes' in args
        assert '--pal' in args and '4' in args


def test_fukui_mode_pre_opt_flag_yes(monkeypatch, tmp_path):
    from delfin.dashboard import local_runner

    _set_env(
        monkeypatch,
        DELFIN_FUKUI_INPUT='C=O',
        DELFIN_FUKUI_PRE_OPT='yes',
    )
    monkeypatch.chdir(tmp_path)

    with mock.patch('delfin.cli_fukui.main', return_value=0) as fukui_main:
        local_runner._run_mode('fukui')
        args = fukui_main.call_args[0][0]
        assert '--pre-opt' in args
        assert '--no-pre-opt' not in args


def test_fukui_mode_pre_opt_flag_no(monkeypatch, tmp_path):
    from delfin.dashboard import local_runner

    _set_env(
        monkeypatch,
        DELFIN_FUKUI_INPUT='/tmp/mol.xyz',
        DELFIN_FUKUI_PRE_OPT='0',
    )
    monkeypatch.chdir(tmp_path)

    with mock.patch('delfin.cli_fukui.main', return_value=0) as fukui_main:
        local_runner._run_mode('fukui')
        args = fukui_main.call_args[0][0]
        assert '--no-pre-opt' in args
        assert '--pre-opt' not in args


def test_fukui_mode_forwards_orca_settings(monkeypatch, tmp_path):
    from delfin.dashboard import local_runner

    _set_env(
        monkeypatch,
        DELFIN_FUKUI_INPUT='mol.xyz',
        DELFIN_FUKUI_FUNCTIONAL='PBE0',
        DELFIN_FUKUI_BASIS='def2-TZVP',
        DELFIN_FUKUI_DISPERSION='D4',
        DELFIN_FUKUI_SOLVATION='CPCM',
        DELFIN_FUKUI_SOLVENT='water',
    )
    monkeypatch.chdir(tmp_path)

    with mock.patch('delfin.cli_fukui.main', return_value=0) as fukui_main:
        local_runner._run_mode('fukui')
        args = fukui_main.call_args[0][0]
        # Verify each setting is forwarded as a flag/value pair
        for flag, value in [
            ('--functional', 'PBE0'),
            ('--basis', 'def2-TZVP'),
            ('--dispersion', 'D4'),
            ('--solvation', 'CPCM'),
            ('--solvent', 'water'),
        ]:
            idx = args.index(flag)
            assert args[idx + 1] == value


def test_pyproject_has_delfin_fukui_entry():
    """Sanity-check the CLI entry point is registered."""
    pyproject = Path(__file__).resolve().parent.parent / 'pyproject.toml'
    text = pyproject.read_text(encoding='utf-8')
    assert 'delfin-fukui = "delfin.cli_fukui:main"' in text
