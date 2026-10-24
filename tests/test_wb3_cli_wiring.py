"""WB3: cli_new_adapter wiring test (phase 2).

Drives the real handler the CLI hook calls (cli.py once the operator
applies .gate/wb3_cli_new_adapter.patch): build_parser + handler
end-to-end on a throwaway tree, with HOME redirected so nothing touches
the real home.  The registry case routes the freshly written module
through the same user-adapter loader production uses
(delfin/tools/_registry.load_user_adapter_file), so it tests the
registration path the scaffold relies on without monkeypatching sys.path.
"""
from __future__ import annotations

import argparse

import delfin.tools._registry as reg
from delfin.agent import cli_new_adapter
from delfin.tools._registry import list_steps


def argparse_Sub():
    return argparse.ArgumentParser().add_subparsers(dest="cmd", required=True)


def _parse(*argv):
    parser = cli_new_adapter.build_parser(argparse_Sub())
    return parser.parse_args(list(argv))


def test_handler_scaffolds_into_throwaway_tree(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(tmp_path / "user_adapters"))
    args = _parse("wb3_demo", "--kind", "parser", "--dir", str(tmp_path))
    rc = cli_new_adapter.handler(args)
    assert rc == 0
    assert (tmp_path / "delfin" / "tools" / "adapters" / "wb3_demo.py").is_file()
    assert (tmp_path / "tests" / "test_adapter_wb3_demo.py").is_file()
    assert (tmp_path / "docs" / "extending_delfin.md").is_file()


def test_scaffolded_module_registers_through_user_adapter_loader(tmp_path, monkeypatch):
    """Same loader production discovery uses — no hand edit, no sys.path game."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(tmp_path / "user_adapters"))
    args = _parse("wb3_demo", "--kind", "parser", "--dir", str(tmp_path))
    assert cli_new_adapter.handler(args) == 0

    reg._DISCOVERED = False
    reg._REGISTRY.clear()
    skeleton = tmp_path / "delfin" / "tools" / "adapters" / "wb3_demo.py"
    reg.load_user_adapter_file(skeleton)
    names = set(list_steps())
    assert "wb3_demo_parser" in names, (
        f"scaffolded adapter missing after load (saw {len(names)} steps)"
    )


def test_handler_refuses_existing_file(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    args = _parse("wb3_demo", "--dir", str(tmp_path))
    assert cli_new_adapter.handler(args) == 0
    rc = cli_new_adapter.handler(args)
    assert rc == 2, "second run must refuse to overwrite (exit 2)"
