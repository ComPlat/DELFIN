"""WB3: the user-adapter destination must produce a file discovery loads.

Red shape (s19 F2a): scaffolding with --dir pointing at a USER adapter
directory (flat *.py — what _load_user_adapters reads) must put the
skeleton at the top level there, so that a real _load_user_adapters()
call — the path production takes, not a direct load_user_adapter_file —
registers the adapter.  Also covers F2b: a SECOND scaffold in the same
tree must succeed (docs entry appended, per-adapter files still refuse
overwrite).
"""
from __future__ import annotations

import argparse

import delfin.tools._registry as reg
from delfin.agent import cli_new_adapter
from delfin.tools._registry import list_steps


def argparse_Sub():
    return argparse.ArgumentParser().add_subparsers(dest="cmd", required=True)


def test_user_adapter_dir_scaffold_is_discovered_by_real_loader(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    user_dir = tmp_path / "adapters"
    user_dir.mkdir()
    # a sibling project file marks this as NOT a repo checkout; the user
    # adapter directory itself stays empty so the flat-layout rule applies
    (tmp_path / "project.txt").write_text("not a repo\n", encoding="utf-8")
    parser = cli_new_adapter.build_parser(argparse_Sub())
    args = parser.parse_args(["wb3_user_demo", "--kind", "parser",
                              "--dir", str(user_dir)])
    assert cli_new_adapter.handler(args) == 0
    assert (user_dir / "wb3_user_demo.py").is_file(), (
        "skeleton must be FLAT in a user adapter directory"
    )

    reg._DISCOVERED = False
    reg._REGISTRY.clear()
    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(user_dir))
    reg._load_user_adapters()
    names = set(list_steps())
    assert "wb3_user_demo_parser" in names, (
        f"user-adapter scaffold not discovered (saw {len(names)} steps)"
    )


def test_second_scaffold_in_same_tree_succeeds(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    parser = cli_new_adapter.build_parser(argparse_Sub())
    for name in ("wb3_two_a", "wb3_two_b"):
        args = parser.parse_args([name, "--kind", "parser", "--dir", str(tmp_path)])
        assert cli_new_adapter.handler(args) == 0, name
    docs = tmp_path / "docs" / "extending_delfin.md"
    text = docs.read_text(encoding="utf-8")
    assert "wb3_two_a_parser" in text and "wb3_two_b_parser" in text, (
        "docs stub must gain a section per adapter, not refuse adapter #2"
    )


def test_scaffold_refuses_to_replace_an_existing_adapter(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    parser = cli_new_adapter.build_parser(argparse_Sub())
    args = parser.parse_args(["wb3_once", "--kind", "parser", "--dir", str(tmp_path)])
    assert cli_new_adapter.handler(args) == 0
    rc = cli_new_adapter.handler(args)
    assert rc == 2, "same adapter twice must be refused (exit 2)"
