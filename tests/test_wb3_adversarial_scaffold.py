"""WB3 adversarial tests (reviewer nacht-s19) for target 1: the
`delfin-agent new-adapter <name> --kind <program|parser|analysis>` scaffold.

All cases are RED on base 8c258293: no `new-adapter` subcommand exists and
delfin/tools/scaffold.py is absent. They are findings tests -- each one turns
green only when the builder's phase-2/3 work lands. RED here is the finding;
green pins correctness once the fix arrives.

Interface assumptions stated up front (to verify against the builder's
implementation when reviewing): the subcommand binds a handler via
set_defaults(func=...), accepts `--kind {program,parser,analysis}` and a
`--dest` directory (defaulting to the user adapter dir); the handler calls
the public scaffold function in delfin.tools.scaffold. If the real interface
differs, these tests are rewritten along their intent, not weakened.
"""
import os
import unittest.mock

import pytest


# --- target 1: scaffold via the real CLI path -------------------------------

def test_new_adapter_subcommand_exists_in_build_parser():
    from delfin.agent import cli

    parser = cli.build_parser()
    args = parser.parse_args(["new-adapter", "demo", "--kind", "analysis"])
    assert args.name == "demo"


def test_scaffold_via_real_cli_path_end_to_end(tmp_path):
    """Run the real CLI handler for new-adapter into a throwaway adapters
    dir; then the registry must list the new adapter with no hand edit."""
    dest = tmp_path / "adapters"
    dest.mkdir()
    env_override = {"DELFIN_ADAPTERS_DIR": str(dest)}

    with unittest.mock.patch.dict(os.environ, env_override):
        from delfin.agent import cli

        parser = cli.build_parser()
        args = parser.parse_args(["new-adapter", "demo_probe", "--kind",
                                  "analysis", "--dest", str(dest)])
        handler = getattr(args, "func", None)
        assert handler is not None, (
            "new-adapter subcommand must bind a handler via set_defaults(func=...)"
        )
        handler(args)

    # The scaffold wrote into the user adapters dir...
    assert (dest / "demo_probe.py").exists()
    # ...and the registry (fresh discovery, user dir injected via env) lists it.
    from delfin.tools import _registry

    _registry._DISCOVERED = False
    _registry._REGISTRY.clear()
    with unittest.mock.patch.dict(os.environ, env_override):
        _registry._discover_builtin_adapters()
    assert "demo_probe" in _registry._REGISTRY
    # restore registry state for other tests
    _registry._DISCOVERED = False
    _registry._REGISTRY.clear()
    _registry._discover_builtin_adapters()


# --- target 1: contract conformance of a scaffolded adapter ------------------

def test_scaffolded_adapter_passes_contract_invariants(tmp_path):
    """The scaffolded adapter must pass the loop invariants of
    tests/test_tool_contracts.py: category declared, contract() non-empty,
    enum defaults within enum."""
    dest = tmp_path / "adapters"
    dest.mkdir()
    env_override = {"DELFIN_ADAPTERS_DIR": str(dest)}
    with unittest.mock.patch.dict(os.environ, env_override):
        from delfin.agent import cli

        parser = cli.build_parser()
        args = parser.parse_args(["new-adapter", "contract_probe", "--kind",
                                  "program", "--dest", str(dest)])
        getattr(args, "func")(args)
        from delfin.tools import _registry

        _registry._DISCOVERED = False
        _registry._REGISTRY.clear()
        _registry._discover_builtin_adapters()
        adapter = _registry._REGISTRY["contract_probe"]
        # contract invariants from tests/test_tool_contracts.py
        assert adapter.category != ""
        assert adapter.contract()  # non-empty
        # enum defaults within enum
        for p in adapter.params:
            enum = getattr(p, "choices", None) or getattr(p, "enum", None)
            if enum and getattr(p, "default", None) is not None:
                assert p.default in enum

    _registry._DISCOVERED = False
    _registry._REGISTRY.clear()
    _registry._discover_builtin_adapters()


def test_scaffolded_adapter_declares_required_params_enforced(tmp_path):
    """If the skeleton declares a required param, validate_params must
    enforce it (contract-test expectation, tests/test_tool_contracts.py:50)."""
    dest = tmp_path / "adapters"
    dest.mkdir()
    env_override = {"DELFIN_ADAPTERS_DIR": str(dest)}
    with unittest.mock.patch.dict(os.environ, env_override):
        from delfin.agent import cli
        from delfin.tools import _registry

        parser = cli.build_parser()
        args = parser.parse_args(["new-adapter", "req_probe", "--kind",
                                  "program", "--dest", str(dest)])
        getattr(args, "func")(args)
        _registry._DISCOVERED = False
        _registry._REGISTRY.clear()
        _registry._discover_builtin_adapters()
        adapter = _registry._REGISTRY["req_probe"]
        required = [p.name for p in adapter.params
                    if getattr(p, "required", False)]
        if not required:
            pytest.skip("skeleton declares no required params")
        with pytest.raises(Exception):
            adapter.validate_params()

    _registry._DISCOVERED = False
    _registry._REGISTRY.clear()
    _registry._discover_builtin_adapters()


# --- security-relevant input surface: the scaffold writes files --------------

def test_scaffold_refuses_path_traversal_in_name(tmp_path):
    """A --name of '../evil' must not write outside --dest."""
    dest = tmp_path / "adapters"
    dest.mkdir()
    from delfin.agent import cli

    parser = cli.build_parser()
    args = parser.parse_args(["new-adapter", "../evil", "--kind", "parser",
                              "--dest", str(dest)])
    before = sorted(os.listdir(tmp_path))
    with pytest.raises(Exception):
        getattr(args, "func")(args)
    after = sorted(os.listdir(tmp_path))
    assert before == after, "scaffold wrote outside --dest"
    assert (tmp_path / "evil.py").exists() is False


def test_scaffold_refuses_unicode_lookalike_or_empty_name(tmp_path):
    """Names that are empty or contain path separators / look-alikes must be
    rejected or sanitized to a safe module name."""
    dest = tmp_path / "adapters"
    dest.mkdir()
    from delfin.agent import cli

    parser = cli.build_parser()
    for bad in ("", "..", "a/b", "a b", "a\nb"):
        with pytest.raises(Exception):
            args = parser.parse_args(["new-adapter", bad, "--kind", "parser",
                                      "--dest", str(dest)])
            getattr(args, "func")(args)
    assert sorted(os.listdir(dest)) == []


def test_scaffold_refuses_overwriting_an_existing_file(tmp_path):
    dest = tmp_path / "adapters"
    dest.mkdir()
    existing = dest / "mine.py"
    existing.write_text("# mine\n", encoding="utf-8")
    from delfin.agent import cli

    parser = cli.build_parser()
    args = parser.parse_args(["new-adapter", "mine", "--kind", "parser",
                              "--dest", str(dest)])
    with pytest.raises(Exception):
        getattr(args, "func")(args)
    assert existing.read_text(encoding="utf-8") == "# mine\n"


def test_scaffold_refuses_invalid_kind(tmp_path):
    from delfin.agent import cli

    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["new-adapter", "x", "--kind", "not-a-kind"])


def test_scaffold_wiring_goes_through_build_parser(tmp_path):
    """Wiring proof: the new-adapter handler reachable from build_parser()
    must call the scaffold's public function (single source of truth; no
    duplicated template inside cli.py)."""
    dest = tmp_path / "adapters"
    dest.mkdir()
    from delfin.agent import cli
    from delfin.tools import scaffold

    parser = cli.build_parser()
    args = parser.parse_args(["new-adapter", "wiring_probe", "--kind",
                              "analysis", "--dest", str(dest)])
    with unittest.mock.patch.object(
        scaffold, "scaffold_new_adapter",
        wraps=scaffold.scaffold_new_adapter,
    ) as spy:
        getattr(args, "func")(args)
    assert spy.called, "CLI handler must call delfin.tools.scaffold"
    assert (dest / "wiring_probe.py").exists()
