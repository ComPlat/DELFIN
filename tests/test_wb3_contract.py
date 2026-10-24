"""WB3 phase 3: a scaffolded adapter passes the real contract tests.

End-to-end: scaffold into a throwaway tree, load the skeleton through the
REAL user-adapter loader, then run the same expectations
``tests/test_tool_contracts.py`` applies to every adapter: contract()
builds, category declared, required params enforced, declared params
accepted.  This is the target-1 acceptance: "a scaffolded adapter passes
the contract tests and is listed by the registry with no hand edit".
"""
from __future__ import annotations

import argparse

import pytest

import delfin.tools._registry as reg
from delfin.agent import cli_new_adapter
from delfin.tools._registry import list_steps


def argparse_Sub():
    return argparse.ArgumentParser().add_subparsers(dest="cmd", required=True)


@pytest.fixture()
def scaffolded(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(tmp_path / "user_adapters"))
    parser = cli_new_adapter.build_parser(argparse_Sub())
    repo = tmp_path / "repo"
    (repo / "delfin" / "tools" / "adapters").mkdir(parents=True)
    for kind in ("program", "parser", "analysis"):
        name = f"wb3_e2e_{kind}"
        args = parser.parse_args([name, "--kind", kind, "--dir", str(repo)])
        assert cli_new_adapter.handler(args) == 0, kind
        skeleton = repo / "delfin" / "tools" / "adapters" / f"{name}.py"
        reg._DISCOVERED = False
        # keep prior registrations visible so all three coexist
        reg.load_user_adapter_file(skeleton)
    yield
    reg._DISCOVERED = False
    reg._REGISTRY.clear()


def test_all_kinds_scaffold_and_register(scaffolded):
    names = set(list_steps())
    for kind in ("program", "parser", "analysis"):
        assert f"wb3_e2e_{kind}" in names, kind


@pytest.mark.parametrize("kind", ["program", "parser", "analysis"])
def test_scaffolded_adapter_passes_contract_expectations(scaffolded, kind):
    """The invariants tests/test_tool_contracts.py runs on every adapter."""
    adapter = list_steps()[f"wb3_e2e_{kind}"]
    contract = adapter.contract()
    assert contract.category == kind, "category must be declared"
    assert contract.params, "params must be declared"
    assert contract.name == f"wb3_e2e_{kind}"

    # required params are enforced by validate_params
    dummy = {"path": "x.xyz", "str": "x"}.get
    required = contract.required_params
    full = {p.name: dummy(p.type, "x") for p in contract.params}
    for req in required:
        kwargs = dict(full)
        kwargs.pop(req, None)
        with pytest.raises(ValueError):
            adapter.validate_params(**kwargs)
    # a fully-specified call passes validation
    adapter.validate_params(**full)
