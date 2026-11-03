"""Licence surface of the agent tool + the refuse-non-commercial option.

Task om1 item 3: mlp_list_backends shows every backend's code and weight
licence (with URL + verified flag), and a filter refuses backends with
non-commercial weight licences when requested.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from delfin.mlp_tools import registry  # noqa: E402


def test_describe_backend_includes_licences():
    d = registry.describe_backend("mace_omol")
    assert d["code_licence"]["name"] == "MIT"
    assert d["weight_licence"]["non_commercial"] is True
    assert d["weight_licence"]["url"].startswith("http")
    uma = registry.describe_backend("fairchem_uma")
    assert uma["weight_licence"]["gated"] is True


def test_mlp_list_backends_tool_shows_licences():
    """The agent-facing tool (delfin.api.mlp_list_backends) carries them."""
    from delfin.api import mlp_list_backends

    out = mlp_list_backends()
    by_name = {b["name"]: b for b in out["backends"]}
    assert by_name["mace_omol"]["weight_licence"]["non_commercial"] is True
    assert by_name["orb"]["weight_licence"]["name"] == "Apache-2.0"


def test_refuse_noncommercial_filter():
    keep, refused = registry.filter_backends_by_licence(
        ["mace_omol", "orb", "fairchem_uma"],
        refuse_noncommercial=True,
    )
    assert "mace_omol" in refused
    assert "mace_omol" not in keep
    assert "orb" in keep          # Apache-2.0 weights
    assert "fairchem_uma" in keep  # gated but not non-commercial


def test_refuse_filter_noop_when_off():
    keep, refused = registry.filter_backends_by_licence(
        ["mace_omol", "orb"], refuse_noncommercial=False
    )
    assert refused == []
    assert set(keep) == {"mace_omol", "orb"}


def test_gate_applies_to_registry_factory():
    """create_calculator honours refuse_noncommercial_weights."""
    import pytest

    try:
        registry.create_calculator(
            "mace_omol", refuse_noncommercial_weights=True
        )
    except ValueError as exc:
        assert "non-commercial" in str(exc).lower()
    else:
        import pytest as _p

        _p.fail("refuse_noncommercial_weights did not refuse mace_omol")
