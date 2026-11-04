"""Phase C: external_cli backend type in the registry.

The orca-mlips ProgExt plugins are NOT ASE-calculator backends: they run
inside ORCA ExtOpt via a ProgExt executable.  These tests pin the
contract: registry answers availability fail-closed, agents can see the
restriction, and the ASE factory refuses to route them.
"""
import sys
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.mlp_tools import registry as reg


def test_orca_mlips_backends_are_registered_and_typed():
    for name in ("orca_mlips_aimnet2", "orca_mlips_uma",
                 "orca_mlips_orb", "orca_mlips_mace"):
        d = reg.describe_backend(name)
        assert d["backend_type"] == "external_cli"
        assert d["executable"].startswith("orca-mlips-")


def test_unavailable_executable_is_fail_closed_with_hint():
    d = reg.describe_backend("orca_mlips_aimnet2")
    assert d["available"] is False
    assert "not found on PATH" in d["reason"]
    assert "pip install orca-mlips" in d["install_hint"]


def test_executable_on_path_makes_aimnet2_available_but_not_uma():
    with mock.patch("shutil.which", return_value="/opt/x/orca-mlips-aimnet2"):
        avail, reason = reg.check_availability("orca_mlips_aimnet2")
        assert avail is True and reason == ""
    # UMA carries a weights-licence note -> restricted even when found
    with mock.patch("shutil.which", return_value="/opt/x/orca-mlips-uma"):
        avail, reason = reg.check_availability("orca_mlips_uma")
        assert avail is True and "licence-restricted" in reason


def test_factory_refuses_external_cli_with_guidance():
    import pytest
    with pytest.raises(NotImplementedError) as ei:
        reg.create_calculator("orca_mlips_aimnet2")
    msg = str(ei.value)
    assert "external_cli" in msg and "ProgExt" in msg and "ExtOpt" in msg


def test_legacy_backends_keep_their_type():
    d = reg.describe_backend("aimnet2")
    assert d["backend_type"] == "ml_potential"
    assert d["executable"] is None
