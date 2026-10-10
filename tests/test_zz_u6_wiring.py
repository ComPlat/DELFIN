"""U6 phase-3 RED TEST (moved to .gate/ after the red control is confirmed).

The refusal message must carry the gate-legal suggestion for the deny shapes
that dominate the wave-13 log (.gate/refusals_w13.txt). This is RED on the
unpatched gate: _run_permission_gate refuses `rm -rf build` with
"command rejected by deny-pattern '...': refusing to run." and names no
route. After the phase-3 wiring (api_client.py:15432 appends
refusal_hints.suggest(command, denied)) the message must name the allowed
form (`mv build .gate/`).
"""
import pytest

from delfin.agent import api_client as A


@pytest.fixture
def perms(tmp_path):
    return A.KitToolPermissions(workspace=tmp_path)


@pytest.mark.parametrize("refused,sought", [
    # rm -rf hits the deny-pattern branch; the suggestion names mv to scratch.
    ("rm -rf build", "mv build .gate/"),
    ("rm -f -r dist", "mv dist .gate/"),
])
def test_deny_refusal_carries_the_suggestion(perms, refused, sought):
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    err = ex._run_permission_gate("bash", {"command": refused}, perms)
    assert err is not None, f"expected a refusal for {refused!r}"
    assert "deny-pattern" in err, err
    assert sought in err, (
        f"refusal for {refused!r} does not name the allowed form:\n{err}"
    )
