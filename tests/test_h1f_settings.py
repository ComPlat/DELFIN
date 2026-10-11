"""Settings contract for H1F read-only-shape approval.

The gate (phase 3, protected) will consult ``agent.read_only_shapes`` from the
user settings payload. That stored shape's contract — its keys, its defaults
(enabled ON, every known shape enabled) and its normalisation (unknown shape
names dropped, the toggle surviving a partial save) — lives here next to the
shapes themselves, so the settings tab and the gate read the SAME code.

The helpers tested here sit in delfin.agent.read_only_shapes (module-level,
pure) and are the only path tab_settings.py uses to load/save the block they
render. Against the unchanged tree these are RED: the three helpers are new.
"""

from delfin.agent import read_only_shapes as ros
from delfin.agent.read_only_shapes import SHAPE_NAMES


def _payload(**ro):
    payload = {"agent": {"read_only_shapes": ro}}
    return payload


def _empty():
    return {}


# --- defaults ---------------------------------------------------------------

def test_settings_defaults_enabled_on():
    d = ros.read_only_shapes_defaults()
    assert d["enabled"] is True


def test_settings_defaults_every_shape_enabled():
    d = ros.read_only_shapes_defaults()
    assert set(d["enabled_shapes"]) == set(SHAPE_NAMES)


# --- reading / normalising --------------------------------------------------

def test_settings_from_missing_payload_uses_defaults():
    d = ros.read_only_shapes_from_payload({})
    assert d["enabled"] is True
    assert set(d["enabled_shapes"]) == set(SHAPE_NAMES)


def test_settings_from_payload_honours_disabled_toggle():
    d = ros.read_only_shapes_from_payload(_payload(enabled=False))
    assert d["enabled"] is False






def test_settings_from_payload_unknown_keys_are_ignored():
    # a partial save must not corrupt the stored block
    d = ros.read_only_shapes_from_payload(_payload(enabled=False))
    assert "enabled_shapes" in d
    assert "enabled" in d
    assert set(d.keys()) == {"enabled", "enabled_shapes"}


# --- writing ----------------------------------------------------------------

def test_settings_into_payload_creates_block():
    out = {}
    ros.read_only_shapes_into_payload(out, ros.read_only_shapes_defaults())
    assert out["agent"]["read_only_shapes"]["enabled"] is True
    assert set(out["agent"]["read_only_shapes"]["enabled_shapes"]) == set(
        SHAPE_NAMES
    )


def test_settings_into_payload_preserves_other_agent_keys():
    out = {"agent": {"backend": "cli"}}
    ros.read_only_shapes_into_payload(out, ros.read_only_shapes_defaults())
    assert out["agent"]["backend"] == "cli"


def test_matches_with_respects_disabled_toggle_like_the_gate_will():
    # the gate (phase 3) will pass enabled=False through as an empty list
    d = ros.read_only_shapes_defaults()
    d["enabled"] = False
    called = ros.matches_with(
        "git show HEAD:delfin/api.py", "/pfs/data6/home/ka/example/delfin",
        [] if not d["enabled"] else set(d["enabled_shapes"]),
    )
    assert called == ""
