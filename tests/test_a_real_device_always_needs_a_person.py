"""Every write to a device that is not the simulator needs a person's
single-use grant, whatever its level and whatever the permission mode;
the simulator keeps its own levels (decided with the maintainer,
2026-11-26: open slowly)."""
from __future__ import annotations

import pytest

from _lab_support import call, lab_fixture, sim  # noqa: F401  (fixture)


def _make_real(device: str) -> None:
    from delfin.lab import manifest, store
    from delfin.lab.safety import pin
    path = store.manifest_path()
    doc = manifest.load(path)
    text = path.read_text(encoding="utf-8")
    lines, out, inside = text.splitlines(), [], False
    for line in lines:
        if line.strip().startswith("- id:"):
            inside = line.split(":", 1)[1].strip() == device
        if inside and line.strip() == "driver: simulator":
            line = line.replace("driver: simulator", "driver: tasmota")
            inside = False
        out.append(line)
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    pin(manifest.load(path), by="alice")
    assert manifest.load(path).devices[device].driver == "tasmota", doc


@pytest.mark.parametrize("mode", ["bypassPermissions", "default"])
def test_a_harmless_step_on_real_hardware_waits_for_a_person(lab, mode):
    _make_real("stirrer")
    out = call(lab, "lab_set", {"device": "stirrer", "channel": "power", "value": "on"}, mode=mode)
    assert out.get("code") == "approval_required", out
    assert "real hardware" in out.get("reason", ""), out


def test_the_simulator_keeps_its_own_levels(lab):
    out = call(lab, "lab_set", {"device": "stirrer", "channel": "power", "value": "on"})
    assert out.get("ok"), out
    assert sim("stirrer")["power"] == "on"
