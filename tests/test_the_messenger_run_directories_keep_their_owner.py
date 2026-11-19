"""The run directories the messenger's units share keep the owner they are made with.

Field defect (Pi, fixed on main in #584 for the bridge's socket directory):
systemd chowns a RuntimeDirectory= to the unit's User=/Group= when it sets
up ExecStart, after an ExecStartPre=+ chgrp, so the gateway lost its group
on the bridge's socket directory on every start. This branch added a second
one, /run/delfin-person-start (the bridge leaves person-start requests, the
root starter reads them), made the same way.

Output: neither directory is a RuntimeDirectory=; the bridge unit creates
both root-side (``install -d``) with the owner and mode they need -- the
socket directory 2750 service user:delfin-gw, the request directory 0700
service user only -- and removes both when it stops.
"""
from __future__ import annotations

from test_the_messenger_provisioning_covers_every_pinned_server import run_script


def _bridge(tmp_path):
    venv = tmp_path / "home-delfin" / "delfin/.venv/bin/python"
    venv.parent.mkdir(parents=True)
    venv.write_text("#!/bin/sh\nexit 0\n")
    venv.chmod(0o755)
    done = run_script(tmp_path)
    assert done.returncode == 0, done.stderr
    unit = (tmp_path / "units" / "delfin-messenger-bridge.service").read_text()
    return [line.strip() for line in unit.splitlines()]


def test_both_run_directories_are_made_root_side_with_their_owner(tmp_path):
    lines = _bridge(tmp_path)
    assert not any(line.startswith("RuntimeDirectory") for line in lines)
    assert ("ExecStartPre=+/usr/bin/install -d -m 2750 -o delfin -g delfin-gw "
            "/run/delfin-messenger-bridge") in lines
    assert ("ExecStartPre=+/usr/bin/install -d -m 0700 -o delfin -g delfin "
            "/run/delfin-person-start") in lines
    for path in ("/run/delfin-messenger-bridge", "/run/delfin-person-start"):
        assert any(line.startswith("ExecStopPost=+") and path in line for line in lines), path
    assert "Environment=DELFIN_MESSENGER_PERSON_STARTS=/run/delfin-person-start" in lines
