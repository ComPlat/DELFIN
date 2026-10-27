"""The Pi scripts run no file they put under a fixed name in /tmp.

Input: scripts/pi/*.sh. Two of them wrote a script to a fixed path in /tmp
and ran it afterwards: deploy_pi.sh --bootstrap uploaded
/tmp/pi_bootstrap.sh and ran it with sudo; pi_bootstrap.sh downloaded the
CLI installer to /tmp/claude-install.sh and ran it as the agent user, who
has passwordless sudo. Any other local user can create that name first
(or swap the file between the write and the run) and have their code run
as root. A name from mktemp is created 0600 by its owner, and in the
sticky /tmp nobody else can replace it.
"""
from __future__ import annotations

import pathlib
import re

_PI = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "pi"


def test_no_pi_script_names_a_fixed_file_in_tmp():
    found = []
    for f in sorted(_PI.glob("*.sh")):
        for no, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            code = line.split("#", 1)[0]
            if re.search(r"/tmp/[A-Za-z0-9_.-]", code):
                found.append(f"{f.name}:{no}: {line.strip()}")
    assert not found, "\n".join(found)
