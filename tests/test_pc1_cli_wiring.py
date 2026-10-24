"""The `delfin prcheck` subcommand wiring (package PC1, phase 3b).

Red control: seen RED through the gate on b391d899 (no prcheck
subcommand; argparse fell through to chat's help) before the operator
built .gate/pc1_cli_prcheck_wiring.patch as c31874c2. Green since.
"""

import subprocess
import sys


def test_delfin_cli_has_a_prcheck_subcommand():
    probe = subprocess.run(
        [sys.executable, "-m", "delfin.agent.cli", "prcheck", "--help"],
        capture_output=True, text=True, timeout=120)
    assert probe.returncode == 0, probe.stderr
    # The subcommand must resolve on its own (the red state fell through
    # to chat's help with "usage: delfin-agent chat ...").
    assert "usage: delfin-agent prcheck" in probe.stdout


def test_prcheck_accepts_the_documented_flags():
    probe = subprocess.run(
        [sys.executable, "-m", "delfin.agent.cli", "prcheck", "--help"],
        capture_output=True, text=True, timeout=120)
    assert "--base" in probe.stdout
    assert "--partner" in probe.stdout
