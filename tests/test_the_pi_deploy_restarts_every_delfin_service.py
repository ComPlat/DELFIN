"""A deploy must leave every running DELFIN service on the new code.

Before: the agent user's sudo rule named only delfin-dashboard and
pi_install.sh restarted only it, so the agent service and the messenger
ran the previous code after every deploy (seen 2026-11-03: started the day
before, three deploys later). The rule now names exactly the four DELFIN
units, each for restart and is-active only.
"""
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
UNITS = ("delfin-dashboard", "delfin-agent-service",
         "delfin-messenger-bridge", "delfin-messenger-gateway")


def _service_rule() -> str:
    text = (ROOT / "scripts/pi/pi_bootstrap.sh").read_text()
    block = text[text.index("    service)"):text.index('> "$SUDOERS.tmp" ;;', text.index("    service)"))]
    script = (block.split("\n", 1)[1].replace('systemctl_bin="$(command -v systemctl)"',
                                              'systemctl_bin=/usr/bin/systemctl')
              + '> /dev/stdout\n')
    out = subprocess.run(["bash", "-c", 'AGENT_USER=delfin\n' + script],
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()


def test_the_rule_names_exactly_the_four_units_for_restart_and_is_active():
    rule = _service_rule()
    assert rule.startswith("delfin ALL=(root) NOPASSWD: ")
    commands = [c.strip() for c in rule.split("NOPASSWD:", 1)[1].split(",")]
    expected = sorted(f"/usr/bin/systemctl {verb} {unit}"
                      for unit in UNITS for verb in ("restart", "is-active"))
    assert sorted(commands) == expected


def test_the_deploy_restarts_each_running_unit_and_names_a_stale_one():
    text = (ROOT / "scripts/pi/pi_install.sh").read_text()
    body = text[text.index("restart_dashboard() {"):]
    body = body[:body.index("\n}\n") + 3]
    for unit in UNITS:
        assert unit in body
    assert re.search(r"is-active --quiet \"\$unit\"", body)
    assert 'sudo -n systemctl restart "$unit"' in body
    assert "still on the previous code until an admin runs" in body
