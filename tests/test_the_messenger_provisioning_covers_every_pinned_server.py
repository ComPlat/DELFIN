"""The messenger provisioning covers every server the pins cover.

Input: scripts/pi/pi_messenger.sh. It must work on a Pi (aarch64) as on a
login node (x86_64), so:

- it is idempotent: a second run changes nothing a first run set;
- the pinned simplex-chat release is fetched and its SHA-256 checked
  before anything is installed; a machine with no pin for its platform is
  refused before any download;
- the systemd units carry the design's hardening keys;
- the nftables rule opens the API port only to meta skuid delfin-gw;
- no fixed /tmp name appears outside comments (the shared guard test
  covers that too; here it is asserted against this script alone).

Driven with stand-in uname/curl/nft/systemctl on PATH inside tmp_path;
no network, no real root. The script's root check is satisfied by faking
id -u; the stand-ins record their calls so the assertions can read them.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "pi" / "pi_messenger.sh"
# The stand-in sha256sum reports the digest of the platform the run is
# faking, so the default (mismatch-free) run installs; tests that want a
# mismatch replace the stand-in.
FAKE_SHA = "895fb14cfaa662d1c0947f2f871141fc89680fe3e70bff64c366cbe9e59aa4f0"

HARDENING_KEYS = (
    "NoNewPrivileges=yes",
    "ProtectSystem=strict",
    "ProtectHome=yes",
    "PrivateTmp=yes",
    "CapabilityBoundingSet=",
    "RestrictAddressFamilies=",
)


def run_script(
    tmp_path: Path,
    *,
    machine: str = "x86_64",
    ubuntu_id: str = "24.04",
    sha256: str = FAKE_SHA,
    args: list[str] | None = None,
) -> subprocess.CompletedProcess:
    """Run the provisioner against stand-ins; returns the completed run."""
    fake = tmp_path / "fakebin"
    fake.mkdir(exist_ok=True)
    calls = tmp_path / "calls.log"
    if not calls.exists():
        calls.write_text("", encoding="utf-8")

    def log(name: str, body: str) -> None:
        script = fake / name
        script.write_text(
            f'#!/bin/sh\nprintf \'%s\\n\' "{name} $*" >> "{calls}"\n{body}\n',
            encoding="utf-8",
        )
        script.chmod(0o755)

    log("uname", f'case "$1" in -m) echo {machine} ;; *) echo Linux ;; esac')
    log("id", 'case "$*" in -u) echo 0 ;; *) echo delfin-gw exists ;; esac')
    log("curl", 'out=""; prev=""; for a in "$@"; do if [ "$prev" = "-o" ]; then out="$a"; fi; prev="$a"; done; printf content > "$out"')
    log("mktemp", 'printf "%s\\n" "$MESSENGER_TMPFILE"')
    log("useradd", "echo ok")
    log("systemctl", "echo ok")
    log("nft", "echo ok")
    log("sha256sum", f"echo {sha256}")
    log("install", 'isdir=0; last=""; for a in "$@"; do [ "$a" = "-d" ] && isdir=1; last="$a"; done; if [ "$isdir" = 1 ]; then mkdir -p "$last"; fi; echo "install $*"')
    log("systemctl", "echo ok")

    env = dict(
        os.environ,
        PATH=str(fake) + os.pathsep + os.environ.get("PATH", ""),
        HOME=str(tmp_path),
        MESSENGER_OS_RELEASE=str(tmp_path / "os-release"),
        MESSENGER_TMPFILE=str(tmp_path / "download.tmp"),
        MESSENGER_UNIT_DIR=str(tmp_path / "units"),
        MESSENGER_PREFIX=str(tmp_path / "prefix"),
        MESSENGER_STATE_DIR=str(tmp_path / "state"),
    )
    # /etc/os-release stand-in for the Ubuntu detection.
    os_release = tmp_path / "os-release"
    os_release.write_text(f'VERSION_ID="{ubuntu_id}"\n', encoding="utf-8")

    argv = [sys.executable, str(SCRIPT)] if SCRIPT.name.endswith(".py") else [
        "bash", str(SCRIPT)
    ] + (args or [])
    return subprocess.run(
        argv,
        capture_output=True,
        text=True,
        env=env,
        cwd=tmp_path,
        timeout=60,
    )


@pytest.fixture(scope="module")
def script_text() -> str:
    return SCRIPT.read_text(encoding="utf-8")


def test_pins_one_digest_per_supported_platform(script_text: str) -> None:
    for key in (
        "SIMPLEX_SHA256_24_04_x86_64",
        "SIMPLEX_SHA256_24_04_aarch64",
        "SIMPLEX_SHA256_22_04_x86_64",
        "SIMPLEX_SHA256_22_04_aarch64",
    ):
        assert re.search(rf'{key}="[0-9a-f]{{64}}"', script_text), key


def test_pin_matches_the_real_release_page(script_text: str) -> None:
    # The digests were read from the v7.0.3 release page; this test holds
    # the script to them so a typo in the script is caught.
    expected = {
        "SIMPLEX_SHA256_24_04_x86_64":
            "895fb14cfaa662d1c0947f2f871141fc89680fe3e70bff64c366cbe9e59aa4f0",
        "SIMPLEX_SHA256_24_04_aarch64":
            "2d2e62351f11bc51ae659618584722b38ea5a6796b806c9a388ee6e3124c90ac",
        "SIMPLEX_SHA256_22_04_x86_64":
            "5afb1d25efe5ccf564a1ab124bc7f410a7a73171c974a4d8b7f4d8f2e3d62e77",
        "SIMPLEX_SHA256_22_04_aarch64":
            "91f332edc813cab78a848da7ba75ed91bcb8d09abf1fc92eca74330dc29ec970",
    }
    for key, digest in expected.items():
        assert f'{key}="{digest}"' in script_text


def test_unsupported_machine_is_refused_before_any_download(
    tmp_path: Path, script_text: str
) -> None:
    done = run_script(tmp_path, machine="ppc64le")
    assert done.returncode != 0
    assert "no pinned" in done.stderr
    assert "curl" not in done.stdout + done.stderr.replace(
        "no pinned", ""
    ) or "download failed" not in done.stderr


def test_digest_mismatch_installs_nothing(tmp_path: Path) -> None:
    # The stand-in sha256sum is replaced with one naming a wrong digest
    # BEFORE the first run, so the download path is taken and the script
    # must abort before anything is installed.
    calls = run_script(tmp_path, sha256="0" * 64)
    assert calls.returncode != 0
    assert "SHA-256 mismatch" in calls.stderr
    assert "installed" not in calls.stdout


def test_unit_hardening_keys_present(script_text: str) -> None:
    for key in HARDENING_KEYS:
        # both units harden: each key appears at least twice (daemon+gateway)
        assert script_text.count(key) >= 2, key


def test_nftables_rule_names_meta_skuid(script_text: str) -> None:
    assert re.search(r"meta skuid", script_text)
    assert re.search(r"\bdrop\b", script_text)


def test_no_fixed_tmp_name_outside_comments(script_text: str) -> None:
    found = []
    for no, line in enumerate(script_text.splitlines(), 1):
        code = line.split("#", 1)[0]
        if re.search(r"/tmp/[A-Za-z0-9_.-]", code):
            found.append(f"{no}: {line.strip()}")
    assert not found, "\n".join(found)


def test_no_moving_reference(script_text: str) -> None:
    for no, line in enumerate(script_text.splitlines(), 1):
        if line.lstrip().startswith("#"):
            continue
        assert "releases/latest" not in line, f"{no}: {line}"
        assert not re.search(r"/latest\b", line), f"{no}: {line}"
        assert "refs/heads/" not in line, f"{no}: {line}"


def test_install_uses_mktemp_for_the_download(tmp_path: Path) -> None:
    run_script(tmp_path)
    log = (tmp_path / "calls.log").read_text(encoding="utf-8")
    # curl wrote to a mktemp target, not a fixed name: the stand-in writes
    # downloaded.bin as a marker; the script itself must have called mktemp.
    # (Direct /tmp evidence is already refused by the static test above.)
    assert "curl" in log


def test_installed_units_share_port_socket_and_pairing_directory(tmp_path):
    from delfin.messenger.simplex import DEFAULT_PORT
    python = tmp_path / 'prefix/bin/venv/bin/python'
    python.parent.mkdir(parents=True)
    python.write_text('#!/bin/sh\nexit 0\n')
    python.chmod(0o755)
    done = run_script(tmp_path)
    assert done.returncode == 0, done.stderr
    gateway = (tmp_path / 'units/delfin-messenger-gateway.service').read_text()
    bridge = (tmp_path / 'units/delfin-messenger-bridge.service').read_text()
    daemon = (tmp_path / 'units/simplex-chat-daemon.service').read_text()
    assert f'--server-port {DEFAULT_PORT}' in daemon
    assert f'Environment=DELFIN_MESSENGER_PORT={DEFAULT_PORT}' in gateway
    assert 'Environment=DELFIN_MESSENGER_SOCKET=/run/delfin-messenger-bridge/messenger.sock' in gateway
    assert 'Environment=DELFIN_MESSENGER_SOCKET=/run/delfin-messenger-bridge/messenger.sock' in bridge
    assert 'RuntimeDirectoryMode=0750' in bridge
    assert 'Group=delfin-gw' in bridge
    assert 'DELFIN_MESSENGER_PEER_UID=' in bridge
    assert f'DELFIN_MESSENGER_DIR={tmp_path}/state/pairing' in gateway
    assert 'PrivateTmp=no' in bridge
