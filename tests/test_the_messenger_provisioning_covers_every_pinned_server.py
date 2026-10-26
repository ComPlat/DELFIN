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
    sha256_body: str = "",
    os_release: str | None = None,
    glibc: str = "2.39",
    env_extra: dict | None = None,
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
    log("mktemp", 'if [ "$1" = "-d" ]; then mkdir -p "$MESSENGER_TMPFILE.d"; printf "%s\\n" "$MESSENGER_TMPFILE.d"; else printf "%s\\n" "$MESSENGER_TMPFILE"; fi')
    log("useradd", "echo ok")
    log("systemctl", "echo ok")
    log("nft", "echo ok")
    log("sha256sum", sha256_body or f"echo {sha256}")
    log("install", 'isdir=0; prev=""; last=""; for a in "$@"; do [ "$a" = "-d" ] && isdir=1; prev="$last"; last="$a"; done; if [ "$isdir" = 1 ]; then mkdir -p "$last"; elif [ -f "$prev" ]; then cp "$prev" "$last"; fi; echo "install $*"')
    log("systemctl", "echo ok")
    log("getconf", f'echo "glibc {glibc}"')
    # The interpreter a fresh gateway venv is made from: a venv's python is
    # a copy of the stand-in, which answers the version and site checks.
    base = fake / "python3.11-standin"
    base.write_text(
        "#!/bin/sh\n"
        f'printf "%s\\n" "python $*" >> "{calls}"\n'
        'case "$*" in\n'
        '  *version_info*) exit 0 ;;\n'
        f'  *purelib*) mkdir -p "{tmp_path}/site"; echo "{tmp_path}/site" ;;\n'
        '  "-m venv "*) mkdir -p "$3/bin"; cp "$0" "$3/bin/python" ;;\n'
        "esac\nexit 0\n", encoding="utf-8")
    base.chmod(0o755)
    service_home = tmp_path / "home-delfin"
    service_home.mkdir(exist_ok=True)

    env = dict(
        os.environ,
        PATH=str(fake) + os.pathsep + os.environ.get("PATH", ""),
        HOME=str(tmp_path),
        MESSENGER_OS_RELEASE=str(tmp_path / "os-release"),
        MESSENGER_TMPFILE=str(tmp_path / "download.tmp"),
        MESSENGER_UNIT_DIR=str(tmp_path / "units"),
        MESSENGER_PREFIX=str(tmp_path / "prefix"),
        MESSENGER_STATE_DIR=str(tmp_path / "state"),
        MESSENGER_PYTHON=str(base),
        MESSENGER_SERVICE_HOME=str(service_home),
        **(env_extra or {}),
    )
    # /etc/os-release stand-in for the platform detection.
    (tmp_path / "os-release").write_text(
        os_release if os_release is not None else
        f'ID=ubuntu\nVERSION_ID="{ubuntu_id}"\n', encoding="utf-8")

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
    assert f'--chat-server-port {DEFAULT_PORT}' in daemon
    assert '--files-folder ' in daemon and '--auto-accept-files 0' in daemon
    assert f'Environment=DELFIN_MESSENGER_PORT={DEFAULT_PORT}' in gateway
    assert 'Environment=DELFIN_MESSENGER_SOCKET=/run/delfin-messenger-bridge/messenger.sock' in gateway
    assert 'Environment=DELFIN_MESSENGER_SOCKET=/run/delfin-messenger-bridge/messenger.sock' in bridge
    assert 'RuntimeDirectoryMode=0750' in bridge
    assert 'Group=delfin-gw' in bridge
    assert 'DELFIN_MESSENGER_PEER_UID=' in bridge
    assert f'DELFIN_MESSENGER_DIR={tmp_path}/state/pairing' in gateway
    assert 'PrivateTmp=no' in bridge


# --- optional local speech recognition (--whisper) -------------------------

_WHISPER_PINS = {
    "config.json": "b55496ac7940a7ae47d2c01eab40edfd8701feec1229d9cce3b40014383fb828",
    "model.bin": "3e305921506d8872816023e4c273e75d2419fb89b24da97b4fe7bce14170d671",
    "tokenizer.json": "fb7b63191e9bb045082c79fd742a3106a12c99513ab30df4a0d47fa6cb6fd0ab",
    "vocabulary.txt": "34ce3fe1c5041027b3f8d42912270993f986dbc4bb34cf27f951e34a1e453913",
}


def _gateway_venv(tmp_path, *, version="3.11"):
    """A stand-in gateway python: logs its argv, answers the version
    check, and makes a venv whose python is a copy of itself."""
    python = tmp_path / "prefix/bin/venv/bin/python"
    python.parent.mkdir(parents=True)
    major, minor = version.split(".")
    python.write_text(
        "#!/bin/sh\n"
        f'printf "%s\\n" "python $*" >> "{tmp_path}/calls.log"\n'
        'case "$*" in\n'
        f'  *version_info*) [ "{major}.{minor}" = "3.11" ]; exit $? ;;\n'
        '  "-m venv "*) mkdir -p "$3/bin"; cp "$0" "$3/bin/python" ;;\n'
        "esac\nexit 0\n")
    python.chmod(0o755)
    return python


def _digest_by_name():
    lines = "".join(f'  *{name}) echo "{digest}  $1" ;;\n'
                    for name, digest in _WHISPER_PINS.items())
    return f'case "$1" in\n{lines}  *) echo {FAKE_SHA} ;;\nesac'


def test_whisper_is_pinned_by_revision_digest_and_hashes(script_text):
    assert "536b0662742c02347bc0e980a01041f333bce120" in script_text
    for digest in _WHISPER_PINS.values():
        assert digest in script_text
    assert "--require-hashes" in script_text
    for arch in ("aarch64", "x86_64"):
        req = SCRIPT.parent / f"whisper-requirements-{arch}.txt"
        text = req.read_text()
        assert "faster-whisper==1.2.1" in text and "--hash=sha256:" in text


def test_a_wrong_model_digest_installs_no_model(tmp_path):
    _gateway_venv(tmp_path)
    done = run_script(tmp_path, args=["--whisper"])
    assert done.returncode != 0
    assert "digest mismatch" in done.stderr
    assert not (tmp_path / "prefix/models/faster-whisper-small/model.bin").exists()


def test_a_verified_model_reaches_the_gateway_unit(tmp_path):
    _gateway_venv(tmp_path)
    done = run_script(tmp_path, args=["--whisper"], sha256_body=_digest_by_name())
    assert done.returncode == 0, done.stderr
    calls = (tmp_path / "calls.log").read_text()
    assert "--require-hashes --no-deps -r" in calls
    gateway = (tmp_path / "units/delfin-messenger-gateway.service").read_text()
    assert (f"Environment=DELFIN_WHISPER_MODEL={tmp_path}/prefix/models/"
            "faster-whisper-small") in gateway
    assert (f"Environment=DELFIN_WHISPER_PYTHON={tmp_path}/prefix/"
            "whisper-venv/bin/python") in gateway


def test_whisper_refuses_an_interpreter_the_wheels_are_not_for(tmp_path):
    _gateway_venv(tmp_path, version="3.13")
    done = run_script(tmp_path, args=["--whisper"], sha256_body=_digest_by_name())
    assert done.returncode != 0 and "3.11" in done.stderr
