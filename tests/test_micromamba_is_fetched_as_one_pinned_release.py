"""Every micromamba bootstrap fetches one pinned release and checks it.

Input: delfin/installers/fetch_micromamba.sh (the pinned version, one
SHA-256 per platform, and the fetch-and-verify routine), the four
installers that bootstrap micromamba, the dashboard's Install micromamba
button, the staging of the installers into the tool directories, and the
CI workflow's own micromamba step.

Output:
- the fetcher downloads the release asset named by version and platform,
  keeps it only when its SHA-256 equals the pinned digest, and removes it
  otherwise; it exits non-zero with both digests named and leaves no file
  at the destination;
- a platform with no pinned digest, or a machine with no SHA-256 tool, is
  refused before anything is installed;
- an explicit MICROMAMBA_URL / MICROMAMBA_SHA256 pair (a mirror) is
  honoured, and the digest still decides;
- no installer names the moving ``micro.mamba.pm/.../latest`` endpoint,
  every one of them reaches the shared fetcher, and no installer pipes a
  download straight into tar or a shell;
- the CI workflow's micromamba pin is the same release and digest as the
  installers' linux-64 entry, so the two cannot drift apart.

Driven with a stand-in ``curl`` and ``uname`` on PATH inside tmp_path; no
network is used.
"""

from __future__ import annotations

import hashlib
import os
import re
import stat
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
FETCHER = REPO / "delfin" / "installers" / "fetch_micromamba.sh"
INSTALLERS = (
    REPO / "delfin" / "installers" / "install_delfin.sh",
    REPO / "delfin" / "qm_tools" / "install_qm_tools.sh",
    REPO / "delfin" / "ai_tools" / "install_ai_tools.sh",
    REPO / "delfin" / "analysis_tools" / "install_analysis_tools.sh",
)
SETTINGS_TAB = REPO / "delfin" / "dashboard" / "tab_settings.py"
CI = REPO / ".github" / "workflows" / "ci.yml"
PLATFORMS = ("linux_64", "linux_aarch64", "osx_64", "osx_arm64")


def _pins() -> dict[str, str]:
    text = FETCHER.read_text(encoding="utf-8")
    return dict(re.findall(r'^DELFIN_MICROMAMBA_([A-Za-z0-9_]+)="([^"]*)"$', text, re.M))


def _fake_tools(tmp_path: Path, payload: bytes, system: str = "Linux",
                machine: str = "x86_64", with_sha_tools: bool = True) -> Path:
    """A PATH holding a recording curl, a fixed uname, and the coreutils the
    fetcher needs -- with or without every SHA-256 tool."""
    fake = tmp_path / "fakebin"
    fake.mkdir()
    (tmp_path / "payload.bin").write_bytes(payload)
    (fake / "curl").write_text(
        "#!/bin/sh\n"
        f'printf "%s\\n" "$*" >> "{tmp_path}/curl.log"\n'
        'out=""\n'
        'while [ $# -gt 0 ]; do\n'
        '  case "$1" in -o) out="$2"; shift 2 ;; *) shift ;; esac\n'
        'done\n'
        f'cat "{tmp_path}/payload.bin" > "$out"\n',
        encoding="utf-8")
    (fake / "uname").write_text(
        "#!/bin/sh\n"
        f'case "$1" in -s) echo {system} ;; -m) echo {machine} ;; *) echo {system} ;; esac\n',
        encoding="utf-8")
    for name in ("curl", "uname"):
        (fake / name).chmod(0o755)
    needed = ["bash", "sh", "cat", "mktemp", "chmod", "mv", "rm", "mkdir", "dirname",
              "awk", "printf", "tr", "cut", "head", "basename"]
    if with_sha_tools:
        needed += ["sha256sum", "shasum", "python3", "openssl"]
    for name in needed:
        for d in os.environ["PATH"].split(os.pathsep):
            cand = Path(d) / name
            if cand.is_file() and os.access(cand, os.X_OK):
                (fake / name).symlink_to(cand)
                break
    return fake


def _run(tmp_path: Path, fake: Path, dest: Path, **env) -> subprocess.CompletedProcess:
    full_env = {"PATH": str(fake), "HOME": str(tmp_path)}
    full_env.update(env)
    return subprocess.run([str(fake / "bash"), str(FETCHER), str(dest)],
                          capture_output=True, text=True, timeout=60, env=full_env)


def test_one_release_and_a_digest_for_every_supported_platform():
    pins = _pins()
    assert re.fullmatch(r"\d+\.\d+\.\d+-\d+", pins.get("VERSION", "")), pins
    for platform in PLATFORMS:
        assert re.fullmatch(r"[0-9a-f]{64}", pins.get("SHA256_" + platform, "")), platform


def test_a_download_that_does_not_match_the_pin_is_removed_and_refused(tmp_path):
    payload = b"not the published micromamba\n"
    fake = _fake_tools(tmp_path, payload)
    dest = tmp_path / "bin" / "micromamba"
    done = _run(tmp_path, fake, dest)
    assert done.returncode != 0
    assert not dest.exists()
    if dest.parent.exists():
        assert list(dest.parent.iterdir()) == []  # no partial download either
    pins = _pins()
    assert pins["SHA256_linux_64"] in done.stderr
    assert hashlib.sha256(payload).hexdigest() in done.stderr
    # The release asset of the pinned version, for this platform.
    asked = (tmp_path / "curl.log").read_text()
    assert f"/download/{pins['VERSION']}/micromamba-linux-64" in asked


@pytest.mark.parametrize("system,machine,asset", [
    ("Linux", "x86_64", "linux-64"),
    ("Linux", "aarch64", "linux-aarch64"),
    ("Darwin", "x86_64", "osx-64"),
    ("Darwin", "arm64", "osx-arm64"),
])
def test_each_platform_asks_for_its_own_asset(tmp_path, system, machine, asset):
    fake = _fake_tools(tmp_path, b"x", system=system, machine=machine)
    done = _run(tmp_path, fake, tmp_path / "micromamba")
    assert done.returncode != 0  # the stand-in bytes match no pin
    assert f"micromamba-{asset}" in (tmp_path / "curl.log").read_text()


def test_a_platform_without_a_pin_is_refused_before_any_download(tmp_path):
    fake = _fake_tools(tmp_path, b"x", system="Linux", machine="sparc64")
    dest = tmp_path / "micromamba"
    done = _run(tmp_path, fake, dest)
    assert done.returncode != 0
    assert "sparc64" in done.stderr
    assert not (tmp_path / "curl.log").exists()
    assert not dest.exists()


def test_no_sha256_tool_means_nothing_is_installed(tmp_path):
    fake = _fake_tools(tmp_path, b"x", with_sha_tools=False)
    dest = tmp_path / "micromamba"
    done = _run(tmp_path, fake, dest)
    assert done.returncode != 0
    assert not dest.exists()
    assert "SHA-256" in done.stderr


def test_a_mirror_with_its_digest_is_accepted_and_made_executable(tmp_path):
    payload = b"#!/bin/sh\necho mirrored\n"
    fake = _fake_tools(tmp_path, payload)
    dest = tmp_path / "bin" / "micromamba"
    done = _run(tmp_path, fake, dest,
                MICROMAMBA_URL="https://mirror.invalid/micromamba",
                MICROMAMBA_SHA256=hashlib.sha256(payload).hexdigest())
    assert done.returncode == 0, done.stderr
    assert dest.read_bytes() == payload
    assert dest.stat().st_mode & stat.S_IXUSR
    assert "https://mirror.invalid/micromamba" in (tmp_path / "curl.log").read_text()
    # Nothing beside it: the partial download was moved, not copied.
    assert sorted(p.name for p in dest.parent.iterdir()) == ["micromamba"]


def test_a_mirror_without_a_matching_digest_is_still_refused(tmp_path):
    fake = _fake_tools(tmp_path, b"other bytes")
    dest = tmp_path / "micromamba"
    done = _run(tmp_path, fake, dest, MICROMAMBA_URL="https://mirror.invalid/micromamba")
    assert done.returncode != 0
    assert not dest.exists()


def test_no_installer_names_the_moving_endpoint_and_all_reach_the_fetcher():
    for path in INSTALLERS:
        text = path.read_text(encoding="utf-8")
        assert "micro.mamba.pm" not in text, path
        assert "fetch_micromamba.sh" in text, path
    tab = SETTINGS_TAB.read_text(encoding="utf-8")
    assert "micro.mamba.pm" not in tab
    assert "fetch_micromamba.sh" in tab
    assert "--insecure" not in tab


_PIPE_INTO = re.compile(r"\b(curl|wget)\b[^\n|]*\|\s*(tar|sh|bash|zsh|python3?)\b")


def test_no_installer_pipes_a_download_into_tar_or_a_shell():
    files = [p for p in (REPO / "delfin").rglob("*.sh")]
    files += [p for p in (REPO / "scripts").rglob("*.sh")]
    files += [p for p in (REPO / "tui" / "scripts").rglob("*.sh")]
    files += [SETTINGS_TAB, REPO / "delfin" / "runtime_setup.py", REPO / "delfin" / "installer.py"]
    hits = []
    for path in files:
        for n, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if line.lstrip().startswith("#"):
                continue
            if _PIPE_INTO.search(line):
                hits.append(f"{path.relative_to(REPO)}:{n}: {line.strip()}")
    assert not hits, hits


def test_the_ci_pin_is_the_installers_pin():
    pins = _pins()
    text = CI.read_text(encoding="utf-8")
    urls = re.findall(r"MICROMAMBA_URL: (\S+)", text)
    digests = re.findall(r"MICROMAMBA_SHA256: ([0-9a-f]{64})", text)
    assert urls and digests
    for url in urls:
        assert url.endswith(f"/download/{pins['VERSION']}/micromamba-linux-64"), url
    assert set(digests) == {pins["SHA256_linux_64"]}


def test_a_staged_installer_finds_the_fetcher(tmp_path):
    from delfin import runtime_setup

    for stage, script in ((runtime_setup.stage_packaged_qm_tools, "install_qm_tools.sh"),
                          (runtime_setup.stage_packaged_ai_tools, "install_ai_tools.sh"),
                          (runtime_setup.stage_packaged_analysis_tools,
                           "install_analysis_tools.sh")):
        target = stage(tmp_path / script.split(".")[0])
        assert (target / script).is_file()
        staged = target / "fetch_micromamba.sh"
        assert staged.read_bytes() == FETCHER.read_bytes(), script


def test_the_settings_button_runs_the_shared_fetcher():
    """The handler hands the destination to the fetcher and nothing else."""
    tab = SETTINGS_TAB.read_text(encoding="utf-8")
    body = tab.split("def _on_install_micromamba(button):", 1)[1].split("\n    def ", 1)[0]
    assert "fetch_micromamba.sh" in body
    assert "curl" not in body and "tar" not in body
