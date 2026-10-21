"""xtb4stda, stda and their parameter files are installed only as published.

Input: delfin/qm_tools/install_qm_tools.sh -- its fetch_verified and
sha256_of functions, driven in bash against file:// sources in tmp_path,
and the xtb4stda bundle's download list.

Output:
- a file whose SHA-256 matches is kept and reported verified;
- a fresh download that does not match exits 1, names both digests, and
  leaves no file behind;
- a kept download that does not match (a truncated earlier run) is fetched
  once more and then accepted when the new copy matches;
- every file of the xtb4stda bundle goes through fetch_verified with a
  digest of its own, and nothing is fetched from a moving branch.

The full install was driven once by hand against local copies of the real
files (17 verified, xtb4stda installed) and once with one byte of xtb4stda
changed (exit 1, "Not installed: xtb4stda", bin/ and downloads/ empty).
"""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "delfin" / "qm_tools" / "install_qm_tools.sh"


def _function(name: str) -> str:
    text = SCRIPT.read_text(encoding="utf-8")
    return name + "() {" + text.split("\n" + name + "() {", 1)[1].split("\n}\n", 1)[0] + "\n}\n"


def _fetch(tmp_path: Path, source: Path, dest: Path, want: str) -> subprocess.CompletedProcess:
    body = ("set -euo pipefail\nFORCE_REDOWNLOAD=0\n"
            + "".join(_function(n) for n in ("log", "warn", "die", "have", "download_file",
                                              "sha256_of", "fetch_verified"))
            + f'fetch_verified "file://{source}" "{dest}" "{want}"\n')
    return subprocess.run(["bash", "-c", body], capture_output=True, text=True, timeout=60,
                          env={"PATH": os.environ["PATH"], "HOME": str(tmp_path)})


def _source(tmp_path: Path, data: bytes) -> tuple[Path, str]:
    src = tmp_path / "upstream.bin"
    src.write_bytes(data)
    return src, hashlib.sha256(data).hexdigest()


def test_a_matching_download_is_kept(tmp_path):
    src, digest = _source(tmp_path, b"published bytes\n")
    dest = tmp_path / "dl" / "tool"
    dest.parent.mkdir()
    done = _fetch(tmp_path, src, dest, digest)
    assert done.returncode == 0, done.stderr
    assert dest.read_bytes() == b"published bytes\n"
    assert "checksum verified: tool" in done.stdout


def test_a_download_that_does_not_match_is_removed_and_refused(tmp_path):
    src, digest = _source(tmp_path, b"published bytes\n")
    src.write_bytes(b"replaced bytes\n")
    dest = tmp_path / "dl" / "tool"
    dest.parent.mkdir()
    done = _fetch(tmp_path, src, dest, digest)
    assert done.returncode == 1
    assert not dest.exists()
    assert f"expected {digest}" in done.stderr
    assert hashlib.sha256(b"replaced bytes\n").hexdigest() in done.stderr


def test_a_kept_download_that_does_not_match_is_fetched_again(tmp_path):
    src, digest = _source(tmp_path, b"published bytes\n")
    dest = tmp_path / "dl" / "tool"
    dest.parent.mkdir()
    dest.write_bytes(b"publ")  # an interrupted earlier run
    done = _fetch(tmp_path, src, dest, digest)
    assert done.returncode == 0, done.stderr
    assert "fetching it again" in done.stdout
    assert dest.read_bytes() == b"published bytes\n"


def test_every_file_of_the_bundle_has_a_digest_and_a_fixed_source():
    text = SCRIPT.read_text(encoding="utf-8")
    bundle = _function("install_xtb4stda_bundle")
    assert "download_file" not in bundle, "a bundle file bypasses the checksum"
    assert bundle.count("fetch_verified") == 3

    names = re.findall(r'^\s+"(\.[^"]+)"$', bundle.split("runtime_files=(", 1)[1].split(")", 1)[0],
                       re.M)
    assert len(names) == 15, names
    body = _function("xtb4stda_runtime_sha256") + "".join(
        f'xtb4stda_runtime_sha256 "{n}"\n' for n in names)
    digests = subprocess.run(["bash", "-c", "set -e\n" + body], capture_output=True,
                             text=True, timeout=30, check=True).stdout.split()
    assert len(digests) == 15 and all(re.fullmatch(r"[0-9a-f]{64}", d) for d in digests)

    for var in ("XTB4STDA_SHA256", "STDA_SHA256"):
        assert re.search(var + r'="\$\{' + var + r':-[0-9a-f]{64}\}"', text), var
    assert re.search(r'XTB4STDA_REF="\$\{XTB4STDA_REF:-[0-9a-f]{40}\}"', text)
    base = re.search(r'^XTB4STDA_RUNTIME_BASE_URL="\$\{XTB4STDA_RUNTIME_BASE_URL:-(.+)\}"$',
                     text, re.M).group(1)
    assert base.endswith("/${XTB4STDA_REF}"), base
