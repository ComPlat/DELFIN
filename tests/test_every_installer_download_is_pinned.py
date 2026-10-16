"""What the installers download beyond micromamba is named exactly and checked.

Input: delfin/qm_tools/install_qm_tools.sh (g-xTB, std2),
delfin/installers/install_delfin.sh (the OpenMPI source),
delfin/analysis_tools/install_analysis_tools.sh (CENSO, its helper
scripts, ANMR), delfin/csp_tools/install_csp_tools.sh (Genarris), and
every installer script under delfin/ for moving references.

Output:
- g-xTB and the std2 source go through fetch_verified with a SHA-256
  written in the installer; there is no path that installs them when no
  SHA-256 tool exists;
- the OpenMPI source tarball has a published SHA-256 per release; a
  download that does not match is removed and not built, and a release
  with no digest is refused before anything is downloaded unless
  OPENMPI_SHA256 names one;
- CENSO is installed from a commit, not a branch; its helper scripts and
  the ANMR binary are fetched from fixed commits/releases and each held
  to its SHA-256; nothing is resolved through ``releases/latest``;
- Genarris is checked out at a commit and the checkout is compared with
  it;
- no installer fetches a branch head (``refs/heads/``) or a ``latest``
  release.

Driven in bash with stand-in downloads in tmp_path; no network.
"""

from __future__ import annotations

import hashlib
import os
import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
QM = REPO / "delfin" / "qm_tools" / "install_qm_tools.sh"
DELFIN_SH = REPO / "delfin" / "installers" / "install_delfin.sh"
ANALYSIS = REPO / "delfin" / "analysis_tools" / "install_analysis_tools.sh"
CSP = REPO / "delfin" / "csp_tools" / "install_csp_tools.sh"
HEX64 = r"[0-9a-f]{64}"
HEX40 = r"[0-9a-f]{40}"


def _function(path: Path, name: str) -> str:
    text = path.read_text(encoding="utf-8")
    return name + "() {" + text.split("\n" + name + "() {", 1)[1].split("\n}\n", 1)[0] + "\n}\n"


def _one_liners(path: Path) -> str:
    return "".join(line + "\n" for line in path.read_text(encoding="utf-8").splitlines()
                   if re.match(r"^(log|warn|die|have)\(\)\s+\{ .* \}$", line))


# --- g-xTB and std2 -------------------------------------------------------------

def test_gxtb_is_held_to_a_pinned_digest_and_never_installed_unchecked():
    text = QM.read_text(encoding="utf-8")
    assert re.search(r'GXTB_SHA256="\$\{GXTB_SHA256:-\$\(gxtb_sha256 "\$\{GXTB_ASSET\}"\)\}"',
                     text) or re.search(r"gxtb_sha256\(\) \{", text)
    body = _function(QM, "install_gxtb")
    assert "fetch_verified" in body
    assert "download_file" not in body
    assert "could not be verified" not in body
    # The linux asset the installer fetches by default has its digest here.
    digests = subprocess.run(
        ["bash", "-c", _function(QM, "gxtb_sha256")
         + 'gxtb_sha256 "xtb-6.7.1-gxtb-140526-linux-x86_64.tar.xz"'],
        capture_output=True, text=True, timeout=30, check=True).stdout.split()
    assert digests == ["cdfb81164f8d3300fc451dbbf2414528b533e71f7a30ec824b564cd30d6d9e5c"]


def test_the_std2_source_is_held_to_a_pinned_digest():
    text = QM.read_text(encoding="utf-8")
    assert re.search(r'STD2_SHA256="\$\{STD2_SHA256:-' + HEX64 + r'\}"', text)
    body = _function(QM, "build_std2_from_source")
    assert 'fetch_verified "${STD2_SRC_URL}"' in body
    assert "download_file" not in body


# --- OpenMPI --------------------------------------------------------------------

def _openmpi_probe(tmp_path: Path, version: str, served: bytes, **env) -> subprocess.CompletedProcess:
    src = tmp_path / "served.tar.gz"
    src.write_bytes(served)
    prefix = tmp_path / "prefix"
    body = ("set -uo pipefail\n" + _one_liners(DELFIN_SH)
            + f'PREFIX="{prefix}"\nDRY_RUN=0\nMAKE_JOBS=1\n'
            + f'fetch() {{ printf "%s\\n" "$1" >> "{tmp_path}/fetched.log"; cp "{src}" "$2"; }}\n'
            + "".join(_function(DELFIN_SH, n) for n in ("sha256_of", "openmpi_sha256",
                                                         "build_openmpi"))
            + f'build_openmpi "{version}"\n')
    full = {"PATH": os.environ["PATH"], "HOME": str(tmp_path)}
    full.update(env)
    return subprocess.run(["bash", "-c", body], capture_output=True, text=True, timeout=60,
                          env=full)


def test_openmpi_releases_have_published_digests():
    table = subprocess.run(
        ["bash", "-c", _function(DELFIN_SH, "openmpi_sha256")
         + "for v in 4.1.1 4.1.6 4.1.8 5.0.8; do openmpi_sha256 $v; done"],
        capture_output=True, text=True, timeout=30, check=True).stdout.split()
    assert table == [
        "d80b9219e80ea1f8bcfe5ad921bd9014285c4948c5965f4156a3831e60776444",
        "44da277b8cdc234e71c62473305a09d63f4dcca292ca40335aab7c4bf0e6a566",
        "fb41086bbed9300baa2f3d7572491facfe5257412fa524ec5a396aa9101d5c62",
        "f891ddf2dab3b604f2521d0569615bc1c07a1ac86a295ac543e059aecb303621",
    ]


def test_an_openmpi_tarball_that_does_not_match_is_removed_and_not_built(tmp_path):
    done = _openmpi_probe(tmp_path, "4.1.8", b"not the release\n")
    assert done.returncode != 0
    assert "fb41086bbed9300baa2f3d7572491facfe5257412fa524ec5a396aa9101d5c62" in done.stderr
    assert not list((tmp_path / "prefix" / ".build").glob("openmpi-4.1.8*"))
    assert not (tmp_path / "prefix" / "openmpi-4.1.8").exists()


def test_an_openmpi_release_without_a_digest_is_refused_before_download(tmp_path):
    done = _openmpi_probe(tmp_path, "9.9.9", b"x")
    assert done.returncode != 0
    assert "OPENMPI_SHA256" in done.stderr
    assert not (tmp_path / "fetched.log").exists()


def test_an_openmpi_digest_given_by_the_user_is_what_decides(tmp_path):
    served = b"a mirror's tarball\n"
    done = _openmpi_probe(tmp_path, "9.9.9", served,
                          OPENMPI_SHA256=hashlib.sha256(b"something else").hexdigest())
    assert done.returncode != 0
    assert (tmp_path / "fetched.log").exists()
    assert not list((tmp_path / "prefix" / ".build").glob("openmpi-9.9.9.tar.gz"))


# --- CENSO, its helper scripts, ANMR ----------------------------------------------

def test_censo_is_installed_from_commits_not_a_branch():
    text = ANALYSIS.read_text(encoding="utf-8")
    for name in ("desired_censo_git_ref", "compatible_censo_git_ref"):
        body = _function(ANALYSIS, name)
        assert '"main"' not in body, name
    assert re.search(r'CENSO_V3_REF="\$\{CENSO_V3_REF:-' + HEX40 + r'\}"', text)
    assert re.search(r'CENSO_V2_REF="\$\{CENSO_V2_REF:-' + HEX40 + r'\}"', text)
    assert "git+https://github.com/grimme-lab/CENSO.git\"" not in text
    assert "refs/heads/" not in text


def test_the_censo_helper_scripts_and_anmr_are_held_to_their_digests():
    text = ANALYSIS.read_text(encoding="utf-8")
    assert "releases/latest" not in text
    assert "api.github.com" not in text
    assert re.search(r'ENSO_ANMR_SHA256="\$\{ENSO_ANMR_SHA256:-' + HEX64 + r'\}"', text)
    anmr = _function(ANALYSIS, "install_anmr")
    assert "fetch_verified" in anmr and "download_to_file" not in anmr
    helpers = _function(ANALYSIS, "install_censo_helper_wrappers_from_repo")
    assert "fetch_verified" in helpers and "download_to_file" not in helpers
    digests = subprocess.run(
        ["bash", "-c", _function(ANALYSIS, "censo_helper_sha256")
         + 'for r in "$1" "$2"; do for h in c2anmr nmrplot; do censo_helper_sha256 "$r" "$h"; done; done',
         "_", "1a82c5ad5e0b46c17cbbc93ecd7890a3b8c8a1e3", "db667bde403c916bd20e8a42544fa3b0bed8f7d7"],
        capture_output=True, text=True, timeout=30, check=True).stdout.split()
    assert len(digests) == 4 and all(re.fullmatch(HEX64, d) for d in digests)


def _analysis_fetch(tmp_path: Path, served: bytes, want: str) -> subprocess.CompletedProcess:
    src = tmp_path / "served"
    src.write_bytes(served)
    dest = tmp_path / "out" / "file"
    body = ("set -uo pipefail\n" + _one_liners(ANALYSIS)
            + "".join(_function(ANALYSIS, n) for n in ("log", "warn", "have", "sha256_of",
                                                        "download_to_file", "fetch_verified")
                      if ("\n" + n + "() {") in ANALYSIS.read_text(encoding="utf-8"))
            + f'mkdir -p "{dest.parent}"\nfetch_verified "file://{src}" "{dest}" "{want}"\n')
    return subprocess.run(["bash", "-c", body], capture_output=True, text=True, timeout=60,
                          env={"PATH": os.environ["PATH"], "HOME": str(tmp_path)})


def test_the_analysis_fetch_removes_a_mismatch(tmp_path):
    refused = _analysis_fetch(tmp_path, b"changed\n", hashlib.sha256(b"published\n").hexdigest())
    assert refused.returncode != 0
    assert not (tmp_path / "out" / "file").exists()
    kept = _analysis_fetch(tmp_path, b"published\n", hashlib.sha256(b"published\n").hexdigest())
    assert kept.returncode == 0, kept.stderr
    assert (tmp_path / "out" / "file").read_bytes() == b"published\n"


# --- Genarris -------------------------------------------------------------------

def test_genarris_is_checked_out_at_a_commit_and_compared_with_it():
    text = CSP.read_text(encoding="utf-8")
    assert re.search(r'GENARRIS_REF="\$\{GENARRIS_REF:-' + HEX40 + r'\}"', text)
    body = _function(CSP, "install_genarris")
    assert 'checkout -q "${GENARRIS_REF}"' in body or 'checkout "${GENARRIS_REF}"' in body
    assert 'rev-parse HEAD' in body
    assert "pull --ff-only" not in body


# --- no moving references anywhere ------------------------------------------------

@pytest.mark.parametrize("path", sorted((REPO / "delfin").rglob("install_*.sh")),
                         ids=lambda p: str(p.relative_to(REPO)))
def test_no_installer_fetches_a_branch_head_or_a_latest_release(path):
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if line.lstrip().startswith("#"):
            continue
        assert "refs/heads/" not in line, f"{path}:{n}"
        assert not re.search(r"releases/latest|/latest\b", line), f"{path}:{n}: {line}"
