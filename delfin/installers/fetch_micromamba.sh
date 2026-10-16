#!/usr/bin/env bash
# ============================================================================
# micromamba for DELFIN's installers: one pinned release, held to its SHA-256.
#
#   bash fetch_micromamba.sh DEST        download, verify, install DEST (0755)
#
# Every installer that needs micromamba and finds none runs this file:
# install_delfin.sh, install_qm_tools.sh, install_ai_tools.sh,
# install_analysis_tools.sh, and the dashboard's Install micromamba button
# (delfin/dashboard/tab_settings.py). It is the only place the release and
# its digests are written; the CI workflow names the same linux-64 pin
# (tests/test_micromamba_is_fetched_as_one_pinned_release.py holds the two
# together).
#
# Source: the standalone binary of the mamba-org/micromamba-releases
# release, not the tarball. Each digest below was read on 2026-10-16 from
# the release's own micromamba-<platform>.sha256 file and agrees with the
# SHA-256 GitHub records for the asset. The binary is ~18 MB against ~7 MB
# for the .tar.bz2, and needs no tar or bzip2 to unpack.
#
# Rejected: micro.mamba.pm/api/micromamba/<platform>/latest, which these
# installers fetched before -- a moving target with no digest to compare,
# piped straight into tar.
#
# A mirror: MICROMAMBA_URL names another location; the pinned digest of
# this platform still applies, or MICROMAMBA_SHA256 names the digest of
# the file the mirror serves. Nothing that fails the check is installed.
#
# Output: progress and errors on stderr only, so a caller can capture its
# own stdout. Exit 0 with DEST in place, non-zero with DEST untouched.
# ============================================================================

DELFIN_MICROMAMBA_VERSION="2.9.0-0"
DELFIN_MICROMAMBA_SHA256_linux_64="366cd9cd8be14df1ab8ed50352a82111082a36686b2d389fdb79a92c3fafb3e3"
DELFIN_MICROMAMBA_SHA256_linux_aarch64="9f93b974adcb4d166996af969b6cd371287d1a3e52733704727884d9b74cb7a7"
DELFIN_MICROMAMBA_SHA256_osx_64="1e71054bb3ac9a076e21f7ec48acfef536f9b3f1408f371a942784bf5ef83d8a"
DELFIN_MICROMAMBA_SHA256_osx_arm64="ec2a072f028e1a7cf20f3e2e74d5a8127cf5a5f27636375b5359811565f4e5be"

_delfin_mm_say() { printf '[micromamba] %s\n' "$*" >&2; }

# The release's platform name for this machine; status 1 when it has none.
delfin_micromamba_platform() {
  local system machine
  system="$(uname -s 2>/dev/null)"
  machine="$(uname -m 2>/dev/null)"
  case "${system}/${machine}" in
    Linux/x86_64|Linux/amd64)   printf '%s\n' linux-64 ;;
    Linux/aarch64|Linux/arm64)  printf '%s\n' linux-aarch64 ;;
    Darwin/x86_64)              printf '%s\n' osx-64 ;;
    Darwin/arm64|Darwin/aarch64) printf '%s\n' osx-arm64 ;;
    *) printf '%s\n' "${system}/${machine}"; return 1 ;;
  esac
}

# The pinned digest of a platform name (linux-64 ...); status 1 when none.
delfin_micromamba_sha256() {
  local var="DELFIN_MICROMAMBA_SHA256_$(printf '%s' "$1" | tr '-' '_')"
  local value="${!var:-}"
  [ -n "${value}" ] || return 1
  printf '%s\n' "${value}"
}

_delfin_mm_sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif command -v python3 >/dev/null 2>&1; then
    python3 -I -c 'import hashlib, sys; print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())' "$1"
  elif command -v openssl >/dev/null 2>&1; then
    openssl dgst -sha256 -r "$1" | awk '{print $1}'
  else
    return 1
  fi
}

delfin_fetch_micromamba() {
  local dest="${1:-}" platform url want got dir partial
  if [ -z "${dest}" ]; then
    _delfin_mm_say "usage: fetch_micromamba.sh DEST"
    return 2
  fi
  if ! platform="$(delfin_micromamba_platform)"; then
    _delfin_mm_say "no pinned micromamba release for ${platform}; nothing was downloaded." \
      "Install micromamba or conda yourself and put it on PATH."
    return 1
  fi
  want="${MICROMAMBA_SHA256:-}"
  if [ -z "${want}" ] && ! want="$(delfin_micromamba_sha256 "${platform}")"; then
    _delfin_mm_say "no pinned digest for ${platform}; nothing was downloaded"
    return 1
  fi
  url="${MICROMAMBA_URL:-https://github.com/mamba-org/micromamba-releases/releases/download/${DELFIN_MICROMAMBA_VERSION}/micromamba-${platform}}"

  dir="$(dirname "${dest}")"
  mkdir -p "${dir}" || return 1
  partial="$(mktemp "${dir}/.micromamba.partial.XXXXXX")" || return 1
  _delfin_mm_say "fetching micromamba ${DELFIN_MICROMAMBA_VERSION} for ${platform} (about 18 MB)"
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL --retry 3 -o "${partial}" "${url}" || { rm -f "${partial}"; _delfin_mm_say "download failed: ${url}"; return 1; }
  elif command -v wget >/dev/null 2>&1; then
    wget -q -O "${partial}" "${url}" || { rm -f "${partial}"; _delfin_mm_say "download failed: ${url}"; return 1; }
  else
    rm -f "${partial}"
    _delfin_mm_say "neither curl nor wget is available to fetch micromamba"
    return 1
  fi
  if ! got="$(_delfin_mm_sha256_of "${partial}")" || [ -z "${got}" ]; then
    rm -f "${partial}"
    _delfin_mm_say "no sha256sum, shasum, python3 or openssl to check the SHA-256 of the download;" \
      "an unchecked binary is not installed"
    return 1
  fi
  if [ "${got}" != "${want}" ]; then
    rm -f "${partial}"
    _delfin_mm_say "checksum mismatch for ${url}: expected ${want}, got ${got}." \
      "The download was removed and nothing was installed."
    return 1
  fi
  chmod 755 "${partial}" && mv -f "${partial}" "${dest}" || { rm -f "${partial}"; return 1; }
  _delfin_mm_say "checksum verified; micromamba is at ${dest}"
}

# Run as a script: fetch into the destination given.
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  delfin_fetch_micromamba "$@"
fi
