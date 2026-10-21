#!/usr/bin/env bash
# Install a pinned Bun release: the archive from Bun's GitHub release,
# checked against a SHA-256 pinned here, into $HOME/.bun/bin -- and onto
# $GITHUB_PATH in a workflow.
#
# Why not oven-sh/setup-bun: a repository may allow GitHub's own actions
# only (DELFIN's other workflows use nothing else), and a workflow naming
# another one then fails before its first step ("startup failure"). One
# checked download needs no action at all.
#
# Usage:  tui/scripts/install-bun.sh [version]    (default 1.4.3)
set -euo pipefail

version=${1:-1.4.3}
case "$(uname -s)-$(uname -m)" in
  Linux-x86_64)  asset=bun-linux-x64 ;;
  Linux-aarch64) asset=bun-linux-aarch64 ;;
  Darwin-arm64)  asset=bun-darwin-aarch64 ;;
  Darwin-x86_64) asset=bun-darwin-x64 ;;
  *) echo "no Bun release for $(uname -s) $(uname -m)" >&2; exit 1 ;;
esac

# From the release's SHASUMS256.txt, checked against GitHub's digest of
# each asset and pinned here rather than fetched from the host that serves
# the archive.
pinned() {
  case "$version/$1" in
    1.4.3/bun-linux-x64)      echo 3e730528d89775f03f87f05c06787c4df4099d0a37f3f252e7b829cb23816e18 ;;
    1.4.3/bun-linux-aarch64)  echo efa9813da5ed72423bf847f916e8d2c47c0d776add972354026a75e10da9aa21 ;;
    1.4.3/bun-darwin-aarch64) echo 80afd4c069b46a6fa8f3ec520bebc750a7843e4a8d5b9583d29ad436c00a3403 ;;
    1.4.3/bun-darwin-x64)     echo 62660e6438eb872278877acca77e1152d5d031324df368d4a81b428211489d2d ;;
    *) return 1 ;;
  esac
}
want=$(pinned "$asset") || { echo "no pinned hash for Bun $version ($asset); add one to this script" >&2; exit 1; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
curl -fsSL --retry 3 -o "$tmp/$asset.zip" \
  "https://github.com/oven-sh/bun/releases/download/bun-v$version/$asset.zip"
if command -v sha256sum >/dev/null; then
  got=$(sha256sum "$tmp/$asset.zip" | cut -d' ' -f1)
else
  got=$(shasum -a 256 "$tmp/$asset.zip" | cut -d' ' -f1)      # macOS
fi
if [ "$got" != "$want" ]; then
  echo "Bun $version ($asset.zip): sha256 $got, expected $want" >&2
  exit 1
fi

unzip -q "$tmp/$asset.zip" -d "$tmp"
bin="$HOME/.bun/bin"
mkdir -p "$bin"
install -m 0755 "$tmp/$asset/bun" "$bin/bun"
ln -sf bun "$bin/bunx"
if [ -n "${GITHUB_PATH:-}" ]; then
  echo "$bin" >> "$GITHUB_PATH"
fi
echo "Bun $("$bin/bun" --version) in $bin"
