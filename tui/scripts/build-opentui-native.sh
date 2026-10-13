#!/usr/bin/env bash
# Build OpenTUI's native library from source with DELFIN's patches applied
# and install it in place of the prebuilt one under node_modules.
#
# Why: opentui#1630. When the split footer grows, the prebuilt library
# scrolls with CSI n S, and xterm.js -- the terminal in VS Code, JupyterLab
# and most browser-based HPC shells -- deletes the rows that leave the top
# instead of keeping them in scrollback. The patch in patches/ scrolls with
# line feeds from the screen's last row instead, which keeps them. The
# measurements are in docs/TUI_REWRITE.md ("Phase 0 findings").
#
# Usage (from tui/):  scripts/build-opentui-native.sh [--check] [target...]
#   targets: linux-x64 (default on Linux x86_64), linux-arm64
#   --check: only verify that the installed library is the patched build.
# Re-run after `bun install`, which puts the prebuilt library back.
set -euo pipefail

here=$(cd "$(dirname "$0")/.." && pwd)
cache=${DELFIN_TUI_CACHE:-$HOME/.cache/delfin-tui}
version=$(sed -n 's/.*"@opentui\/core": "\([0-9.]*\)".*/\1/p' "$here/package.json")

# The exact source each OpenTUI version is built from: a tag can be moved,
# a commit cannot.
case "$version" in
  0.6.4) commit=19ac2963a41f182122c6e1a88f3c63293de29bbd ;;
  *) echo "no pinned OpenTUI commit for $version; add one to this script"; exit 1 ;;
esac

# Zig 0.16.0, checked against hashes pinned here rather than against an
# index served by the same host as the tarball.
zig_version=0.16.0
zig_sha() {
  case "$1" in
    x86_64-linux)  echo 70e49664a74374b48b51e6f3fdfbf437f6395d42509050588bd49abe52ba3d00 ;;
    aarch64-linux) echo ea4b09bfb22ec6f6c6ceac57ab63efb6b46e17ab08d21f69f3a48b38e1534f17 ;;
    x86_64-macos)  echo 0387557ed1877bc6a2e1802c8391953baddba76081876301c522f52977b52ba7 ;;
    aarch64-macos) echo b23d70deaa879b5c2d486ed3316f7eaa53e84acf6fc9cc747de152450d401489 ;;
    *) return 1 ;;
  esac
}

shopt -s nullglob
patches=("$here"/patches/opentui-v"$version"-*.patch)
shopt -u nullglob
[ ${#patches[@]} -gt 0 ] || { echo "no patches for OpenTUI $version in patches/"; exit 1; }
patch_hash=$(cat "${patches[@]}" | sha256sum | cut -c1-12)

check_only=0
targets=()
for arg in "$@"; do
  case "$arg" in
    --check) check_only=1 ;;
    *) targets+=("$arg") ;;
  esac
done
if [ ${#targets[@]} -eq 0 ]; then
  case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) targets=(linux-x64) ;;
    Linux-aarch64) targets=(linux-arm64) ;;
    *) echo "no default target for $(uname -s)-$(uname -m); name one"; exit 1 ;;
  esac
fi

target_info() {   # target -> "zig-target out-dir lib-name"
  case "$1" in
    linux-x64)   echo "x86_64-linux-gnu.2.17 x86_64-linux libopentui.so" ;;
    linux-arm64) echo "aarch64-linux-gnu.2.17 aarch64-linux libopentui.so" ;;
    *) return 1 ;;
  esac
}

# --- verify only ------------------------------------------------------------
if [ "$check_only" = 1 ]; then
  status=0
  for target in "${targets[@]}"; do
    read -r _zt _out lib < <(target_info "$target") || { echo "unknown target $target"; exit 1; }
    pkg=$here/node_modules/@opentui/core-$target
    built=$cache/build-v$version-$patch_hash-$target/$lib
    if [ -f "$built" ] && [ -f "$pkg/$lib" ] && cmp -s "$built" "$pkg/$lib"; then
      echo "$target: patched build $patch_hash installed"
    else
      echo "$target: NOT the patched build (run scripts/build-opentui-native.sh)"
      status=1
    fi
  done
  exit $status
fi

# --- Zig ----------------------------------------------------------------------
zig_dir=$cache/zig-$zig_version
if [ ! -x "$zig_dir/zig" ]; then
  host="$(uname -m)-$(uname -s | tr 'A-Z' 'a-z')"
  host=${host/darwin/macos}
  host=${host/arm64-macos/aarch64-macos}
  sha=$(zig_sha "$host") || { echo "no pinned Zig hash for $host"; exit 1; }
  mkdir -p "$cache"
  tmp=$(mktemp -d "$cache/zig.XXXXXX")
  curl -fsSL "https://ziglang.org/download/$zig_version/zig-$host-$zig_version.tar.xz" -o "$tmp/zig.tar.xz"
  echo "$sha  $tmp/zig.tar.xz" | sha256sum -c - >/dev/null
  mkdir -p "$tmp/zig" && tar -xf "$tmp/zig.tar.xz" -C "$tmp/zig" --strip-components=1
  rm -rf "$zig_dir" && mv "$tmp/zig" "$zig_dir" && rm -rf "$tmp"
fi
export PATH="$zig_dir:$PATH"

# --- source at the pinned commit, patched, kept only when every patch applied --
src=$cache/opentui-v$version-$patch_hash
if [ ! -f "$src/.delfin-patched" ]; then
  rm -rf "$src"
  tmp=$(mktemp -d "$cache/opentui.XXXXXX")
  git clone -q --branch "v$version" https://github.com/anomalyco/opentui.git "$tmp/src"
  git -C "$tmp/src" -c advice.detachedHead=false checkout -q "$commit"
  [ "$(git -C "$tmp/src" rev-parse HEAD)" = "$commit" ] || { echo "OpenTUI commit mismatch"; exit 1; }
  for p in "${patches[@]}"; do git -C "$tmp/src" apply "$p"; done
  for p in "${patches[@]}"; do git -C "$tmp/src" apply -R --check "$p"; done
  echo "$patch_hash" > "$tmp/src/.delfin-patched"
  mv "$tmp/src" "$src" && rm -rf "$tmp"
fi
(cd "$src/packages/native" && sh scripts/prepare-zig-deps.sh)

for target in "${targets[@]}"; do
  read -r zig_target out lib < <(target_info "$target") || { echo "unknown target $target"; exit 1; }
  pkg=$here/node_modules/@opentui/core-$target
  [ -d "$pkg" ] || { echo "skip $target: $pkg not installed"; continue; }
  built=$cache/build-v$version-$patch_hash-$target/$lib
  if [ ! -f "$built" ]; then
    echo "building OpenTUI $version ($target, patches $patch_hash); a clean build takes minutes"
    # The renderer's own split-footer tests, including the one the patch
    # updates, must pass before anything is installed.
    (cd "$src/packages/native" && zig build test -Dtest-filter="split" --summary none)
    (cd "$src/packages/native" && zig build -Doptimize=ReleaseFast -Dlibrary-target="$zig_target")
    mkdir -p "$(dirname "$built")"
    strip --strip-unneeded -o "$built.tmp" "$src/packages/native/lib/$out/$lib"
    mv -f "$built.tmp" "$built"
  fi
  # Bun hardlinks node_modules files to its global cache, so writing into
  # the installed file would change the cache -- and every project on the
  # machine. Copy to a new file and rename it over the old name instead.
  [ -f "$pkg/$lib.prebuilt" ] || cp "$pkg/$lib" "$pkg/$lib.prebuilt"
  cp "$built" "$pkg/$lib.delfin-new"
  mv -f "$pkg/$lib.delfin-new" "$pkg/$lib"
  printf '%s %s\n' "$patch_hash" "$(sha256sum "$built" | cut -d' ' -f1)" > "$pkg/.delfin-patched"
  echo "installed patched $lib for $target ($patch_hash)"
done
