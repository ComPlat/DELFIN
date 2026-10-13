#!/usr/bin/env bash
# Build OpenTUI's native library from source with DELFIN's patches applied
# and install it in place of the prebuilt one under node_modules.
#
# Why: opentui#1630. When the split footer grows, the prebuilt library
# scrolls with CSI n S, and xterm.js -- the terminal in VS Code, JupyterLab
# and most browser-based HPC shells -- deletes the rows that leave the top
# instead of keeping them in scrollback. patches/ makes it scroll with line
# feeds, which every terminal keeps. Measured in docs/TUI_REWRITE.md.
#
# Usage (from tui/): scripts/build-opentui-native.sh [target...]
#   targets: linux-x64 (default on Linux x86_64), linux-arm64
# Re-run after `bun install`, which restores the prebuilt library.
set -euo pipefail

here=$(cd "$(dirname "$0")/.." && pwd)
cache=${DELFIN_TUI_CACHE:-$HOME/.cache/delfin-tui}
version=$(sed -n 's/.*"@opentui\/core": "\([0-9.]*\)".*/\1/p' "$here/package.json")
zig_version=0.16.0
patches=("$here"/patches/opentui-v"$version"-*.patch)
[ -e "${patches[0]}" ] || { echo "no patches for OpenTUI $version in patches/"; exit 1; }
patch_hash=$(cat "${patches[@]}" | sha256sum | cut -c1-12)

targets=("$@")
if [ ${#targets[@]} -eq 0 ]; then
  case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) targets=(linux-x64) ;;
    Linux-aarch64) targets=(linux-arm64) ;;
    *) echo "no default target for $(uname -s)-$(uname -m); name one"; exit 1 ;;
  esac
fi

# --- Zig, verified against ziglang.org's release index -------------------
zig_dir=$cache/zig-$zig_version
if [ ! -x "$zig_dir/zig" ]; then
  mkdir -p "$cache"
  host="$(uname -m)-$(uname -s | tr 'A-Z' 'a-z')"
  curl -fsS https://ziglang.org/download/index.json -o "$cache/zig-index.json"
  read -r url sha < <(python3 -I -c 'import json,sys; e=json.load(open(sys.argv[1]))[sys.argv[2]][sys.argv[3]]; print(e["tarball"], e["shasum"])' \
    "$cache/zig-index.json" "$zig_version" "$host")
  curl -fsSL "$url" -o "$cache/zig.tar.xz"
  echo "$sha  $cache/zig.tar.xz" | sha256sum -c - >/dev/null
  mkdir -p "$zig_dir" && tar -xf "$cache/zig.tar.xz" -C "$zig_dir" --strip-components=1
  rm "$cache/zig.tar.xz"
fi
export PATH="$zig_dir:$PATH"

# --- OpenTUI source at the pinned tag, with the patches applied ----------
src=$cache/opentui-v$version-$patch_hash
if [ ! -d "$src" ]; then
  git clone -q --depth 1 --branch "v$version" https://github.com/anomalyco/opentui.git "$src"
  for p in "${patches[@]}"; do git -C "$src" apply "$p"; done
fi
(cd "$src/packages/native" && sh scripts/prepare-zig-deps.sh)

for target in "${targets[@]}"; do
  case "$target" in
    linux-x64)   zig_target=x86_64-linux-gnu.2.17;  out=x86_64-linux;  lib=libopentui.so ;;
    linux-arm64) zig_target=aarch64-linux-gnu.2.17; out=aarch64-linux; lib=libopentui.so ;;
    *) echo "unsupported target: $target"; exit 1 ;;
  esac
  pkg=$here/node_modules/@opentui/core-$target
  [ -d "$pkg" ] || { echo "skip $target: $pkg not installed"; continue; }
  built=$cache/build-v$version-$patch_hash-$target/$lib
  if [ ! -f "$built" ]; then
    echo "building OpenTUI $version ($target, patches $patch_hash); a clean build takes minutes"
    (cd "$src/packages/native" && zig build -Doptimize=ReleaseFast -Dlibrary-target="$zig_target")
    mkdir -p "$(dirname "$built")"
    strip --strip-unneeded -o "$built" "$src/packages/native/lib/$out/$lib"
  fi
  [ -f "$pkg/$lib.prebuilt" ] || cp "$pkg/$lib" "$pkg/$lib.prebuilt"
  cp "$built" "$pkg/$lib"
  printf '%s\n' "$patch_hash" > "$pkg/.delfin-patched"
  echo "installed patched $lib for $target ($patch_hash)"
done
