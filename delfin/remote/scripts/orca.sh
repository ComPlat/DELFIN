#!/usr/bin/env bash
# Where DELFIN's installer would find ORCA on this cluster
# (delfin/remote/setup.py, orca_status). Fed over the tunnel as
#
#   bash -s -- <software_dir>
#
# and prints one line, in the order install_delfin.sh's find_orca looks:
#
#   dir <path>       an ORCA directory (orca plus orca_plot or orca_mdci):
#                    the one `orca` on the PATH belongs to, else the newest
#                    orca_*/ under the software dir or the home
#   tarball <path>   no directory, but an ORCA tarball the installer unpacks
#   none             neither: the user puts their ORCA into the home
#
# Read-only: nothing is unpacked, copied or written. /usr/bin/orca alone is
# the GNOME screen reader, not ORCA, and is not taken for it.
set -u

home_path() {
  case "$1" in
    "~") printf '%s' "$HOME" ;;
    "~/"*) printf '%s/%s' "$HOME" "${1#"~/"}" ;;
    *) printf '%s' "$1" ;;
  esac
}

is_orca_dir() {
  [ -x "$1/orca" ] && { [ -x "$1/orca_plot" ] || [ -x "$1/orca_mdci" ]; }
}

prefix="$(home_path "${1:-~/delfin/software}")"

candidate="$(command -v orca 2>/dev/null || true)"
if [ -n "$candidate" ]; then
  candidate="$(dirname "$(readlink -f "$candidate")")"
  if is_orca_dir "$candidate"; then
    printf 'dir %s\n' "$candidate"
    exit 0
  fi
fi
while IFS= read -r candidate; do
  if is_orca_dir "$candidate"; then
    printf 'dir %s\n' "$candidate"
    exit 0
  fi
done < <(ls -d "$prefix"/orca_*/ "$HOME"/orca_*/ 2>/dev/null | sed 's:/$::' | sort -V -r)
tarball="$(ls -1 "$HOME"/orca_*.tar.xz "$prefix"/orca_*.tar.xz 2>/dev/null | sort -V -r | head -n 1)"
if [ -n "$tarball" ]; then
  printf 'tarball %s\n' "$tarball"
  exit 0
fi
printf 'none\n'
