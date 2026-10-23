#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# The Raspberry Pi OS kernel with Landlock, beside the stock kernel.
#
# Raspberry Pi OS builds its kernels without CONFIG_SECURITY_LANDLOCK. The
# agent's command cage is bubblewrap PLUS Landlock (agent.bash_isolation);
# only with both does a command off the auto-allow list run in the cage
# without a question (_caged_command_may_run_unasked: Landlock confines its
# reads). Without Landlock each of them is asked about instead.
#
# This builds the INSTALLED stock kernel again -- the same Debian source
# package version (linux-source, from the signed Raspberry Pi archive) and
# the same /boot/config-<release> -- with exactly two changes:
# CONFIG_SECURITY_LANDLOCK=y and CONFIG_LSM="landlock". Nothing of the stock
# kernel is touched:
#   /boot/firmware/kernel_<flavour>_landlock.img, initramfs_<flavour>_landlock,
#   /lib/modules/<release>-landlock
# are files of their own, which neither apt nor the raspi-firmware hooks
# write (their flavour pattern does not match "-landlock").
#
# The first boot is a tryboot: tryboot.txt selects the new kernel and a
# cmdline with panic=10, so a kernel that does not come up reboots into the
# stock one by itself. A boot unit then checks that Landlock is active and
# only then adds the kernel line to config.txt, inside the board's own
# section ([pi5] or [pi4]): the card moved to another model boots that
# model's stock kernel.
#
# Usage (on the Pi, as root):
#   sudo bash pi_kernel_landlock.sh build    # ~1 h on a Pi 5; idempotent
#   sudo bash pi_kernel_landlock.sh try      # reboot once into it (tryboot)
#   sudo bash pi_kernel_landlock.sh status   # what runs, what is installed
#   sudo bash pi_kernel_landlock.sh remove   # back to the stock kernel
# After an apt kernel upgrade `status` says so; run build and try again.
# ============================================================================

FIRMWARE=/boot/firmware
CONFIG_TXT="$FIRMWARE/config.txt"
# The build runs as root: in a directory of its own, made fresh (mktemp, 0700),
# never at a fixed name another user could have created first.
BUILD_PARENT="${DELFIN_KERNEL_BUILD_PARENT:-/var/tmp}"
SUFFIX="-landlock"
BEGIN_MARK="# BEGIN delfin-landlock (pi_kernel_landlock.sh; remove this block for the stock kernel)"
END_MARK="# END delfin-landlock"
SELF_COPY=/usr/local/sbin/delfin-landlock-kernel
UNIT=/etc/systemd/system/delfin-landlock-commit.service

die() { echo "error: $*" >&2; exit 1; }
landlock_active() { grep -qw landlock /sys/kernel/security/lsm 2>/dev/null; }

# The flavour of the running stock kernel line: 2712 (Pi 5) or v8 (Pi 4).
flavour() {
    local rel
    rel="$(uname -r)"
    rel="${rel%"$SUFFIX"}"
    case "$rel" in
        *-rpi-2712) echo 2712 ;;
        *-rpi-v8) echo v8 ;;
        *) return 1 ;;
    esac
}

# The config.txt section that holds only this board family.
board_section() {
    local model
    model="$(tr -d '\0' < /proc/device-tree/model 2>/dev/null || true)"
    case "$1:$model" in
        2712:"Raspberry Pi 5"*|2712:"Raspberry Pi Compute Module 5"*|2712:"Raspberry Pi 500"*) echo "[pi5]" ;;
        v8:"Raspberry Pi 4"*|v8:"Raspberry Pi 400"*|v8:"Raspberry Pi Compute Module 4"*) echo "[pi4]" ;;
        *) return 1 ;;
    esac
}

# Stock file names of the flavour, and ours beside them.
stock_kernel() { [[ "$1" == 2712 ]] && echo kernel_2712.img || echo kernel8.img; }
our_kernel() { [[ "$1" == 2712 ]] && echo kernel_2712_landlock.img || echo kernel8_landlock.img; }
# auto_initramfs: "kernel" -> "initramfs", extension dropped.
our_initramfs() { local k; k="$(our_kernel "$1")"; k="${k%.img}"; echo "initramfs${k#kernel}"; }

# The release the stock metapackage points at now (follows apt upgrades).
stock_release() {
    dpkg-query -W -f='${Depends}' "linux-image-rpi-$1" 2>/dev/null \
        | grep -o "linux-image-[^ ,]*-rpi-$1" | head -1 | sed 's/^linux-image-//'
}

block() {
    printf '%s\n%s\nkernel=%s\n[all]\n%s\n' "$BEGIN_MARK" "$1" "$2" "$END_MARK"
}

strip_block() {   # stdin -> stdout without our block
    awk -v b="$BEGIN_MARK" -v e="$END_MARK" '$0==b{skip=1;next} $0==e{skip=0;next} !skip'
}

require_pi() {
    [[ $EUID -eq 0 ]] || die "run with sudo"
    [[ -f "$CONFIG_TXT" ]] || die "$CONFIG_TXT not found: not Raspberry Pi OS"
    FLAVOUR="$(flavour)" || die "running kernel $(uname -r) is not a Raspberry Pi OS kernel (rpi-2712 / rpi-v8)"
    SECTION="$(board_section "$FLAVOUR")" \
        || die "this board ($(tr -d '\0' < /proc/device-tree/model 2>/dev/null)) is not one this script was made for"
    STOCK_REL="$(stock_release "$FLAVOUR")"
    [[ -n "$STOCK_REL" ]] || die "linux-image-rpi-$FLAVOUR is not installed"
    OUR_REL="$STOCK_REL$SUFFIX"
}

cmd_build() {
    require_pi
    if [[ -f "$FIRMWARE/$(our_kernel "$FLAVOUR")" && -d "/lib/modules/$OUR_REL" ]]; then
        echo "==> $OUR_REL is built and installed already"
        install_commit_unit
        return 0
    fi
    local cfg="/boot/config-$STOCK_REL"
    [[ -f "$cfg" ]] || die "$cfg missing"
    local pkgver
    pkgver="$(dpkg-query -W -f='${Version}' "linux-image-$STOCK_REL")"

    echo "==> build tools and linux-source $pkgver (signed Raspberry Pi archive)"
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -q
    apt-get install -y -q bc bison flex libssl-dev libelf-dev make gcc kmod cpio \
        xz-utils initramfs-tools "linux-source=$pkgver"
    # linux-source depends on linux-source-<major.minor>, which ships the tree
    # as /usr/src/linux-source-<major.minor>.tar.xz.
    local upstream="${pkgver#*:}"; upstream="${upstream%%-*}"
    local tarball
    tarball="$(dpkg -L "linux-source-${upstream%.*}" | grep -m1 '^/usr/src/linux-source-.*\.tar\.xz$')" \
        || die "no source tarball in linux-source-${upstream%.*}"

    BUILD_ROOT="$(mktemp -d "$BUILD_PARENT/delfin-landlock-kernel.XXXXXXXX")"
    trap 'rm -rf "$BUILD_ROOT"' EXIT
    tar -xJf "$tarball" -C "$BUILD_ROOT"
    local src
    src="$BUILD_ROOT/$(basename "$tarball" .tar.xz)"
    [[ -d "$src" ]] || die "$tarball did not unpack to $src"
    cp "$cfg" "$src/.config"
    ( cd "$src"
      scripts/config --enable SECURITY_LANDLOCK --set-str LSM "landlock"
      make -s olddefconfig
      grep -q '^CONFIG_SECURITY_LANDLOCK=y' .config || die "Landlock did not stay enabled"
      [[ "$(make -s KERNELRELEASE="$OUR_REL" kernelrelease)" == "$OUR_REL" ]] \
          || die "kernel release is not $OUR_REL"
      echo "==> building $OUR_REL with $(nproc) jobs (about an hour on a Pi 5)"
      nice -n 10 make -j"$(nproc)" KERNELRELEASE="$OUR_REL" Image modules
      make -s KERNELRELEASE="$OUR_REL" INSTALL_MOD_STRIP=1 modules_install
      install -m 0644 arch/arm64/boot/Image "$FIRMWARE/$(our_kernel "$FLAVOUR").new"
      install -m 0644 .config "/boot/config-$OUR_REL" )
    depmod "$OUR_REL"
    update-initramfs -c -k "$OUR_REL"
    install -m 0644 "/boot/initrd.img-$OUR_REL" "$FIRMWARE/$(our_initramfs "$FLAVOUR")"
    mv "$FIRMWARE/$(our_kernel "$FLAVOUR").new" "$FIRMWARE/$(our_kernel "$FLAVOUR")"
    rm -rf "$BUILD_ROOT"
    install_commit_unit
    echo "==> built. Boot it once with: sudo bash $0 try"
}

install_commit_unit() {
    install -m 0755 "$(readlink -f "$0")" "$SELF_COPY"
    cat > "$UNIT" <<EOF
[Unit]
Description=Keep the Landlock kernel once a tryboot of it came up with Landlock
After=local-fs.target

[Service]
Type=oneshot
ExecStart=$SELF_COPY commit

[Install]
WantedBy=multi-user.target
EOF
    systemctl daemon-reload
    systemctl enable -q delfin-landlock-commit.service
}

cmd_try() {
    require_pi
    [[ -f "$FIRMWARE/$(our_kernel "$FLAVOUR")" ]] || die "not built yet: sudo bash $0 build"
    # The stock cmdline with panic=10: a kernel that does not come up
    # reboots, and a reboot without the tryboot flag is the stock kernel.
    printf '%s panic=10\n' "$(tr -d '\n' < "$FIRMWARE/cmdline.txt")" > "$FIRMWARE/cmdline_landlock_try.txt"
    { strip_block < "$CONFIG_TXT"
      printf '%s\n%s\nkernel=%s\ncmdline=cmdline_landlock_try.txt\n[all]\n%s\n' \
          "$BEGIN_MARK" "$SECTION" "$(our_kernel "$FLAVOUR")" "$END_MARK"; } > "$FIRMWARE/tryboot.txt"
    sync
    echo "==> rebooting once into $(our_kernel "$FLAVOUR") (tryboot)"
    reboot '0 tryboot'
}

booted_by_tryboot() {
    local f=/proc/device-tree/chosen/bootloader/tryboot
    [[ -r "$f" ]] && [[ "$(od -An -tu4 --endian=big "$f" | tr -d ' ')" == "1" ]]
}

cmd_commit() {
    FLAVOUR="$(flavour)" || exit 0
    SECTION="$(board_section "$FLAVOUR")" || exit 0
    if landlock_active && [[ "$(uname -r)" == *"$SUFFIX" ]]; then
        if ! grep -qxF "$BEGIN_MARK" "$CONFIG_TXT"; then
            cp "$CONFIG_TXT" "$CONFIG_TXT.before-landlock"
            { strip_block < "$CONFIG_TXT"; block "$SECTION" "$(our_kernel "$FLAVOUR")"; } > "$CONFIG_TXT.new"
            mv "$CONFIG_TXT.new" "$CONFIG_TXT"
            echo "delfin-landlock: $(uname -r) runs with Landlock; kept in config.txt"
        fi
        rm -f "$FIRMWARE/tryboot.txt" "$FIRMWARE/cmdline_landlock_try.txt"
    elif booted_by_tryboot; then
        echo "delfin-landlock: tryboot came up without Landlock; config.txt left as it was" >&2
    fi
}

cmd_status() {
    local fl rel
    echo "running:  $(uname -r)  lsm=$(cat /sys/kernel/security/lsm 2>/dev/null || echo '?')"
    landlock_active && echo "landlock: active" || echo "landlock: NOT active"
    fl="$(flavour)" || { echo "not a Raspberry Pi OS kernel"; return 0; }
    rel="$(stock_release "$fl")"
    echo "stock:    $rel"
    if [[ -f "$FIRMWARE/$(our_kernel "$fl")" ]]; then
        echo "built:    $(our_kernel "$fl") ($(ls /boot/config-*"$SUFFIX" 2>/dev/null | sed 's#.*/config-##' | tr '\n' ' '))"
        [[ -d "/lib/modules/$rel$SUFFIX" ]] \
            || echo "outdated: the stock kernel is $rel now -- run build and try again"
    else
        echo "built:    no"
    fi
    grep -qxF "$BEGIN_MARK" "$CONFIG_TXT" 2>/dev/null && echo "boot:     config.txt boots the Landlock kernel" \
        || echo "boot:     config.txt boots the stock kernel"
}

cmd_remove() {
    [[ $EUID -eq 0 ]] || die "run with sudo"
    if [[ -f "$CONFIG_TXT" ]]; then
        strip_block < "$CONFIG_TXT" > "$CONFIG_TXT.new" && mv "$CONFIG_TXT.new" "$CONFIG_TXT"
    fi
    rm -f "$FIRMWARE/tryboot.txt" "$FIRMWARE/cmdline_landlock_try.txt" \
        "$FIRMWARE"/kernel*_landlock.img "$FIRMWARE"/initramfs*_landlock
    systemctl disable -q delfin-landlock-commit.service 2>/dev/null || true
    rm -f "$UNIT" "$SELF_COPY"
    systemctl daemon-reload
    echo "==> config.txt boots the stock kernel again (modules under /lib/modules/*$SUFFIX stay until removed by hand)"
}

case "${1:-}" in
    build) cmd_build ;;
    try) cmd_try ;;
    commit) cmd_commit ;;
    status) cmd_status ;;
    remove) cmd_remove ;;
    *) echo "usage: $0 {build|try|status|remove}" >&2; exit 2 ;;
esac
