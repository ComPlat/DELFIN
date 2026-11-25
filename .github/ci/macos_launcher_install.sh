#!/bin/bash
# The macOS launcher as a person gets it: the release ZIP, marked as a
# browser download, extracted the way Finder does, Install.command answered
# both ways at the download-mark question, then uninstall keeping and
# removing. Runs on macos-latest (.github/workflows/macos-launcher.yml)
# with the Mac's own /bin/bash 3.2.
set -eo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd -P)
work=$(mktemp -d)
home="$work/home"
mkdir -p "$home"
export DELFIN_LAUNCHER_STORE="$home/store" DELFIN_LAUNCHER_APPS="$home/Applications" DELFIN_LAUNCHER_SSH_DIR="$home/.ssh"
fail() { echo "FAIL: $*" >&2; exit 1; }
mark='0083;67000000;Safari;'

# The committed ZIP, as the release publishes it (that it matches this tree
# is tests/test_the_macos_launcher_package_is_current.py).
cp "$repo/tools/delfin-macos-launcher.zip" "$work/"
xattr -w com.apple.quarantine "$mark" "$work/delfin-macos-launcher.zip"
# ditto is Archive Utility's engine; Archive Utility also passes the mark on
# to every extracted file, which is done here by hand.
ditto -x -k "$work/delfin-macos-launcher.zip" "$work/x"
pkg="$work/x/macos-launcher"
for file in "$pkg"/*; do xattr -w com.apple.quarantine "$mark" "$file"; done
[ -x "$pkg/Install.command" ] && [ -x "$pkg/Uninstall.command" ] || fail 'the .command files are not executable after extraction'

echo '== answer n at the download-mark question'
if printf 'n\n\n' | /bin/bash "$pkg/Install.command" > "$work/no.log" 2>&1; then cat "$work/no.log"; fail 'n installed anyway'; fi
cat "$work/no.log"
grep -q 'Remove the mark and continue? (y/N)' "$work/no.log" || fail 'no question'
grep -q 'Install.command' "$work/no.log" || fail 'the marked files are not listed'
[ ! -e "$home/Applications/DELFIN.app" ] || fail 'n left an app behind'
xattr -p com.apple.quarantine "$pkg/DELFIN.sh" >/dev/null || fail 'n removed the mark'

echo '== answer y, then do not start'
printf 'y\nn\n\n' | /bin/bash "$pkg/Install.command" > "$work/yes.log" 2>&1 || { cat "$work/yes.log"; fail 'install failed'; }
cat "$work/yes.log"
for file in "$pkg"/*; do
    if xattr -p com.apple.quarantine "$file" >/dev/null 2>&1; then fail "mark still on $file"; fi
done
app="$home/Applications/DELFIN.app/Contents"
plutil -lint "$app/Info.plist"
[ "$(plutil -extract CFBundleExecutable raw "$app/Info.plist")" = DELFIN ] || fail 'CFBundleExecutable'
[ -x "$app/MacOS/DELFIN" ] || fail 'the app has no executable'
[ -s "$app/Resources/DELFIN.icns" ] || fail 'no icon (sips/iconutil)'
file "$app/Resources/DELFIN.icns"
/bin/bash "$app/Resources/DELFIN.sh" validate "$app/Resources"
[ "$(stat -f %Lp "$home/store")" = 700 ] || fail 'the settings folder is not 0700'
if xattr -p com.apple.quarantine "$app/Resources/DELFIN.sh" >/dev/null 2>&1; then fail 'the installed app carries a mark the person removed'; fi

echo '== the menu lists what it should (osascript replaced by a recorder)'
# The dialogs themselves need a person; what the app offers is checked here.
stub="$work/stub"
mkdir -p "$stub"
printf '#!/bin/bash\ncat > /dev/null\nprintf "%%s\\n" "$@" > %q\n' "$work/menu.txt" > "$stub/osascript"
chmod 755 "$stub/osascript"
sed "s#^PATH=/usr/bin:/bin:/usr/sbin:/sbin#PATH=$stub:/usr/bin:/bin:/usr/sbin:/sbin#" "$app/Resources/DELFIN.sh" > "$work/menu.sh"
cp "$app/Resources/"{remote_bootstrap.sh,remote_launcher.py} "$work/"
/bin/bash "$work/menu.sh" app
cat "$work/menu.txt"
grep -qx 'Set up DELFIN on this Mac...' "$work/menu.txt" || fail 'the menu does not offer DELFIN on this Mac'
grep -qx 'New server connection...' "$work/menu.txt" || fail 'the menu does not offer a server'

echo '== uninstall keeping the connections'
mkdir -p "$home/store/profiles" && echo 'target=delfin@192.0.2.10' > "$home/store/profiles/0b8f1f7e-8b1c-4a4f-9f1e-3c5d7e9a1b2c.conf"
printf '\n\n' | /bin/bash "$pkg/Uninstall.command" > "$work/u1.log" 2>&1 || { cat "$work/u1.log"; fail 'uninstall'; }
cat "$work/u1.log"
[ ! -e "$home/Applications/DELFIN.app" ] || fail 'the app is still there'
[ -f "$home/store/profiles/0b8f1f7e-8b1c-4a4f-9f1e-3c5d7e9a1b2c.conf" ] || fail 'Y (the default) did not keep the connections'

echo '== install again, uninstall removing everything'
printf 'n\n\n' | /bin/bash "$pkg/Install.command" > "$work/again.log" 2>&1 || { cat "$work/again.log"; fail 'reinstall'; }
printf 'n\n\n' | /bin/bash "$pkg/Uninstall.command" > "$work/u2.log" 2>&1 || { cat "$work/u2.log"; fail 'uninstall 2'; }
cat "$work/u2.log"
[ ! -e "$home/store" ] || { ls -la "$home/store"; fail 'n did not remove the connections'; }
echo 'OK: install and uninstall on macOS'
