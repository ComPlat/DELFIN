#!/bin/bash
# The macOS launcher as a person gets it, on macos-latest with the Mac's own
# /bin/bash 3.2 (.github/workflows/macos-launcher.yml):
#  - the one-line install from README.md, with the committed ZIP in place of
#    the release download (the asset exists only once a release carries it):
#    its SHA-256 check, extraction and Install.command;
#  - the ZIP marked as a browser download, extracted as Finder does, and
#    Install.command answered n and y at the download-mark question;
#  - the app bundle: plutil, icon, ad-hoc signature verified strictly (what
#    a Developer ID signature needs from the bundle's structure);
#  - uninstall keeping and removing.
set -eo pipefail
repo=$(cd "$(dirname "$0")/../.." && pwd -P)
work=$(mktemp -d)
home="$work/home"
mkdir -p "$home"
export DELFIN_LAUNCHER_STORE="$home/store" DELFIN_LAUNCHER_APPS="$home/Applications" DELFIN_LAUNCHER_SSH_DIR="$home/.ssh"
fail() { echo "FAIL: $*" >&2; exit 1; }
mark='0083;67000000;Safari;'
zip="$repo/tools/delfin-macos-launcher.zip"

echo '== the one-line install (download replaced by the committed ZIP)'
line=$(sed -n '/^<!-- macos-install-line/,/^<!-- \/macos-install-line/p' "$repo/README.md" | grep '^d=\$(mktemp -d)')
[ -n "$line" ] || fail 'no one-line install in README.md'
case $line in *'curl -fL --proto =https -o "$d/delfin.zip" https://github.com/ComPlat/DELFIN/releases/download/v'*) ;; *) fail "unexpected line: $line" ;; esac
local_line=$(printf '%s' "$line" | sed -E "s#curl -fL --proto =https -o \"\\\$d/delfin.zip\" [^ ]+#cp '$zip' \"\$d/delfin.zip\"#")
echo "$local_line"
printf 'n\n\n' | /bin/bash -c "$local_line" > "$work/line.log" 2>&1 || { cat "$work/line.log"; fail 'the one-line install failed'; }
cat "$work/line.log"
grep -q ': OK$' "$work/line.log" || fail 'the SHA-256 check did not run'
grep -q 'download mark' "$work/line.log" && fail 'curl sets no mark; nothing should be asked'
[ -d "$home/Applications/DELFIN.app" ] || fail 'the one-line install installed no app'
printf '\n\n' | /bin/bash "$home/Applications/DELFIN.app/Contents/Resources/DELFIN.sh" uninstall --remove-profiles --remove-local > /dev/null
[ ! -e "$home/Applications/DELFIN.app" ] || fail 'uninstall after the one-line install'
broken=$(printf '%s' "$local_line" | sed -E 's/echo "[0-9a-f]{64}/echo "0000000000000000000000000000000000000000000000000000000000000000/')
if printf 'n\n\n' | /bin/bash -c "$broken" > "$work/broken.log" 2>&1; then fail 'a wrong digest did not stop the install'; fi
[ ! -e "$home/Applications/DELFIN.app" ] || fail 'a wrong digest installed anyway'
echo 'one-line install: digest checked before anything ran; a wrong digest stops it'

echo '== downloaded with a browser: the mark, and Finder extracting it'
cp "$zip" "$work/"
xattr -w com.apple.quarantine "$mark" "$work/delfin-macos-launcher.zip"
# ditto is Archive Utility's engine; Archive Utility also passes the mark on
# to every extracted file, which is done here by hand.
ditto -x -k "$work/delfin-macos-launcher.zip" "$work/x"
pkg="$work/x/macos-launcher"
find "$pkg" -type f -exec xattr -w com.apple.quarantine "$mark" {} \;
xattr -w com.apple.quarantine "$mark" "$pkg/DELFIN.app"
[ -x "$pkg/Install.command" ] && [ -x "$pkg/Uninstall.command" ] && [ -x "$pkg/DELFIN.app/Contents/MacOS/DELFIN" ] \
    || fail 'executables lost their mode in extraction'
if spctl --assess --type execute "$pkg/DELFIN.app" 2>/dev/null; then fail 'an unsigned app passed Gatekeeper?'; fi

echo '== answer n at the download-mark question'
if printf 'n\n\n' | /bin/bash "$pkg/Install.command" > "$work/no.log" 2>&1; then cat "$work/no.log"; fail 'n installed anyway'; fi
cat "$work/no.log"
grep -q 'Remove the mark and continue? (y/N)' "$work/no.log" || fail 'no question'
grep -q '   DELFIN.app' "$work/no.log" || fail 'the marked app is not listed'
[ ! -e "$home/Applications/DELFIN.app" ] || fail 'n left an app behind'
xattr -p com.apple.quarantine "$pkg/DELFIN.app/Contents/Resources/DELFIN.sh" >/dev/null || fail 'n removed the mark'

echo '== answer y, then do not start'
printf 'y\nn\n\n' | /bin/bash "$pkg/Install.command" > "$work/yes.log" 2>&1 || { cat "$work/yes.log"; fail 'install failed'; }
cat "$work/yes.log"
[ -z "$(find "$pkg" -xattrname com.apple.quarantine)" ] || fail 'mark still on some files'
app="$home/Applications/DELFIN.app"
plutil -lint "$app/Contents/Info.plist"
[ "$(plutil -extract CFBundleExecutable raw "$app/Contents/Info.plist")" = DELFIN ] || fail 'CFBundleExecutable'
[ -x "$app/Contents/MacOS/DELFIN" ] || fail 'the app has no executable'
file "$app/Contents/Resources/DELFIN.icns" | grep -qi 'icon' || fail 'DELFIN.icns is no Mac icon'
diff -r "$pkg/DELFIN.app" "$app" || fail 'the installed app differs from the package'
[ "$(stat -f %Lp "$home/store")" = 700 ] || fail 'the settings folder is not 0700'

echo '== the bundle can be signed as it is'
ditto "$app" "$work/signed/DELFIN.app"
codesign --force --sign - --options runtime "$work/signed/DELFIN.app"
codesign --verify --deep --strict --verbose=2 "$work/signed/DELFIN.app"
(cd "$pkg" && shasum -a 256 -c SHA256SUMS.txt >/dev/null) || fail 'checksums'
for name in $(awk '{print $2}' "$pkg/SHA256SUMS.txt" | grep '^DELFIN.app/'); do
    cmp -s "$pkg/$name" "$work/signed/$name" || fail "signing changed $name (Install.command's check would fail)"
done
echo 'ad-hoc signature valid; no checked file changed by signing'

echo '== the menu lists what it should (osascript replaced by a recorder)'
stub="$work/stub"
mkdir -p "$stub"
printf '#!/bin/bash\ncat > /dev/null\nprintf "%%s\\n" "$@" > %q\n' "$work/menu.txt" > "$stub/osascript"
chmod 755 "$stub/osascript"
mkdir -p "$work/menu/Contents/Resources"
sed "s#^PATH=/usr/bin:/bin:/usr/sbin:/sbin#PATH=$stub:/usr/bin:/bin:/usr/sbin:/sbin#" "$app/Contents/Resources/DELFIN.sh" > "$work/menu/Contents/Resources/DELFIN.sh"
/bin/bash "$work/menu/Contents/Resources/DELFIN.sh" app
cat "$work/menu.txt"
grep -qx 'Set up DELFIN on this Mac...' "$work/menu.txt" || fail 'the menu does not offer DELFIN on this Mac'
grep -qx 'New server connection...' "$work/menu.txt" || fail 'the menu does not offer a server'

echo '== uninstall keeping the connections'
mkdir -p "$home/store/profiles" && echo 'target=delfin@192.0.2.10' > "$home/store/profiles/0b8f1f7e-8b1c-4a4f-9f1e-3c5d7e9a1b2c.conf"
printf '\n\n' | /bin/bash "$pkg/Uninstall.command" > "$work/u1.log" 2>&1 || { cat "$work/u1.log"; fail 'uninstall'; }
cat "$work/u1.log"
[ ! -e "$app" ] || fail 'the app is still there'
[ -f "$home/store/profiles/0b8f1f7e-8b1c-4a4f-9f1e-3c5d7e9a1b2c.conf" ] || fail 'Y (the default) did not keep the connections'

echo '== install again, uninstall removing everything'
printf 'n\n\n' | /bin/bash "$pkg/Install.command" > "$work/again.log" 2>&1 || { cat "$work/again.log"; fail 'reinstall'; }
printf 'n\n\n' | /bin/bash "$pkg/Uninstall.command" > "$work/u2.log" 2>&1 || { cat "$work/u2.log"; fail 'uninstall 2'; }
cat "$work/u2.log"
[ ! -e "$home/store" ] || { ls -la "$home/store"; fail 'n did not remove the connections'; }
ls "$home/Applications" | grep -q . && fail 'something is left in Applications'
echo 'OK: install and uninstall on macOS'
