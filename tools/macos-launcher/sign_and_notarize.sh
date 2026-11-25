#!/bin/bash
# Sign and notarize the macOS launcher for a release (maintainer-side, on a Mac
# with Xcode's command line tools and an Apple Developer ID). Unset, nothing
# here runs and the release ships the unsigned package, which works the same
# but asks the person to allow it once (Install.command, START-HERE.txt).
#
#   DELFIN_SIGN_IDENTITY   "Developer ID Application: <Team name> (<TEAMID>)",
#                          a certificate in the login keychain
#   DELFIN_NOTARY_PROFILE  a notarytool keychain profile, made once with
#                          xcrun notarytool store-credentials <profile>
#                          (Apple ID + app-specific password, or an API key);
#                          nothing of it is in the repository
#
#   bash tools/macos-launcher/sign_and_notarize.sh
#
# Steps, all Apple's own tools:
#   1. build_package.py's ZIP is extracted with ditto.
#   2. DELFIN.app is signed with the hardened runtime and a secure timestamp
#      (codesign --options runtime --timestamp), then verified.
#   3. The app is submitted to Apple's notary service (notarytool --wait) and
#      the ticket stapled to it (stapler), so Gatekeeper accepts it offline.
#   4. The ZIP is rebuilt from the signed tree with ditto, and spctl must
#      accept the app. The signature adds Contents/_CodeSignature and the
#      ticket Contents/CodeResources; neither is in SHA256SUMS.txt, and no
#      listed file changes, so Install.command's check still holds.
#   5. tools/delfin-macos-launcher.zip is replaced. Then update the
#      one-line install's digest (README.md, between the macos-install-line
#      markers) to the printed SHA-256 and commit both before tagging;
#      tests/test_the_macos_launcher_package_is_current.py checks the ZIP's
#      files and the line's digest.
#
# Install.command and Uninstall.command stay outside the app: a plain script
# cannot carry a stapled ticket. With a notarized app, Install.command does not
# ask about the download mark (Gatekeeper accepts the app with it), and the
# one-line install never sets the mark (curl does not). A signed installer
# package (pkgbuild/productbuild + productsign, "Developer ID Installer")
# would remove the remaining first-open question for Install.command in
# Finder; it is not built here.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd -P)
zip="$here/../delfin-macos-launcher.zip"
: "${DELFIN_SIGN_IDENTITY:?set DELFIN_SIGN_IDENTITY (Developer ID Application)}"
: "${DELFIN_NOTARY_PROFILE:?set DELFIN_NOTARY_PROFILE (xcrun notarytool store-credentials)}"
[ "$(uname -s)" = Darwin ] || { echo 'Signing needs macOS.' >&2; exit 1; }
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
ditto -x -k "$zip" "$work"
app="$work/macos-launcher/DELFIN.app"
codesign --force --options runtime --timestamp --sign "$DELFIN_SIGN_IDENTITY" "$app"
codesign --verify --deep --strict --verbose=2 "$app"
ditto -c -k --keepParent "$app" "$work/notarize.zip"
xcrun notarytool submit "$work/notarize.zip" --keychain-profile "$DELFIN_NOTARY_PROFILE" --wait
xcrun stapler staple "$app"
xcrun stapler validate "$app"
spctl --assess --type execute --verbose=2 "$app"
(cd "$work/macos-launcher" && shasum -a 256 -c SHA256SUMS.txt >/dev/null)
rm -f "$work/signed.zip"
(cd "$work" && ditto -c -k --keepParent macos-launcher signed.zip)
cp "$work/signed.zip" "$zip"
echo "Signed and notarized: $zip"
echo "SHA-256 for the one-line install: $(shasum -a 256 "$zip" | awk '{print $1}')"
