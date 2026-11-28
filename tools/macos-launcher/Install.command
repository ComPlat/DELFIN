#!/bin/bash
# DELFIN - installation. Double-click in Finder (Terminal opens and runs
# this), or run it with bash: the one-line install in README.md does.
cd "$(dirname "$0")" || exit 1
here=$(pwd -P)
app="$here/DELFIN.app"
script="$app/Contents/Resources/DELFIN.sh"
echo 'DELFIN - installation'
echo
finish() { echo; printf 'Press Enter to close. '; read -r _; exit "$1"; }
if [ ! -f "$script" ] || [ ! -f "$here/SHA256SUMS.txt" ]; then
    echo 'The files of this package are not next to Install.command.'
    echo 'Double-click the ZIP file in Finder first: it extracts into a folder.'
    echo 'Then double-click Install.command inside that folder. Nothing was installed.'
    finish 1
fi
# macOS marks every file that came from the internet (the quarantine
# attribute), and Gatekeeper then refuses to open software that Apple has
# not notarized. A notarized DELFIN.app opens with the mark on, and the mark
# stays. For an unsigned build, removing the mark is offered, never done
# silently: it is a protection, and lifting it for these files is the
# person's decision, not the installer's. Gatekeeper itself stays on.
if xattr -p com.apple.quarantine "$app" >/dev/null 2>&1 && spctl --assess --type execute "$app" >/dev/null 2>&1; then
    echo 'DELFIN.app is signed and notarized; macOS opens it with the download mark on.'
else
    marked=''
    for name in *; do
        if [ -n "$(find "$name" -xattrname com.apple.quarantine -print -quit 2>/dev/null)" ]; then marked="$marked $name"; fi
    done
    if [ -n "$marked" ]; then
        echo 'macOS marked these files as downloaded from the internet:'
        for name in $marked; do echo "   $name"; done
        echo
        echo 'This build of DELFIN is not notarized by Apple. With the mark, macOS'
        echo 'refuses to open DELFIN.app until you allow it in System Settings.'
        echo 'Removing the mark clears it on THESE files and nothing else. No setting'
        echo 'is changed and Gatekeeper stays on. Do it only if you trust where this came from.'
        echo
        printf 'Remove the mark and continue? (y/N) '
        read -r answer
        case $answer in
            y|Y) for name in $marked; do xattr -dr com.apple.quarantine "$name"; done; echo 'Removed.' ;;
            *) echo 'Nothing was changed. START-HERE.txt describes the other way'
               echo '(System Settings > Privacy & Security > Open Anyway).'
               finish 1 ;;
        esac
    fi
fi
echo 'Checking the package...'
/bin/bash "$script" validate "$here" || { echo 'The package is incomplete or changed. Nothing was installed.'; finish 1; }
/bin/bash "$script" install || { echo; echo 'Installation failed. Read the error above, then START-HERE.txt.'; finish 1; }
apps=${DELFIN_LAUNCHER_APPS:-"$HOME/Applications"}
echo
printf 'Start DELFIN now? (Y/n) '
read -r answer
case $answer in n|N) ;; *) open "$apps/DELFIN.app" ;; esac
finish 0
