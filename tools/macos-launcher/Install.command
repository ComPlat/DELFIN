#!/bin/bash
# DELFIN - installation. Double-click in Finder; Terminal opens and runs this.
cd "$(dirname "$0")" || exit 1
here=$(pwd -P)
echo 'DELFIN - installation'
echo
finish() { echo; printf 'Press Enter to close. '; read -r _; exit "$1"; }
if [ ! -f "$here/DELFIN.sh" ] || [ ! -f "$here/remote_bootstrap.sh" ]; then
    echo 'The files of this package are not next to Install.command.'
    echo 'Double-click the ZIP file in Finder first: it extracts into a folder.'
    echo 'Then double-click Install.command inside that folder. Nothing was installed.'
    finish 1
fi
# macOS marks every file that came from the internet (the quarantine
# attribute), and Gatekeeper then refuses to open an unsigned app or script
# with it -- each one, one at a time. Removing the mark is offered, never
# done silently: it is a protection, and lifting it for these files is the
# person's decision, not the installer's.
marked=''
for name in *; do
    if [ -f "$name" ] && xattr -p com.apple.quarantine "$name" >/dev/null 2>&1; then marked="$marked $name"; fi
done
if [ -n "$marked" ]; then
    echo 'macOS marked these files as downloaded from the internet:'
    for name in $marked; do echo "   $name"; done
    echo
    echo 'With the mark, macOS asks about DELFIN again each time it opens one of'
    echo 'them, and an unsigned one is only opened through System Settings.'
    echo 'Removing it clears the mark on THESE files and nothing else. No setting'
    echo 'is changed and Gatekeeper stays on. Do it only if you trust where this came from.'
    echo
    printf 'Remove the mark and continue? (y/N) '
    read -r answer
    case $answer in
        y|Y) for name in $marked; do xattr -d com.apple.quarantine "$name"; done; echo 'Removed.' ;;
        *) echo 'Nothing was changed. START-HERE.txt describes the other way'
           echo '(System Settings > Privacy & Security > Open Anyway, for each file).'
           finish 1 ;;
    esac
fi
echo 'Checking the package...'
/bin/bash "$here/DELFIN.sh" validate || { echo 'The package is incomplete or changed. Nothing was installed.'; finish 1; }
/bin/bash "$here/DELFIN.sh" install || { echo; echo 'Installation failed. Read the error above, then START-HERE.txt.'; finish 1; }
apps=${DELFIN_LAUNCHER_APPS:-"$HOME/Applications"}
echo
printf 'Start DELFIN now? (Y/n) '
read -r answer
case $answer in n|N) ;; *) open "$apps/DELFIN.app" ;; esac
finish 0
