#!/bin/bash
# DELFIN - uninstall. Double-click in Finder; Terminal opens and runs this.
# Asks whether to keep the saved connections, and whether to keep DELFIN
# set up on this Mac. Nothing on any server is touched.
apps=${DELFIN_LAUNCHER_APPS:-"$HOME/Applications"}
script="$(cd "$(dirname "$0")" && pwd -P)/DELFIN.app/Contents/Resources/DELFIN.sh"
[ -f "$script" ] || script="$apps/DELFIN.app/Contents/Resources/DELFIN.sh"
cd / || exit 1
echo 'DELFIN - uninstall'
echo
if [ ! -f "$script" ]; then
    echo 'DELFIN is not installed for this account.'
else
    /bin/bash "$script" uninstall "$@"
fi
echo
printf 'Press Enter to close. '
read -r _
