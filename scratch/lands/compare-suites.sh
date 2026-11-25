#!/usr/bin/env bash
# Compare the bundle's full-suite failures with the origin/main control.
#
# A node that fails in BOTH logs was already red before this bundle; only a
# node that PASSES on the control and FAILS on the bundle is a regression of
# it. No interpreter path is hardcoded: pass it as $3 if the default is wrong,
# because a committed script must not carry a home path.
#
#   ./compare-suites.sh scratch/lands/full-suite-bundle.log \
#                       scratch/lands/full-suite-control.log
set -uo pipefail

bundle="${1:-scratch/lands/full-suite-bundle.log}"
control="${2:-scratch/lands/full-suite-control.log}"

for f in "$bundle" "$control"; do
    if [ ! -f "$f" ]; then
        echo "missing log: $f" >&2
        exit 2
    fi
done

# pytest -q writes "FAILED <nodeid> - <reason>" lines in its short summary.
# ERROR lines are collected too: a fixture that blows up is a failure.
#
# The scratch files sit beside the logs rather than in the system temporary
# folder: this harness confines writes to the workspace, and a path that a
# previous run cannot open is a comparison that silently produces nothing.
scratch_dir="$(dirname "$bundle")"
[ -d "$scratch_dir" ] || { echo "no such directory: $scratch_dir" >&2; exit 2; }
failures() {
    grep -oE '^(FAILED|ERROR) [^ ]+' "$1" | awk '{print $2}' | sort -u
}

failures "$bundle"  > "$scratch_dir/.compare-bundle.fails"
failures "$control" > "$scratch_dir/.compare-control.fails"
bundle_fails="$scratch_dir/.compare-bundle.fails"
control_fails="$scratch_dir/.compare-control.fails"

nb=$(wc -l < "$bundle_fails")
nc=$(wc -l < "$control_fails")

echo "bundle failures:  $nb"
echo "control failures: $nc"
echo
echo "REGRESSIONS (fail on the bundle, pass on the control):"
comm -23 "$bundle_fails" "$control_fails"
echo
echo "ALREADY RED (fail on both -- not this bundle's doing):"
comm -12 "$bundle_fails" "$control_fails" | head -50
comm -12 "$bundle_fails" "$control_fails" | wc -l \
    | xargs -r echo "already-red total (untruncated):"
echo
echo "FIXED BY THE BUNDLE (fail on the control, pass here):"
comm -13 "$bundle_fails" "$control_fails"

rm -f "$bundle_fails" "$control_fails"

# The exit codes the two runs reported, so a truncated summary cannot hide a
# run that died rather than finished.
echo
for f in "$bundle" "$control"; do
    printf '%s: %s\n' "$(basename "$f")" \
        "$(grep -o 'PYTEST_EXIT=[0-9]*' "$f" | tail -1 || echo 'no exit line yet')"
done
